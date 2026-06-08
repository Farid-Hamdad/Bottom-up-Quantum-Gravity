#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS — add dimensional sector v1

Adds to the SLACS CSV:
    ds_mean
    dw_mean
    alpha_eff_mean

with:
    alpha_eff = 2 ds / dw + dw - 4

Graph:
    polar grid
    Sersic profile
    W_ij = sqrt(sigma_i sigma_j) exp(-d_ij^2 / 2 xi^2)

Spectral dimension:
    P(t) = (1/N) Tr exp(-t L_norm)
    ds = -2 d log P / d log t

Walk dimension:
    <r^2(t)> ~ t^(2/dw)
    dw = 2 / slope(log MSD vs log t)

This script does NOT force alpha_eff = 1.
It only computes the graph-derived dimensional sector.
"""

import argparse
import json
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


def read_csv_auto(path):
    return pd.read_csv(path, sep=None, engine="python")


def sersic_b_n(n):
    n = float(max(n, 0.3))
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n * n)


def sersic_sigma(R, Re, n):
    R = np.asarray(R, dtype=float)
    Re = max(float(Re), 1e-12)
    n = max(float(n), 0.3)
    b = sersic_b_n(n)
    x = np.maximum(R / Re, 1e-12)
    return np.exp(-b * (x ** (1.0 / n) - 1.0))


def get_sersic_n(row, default=4.0):
    for c in ["sersic_n_final", "sersic_n_measured", "n_sersic", "sersic_n", "n"]:
        if c in row and pd.notna(row[c]):
            try:
                val = float(row[c])
                if np.isfinite(val) and val > 0:
                    return val, c
            except Exception:
                pass
    return float(default), "default"


def get_q(row):
    for c in ["q_SIE", "q", "axis_ratio"]:
        if c in row and pd.notna(row[c]):
            try:
                q = float(row[c])
                if np.isfinite(q) and 0.1 <= q <= 1.5:
                    return q
            except Exception:
                pass
    return 1.0


def get_mstar_rel(row):
    if "Mstar_rel" in row and pd.notna(row["Mstar_rel"]):
        return float(row["Mstar_rel"])
    if "Mstar_rel_used" in row and pd.notna(row["Mstar_rel_used"]):
        return float(row["Mstar_rel_used"])
    if "logMs" in row and pd.notna(row["logMs"]):
        return float(10 ** (float(row["logMs"]) - 11.0))
    return 1.0


def make_polar_nodes(Re, q, n_rings, n_theta, rmin_factor, rmax_factor):
    r = np.logspace(np.log10(rmin_factor * Re), np.log10(rmax_factor * Re), n_rings)
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)

    nodes = []
    for rr in r:
        for tt in th:
            x = rr * np.cos(tt)
            y = rr * q * np.sin(tt)
            Rell = np.sqrt(x * x + (y / max(q, 1e-6)) ** 2)
            nodes.append((x, y, Rell, rr, tt))

    return np.array(nodes, dtype=float)


def build_W(nodes, sigma_amp, xi):
    x = nodes[:, 0]
    y = nodes[:, 1]

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    d2 = dx * dx + dy * dy

    sigma_amp = np.maximum(np.asarray(sigma_amp, dtype=float), 1e-300)

    W = np.sqrt(sigma_amp[:, None] * sigma_amp[None, :]) * np.exp(
        -d2 / max(2.0 * xi * xi, 1e-12)
    )

    np.fill_diagonal(W, 0.0)
    return W


def laplacian_norm(W):
    d = np.sum(W, axis=1)
    invsqrt = 1.0 / np.sqrt(np.maximum(d, 1e-12))
    return np.eye(W.shape[0]) - invsqrt[:, None] * W * invsqrt[None, :]


def transition_matrix(W):
    d = np.sum(W, axis=1)
    return W / np.maximum(d[:, None], 1e-12)


def pairwise_dist2(nodes):
    x = nodes[:, 0]
    y = nodes[:, 1]
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return dx * dx + dy * dy


def fit_slope_loglog(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.sum(m) < 4:
        return np.nan, np.nan

    X = np.log(x[m]).reshape(-1, 1)
    Y = np.log(y[m])

    model = LinearRegression()
    model.fit(X, Y)

    slope = float(model.coef_[0])
    r2 = float(model.score(X, Y))

    return slope, r2


def compute_spectral_dimension(L, t_min=0.05, t_max=20.0, n_t=80, p_low=0.05, p_high=0.85):
    vals = np.linalg.eigvalsh(L)
    vals = np.maximum(vals, 0.0)

    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)
    P = np.array([np.mean(np.exp(-tt * vals)) for tt in t], dtype=float)

    mask = np.isfinite(P) & (P > p_low) & (P < p_high)

    if np.sum(mask) < 6:
        # fallback middle range
        i0 = int(0.20 * len(t))
        i1 = int(0.75 * len(t))
        mask = np.zeros_like(P, dtype=bool)
        mask[i0:i1] = True

    slope, r2 = fit_slope_loglog(t[mask], P[mask])

    ds = -2.0 * slope if np.isfinite(slope) else np.nan

    return {
        "ds_mean": float(ds) if np.isfinite(ds) else np.nan,
        "ds_slope": float(slope) if np.isfinite(slope) else np.nan,
        "ds_fit_r2": float(r2) if np.isfinite(r2) else np.nan,
        "ds_n_fit": int(np.sum(mask)),
    }


def compute_walk_dimension(W, nodes, t_min_step=1, t_max_step=25):
    P = transition_matrix(W)
    D2 = pairwise_dist2(nodes)

    N = W.shape[0]
    Pt = np.eye(N)

    steps = []
    msd = []

    for step in range(1, t_max_step + 1):
        Pt = Pt @ P

        if step >= t_min_step:
            # Uniform start over all nodes
            val = np.mean(np.sum(Pt * D2, axis=1))
            if np.isfinite(val) and val > 0:
                steps.append(step)
                msd.append(val)

    steps = np.asarray(steps, dtype=float)
    msd = np.asarray(msd, dtype=float)

    slope, r2 = fit_slope_loglog(steps, msd)

    if np.isfinite(slope) and abs(slope) > 1e-12:
        dw = 2.0 / slope
    else:
        dw = np.nan

    return {
        "dw_mean": float(dw) if np.isfinite(dw) else np.nan,
        "dw_slope": float(slope) if np.isfinite(slope) else np.nan,
        "dw_fit_r2": float(r2) if np.isfinite(r2) else np.nan,
        "dw_n_fit": int(len(steps)),
    }


def compute_one(row, args):
    Re = float(row["Re_arcsec"])
    q = get_q(row)
    nser, n_source = get_sersic_n(row, args.default_sersic_n)

    nodes = make_polar_nodes(
        Re=Re,
        q=q,
        n_rings=args.n_rings,
        n_theta=args.n_theta,
        rmin_factor=args.rmin_factor,
        rmax_factor=args.rmax_factor,
    )

    R = nodes[:, 2]

    sigma = sersic_sigma(R, Re, nser)
    sigma = sigma / np.maximum(np.sum(sigma), 1e-12)

    if args.amp_mode == "shape":
        amp = 1.0
    elif args.amp_mode == "sqrtM":
        amp = np.sqrt(max(get_mstar_rel(row), 1e-12))
    elif args.amp_mode == "Mstar":
        amp = max(get_mstar_rel(row), 1e-12)
    else:
        raise ValueError(f"Unknown amp_mode: {args.amp_mode}")

    sigma_amp = sigma * amp

    xi = args.xi_factor * Re

    W = build_W(nodes, sigma_amp, xi)

    L = laplacian_norm(W)

    spec = compute_spectral_dimension(
        L,
        t_min=args.t_min,
        t_max=args.t_max,
        n_t=args.n_t,
        p_low=args.p_low,
        p_high=args.p_high,
    )

    walk = compute_walk_dimension(
        W,
        nodes,
        t_min_step=args.rw_t_min,
        t_max_step=args.rw_t_max,
    )

    ds = spec["ds_mean"]
    dw = walk["dw_mean"]

    if np.isfinite(ds) and np.isfinite(dw) and abs(dw) > 1e-12:
        alpha = 2.0 * ds / dw + dw - 4.0
    else:
        alpha = np.nan

    out = {
        "ds_mean": ds,
        "dw_mean": dw,
        "alpha_eff_mean": float(alpha) if np.isfinite(alpha) else np.nan,

        "ds_fit_r2": spec["ds_fit_r2"],
        "dw_fit_r2": walk["dw_fit_r2"],
        "ds_n_fit": spec["ds_n_fit"],
        "dw_n_fit": walk["dw_n_fit"],

        "sersic_n_used_dim": nser,
        "sersic_n_source_dim": n_source,
        "q_used_dim": q,
        "xi_used_dim": xi,
        "amp_mode_dim": args.amp_mode,
        "amp_used_dim": amp,
    }

    return out


def summarize(df):
    cols = ["ds_mean", "dw_mean", "alpha_eff_mean", "ds_fit_r2", "dw_fit_r2"]
    s = {}
    for c in cols:
        if c in df.columns:
            x = pd.to_numeric(df[c], errors="coerce")
            s[c] = {
                "n": int(x.notna().sum()),
                "mean": float(x.mean()),
                "median": float(x.median()),
                "std": float(x.std()),
                "min": float(x.min()),
                "max": float(x.max()),
            }

    if "logMs" in df.columns and "alpha_eff_mean" in df.columns:
        x = pd.to_numeric(df["logMs"], errors="coerce")
        y = pd.to_numeric(df["alpha_eff_mean"], errors="coerce")
        m = x.notna() & y.notna()
        if m.sum() > 3:
            s["alpha_vs_logMs"] = {
                "pearson": float(np.corrcoef(x[m], y[m])[0, 1]),
                "n": int(m.sum()),
            }

    return s


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--output-summary", required=True)

    ap.add_argument("--n-rings", type=int, default=12)
    ap.add_argument("--n-theta", type=int, default=16)
    ap.add_argument("--rmin-factor", type=float, default=0.03)
    ap.add_argument("--rmax-factor", type=float, default=5.0)
    ap.add_argument("--xi-factor", type=float, default=0.35)
    ap.add_argument("--default-sersic-n", type=float, default=4.0)

    ap.add_argument("--amp-mode", choices=["shape", "sqrtM", "Mstar"], default="shape")

    ap.add_argument("--t-min", type=float, default=0.05)
    ap.add_argument("--t-max", type=float, default=20.0)
    ap.add_argument("--n-t", type=int, default=80)
    ap.add_argument("--p-low", type=float, default=0.05)
    ap.add_argument("--p-high", type=float, default=0.85)

    ap.add_argument("--rw-t-min", type=int, default=1)
    ap.add_argument("--rw-t-max", type=int, default=25)

    args = ap.parse_args()

    df = read_csv_auto(args.input_csv)

    required = ["name", "Re_arcsec", "logMs"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    rows = []

    print("Computing ds_mean, dw_mean, alpha_eff_mean...")
    for k, (_, row) in enumerate(df.iterrows(), start=1):
        name = row["name"] if "name" in row else f"gal_{k:03d}"
        print(f"[{k}/{len(df)}] {name}")

        try:
            out = compute_one(row, args)
        except Exception as e:
            print(f"  FAILED: {e}")
            out = {
                "ds_mean": np.nan,
                "dw_mean": np.nan,
                "alpha_eff_mean": np.nan,
                "ds_fit_r2": np.nan,
                "dw_fit_r2": np.nan,
                "ds_n_fit": 0,
                "dw_n_fit": 0,
                "sersic_n_used_dim": np.nan,
                "sersic_n_source_dim": "failed",
                "q_used_dim": np.nan,
                "xi_used_dim": np.nan,
                "amp_mode_dim": args.amp_mode,
                "amp_used_dim": np.nan,
            }

        rows.append(out)

    dim = pd.DataFrame(rows)

    # Remove old dimensional columns if present, then replace cleanly
    replace_cols = list(dim.columns)
    base = df.drop(columns=[c for c in replace_cols if c in df.columns], errors="ignore")

    outdf = pd.concat([base.reset_index(drop=True), dim.reset_index(drop=True)], axis=1)

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_summary), exist_ok=True)

    outdf.to_csv(args.output_csv, index=False)

    summary = {
        "title": "BuP SLACS dimensional sector added v1",
        "input_csv": args.input_csv,
        "output_csv": args.output_csv,
        "n": int(len(outdf)),
        "settings": vars(args),
        "formula": "alpha_eff_mean = 2*ds_mean/dw_mean + dw_mean - 4",
        "summary": summarize(outdf),
    }

    with open(args.output_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 100)
    print("DIMENSIONAL SECTOR SUMMARY")
    print("=" * 100)
    print(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(args.output_csv)
    print(args.output_summary)


if __name__ == "__main__":
    main()

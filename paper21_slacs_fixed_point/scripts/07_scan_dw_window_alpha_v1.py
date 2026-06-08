#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS — scan diffusion windows for alpha_eff = 1

Goal:
    Scan many MSD fit windows [t_min, t_max] on the same graph and test whether
    alpha_eff = 2 ds / dw + dw - 4 can cross 1.

Key idea:
    Previous estimators bracket alpha_eff:
      - fixed late window gives alpha_eff ~ 1.48
      - adaptive early window gives alpha_eff < 0

    Therefore alpha_eff = 1 may occur at an intermediate diffusion window.

Outputs:
    scan_window_results.csv
    scan_window_summary.csv
    best_alpha1_windows.csv
    summary.json
"""

import argparse
import json
import os
import numpy as np
import pandas as pd

from scipy.linalg import eigh
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
                v = float(row[c])
                if np.isfinite(v) and 0.3 < v < 12:
                    return v, c
            except Exception:
                pass
    return float(default), "default"


def get_q(row):
    for c in ["q_SIE", "q", "axis_ratio"]:
        if c in row and pd.notna(row[c]):
            try:
                v = float(row[c])
                if np.isfinite(v) and 0.1 <= v <= 1.5:
                    return v
            except Exception:
                pass
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


def transition_matrix(W):
    d = np.sum(W, axis=1)
    return W / np.maximum(d[:, None], 1e-15)


def laplacian_norm(W):
    d = np.sum(W, axis=1)
    invsqrt = 1.0 / np.sqrt(np.maximum(d, 1e-15))
    return np.eye(W.shape[0]) - invsqrt[:, None] * W * invsqrt[None, :]


def pairwise_dist2(nodes, Re):
    x = nodes[:, 0] / max(Re, 1e-12)
    y = nodes[:, 1] / max(Re, 1e-12)

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return dx * dx + dy * dy


def fit_loglog(t, y):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(t) & np.isfinite(y) & (t > 0) & (y > 0)

    if np.sum(m) < 4:
        return np.nan, np.nan, int(np.sum(m))

    X = np.log(t[m]).reshape(-1, 1)
    Y = np.log(y[m])

    model = LinearRegression()
    model.fit(X, Y)

    slope = float(model.coef_[0])
    r2 = float(model.score(X, Y))

    return slope, r2, int(np.sum(m))


def compute_ds(Lnorm, t_min=0.05, t_max=20.0, n_t=80, p_low=0.05, p_high=0.85):
    vals = np.linalg.eigvalsh(Lnorm)
    vals = np.maximum(vals, 0.0)

    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)
    P = np.array([np.mean(np.exp(-tt * vals)) for tt in t])

    mask = np.isfinite(P) & (P > p_low) & (P < p_high)

    if np.sum(mask) < 6:
        i0 = int(0.20 * len(t))
        i1 = int(0.75 * len(t))
        mask = np.zeros_like(P, dtype=bool)
        mask[i0:i1] = True

    slope, r2, nfit = fit_loglog(t[mask], P[mask])
    ds = -2.0 * slope if np.isfinite(slope) else np.nan

    return float(ds), float(slope), float(r2), int(nfit)


def msd_discrete(P, D2, t_max):
    Pt = np.eye(P.shape[0])
    times = []
    msd = []

    for step in range(1, t_max + 1):
        Pt = Pt @ P
        val = np.mean(np.sum(Pt * D2, axis=1))

        if np.isfinite(val) and val > 0:
            times.append(float(step))
            msd.append(float(val))

    return np.asarray(times), np.asarray(msd)


def dw_from_window(times, msd, t_low, t_high):
    mask = (times >= t_low) & (times <= t_high)

    slope, r2, nfit = fit_loglog(times[mask], msd[mask])

    if not np.isfinite(slope) or abs(slope) < 1e-12:
        return np.nan, np.nan, r2, nfit

    dw = 2.0 / slope
    return float(dw), float(slope), float(r2), int(nfit)


def alpha_eff(ds, dw):
    if not np.isfinite(ds) or not np.isfinite(dw) or abs(dw) < 1e-12:
        return np.nan
    return float(2.0 * ds / dw + dw - 4.0)


def compute_graph_for_row(row, args):
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
    sigma = sigma / np.maximum(np.sum(sigma), 1e-15)

    xi = args.xi_factor * Re
    W = build_W(nodes, sigma, xi)

    P = transition_matrix(W)
    Lnorm = laplacian_norm(W)
    D2 = pairwise_dist2(nodes, Re)

    ds, ds_slope, ds_r2, ds_nfit = compute_ds(
        Lnorm,
        t_min=args.ds_t_min,
        t_max=args.ds_t_max,
        n_t=args.ds_n_t,
        p_low=args.ds_p_low,
        p_high=args.ds_p_high,
    )

    times, msd = msd_discrete(P, D2, t_max=args.t_max)

    return {
        "Re": Re,
        "q": q,
        "sersic_n": nser,
        "sersic_n_source": n_source,
        "ds": ds,
        "ds_slope": ds_slope,
        "ds_r2": ds_r2,
        "ds_nfit": ds_nfit,
        "times": times,
        "msd": msd,
    }


def scan_windows_for_dataset(df, args):
    windows = []

    t_lows = list(range(args.t_low_min, args.t_low_max + 1))
    t_highs = list(range(args.t_high_min, args.t_high_max + 1))

    valid_windows = []
    for lo in t_lows:
        for hi in t_highs:
            if hi - lo + 1 >= args.min_fit_points and hi > lo:
                valid_windows.append((lo, hi))

    print(f"Number of diffusion windows tested: {len(valid_windows)}")

    all_rows = []

    # Precompute graphs once per galaxy
    graph_cache = []

    for k, (_, row) in enumerate(df.iterrows(), start=1):
        name = row["name"]
        print(f"[{k}/{len(df)}] building graph: {name}")

        try:
            g = compute_graph_for_row(row, args)
            graph_cache.append((row, g))
        except Exception as e:
            print("  FAILED:", e)

    for lo, hi in valid_windows:
        alpha_vals = []
        dw_vals = []
        ds_vals = []
        r2_vals = []
        logM_vals = []
        names = []

        for row, g in graph_cache:
            dw, slope, r2, nfit = dw_from_window(g["times"], g["msd"], lo, hi)
            a = alpha_eff(g["ds"], dw)

            if np.isfinite(a) and np.isfinite(dw):
                alpha_vals.append(a)
                dw_vals.append(dw)
                ds_vals.append(g["ds"])
                r2_vals.append(r2)
                logM_vals.append(float(row["logMs"]))
                names.append(row["name"])

        if len(alpha_vals) < 8:
            continue

        alpha_vals = np.asarray(alpha_vals)
        dw_vals = np.asarray(dw_vals)
        ds_vals = np.asarray(ds_vals)
        r2_vals = np.asarray(r2_vals)
        logM_vals = np.asarray(logM_vals)

        # Core metrics
        mean_alpha = float(np.mean(alpha_vals))
        med_alpha = float(np.median(alpha_vals))
        std_alpha = float(np.std(alpha_vals))
        mean_dw = float(np.mean(dw_vals))
        mean_ds = float(np.mean(ds_vals))
        mean_r2 = float(np.nanmean(r2_vals))

        # Does this diffusion window globally approach alpha=1?
        abs_mean_minus_1 = abs(mean_alpha - 1.0)
        abs_med_minus_1 = abs(med_alpha - 1.0)

        # Focus around observed fixed-point window
        mfix = (logM_vals >= args.fixed_low) & (logM_vals <= args.fixed_high)

        if np.sum(mfix) >= 4:
            alpha_fix = alpha_vals[mfix]
            dw_fix = dw_vals[mfix]
            ds_fix = ds_vals[mfix]
            mean_alpha_fix = float(np.mean(alpha_fix))
            med_alpha_fix = float(np.median(alpha_fix))
            abs_fix_minus_1 = abs(mean_alpha_fix - 1.0)
            n_fix = int(np.sum(mfix))
            mean_dw_fix = float(np.mean(dw_fix))
            mean_ds_fix = float(np.mean(ds_fix))
        else:
            mean_alpha_fix = np.nan
            med_alpha_fix = np.nan
            abs_fix_minus_1 = np.nan
            n_fix = int(np.sum(mfix))
            mean_dw_fix = np.nan
            mean_ds_fix = np.nan

        # Correlation alpha vs logM
        if np.std(alpha_vals) > 1e-12 and np.std(logM_vals) > 1e-12:
            corr_logM = float(np.corrcoef(logM_vals, alpha_vals)[0, 1])
        else:
            corr_logM = np.nan

        all_rows.append({
            "t_low": lo,
            "t_high": hi,
            "n": int(len(alpha_vals)),
            "mean_alpha": mean_alpha,
            "median_alpha": med_alpha,
            "std_alpha": std_alpha,
            "mean_dw": mean_dw,
            "mean_ds": mean_ds,
            "mean_r2": mean_r2,
            "abs_mean_alpha_minus_1": abs_mean_minus_1,
            "abs_median_alpha_minus_1": abs_med_minus_1,
            "n_fixed_window": n_fix,
            "mean_alpha_fixed_window": mean_alpha_fix,
            "median_alpha_fixed_window": med_alpha_fix,
            "mean_dw_fixed_window": mean_dw_fix,
            "mean_ds_fixed_window": mean_ds_fix,
            "abs_fixed_window_alpha_minus_1": abs_fix_minus_1,
            "alpha_corr_logM": corr_logM,
        })

    return pd.DataFrame(all_rows)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--tag", required=True)

    ap.add_argument("--n-rings", type=int, default=12)
    ap.add_argument("--n-theta", type=int, default=16)
    ap.add_argument("--rmin-factor", type=float, default=0.03)
    ap.add_argument("--rmax-factor", type=float, default=5.0)
    ap.add_argument("--xi-factor", type=float, default=0.35)
    ap.add_argument("--default-sersic-n", type=float, default=4.0)

    ap.add_argument("--ds-t-min", type=float, default=0.05)
    ap.add_argument("--ds-t-max", type=float, default=20.0)
    ap.add_argument("--ds-n-t", type=int, default=80)
    ap.add_argument("--ds-p-low", type=float, default=0.05)
    ap.add_argument("--ds-p-high", type=float, default=0.85)

    ap.add_argument("--t-max", type=int, default=60)
    ap.add_argument("--t-low-min", type=int, default=1)
    ap.add_argument("--t-low-max", type=int, default=25)
    ap.add_argument("--t-high-min", type=int, default=4)
    ap.add_argument("--t-high-max", type=int, default=60)
    ap.add_argument("--min-fit-points", type=int, default=5)

    ap.add_argument("--fixed-low", type=float, default=11.545)
    ap.add_argument("--fixed-high", type=float, default=11.645)

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = read_csv_auto(args.input_csv)

    required = ["name", "Re_arcsec", "logMs"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    result = scan_windows_for_dataset(df, args)

    out_csv = os.path.join(args.output_dir, f"{args.tag}_scan_window_results.csv")
    result.to_csv(out_csv, index=False)

    # Best windows globally and inside fixed-point mass window.
    best_global = result.sort_values("abs_mean_alpha_minus_1").head(30)
    best_fixed = result.dropna(subset=["abs_fixed_window_alpha_minus_1"]).sort_values(
        "abs_fixed_window_alpha_minus_1"
    ).head(30)

    best_global_csv = os.path.join(args.output_dir, f"{args.tag}_best_global_alpha1_windows.csv")
    best_fixed_csv = os.path.join(args.output_dir, f"{args.tag}_best_fixedmass_alpha1_windows.csv")

    best_global.to_csv(best_global_csv, index=False)
    best_fixed.to_csv(best_fixed_csv, index=False)

    summary = {
        "title": "BuP SLACS diffusion-window alpha scan v1",
        "tag": args.tag,
        "input_csv": args.input_csv,
        "n": int(len(df)),
        "fixed_mass_window": [args.fixed_low, args.fixed_high],
        "best_global": best_global.iloc[0].to_dict() if len(best_global) else None,
        "best_fixed_mass_window": best_fixed.iloc[0].to_dict() if len(best_fixed) else None,
        "files": {
            "all_windows": os.path.basename(out_csv),
            "best_global": os.path.basename(best_global_csv),
            "best_fixed_mass_window": os.path.basename(best_fixed_csv),
        }
    }

    with open(os.path.join(args.output_dir, f"{args.tag}_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\nBEST GLOBAL WINDOWS")
    print(best_global.head(15).to_string(index=False))

    print("\nBEST FIXED-MASS WINDOWS")
    print(best_fixed.head(15).to_string(index=False))

    print("\nWrote:")
    print(out_csv)
    print(best_global_csv)
    print(best_fixed_csv)


if __name__ == "__main__":
    main()

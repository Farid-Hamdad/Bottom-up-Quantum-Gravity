#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS — alternative d_w estimators v1

Goal:
    Test whether the current high d_w ~ 4.77 is intrinsic to the graph
    or caused by the MSD fitting method.

Inputs:
    CSV with Re_arcsec, logMs, Sérsic n columns.

Outputs:
    - CSV with ds_mean and several d_w estimators
    - alpha_eff for each d_w estimator
    - summary table
    - alpha=1 crossing test for each estimator

Definitions:
    alpha_eff = 2 ds / dw + dw - 4
"""

import argparse
import json
import os
import numpy as np
import pandas as pd

from scipy.linalg import eigh
from scipy.optimize import brentq
from sklearn.linear_model import LinearRegression, HuberRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


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
                if np.isfinite(val) and 0.3 < val < 12:
                    return val, c
            except Exception:
                pass
    return float(default), "default"


def get_q(row):
    for c in ["q_SIE", "q", "axis_ratio"]:
        if c in row and pd.notna(row[c]):
            try:
                val = float(row[c])
                if np.isfinite(val) and 0.1 <= val <= 1.5:
                    return val
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


def pairwise_dist2(nodes, normalize_re=None):
    x = nodes[:, 0]
    y = nodes[:, 1]

    if normalize_re is not None and normalize_re > 0:
        x = x / normalize_re
        y = y / normalize_re

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return dx * dx + dy * dy


def transition_matrix(W):
    d = np.sum(W, axis=1)
    return W / np.maximum(d[:, None], 1e-15)


def laplacian_raw(W):
    d = np.sum(W, axis=1)
    return np.diag(d) - W


def laplacian_norm(W):
    d = np.sum(W, axis=1)
    invsqrt = 1.0 / np.sqrt(np.maximum(d, 1e-15))
    return np.eye(W.shape[0]) - invsqrt[:, None] * W * invsqrt[None, :]


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


def best_adaptive_window(t, msd, min_points=6, frac_low=0.02, frac_high=0.45):
    t = np.asarray(t, dtype=float)
    msd = np.asarray(msd, dtype=float)

    m = np.isfinite(t) & np.isfinite(msd) & (t > 0) & (msd > 0)
    t = t[m]
    msd = msd[m]

    if len(t) < min_points:
        return np.nan, np.nan, 0, "too_few"

    sat = np.nanmax(msd)
    if not np.isfinite(sat) or sat <= 0:
        return np.nan, np.nan, 0, "bad_sat"

    frac = msd / sat

    mask = (frac >= frac_low) & (frac <= frac_high)

    if np.sum(mask) >= min_points:
        slope, r2, nfit = fit_loglog(t[mask], msd[mask])
        return slope, r2, nfit, "fraction_window"

    # fallback: search best contiguous log-log window
    best = None
    n = len(t)

    for w in range(min_points, min(20, n) + 1):
        for i in range(0, n - w + 1):
            tt = t[i:i + w]
            yy = msd[i:i + w]
            slope, r2, nfit = fit_loglog(tt, yy)

            if not np.isfinite(slope) or not np.isfinite(r2):
                continue
            if slope <= 0:
                continue

            score = r2 - 0.05 * abs(slope - 1.0)

            if best is None or score > best[0]:
                best = (score, slope, r2, nfit)

    if best is None:
        slope, r2, nfit = fit_loglog(t, msd)
        return slope, r2, nfit, "fallback_all"

    _, slope, r2, nfit = best
    return slope, r2, nfit, "best_contiguous"


def dw_from_slope(slope):
    if not np.isfinite(slope) or abs(slope) < 1e-12:
        return np.nan
    return 2.0 / slope


def compute_ds_from_Lnorm(Lnorm, t_min=0.05, t_max=20.0, n_t=80, p_low=0.05, p_high=0.85):
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

    return ds, slope, r2, nfit


def msd_discrete(P, D2, t_max=60):
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


def heat_kernel_msd(L, D2, t_grid):
    # L must be symmetric.
    vals, vecs = eigh(L)
    vals = np.maximum(vals, 0.0)

    msd = []

    for t in t_grid:
        K = (vecs * np.exp(-t * vals)[None, :]) @ vecs.T

        # Turn kernel into row-wise transition-like kernel.
        K = np.maximum(K, 0.0)
        rowsum = np.sum(K, axis=1)
        K = K / np.maximum(rowsum[:, None], 1e-15)

        val = np.mean(np.sum(K * D2, axis=1))
        msd.append(float(val))

    return np.asarray(t_grid, dtype=float), np.asarray(msd, dtype=float)


def compute_dw_estimators(W, nodes, Re, args):
    D2 = pairwise_dist2(nodes, normalize_re=Re)

    P = transition_matrix(W)
    Lraw = laplacian_raw(W)
    Lnorm = laplacian_norm(W)

    out = {}

    # 1. Old-ish fixed discrete estimator.
    t_disc, msd_disc = msd_discrete(P, D2, t_max=args.discrete_t_max)

    fixed_mask = (t_disc >= args.fixed_t_min) & (t_disc <= args.fixed_t_max)
    slope, r2, nfit = fit_loglog(t_disc[fixed_mask], msd_disc[fixed_mask])
    out["dw_discrete_fixed"] = dw_from_slope(slope)
    out["dw_discrete_fixed_slope"] = slope
    out["dw_discrete_fixed_r2"] = r2
    out["dw_discrete_fixed_nfit"] = nfit

    # 2. Adaptive discrete estimator.
    slope, r2, nfit, method = best_adaptive_window(
        t_disc,
        msd_disc,
        min_points=args.min_fit_points,
        frac_low=args.msd_frac_low,
        frac_high=args.msd_frac_high,
    )
    out["dw_discrete_adaptive"] = dw_from_slope(slope)
    out["dw_discrete_adaptive_slope"] = slope
    out["dw_discrete_adaptive_r2"] = r2
    out["dw_discrete_adaptive_nfit"] = nfit
    out["dw_discrete_adaptive_method"] = method

    # 3. Heat kernel with normalized Laplacian.
    t_heat = np.logspace(np.log10(args.heat_t_min), np.log10(args.heat_t_max), args.heat_n_t)
    th, mh = heat_kernel_msd(Lnorm, D2, t_heat)
    slope, r2, nfit, method = best_adaptive_window(
        th,
        mh,
        min_points=args.min_fit_points,
        frac_low=args.msd_frac_low,
        frac_high=args.msd_frac_high,
    )
    out["dw_heat_norm_adaptive"] = dw_from_slope(slope)
    out["dw_heat_norm_slope"] = slope
    out["dw_heat_norm_r2"] = r2
    out["dw_heat_norm_nfit"] = nfit
    out["dw_heat_norm_method"] = method

    # 4. Heat kernel with rescaled raw Laplacian.
    # Rescale Lraw by its median positive eigenvalue to avoid time-scale problems.
    vals = np.linalg.eigvalsh(Lraw)
    pos = vals[vals > 1e-12]
    scale = float(np.median(pos)) if len(pos) else 1.0
    Lraw_scaled = Lraw / max(scale, 1e-15)

    th, mh = heat_kernel_msd(Lraw_scaled, D2, t_heat)
    slope, r2, nfit, method = best_adaptive_window(
        th,
        mh,
        min_points=args.min_fit_points,
        frac_low=args.msd_frac_low,
        frac_high=args.msd_frac_high,
    )
    out["dw_heat_raw_adaptive"] = dw_from_slope(slope)
    out["dw_heat_raw_slope"] = slope
    out["dw_heat_raw_r2"] = r2
    out["dw_heat_raw_nfit"] = nfit
    out["dw_heat_raw_method"] = method

    return out


def alpha_from(ds, dw):
    if not np.isfinite(ds) or not np.isfinite(dw) or abs(dw) < 1e-12:
        return np.nan
    return 2.0 * ds / dw + dw - 4.0


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
    sigma = sigma / np.maximum(np.sum(sigma), 1e-15)

    if args.amp_mode == "shape":
        amp = 1.0
    elif args.amp_mode == "sqrtM":
        amp = np.sqrt(max(get_mstar_rel(row), 1e-12))
    elif args.amp_mode == "Mstar":
        amp = max(get_mstar_rel(row), 1e-12)
    else:
        raise ValueError(args.amp_mode)

    sigma_amp = sigma * amp
    xi = args.xi_factor * Re

    W = build_W(nodes, sigma_amp, xi)
    Lnorm = laplacian_norm(W)

    ds, ds_slope, ds_r2, ds_nfit = compute_ds_from_Lnorm(
        Lnorm,
        t_min=args.ds_t_min,
        t_max=args.ds_t_max,
        n_t=args.ds_n_t,
        p_low=args.ds_p_low,
        p_high=args.ds_p_high,
    )

    dws = compute_dw_estimators(W, nodes, Re, args)

    out = {
        "ds_mean": ds,
        "ds_slope": ds_slope,
        "ds_fit_r2": ds_r2,
        "ds_nfit": ds_nfit,
        "sersic_n_used": nser,
        "sersic_n_source": n_source,
        "q_used": q,
        "xi_used": xi,
        "amp_used": amp,
    }

    out.update(dws)

    for key in [
        "dw_discrete_fixed",
        "dw_discrete_adaptive",
        "dw_heat_norm_adaptive",
        "dw_heat_raw_adaptive",
    ]:
        out[f"alpha_eff_{key}"] = alpha_from(ds, out.get(key, np.nan))

    return out


def fit_curve(x, y, model_type):
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    y = np.asarray(y, dtype=float)

    if model_type == "linear":
        model = LinearRegression()
    elif model_type == "quadratic":
        model = make_pipeline(PolynomialFeatures(2, include_bias=False), LinearRegression())
    elif model_type == "cubic":
        model = make_pipeline(PolynomialFeatures(3, include_bias=False), LinearRegression())
    elif model_type == "huber":
        model = HuberRegressor(epsilon=1.35, alpha=1e-6, max_iter=1000)
    else:
        raise ValueError(model_type)

    model.fit(x, y)
    return model


def predict_curve(model, x):
    return np.asarray(model.predict(np.asarray(x).reshape(-1, 1)), dtype=float)


def crossing_alpha_one(df, alpha_col):
    rows = []

    work = df.dropna(subset=["logMs", alpha_col]).copy()

    if len(work) < 8:
        return pd.DataFrame()

    xmin = float(work["logMs"].min())
    xmax = float(work["logMs"].max())

    for mt in ["linear", "quadratic", "cubic", "huber"]:
        try:
            model = fit_curve(work["logMs"].values, work[alpha_col].values, mt)
            grid = np.linspace(xmin, xmax, 1000)
            vals = predict_curve(model, grid) - 1.0

            roots = []
            for i in range(len(grid) - 1):
                a, b = grid[i], grid[i + 1]
                fa, fb = vals[i], vals[i + 1]

                if not np.isfinite(fa) or not np.isfinite(fb):
                    continue

                if fa == 0:
                    roots.append(a)
                elif fa * fb < 0:
                    try:
                        roots.append(brentq(lambda z: float(predict_curve(model, [z])[0] - 1.0), a, b))
                    except Exception:
                        pass

            if roots:
                root = min(roots, key=lambda r: abs(r - 11.58))
                method = "sign_change"
                val = float(predict_curve(model, [root])[0])
            else:
                idx = int(np.argmin(np.abs(vals)))
                root = float(grid[idx])
                method = "closest_no_sign_change"
                val = float(vals[idx] + 1.0)

            rows.append({
                "alpha_col": alpha_col,
                "model_type": mt,
                "root_logM": root,
                "method": method,
                "value_at_root": val,
                "n": int(len(work)),
            })

        except Exception as e:
            rows.append({
                "alpha_col": alpha_col,
                "model_type": mt,
                "root_logM": np.nan,
                "method": f"failed:{e}",
                "value_at_root": np.nan,
                "n": int(len(work)),
            })

    return pd.DataFrame(rows)


def summarize(df, tag):
    cols = [
        "ds_mean",
        "dw_discrete_fixed",
        "dw_discrete_adaptive",
        "dw_heat_norm_adaptive",
        "dw_heat_raw_adaptive",
        "alpha_eff_dw_discrete_fixed",
        "alpha_eff_dw_discrete_adaptive",
        "alpha_eff_dw_heat_norm_adaptive",
        "alpha_eff_dw_heat_raw_adaptive",
    ]

    rows = []

    for c in cols:
        if c not in df.columns:
            continue

        x = pd.to_numeric(df[c], errors="coerce")
        rows.append({
            "tag": tag,
            "quantity": c,
            "n": int(x.notna().sum()),
            "mean": float(x.mean()),
            "median": float(x.median()),
            "std": float(x.std()),
            "min": float(x.min()),
            "max": float(x.max()),
        })

    return pd.DataFrame(rows)


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
    ap.add_argument("--amp-mode", choices=["shape", "sqrtM", "Mstar"], default="shape")

    ap.add_argument("--ds-t-min", type=float, default=0.05)
    ap.add_argument("--ds-t-max", type=float, default=20.0)
    ap.add_argument("--ds-n-t", type=int, default=80)
    ap.add_argument("--ds-p-low", type=float, default=0.05)
    ap.add_argument("--ds-p-high", type=float, default=0.85)

    ap.add_argument("--discrete-t-max", type=int, default=60)
    ap.add_argument("--fixed-t-min", type=int, default=1)
    ap.add_argument("--fixed-t-max", type=int, default=25)

    ap.add_argument("--heat-t-min", type=float, default=0.01)
    ap.add_argument("--heat-t-max", type=float, default=30.0)
    ap.add_argument("--heat-n-t", type=int, default=80)

    ap.add_argument("--min-fit-points", type=int, default=6)
    ap.add_argument("--msd-frac-low", type=float, default=0.02)
    ap.add_argument("--msd-frac-high", type=float, default=0.45)

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = read_csv_auto(args.input_csv)

    required = ["name", "Re_arcsec", "logMs"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    rows = []

    print(f"Running d_w estimator test for {args.tag}")
    print("Input:", args.input_csv)

    for k, (_, row) in enumerate(df.iterrows(), start=1):
        name = row["name"]
        print(f"[{k}/{len(df)}] {name}")

        try:
            out = compute_one(row, args)
        except Exception as e:
            print("  FAILED:", e)
            out = {
                "ds_mean": np.nan,
                "dw_discrete_fixed": np.nan,
                "dw_discrete_adaptive": np.nan,
                "dw_heat_norm_adaptive": np.nan,
                "dw_heat_raw_adaptive": np.nan,
                "alpha_eff_dw_discrete_fixed": np.nan,
                "alpha_eff_dw_discrete_adaptive": np.nan,
                "alpha_eff_dw_heat_norm_adaptive": np.nan,
                "alpha_eff_dw_heat_raw_adaptive": np.nan,
                "sersic_n_used": np.nan,
                "sersic_n_source": "failed",
            }

        rows.append(out)

    est = pd.DataFrame(rows)

    full = pd.concat([df.reset_index(drop=True), est.reset_index(drop=True)], axis=1)

    out_csv = os.path.join(args.output_dir, f"{args.tag}_dw_estimators.csv")
    full.to_csv(out_csv, index=False)

    summary = summarize(full, args.tag)
    summary_csv = os.path.join(args.output_dir, f"{args.tag}_dw_estimators_summary.csv")
    summary.to_csv(summary_csv, index=False)

    alpha_cols = [
        "alpha_eff_dw_discrete_fixed",
        "alpha_eff_dw_discrete_adaptive",
        "alpha_eff_dw_heat_norm_adaptive",
        "alpha_eff_dw_heat_raw_adaptive",
    ]

    cross_tables = []
    for ac in alpha_cols:
        cross = crossing_alpha_one(full, ac)
        if len(cross):
            cross["tag"] = args.tag
            cross_tables.append(cross)

    if cross_tables:
        crossings = pd.concat(cross_tables, ignore_index=True)
    else:
        crossings = pd.DataFrame()

    crossing_csv = os.path.join(args.output_dir, f"{args.tag}_alpha_one_crossings.csv")
    crossings.to_csv(crossing_csv, index=False)

    result = {
        "title": "BuP SLACS alternative d_w estimators v1",
        "tag": args.tag,
        "input_csv": args.input_csv,
        "n": int(len(full)),
        "settings": vars(args),
        "files": {
            "estimators": os.path.basename(out_csv),
            "summary": os.path.basename(summary_csv),
            "alpha_one_crossings": os.path.basename(crossing_csv),
        },
    }

    with open(os.path.join(args.output_dir, f"{args.tag}_summary.json"), "w") as f:
        json.dump(result, f, indent=2)

    print("\nSUMMARY")
    print(summary.to_string(index=False))

    print("\nALPHA=1 CROSSINGS")
    print(crossings.to_string(index=False))

    print("\nWrote:")
    print(out_csv)
    print(summary_csv)
    print(crossing_csv)


if __name__ == "__main__":
    main()

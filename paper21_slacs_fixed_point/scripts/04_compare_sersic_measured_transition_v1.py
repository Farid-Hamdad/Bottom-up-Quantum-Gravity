#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS — compare measured Sérsic n vs forced n=4 transition v1

Goal:
    For each dataset:
      - measured_n_only
      - matched_forced_n4
      - measured_plus_fallback

    Use already computed absolute_amplitude_v2_features.csv,
    compare Phi_BuP vs logM galaxy-by-galaxy,
    find the sliding-window transition where Phi_BuP beats logM,
    and estimate C_obs=0 fixed point.

Outputs:
    compare_sersic_transition_summary.csv
    compare_sersic_transition_windows.csv
    compare_sersic_transition_fixed_points.csv
    summary.json
"""

import argparse
import json
import os
import numpy as np
import pandas as pd

from scipy.optimize import brentq
from sklearn.linear_model import LinearRegression, HuberRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


def load_features(path):
    df = pd.read_csv(path)
    required = [
        "name",
        "logMs",
        "theta_E_obs_arcsec",
        "theta_E_baryon_arcsec",
    ]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column {c} in {path}")

    df = df.dropna(subset=required).copy()
    df = df[(df["theta_E_obs_arcsec"] > 0) & (df["theta_E_baryon_arcsec"] > 0)].copy()
    df["C_obs"] = np.log(df["theta_E_obs_arcsec"] / df["theta_E_baryon_arcsec"])
    return df


def standardize(train, test, features):
    Xtr = [np.ones(len(train))]
    Xte = [np.ones(len(test))]

    for f in features:
        xtr = train[f].values.astype(float)
        xte = test[f].values.astype(float)

        mu = np.nanmean(xtr)
        sd = np.nanstd(xtr)

        if not np.isfinite(sd) or sd < 1e-12:
            sd = 1.0

        xtr = np.where(np.isfinite(xtr), xtr, mu)
        xte = np.where(np.isfinite(xte), xte, mu)

        Xtr.append((xtr - mu) / sd)
        Xte.append((xte - mu) / sd)

    return np.column_stack(Xtr), np.column_stack(Xte)


def fit_predict(train, test, features):
    y = train["C_obs"].values.astype(float)
    Xtr, Xte = standardize(train, test, features)
    beta = np.linalg.lstsq(Xtr, y, rcond=None)[0]
    return Xte @ beta


def loo_pred(df, features):
    pred = np.zeros(len(df))
    for i in range(len(df)):
        train = df.drop(df.index[i])
        test = df.iloc[[i]]
        pred[i] = fit_predict(train, test, features)[0]
    return pred


def score_logtheta(df, C_pred):
    log_obs = np.log(df["theta_E_obs_arcsec"].values)
    log_bar = np.log(df["theta_E_baryon_arcsec"].values)
    log_pred = log_bar + C_pred

    rmse_bar = np.sqrt(np.mean((log_obs - log_bar) ** 2))
    rmse = np.sqrt(np.mean((log_obs - log_pred) ** 2))
    improvement = 100.0 * (rmse_bar - rmse) / rmse_bar

    return rmse, improvement


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


def zero_crossing(df, model_type):
    x = df["logMs"].values.astype(float)
    y = df["C_obs"].values.astype(float)

    model = fit_curve(x, y, model_type)

    xmin, xmax = float(np.min(x)), float(np.max(x))
    grid = np.linspace(xmin, xmax, 1000)
    vals = predict_curve(model, grid)

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
                roots.append(brentq(lambda z: float(predict_curve(model, [z])[0]), a, b))
            except Exception:
                pass

    if roots:
        root = min(roots, key=lambda r: abs(r - 11.58))
        method = "sign_change"
        val = float(predict_curve(model, [root])[0])
        n_roots = len(roots)
    else:
        idx = int(np.argmin(np.abs(vals)))
        root = float(grid[idx])
        method = "closest_no_sign_change"
        val = float(vals[idx])
        n_roots = 0

    return {
        "model_type": model_type,
        "root_logM": float(root),
        "method": method,
        "value_at_root": val,
        "n_roots": int(n_roots),
    }


def bootstrap_zero(df, model_type, n_boot, seed):
    rng = np.random.default_rng(seed)
    roots = []
    n = len(df)

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = df.iloc[idx].copy()
        try:
            z = zero_crossing(sample, model_type)
            roots.append(z["root_logM"])
        except Exception:
            pass

    roots = np.asarray(roots, dtype=float)
    roots = roots[np.isfinite(roots)]

    if len(roots) == 0:
        return {}

    return {
        "n_ok": int(len(roots)),
        "median": float(np.quantile(roots, 0.50)),
        "q05": float(np.quantile(roots, 0.05)),
        "q16": float(np.quantile(roots, 0.16)),
        "q84": float(np.quantile(roots, 0.84)),
        "q95": float(np.quantile(roots, 0.95)),
        "std": float(np.std(roots)),
        "frac_in_11p50_11p65": float(np.mean((roots >= 11.50) & (roots <= 11.65))),
        "frac_in_11p55_11p60": float(np.mean((roots >= 11.55) & (roots <= 11.60))),
    }


def sliding_windows(df, phi_feature, widths, step, n_perm, seed):
    rng = np.random.default_rng(seed)

    df = df.copy()

    df["C_pred_logM"] = loo_pred(df, ["logMs"])
    df["C_pred_Phi"] = loo_pred(df, [phi_feature])

    df["abs_err_logM"] = np.abs(df["C_obs"] - df["C_pred_logM"])
    df["abs_err_Phi"] = np.abs(df["C_obs"] - df["C_pred_Phi"])
    df["delta_abs_logM_minus_Phi"] = df["abs_err_logM"] - df["abs_err_Phi"]
    df["Phi_better_than_logM"] = df["delta_abs_logM_minus_Phi"] > 0

    rows = []

    xmin = float(df["logMs"].min())
    xmax = float(df["logMs"].max())

    all_delta = df["delta_abs_logM_minus_Phi"].values.copy()

    for width in widths:
        centers = np.arange(xmin + width / 2, xmax - width / 2 + 1e-9, step)

        for center in centers:
            lo = center - width / 2
            hi = center + width / 2
            g = df[(df["logMs"] >= lo) & (df["logMs"] < hi)]

            if len(g) < 6:
                continue

            obs_mean = float(g["delta_abs_logM_minus_Phi"].mean())
            obs_frac = float(g["Phi_better_than_logM"].mean())
            obs_median = float(g["delta_abs_logM_minus_Phi"].median())

            perm_means = []
            perm_fracs = []

            idx = g.index.to_numpy()

            for _ in range(n_perm):
                sh = rng.permutation(all_delta)
                win = sh[idx]
                perm_means.append(np.mean(win))
                perm_fracs.append(np.mean(win > 0))

            perm_means = np.asarray(perm_means)
            perm_fracs = np.asarray(perm_fracs)

            p_mean = (1 + np.sum(perm_means >= obs_mean)) / (n_perm + 1)
            p_frac = (1 + np.sum(perm_fracs >= obs_frac)) / (n_perm + 1)

            rows.append({
                "width": float(width),
                "center": float(center),
                "logM_low": float(lo),
                "logM_high": float(hi),
                "n": int(len(g)),
                "frac_Phi_better_logM": obs_frac,
                "mean_delta_logM_minus_Phi": obs_mean,
                "median_delta_logM_minus_Phi": obs_median,
                "p_mean_delta_one_sided": float(p_mean),
                "p_frac_one_sided": float(p_frac),
                "mean_abs_err_logM": float(g["abs_err_logM"].mean()),
                "mean_abs_err_Phi": float(g["abs_err_Phi"].mean()),
                "names": ";".join(g["name"].astype(str).tolist()),
            })

    win = pd.DataFrame(rows)
    win = win.sort_values(["mean_delta_logM_minus_Phi", "frac_Phi_better_logM"], ascending=False)

    return df, win


def analyze_dataset(tag, feature_path, phi_feature, args):
    df = load_features(feature_path)

    if phi_feature not in df.columns:
        raise ValueError(f"Missing phi feature {phi_feature} in {feature_path}")

    df = df.dropna(subset=[phi_feature, "logMs", "C_obs"]).copy()

    # LOO global
    pred_logM = loo_pred(df, ["logMs"])
    pred_phi = loo_pred(df, [phi_feature])

    rmse_logM, imp_logM = score_logtheta(df, pred_logM)
    rmse_phi, imp_phi = score_logtheta(df, pred_phi)

    df_gal, win = sliding_windows(
        df,
        phi_feature=phi_feature,
        widths=args.widths,
        step=args.step,
        n_perm=args.n_perm,
        seed=args.seed,
    )

    best_win = win.iloc[0].to_dict() if len(win) else {}

    # Fixed points
    fixed_rows = []
    for mt in ["linear", "quadratic", "cubic", "huber"]:
        z = zero_crossing(df, mt)
        b = bootstrap_zero(df, mt, args.n_boot, args.seed)

        row = {
            "tag": tag,
            **z,
            **{f"boot_{k}": v for k, v in b.items()},
        }
        fixed_rows.append(row)

    fixed = pd.DataFrame(fixed_rows)

    summary = {
        "tag": tag,
        "n": int(len(df)),
        "phi_feature": phi_feature,
        "logM_LOO_improvement": float(imp_logM),
        "Phi_LOO_improvement": float(imp_phi),
        "Phi_minus_logM_LOO": float(imp_phi - imp_logM),
        "frac_Phi_better_logM_global": float(df_gal["Phi_better_than_logM"].mean()),
        "mean_delta_global": float(df_gal["delta_abs_logM_minus_Phi"].mean()),
        "best_window": best_win,
        "huber_Cobs_zero": float(fixed[fixed["model_type"] == "huber"]["root_logM"].iloc[0]),
        "linear_Cobs_zero": float(fixed[fixed["model_type"] == "linear"]["root_logM"].iloc[0]),
    }

    return summary, df_gal, win, fixed


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--output-dir", required=True)

    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--step", type=float, default=0.025)

    args = ap.parse_args()

    args.widths = [0.10, 0.15, 0.20, 0.25, 0.30]

    os.makedirs(args.output_dir, exist_ok=True)

    configs = {
        "measured_n_only": {
            "path": f"{args.base_dir}/results/paper9_absolute_amplitude_v2_measured_n_only/absolute_amplitude_v2_features.csv",
            "phi_feature": "rawL_zS_amp_sqrtMstar_phi_mid_mean",
        },
        "matched_forced_n4": {
            "path": f"{args.base_dir}/results/paper9_absolute_amplitude_v2_matched_forced_n4/absolute_amplitude_v2_features.csv",
            "phi_feature": "norm_rawS_amp_sqrtMstar_phi_inner_minus_outer",
        },
        "measured_plus_fallback": {
            "path": f"{args.base_dir}/results/paper9_absolute_amplitude_v2_measured_plus_fallback/absolute_amplitude_v2_features.csv",
            "phi_feature": "rawL_zS_amp_sqrtMstar_phi_mid_mean",
        },
    }

    summaries = []
    all_windows = []
    all_fixed = []

    for tag, cfg in configs.items():
        print("\n" + "=" * 100)
        print("Analyzing", tag)
        print("=" * 100)

        summary, pergal, windows, fixed = analyze_dataset(
            tag,
            cfg["path"],
            cfg["phi_feature"],
            args,
        )

        summaries.append(summary)

        pergal.to_csv(os.path.join(args.output_dir, f"{tag}_per_galaxy_phi_vs_logM.csv"), index=False)

        windows["tag"] = tag
        fixed["tag"] = tag

        all_windows.append(windows)
        all_fixed.append(fixed)

        print(json.dumps(summary, indent=2))

    sumdf = pd.DataFrame(summaries)

    windf = pd.concat(all_windows, ignore_index=True)
    fixdf = pd.concat(all_fixed, ignore_index=True)

    sumdf.to_csv(os.path.join(args.output_dir, "compare_sersic_transition_summary.csv"), index=False)
    windf.to_csv(os.path.join(args.output_dir, "compare_sersic_transition_windows.csv"), index=False)
    fixdf.to_csv(os.path.join(args.output_dir, "compare_sersic_transition_fixed_points.csv"), index=False)

    # top windows by tag
    top = (
        windf.sort_values("mean_delta_logM_minus_Phi", ascending=False)
             .groupby("tag")
             .head(10)
             .sort_values(["tag", "mean_delta_logM_minus_Phi"], ascending=[True, False])
    )
    top.to_csv(os.path.join(args.output_dir, "compare_sersic_transition_top_windows.csv"), index=False)

    result = {
        "title": "BuP SLACS measured Sérsic n transition comparison v1",
        "summaries": summaries,
        "files": {
            "summary": "compare_sersic_transition_summary.csv",
            "windows": "compare_sersic_transition_windows.csv",
            "top_windows": "compare_sersic_transition_top_windows.csv",
            "fixed_points": "compare_sersic_transition_fixed_points.csv",
        },
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(result, f, indent=2)

    print("\n" + "=" * 100)
    print("FINAL COMPARISON")
    print("=" * 100)
    print(sumdf.to_string(index=False))

    print("\nTOP WINDOWS")
    print(top[[
        "tag", "width", "center", "logM_low", "logM_high", "n",
        "frac_Phi_better_logM", "mean_delta_logM_minus_Phi",
        "p_mean_delta_one_sided", "p_frac_one_sided"
    ]].head(40).to_string(index=False))

    print("\nFIXED POINTS")
    print(fixdf.to_string(index=False))


if __name__ == "__main__":
    main()

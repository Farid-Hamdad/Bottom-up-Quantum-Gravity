#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP vs ΛCDM/NFW — fresh-fit comparison v2

This script recomputes BOTH models from the same SPARC rotmod files.
It is designed to remove the asymmetry of comparing precomputed BuP chi2
against freshly fitted NFW chi2.

Models
------
NFW:
    V_model^2 = V_bar^2 + V_NFW^2

Fresh BuP phenomenological bridge:
    V_model^2 = V_bar^2 + A * S(r; r_t, w)

where S is a smooth transition. In predictive mode, the transition scale
is fixed from Paper 14:
    r_t = lambda_corr = f_lambda_pred * R_d

The only BuP fitted parameters in predictive mode are:
    A       transition amplitude in velocity^2
    w       transition width

In semi_flexible mode, r_t is also fitted.

Scientific status
-----------------
This is a symmetric fresh-fit comparison, but still a phenomenological
BuP rotation implementation, not the final microscopic BuP solver.

Recommended first run
---------------------
cd ~/bottomup
python3 papers/paper14_modular_spectral_action/scripts/bup_vs_nfw_publication_grade_v2.py \
  --rotmod-dir "/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2" \
  --galaxy-list-csv sparc_fit_quality_category_175.csv \
  --lambda-pred-csv papers/paper14_modular_spectral_action/results/predictive_lambda_LOO_v1/loo_predicted_lambda_results.csv \
  --lambda-model rf_classifier \
  --max-galaxies 32 \
  --select-mode best_bup \
  --bup-mode predictive \
  --output-dir papers/paper14_modular_spectral_action/results/bup_vs_nfw_publication_grade_v2_excellent32
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

G = 4.30091e-6  # kpc (km/s)^2 / Msun


def parse_args():
    p = argparse.ArgumentParser(description="Fresh-fit BuP vs NFW comparison on SPARC rotmod files.")
    p.add_argument("--rotmod-dir", required=True)
    p.add_argument("--galaxy-list-csv", required=True)
    p.add_argument("--lambda-pred-csv", default=None)
    p.add_argument("--lambda-model", default="rf_classifier")
    p.add_argument("--galaxies", nargs="*", default=None)
    p.add_argument("--max-galaxies", type=int, default=32)
    p.add_argument("--select-mode", choices=["all", "best_bup", "quality_balanced", "representative"], default="best_bup")
    p.add_argument("--bup-mode", choices=["predictive", "fixed", "semi_flexible"], default="predictive")
    p.add_argument("--fixed-f-lambda", type=float, default=1.0)
    p.add_argument("--ups-disk", type=float, default=0.5)
    p.add_argument("--ups-bul", type=float, default=0.7)
    p.add_argument("--fixed-c", type=float, default=10.0)
    p.add_argument("--h0", type=float, default=70.0)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--plot-limit", type=int, default=16)
    p.add_argument("--min-points", type=int, default=4)
    return p.parse_args()


def find_rotmod(galaxy: str, rotmod_dir: str):
    d = Path(rotmod_dir)
    for name in [f"{galaxy}_rotmod.dat", f"{galaxy}.dat", f"{galaxy}_rotmod.txt"]:
        p = d / name
        if p.exists():
            return p
    g = galaxy.lower()
    for p in d.glob("*"):
        if p.is_file() and g in p.name.lower() and "rotmod" in p.name.lower():
            return p
    return None


def load_rotmod(path: Path):
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith(";"):
                continue
            nums = []
            for z in s.replace(",", " ").split():
                try:
                    nums.append(float(z))
                except Exception:
                    pass
            if len(nums) >= 6:
                rows.append(nums)

    if not rows:
        raise ValueError(f"No numerical rows found in {path}")

    maxlen = max(len(r) for r in rows)
    arr = np.full((len(rows), maxlen), np.nan)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r

    data = {
        "r": arr[:, 0],
        "vobs": arr[:, 1],
        "verr": arr[:, 2],
        "vgas": arr[:, 3],
        "vdisk": arr[:, 4],
        "vbul": arr[:, 5] if maxlen > 5 else np.zeros(len(arr)),
        "sbdisk": arr[:, 6] if maxlen > 6 else np.full(len(arr), np.nan),
        "sbbul": arr[:, 7] if maxlen > 7 else np.full(len(arr), np.nan),
    }

    mask = (
        np.isfinite(data["r"]) & np.isfinite(data["vobs"]) & np.isfinite(data["verr"]) &
        np.isfinite(data["vgas"]) & np.isfinite(data["vdisk"]) & np.isfinite(data["vbul"]) &
        (data["r"] > 0) & (data["verr"] > 0)
    )
    for k in data:
        data[k] = np.asarray(data[k][mask], dtype=float)
    order = np.argsort(data["r"])
    for k in data:
        data[k] = data[k][order]
    return data


def estimate_rd(data):
    r = data["r"]
    sb = data.get("sbdisk", np.full_like(r, np.nan))
    mask = np.isfinite(sb) & (sb > 0) & np.isfinite(r) & (r > 0)
    if mask.sum() >= 4:
        x = r[mask]
        y = np.log(sb[mask])
        try:
            slope, _ = np.polyfit(x, y, 1)
            if slope < 0:
                rd = -1.0 / slope
                if np.isfinite(rd) and 0.02 < rd < 100:
                    return float(rd), "SBdisk_expfit"
        except Exception:
            pass
    return max(float(np.nanmax(r)) / 3.0, 0.1), "rmax_over_3"


def vbar_squared(data, ups_disk=0.5, ups_bul=0.7):
    vgas = data["vgas"]
    vb2 = vgas * np.abs(vgas) + ups_disk * data["vdisk"]**2 + ups_bul * data["vbul"]**2
    return np.maximum(vb2, 0.0)


def rho_crit_msun_kpc3(h0):
    H = h0 / 1000.0
    return 3.0 * H**2 / (8.0 * np.pi * G)


def r200_from_m200(m200, h0):
    rhoc = rho_crit_msun_kpc3(h0)
    return (3.0 * m200 / (4.0 * np.pi * 200.0 * rhoc)) ** (1.0 / 3.0)


def vnfw_squared(r, log10_m200, log10_c, h0):
    m200 = 10.0 ** log10_m200
    c = 10.0 ** log10_c
    r200 = r200_from_m200(m200, h0)
    x = np.maximum(r / r200, 1e-12)
    fc = np.log1p(c) - c / (1.0 + c)
    y = c * x
    fy = np.log1p(y) - y / (1.0 + y)
    v2002 = G * m200 / r200
    return np.maximum(v2002 * fy / (x * fc), 0.0)


def nfw_velocity(data, params, h0, ups_disk, ups_bul, fixed_c=None):
    vb2 = vbar_squared(data, ups_disk, ups_bul)
    if fixed_c is None:
        logm, logc = params
    else:
        logm = params[0]
        logc = np.log10(fixed_c)
    vh2 = vnfw_squared(data["r"], logm, logc, h0)
    return np.sqrt(np.maximum(vb2 + vh2, 0.0))


def fit_nfw(data, h0, ups_disk, ups_bul, fixed_c=None):
    vobs, verr = data["vobs"], data["verr"]
    if fixed_c is None:
        x0 = np.array([11.0, np.log10(10.0)])
        bounds = ([8.0, np.log10(1.0)], [13.8, np.log10(40.0)])
        k = 2
    else:
        x0 = np.array([11.0])
        bounds = ([8.0], [13.8])
        k = 1

    def residuals(x):
        vm = nfw_velocity(data, x, h0, ups_disk, ups_bul, fixed_c=fixed_c)
        return (vobs - vm) / verr

    res = least_squares(residuals, x0=x0, bounds=bounds, max_nfev=6000)
    vm = nfw_velocity(data, res.x, h0, ups_disk, ups_bul, fixed_c=fixed_c)
    chi2 = float(np.sum(((vobs - vm) / verr) ** 2))
    dof = max(len(vobs) - k, 1)
    if fixed_c is None:
        logm = float(res.x[0])
        c = float(10.0 ** res.x[1])
    else:
        logm = float(res.x[0])
        c = float(fixed_c)
    return {
        "k": k, "chi2": chi2, "chi2_red": chi2 / dof, "dof": dof,
        "log10_M200": logm, "M200": float(10.0 ** logm), "c": c,
        "vmodel": vm,
        "rmse": float(np.sqrt(np.mean((vobs - vm) ** 2))),
        "mae": float(np.mean(np.abs(vobs - vm))),
    }


def bup_shape(r, rt, width):
    width = max(float(width), 1e-3)
    rt = max(float(rt), 1e-3)
    sigmoid = 1.0 / (1.0 + np.exp(-(r - rt) / width))
    inner = 1.0 - np.exp(-r / rt)
    return sigmoid * inner


def bup_velocity(data, params, ups_disk, ups_bul, mode, rt_fixed=None):
    vb2 = vbar_squared(data, ups_disk, ups_bul)
    r = data["r"]
    if mode in ["predictive", "fixed"]:
        logA, logw = params
        rt = float(rt_fixed)
    elif mode == "semi_flexible":
        logA, logrt, logw = params
        rt = float(np.exp(logrt))
    else:
        raise ValueError(mode)
    A = float(np.exp(logA))
    width = float(np.exp(logw))
    return np.sqrt(np.maximum(vb2 + A * bup_shape(r, rt, width), 0.0))


def fit_bup(data, ups_disk, ups_bul, mode, rt_fixed):
    vobs, verr = data["vobs"], data["verr"]
    vb = np.sqrt(vbar_squared(data, ups_disk, ups_bul))
    excess = np.maximum(np.nanpercentile(vobs**2 - vb**2, 70), 10.0)

    if mode in ["predictive", "fixed"]:
        rt = max(float(rt_fixed), 0.05)
        x0 = np.array([np.log(excess), np.log(max(rt / 2.0, 0.05))])
        bounds = ([np.log(1e-3), np.log(0.02)], [np.log(1e6), np.log(100.0)])
        k = 2
    else:
        r = data["r"]
        rt0 = max(float(rt_fixed), float(np.nanmedian(r)), 0.05)
        x0 = np.array([np.log(excess), np.log(rt0), np.log(max(rt0 / 2.0, 0.05))])
        bounds = ([np.log(1e-3), np.log(0.02), np.log(0.02)],
                  [np.log(1e6), np.log(200.0), np.log(100.0)])
        k = 3

    def residuals(x):
        vm = bup_velocity(data, x, ups_disk, ups_bul, mode, rt_fixed=rt_fixed)
        return (vobs - vm) / verr

    res = least_squares(residuals, x0=x0, bounds=bounds, max_nfev=6000)
    vm = bup_velocity(data, res.x, ups_disk, ups_bul, mode, rt_fixed=rt_fixed)
    chi2 = float(np.sum(((vobs - vm) / verr) ** 2))
    dof = max(len(vobs) - k, 1)
    if mode in ["predictive", "fixed"]:
        logA, logw = res.x
        rt = float(rt_fixed)
    else:
        logA, logrt, logw = res.x
        rt = float(np.exp(logrt))
    return {
        "k": k, "chi2": chi2, "chi2_red": chi2 / dof, "dof": dof,
        "A_v2": float(np.exp(logA)), "r_trans_BuP": rt, "width_BuP": float(np.exp(logw)),
        "vmodel": vm,
        "rmse": float(np.sqrt(np.mean((vobs - vm) ** 2))),
        "mae": float(np.mean(np.abs(vobs - vm))),
    }


def aic(chi2, k):
    return chi2 + 2.0 * k


def bic(chi2, k, n):
    return chi2 + k * np.log(max(n, 1))


def select_galaxies(df, mode, max_galaxies, explicit=None):
    if explicit:
        return [str(g).strip() for g in explicit]
    d = df.copy()
    if mode == "all":
        return d["galaxy"].astype(str).head(max_galaxies).tolist()
    if mode == "best_bup" and "chi2_red" in d.columns:
        return d.sort_values("chi2_red")["galaxy"].astype(str).head(max_galaxies).tolist()
    if mode == "quality_balanced" and "fit_quality_category" in d.columns:
        names = []
        cats = ["excellent", "good", "medium", "poor"]
        per = max(1, max_galaxies // len(cats))
        for cat in cats:
            sub = d[d["fit_quality_category"] == cat]
            if "chi2_red" in sub.columns:
                sub = sub.sort_values("chi2_red")
            names += sub["galaxy"].astype(str).head(per).tolist()
        if len(names) < max_galaxies:
            extra = d[~d["galaxy"].isin(names)]
            if "chi2_red" in extra.columns:
                extra = extra.sort_values("chi2_red")
            names += extra["galaxy"].astype(str).head(max_galaxies - len(names)).tolist()
        return names[:max_galaxies]
    if "chi2_red" in d.columns:
        d = d.sort_values("chi2_red").reset_index(drop=True)
        idx = np.linspace(0, len(d) - 1, min(max_galaxies, len(d))).round().astype(int)
        return d.iloc[idx]["galaxy"].astype(str).tolist()
    return d["galaxy"].astype(str).head(max_galaxies).tolist()


def load_lambda_predictions(path, model_name):
    if path is None:
        return {}
    df = pd.read_csv(path)
    if "model" in df.columns:
        df = df[df["model"] == model_name].copy()
    if "galaxy" not in df.columns or "pred_f_lambda" not in df.columns:
        raise ValueError("lambda-pred-csv must contain galaxy and pred_f_lambda columns")
    return dict(zip(df["galaxy"].astype(str).str.strip(), df["pred_f_lambda"].astype(float)))


def main():
    args = parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    plots_dir = out / "plots"
    plots_dir.mkdir(exist_ok=True)

    galdf = pd.read_csv(args.galaxy_list_csv)
    galdf["galaxy"] = galdf["galaxy"].astype(str).str.strip()
    galaxies = select_galaxies(galdf, args.select_mode, args.max_galaxies, args.galaxies)
    pred_map = load_lambda_predictions(args.lambda_pred_csv, args.lambda_model)

    rows = []
    plot_count = 0
    for gal in galaxies:
        rot = find_rotmod(gal, args.rotmod_dir)
        if rot is None:
            print(f"[skip] {gal}: rotmod not found")
            continue
        try:
            data = load_rotmod(rot)
            if len(data["r"]) < args.min_points:
                print(f"[skip] {gal}: too few points")
                continue
            rd, rd_source = estimate_rd(data)
            f_lam = float(pred_map.get(gal, args.fixed_f_lambda)) if args.bup_mode == "predictive" else float(args.fixed_f_lambda)
            lambda_corr = max(f_lam * rd, 0.05)

            bup = fit_bup(data, args.ups_disk, args.ups_bul, args.bup_mode, rt_fixed=lambda_corr)
            nfw2 = fit_nfw(data, args.h0, args.ups_disk, args.ups_bul, fixed_c=None)
            nfwc = fit_nfw(data, args.h0, args.ups_disk, args.ups_bul, fixed_c=args.fixed_c)
            n = len(data["r"])
            meta_df = galdf[galdf["galaxy"] == gal]
            meta = meta_df.iloc[0].to_dict() if len(meta_df) else {}

            row = {
                "galaxy": gal,
                "fit_quality_category": meta.get("fit_quality_category", np.nan),
                "legacy_chi2_red_BuP": meta.get("chi2_red", np.nan),
                "n_points": n,
                "rotmod": str(rot),
                "R_d_est": rd,
                "R_d_source": rd_source,
                "f_lambda_used": f_lam,
                "lambda_corr": lambda_corr,
                "bup_mode": args.bup_mode,
                "chi2_red_BuP_fresh": bup["chi2_red"],
                "chi2_BuP_fresh": bup["chi2"],
                "AIC_BuP_fresh": aic(bup["chi2"], bup["k"]),
                "BIC_BuP_fresh": bic(bup["chi2"], bup["k"], n),
                "k_BuP": bup["k"],
                "A_v2_BuP": bup["A_v2"],
                "r_trans_BuP": bup["r_trans_BuP"],
                "width_BuP": bup["width_BuP"],
                "rmse_BuP": bup["rmse"],
                "mae_BuP": bup["mae"],
                "chi2_red_NFW_2param": nfw2["chi2_red"],
                "chi2_NFW_2param": nfw2["chi2"],
                "AIC_NFW_2param": aic(nfw2["chi2"], nfw2["k"]),
                "BIC_NFW_2param": bic(nfw2["chi2"], nfw2["k"], n),
                "log10_M200_NFW_2param": nfw2["log10_M200"],
                "c_NFW_2param": nfw2["c"],
                "rmse_NFW_2param": nfw2["rmse"],
                "mae_NFW_2param": nfw2["mae"],
                "chi2_red_NFW_cfixed": nfwc["chi2_red"],
                "chi2_NFW_cfixed": nfwc["chi2"],
                "AIC_NFW_cfixed": aic(nfwc["chi2"], nfwc["k"]),
                "BIC_NFW_cfixed": bic(nfwc["chi2"], nfwc["k"], n),
                "log10_M200_NFW_cfixed": nfwc["log10_M200"],
                "c_NFW_cfixed": nfwc["c"],
                "rmse_NFW_cfixed": nfwc["rmse"],
                "mae_NFW_cfixed": nfwc["mae"],
            }
            row["winner_chi2red_vs_NFW2"] = "BuP_fresh" if row["chi2_red_BuP_fresh"] < row["chi2_red_NFW_2param"] else "NFW_2param"
            row["winner_chi2red_vs_NFWcfixed"] = "BuP_fresh" if row["chi2_red_BuP_fresh"] < row["chi2_red_NFW_cfixed"] else "NFW_cfixed"
            row["delta_chi2red_BuP_minus_NFW2"] = row["chi2_red_BuP_fresh"] - row["chi2_red_NFW_2param"]
            row["delta_chi2red_BuP_minus_NFWcfixed"] = row["chi2_red_BuP_fresh"] - row["chi2_red_NFW_cfixed"]
            row["delta_AIC_BuP_minus_NFW2"] = row["AIC_BuP_fresh"] - row["AIC_NFW_2param"]
            row["delta_AIC_BuP_minus_NFWcfixed"] = row["AIC_BuP_fresh"] - row["AIC_NFW_cfixed"]
            rows.append(row)

            if plot_count < args.plot_limit:
                r = data["r"]
                vb = np.sqrt(vbar_squared(data, args.ups_disk, args.ups_bul))
                plt.figure(figsize=(7, 5))
                plt.errorbar(r, data["vobs"], yerr=data["verr"], fmt="o", label="SPARC Vobs")
                plt.plot(r, vb, "--", label="Baryons")
                plt.plot(r, bup["vmodel"], label=f"BuP fresh chi2red={bup['chi2_red']:.2f}")
                plt.plot(r, nfw2["vmodel"], label=f"NFW 2p chi2red={nfw2['chi2_red']:.2f}")
                plt.plot(r, nfwc["vmodel"], label=f"NFW c={args.fixed_c:g} chi2red={nfwc['chi2_red']:.2f}")
                plt.title(f"{gal} | f_lambda={f_lam:.2g}, Rd={rd:.2g} kpc")
                plt.xlabel("r [kpc]")
                plt.ylabel("V [km/s]")
                plt.legend(fontsize=8)
                plt.tight_layout()
                plt.savefig(plots_dir / f"{gal}_fresh_bup_vs_nfw.png", dpi=180)
                plt.close()
                plot_count += 1
            print(f"[ok] {gal}: BuPfresh={bup['chi2_red']:.3g} NFW2={nfw2['chi2_red']:.3g} NFWc={nfwc['chi2_red']:.3g} f={f_lam:.2g} Rd={rd:.2g}")
        except Exception as exc:
            print(f"[failed] {gal}: {exc}")

    res = pd.DataFrame(rows)
    if res.empty:
        raise RuntimeError("No successful fits")
    res.to_csv(out / "bup_vs_nfw_publication_grade_results.csv", index=False)

    summary = {
        "experiment": "BuP vs LambdaCDM/NFW fresh-fit prototype v2",
        "n_success": int(len(res)),
        "selection_mode": args.select_mode,
        "max_galaxies": args.max_galaxies,
        "bup_mode": args.bup_mode,
        "lambda_model": args.lambda_model,
        "ups_disk": args.ups_disk,
        "ups_bul": args.ups_bul,
        "fixed_c": args.fixed_c,
        "median_chi2red_BuP_fresh": float(res["chi2_red_BuP_fresh"].median()),
        "median_chi2red_NFW_2param": float(res["chi2_red_NFW_2param"].median()),
        "median_chi2red_NFW_cfixed": float(res["chi2_red_NFW_cfixed"].median()),
        "mean_chi2red_BuP_fresh": float(res["chi2_red_BuP_fresh"].mean()),
        "mean_chi2red_NFW_2param": float(res["chi2_red_NFW_2param"].mean()),
        "mean_chi2red_NFW_cfixed": float(res["chi2_red_NFW_cfixed"].mean()),
        "BuP_wins_vs_NFW2_fraction": float((res["winner_chi2red_vs_NFW2"] == "BuP_fresh").mean()),
        "BuP_wins_vs_NFWcfixed_fraction": float((res["winner_chi2red_vs_NFWcfixed"] == "BuP_fresh").mean()),
        "median_delta_chi2red_BuP_minus_NFW2": float(res["delta_chi2red_BuP_minus_NFW2"].median()),
        "median_delta_chi2red_BuP_minus_NFWcfixed": float(res["delta_chi2red_BuP_minus_NFWcfixed"].median()),
        "median_delta_AIC_BuP_minus_NFW2": float(res["delta_AIC_BuP_minus_NFW2"].median()),
        "median_delta_AIC_BuP_minus_NFWcfixed": float(res["delta_AIC_BuP_minus_NFWcfixed"].median()),
    }
    with open(out / "bup_vs_nfw_publication_grade_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    by_quality = res.groupby("fit_quality_category").agg(
        n=("galaxy", "count"),
        median_BuPfresh=("chi2_red_BuP_fresh", "median"),
        median_NFW2=("chi2_red_NFW_2param", "median"),
        median_NFWcfixed=("chi2_red_NFW_cfixed", "median"),
        BuP_win_NFW2=("winner_chi2red_vs_NFW2", lambda x: (x == "BuP_fresh").mean()),
        BuP_win_NFWcfixed=("winner_chi2red_vs_NFWcfixed", lambda x: (x == "BuP_fresh").mean()),
    ).reset_index()
    by_quality.to_csv(out / "bup_vs_nfw_publication_grade_by_quality.csv", index=False)

    print("\n" + "=" * 110)
    print("BuP vs LambdaCDM/NFW — fresh-fit prototype v2")
    print("=" * 110)
    print(json.dumps(summary, indent=2))
    print("\nTop table:")
    cols = ["galaxy", "fit_quality_category", "legacy_chi2_red_BuP", "chi2_red_BuP_fresh", "chi2_red_NFW_2param", "chi2_red_NFW_cfixed", "winner_chi2red_vs_NFW2", "winner_chi2red_vs_NFWcfixed", "f_lambda_used", "R_d_est", "lambda_corr", "r_trans_BuP", "width_BuP"]
    print(res[cols].sort_values("chi2_red_BuP_fresh").head(40).to_string(index=False))
    print("\nBy quality:")
    print(by_quality.to_string(index=False))
    print("\nFiles written:")
    print(out / "bup_vs_nfw_publication_grade_results.csv")
    print(out / "bup_vs_nfw_publication_grade_summary.json")
    print(out / "bup_vs_nfw_publication_grade_by_quality.csv")
    print(plots_dir)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()

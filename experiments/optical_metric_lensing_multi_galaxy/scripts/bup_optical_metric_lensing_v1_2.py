#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

G_KPC = 4.30091e-6
C_KMS = 299792.458
RAD_TO_ARCSEC = 206264.80624709636
ARCSEC_TO_RAD = 1.0 / RAD_TO_ARCSEC


def load_report(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_points(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"name", "r_kpc", "rd_kpc", "vbar_sq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {sorted(missing)}")
    return df.copy()


def get_global_bump_params(report: dict):
    if "bump_fit" in report:
        bf = report["bump_fit"]
        if all(k in bf for k in ("delta_d", "x0", "sigma")):
            return float(bf["delta_d"]), float(bf["x0"]), float(bf["sigma"])
    metrics = report.get("metrics", {})
    bump = metrics.get("bump_rnorm", {})
    glob = bump.get("global", {})
    if all(k in glob for k in ("delta_d", "x0", "sigma")):
        return float(glob["delta_d"]), float(glob["x0"]), float(glob["sigma"])
    return None


def get_x_ref(report: dict, default: float = 1.0) -> float:
    if "config" in report and "x_ref" in report["config"]:
        return float(report["config"]["x_ref"])
    if "x_ref" in report:
        return float(report["x_ref"])
    return default


def d_profile_x(x: np.ndarray, delta_d: float, x0: float, sigma: float) -> np.ndarray:
    x = np.maximum(np.asarray(x, dtype=float), 1e-10)
    return 3.0 - delta_d * np.exp(-(np.log(x / x0) ** 2) / (2.0 * sigma ** 2))


def geff_factor_x(x: np.ndarray, delta_d: float, x0: float, sigma: float, x_ref: float) -> np.ndarray:
    d = d_profile_x(x, delta_d, x0, sigma)
    x = np.maximum(np.asarray(x, dtype=float), 1e-10)
    return np.power(x / x_ref, 3.0 - d)


def ensure_monotonic_radius(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values("r_kpc").reset_index(drop=True)


def E_z(z: np.ndarray, omega_m: float, omega_l: float) -> np.ndarray:
    return np.sqrt(omega_m * (1.0 + z) ** 3 + omega_l)


def comoving_distance_mpc(z: float, H0: float, omega_m: float) -> float:
    if z <= 0:
        return 0.0
    omega_l = 1.0 - omega_m
    zz = np.linspace(0.0, z, 4000)
    integrand = 1.0 / E_z(zz, omega_m, omega_l)
    c_over_H0 = C_KMS / H0
    return float(c_over_H0 * np.trapezoid(integrand, zz))


def angular_diameter_distance_mpc(z: float, H0: float, omega_m: float) -> float:
    return comoving_distance_mpc(z, H0, omega_m) / (1.0 + z)


def angular_diameter_distance_between_mpc(z1: float, z2: float, H0: float, omega_m: float) -> float:
    if z2 <= z1:
        return 0.0
    dc1 = comoving_distance_mpc(z1, H0, omega_m)
    dc2 = comoving_distance_mpc(z2, H0, omega_m)
    return (dc2 - dc1) / (1.0 + z2)


def theta_arcsec_from_R_kpc(R_kpc: np.ndarray, D_l_mpc: float) -> np.ndarray:
    D_l_kpc = D_l_mpc * 1e3
    return np.asarray(R_kpc, dtype=float) / D_l_kpc * RAD_TO_ARCSEC


def find_crossing(x: np.ndarray, y: np.ndarray, target) -> Optional[float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(target, dtype=float)
    if t.ndim == 0:
        t = np.full_like(y, float(t))
    diff = y - t
    for i in range(len(x) - 1):
        if np.isnan(diff[i]) or np.isnan(diff[i + 1]):
            continue
        if diff[i] == 0:
            return float(x[i])
        if diff[i] * diff[i + 1] < 0:
            w = diff[i] / (diff[i] - diff[i + 1])
            return float(x[i] + w * (x[i + 1] - x[i]))
    return np.nan


def build_baryonic_speed_profile(g: pd.DataFrame):
    R_kpc = g["r_kpc"].to_numpy(dtype=float)
    Rd_kpc = float(np.nanmedian(g["rd_kpc"].to_numpy(dtype=float)))
    x = g["x"].to_numpy(dtype=float) if "x" in g.columns else (R_kpc / max(Rd_kpc, 1e-12))
    vbar_sq = np.maximum(g["vbar_sq"].to_numpy(dtype=float), 0.0)
    v_bary = np.sqrt(vbar_sq)
    return R_kpc, Rd_kpc, x, v_bary


def interpolate_speed(R_eval: np.ndarray, R_grid: np.ndarray, v_grid: np.ndarray) -> np.ndarray:
    return np.interp(R_eval, R_grid, v_grid, left=v_grid[0], right=v_grid[-1])


def g_bary_from_v(r_kpc: np.ndarray, R_grid: np.ndarray, v_bary_grid: np.ndarray) -> np.ndarray:
    r = np.maximum(np.asarray(r_kpc, dtype=float), 1e-10)
    v = interpolate_speed(r, R_grid, v_bary_grid)
    return np.square(v) / r


def alpha_from_geff(R_grid_kpc: np.ndarray, g_eff_func, zmax_kpc: float, nz: int = 3001) -> np.ndarray:
    z = np.linspace(-zmax_kpc, zmax_kpc, nz)
    out = np.zeros_like(R_grid_kpc, dtype=float)
    for i, R in enumerate(R_grid_kpc):
        r = np.sqrt(R**2 + z**2)
        g_eff = g_eff_func(r)
        g_perp = g_eff * (R / np.maximum(r, 1e-10))
        alpha_rad = (2.0 / (C_KMS ** 2)) * np.trapezoid(g_perp, z)
        out[i] = alpha_rad * RAD_TO_ARCSEC
    return out


def local_weight_x(x: np.ndarray, x0: float, width: float) -> np.ndarray:
    x = np.maximum(np.asarray(x, dtype=float), 1e-10)
    return np.exp(-(np.log(x / x0) ** 2) / (2.0 * width ** 2))


def robust_ratio_stats(alpha_bup: np.ndarray, alpha_bary: np.ndarray, theta: np.ndarray, frac_threshold: float = 0.1):
    a_bary_max = float(np.nanmax(alpha_bary))
    thresh = frac_threshold * a_bary_max
    mask = np.asarray(alpha_bary) >= thresh
    ratio = np.full_like(alpha_bary, np.nan, dtype=float)
    ratio[mask] = np.asarray(alpha_bup)[mask] / np.maximum(np.asarray(alpha_bary)[mask], 1e-12)

    if np.any(mask):
        mean_ratio = float(np.nanmean(ratio[mask]))
        median_ratio = float(np.nanmedian(ratio[mask]))
        max_ratio = float(np.nanmax(ratio[mask]))
        theta_min = float(np.nanmin(theta[mask]))
        theta_max = float(np.nanmax(theta[mask]))
    else:
        mean_ratio = np.nan
        median_ratio = np.nan
        max_ratio = np.nan
        theta_min = np.nan
        theta_max = np.nan

    return {
        "mask": mask,
        "ratio": ratio,
        "threshold": thresh,
        "mean_ratio": mean_ratio,
        "median_ratio": median_ratio,
        "max_ratio": max_ratio,
        "theta_window_min_arcsec": theta_min,
        "theta_window_max_arcsec": theta_max,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points-csv", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--galaxy", required=True)
    ap.add_argument("--z-lens", type=float, default=0.3)
    ap.add_argument("--z-source", type=float, default=1.0)
    ap.add_argument("--H0", type=float, default=70.0)
    ap.add_argument("--omega-m", type=float, default=0.3)
    ap.add_argument("--zmax-factor", type=float, default=20.0)
    ap.add_argument("--nz", type=int, default=3001)
    ap.add_argument("--ratio-threshold-frac", type=float, default=0.10)
    ap.add_argument("--lambda-bup", type=float, default=1.0)
    ap.add_argument("--weight-width", type=float, default=0.6, help="log-space width of radial localization around x0")
    args = ap.parse_args()

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    report = load_report(Path(args.report_json).expanduser())
    df = load_points(Path(args.points_csv).expanduser())
    g = df[df["name"] == args.galaxy].copy()
    if g.empty:
        raise ValueError(f"Galaxy not found: {args.galaxy}")

    global_bump = get_global_bump_params(report)
    if global_bump is None:
        raise ValueError("Global bump parameters not found in report.")
    delta_d, x0, sigma = global_bump
    x_ref = get_x_ref(report, default=1.0)

    g = ensure_monotonic_radius(g)
    R_kpc, Rd_kpc, x_grid, v_bary = build_baryonic_speed_profile(g)
    d_used = d_profile_x(x_grid, delta_d, x0, sigma)
    f_d_grid = geff_factor_x(x_grid, delta_d, x0, sigma, x_ref)
    w_grid = local_weight_x(x_grid, x0=x0, width=args.weight_width)

    def g_bary_func(r):
        return g_bary_from_v(r, R_kpc, v_bary)

    def f_d_func(r):
        x = np.asarray(r, dtype=float) / max(Rd_kpc, 1e-12)
        return geff_factor_x(x, delta_d, x0, sigma, x_ref)

    def w_func(r):
        x = np.asarray(r, dtype=float) / max(Rd_kpc, 1e-12)
        return local_weight_x(x, x0=x0, width=args.weight_width)

    def g_eff_func(r):
        fd = f_d_func(r)
        ww = w_func(r)
        return g_bary_func(r) * (1.0 + args.lambda_bup * ww * (fd - 1.0))

    zmax_kpc = float(args.zmax_factor * np.max(R_kpc))
    alpha_bary = alpha_from_geff(R_kpc, g_bary_func, zmax_kpc=zmax_kpc, nz=args.nz)
    alpha_bup = alpha_from_geff(R_kpc, g_eff_func, zmax_kpc=zmax_kpc, nz=args.nz)

    D_l = angular_diameter_distance_mpc(args.z_lens, args.H0, args.omega_m)
    D_s = angular_diameter_distance_mpc(args.z_source, args.H0, args.omega_m)
    D_ls = angular_diameter_distance_between_mpc(args.z_lens, args.z_source, args.H0, args.omega_m)
    theta_arcsec = theta_arcsec_from_R_kpc(R_kpc, D_l)

    thetaE_bary = find_crossing(theta_arcsec, alpha_bary, theta_arcsec)
    thetaE_bup = find_crossing(theta_arcsec, alpha_bup, theta_arcsec)

    enhancement = alpha_bup / np.maximum(alpha_bary, 1e-12)
    delta_alpha = alpha_bup - alpha_bary
    rr = robust_ratio_stats(alpha_bup, alpha_bary, theta_arcsec, frac_threshold=args.ratio_threshold_frac)

    i_alpha_bary = int(np.nanargmax(alpha_bary))
    i_alpha_bup = int(np.nanargmax(alpha_bup))
    i_delta = int(np.nanargmax(delta_alpha))
    i_ratio = int(np.nanargmax(rr["ratio"])) if np.any(np.isfinite(rr["ratio"])) else None

    out_df = pd.DataFrame({
        "name": args.galaxy,
        "R_kpc": R_kpc,
        "theta_arcsec": theta_arcsec,
        "x_over_Rd": x_grid,
        "d_used": d_used,
        "f_d": f_d_grid,
        "w_local": w_grid,
        "v_bary_kms": v_bary,
        "alpha_bary_arcsec": alpha_bary,
        "alpha_bup_arcsec": alpha_bup,
        "delta_alpha_arcsec": delta_alpha,
        "alpha_ratio_bup_over_bary_raw": enhancement,
        "alpha_ratio_bup_over_bary_robust": rr["ratio"],
        "robust_mask": rr["mask"].astype(int),
    })
    out_df.to_csv(outdir / f"{args.galaxy}_optical_metric_v1_2_profile.csv", index=False)

    report_out = {
        "model": "bup_optical_metric_lensing_v1_2",
        "interpretation": "localized effective acceleration ray-tracing",
        "galaxy": args.galaxy,
        "cosmology": {
            "z_lens": args.z_lens,
            "z_source": args.z_source,
            "H0": args.H0,
            "omega_m": args.omega_m,
            "D_l_Mpc": D_l,
            "D_s_Mpc": D_s,
            "D_ls_Mpc": D_ls,
        },
        "bup_profile": {
            "delta_d": delta_d,
            "x0": x0,
            "sigma": sigma,
            "x_ref": x_ref,
            "d_min": 3.0 - delta_d,
            "lambda_bup": args.lambda_bup,
            "weight_width": args.weight_width,
        },
        "numerics": {
            "zmax_kpc": zmax_kpc,
            "nz": args.nz,
            "ratio_threshold_frac": args.ratio_threshold_frac,
        },
        "diagnostics": {
            "alpha_bary_max_arcsec": float(np.nanmax(alpha_bary)),
            "alpha_bup_max_arcsec": float(np.nanmax(alpha_bup)),
            "delta_alpha_max_arcsec": float(np.nanmax(delta_alpha)),
            "alpha_enhancement_max_raw": float(np.nanmax(enhancement)),
            "alpha_enhancement_mean_robust": rr["mean_ratio"],
            "alpha_enhancement_median_robust": rr["median_ratio"],
            "alpha_enhancement_max_robust": rr["max_ratio"],
            "theta_window_min_arcsec": rr["theta_window_min_arcsec"],
            "theta_window_max_arcsec": rr["theta_window_max_arcsec"],
            "thetaE_bary_arcsec": None if np.isnan(thetaE_bary) else float(thetaE_bary),
            "thetaE_bup_arcsec": None if np.isnan(thetaE_bup) else float(thetaE_bup),
            "R_alpha_bary_max_kpc": float(R_kpc[i_alpha_bary]),
            "R_alpha_bup_max_kpc": float(R_kpc[i_alpha_bup]),
            "R_delta_alpha_max_kpc": float(R_kpc[i_delta]),
            "R_ratio_max_kpc": None if i_ratio is None else float(R_kpc[i_ratio]),
            "theta_ratio_max_arcsec": None if i_ratio is None else float(theta_arcsec[i_ratio]),
        }
    }

    with open(outdir / "report.json", "w", encoding="utf-8") as f:
        json.dump(report_out, f, indent=2)

    plt.figure(figsize=(8, 5))
    plt.plot(theta_arcsec, alpha_bary, label="alpha_bary")
    plt.plot(theta_arcsec, alpha_bup, label="alpha_BuP localized")
    plt.plot(theta_arcsec, theta_arcsec, "--", label="alpha=theta")
    plt.xlabel("theta [arcsec]")
    plt.ylabel("alpha(theta) [arcsec]")
    plt.title(f"{args.galaxy} — BuP optical metric v1.2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{args.galaxy}_alpha_compare.png", dpi=180, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(R_kpc, delta_alpha)
    plt.xlabel("R [kpc]")
    plt.ylabel("alpha_BuP - alpha_bary [arcsec]")
    plt.title(f"{args.galaxy} — excess deflection v1.2")
    plt.tight_layout()
    plt.savefig(outdir / f"{args.galaxy}_delta_alpha.png", dpi=180, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(R_kpc, enhancement, label="raw ratio")
    if np.any(np.isfinite(rr["ratio"])):
        plt.plot(R_kpc, rr["ratio"], label="robust ratio")
    plt.xlabel("R [kpc]")
    plt.ylabel("alpha_BuP / alpha_bary")
    plt.title(f"{args.galaxy} — deflection enhancement v1.2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{args.galaxy}_enhancement.png", dpi=180, bbox_inches="tight")
    plt.close()

    print("=" * 72)
    print("BuP OPTICAL METRIC LENSING v1.2 — SUMMARY")
    print("=" * 72)
    print(f"Galaxy                        : {args.galaxy}")
    print(f"z_lens / z_source             : {args.z_lens:.3f} / {args.z_source:.3f}")
    print(f"d_min                         : {3.0 - delta_d:.4f}")
    print(f"lambda_bup / weight_width     : {args.lambda_bup:.3f} / {args.weight_width:.3f}")
    print(f"zmax_kpc / nz                 : {zmax_kpc:.2f} / {args.nz}")
    print(f"alpha_bary,max [arcsec]       : {np.nanmax(alpha_bary):.4f}")
    print(f"alpha_BuP,max  [arcsec]       : {np.nanmax(alpha_bup):.4f}")
    print(f"delta_alpha,max [arcsec]      : {np.nanmax(delta_alpha):.4f}")
    print(f"enhancement max raw           : {np.nanmax(enhancement):.4f}")
    print(f"enhancement mean robust       : {rr['mean_ratio']:.4f}")
    print(f"enhancement median robust     : {rr['median_ratio']:.4f}")
    print(f"theta robust window [arcsec]  : {rr['theta_window_min_arcsec']:.4f} -> {rr['theta_window_max_arcsec']:.4f}")
    print(f"thetaE_bary [arcsec]          : {'nan' if np.isnan(thetaE_bary) else f'{thetaE_bary:.4f}'}")
    print(f"thetaE_BuP  [arcsec]          : {'nan' if np.isnan(thetaE_bup) else f'{thetaE_bup:.4f}'}")
    print(f"Output dir                    : {outdir}")


if __name__ == "__main__":
    main()

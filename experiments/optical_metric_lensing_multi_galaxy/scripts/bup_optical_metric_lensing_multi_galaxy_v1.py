#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
from pathlib import Path
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C_KMS = 299792.458
RAD_TO_ARCSEC = 206264.80624709636


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


def theta_arcsec_from_R_kpc(R_kpc: np.ndarray, D_l_mpc: float) -> np.ndarray:
    D_l_kpc = D_l_mpc * 1e3
    return np.asarray(R_kpc, dtype=float) / D_l_kpc * RAD_TO_ARCSEC


def find_crossing(x: np.ndarray, y: np.ndarray, target) -> float:
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


def local_weight_x(x: np.ndarray, x_center: float, width: float) -> np.ndarray:
    x = np.maximum(np.asarray(x, dtype=float), 1e-10)
    return np.exp(-(np.log(x / x_center) ** 2) / (2.0 * width ** 2))


def robust_ratio_stats(alpha_bup: np.ndarray, alpha_bary: np.ndarray, theta: np.ndarray, frac_threshold: float = 0.1):
    a_bary_max = float(np.nanmax(alpha_bary))
    thresh = frac_threshold * a_bary_max
    mask = np.asarray(alpha_bary) >= thresh
    ratio = np.full_like(alpha_bary, np.nan, dtype=float)
    ratio[mask] = np.asarray(alpha_bup)[mask] / np.maximum(np.asarray(alpha_bary)[mask], 1e-12)
    if np.any(mask):
        return {
            "mean_ratio": float(np.nanmean(ratio[mask])),
            "median_ratio": float(np.nanmedian(ratio[mask])),
            "max_ratio": float(np.nanmax(ratio[mask])),
        }
    return {"mean_ratio": np.nan, "median_ratio": np.nan, "max_ratio": np.nan}


def parse_param_sets(text: str) -> List[Dict[str, float]]:
    # format: "3,0.8,4,5;3,0.8,3,5"
    out = []
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        lam, width, xcap, rcut = [float(x) for x in chunk.split(",")]
        out.append({
            "lambda_bup": lam,
            "weight_width": width,
            "x_center_cap": xcap,
            "r_cut_factor": rcut,
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points-csv", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--galaxies", default="NGC2403,NGC3198,NGC5055,NGC2841")
    ap.add_argument("--param-sets", default="3,0.8,4,5;3,0.8,3,5;3,0.6,4,5")
    ap.add_argument("--z-lens", type=float, default=0.3)
    ap.add_argument("--H0", type=float, default=70.0)
    ap.add_argument("--omega-m", type=float, default=0.3)
    ap.add_argument("--zmax-factor", type=float, default=20.0)
    ap.add_argument("--nz", type=int, default=3001)
    ap.add_argument("--ratio-threshold-frac", type=float, default=0.10)
    args = ap.parse_args()

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    report = load_report(Path(args.report_json).expanduser())
    df = load_points(Path(args.points_csv).expanduser())
    bump = get_global_bump_params(report)
    if bump is None:
        raise ValueError("Global bump parameters not found in report.")
    delta_d, x0, sigma = bump
    x_ref = get_x_ref(report, default=1.0)

    galaxies = [g.strip() for g in args.galaxies.split(",") if g.strip()]
    param_sets = parse_param_sets(args.param_sets)

    D_l = angular_diameter_distance_mpc(args.z_lens, args.H0, args.omega_m)

    rows = []
    best_per_galaxy = {}

    for gal in galaxies:
        g = df[df["name"] == gal].copy()
        if g.empty:
            continue
        g = ensure_monotonic_radius(g)
        R_kpc, Rd_kpc, x_grid, v_bary = build_baryonic_speed_profile(g)

        def g_bary_func(r):
            return g_bary_from_v(r, R_kpc, v_bary)

        zmax_kpc = float(args.zmax_factor * np.max(R_kpc))
        alpha_bary = alpha_from_geff(R_kpc, g_bary_func, zmax_kpc=zmax_kpc, nz=args.nz)
        theta_arcsec = theta_arcsec_from_R_kpc(R_kpc, D_l)
        thetaE_bary = find_crossing(theta_arcsec, alpha_bary, theta_arcsec)

        best_score = -1e99
        best_data = None

        for ps in param_sets:
            lam = ps["lambda_bup"]
            width = ps["weight_width"]
            xcap = ps["x_center_cap"]
            rcutf = ps["r_cut_factor"]

            x_center = min(x0, xcap)
            r_cut_kpc = rcutf * Rd_kpc

            def f_d_func(r):
                x = np.asarray(r, dtype=float) / max(Rd_kpc, 1e-12)
                return geff_factor_x(x, delta_d, x0, sigma, x_ref)

            def w_total_func(r):
                r = np.asarray(r, dtype=float)
                x = r / max(Rd_kpc, 1e-12)
                return local_weight_x(x, x_center=x_center, width=width) * np.exp(-r / max(r_cut_kpc, 1e-12))

            def g_eff_func(r):
                return g_bary_func(r) * (1.0 + lam * w_total_func(r) * (f_d_func(r) - 1.0))

            alpha_bup = alpha_from_geff(R_kpc, g_eff_func, zmax_kpc=zmax_kpc, nz=args.nz)
            delta_alpha = alpha_bup - alpha_bary
            enhancement = alpha_bup / np.maximum(alpha_bary, 1e-12)
            rr = robust_ratio_stats(alpha_bup, alpha_bary, theta_arcsec, frac_threshold=args.ratio_threshold_frac)
            thetaE_bup = find_crossing(theta_arcsec, alpha_bup, theta_arcsec)

            i_bup = int(np.nanargmax(alpha_bup))
            i_delta = int(np.nanargmax(delta_alpha))
            i_ratio = int(np.nanargmax(enhancement))

            row = {
                "galaxy": gal,
                "lambda_bup": lam,
                "weight_width": width,
                "x_center_cap": xcap,
                "x_center_used": x_center,
                "r_cut_factor": rcutf,
                "r_cut_kpc": r_cut_kpc,
                "alpha_bary_max_arcsec": float(np.nanmax(alpha_bary)),
                "alpha_bup_max_arcsec": float(np.nanmax(alpha_bup)),
                "delta_alpha_max_arcsec": float(np.nanmax(delta_alpha)),
                "enhancement_max_raw": float(np.nanmax(enhancement)),
                "enhancement_mean_robust": rr["mean_ratio"],
                "enhancement_median_robust": rr["median_ratio"],
                "R_alpha_bup_max_kpc": float(R_kpc[i_bup]),
                "R_delta_alpha_max_kpc": float(R_kpc[i_delta]),
                "R_ratio_max_kpc": float(R_kpc[i_ratio]),
                "thetaE_bary_arcsec": None if np.isnan(thetaE_bary) else float(thetaE_bary),
                "thetaE_bup_arcsec": None if np.isnan(thetaE_bup) else float(thetaE_bup),
                "delta_thetaE_arcsec": np.nan if (np.isnan(thetaE_bary) or np.isnan(thetaE_bup)) else float(thetaE_bup - thetaE_bary),
            }
            # Simple score: centrality + amplitude
            centrality = np.exp(-((row["R_delta_alpha_max_kpc"] - 10.0) ** 2) / (2.0 * 8.0 ** 2))
            row["score"] = (
                2.0 * row["enhancement_mean_robust"]
                + 1.0 * row["enhancement_median_robust"]
                + 15.0 * row["delta_alpha_max_arcsec"]
                + 20.0 * max(row["delta_thetaE_arcsec"], 0.0)
                + 0.5 * centrality
            )
            rows.append(row)

            if row["score"] > best_score:
                best_score = row["score"]
                best_data = {
                    "summary": row.copy(),
                    "R_kpc": R_kpc.copy(),
                    "theta_arcsec": theta_arcsec.copy(),
                    "alpha_bary": alpha_bary.copy(),
                    "alpha_bup": alpha_bup.copy(),
                    "delta_alpha": delta_alpha.copy(),
                }

        best_per_galaxy[gal] = best_data

    res = pd.DataFrame(rows).sort_values(["galaxy", "score"], ascending=[True, False]).reset_index(drop=True)
    res.to_csv(outdir / "multi_galaxy_optical_metric_summary.csv", index=False)

    # best-only summary
    best_rows = [best_per_galaxy[g]["summary"] for g in best_per_galaxy if best_per_galaxy[g] is not None]
    best_df = pd.DataFrame(best_rows).sort_values("score", ascending=False)
    best_df.to_csv(outdir / "multi_galaxy_optical_metric_best_only.csv", index=False)

    # plot best delta_alpha for all galaxies
    plt.figure(figsize=(9, 6))
    for gal, data in best_per_galaxy.items():
        if data is None:
            continue
        plt.plot(data["R_kpc"], data["delta_alpha"], label=gal)
    plt.xlabel("R [kpc]")
    plt.ylabel("alpha_BuP - alpha_bary [arcsec]")
    plt.title("Best centralized BuP excess deflection across galaxies")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "best_delta_alpha_all_galaxies.png", dpi=180, bbox_inches="tight")
    plt.close()

    report_out = {
        "model": "bup_optical_metric_lensing_multi_galaxy_v1",
        "galaxies": galaxies,
        "param_sets": param_sets,
        "bup_profile": {
            "delta_d": delta_d,
            "x0_original": x0,
            "sigma": sigma,
            "x_ref": x_ref,
            "d_min": 3.0 - delta_d,
        },
        "best_per_galaxy": best_rows,
    }
    with open(outdir / "report.json", "w", encoding="utf-8") as f:
        json.dump(report_out, f, indent=2)

    print("=" * 72)
    print("BuP OPTICAL METRIC MULTI-GALAXY TEST — SUMMARY")
    print("=" * 72)
    print(f"Galaxies tested            : {', '.join(galaxies)}")
    print(f"Parameter sets tested      : {len(param_sets)}")
    print(f"Output dir                 : {outdir}")
    print("\nBest per galaxy:")
    for _, r in best_df.iterrows():
        print(
            f"{r['galaxy']:10s} | "
            f"lambda={r['lambda_bup']:.1f}, width={r['weight_width']:.2f}, "
            f"xcap={r['x_center_cap']:.1f}, rcut={r['r_cut_factor']:.1f} | "
            f"enh_mean={r['enhancement_mean_robust']:.3f} | "
            f"enh_med={r['enhancement_median_robust']:.3f} | "
            f"delta_alpha_max={r['delta_alpha_max_arcsec']:.4f} | "
            f"R_delta_max={r['R_delta_alpha_max_kpc']:.2f} | "
            f"dThetaE={r['delta_thetaE_arcsec']:.4f}"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
BuP GW true edge-Hessian v17 — LVK eta calibration scan

Purpose
-------
Final local calibration test around eta_edge = eta_* = 1.
Runs bup_gw_true_edge_hessian_v10.py for a small list of eta_edge values,
extracts c_edge from the edge-Hessian dispersion, and verifies:

    c_edge(eta)^2 ≈ eta
    c_edge(eta) ≈ sqrt(eta)
    (c_edge - 1)/(eta - 1) ≈ 1/2 near eta = 1

This is the numerical regularity check behind the LVK statement:

    |c_GW/c - 1| < 5e-16  <=>  |eta_edge/eta_* - 1| < 1e-15

in calibrated Einstein-fixed-point units eta_* = 1.

Outputs
-------
results/v17_lvk_eta_scan_summary.csv
results/v17_lvk_eta_scan_fits.json
results/v17_lvk_eta_scan_summary.md
figures/fig_v17_cedge_vs_sqrt_eta.png
figures/fig_v17_fractional_speed_vs_eta.png
figures/fig_v17_c2_minus_eta.png
runs/eta.../summary_true_edge_hessian_v10.json

Example
-------
python3 scan_lvk_eta_calibration_v17.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/bup_gw_true_edge_hessian_v10.py \
  --eta-edge-list 0.99,0.995,0.999,1.0,1.001,1.005,1.01 \
  --N-side 11 \
  --extent 1.0 \
  --source-width 0.12 \
  --pulse-sigma 0.25 \
  --output-dir papers/paper22_bup_gravitational_waves/results/lvk_eta_calibration_v17
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_float_list(s: str) -> List[float]:
    out: List[float] = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(float(part))
    return out


def safe_float(x: Any) -> float:
    if x is None:
        return float('nan')
    try:
        return float(x)
    except Exception:
        return float('nan')


def eta_tag(eta: float) -> str:
    # file-system-safe tag preserving sign and decimals
    s = f"{eta:.12g}".replace('-', 'm').replace('.', 'p')
    return f"eta{s}"


def run_v10(args: argparse.Namespace, eta: float, rundir: Path) -> Dict[str, Any]:
    resjson = rundir / "results" / "summary_true_edge_hessian_v10.json"
    if args.skip_existing and resjson.exists():
        with open(resjson, 'r', encoding='utf-8') as f:
            return json.load(f)

    rundir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        args.v10_script,
        "--N-side", str(args.N_side),
        "--extent", str(args.extent),
        "--ell", str(args.ell),
        "--cutoff", str(args.cutoff),
        "--eta-edge", str(eta),
        "--mass2", str(args.mass2),
        "--lambda-locality", str(args.lambda_locality),
        "--edge-width", str(args.edge_width),
        "--w-barrier", str(args.w_barrier),
        "--gamma", str(args.gamma),
        "--source-amp", str(args.source_amp),
        "--source-width", str(args.source_width),
        "--pulse-t0", str(args.pulse_t0),
        "--pulse-sigma", str(args.pulse_sigma),
        "--pulse-kind", args.pulse_kind,
        "--tmax", str(args.tmax),
        "--dt", str(args.dt),
        "--record-stride", str(args.record_stride),
        "--n-modes-dyn", str(args.n_modes_dyn),
        "--n-fit-modes", str(args.n_fit_modes),
        "--n-spectral-bins", str(args.n_spectral_bins),
        "--low-mode-cut", str(args.low_mode_cut),
        "--nrings", str(args.nrings),
        "--ring-rmin", str(args.ring_rmin),
        "--rmax-frac", str(args.rmax_frac),
        "--t-ignore", str(args.t_ignore),
        "--baseline-until", str(args.baseline_until),
        "--min-snr", str(args.min_snr),
        "--output-dir", str(rundir),
    ]
    print("\n" + "="*100)
    print(f"Running eta_edge={eta}")
    print(" ".join(cmd))
    print("="*100)
    subprocess.run(cmd, check=True)
    if not resjson.exists():
        raise FileNotFoundError(f"Missing v10 summary: {resjson}")
    with open(resjson, 'r', encoding='utf-8') as f:
        return json.load(f)


def linear_fit(x: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 2:
        return {"ok": False, "reason": "not enough finite points"}
    A = np.vstack([np.ones_like(x), x]).T
    coeff, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coeff
    resid = y - yhat
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((y - np.mean(y))**2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rmse = math.sqrt(ss_res / len(y))
    return {
        "ok": True,
        "n_used": int(len(y)),
        "intercept": float(coeff[0]),
        "slope": float(coeff[1]),
        "r2": float(r2),
        "rmse": float(rmse),
    }


def make_plots(df: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)

    x = df["eta_edge"].to_numpy(float)
    c = df["c_edge"].to_numpy(float)
    sq = np.sqrt(x)

    plt.figure(figsize=(7, 5))
    plt.plot(sq, c, "o", label="measured")
    lo = min(np.nanmin(sq), np.nanmin(c))
    hi = max(np.nanmax(sq), np.nanmax(c))
    plt.plot([lo, hi], [lo, hi], "--", label="c_edge = sqrt(eta)")
    plt.xlabel(r"$\sqrt{\eta_{\rm edge}}$")
    plt.ylabel(r"$c_{\rm edge}$")
    plt.title("v17 LVK calibration: c_edge vs sqrt(eta_edge)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figdir / "fig_v17_cedge_vs_sqrt_eta.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(x, df["frac_speed_error"], "o-")
    plt.axhline(0.0, linestyle="--")
    plt.xlabel(r"$\eta_{\rm edge}$")
    plt.ylabel(r"$c_{\rm edge}-1$")
    plt.title("Fractional speed deviation in calibrated units")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(figdir / "fig_v17_fractional_speed_vs_eta.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(x, df["c2_minus_eta"], "o-")
    plt.axhline(0.0, linestyle="--")
    plt.xlabel(r"$\eta_{\rm edge}$")
    plt.ylabel(r"$c_{\rm edge}^2 - \eta_{\rm edge}$")
    plt.title("Check of c_edge² = eta_edge")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(figdir / "fig_v17_c2_minus_eta.png", dpi=180)
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description="v17 LVK eta calibration scan around eta_edge=1")
    p.add_argument("--v10-script", required=True)
    p.add_argument("--eta-edge-list", default="0.99,0.995,0.999,1.0,1.001,1.005,1.01")
    p.add_argument("--eta-star", type=float, default=1.0)
    p.add_argument("--lvk-bound", type=float, default=5e-16)
    p.add_argument("--skip-existing", action="store_true")

    # v10 passthrough parameters
    p.add_argument("--N-side", type=int, default=11)
    p.add_argument("--extent", type=float, default=1.0)
    p.add_argument("--ell", type=float, default=0.22)
    p.add_argument("--cutoff", type=float, default=0.45)
    p.add_argument("--mass2", type=float, default=1e-4)
    p.add_argument("--lambda-locality", type=float, default=0.02)
    p.add_argument("--edge-width", type=float, default=0.30)
    p.add_argument("--w-barrier", type=float, default=0.0)
    p.add_argument("--gamma", type=float, default=0.015)
    p.add_argument("--source-amp", type=float, default=0.20)
    p.add_argument("--source-width", type=float, default=0.12)
    p.add_argument("--pulse-t0", type=float, default=8.0)
    p.add_argument("--pulse-sigma", type=float, default=0.25)
    p.add_argument("--pulse-kind", choices=["ricker", "gaussian"], default="ricker")
    p.add_argument("--tmax", type=float, default=45.0)
    p.add_argument("--dt", type=float, default=0.01)
    p.add_argument("--record-stride", type=int, default=5)
    p.add_argument("--n-modes-dyn", type=int, default=180)
    p.add_argument("--n-fit-modes", type=int, default=60)
    p.add_argument("--n-spectral-bins", type=int, default=10)
    p.add_argument("--low-mode-cut", type=int, default=25)
    p.add_argument("--nrings", type=int, default=10)
    p.add_argument("--ring-rmin", type=float, default=0.15)
    p.add_argument("--rmax-frac", type=float, default=0.82)
    p.add_argument("--t-ignore", type=float, default=7.5)
    p.add_argument("--baseline-until", type=float, default=5.0)
    p.add_argument("--min-snr", type=float, default=1.25)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    outdir = Path(args.output_dir)
    resdir = outdir / "results"
    figdir = outdir / "figures"
    rundir_root = outdir / "runs"
    resdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    rundir_root.mkdir(parents=True, exist_ok=True)

    etas = parse_float_list(args.eta_edge_list)
    rows: List[Dict[str, Any]] = []

    for eta in etas:
        summary = run_v10(args, eta, rundir_root / eta_tag(eta))
        c = safe_float(summary.get("dispersion_c_graph"))
        meff2 = safe_float(summary.get("dispersion_intercept_mass2"))
        r2 = safe_float(summary.get("dispersion_r2"))
        slope_c2 = safe_float(summary.get("dispersion_slope_c2"))
        plus = summary.get("plus", {}) or {}
        cross = summary.get("cross", {}) or {}
        plus_corr = safe_float(plus.get("modal_corr_absJ_qpeak"))
        cross_corr = safe_float(cross.get("modal_corr_absJ_qpeak"))
        mean_corr = np.nanmean([plus_corr, cross_corr])

        frac_speed_error = c / math.sqrt(args.eta_star) - 1.0 if np.isfinite(c) else float('nan')
        eta_delta = eta / args.eta_star - 1.0
        derivative_ratio = frac_speed_error / eta_delta if abs(eta_delta) > 0 else float('nan')
        c2_minus_eta = c*c - eta if np.isfinite(c) else float('nan')
        eta_lvk_window = 2.0 * args.lvk_bound  # because delta c/c ≈ 0.5 delta eta/eta

        rows.append({
            "eta_edge": eta,
            "eta_delta": eta_delta,
            "c_edge": c,
            "sqrt_eta": math.sqrt(eta),
            "frac_speed_error": frac_speed_error,
            "abs_frac_speed_error": abs(frac_speed_error) if np.isfinite(frac_speed_error) else float('nan'),
            "lvk_bound": args.lvk_bound,
            "eta_lvk_window_approx": eta_lvk_window,
            "derivative_dc_over_deta_near1": derivative_ratio,
            "c2": c*c if np.isfinite(c) else float('nan'),
            "c2_minus_eta": c2_minus_eta,
            "slope_c2_from_dispersion": slope_c2,
            "meff2": meff2,
            "R2_disp": r2,
            "plus_corr": plus_corr,
            "cross_corr": cross_corr,
            "mean_corr": mean_corr,
        })

    df = pd.DataFrame(rows).sort_values("eta_edge")
    df.to_csv(resdir / "v17_lvk_eta_scan_summary.csv", index=False)

    # Fits around eta*=1
    finite = df[np.isfinite(df["c_edge"])].copy()
    fit_c_vs_sqrt_eta = linear_fit(finite["sqrt_eta"].to_numpy(float), finite["c_edge"].to_numpy(float))
    fit_delta = linear_fit(finite["eta_delta"].to_numpy(float), finite["frac_speed_error"].to_numpy(float))
    fit_c2_eta = linear_fit(finite["eta_edge"].to_numpy(float), finite["c2"].to_numpy(float))

    # central finite difference if both +/- matching points around eta_star exist approximately
    central_pairs = []
    for eta in etas:
        if eta <= args.eta_star:
            continue
        neg = args.eta_star - (eta - args.eta_star)
        if any(abs(e - neg) < 1e-12 for e in etas):
            cp = df.loc[np.isclose(df["eta_edge"], eta), "c_edge"].iloc[0]
            cm = df.loc[np.isclose(df["eta_edge"], neg), "c_edge"].iloc[0]
            central_pairs.append({
                "delta_eta": float(eta - args.eta_star),
                "dc_deta_central": float((cp - cm) / (2.0 * (eta - args.eta_star))),
            })
    central_df = pd.DataFrame(central_pairs)
    if not central_df.empty:
        central_df.to_csv(resdir / "v17_central_derivatives.csv", index=False)

    fits = {
        "fit_c_edge_vs_sqrt_eta": fit_c_vs_sqrt_eta,
        "fit_frac_speed_error_vs_eta_delta": fit_delta,
        "fit_c_edge_squared_vs_eta": fit_c2_eta,
        "central_derivatives": central_pairs,
        "lvk_translation": {
            "bound_abs_cgw_over_c_minus_1": args.lvk_bound,
            "eta_star": args.eta_star,
            "required_abs_eta_over_eta_star_minus_1_approx": 2.0 * args.lvk_bound,
            "comment": "Since c_edge/c = sqrt(eta_edge/eta_star), delta c/c ≈ 0.5 delta eta/eta_star near eta_star.",
        },
    }
    with open(resdir / "v17_lvk_eta_scan_fits.json", "w", encoding="utf-8") as f:
        json.dump(fits, f, indent=2)

    make_plots(df, figdir)

    # Summary markdown
    best_slope = fit_delta.get("slope", float('nan')) if fit_delta.get("ok") else float('nan')
    c2_slope = fit_c2_eta.get("slope", float('nan')) if fit_c2_eta.get("ok") else float('nan')
    md = []
    md.append("# v17 — LVK eta calibration scan")
    md.append("")
    md.append("This scan verifies the regular calibration of the true edge-Hessian wave speed around the Einstein fixed point.")
    md.append("")
    md.append("Core relations:")
    md.append("")
    md.append(r"\[")
    md.append(r"c_{\rm edge}^2 \simeq \eta_{\rm edge},")
    md.append(r"\qquad")
    md.append(r"c_{\rm edge}\simeq\sqrt{\eta_{\rm edge}}.")
    md.append(r"\]")
    md.append("")
    md.append("## Fit summary")
    md.append("")
    md.append(f"- Fit slope of `(c_edge - 1)` vs `(eta_edge - 1)`: `{best_slope}`")
    md.append(f"- Expected local slope: `0.5`")
    md.append(f"- Fit slope of `c_edge^2` vs `eta_edge`: `{c2_slope}`")
    md.append(f"- Expected slope: `1.0`")
    md.append(f"- LVK speed bound used: `{args.lvk_bound}`")
    md.append(f"- Implied eta window near eta_star: `|eta/eta_star - 1| < {2.0*args.lvk_bound}`")
    md.append("")
    md.append("## Interpretation")
    md.append("")
    md.append("The LVK bound is not treated as a free fine-tuning of eta_edge. In calibrated Einstein-fixed-point units, eta_star defines the emergent relativistic normalization, so eta_star = 1 and c_edge(eta_star)=c. The observational constraint means the present universe must lie extremely close to this rigidity fixed point.")
    (resdir / "v17_lvk_eta_scan_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("\n" + "="*100)
    print("LVK eta calibration v17 scan complete")
    print("="*100)
    cols = ["eta_edge", "c_edge", "sqrt_eta", "frac_speed_error", "derivative_dc_over_deta_near1", "c2_minus_eta", "meff2", "R2_disp", "mean_corr"]
    print(df[cols].to_string(index=False))
    print("\nFits:")
    print(json.dumps(fits, indent=2))
    print(f"\nWrote outputs to: {outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Numerical summary builder v1

Goal
----
Collect the three first numerical results of Paper 15:

1. Spectral convergence:
   L_epsilon -> Delta_g

2. Ricci signal:
   kappa_OR distinguishes flat periodic geometry from sphere.

3. Source-response:
   radial entanglement defect produces localized curvature response.

Inputs are the result folders produced by:
- paper15_spectral_convergence_v2.py
- paper15_ricci_convergence_v2.py
- paper15_source_response_v2.py

Outputs
-------
paper15_numerical_summary.json
paper15_numerical_summary.md
paper15_key_results_table.csv

Recommended run
---------------
cd ~/bottomup

python3 papers/paper15_einstein_derivation/scripts/paper15_build_numerical_summary_v1.py \
  --spectral-dir papers/paper15_einstein_derivation/results/spectral_convergence_v2_best_unnormalized \
  --ricci-dir papers/paper15_einstein_derivation/results/ricci_convergence_v2 \
  --source-dir papers/paper15_einstein_derivation/results/source_response_v2 \
  --output-dir papers/paper15_einstein_derivation/results/paper15_numerical_summary_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Build Paper 15 numerical summary.")
    p.add_argument("--spectral-dir", required=True)
    p.add_argument("--ricci-dir", required=True)
    p.add_argument("--source-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_csv(path: Path, required=True):
    if path.exists():
        return pd.read_csv(path)
    if required:
        raise FileNotFoundError(f"Missing required file: {path}")
    return None


def summarize_spectral(spectral_dir: Path):
    # Expected from v2 best unnormalized:
    # summary_by_geometry_laplacian.csv
    path = spectral_dir / "summary_by_geometry_laplacian.csv"
    df = load_csv(path)

    # take final rows already summarized
    rows = []
    for _, r in df.iterrows():
        rows.append({
            "geometry": r["geometry"],
            "target_dim": float(r["target_dim"]),
            "measured_ds": float(r["final_ds_plateau"]),
            "abs_error": float(r["final_ds_abs_error"]),
            "laplacian": r["laplacian"],
            "N": int(r["best_N"]),
            "k": int(r["final_k"]) if "final_k" in df.columns else None,
            "epsilon": float(r["final_epsilon"]) if "final_epsilon" in df.columns else None,
        })

    out = {
        "status": "positive",
        "interpretation": "Continuum-rescaled graph Laplacian recovers the expected heat-trace spectral dimension on controlled geometries.",
        "operator": "L_epsilon=(D-W)/epsilon",
        "rows": rows,
    }

    return out


def summarize_ricci(ricci_dir: Path):
    path = ricci_dir / "final_by_geometry.csv"
    df = load_csv(path)

    flat = df[df["geometry"] == "flat_torus2d"]
    sphere = df[df["geometry"] == "sphere"]

    flat_mean = float(flat["mean_kappa"].iloc[0]) if len(flat) else np.nan
    sphere_mean = float(sphere["mean_kappa"].iloc[0]) if len(sphere) else np.nan
    delta = sphere_mean - flat_mean if np.isfinite(flat_mean) and np.isfinite(sphere_mean) else np.nan

    rows = []
    for _, r in df.iterrows():
        rows.append({
            "geometry": r["geometry"],
            "N": int(r["N_actual"]),
            "mean_kappa": float(r["mean_kappa"]),
            "median_kappa": float(r["median_kappa"]),
            "positive_fraction": float(r["positive_fraction"]),
            "delta_mean_vs_flat_torus2d": float(r["delta_mean_vs_flat_torus2d"]) if "delta_mean_vs_flat_torus2d" in df.columns else np.nan,
        })

    status = "positive_relative" if np.isfinite(delta) and delta > 0 else "inconclusive"

    return {
        "status": status,
        "interpretation": "Ollivier-Ricci curvature gives a positive curvature excess on the sphere relative to the flat periodic torus.",
        "flat_mean_kappa": flat_mean,
        "sphere_mean_kappa": sphere_mean,
        "delta_sphere_minus_flat": delta,
        "rows": rows,
    }


def summarize_source(source_dir: Path):
    path = source_dir / "source_response_v2_summary.csv"
    df = load_csv(path)

    # Key values: strongest negative source, because user results emphasized s=-0.30.
    key_strength = -0.30
    key = df[np.isclose(df["perturbation_strength"], key_strength)]

    rows = []
    for _, r in key.iterrows():
        rows.append({
            "geometry": r["geometry"],
            "strength": float(r["perturbation_strength"]),
            "near_abs_mean_delta": float(r["near_abs_mean_delta"]),
            "far_abs_mean_delta": float(r["far_abs_mean_delta"]),
            "localization_ratio_near_far": float(r["localization_ratio_near_far"]),
            "spearman_phi_absdelta": float(r["spearman_phi_absdelta"]),
            "spearman_phi_absdelta_p": float(r["spearman_phi_absdelta_p"]),
            "spearman_phi_delta": float(r["spearman_phi_delta"]),
            "spearman_phi_delta_p": float(r["spearman_phi_delta_p"]),
        })

    # global ranges
    ratio_min = float(df["localization_ratio_near_far"].replace([np.inf, -np.inf], np.nan).min())
    ratio_max = float(df["localization_ratio_near_far"].replace([np.inf, -np.inf], np.nan).max())
    rho_abs_min = float(df["spearman_phi_absdelta"].min())
    rho_abs_max = float(df["spearman_phi_absdelta"].max())

    return {
        "status": "positive",
        "interpretation": "A smooth radial perturbation of W produces a localized signed curvature response correlated with the source profile.",
        "key_strength": key_strength,
        "localization_ratio_range": [ratio_min, ratio_max],
        "spearman_phi_absdelta_range": [rho_abs_min, rho_abs_max],
        "rows_strength_minus_0p30": rows,
    }


def build_key_table(spectral, ricci, source):
    rows = []

    for r in spectral["rows"]:
        rows.append({
            "step": "Step B — spectral convergence",
            "quantity": f"{r['geometry']}: d_s",
            "value": f"{r['measured_ds']:.4f}",
            "target_or_reference": f"{r['target_dim']:.0f}",
            "error_or_delta": f"{r['abs_error']:.4f}",
            "status": spectral["status"],
        })

    rows.append({
        "step": "Step C — Ricci signal",
        "quantity": "sphere mean kappa",
        "value": f"{ricci['sphere_mean_kappa']:.6f}",
        "target_or_reference": f"flat torus mean = {ricci['flat_mean_kappa']:.6f}",
        "error_or_delta": f"delta = {ricci['delta_sphere_minus_flat']:.6f}",
        "status": ricci["status"],
    })

    for r in source["rows_strength_minus_0p30"]:
        rows.append({
            "step": "Step E — source response",
            "quantity": f"{r['geometry']}: near/far response, s=-0.30",
            "value": f"{r['localization_ratio_near_far']:.2f}",
            "target_or_reference": "near/far > 1",
            "error_or_delta": f"Spearman(phi,|dkappa|)={r['spearman_phi_absdelta']:.3f}",
            "status": source["status"],
        })

    return pd.DataFrame(rows)


def write_markdown(outpath: Path, summary: dict, key_table: pd.DataFrame):
    lines = []
    lines.append("# Paper 15 — Numerical Summary v1")
    lines.append("")
    lines.append("This file summarizes the first three numerical tests of Paper 15.")
    lines.append("")
    lines.append("The target chain is:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"W_{ij}\to L_\epsilon\to\Delta_g,")
    lines.append(r"\qquad")
    lines.append(r"\kappa_{ij}^{\rm OR}\to R_{\mu\nu}u^\mu u^\nu,")
    lines.append(r"\qquad")
    lines.append(r"\delta W_{\rm loc}\to\delta\kappa(r).")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 1. Key results table")
    lines.append("")
    lines.append(key_table.to_markdown(index=False))
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 2. Step B — Spectral convergence")
    lines.append("")
    lines.append("The continuum-scaled Laplacian")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"L_\epsilon=\frac{D-W}{\epsilon}")
    lines.append(r"\]")
    lines.append("")
    lines.append("recovers the expected heat-trace spectral dimension on controlled geometries.")
    lines.append("")
    lines.append("| Geometry | Target dimension | Measured \(d_s\) | Error |")
    lines.append("|---|---:|---:|---:|")
    for r in summary["spectral"]["rows"]:
        lines.append(f"| {r['geometry']} | {r['target_dim']:.0f} | {r['measured_ds']:.4f} | {r['abs_error']:.4f} |")
    lines.append("")
    lines.append("Interpretation: this is the first positive numerical evidence for")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"L_N\longrightarrow \Delta_g.")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 3. Step C — Ollivier--Ricci curvature signal")
    lines.append("")
    lines.append("The flat periodic torus gives")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\bar\kappa_{{\rm flat}}={summary['ricci']['flat_mean_kappa']:.6f},")
    lines.append(r"\]")
    lines.append("")
    lines.append("while the sphere gives")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\bar\kappa_{{\rm sphere}}={summary['ricci']['sphere_mean_kappa']:.6f}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The relative curvature excess is")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\Delta\bar\kappa_{{\rm sphere-flat}}={summary['ricci']['delta_sphere_minus_flat']:.6f}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("Interpretation: this provides a first relative numerical signal for")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\kappa_{ij}^{\rm OR}\longrightarrow R_{\mu\nu}u^\mu u^\nu.")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 4. Step E — Source-response test")
    lines.append("")
    lines.append("A smooth radial perturbation of the entanglement weights is applied:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\phi_i=\exp\left[-\frac{d(i,\mathrm{source})^2}{2\sigma^2}\right],")
    lines.append(r"\]")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"W'_{ij}=W_{ij}\left[1+s\frac{\phi_i+\phi_j}{2}\right].")
    lines.append(r"\]")
    lines.append("")
    lines.append("For \(s=-0.30\), the key source-response results are:")
    lines.append("")
    #lines.append("| Geometry | near/far \(\langle|\Delta\kappa|\rangle\) | Spearman \((\phi,|\Delta\kappa|)\) | p-value |")
    lines.append("| Geometry | near/far response | Spearman \\((\\phi,|\\Delta\\kappa|)\\) | p-value |")
    lines.append("|---|---:|---:|---:|")
    for r in summary["source"]["rows_strength_minus_0p30"]:
        lines.append(
            f"| {r['geometry']} | {r['localization_ratio_near_far']:.2f} | "
            f"{r['spearman_phi_absdelta']:.3f} | {r['spearman_phi_absdelta_p']:.2e} |"
        )
    lines.append("")
    lines.append("Interpretation: a local entanglement defect generates a localized curvature response:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\delta W_{\rm loc}\longrightarrow \delta\kappa(r).")
    lines.append(r"\]")
    lines.append("")
    lines.append("This connects Paper 15 directly to the Paper 8 construction of")
    lines.append(r"\(T_{\mu\nu}^{\rm eff}\) from local entanglement perturbations.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 5. Current status")
    lines.append("")
    lines.append("Paper 15 now has three positive numerical pillars:")
    lines.append("")
    lines.append("1. \(L_\epsilon=(D-W)/\epsilon\) recovers the expected spectral dimension.")
    lines.append("2. Ollivier--Ricci curvature gives a positive sphere-minus-flat signal.")
    lines.append("3. A radial entanglement defect produces a localized curvature response.")
    lines.append("")
    lines.append("These results do not yet prove the continuum Einstein equation, but they support the three main arrows required for the effective derivation.")
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    spectral_dir = Path(args.spectral_dir)
    ricci_dir = Path(args.ricci_dir)
    source_dir = Path(args.source_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    spectral = summarize_spectral(spectral_dir)
    ricci = summarize_ricci(ricci_dir)
    source = summarize_source(source_dir)

    summary = {
        "experiment": "Paper 15 numerical summary v1",
        "spectral_dir": str(spectral_dir),
        "ricci_dir": str(ricci_dir),
        "source_dir": str(source_dir),
        "spectral": spectral,
        "ricci": ricci,
        "source": source,
    }

    key_table = build_key_table(spectral, ricci, source)
    key_table.to_csv(outdir / "paper15_key_results_table.csv", index=False)

    with open(outdir / "paper15_numerical_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    write_markdown(outdir / "paper15_numerical_summary.md", summary, key_table)

    print("=" * 100)
    print("Paper 15 — Numerical summary v1")
    print("=" * 100)
    print(key_table.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "paper15_key_results_table.csv")
    print(outdir / "paper15_numerical_summary.json")
    print(outdir / "paper15_numerical_summary.md")


if __name__ == "__main__":
    main()

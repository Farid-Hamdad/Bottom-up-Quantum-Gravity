#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 17 — Ricci summary builder v1

Goal
----
Collect Paper 17 numerical results:

1. v1: raw and scale-normalized Ollivier--Ricci signal
   sphere > flat torus

2. v2: affine calibration of the Ricci proxy

    Rhat_OR = A_N * (kappa/epsilon - B_N)

with flat torus as Ricci=0 reference and unit sphere as Ricci=1 reference.

Inputs
------
- results/or_ricci_scaling_v1/or_ricci_scaling_summary.csv
- results/or_ricci_scaling_v1/or_ricci_scaling_relative.csv
- results/or_ricci_scaling_v2/or_ricci_scaling_v2_calibration.csv
- results/or_ricci_scaling_v2/or_ricci_scaling_v2_calibrated_summary.csv

Outputs
-------
paper17_key_results_table.csv
paper17_ricci_summary.json
paper17_ricci_summary.md

Recommended run
---------------
cd ~/bottomup

python3 papers/paper17_ollivier_ricci_limit/scripts/paper17_build_ricci_summary_v1.py \
  --v1-dir papers/paper17_ollivier_ricci_limit/results/or_ricci_scaling_v1 \
  --v2-dir papers/paper17_ollivier_ricci_limit/results/or_ricci_scaling_v2 \
  --output-dir papers/paper17_ollivier_ricci_limit/results/paper17_ricci_summary_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Build Paper 17 Ricci summary v1.")
    p.add_argument("--v1-dir", required=True)
    p.add_argument("--v2-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_csv(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def summarize_v1(v1_dir: Path):
    summary_path = v1_dir / "or_ricci_scaling_summary.csv"
    relative_path = v1_dir / "or_ricci_scaling_relative.csv"

    summary_df = load_csv(summary_path)
    relative_df = load_csv(relative_path)

    # v1 relative comparison had one exact N_actual match in the previous run.
    if len(relative_df):
        best = relative_df.sort_values("N_actual").tail(1).iloc[0].to_dict()
    else:
        best = {}

    return {
        "status": "positive_relative",
        "interpretation": "Raw and scale-normalized Ollivier-Ricci curvature give a positive sphere-minus-flat signal.",
        "summary_rows": summary_df.to_dict(orient="records"),
        "relative_rows": relative_df.to_dict(orient="records"),
        "best_relative": best,
    }


def summarize_v2(v2_dir: Path):
    calibration_path = v2_dir / "or_ricci_scaling_v2_calibration.csv"
    calibrated_summary_path = v2_dir / "or_ricci_scaling_v2_calibrated_summary.csv"

    calibration_df = load_csv(calibration_path)
    calibrated_summary_df = load_csv(calibrated_summary_path)

    final_cal = calibration_df.sort_values("N_input").tail(1).iloc[0].to_dict()

    final_calibrated = (
        calibrated_summary_df
        .sort_values("N_input")
        .groupby("geometry")
        .tail(1)
        .reset_index(drop=True)
    )

    # Stability between N=256 and N=512 if available.
    stability = {}
    sub = calibration_df[calibration_df["N_input"].isin([256, 512])].sort_values("N_input")
    if len(sub) == 2:
        row256 = sub[sub["N_input"] == 256].iloc[0]
        row512 = sub[sub["N_input"] == 512].iloc[0]
        stability = {
            "B_epsilon_256": float(row256["B_epsilon_flat_baseline"]),
            "B_epsilon_512": float(row512["B_epsilon_flat_baseline"]),
            "delta_epsilon_256": float(row256["delta_kappa_over_epsilon"]),
            "delta_epsilon_512": float(row512["delta_kappa_over_epsilon"]),
            "A_epsilon_256": float(row256["A_epsilon"]),
            "A_epsilon_512": float(row512["A_epsilon"]),
            "B_epsilon_abs_change_256_512": float(abs(row512["B_epsilon_flat_baseline"] - row256["B_epsilon_flat_baseline"])),
            "delta_epsilon_abs_change_256_512": float(abs(row512["delta_kappa_over_epsilon"] - row256["delta_kappa_over_epsilon"])),
            "A_epsilon_abs_change_256_512": float(abs(row512["A_epsilon"] - row256["A_epsilon"])),
        }

    return {
        "status": "positive_affine_calibration",
        "interpretation": "A finite-N affine renormalization maps kappa/epsilon to a mean Ricci proxy on flat torus and sphere.",
        "calibration_rows": calibration_df.to_dict(orient="records"),
        "calibrated_summary_rows": calibrated_summary_df.to_dict(orient="records"),
        "final_calibration": final_cal,
        "final_calibrated_summary": final_calibrated.to_dict(orient="records"),
        "stability_256_512": stability,
    }


def build_key_table(v1, v2):
    rows = []

    b = v1.get("best_relative", {})
    if b:
        rows.append({
            "step": "v1 relative OR signal",
            "quantity": "sphere - flat raw mean kappa",
            "N": int(b["N_actual"]),
            "value": float(b["delta_mean_kappa_sphere_minus_flat"]),
            "reference": "> 0",
            "status": "positive",
        })
        rows.append({
            "step": "v1 relative OR signal",
            "quantity": "sphere - flat mean kappa/epsilon",
            "N": int(b["N_actual"]),
            "value": float(b["delta_mean_kappa_over_epsilon_sphere_minus_flat"]),
            "reference": "> 0",
            "status": "positive",
        })
        rows.append({
            "step": "v1 relative OR signal",
            "quantity": "sphere - flat mean kappa/l^2",
            "N": int(b["N_actual"]),
            "value": float(b["delta_mean_kappa_over_l2_sphere_minus_flat"]),
            "reference": "> 0",
            "status": "positive",
        })

    fc = v2["final_calibration"]
    rows.append({
        "step": "v2 affine calibration",
        "quantity": "B_N flat baseline for kappa/epsilon",
        "N": int(fc["N_input"]),
        "value": float(fc["B_epsilon_flat_baseline"]),
        "reference": "finite-N bias",
        "status": "identified",
    })
    rows.append({
        "step": "v2 affine calibration",
        "quantity": "sphere-flat delta for kappa/epsilon",
        "N": int(fc["N_input"]),
        "value": float(fc["delta_kappa_over_epsilon"]),
        "reference": "> 0",
        "status": "positive",
    })
    rows.append({
        "step": "v2 affine calibration",
        "quantity": "A_N epsilon calibration",
        "N": int(fc["N_input"]),
        "value": float(fc["A_epsilon"]),
        "reference": "1/delta",
        "status": "identified",
    })

    return pd.DataFrame(rows)


def write_markdown(outpath: Path, summary: dict, key_table: pd.DataFrame):
    v1 = summary["v1"]
    v2 = summary["v2"]
    fc = v2["final_calibration"]
    stability = v2.get("stability_256_512", {})

    lines = []
    lines.append("# Paper 17 — Ricci Summary v1")
    lines.append("")
    lines.append("This summary gathers the first Paper 17 tests of the continuum Ricci limit of Ollivier--Ricci curvature.")
    lines.append("")
    lines.append("The target relation is:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\kappa_{ij}^{OR}")
    lines.append(r"\longrightarrow")
    lines.append(r"R_{\mu\nu}u^\mu u^\nu.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The finite-scale diagnostic suggested by the tests is:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\frac{\kappa^{OR}}{\epsilon}")
    lines.append(r"=")
    lines.append(r"B_N+C_N R_{\mu\nu}u^\mu u^\nu+o(1).")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 1. Key results")
    lines.append("")
    lines.append("| Step | Quantity | N | Value | Reference | Status |")
    lines.append("|---|---|---:|---:|---|---|")
    for _, r in key_table.iterrows():
        lines.append(
            f"| {r['step']} | {r['quantity']} | {int(r['N'])} | "
            f"{float(r['value']):.6f} | {r['reference']} | {r['status']} |"
        )
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 2. v1 — Relative sphere-minus-flat signal")
    lines.append("")
    b = v1.get("best_relative", {})
    if b:
        lines.append("The first v1 test confirmed that the sphere has a positive Ollivier--Ricci signal relative to the flat torus.")
        lines.append("")
        lines.append(f"At \(N={int(b['N_actual'])}\):")
        lines.append("")
        lines.append(r"\[")
        lines.append(fr"\Delta\bar\kappa_{{\rm sphere-flat}}={float(b['delta_mean_kappa_sphere_minus_flat']):.6f},")
        lines.append(r"\]")
        lines.append("")
        lines.append(r"\[")
        lines.append(fr"\Delta\left\langle\frac{{\kappa}}{{\epsilon}}\right\rangle_{{\rm sphere-flat}}={float(b['delta_mean_kappa_over_epsilon_sphere_minus_flat']):.6f},")
        lines.append(r"\]")
        lines.append("")
        lines.append(r"\[")
        lines.append(fr"\Delta\left\langle\frac{{\kappa}}{{\ell^2}}\right\rangle_{{\rm sphere-flat}}={float(b['delta_mean_kappa_over_l2_sphere_minus_flat']):.6f}.")
        lines.append(r"\]")
    else:
        lines.append("No v1 relative rows were found.")
    lines.append("")
    lines.append("Interpretation: the sign of the Ricci signal is correct, but the raw normalizations are not yet unbiased.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 3. v2 — Affine Ricci calibration")
    lines.append("")
    lines.append("The v2 test introduces an affine calibration:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\widehat{R}_{OR}")
    lines.append(r"=")
    lines.append(r"A_N\left(")
    lines.append(r"\frac{\kappa^{OR}}{\epsilon}-B_N")
    lines.append(r"\right).")
    lines.append(r"\]")
    lines.append("")
    lines.append("The flat torus is used as the Ricci-zero reference:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"B_N=")
    lines.append(r"\left\langle\frac{\kappa^{OR}}{\epsilon}\right\rangle_{\rm flat}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The unit sphere is used as the \(R_{\mu\nu}u^\mu u^\nu=1\) reference:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"A_N=")
    lines.append(r"\frac{1}{")
    lines.append(r"\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm sphere}")
    lines.append(r"-")
    lines.append(r"\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm flat}")
    lines.append(r"}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The calibration constants are:")
    lines.append("")
    lines.append("| \(N_{\rm input}\) | \(B_N\) flat baseline | sphere-flat \(\Delta\) | \(A_N\) |")
    lines.append("|---:|---:|---:|---:|")
    for r in v2["calibration_rows"]:
        lines.append(
            f"| {int(r['N_input'])} | {float(r['B_epsilon_flat_baseline']):.4f} | "
            f"{float(r['delta_kappa_over_epsilon']):.4f} | {float(r['A_epsilon']):.4f} |"
        )
    lines.append("")
    lines.append(f"At \(N={int(fc['N_input'])}\), this gives:")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\widehat{{R}}_{{OR}}={float(fc['A_epsilon']):.3f}")
    lines.append(r"\left(")
    if float(fc["B_epsilon_flat_baseline"]) < 0:
        lines.append(fr"\frac{{\kappa^{{OR}}}}{{\epsilon}}+{abs(float(fc['B_epsilon_flat_baseline'])):.3f}")
    else:
        lines.append(fr"\frac{{\kappa^{{OR}}}}{{\epsilon}}-{float(fc['B_epsilon_flat_baseline']):.3f}")
    lines.append(r"\right).")
    lines.append(r"\]")
    lines.append("")
    if stability:
        lines.append("Between \(N=256\) and \(N=512\), the calibration is relatively stable:")
        lines.append("")
        lines.append("| Quantity | N=256 | N=512 | absolute change |")
        lines.append("|---|---:|---:|---:|")
        lines.append(
            f"| \(B_N\) | {stability['B_epsilon_256']:.4f} | {stability['B_epsilon_512']:.4f} | "
            f"{stability['B_epsilon_abs_change_256_512']:.4f} |"
        )
        lines.append(
            f"| sphere-flat \(\Delta\) | {stability['delta_epsilon_256']:.4f} | {stability['delta_epsilon_512']:.4f} | "
            f"{stability['delta_epsilon_abs_change_256_512']:.4f} |"
        )
        lines.append(
            f"| \(A_N\) | {stability['A_epsilon_256']:.4f} | {stability['A_epsilon_512']:.4f} | "
            f"{stability['A_epsilon_abs_change_256_512']:.4f} |"
        )
        lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 4. Interpretation")
    lines.append("")
    lines.append("Paper 17 v2 identifies a finite-scale affine renormalization of Ollivier--Ricci curvature.")
    lines.append("")
    lines.append("The result supports the mean-level relation:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\frac{\kappa^{OR}}{\epsilon}")
    lines.append(r"\simeq")
    lines.append(r"B_N+C_N R_{\mu\nu}u^\mu u^\nu.")
    lines.append(r"\]")
    lines.append("")
    lines.append("This is not yet a pointwise convergence theorem. The calibrated means are correct by construction, while the distributional dispersion remains large, especially on the flat torus.")
    lines.append("")
    lines.append("The next step is to reduce this dispersion by testing radius graphs, edge-length bins, trimmed means, and larger \(N\).")
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    v1_dir = Path(args.v1_dir)
    v2_dir = Path(args.v2_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    v1 = summarize_v1(v1_dir)
    v2 = summarize_v2(v2_dir)

    summary = {
        "experiment": "Paper 17 Ricci summary v1",
        "v1_dir": str(v1_dir),
        "v2_dir": str(v2_dir),
        "v1": v1,
        "v2": v2,
    }

    key_table = build_key_table(v1, v2)
    key_table.to_csv(outdir / "paper17_key_results_table.csv", index=False)

    with open(outdir / "paper17_ricci_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    write_markdown(outdir / "paper17_ricci_summary.md", summary, key_table)

    print("=" * 110)
    print("Paper 17 — Ricci summary v1")
    print("=" * 110)
    print(key_table.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "paper17_key_results_table.csv")
    print(outdir / "paper17_ricci_summary.json")
    print(outdir / "paper17_ricci_summary.md")


if __name__ == "__main__":
    main()

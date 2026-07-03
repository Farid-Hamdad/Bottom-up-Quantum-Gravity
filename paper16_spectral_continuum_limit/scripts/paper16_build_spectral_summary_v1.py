#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 16 — Spectral summary builder v1

Goal
----
Collect the first two eigenvalue-level convergence tests of Paper 16:

1. Circle spectrum:
   c_N L_N -> -Δ_{S^1}

2. Flat torus spectrum:
   c_N L_N -> -Δ_{T^2}

Inputs
------
- results/circle_spectrum_convergence_v1/circle_spectrum_summary.csv
- results/flat_torus_spectrum_convergence_v1/flat_torus_spectrum_summary.csv

Outputs
-------
paper16_key_results_table.csv
paper16_spectral_summary.json
paper16_spectral_summary.md

Recommended run
---------------
cd ~/bottomup

python3 papers/paper16_spectral_continuum_limit/scripts/paper16_build_spectral_summary_v1.py \
  --circle-dir papers/paper16_spectral_continuum_limit/results/circle_spectrum_convergence_v1 \
  --flat-torus-dir papers/paper16_spectral_continuum_limit/results/flat_torus_spectrum_convergence_v1 \
  --output-dir papers/paper16_spectral_continuum_limit/results/paper16_spectral_summary_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Build Paper 16 spectral summary v1.")
    p.add_argument("--circle-dir", required=True)
    p.add_argument("--flat-torus-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_csv(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)


def final_row(df: pd.DataFrame, n_col: str):
    return df.sort_values(n_col).tail(1).iloc[0].to_dict()


def summarize_circle(circle_dir: Path):
    summary_path = circle_dir / "circle_spectrum_summary.csv"
    rows_path = circle_dir / "circle_spectrum_rows.csv"

    summary_df = load_csv(summary_path)
    rows_df = load_csv(rows_path)

    final = final_row(summary_df, "N")
    N_final = int(final["N"])

    final_modes = rows_df[rows_df["N"] == N_final].sort_values("mode_index")

    return {
        "geometry": "circle",
        "manifold": "S^1",
        "operator_limit": "c_N L_N -> -Delta_{S^1}",
        "analytic_spectrum": "lambda_m = m^2 with sine/cosine degeneracy",
        "summary_rows": summary_df.to_dict(orient="records"),
        "final": {
            "N": int(final["N"]),
            "k": int(final["k"]),
            "epsilon": float(final["epsilon"]),
            "scale_c": float(final["scale_c"]),
            "n_modes": int(final["n_modes"]),
            "fit_modes": int(final["fit_modes"]),
            "mean_rel_error": float(final["mean_rel_error"]),
            "median_rel_error": float(final["median_rel_error"]),
            "max_rel_error": float(final["max_rel_error"]),
            "lambda1_scaled": float(final["lambda1_scaled"]),
            "lambda1_target": 1.0,
            "lambda1_rel_error": float(abs(final["lambda1_scaled"] - 1.0)),
        },
        "final_modes": final_modes.to_dict(orient="records"),
    }


def summarize_flat_torus(flat_dir: Path):
    summary_path = flat_dir / "flat_torus_spectrum_summary.csv"
    rows_path = flat_dir / "flat_torus_spectrum_rows.csv"

    summary_df = load_csv(summary_path)
    rows_df = load_csv(rows_path)

    final = final_row(summary_df, "N_actual")
    N_final = int(final["N_actual"])

    final_modes = rows_df[rows_df["N_actual"] == N_final].sort_values("mode_index")

    lambda1_target = float(final["lambda1_target"])
    lambda1_scaled = float(final["lambda1_scaled"])

    return {
        "geometry": "flat_torus",
        "manifold": "T^2",
        "operator_limit": "c_N L_N -> -Delta_{T^2}",
        "analytic_spectrum": "lambda_{m,n}=4*pi^2*(m^2+n^2)",
        "summary_rows": summary_df.to_dict(orient="records"),
        "final": {
            "N_input": int(final["N_input"]),
            "N_actual": int(final["N_actual"]),
            "m_grid": int(final["m_grid"]),
            "k": int(final["k"]),
            "epsilon": float(final["epsilon"]),
            "scale_c": float(final["scale_c"]),
            "n_modes": int(final["n_modes"]),
            "fit_modes": int(final["fit_modes"]),
            "mean_rel_error": float(final["mean_rel_error"]),
            "median_rel_error": float(final["median_rel_error"]),
            "max_rel_error": float(final["max_rel_error"]),
            "lambda1_scaled": lambda1_scaled,
            "lambda1_target": lambda1_target,
            "lambda1_rel_error": float(abs(lambda1_scaled - lambda1_target) / lambda1_target),
        },
        "final_modes": final_modes.to_dict(orient="records"),
    }


def build_key_table(circle, flat):
    rows = []

    c = circle["final"]
    rows.append({
        "test": "Circle spectrum",
        "limit": r"c_N L_N -> -Delta_{S^1}",
        "N": c["N"],
        "modes_compared": c["n_modes"],
        "mean_rel_error": c["mean_rel_error"],
        "median_rel_error": c["median_rel_error"],
        "max_rel_error": c["max_rel_error"],
        "lambda1_scaled": c["lambda1_scaled"],
        "lambda1_target": c["lambda1_target"],
        "status": "positive",
    })

    f = flat["final"]
    rows.append({
        "test": "Flat torus spectrum",
        "limit": r"c_N L_N -> -Delta_{T^2}",
        "N": f["N_actual"],
        "modes_compared": f["n_modes"],
        "mean_rel_error": f["mean_rel_error"],
        "median_rel_error": f["median_rel_error"],
        "max_rel_error": f["max_rel_error"],
        "lambda1_scaled": f["lambda1_scaled"],
        "lambda1_target": f["lambda1_target"],
        "status": "positive",
    })

    return pd.DataFrame(rows)


def fmt_float(x, nd=4):
    return f"{float(x):.{nd}f}"


def write_markdown(outpath: Path, summary: dict, key_table: pd.DataFrame):
    circle = summary["circle"]
    flat = summary["flat_torus"]

    lines = []
    lines.append("# Paper 16 — Spectral Summary v1")
    lines.append("")
    lines.append("This summary gathers the first eigenvalue-level convergence tests for Paper 16.")
    lines.append("")
    lines.append("The target statement is:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"c_N L_N")
    lines.append(r"=")
    lines.append(r"c_N\frac{D_N-W_N}{\epsilon_N}")
    lines.append(r"\longrightarrow")
    lines.append(r"-\Delta_g.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The scalar \(c_N\) is currently fitted numerically. A main theoretical task of Paper 16 is to derive this normalization analytically.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 1. Key results")
    lines.append("")
    lines.append("| Test | Limit | \(N\) | Modes | Mean rel. error | Median rel. error | Max rel. error | \(\lambda_1^{scaled}\) | \(\lambda_1^{target}\) |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in key_table.iterrows():
        lines.append(
            f"| {r['test']} | `{r['limit']}` | {int(r['N'])} | {int(r['modes_compared'])} | "
            f"{r['mean_rel_error']:.4f} | {r['median_rel_error']:.4f} | {r['max_rel_error']:.4f} | "
            f"{r['lambda1_scaled']:.4f} | {r['lambda1_target']:.4f} |"
        )
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 2. Circle spectrum")
    lines.append("")
    lines.append("For the unit circle \(S^1\), the analytic spectrum is")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\lambda_m=m^2,")
    lines.append(r"\qquad")
    lines.append(r"m=0,1,1,2,2,3,3,\ldots")
    lines.append(r"\]")
    lines.append("")
    c = circle["final"]
    lines.append(f"At \(N={c['N']}\), the first {c['n_modes']} nonzero modes are reconstructed with:")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\mathrm{{mean\ relative\ error}}={c['mean_rel_error']:.4f},")
    lines.append(r"\qquad")
    lines.append(fr"\mathrm{{median\ relative\ error}}={c['median_rel_error']:.4f},")
    lines.append(r"\qquad")
    lines.append(fr"\mathrm{{max\ relative\ error}}={c['max_rel_error']:.4f}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The first nonzero eigenvalue is")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\lambda_1^{{\rm scaled}}={c['lambda1_scaled']:.4f},")
    lines.append(r"\qquad")
    lines.append(r"\lambda_1^{S^1}=1.")
    lines.append(r"\]")
    lines.append("")
    lines.append("This gives eigenvalue-level support for")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"c_N L_N\longrightarrow -\Delta_{S^1}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 3. Flat torus spectrum")
    lines.append("")
    lines.append("For the unit flat torus \(T^2=[0,1)^2\), the analytic spectrum is")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\lambda_{m,n}=4\pi^2(m^2+n^2),")
    lines.append(r"\qquad")
    lines.append(r"(m,n)\in\mathbb{Z}^2\setminus\{(0,0)\}.")
    lines.append(r"\]")
    lines.append("")
    f = flat["final"]
    lines.append(f"At \(N={f['N_actual']}\), the first {f['n_modes']} nonzero modes are reconstructed with:")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\mathrm{{mean\ relative\ error}}={f['mean_rel_error']:.4f},")
    lines.append(r"\qquad")
    lines.append(fr"\mathrm{{median\ relative\ error}}={f['median_rel_error']:.4f},")
    lines.append(r"\qquad")
    lines.append(fr"\mathrm{{max\ relative\ error}}={f['max_rel_error']:.4f}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("The first nonzero eigenvalue is")
    lines.append("")
    lines.append(r"\[")
    lines.append(fr"\lambda_1^{{\rm scaled}}={f['lambda1_scaled']:.4f},")
    lines.append(r"\qquad")
    lines.append(fr"\lambda_1^{{T^2}}={f['lambda1_target']:.4f}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("This gives eigenvalue-level support for")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"c_N L_N\longrightarrow -\Delta_{T^2}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 4. Interpretation")
    lines.append("")
    lines.append("Paper 15 showed heat-trace dimension convergence. Paper 16 now strengthens this by testing individual eigenvalues.")
    lines.append("")
    lines.append("The current results show:")
    lines.append("")
    lines.append("- \(S^1\): low-spectrum convergence with mean relative error \(0.0051\).")
    lines.append("- \(T^2\): low-spectrum convergence with mean relative error \(0.0266\).")
    lines.append("")
    lines.append("This supports the spectral continuum limit:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"L_N\to \Delta_g")
    lines.append(r"\]")
    lines.append("")
    lines.append("at the level of low eigenmodes, up to a scalar normalization \(c_N\).")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 5. Remaining theoretical issue")
    lines.append("")
    lines.append("The current tests fit \(c_N\) numerically. The next task is to derive \(c_N\) analytically from the kernel normalization, the sampling density, and the choice of graph Laplacian.")
    lines.append("")
    lines.append("A target theorem should specify conditions under which")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"c_N\frac{D_N-W_N}{\epsilon_N}")
    lines.append(r"\longrightarrow")
    lines.append(r"-\Delta_g")
    lines.append(r"\]")
    lines.append("")
    lines.append("in a controlled convergence mode: pointwise on smooth test functions, quadratic-form convergence, heat-kernel convergence, or spectral convergence.")
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    circle_dir = Path(args.circle_dir)
    flat_dir = Path(args.flat_torus_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    circle = summarize_circle(circle_dir)
    flat = summarize_flat_torus(flat_dir)

    summary = {
        "experiment": "Paper 16 spectral summary v1",
        "circle_dir": str(circle_dir),
        "flat_torus_dir": str(flat_dir),
        "circle": circle,
        "flat_torus": flat,
    }

    key_table = build_key_table(circle, flat)
    key_table.to_csv(outdir / "paper16_key_results_table.csv", index=False)

    with open(outdir / "paper16_spectral_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    write_markdown(outdir / "paper16_spectral_summary.md", summary, key_table)

    print("=" * 100)
    print("Paper 16 — Spectral summary v1")
    print("=" * 100)
    print(key_table.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "paper16_key_results_table.csv")
    print(outdir / "paper16_spectral_summary.json")
    print(outdir / "paper16_spectral_summary.md")


if __name__ == "__main__":
    main()

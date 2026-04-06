#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_bup_v3_results.py

Analyse complète des résultats de la v3 :
- lit automatiquement aggregate_rows.csv
- calcule les stats par type d'état
- calcule les stats par anisotropie eta_micro
- calcule la fraction de cas avec amplification > 1
- écrit des CSV de synthèse
- génère des figures PNG

Usage recommandé depuis ~/bottomup/jauge :

python analyze_bup_v3_results.py \
  --csv results_bup_modular_referee_v3/aggregate_rows.csv \
  --output-dir results_bup_modular_referee_v3_analysis

Si tu veux, tu peux aussi juste lancer :
python analyze_bup_v3_results.py

Dans ce cas, le script cherchera automatiquement :
- results_bup_modular_referee_v3/aggregate_rows.csv
- ./aggregate_rows.csv
"""

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def find_default_csv() -> Path:
    candidates = [
        Path("results_bup_modular_referee_v3/aggregate_rows.csv"),
        Path("./aggregate_rows.csv"),
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise FileNotFoundError(
        "Impossible de trouver aggregate_rows.csv automatiquement.\n"
        "Essaye avec : --csv results_bup_modular_referee_v3/aggregate_rows.csv"
    )


def make_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_mean(x):
    x = pd.Series(x).replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.mean()) if len(x) else np.nan


def safe_std(x):
    x = pd.Series(x).replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.std(ddof=1)) if len(x) > 1 else np.nan


def check_required_columns(df: pd.DataFrame):
    required = {
        "state_kind",
        "eta_micro",
        "amplification_baseline",
        "fit_r2_baseline",
        "locality_score",
        "delta_far",
        "delta_near",
        "N",
        "A_size",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Colonnes manquantes dans le CSV : {sorted(missing)}\n"
            f"Colonnes disponibles : {list(df.columns)}"
        )


def main():
    parser = argparse.ArgumentParser(description="Analyse des résultats BUP modular referee v3")
    parser.add_argument(
        "--csv",
        type=str,
        default=None,
        help="Chemin vers aggregate_rows.csv "
             "(par défaut le script cherche results_bup_modular_referee_v3/aggregate_rows.csv)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results_bup_modular_referee_v3_analysis",
        help="Dossier de sortie pour les tableaux et figures",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv).resolve() if args.csv else find_default_csv()
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV introuvable : {csv_path}")

    outdir = Path(args.output_dir).resolve()
    make_output_dir(outdir)

    df = pd.read_csv(csv_path)
    check_required_columns(df)

    # Nettoyage léger
    df = df.copy()
    df["amplification_baseline_clean"] = pd.to_numeric(df["amplification_baseline"], errors="coerce")
    df["fit_r2_baseline_clean"] = pd.to_numeric(df["fit_r2_baseline"], errors="coerce")
    df["locality_score_clean"] = pd.to_numeric(df["locality_score"], errors="coerce")
    df["delta_far_clean"] = pd.to_numeric(df["delta_far"], errors="coerce")
    df["delta_near_clean"] = pd.to_numeric(df["delta_near"], errors="coerce")
    df["eta_micro"] = pd.to_numeric(df["eta_micro"], errors="coerce")
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df["A_size"] = pd.to_numeric(df["A_size"], errors="coerce")

    # 1) Stats par type d'état
    rows = []
    for state, sub in df.groupby("state_kind"):
        amp = sub["amplification_baseline_clean"].replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({
            "state_kind": state,
            "n_cases": int(len(sub)),
            "mean_amplification": safe_mean(amp),
            "std_amplification": safe_std(amp),
            "mean_fit_r2": safe_mean(sub["fit_r2_baseline_clean"]),
            "mean_locality_score": safe_mean(sub["locality_score_clean"]),
            "mean_delta_near": safe_mean(sub["delta_near_clean"]),
            "mean_delta_far": safe_mean(sub["delta_far_clean"]),
            "fraction_amplification_gt_1": float((amp > 1.0).mean()) if len(amp) else np.nan,
        })
    by_state = pd.DataFrame(rows).sort_values("state_kind")
    by_state.to_csv(outdir / "summary_by_state_kind.csv", index=False)

    # 2) Stats par type d'état et anisotropie
    rows = []
    for (state, eta), sub in df.groupby(["state_kind", "eta_micro"]):
        amp = sub["amplification_baseline_clean"].replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({
            "state_kind": state,
            "eta_micro": float(eta),
            "n_cases": int(len(sub)),
            "mean_amplification": safe_mean(amp),
            "std_amplification": safe_std(amp),
            "mean_fit_r2": safe_mean(sub["fit_r2_baseline_clean"]),
            "mean_locality_score": safe_mean(sub["locality_score_clean"]),
            "fraction_amplification_gt_1": float((amp > 1.0).mean()) if len(amp) else np.nan,
        })
    by_state_eta = pd.DataFrame(rows).sort_values(["state_kind", "eta_micro"])
    by_state_eta.to_csv(outdir / "summary_by_state_and_eta.csv", index=False)

    # 3) Stats par N et A_size
    rows = []
    for (state, N, A_size), sub in df.groupby(["state_kind", "N", "A_size"]):
        amp = sub["amplification_baseline_clean"].replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({
            "state_kind": state,
            "N": int(N),
            "A_size": int(A_size),
            "n_cases": int(len(sub)),
            "mean_amplification": safe_mean(amp),
            "mean_fit_r2": safe_mean(sub["fit_r2_baseline_clean"]),
            "mean_locality_score": safe_mean(sub["locality_score_clean"]),
        })
    by_state_scale = pd.DataFrame(rows).sort_values(["state_kind", "N", "A_size"])
    by_state_scale.to_csv(outdir / "summary_by_state_N_A.csv", index=False)

    # Console summary
    print("==========================================")
    print("BUP V3 ANALYSIS")
    print("==========================================")
    print(f"CSV lu           : {csv_path}")
    print(f"Dossier de sortie: {outdir}")
    print()
    print("=== Stats par type d'état ===")
    print(by_state.to_string(index=False))
    print()
    print("=== Stats par type d'état et eta_micro ===")
    print(by_state_eta.to_string(index=False))
    print()

    # -------------------------
    # PLOTS
    # -------------------------

    # Plot 1: amplification vs eta_micro par type d'état
    plt.figure(figsize=(7, 5))
    for state, sub in by_state_eta.groupby("state_kind"):
        plt.plot(sub["eta_micro"], sub["mean_amplification"], marker="o", label=state)
    plt.axhline(1.0, linestyle="--")
    plt.xlabel("eta_micro = Jz / Jxy")
    plt.ylabel("mean amplification")
    plt.title("Amplification vs microscopic anisotropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "amplification_vs_eta_micro_by_state.png", dpi=180)
    plt.close()

    # Plot 2: locality score vs eta_micro
    plt.figure(figsize=(7, 5))
    for state, sub in by_state_eta.groupby("state_kind"):
        plt.plot(sub["eta_micro"], sub["mean_locality_score"], marker="o", label=state)
    plt.xlabel("eta_micro = Jz / Jxy")
    plt.ylabel("mean locality score")
    plt.title("Locality score vs microscopic anisotropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "locality_vs_eta_micro_by_state.png", dpi=180)
    plt.close()

    # Plot 3: fit R² vs eta_micro
    plt.figure(figsize=(7, 5))
    for state, sub in by_state_eta.groupby("state_kind"):
        plt.plot(sub["eta_micro"], sub["mean_fit_r2"], marker="o", label=state)
    plt.xlabel("eta_micro = Jz / Jxy")
    plt.ylabel("mean fit R²")
    plt.title("Fit quality vs microscopic anisotropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fit_r2_vs_eta_micro_by_state.png", dpi=180)
    plt.close()

    # Plot 4: mean amplification by state kind
    plt.figure(figsize=(7, 5))
    plt.bar(by_state["state_kind"], by_state["mean_amplification"])
    plt.axhline(1.0, linestyle="--")
    plt.ylabel("mean amplification")
    plt.title("Mean amplification by state kind")
    plt.tight_layout()
    plt.savefig(outdir / "amplification_by_state_kind.png", dpi=180)
    plt.close()

    # Plot 5: fraction amplification > 1 by state kind
    plt.figure(figsize=(7, 5))
    plt.bar(by_state["state_kind"], by_state["fraction_amplification_gt_1"])
    plt.ylabel("fraction amplification > 1")
    plt.title("How often amplification exceeds 1")
    plt.tight_layout()
    plt.savefig(outdir / "fraction_amplification_gt_1_by_state.png", dpi=180)
    plt.close()

    # README / report
    with open(outdir / "REPORT.md", "w", encoding="utf-8") as f:
        f.write("# BUP v3 analysis report\n\n")
        f.write(f"- Input CSV: `{csv_path}`\n")
        f.write(f"- Number of rows: **{len(df)}**\n\n")
        f.write("## Main outputs\n\n")
        f.write("- `summary_by_state_kind.csv`\n")
        f.write("- `summary_by_state_and_eta.csv`\n")
        f.write("- `summary_by_state_N_A.csv`\n")
        f.write("- `amplification_vs_eta_micro_by_state.png`\n")
        f.write("- `locality_vs_eta_micro_by_state.png`\n")
        f.write("- `fit_r2_vs_eta_micro_by_state.png`\n")
        f.write("- `amplification_by_state_kind.png`\n")
        f.write("- `fraction_amplification_gt_1_by_state.png`\n\n")
        f.write("## Quick interpretation guide\n\n")
        f.write("- amplification > 1 : anisotropy enhancement\n")
        f.write("- amplification ~ 1 : anisotropy inheritance only\n")
        f.write("- amplification < 1 : anisotropy smoothing\n")
        f.write("- locality_score > 0 : stronger response on local-support terms\n")
        f.write("- fit R² close to 1 : effective 1+2-body description is accurate\n")

    print(f"Fichiers écrits dans : {outdir}")
    print("Terminé.")


if __name__ == "__main__":
    main()

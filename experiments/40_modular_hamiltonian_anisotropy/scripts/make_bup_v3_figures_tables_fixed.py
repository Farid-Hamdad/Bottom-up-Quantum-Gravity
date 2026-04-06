#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_default_csv() -> Path:
    candidates = [
        Path("results_bup_modular_referee_v3_scaling/aggregate_rows.csv"),
        Path("results_bup_modular_referee_v3/aggregate_rows.csv"),
        Path("aggregate_rows.csv"),
    ]
    for p in candidates:
        if p.exists():
            return p.resolve()
    raise FileNotFoundError(
        "aggregate_rows.csv introuvable. Donne-le explicitement avec --csv "
        "ou place-toi dans ~/bottomup/jauge."
    )


def safe_series(s):
    return pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()


def mean_std_df(df, group_cols, value_cols):
    rows = []
    for keys, sub in df.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: val for col, val in zip(group_cols, keys)}
        row["n_cases"] = int(len(sub))
        for v in value_cols:
            x = safe_series(sub[v])
            row[f"{v}_mean"] = float(x.mean()) if len(x) else np.nan
            row[f"{v}_std"] = float(x.std(ddof=1)) if len(x) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot_amplification_by_state(df, outpath: Path):
    stat = mean_std_df(df, ["state_kind"], ["amplification_baseline"]).sort_values("state_kind")
    x = np.arange(len(stat))
    y = stat["amplification_baseline_mean"].to_numpy(dtype=float)
    yerr = stat["amplification_baseline_std"].fillna(0).to_numpy(dtype=float)

    plt.figure(figsize=(7, 5))
    plt.bar(x, y, yerr=yerr, capsize=5)
    plt.axhline(1.0, linestyle="--")
    plt.xticks(x, stat["state_kind"])
    plt.ylabel("Amplification")
    plt.title("Anisotropy amplification by state class")
    plt.tight_layout()
    plt.savefig(outpath, dpi=220)
    plt.close()
    return stat


def plot_amplification_vs_eta(df, outpath: Path):
    stat = mean_std_df(df, ["state_kind", "eta_micro"], ["amplification_baseline"]).sort_values(["state_kind", "eta_micro"])
    plt.figure(figsize=(7, 5))
    for state, sub in stat.groupby("state_kind"):
        plt.errorbar(
            sub["eta_micro"].to_numpy(dtype=float),
            sub["amplification_baseline_mean"].to_numpy(dtype=float),
            yerr=sub["amplification_baseline_std"].fillna(0).to_numpy(dtype=float),
            marker="o",
            capsize=4,
            label=str(state),
        )
    plt.axhline(1.0, linestyle="--")
    plt.xlabel(r"$\eta_{\rm micro}=J_z/J_{xy}$")
    plt.ylabel(r"Amplification $A=\eta_{\rm eff}/\eta_{\rm micro}$")
    plt.title("Amplification vs microscopic anisotropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=220)
    plt.close()
    return stat


def plot_locality_vs_eta(df, outpath: Path):
    stat = mean_std_df(df, ["state_kind", "eta_micro"], ["locality_score"]).sort_values(["state_kind", "eta_micro"])
    plt.figure(figsize=(7, 5))
    for state, sub in stat.groupby("state_kind"):
        plt.errorbar(
            sub["eta_micro"].to_numpy(dtype=float),
            sub["locality_score_mean"].to_numpy(dtype=float),
            yerr=sub["locality_score_std"].fillna(0).to_numpy(dtype=float),
            marker="o",
            capsize=4,
            label=str(state),
        )
    plt.xlabel(r"$\eta_{\rm micro}=J_z/J_{xy}$")
    plt.ylabel("Locality score")
    plt.title("Quasi-local response vs microscopic anisotropy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=220)
    plt.close()
    return stat


def plot_fit_r2_by_state(df, outpath: Path):
    stat = mean_std_df(df, ["state_kind"], ["fit_r2_baseline"]).sort_values("state_kind")
    x = np.arange(len(stat))
    y = stat["fit_r2_baseline_mean"].to_numpy(dtype=float)
    yerr = stat["fit_r2_baseline_std"].fillna(0).to_numpy(dtype=float)

    plt.figure(figsize=(7, 5))
    plt.bar(x, y, yerr=yerr, capsize=5)
    plt.xticks(x, stat["state_kind"])
    plt.ylabel(r"Fit $R^2$")
    plt.title("Quality of the 1+2-body modular fit")
    plt.tight_layout()
    plt.savefig(outpath, dpi=220)
    plt.close()
    return stat


def make_typical_coeff_table(df, outdir: Path):
    sub = df[df["state_kind"] == "ground"].copy()
    if sub.empty:
        sub = df.copy()
    sub["eta_dist"] = (pd.to_numeric(sub["eta_micro"], errors="coerce") - 1.5).abs()
    sub = sub.sort_values(["eta_dist", "N", "A_size"], ascending=[True, False, False])
    row = sub.iloc[0]

    table = pd.DataFrame([{
        "state_kind": row["state_kind"],
        "N": int(row["N"]),
        "A_size": int(row["A_size"]),
        "eta_micro": float(row["eta_micro"]),
        "mean_abs_J_xx": float(row["mean_abs_xx_baseline"]),
        "mean_abs_J_yy": float(row["mean_abs_yy_baseline"]),
        "mean_abs_J_zz": float(row["mean_abs_zz_baseline"]),
        "mean_abs_J_cross": float(row["mean_abs_cross_baseline"]),
        "fit_r2": float(row["fit_r2_baseline"]),
        "amplification": float(row["amplification_baseline"]),
        "xxz_regime": row["xxz_regime_baseline"],
    }])
    table.to_csv(outdir / "table_typical_coefficients.csv", index=False)
    return table


def dataframe_to_plain_text(df: pd.DataFrame) -> str:
    return df.to_string(index=False)


def main():
    parser = argparse.ArgumentParser(description="Generate figures and tables for the BUP v3 paper")
    parser.add_argument("--csv", type=str, default=None, help="Path to aggregate_rows.csv")
    parser.add_argument("--output-dir", type=str, default="results_bup_v3_paper_assets", help="Output directory")
    args = parser.parse_args()

    csv_path = Path(args.csv).resolve() if args.csv else find_default_csv()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)

    required_cols = {
        "state_kind", "eta_micro", "amplification_baseline", "fit_r2_baseline",
        "locality_score", "N", "A_size",
        "mean_abs_xx_baseline", "mean_abs_yy_baseline", "mean_abs_zz_baseline",
        "mean_abs_cross_baseline", "xxz_regime_baseline"
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Colonnes manquantes : {sorted(missing)}")

    table_state = mean_std_df(
        df,
        ["state_kind"],
        ["amplification_baseline", "fit_r2_baseline", "locality_score"]
    ).sort_values("state_kind")
    table_state.to_csv(outdir / "table_state_summary.csv", index=False)

    table_state_eta = mean_std_df(
        df,
        ["state_kind", "eta_micro"],
        ["amplification_baseline", "fit_r2_baseline", "locality_score"]
    ).sort_values(["state_kind", "eta_micro"])
    table_state_eta.to_csv(outdir / "table_state_eta_summary.csv", index=False)

    table_state_scale = mean_std_df(
        df,
        ["state_kind", "N", "A_size"],
        ["amplification_baseline", "fit_r2_baseline", "locality_score"]
    ).sort_values(["state_kind", "N", "A_size"])
    table_state_scale.to_csv(outdir / "table_state_scaling_summary.csv", index=False)

    plot_amplification_by_state(df, outdir / "fig_amplification_by_state.png")
    plot_amplification_vs_eta(df, outdir / "fig_amplification_vs_eta.png")
    plot_locality_vs_eta(df, outdir / "fig_locality_vs_eta.png")
    plot_fit_r2_by_state(df, outdir / "fig_fit_r2_by_state.png")

    typical = make_typical_coeff_table(df, outdir)

    with open(outdir / "REPORT.md", "w", encoding="utf-8") as f:
        f.write("# BUP v3 paper assets\n\n")
        f.write(f"- Input CSV: `{csv_path}`\n")
        f.write(f"- Number of rows: **{len(df)}**\n\n")
        f.write("## Figures\n\n")
        f.write("- `fig_amplification_by_state.png`\n")
        f.write("- `fig_amplification_vs_eta.png`\n")
        f.write("- `fig_locality_vs_eta.png`\n")
        f.write("- `fig_fit_r2_by_state.png`\n\n")
        f.write("## Tables\n\n")
        f.write("- `table_state_summary.csv`\n")
        f.write("- `table_state_eta_summary.csv`\n")
        f.write("- `table_state_scaling_summary.csv`\n")
        f.write("- `table_typical_coefficients.csv`\n\n")
        f.write("## Typical coefficient example\n\n")
        f.write(dataframe_to_plain_text(typical))
        f.write("\n")

    print("======================================")
    print("BUP V3 FIGURES + TABLES GENERATED")
    print("======================================")
    print(f"Input CSV : {csv_path}")
    print(f"Output dir: {outdir}")
    for name in [
        "fig_amplification_by_state.png",
        "fig_amplification_vs_eta.png",
        "fig_locality_vs_eta.png",
        "fig_fit_r2_by_state.png",
        "table_state_summary.csv",
        "table_state_eta_summary.csv",
        "table_state_scaling_summary.csv",
        "table_typical_coefficients.csv",
        "REPORT.md",
    ]:
        print(" -", outdir / name)
    print("======================================")


if __name__ == "__main__":
    main()

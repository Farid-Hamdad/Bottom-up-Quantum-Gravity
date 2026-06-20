#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP gravitational waves — true edge Hessian v10 parameter scan.

This script launches bup_gw_true_edge_hessian_v10.py over a parameter grid and
collects the spectral diagnostics:
  - edge-Hessian dispersion: c_edge, m_eff^2, R2
  - modal response: corr(|J_n|, q_peak) for plus/cross
  - spectral concentration: low-mode energy, k_weighted, dominant mode/fraction

It intentionally treats packet/ring speed as secondary because high-k edge modes
need not form a clean radial packet.

Example:
python3 scripts/scan_true_edge_hessian_v10.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/bup_gw_true_edge_hessian_v10.py \
  --output-dir papers/paper22_bup_gravitational_waves/results/scan_true_edge_hessian_v10
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def safe_float(x):
    if x is None:
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan


def run_one(cmd: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    return int(proc.returncode)


def get_nested(d: dict, path: list[str], default=np.nan):
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def read_summary(run_dir: Path) -> dict | None:
    p = run_dir / "results" / "summary_true_edge_hessian_v10.json"
    if not p.exists():
        return None
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def scatter_plot(df: pd.DataFrame, x: str, y: str, outpath: Path, title: str, xlabel: str, ylabel: str):
    ok = np.isfinite(df[x].to_numpy(dtype=float)) & np.isfinite(df[y].to_numpy(dtype=float))
    plt.figure(figsize=(6.2, 4.6))
    if ok.sum() > 0:
        plt.scatter(df.loc[ok, x], df.loc[ok, y], s=50)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def heatmap_pivot(df: pd.DataFrame, row: str, col: str, val: str, outpath: Path, title: str):
    tmp = df.groupby([row, col])[val].mean().reset_index()
    piv = tmp.pivot(index=row, columns=col, values=val)
    plt.figure(figsize=(7.0, 4.8))
    im = plt.imshow(piv.to_numpy(dtype=float), aspect="auto", origin="lower")
    plt.colorbar(im, label=val)
    plt.xticks(np.arange(len(piv.columns)), [str(x) for x in piv.columns], rotation=45)
    plt.yticks(np.arange(len(piv.index)), [str(x) for x in piv.index])
    plt.xlabel(col)
    plt.ylabel(row)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--v10-script", required=True, help="Path to bup_gw_true_edge_hessian_v10.py")
    p.add_argument("--output-dir", required=True)

    # Scan parameters. Defaults centered around the successful v10_highk run.
    p.add_argument("--eta-edge-list", default="0.75,1.0,1.25")
    p.add_argument("--source-width-list", default="0.12,0.16,0.22")
    p.add_argument("--pulse-sigma-list", default="0.25,0.35,0.50")

    # Fixed physical/numerical parameters.
    p.add_argument("--N-side", type=int, default=11)
    p.add_argument("--extent", type=float, default=1.0)
    p.add_argument("--ell", type=float, default=0.22)
    p.add_argument("--cutoff", type=float, default=0.45)
    p.add_argument("--mass2", type=float, default=1e-4)
    p.add_argument("--lambda-locality", type=float, default=0.02)
    p.add_argument("--edge-width", type=float, default=0.30)
    p.add_argument("--gamma", type=float, default=0.010)
    p.add_argument("--source-amp", type=float, default=0.35)
    p.add_argument("--pulse-t0", type=float, default=6.0)
    p.add_argument("--pulse-kind", default="ricker")
    p.add_argument("--tmax", type=float, default=28.0)
    p.add_argument("--dt", type=float, default=0.006)
    p.add_argument("--record-stride", type=int, default=3)
    p.add_argument("--n-modes-dyn", type=int, default=500)
    p.add_argument("--n-fit-modes", type=int, default=80)
    p.add_argument("--n-spectral-bins", type=int, default=10)
    p.add_argument("--low-mode-cut", type=int, default=25)
    p.add_argument("--nrings", type=int, default=10)
    p.add_argument("--ring-rmin", type=float, default=0.16)
    p.add_argument("--rmax-frac", type=float, default=0.82)
    p.add_argument("--t-ignore", type=float, default=5.8)
    p.add_argument("--baseline-until", type=float, default=4.0)
    p.add_argument("--min-snr", type=float, default=1.20)

    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    runs_dir = outdir / "runs"
    results_dir = outdir / "results"
    figures_dir = outdir / "figures"
    logs_dir = outdir / "logs"
    for d in [runs_dir, results_dir, figures_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    eta_list = parse_float_list(args.eta_edge_list)
    sw_list = parse_float_list(args.source_width_list)
    ps_list = parse_float_list(args.pulse_sigma_list)
    grid = list(itertools.product(eta_list, sw_list, ps_list))

    print("=" * 100)
    print("BuP true edge-Hessian v10 scan")
    print("=" * 100)
    print(f"Runs: {len(grid)}")
    print(f"v10 script: {args.v10_script}")
    print(f"output dir: {outdir}")

    rows = []
    for run_idx, (eta, sw, ps) in enumerate(grid, start=1):
        tag = f"eta{eta:.3g}_sw{sw:.3g}_ps{ps:.3g}".replace(".", "p")
        run_dir = runs_dir / tag
        log_path = logs_dir / f"{tag}.log"
        summary_path = run_dir / "results" / "summary_true_edge_hessian_v10.json"

        cmd = [
            sys.executable, str(args.v10_script),
            "--N-side", str(args.N_side),
            "--extent", str(args.extent),
            "--ell", str(args.ell),
            "--cutoff", str(args.cutoff),
            "--eta-edge", str(eta),
            "--mass2", str(args.mass2),
            "--lambda-locality", str(args.lambda_locality),
            "--edge-width", str(args.edge_width),
            "--gamma", str(args.gamma),
            "--source-amp", str(args.source_amp),
            "--source-width", str(sw),
            "--pulse-t0", str(args.pulse_t0),
            "--pulse-sigma", str(ps),
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
            "--output-dir", str(run_dir),
        ]

        print(f"\n[{run_idx}/{len(grid)}] {tag}")
        if args.dry_run:
            print(" ".join(cmd))
            continue

        if args.skip_existing and summary_path.exists():
            code = 0
            print("existing summary found; skipping")
        else:
            code = run_one(cmd, log_path)
            print(f"return code: {code}")

        s = read_summary(run_dir)
        row = {
            "run": tag,
            "return_code": code,
            "eta_edge": eta,
            "source_width": sw,
            "pulse_sigma": ps,
            "run_dir": str(run_dir),
            "log_path": str(log_path),
        }
        if s is not None:
            row.update({
                "N_nodes": s.get("N_nodes"),
                "N_edges": s.get("N_edges"),
                "c_edge": safe_float(s.get("dispersion_c_graph")),
                "meff2": safe_float(s.get("dispersion_intercept_mass2")),
                "R2_disp": safe_float(s.get("dispersion_r2")),
                "plus_corr": safe_float(get_nested(s, ["plus", "modal_corr_absJ_qpeak"])),
                "cross_corr": safe_float(get_nested(s, ["cross", "modal_corr_absJ_qpeak"])),
                "plus_low25_energy": safe_float(get_nested(s, ["plus", "low25_energy"])),
                "cross_low25_energy": safe_float(get_nested(s, ["cross", "low25_energy"])),
                "plus_k_weighted": safe_float(get_nested(s, ["plus", "spectral_stats", "energy_weighted_k"])),
                "cross_k_weighted": safe_float(get_nested(s, ["cross", "spectral_stats", "energy_weighted_k"])),
                "plus_dominant_mode": get_nested(s, ["plus", "spectral_stats", "dominant_mode"], None),
                "cross_dominant_mode": get_nested(s, ["cross", "spectral_stats", "dominant_mode"], None),
                "plus_dominant_frac": safe_float(get_nested(s, ["plus", "spectral_stats", "dominant_mode_energy_fraction"])),
                "cross_dominant_frac": safe_float(get_nested(s, ["cross", "spectral_stats", "dominant_mode_energy_fraction"])),
                "plus_packet_v": safe_float(get_nested(s, ["plus", "packet_v"])),
                "cross_packet_v": safe_float(get_nested(s, ["cross", "packet_v"])),
                "plus_packet_r2": safe_float(get_nested(s, ["plus", "packet_r2"])),
                "cross_packet_r2": safe_float(get_nested(s, ["cross", "packet_r2"])),
            })
            row["mean_corr"] = np.nanmean([row["plus_corr"], row["cross_corr"]])
            row["mean_low25_energy"] = np.nanmean([row["plus_low25_energy"], row["cross_low25_energy"]])
            row["mean_k_weighted"] = np.nanmean([row["plus_k_weighted"], row["cross_k_weighted"]])
            row["score"] = row["mean_corr"] - 0.15 * row["mean_low25_energy"]
            print(
                f"c={row['c_edge']:.6f}, R2={row['R2_disp']:.6f}, "
                f"corr+= {row['plus_corr']:.4f}, corrx= {row['cross_corr']:.4f}, "
                f"low25 mean={row['mean_low25_energy']:.4f}"
            )
        else:
            print("missing summary")
        rows.append(row)

    if args.dry_run:
        return

    df = pd.DataFrame(rows)
    df.to_csv(results_dir / "scan_v10_summary.csv", index=False)
    with open(results_dir / "scan_v10_summary.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    ok = df[(df["return_code"] == 0) & np.isfinite(df.get("mean_corr", np.nan))].copy()
    if len(ok) > 0:
        ok = ok.sort_values("score", ascending=False)
        ok.to_csv(results_dir / "scan_v10_ranked.csv", index=False)
        best = ok.iloc[0].to_dict()
        with open(results_dir / "scan_v10_best.json", "w", encoding="utf-8") as f:
            json.dump(best, f, indent=2)

        scatter_plot(ok, "source_width", "mean_corr", figures_dir / "fig_scan_corr_vs_source_width.png", "Mean modal response vs source width", "source width", "mean corr")
        scatter_plot(ok, "pulse_sigma", "mean_corr", figures_dir / "fig_scan_corr_vs_pulse_sigma.png", "Mean modal response vs pulse duration", "pulse sigma", "mean corr")
        scatter_plot(ok, "mean_k_weighted", "mean_corr", figures_dir / "fig_scan_corr_vs_kweighted.png", "Mean modal response vs weighted k", "weighted k", "mean corr")
        scatter_plot(ok, "mean_low25_energy", "mean_corr", figures_dir / "fig_scan_corr_vs_low25.png", "Mean modal response vs low-mode energy", "low25 energy fraction", "mean corr")
        heatmap_pivot(ok, "source_width", "pulse_sigma", "mean_corr", figures_dir / "fig_scan_heatmap_corr_sw_ps.png", "Mean modal response")
        heatmap_pivot(ok, "source_width", "pulse_sigma", "mean_k_weighted", figures_dir / "fig_scan_heatmap_k_sw_ps.png", "Mean weighted k")

        md = [
            "# BuP true edge-Hessian v10 scan",
            "",
            "## Best run by score",
            "",
        ]
        for key in [
            "run", "eta_edge", "source_width", "pulse_sigma", "c_edge", "meff2", "R2_disp",
            "plus_corr", "cross_corr", "mean_corr", "plus_low25_energy", "cross_low25_energy",
            "mean_low25_energy", "plus_k_weighted", "cross_k_weighted", "mean_k_weighted",
            "plus_dominant_mode", "cross_dominant_mode", "score",
        ]:
            if key in best:
                md.append(f"- `{key}`: `{best[key]}`")
        md += [
            "",
            "## Interpretation",
            "",
            "The scan tests whether the true edge-Hessian dynamic generation result is stable under variations of source localization, pulse duration, and edge rigidity.",
            "The packet/ring diagnostic is secondary; the main validation criteria are dispersion R2, modal response correlation, and spectral localization.",
        ]
        (results_dir / "scan_v10_summary.md").write_text("\n".join(md), encoding="utf-8")

        print("\nBest run:")
        for key in ["run", "score", "mean_corr", "mean_low25_energy", "mean_k_weighted", "c_edge", "R2_disp"]:
            print(f"  {key}: {best.get(key)}")

    print("\nFiles written:")
    for pth in [results_dir / "scan_v10_summary.csv", results_dir / "scan_v10_ranked.csv", results_dir / "scan_v10_best.json", results_dir / "scan_v10_summary.md", figures_dir]:
        print(pth)


if __name__ == "__main__":
    main()

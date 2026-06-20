#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP gravitational-wave dynamic generation
v11 finite-size scan for true edge-Hessian dynamics.

Goal
----
Repeat the best v10 regime over increasing graph sizes:
    eta_edge = 0.75
    source_width = 0.12
    pulse_sigma = 0.25

Default N-side values:
    11, 16, 20, 25
corresponding to N = 121, 256, 400, 625.

For each size, this wrapper runs bup_gw_true_edge_hessian_v10.py and extracts:
    - c_edge, meff2, R2_disp
    - plus/cross modal correlations
    - plus/cross k_weighted
    - plus/cross packet velocity and R2
    - group velocity estimated as d omega / d k near k_weighted

Usage
-----
python3 scan_finite_size_true_edge_hessian_v11.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/bup_gw_true_edge_hessian_v10.py \
  --output-dir papers/paper22_bup_gravitational_waves/results/finite_size_true_edge_hessian_v11

Notes
-----
N=625 may be heavy depending on how v10 diagonalizes K_ef. If it becomes slow,
first run with --N-side-list 11,16,20, then run 25 separately.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_float_list(s: str) -> list[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def safe_float(x: Any, default: float = float("nan")) -> float:
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def find_file(root: Path, candidates: Iterable[str]) -> Optional[Path]:
    for c in candidates:
        p = root / c
        if p.exists():
            return p
    for c in candidates:
        hits = list(root.rglob(Path(c).name))
        if hits:
            return hits[0]
    return None


def load_summary(run_dir: Path) -> Dict[str, Any]:
    p = find_file(run_dir, [
        "results/summary_true_edge_hessian_v10.json",
        "summary_true_edge_hessian_v10.json",
        "results/summary.json",
        "summary.json",
    ])
    if p is None:
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}


def load_modes(run_dir: Path) -> Optional[pd.DataFrame]:
    p = find_file(run_dir, [
        "results/v10_edge_modes.csv",
        "v10_edge_modes.csv",
        "results/edge_modes.csv",
        "edge_modes.csv",
    ])
    if p is None:
        return None
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def pick_col(df: pd.DataFrame, names: list[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    # fuzzy fallback
    for c in df.columns:
        lc = c.lower()
        if any(n.lower() in lc for n in names):
            return c
    return None


def estimate_group_velocity(modes: Optional[pd.DataFrame], k_target: float, window: int = 14) -> Tuple[float, float, int]:
    """Estimate d omega / d k near k_target from the modes table.

    Returns: (v_group, local_R2, n_used)
    """
    if modes is None or modes.empty or not np.isfinite(k_target):
        return float("nan"), float("nan"), 0

    k_col = pick_col(modes, ["k_eff", "k", "k_edge", "k_graph"])
    omega_col = pick_col(modes, ["omega", "omega_hessian", "omega_edge"])
    omega2_col = pick_col(modes, ["omega2", "omega2_hessian", "omega2_edge", "eigenvalue"])
    k2_col = pick_col(modes, ["k2_eff", "k2", "k2_edge", "k2_graph"])

    df = modes.copy()
    if k_col is None and k2_col is not None:
        vals = pd.to_numeric(df[k2_col], errors="coerce").to_numpy(dtype=float)
        df["__k__"] = np.sqrt(np.maximum(vals, 0.0))
        k_col = "__k__"
    if omega_col is None and omega2_col is not None:
        vals = pd.to_numeric(df[omega2_col], errors="coerce").to_numpy(dtype=float)
        df["__omega__"] = np.sqrt(np.maximum(vals, 0.0))
        omega_col = "__omega__"

    if k_col is None or omega_col is None:
        return float("nan"), float("nan"), 0

    k = pd.to_numeric(df[k_col], errors="coerce").to_numpy(dtype=float)
    omega = pd.to_numeric(df[omega_col], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(k) & np.isfinite(omega) & (k > 0) & (omega > 0)
    k = k[ok]
    omega = omega[ok]
    if len(k) < 5:
        return float("nan"), float("nan"), int(len(k))

    order = np.argsort(np.abs(k - k_target))
    n = min(max(5, window), len(order))
    sel = order[:n]
    x = k[sel]
    y = omega[sel]
    # Guard against duplicate k values or tiny spread.
    if np.nanmax(x) - np.nanmin(x) < 1e-10:
        return float("nan"), float("nan"), int(n)

    X = np.column_stack([x, np.ones_like(x)])
    slope, intercept = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ np.array([slope, intercept])
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(r2), int(n)




def analytic_group_velocity(c_edge: float, meff2: float, k_target: float) -> float:
    """Analytic group velocity from omega^2 = c^2 k^2 + m_eff^2.

    v_g = d omega/dk = c^2 k / sqrt(c^2 k^2 + m_eff^2).
    This is the physically meaningful group velocity for the fitted dispersion branch.
    The earlier local finite-difference estimate can be biased by branch mixing and
    degeneracies in the discrete graph spectrum, especially in high-k bands.
    """
    try:
        c = float(c_edge)
        m2 = float(meff2)
        k = float(k_target)
    except Exception:
        return float("nan")
    if not (np.isfinite(c) and np.isfinite(m2) and np.isfinite(k)) or c <= 0 or k <= 0:
        return float("nan")
    omega = np.sqrt(max(c*c*k*k + m2, 0.0))
    if omega <= 0:
        return float("nan")
    return float((c*c*k) / omega)

def get_first(summary: Dict[str, Any], keys: list[str], default: float = float("nan")) -> float:
    """Extract a numeric field from a v10 summary JSON.

    v10 has both flat fields such as ``dispersion_c_graph`` and nested fields:
        plus.modal_corr_absJ_qpeak
        plus.spectral_stats.energy_weighted_k
    This helper supports top-level, dotted paths and recursive fallback.
    """
    def walk(obj: Any, parts: list[str]) -> Any:
        cur = obj
        for part in parts:
            if not isinstance(cur, dict) or part not in cur:
                return None
            cur = cur[part]
        return cur

    for k in keys:
        if "." in k:
            v = walk(summary, k.split("."))
            if v is not None:
                return safe_float(v, default)
        elif k in summary:
            return safe_float(summary[k], default)

    def recursive_find(obj: Any, target: str) -> Any:
        if not isinstance(obj, dict):
            return None
        if target in obj:
            return obj[target]
        for v in obj.values():
            hit = recursive_find(v, target)
            if hit is not None:
                return hit
        return None

    for k in keys:
        if "." not in k:
            v = recursive_find(summary, k)
            if v is not None:
                return safe_float(v, default)
    return default


def pol(summary: Dict[str, Any], pol_name: str, keys: list[str], default: float = float("nan")) -> float:
    """Extract a field under summary['plus'] or summary['cross']."""
    block = summary.get(pol_name, {})
    if not isinstance(block, dict):
        return default
    return get_first(block, keys, default)


def build_v10_command(args: argparse.Namespace, n_side: int, run_dir: Path) -> list[str]:
    return [
        sys.executable,
        str(args.v10_script),
        "--N-side", str(n_side),
        "--extent", str(args.extent),
        "--ell", str(args.ell),
        "--cutoff", str(args.cutoff),
        "--eta-edge", str(args.eta_edge),
        "--mass2", str(args.mass2),
        "--lambda-locality", str(args.lambda_locality),
        "--edge-width", str(args.edge_width),
        "--gamma", str(args.gamma),
        "--source-amp", str(args.source_amp),
        "--source-width", str(args.source_width),
        "--pulse-t0", str(args.pulse_t0),
        "--pulse-sigma", str(args.pulse_sigma),
        "--pulse-kind", str(args.pulse_kind),
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


def plot_metric(df: pd.DataFrame, outpath: Path, ycols: list[str], title: str, ylabel: str) -> None:
    plt.figure(figsize=(7, 4.5))
    for y in ycols:
        if y in df.columns:
            plt.plot(df["N_nodes"], df[y], marker="o", label=y)
    plt.xlabel("N nodes")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--v10-script", type=Path, required=True)
    p.add_argument("--N-side-list", default="11,16,20,25")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--skip-existing", action="store_true")

    # Best v10 scan point defaults.
    p.add_argument("--eta-edge", type=float, default=0.75)
    p.add_argument("--source-width", type=float, default=0.12)
    p.add_argument("--pulse-sigma", type=float, default=0.25)

    # v10 graph and dynamics defaults.
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
    p.add_argument("--group-window", type=int, default=14)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.output_dir
    runs_dir = outdir / "runs"
    results_dir = outdir / "results"
    figures_dir = outdir / "figures"
    runs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    n_sides = parse_int_list(args.N_side_list)
    rows = []

    for n_side in n_sides:
        n_nodes = n_side * n_side
        run_name = f"N{n_nodes}_side{n_side}_eta{str(args.eta_edge).replace('.', 'p')}_sw{str(args.source_width).replace('.', 'p')}_ps{str(args.pulse_sigma).replace('.', 'p')}"
        run_dir = runs_dir / run_name
        summary_file = find_file(run_dir, ["results/summary_true_edge_hessian_v10.json", "summary_true_edge_hessian_v10.json"])

        if args.skip_existing and summary_file is not None:
            print(f"[skip] {run_name}")
        else:
            print("=" * 100)
            print(f"Running finite-size v11: {run_name}")
            print("=" * 100)
            run_dir.mkdir(parents=True, exist_ok=True)
            cmd = build_v10_command(args, n_side, run_dir)
            print(" ".join(cmd))
            proc = subprocess.run(cmd, text=True)
            if proc.returncode != 0:
                print(f"[warning] run failed: {run_name}, returncode={proc.returncode}")

        summary = load_summary(run_dir)
        modes = load_modes(run_dir)

        plus_k = get_first(summary, [
            "plus_energy_weighted_k",
            "plus.spectral_stats.energy_weighted_k",
            "plus_k_weighted", "k_weighted_plus", "plus_spectral_k_weighted",
        ])
        cross_k = get_first(summary, [
            "cross_energy_weighted_k",
            "cross.spectral_stats.energy_weighted_k",
            "cross_k_weighted", "k_weighted_cross", "cross_spectral_k_weighted",
        ])
        # Numeric local slope is retained as a diagnostic only. On discrete graph spectra,
        # nearby modes can belong to different branches/degenerate shells, so this value
        # can underestimate the physical group speed. The primary group velocity reported
        # below is the analytic derivative of the fitted dispersion branch.
        c_edge = get_first(summary, ["c_edge", "dispersion_c_graph", "c_graph", "dispersion_c", "c_graph_from_hessian"])
        meff2 = get_first(summary, ["meff2", "dispersion_intercept_mass2", "dispersion_meff2", "m_eff2", "mass2_eff"])
        r2_disp = get_first(summary, ["R2_disp", "dispersion_r2", "r2_disp", "dispersion_R2"])

        plus_vg_numeric, plus_vg_r2, plus_vg_n = estimate_group_velocity(modes, plus_k, args.group_window)
        cross_vg_numeric, cross_vg_r2, cross_vg_n = estimate_group_velocity(modes, cross_k, args.group_window)
        plus_vg = analytic_group_velocity(c_edge, meff2, plus_k)
        cross_vg = analytic_group_velocity(c_edge, meff2, cross_k)

        plus_corr = get_first(summary, ["plus_corr", "plus.modal_corr_absJ_qpeak", "plus_modal_corr", "corr_plus"])
        cross_corr = get_first(summary, ["cross_corr", "cross.modal_corr_absJ_qpeak", "cross_modal_corr", "corr_cross"])
        plus_low25 = get_first(summary, ["plus_low25_energy", "plus.low25_energy", "low25_energy_plus"])
        cross_low25 = get_first(summary, ["cross_low25_energy", "cross.low25_energy", "low25_energy_cross"])

        row = {
            "run": run_name,
            "N_side": n_side,
            "N_nodes": n_nodes,
            "eta_edge": args.eta_edge,
            "source_width": args.source_width,
            "pulse_sigma": args.pulse_sigma,
            "c_edge": get_first(summary, ["c_edge", "dispersion_c_graph", "c_graph", "dispersion_c", "c_graph_from_hessian"]),
            "meff2": get_first(summary, ["meff2", "mass2_eff", "dispersion_intercept_mass2", "meff2_disp"]),
            "R2_disp": get_first(summary, ["R2_disp", "dispersion_r2", "r2_disp"]),
            "plus_corr": plus_corr,
            "cross_corr": cross_corr,
            "mean_corr": np.nanmean([plus_corr, cross_corr]),
            "plus_low25_energy": plus_low25,
            "cross_low25_energy": cross_low25,
            "mean_low25_energy": np.nanmean([plus_low25, cross_low25]),
            "plus_k_weighted": plus_k,
            "cross_k_weighted": cross_k,
            "mean_k_weighted": np.nanmean([plus_k, cross_k]),
            "plus_dominant_mode": get_first(summary, ["plus_dominant_mode", "plus.spectral_stats.dominant_mode", "dominant_mode_plus"]),
            "cross_dominant_mode": get_first(summary, ["cross_dominant_mode", "cross.spectral_stats.dominant_mode", "dominant_mode_cross"]),
            "plus_dominant_frac": get_first(summary, ["plus_dominant_frac", "plus.spectral_stats.dominant_mode_energy_fraction", "dominant_frac_plus"]),
            "cross_dominant_frac": get_first(summary, ["cross_dominant_frac", "cross.spectral_stats.dominant_mode_energy_fraction", "dominant_frac_cross"]),
            "plus_packet_v": get_first(summary, ["plus_packet_v", "plus.packet_v", "plus_packet_velocity", "packet_v_plus"]),
            "cross_packet_v": get_first(summary, ["cross_packet_v", "cross.packet_v", "cross_packet_velocity", "packet_v_cross"]),
            "plus_packet_R2": get_first(summary, ["plus_packet_R2", "plus.packet_r2", "plus_packet_r2", "packet_R2_plus"]),
            "cross_packet_R2": get_first(summary, ["cross_packet_R2", "cross.packet_r2", "cross_packet_r2", "packet_R2_cross"]),
            "plus_group_velocity": plus_vg,
            "plus_group_velocity_numeric": plus_vg_numeric,
            "plus_group_R2": plus_vg_r2,
            "plus_group_points": plus_vg_n,
            "cross_group_velocity": cross_vg,
            "cross_group_velocity_numeric": cross_vg_numeric,
            "cross_group_R2": cross_vg_r2,
            "cross_group_points": cross_vg_n,
            "run_dir": str(run_dir),
        }
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("N_nodes")
    df.to_csv(results_dir / "finite_size_v12_summary.csv", index=False)

    stable_cols = [
        "N_nodes", "c_edge", "meff2", "R2_disp",
        "plus_corr", "cross_corr", "mean_corr",
        "plus_low25_energy", "cross_low25_energy", "mean_low25_energy",
        "plus_k_weighted", "cross_k_weighted", "mean_k_weighted",
        "plus_group_velocity", "cross_group_velocity",
        "plus_packet_v", "cross_packet_v", "plus_packet_R2", "cross_packet_R2",
    ]
    df[[c for c in stable_cols if c in df.columns]].to_csv(results_dir / "finite_size_v12_key_metrics.csv", index=False)

    plot_metric(df, figures_dir / "fig_v12_corr_vs_N.png", ["plus_corr", "cross_corr", "mean_corr"], "Modal response vs graph size", "corr(|J_n|, q_peak)")
    plot_metric(df, figures_dir / "fig_v12_low25_vs_N.png", ["plus_low25_energy", "cross_low25_energy", "mean_low25_energy"], "Low-mode energy vs graph size", "Low25 energy fraction")
    plot_metric(df, figures_dir / "fig_v12_dispersion_vs_N.png", ["c_edge", "R2_disp"], "Dispersion diagnostics vs graph size", "value")
    plot_metric(df, figures_dir / "fig_v12_group_velocity_vs_N.png", ["plus_group_velocity", "cross_group_velocity"], "Group velocity near excited k band", "d omega / d k")
    plot_metric(df, figures_dir / "fig_v12_packet_vs_N.png", ["plus_packet_v", "cross_packet_v"], "Packet diagnostic velocity vs graph size", "packet v")

    summary = {
        "N_side_list": n_sides,
        "eta_edge": args.eta_edge,
        "source_width": args.source_width,
        "pulse_sigma": args.pulse_sigma,
        "n_runs": int(len(df)),
        "mean_corr_min": safe_float(df["mean_corr"].min()) if "mean_corr" in df else float("nan"),
        "mean_corr_mean": safe_float(df["mean_corr"].mean()) if "mean_corr" in df else float("nan"),
        "R2_disp_min": safe_float(df["R2_disp"].min()) if "R2_disp" in df else float("nan"),
        "meff2_mean": safe_float(df["meff2"].mean()) if "meff2" in df else float("nan"),
        "interpretation": "Finite-size scan for true edge-Hessian dynamic generation at the best v10 source point.",
    }
    (results_dir / "finite_size_v12_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        "# True edge-Hessian v11 finite-size scan",
        "",
        f"Best v10 source point repeated over `N_side={n_sides}`.",
        "",
        "## Fixed parameters",
        "",
        f"- `eta_edge`: `{args.eta_edge}`",
        f"- `source_width`: `{args.source_width}`",
        f"- `pulse_sigma`: `{args.pulse_sigma}`",
        "",
        "## Aggregate summary",
        "",
    ]
    for k, v in summary.items():
        md.append(f"- `{k}`: `{v}`")
    md.append("")
    md.append("## Key metrics")
    md.append("")
    md.append(df[[c for c in stable_cols if c in df.columns]].to_markdown(index=False))
    (results_dir / "finite_size_v12_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Finite-size v12 scan complete")
    print("=" * 100)
    print(df[[c for c in stable_cols if c in df.columns]].to_string(index=False))
    print(f"\nWrote outputs to: {outdir}")


if __name__ == "__main__":
    main()

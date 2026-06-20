#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP GW true edge-Hessian fixed-density mass-gap scan v14

Goal
----
Disentangle the small fitted mass gap m_eff^2 observed in the true edge-Hessian
wave dispersion:

    omega^2 = c_edge^2 k^2 + m_eff^2

from three possible origins:
  1. UV/discretization gap: m_eff^2 -> 0 as lattice spacing a -> 0.
  2. IR/finite-volume gap: m_eff^2 -> 0 as physical size L -> infinity.
  3. Residual physical BuP gap: m_eff^2 -> m0^2 > 0.

This script orchestrates repeated runs of bup_gw_true_edge_hessian_v10.py with
controlled geometry families, extracts dispersion and modal-response metrics, and
fits simple scaling laws.

Scan families
-------------
fixed-volume:
    extent fixed, N_side varied. The lattice spacing a changes.
fixed-density:
    lattice spacing a fixed approximately, N_side varied, extent varied.
fixed-N-density:
    N_side fixed, extent varied. Density/lattice spacing changes at fixed node count.
combined:
    all of the above.

Examples
--------
python3 scan_fixed_density_mass_gap_v14.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/bup_gw_true_edge_hessian_v10.py \
  --scan-type fixed-density \
  --fixed-density-spacing 0.1052631579 \
  --fixed-density-N-side-list 11,16,20,25 \
  --output-dir papers/paper22_bup_gravitational_waves/results/mass_gap_v14_fixed_density
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


# ----------------------------- helpers -------------------------------------

def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def safe_name(x: Any) -> str:
    return str(x).replace("-", "m").replace(".", "p").replace(",", "_")


def get_nested(d: Dict[str, Any], path: str, default=np.nan):
    cur: Any = d
    for part in path.split("."):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return default
    return cur


def first_existing(d: Dict[str, Any], paths: Sequence[str], default=np.nan):
    for p in paths:
        val = get_nested(d, p, None)
        if val is not None:
            return val
    return default


def to_float(x: Any, default=np.nan) -> float:
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def read_summary_json(run_dir: Path) -> Optional[Dict[str, Any]]:
    candidates = [
        run_dir / "results" / "summary_true_edge_hessian_v10.json",
        run_dir / "summary_true_edge_hessian_v10.json",
        run_dir / "results" / "summary.json",
        run_dir / "summary.json",
    ]
    for path in candidates:
        if path.exists():
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    return None


def extract_v10_metrics(summary: Dict[str, Any]) -> Dict[str, float]:
    c_edge = to_float(first_existing(summary, [
        "dispersion_c_graph", "c_edge", "c_graph", "dispersion.c_edge", "dispersion.c_graph"
    ]))
    meff2 = to_float(first_existing(summary, [
        "dispersion_intercept_mass2", "meff2", "m_eff2", "mass2", "dispersion.meff2", "dispersion.mass2"
    ]))
    r2 = to_float(first_existing(summary, [
        "dispersion_r2", "R2_disp", "r2_disp", "dispersion.R2", "dispersion.r2"
    ]))
    plus_corr = to_float(first_existing(summary, [
        "plus.modal_corr_absJ_qpeak", "plus_corr", "plus.modal_corr", "modal_plus_corr"
    ]))
    cross_corr = to_float(first_existing(summary, [
        "cross.modal_corr_absJ_qpeak", "cross_corr", "cross.modal_corr", "modal_cross_corr"
    ]))
    plus_low25 = to_float(first_existing(summary, [
        "plus.low_mode_energy", "plus.low25_energy", "plus_low25_energy", "plus.spectral_stats.low25_energy"
    ]))
    cross_low25 = to_float(first_existing(summary, [
        "cross.low_mode_energy", "cross.low25_energy", "cross_low25_energy", "cross.spectral_stats.low25_energy"
    ]))
    plus_k = to_float(first_existing(summary, [
        "plus.spectral_stats.energy_weighted_k", "plus.k_weighted", "plus_k_weighted"
    ]))
    cross_k = to_float(first_existing(summary, [
        "cross.spectral_stats.energy_weighted_k", "cross.k_weighted", "cross_k_weighted"
    ]))
    plus_packet_v = to_float(first_existing(summary, [
        "plus.packet_v", "plus_packet_v", "plus.packet.velocity", "plus.packet.fit_v"
    ]))
    cross_packet_v = to_float(first_existing(summary, [
        "cross.packet_v", "cross_packet_v", "cross.packet.velocity", "cross.packet.fit_v"
    ]))
    plus_packet_r2 = to_float(first_existing(summary, [
        "plus.packet_r2", "plus_packet_R2", "plus_packet_r2", "plus.packet.r2"
    ]))
    cross_packet_r2 = to_float(first_existing(summary, [
        "cross.packet_r2", "cross_packet_R2", "cross_packet_r2", "cross.packet.r2"
    ]))
    return {
        "c_edge": c_edge,
        "meff2": meff2,
        "R2_disp": r2,
        "plus_corr": plus_corr,
        "cross_corr": cross_corr,
        "mean_corr": np.nanmean([plus_corr, cross_corr]),
        "plus_low25_energy": plus_low25,
        "cross_low25_energy": cross_low25,
        "mean_low25_energy": np.nanmean([plus_low25, cross_low25]),
        "plus_k_weighted": plus_k,
        "cross_k_weighted": cross_k,
        "mean_k_weighted": np.nanmean([plus_k, cross_k]),
        "plus_packet_v": plus_packet_v,
        "cross_packet_v": cross_packet_v,
        "plus_packet_R2": plus_packet_r2,
        "cross_packet_R2": cross_packet_r2,
    }


def analytic_group_velocity(c_edge: float, meff2: float, k: float) -> float:
    if not np.isfinite(c_edge) or not np.isfinite(meff2) or not np.isfinite(k):
        return np.nan
    omega = math.sqrt(max(c_edge * c_edge * k * k + meff2, 0.0))
    if omega <= 0:
        return np.nan
    return (c_edge * c_edge * k) / omega


@dataclass
class Case:
    family: str
    N_side: int
    extent: float
    ell: float
    cutoff: float

    @property
    def N_nodes(self) -> int:
        return self.N_side * self.N_side

    @property
    def L(self) -> float:
        return 2.0 * self.extent

    @property
    def a(self) -> float:
        return self.L / max(self.N_side - 1, 1)

    @property
    def invL2(self) -> float:
        return 1.0 / (self.L * self.L) if self.L > 0 else np.nan

    def run_name(self, eta_edge: float, source_width: float, pulse_sigma: float) -> str:
        return (
            f"{self.family}_N{self.N_nodes}_side{self.N_side}_"
            f"extent{safe_name(self.extent)}_a{safe_name(round(self.a, 6))}_"
            f"eta{safe_name(eta_edge)}_sw{safe_name(source_width)}_ps{safe_name(pulse_sigma)}"
        )


def build_cases(args) -> List[Case]:
    cases: List[Case] = []
    scan_types = [s.strip() for s in args.scan_type.split(",") if s.strip()]
    if "combined" in scan_types:
        scan_types = ["fixed-volume", "fixed-density", "fixed-N-density"]

    if "fixed-volume" in scan_types:
        for ns in parse_int_list(args.fixed_volume_N_side_list):
            cases.append(Case(
                family="fixed_volume",
                N_side=ns,
                extent=args.fixed_volume_extent,
                ell=args.ell,
                cutoff=args.cutoff,
            ))

    if "fixed-density" in scan_types:
        # Keep a approximately fixed by setting extent = a*(N_side-1)/2.
        a0 = args.fixed_density_spacing
        for ns in parse_int_list(args.fixed_density_N_side_list):
            extent = 0.5 * a0 * (ns - 1)
            cases.append(Case(
                family="fixed_density",
                N_side=ns,
                extent=extent,
                ell=args.ell,
                cutoff=args.cutoff,
            ))

    if "fixed-N-density" in scan_types:
        ns = args.fixed_N_side
        for ext in parse_float_list(args.fixed_N_extent_list):
            cases.append(Case(
                family="fixed_N_density",
                N_side=ns,
                extent=ext,
                ell=args.ell,
                cutoff=args.cutoff,
            ))

    # remove exact duplicates while preserving order
    seen = set()
    unique: List[Case] = []
    for c in cases:
        key = (c.family, c.N_side, round(c.extent, 12), round(c.ell, 12), round(c.cutoff, 12))
        if key not in seen:
            seen.add(key)
            unique.append(c)
    return unique


def run_v10(args, case: Case, run_dir: Path) -> int:
    cmd = [
        sys.executable, str(args.v10_script),
        "--N-side", str(case.N_side),
        "--extent", str(case.extent),
        "--ell", str(case.ell),
        "--cutoff", str(case.cutoff),
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
    print("\n" + "=" * 100)
    print(f"Running {run_dir.name}")
    print(" ".join(cmd))
    print("=" * 100)
    return subprocess.call(cmd)


# ----------------------------- fitting -------------------------------------

def fit_linear_model(y: np.ndarray, X: np.ndarray, labels: List[str]) -> Dict[str, Any]:
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < X.shape[1]:
        return {"ok": False, "reason": "not_enough_points", "n_used": int(ok.sum())}
    beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    pred = X[ok] @ beta
    ss_res = float(np.sum((y[ok] - pred) ** 2))
    ss_tot = float(np.sum((y[ok] - y[ok].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {
        "ok": True,
        "n_used": int(ok.sum()),
        "labels": labels,
        "coefficients": {lab: float(val) for lab, val in zip(labels, beta)},
        "r2": float(r2),
        "rmse": float(np.sqrt(np.mean((y[ok] - pred) ** 2))),
    }


def fit_power_law_with_offset(a: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
    """Grid over p, fit y = m0 + A a^p."""
    ok = np.isfinite(a) & np.isfinite(y) & (a > 0)
    if ok.sum() < 3:
        return {"ok": False, "reason": "not_enough_points", "n_used": int(ok.sum())}
    aa = a[ok]
    yy = y[ok]
    best = None
    for p in np.linspace(0.25, 4.0, 151):
        X = np.column_stack([np.ones_like(aa), aa ** p])
        beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
        pred = X @ beta
        rmse = float(np.sqrt(np.mean((yy - pred) ** 2)))
        if best is None or rmse < best["rmse"]:
            ss_res = float(np.sum((yy - pred) ** 2))
            ss_tot = float(np.sum((yy - yy.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            best = {
                "ok": True,
                "n_used": int(ok.sum()),
                "m0_squared": float(beta[0]),
                "A": float(beta[1]),
                "p": float(p),
                "rmse": rmse,
                "r2": float(r2),
            }
    return best or {"ok": False, "reason": "fit_failed"}


def compute_fits(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    y = np.array([to_float(r.get("meff2")) for r in rows], dtype=float)
    a = np.array([to_float(r.get("a")) for r in rows], dtype=float)
    invL2 = np.array([to_float(r.get("invL2")) for r in rows], dtype=float)
    family = np.array([str(r.get("family")) for r in rows])

    fits: Dict[str, Any] = {}
    # global simple decompositions
    X1 = np.column_stack([np.ones_like(y), a])
    fits["global_linear_a"] = fit_linear_model(y, X1, ["m0_squared", "A_a"])
    X2 = np.column_stack([np.ones_like(y), invL2])
    fits["global_linear_invL2"] = fit_linear_model(y, X2, ["m0_squared", "B_invL2"])
    X3 = np.column_stack([np.ones_like(y), a, invL2])
    fits["global_linear_a_invL2"] = fit_linear_model(y, X3, ["m0_squared", "A_a", "B_invL2"])
    fits["global_power_a_offset"] = fit_power_law_with_offset(a, y)

    for fam in sorted(set(family)):
        idx = family == fam
        if idx.sum() >= 3:
            fits[f"{fam}_linear_a"] = fit_linear_model(
                y[idx], np.column_stack([np.ones(idx.sum()), a[idx]]), ["m0_squared", "A_a"]
            )
            fits[f"{fam}_linear_invL2"] = fit_linear_model(
                y[idx], np.column_stack([np.ones(idx.sum()), invL2[idx]]), ["m0_squared", "B_invL2"]
            )
            fits[f"{fam}_power_a_offset"] = fit_power_law_with_offset(a[idx], y[idx])
    return fits


# ----------------------------- outputs -------------------------------------

def write_csv(path: Path, rows: List[Dict[str, Any]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def make_plots(outdir: Path, rows: List[Dict[str, Any]]):
    if plt is None or not rows:
        return
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    a = np.array([to_float(r.get("a")) for r in rows])
    invL2 = np.array([to_float(r.get("invL2")) for r in rows])
    meff2 = np.array([to_float(r.get("meff2")) for r in rows])
    N = np.array([to_float(r.get("N_nodes")) for r in rows])
    corr = np.array([to_float(r.get("mean_corr")) for r in rows])
    c = np.array([to_float(r.get("c_edge")) for r in rows])
    fam = [str(r.get("family")) for r in rows]

    def scatter_by_family(x, y, xlabel, ylabel, title, filename):
        plt.figure(figsize=(6.2, 4.5))
        for f in sorted(set(fam)):
            idx = np.array([ff == f for ff in fam])
            plt.scatter(x[idx], y[idx], label=f)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(figdir / filename, dpi=170)
        plt.close()

    scatter_by_family(a, meff2, "lattice spacing a", r"$m_{eff}^2$", "Mass gap vs lattice spacing", "fig_v14_meff2_vs_a.png")
    scatter_by_family(invL2, meff2, r"$1/L^2$", r"$m_{eff}^2$", "Mass gap vs inverse volume scale", "fig_v14_meff2_vs_invL2.png")
    scatter_by_family(N, meff2, "N nodes", r"$m_{eff}^2$", "Mass gap vs graph size", "fig_v14_meff2_vs_N.png")
    scatter_by_family(N, corr, "N nodes", "mean modal correlation", "Modal response vs graph size", "fig_v14_corr_vs_N.png")
    scatter_by_family(N, c, "N nodes", r"$c_{edge}$", "Edge wave speed vs graph size", "fig_v14_cedge_vs_N.png")


def write_markdown(outdir: Path, rows: List[Dict[str, Any]], fits: Dict[str, Any]):
    lines: List[str] = []
    lines.append("# BuP true edge-Hessian mass-gap scan v13")
    lines.append("")
    lines.append("This scan tests whether the fitted mass gap is a discretization artifact, a finite-volume artifact, or a residual physical gap.")
    lines.append("")
    lines.append("## Key columns")
    lines.append("")
    lines.append("- `a`: lattice spacing, `a = 2 extent / (N_side - 1)`. ")
    lines.append("- `L`: physical box size, `L = 2 extent`. ")
    lines.append("- `meff2`: intercept of the edge-Hessian dispersion fit. ")
    lines.append("- `mean_corr`: mean of plus/cross modal response correlations. ")
    lines.append("")
    if rows:
        lines.append("## Runs")
        lines.append("")
        lines.append("| family | N | extent | a | L | meff2 | c_edge | R2_disp | mean_corr |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in rows:
            lines.append(
                f"| {r.get('family')} | {r.get('N_nodes')} | {float(r.get('extent')):.4g} | "
                f"{float(r.get('a')):.4g} | {float(r.get('L')):.4g} | {float(r.get('meff2')):.6g} | "
                f"{float(r.get('c_edge')):.6g} | {float(r.get('R2_disp')):.6g} | {float(r.get('mean_corr')):.6g} |"
            )
        lines.append("")
    lines.append("## Fits")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps(fits, indent=2))
    lines.append("```")
    lines.append("")
    lines.append("## Interpretation rule")
    lines.append("")
    lines.append("- If fitted `m0_squared` is consistent with zero and the `a` term dominates, the gap is likely a discretization artifact.")
    lines.append("- If fitted `B_invL2` dominates, the gap is likely an IR finite-volume gap.")
    lines.append("- If fitted `m0_squared` remains positive across scan families, the gap may be a residual physical BuP mass scale.")
    (outdir / "results" / "mass_gap_v14_summary.md").write_text("\n".join(lines), encoding="utf-8")


# ----------------------------- CLI -----------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--v10-script", required=True, type=Path)
    p.add_argument("--scan-type", default="fixed-density", help="fixed-density by default; also accepts combined or comma list of fixed-volume,fixed-density,fixed-N-density")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--skip-existing", action="store_true")

    # geometry families
    p.add_argument("--fixed-volume-N-side-list", default="11,16,20,25")
    p.add_argument("--fixed-volume-extent", type=float, default=1.0)
    p.add_argument("--fixed-density-N-side-list", default="11,16,20,25")
    p.add_argument("--fixed-density-spacing", type=float, default=0.1052631579)
    p.add_argument("--fixed-N-side", type=int, default=20)
    p.add_argument("--fixed-N-extent-list", default="0.7,1.0,1.4,2.0")
    p.add_argument("--ell", type=float, default=0.22)
    p.add_argument("--cutoff", type=float, default=0.45)

    # v10 physics / dynamics parameters
    p.add_argument("--eta-edge", type=float, default=1.0)
    p.add_argument("--mass2", type=float, default=1e-4)
    p.add_argument("--lambda-locality", type=float, default=0.02)
    p.add_argument("--edge-width", type=float, default=0.30)
    p.add_argument("--gamma", type=float, default=0.010)
    p.add_argument("--source-amp", type=float, default=0.35)
    p.add_argument("--source-width", type=float, default=0.12)
    p.add_argument("--pulse-t0", type=float, default=6.0)
    p.add_argument("--pulse-sigma", type=float, default=0.25)
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
    return p.parse_args()


def main():
    args = parse_args()
    outdir: Path = args.output_dir
    (outdir / "runs").mkdir(parents=True, exist_ok=True)
    (outdir / "results").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    cases = build_cases(args)
    rows: List[Dict[str, Any]] = []

    for case in cases:
        name = case.run_name(args.eta_edge, args.source_width, args.pulse_sigma)
        run_dir = outdir / "runs" / name
        summary = read_summary_json(run_dir) if args.skip_existing else None
        if summary is None:
            rc = run_v10(args, case, run_dir)
            if rc != 0:
                print(f"WARNING: run failed with code {rc}: {name}")
            summary = read_summary_json(run_dir)
        if summary is None:
            print(f"WARNING: no summary JSON found for {name}")
            row = {
                "run": name,
                "family": case.family,
                "N_side": case.N_side,
                "N_nodes": case.N_nodes,
                "extent": case.extent,
                "L": case.L,
                "a": case.a,
                "invL2": case.invL2,
                "status": "missing_summary",
            }
            rows.append(row)
            continue

        metrics = extract_v10_metrics(summary)
        plus_vg = analytic_group_velocity(metrics["c_edge"], metrics["meff2"], metrics["plus_k_weighted"])
        cross_vg = analytic_group_velocity(metrics["c_edge"], metrics["meff2"], metrics["cross_k_weighted"])
        row = {
            "run": name,
            "family": case.family,
            "N_side": case.N_side,
            "N_nodes": case.N_nodes,
            "extent": case.extent,
            "L": case.L,
            "a": case.a,
            "invL2": case.invL2,
            "eta_edge": args.eta_edge,
            "source_width": args.source_width,
            "pulse_sigma": args.pulse_sigma,
            "status": "ok",
            **metrics,
            "plus_group_velocity": plus_vg,
            "cross_group_velocity": cross_vg,
        }
        rows.append(row)

    fits = compute_fits([r for r in rows if r.get("status") == "ok"])

    write_csv(outdir / "results" / "mass_gap_v14_summary.csv", rows)
    with open(outdir / "results" / "mass_gap_v14_fits.json", "w", encoding="utf-8") as f:
        json.dump(fits, f, indent=2)
    with open(outdir / "results" / "mass_gap_v14_summary.json", "w", encoding="utf-8") as f:
        json.dump({"rows": rows, "fits": fits}, f, indent=2)
    write_markdown(outdir, rows, fits)
    make_plots(outdir, rows)

    print("\n" + "=" * 100)
    print("Mass-gap v14 scan complete")
    print("=" * 100)
    if pd is not None:
        cols = [
            "family", "N_nodes", "extent", "a", "L", "c_edge", "meff2", "R2_disp",
            "mean_corr", "mean_low25_energy", "mean_k_weighted", "plus_group_velocity", "cross_group_velocity",
        ]
        df = pd.DataFrame(rows)
        print(df[[c for c in cols if c in df.columns]].to_string(index=False))
    else:
        for r in rows:
            print(r)
    print("\nFits:")
    print(json.dumps(fits, indent=2))
    print(f"\nWrote outputs to: {outdir}")


if __name__ == "__main__":
    main()

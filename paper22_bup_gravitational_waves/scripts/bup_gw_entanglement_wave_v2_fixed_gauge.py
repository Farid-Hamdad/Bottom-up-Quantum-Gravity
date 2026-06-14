#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 21 — BuP gravitational-wave prototype v2

Fix compared to v1
------------------
v1 measured the ring deformation after an unconstrained MDS reconstruction.
MDS has a rotational gauge freedom: for a cross-polarized ellipse, it rotates
the reconstructed axes onto the principal axes, so the cross mode is converted
into a plus-like mode. This erases the polarization information.

v2 measures the deformation in a fixed physical/emergent gauge using the
original graph coordinates. The perturbation is still applied to the
entanglement graph W_ij(t), but the quadrupole response is extracted from the
induced change of pairwise entanglement distances on the ring.

Expected result
---------------
plus:
  q_plus correlates with input h(t), q_cross ~ 0

cross:
  q_cross correlates with input h(t), q_plus ~ 0

Run
---
cd ~/bottomup

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_entanglement_wave_v2_fixed_gauge.py \
  --polarization plus \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_plus_v2

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_entanglement_wave_v2_fixed_gauge.py \
  --polarization cross \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_cross_v2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def build_grid(n_side: int):
    xs = np.linspace(-1.0, 1.0, n_side)
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def pairwise_geometry(coords: np.ndarray):
    diff = coords[:, None, :] - coords[None, :, :]
    r = np.sqrt(np.sum(diff * diff, axis=-1))
    nvec = np.zeros_like(diff)
    mask = r > 0
    nvec[mask] = diff[mask] / r[mask, None]
    return diff, r, nvec


def base_weights(coords: np.ndarray, ell: float, cutoff: float):
    _, r, _ = pairwise_geometry(coords)
    W = np.exp(-r / ell)
    W[r > cutoff] = 0.0
    np.fill_diagonal(W, 0.0)
    return W


def strain_tensor(polarization: str, phase: float, amplitude: float):
    h = amplitude * np.cos(phase)
    if polarization == "plus":
        return np.array([[h, 0.0], [0.0, -h]])
    if polarization == "cross":
        return np.array([[0.0, h], [h, 0.0]])
    if polarization == "mixed":
        hp = amplitude * np.cos(phase)
        hx = amplitude * np.sin(phase)
        return np.array([[hp, hx], [hx, -hp]])
    raise ValueError("polarization must be plus, cross or mixed")


def perturb_weights(coords: np.ndarray, W0: np.ndarray, H: np.ndarray):
    _, r, nvec = pairwise_geometry(coords)
    q = np.einsum("...a,ab,...b->...", nvec, H, nvec)
    # Length perturbation: dl/l ~= 0.5 n^T H n.
    # Entanglement weight decreases when effective length increases.
    W = W0 * np.exp(-0.5 * q)
    np.fill_diagonal(W, 0.0)
    return W


def entanglement_distance(W: np.ndarray, l0: float, eps: float = 1e-12):
    Wpos = W[W > 0]
    Wref = np.max(Wpos)
    D = np.full_like(W, np.nan, dtype=float)
    mask = W > 0
    ratio = np.clip((W[mask] + eps) / (Wref + eps), 1e-300, 1.0)
    D[mask] = -l0 * np.log(ratio)
    np.fill_diagonal(D, 0.0)
    return D


def select_ring(coords: np.ndarray, radius: float, width: float):
    r = np.sqrt(np.sum(coords * coords, axis=1))
    idx = np.where(np.abs(r - radius) <= width)[0]
    theta = np.arctan2(coords[idx, 1], coords[idx, 0])
    return idx[np.argsort(theta)]


def ring_quadrupole_fixed_gauge(coords_ring: np.ndarray, D0: np.ndarray, Dt: np.ndarray):
    """
    Extract signed plus/cross response from distance changes in a fixed gauge.

    For a small strain H, the fractional squared-distance change satisfies
        delta(d^2) / d^2 ~= n^T H n
    with
        n^T H_plus n  = h cos(2theta)
        n^T H_cross n = h sin(2theta)

    We regress the measured fractional change against cos(2theta), sin(2theta).
    """
    n = len(coords_ring)
    diff = coords_ring[:, None, :] - coords_ring[None, :, :]
    r2 = np.sum(diff * diff, axis=-1)

    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    mask &= np.isfinite(D0) & np.isfinite(Dt) & (D0 > 0) & (r2 > 0)

    dx = diff[..., 0][mask]
    dy = diff[..., 1][mask]
    theta = np.arctan2(dy, dx)

    y = (Dt[mask] ** 2 - D0[mask] ** 2) / (D0[mask] ** 2)

    X = np.column_stack([np.cos(2 * theta), np.sin(2 * theta), np.ones_like(theta)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    q_plus, q_cross, offset = coef

    pred = X @ coef
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2_fit = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return float(q_plus), float(q_cross), float(offset), float(r2_fit)


def synthetic_deformed_ring(coords_ring: np.ndarray, H: np.ndarray):
    # Visual only: apply half-strain displacement x -> (I + H/2)x.
    return coords_ring @ (np.eye(2) + 0.5 * H).T


def plot_ring(coords_ring, coords_def, outpath, title):
    plt.figure(figsize=(5, 5))
    plt.scatter(coords_ring[:, 0], coords_ring[:, 1], s=14, label="reference")
    plt.scatter(coords_def[:, 0], coords_def[:, 1], s=14, label="perturbed")
    plt.plot(np.r_[coords_def[:, 0], coords_def[0, 0]], np.r_[coords_def[:, 1], coords_def[0, 1]], alpha=0.6)
    lim = max(np.max(np.abs(coords_ring)), np.max(np.abs(coords_def))) * 1.25
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.legend()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-side", type=int, default=31)
    p.add_argument("--ell", type=float, default=0.18)
    p.add_argument("--cutoff", type=float, default=0.45)
    p.add_argument("--l0", type=float, default=0.18)
    p.add_argument("--amplitude", type=float, default=0.03)
    p.add_argument("--omega", type=float, default=1.0)
    p.add_argument("--n-frames", type=int, default=80)
    p.add_argument("--polarization", choices=["plus", "cross", "mixed"], default="plus")
    p.add_argument("--ring-radius", type=float, default=0.55)
    p.add_argument("--ring-width", type=float, default=0.04)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    frames = outdir / "frames"
    outdir.mkdir(parents=True, exist_ok=True)
    frames.mkdir(parents=True, exist_ok=True)

    coords = build_grid(args.N_side)
    W0 = base_weights(coords, args.ell, args.cutoff)
    D0_all = entanglement_distance(W0, args.l0)

    ring_idx = select_ring(coords, args.ring_radius, args.ring_width)
    if len(ring_idx) < 12:
        raise ValueError(f"Ring has too few nodes: {len(ring_idx)}. Increase --ring-width.")

    coords_ring = coords[ring_idx]
    D0 = D0_all[np.ix_(ring_idx, ring_idx)]

    rows = []

    for t in range(args.n_frames):
        phase = 2 * np.pi * args.omega * t / args.n_frames
        h_input = args.amplitude * np.cos(phase)
        H = strain_tensor(args.polarization, phase, args.amplitude)

        Wt = perturb_weights(coords, W0, H)
        Dt_all = entanglement_distance(Wt, args.l0)
        Dt = Dt_all[np.ix_(ring_idx, ring_idx)]

        q_plus, q_cross, offset, fit_r2 = ring_quadrupole_fixed_gauge(coords_ring, D0, Dt)
        ellip = float(np.sqrt(q_plus * q_plus + q_cross * q_cross))

        rows.append({
            "frame": t,
            "phase": phase,
            "h_input": h_input,
            "q_plus": q_plus,
            "q_cross": q_cross,
            "q_amp": ellip,
            "offset": offset,
            "quadrupole_fit_r2": fit_r2,
            "n_ring": int(len(ring_idx)),
        })

        if t in {0, args.n_frames // 4, args.n_frames // 2, 3 * args.n_frames // 4, args.n_frames - 1}:
            coords_def = synthetic_deformed_ring(coords_ring, H)
            plot_ring(coords_ring, coords_def, frames / f"frame_{t:04d}.png",
                      f"{args.polarization}, frame {t}")

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "ring_deformation_timeseries.csv", index=False)

    corr_plus = np.corrcoef(df["h_input"], df["q_plus"])[0, 1] if df["q_plus"].std() > 0 else np.nan
    corr_cross = np.corrcoef(df["h_input"], df["q_cross"])[0, 1] if df["q_cross"].std() > 0 else np.nan

    summary = {
        "N_side": args.N_side,
        "N_nodes": int(len(coords)),
        "n_ring": int(len(ring_idx)),
        "polarization": args.polarization,
        "input_amplitude": float(args.amplitude),
        "q_plus_amplitude": float(0.5 * (df["q_plus"].max() - df["q_plus"].min())),
        "q_cross_amplitude": float(0.5 * (df["q_cross"].max() - df["q_cross"].min())),
        "q_total_amplitude": float(0.5 * (df["q_amp"].max() - df["q_amp"].min())),
        "corr_input_with_q_plus": float(corr_plus),
        "corr_input_with_q_cross": float(corr_cross),
        "mean_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].mean()),
        "median_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].median()),
        "interpretation": "fixed-gauge extraction: plus -> q_plus; cross -> q_cross.",
    }

    pd.DataFrame([summary]).to_csv(outdir / "summary.csv", index=False)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    plt.figure(figsize=(7, 4))
    plt.plot(df["frame"], df["h_input"], label="input h(t)")
    plt.plot(df["frame"], df["q_plus"], label="q_plus")
    plt.plot(df["frame"], df["q_cross"], label="q_cross")
    plt.xlabel("frame")
    plt.ylabel("strain / quadrupole")
    plt.title(f"BuP fixed-gauge entanglement-wave response ({args.polarization})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_strain_timeseries.png", dpi=170)
    plt.close()

    first = frames / "frame_0000.png"
    last = frames / f"frame_{args.n_frames-1:04d}.png"
    if first.exists():
        (outdir / "fig_ring_initial.png").write_bytes(first.read_bytes())
    if last.exists():
        (outdir / "fig_ring_final.png").write_bytes(last.read_bytes())

    md = [
        "# Paper 21 — BuP gravitational-wave prototype v2",
        "",
        "v2 uses a fixed-gauge quadrupole extraction on entanglement-distance changes.",
        "",
        r"\[",
        r"\frac{\delta d_{ij}^2}{d_{ij}^2} \simeq q_+\cos(2\theta_{ij}) + q_\times\sin(2\theta_{ij}).",
        r"\]",
        "",
        "## Summary",
        "",
    ]
    for k, v in summary.items():
        md.append(f"- `{k}`: `{v}`")
    (outdir / "summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Paper 21 — BuP gravitational-wave prototype v2 fixed gauge")
    print("=" * 100)
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nFiles written:")
    for name in [
        "summary.csv",
        "summary.json",
        "summary.md",
        "ring_deformation_timeseries.csv",
        "fig_strain_timeseries.png",
        "fig_ring_initial.png",
        "fig_ring_final.png",
        "frames/",
    ]:
        print(outdir / name)


if __name__ == "__main__":
    main()

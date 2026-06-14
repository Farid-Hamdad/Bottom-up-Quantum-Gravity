#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 21 — BuP gravitational-wave prototype v3

Goal
----
v2 showed that the two transverse polarizations can be extracted from a
fixed-gauge perturbation of the entanglement graph.

v3 tests propagation.

A BuP gravitational wave is modeled as a traveling collective perturbation
of the entanglement graph:

    W_ij(t) = W_ij^0 + delta W_ij(x,t)

with

    h(x,t) = A cos(k x - omega t).

The response is measured on several test rings placed at different x positions.
For plus polarization:

    q_+(x,t) ~ cos(k x - omega t)

For cross polarization:

    q_x(x,t) ~ cos(k x - omega t)

Outputs
-------
summary.csv
summary.json
ring_timeseries.csv
fig_mode_timeseries_by_ring.png
fig_phase_fit.png
fig_spacetime_mode.png
frames/

Run
---
cd ~/bottomup

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_entanglement_wave_v3_propagation.py \
  --polarization plus \
  --N-side 41 \
  --amplitude 0.03 \
  --wavelength 1.6 \
  --n-frames 100 \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_plus_v3_propagation

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_entanglement_wave_v3_propagation.py \
  --polarization cross \
  --N-side 41 \
  --amplitude 0.03 \
  --wavelength 1.6 \
  --n-frames 100 \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_cross_v3_propagation
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


def pairwise_geometry(coords):
    diff = coords[:, None, :] - coords[None, :, :]
    r = np.sqrt(np.sum(diff * diff, axis=-1))
    nvec = np.zeros_like(diff)
    mask = r > 0
    nvec[mask] = diff[mask] / r[mask, None]
    mid = 0.5 * (coords[:, None, :] + coords[None, :, :])
    return diff, r, nvec, mid


def base_weights(coords, ell, cutoff):
    _, r, _, _ = pairwise_geometry(coords)
    W = np.exp(-r / ell)
    W[r > cutoff] = 0.0
    np.fill_diagonal(W, 0.0)
    return W


def traveling_strain(coords, polarization, amplitude, k_wave, omega, t, envelope_sigma):
    """
    Node-local traveling wave H_ab(x,t).
    The pair strain uses the midpoint x coordinate.
    This function returns pairwise H evaluated at each edge midpoint.
    """
    _, _, _, mid = pairwise_geometry(coords)
    xmid = mid[..., 0]
    ymid = mid[..., 1]
    envelope = np.exp(-(ymid ** 2) / (2 * envelope_sigma ** 2))
    h = amplitude * envelope * np.cos(k_wave * xmid - omega * t)

    H = np.zeros(xmid.shape + (2, 2), dtype=float)
    if polarization == "plus":
        H[..., 0, 0] = h
        H[..., 1, 1] = -h
    elif polarization == "cross":
        H[..., 0, 1] = h
        H[..., 1, 0] = h
    elif polarization == "mixed":
        hp = amplitude * envelope * np.cos(k_wave * xmid - omega * t)
        hx = amplitude * envelope * np.sin(k_wave * xmid - omega * t)
        H[..., 0, 0] = hp
        H[..., 1, 1] = -hp
        H[..., 0, 1] = hx
        H[..., 1, 0] = hx
    else:
        raise ValueError("polarization must be plus, cross, or mixed")
    return H


def perturb_weights_traveling(coords, W0, polarization, amplitude, k_wave, omega, t, envelope_sigma):
    _, r, nvec, _ = pairwise_geometry(coords)
    H = traveling_strain(coords, polarization, amplitude, k_wave, omega, t, envelope_sigma)
    q = np.einsum("...a,...ab,...b->...", nvec, H, nvec)
    W = W0 * np.exp(-0.5 * q)
    np.fill_diagonal(W, 0.0)
    return W


def entanglement_distance(W, l0, eps=1e-12):
    Wpos = W[W > 0]
    Wref = np.max(Wpos)
    D = np.full_like(W, np.nan, dtype=float)
    mask = W > 0
    ratio = np.clip((W[mask] + eps) / (Wref + eps), 1e-300, 1.0)
    D[mask] = -l0 * np.log(ratio)
    np.fill_diagonal(D, 0.0)
    return D


def select_ring(coords, center_x, center_y, radius, width):
    rel = coords - np.array([center_x, center_y])
    r = np.sqrt(np.sum(rel * rel, axis=1))
    idx = np.where(np.abs(r - radius) <= width)[0]
    theta = np.arctan2(rel[idx, 1], rel[idx, 0])
    return idx[np.argsort(theta)]


def ring_quadrupole_fixed_gauge(coords_ring, D0, Dt):
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


def fit_phase_and_amplitude(times, signal, omega):
    """
    Fit signal(t) = A cos(omega t) + B sin(omega t) + C.
    Return amplitude, phase phi in A0 cos(omega t - phi), offset, r2.
    """
    X = np.column_stack([np.cos(omega * times), np.sin(omega * times), np.ones_like(times)])
    coef, *_ = np.linalg.lstsq(X, signal, rcond=None)
    A, B, C = coef
    pred = X @ coef
    ss_res = np.sum((signal - pred) ** 2)
    ss_tot = np.sum((signal - signal.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    amp = np.sqrt(A * A + B * B)
    # A cos(wt) + B sin(wt) = amp cos(wt - phi), with phi=atan2(B,A)
    phi = np.arctan2(B, A)
    return float(amp), float(phi), float(C), float(r2)


def unwrap_phase(phases):
    return np.unwrap(np.asarray(phases))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-side", type=int, default=41)
    p.add_argument("--ell", type=float, default=0.14)
    p.add_argument("--cutoff", type=float, default=0.36)
    p.add_argument("--l0", type=float, default=0.14)
    p.add_argument("--amplitude", type=float, default=0.03)
    p.add_argument("--wavelength", type=float, default=1.6)
    p.add_argument("--period", type=float, default=1.0)
    p.add_argument("--n-frames", type=int, default=100)
    p.add_argument("--polarization", choices=["plus", "cross", "mixed"], default="plus")
    p.add_argument("--ring-radius", type=float, default=0.22)
    p.add_argument("--ring-width", type=float, default=0.035)
    p.add_argument("--ring-centers", nargs="+", type=float, default=[-0.6, -0.3, 0.0, 0.3, 0.6])
    p.add_argument("--envelope-sigma", type=float, default=0.8)
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

    k_wave = 2 * np.pi / args.wavelength
    omega = 2 * np.pi / args.period
    times = np.linspace(0.0, args.period, args.n_frames, endpoint=False)

    rings = []
    for cx in args.ring_centers:
        idx = select_ring(coords, cx, 0.0, args.ring_radius, args.ring_width)
        if len(idx) < 10:
            raise ValueError(f"Ring at x={cx} has too few nodes ({len(idx)}). Increase --ring-width.")
        rings.append((cx, idx, coords[idx], D0_all[np.ix_(idx, idx)]))

    rows = []

    for frame, t in enumerate(times):
        Wt = perturb_weights_traveling(
            coords, W0, args.polarization, args.amplitude, k_wave, omega, t, args.envelope_sigma
        )
        Dt_all = entanglement_distance(Wt, args.l0)

        for cx, idx, coords_ring, D0 in rings:
            Dt = Dt_all[np.ix_(idx, idx)]
            q_plus, q_cross, offset, fit_r2 = ring_quadrupole_fixed_gauge(coords_ring, D0, Dt)
            h_local = args.amplitude * np.cos(k_wave * cx - omega * t)
            rows.append({
                "frame": frame,
                "t": t,
                "center_x": cx,
                "h_local_expected": h_local,
                "q_plus": q_plus,
                "q_cross": q_cross,
                "q_amp": float(np.sqrt(q_plus*q_plus + q_cross*q_cross)),
                "offset": offset,
                "quadrupole_fit_r2": fit_r2,
                "n_ring": int(len(idx)),
            })

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "ring_timeseries.csv", index=False)

    # Select active mode.
    mode_col = "q_plus" if args.polarization == "plus" else "q_cross"
    if args.polarization == "mixed":
        mode_col = "q_amp"

    phase_rows = []
    for cx in args.ring_centers:
        sub = df[df["center_x"] == cx].sort_values("t")
        amp, phi, off, r2 = fit_phase_and_amplitude(sub["t"].values, sub[mode_col].values, omega)
        expected_phi = k_wave * cx
        phase_rows.append({
            "center_x": cx,
            "mode": mode_col,
            "amplitude_fit": amp,
            "phase_fit": phi,
            "phase_expected_modulo": np.arctan2(np.sin(expected_phi), np.cos(expected_phi)),
            "offset_fit": off,
            "temporal_fit_r2": r2,
            "corr_with_local_expected": float(np.corrcoef(sub["h_local_expected"], sub[mode_col])[0, 1]),
            "mean_quadrupole_fit_r2": float(sub["quadrupole_fit_r2"].mean()),
        })

    ph = pd.DataFrame(phase_rows)
    ph["phase_fit_unwrapped"] = unwrap_phase(ph["phase_fit"].values)
    ph.to_csv(outdir / "phase_fit.csv", index=False)

    # Fit phase vs x. For h=A cos(kx - wt), extracted q(t)=amp cos(wt - kx) -> phase ~= kx up to sign convention.
    X = np.column_stack([ph["center_x"].values, np.ones(len(ph))])
    y = ph["phase_fit_unwrapped"].values
    slope, intercept = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ np.array([slope, intercept])
    ss_res = np.sum((y - pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    phase_r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan

    # Global correlations.
    corr_active = np.corrcoef(df["h_local_expected"], df[mode_col])[0, 1]
    corr_plus = np.corrcoef(df["h_local_expected"], df["q_plus"])[0, 1]
    corr_cross = np.corrcoef(df["h_local_expected"], df["q_cross"])[0, 1]

    summary = {
        "N_side": args.N_side,
        "N_nodes": int(len(coords)),
        "polarization": args.polarization,
        "active_mode": mode_col,
        "input_amplitude": args.amplitude,
        "wavelength": args.wavelength,
        "k_wave_expected": float(k_wave),
        "omega": float(omega),
        "n_frames": args.n_frames,
        "n_rings": len(rings),
        "ring_centers": args.ring_centers,
        "global_corr_expected_with_active_mode": float(corr_active),
        "global_corr_expected_with_q_plus": float(corr_plus),
        "global_corr_expected_with_q_cross": float(corr_cross),
        "mean_temporal_fit_r2": float(ph["temporal_fit_r2"].mean()),
        "median_temporal_fit_r2": float(ph["temporal_fit_r2"].median()),
        "phase_slope_fit": float(slope),
        "phase_intercept_fit": float(intercept),
        "phase_vs_x_r2": float(phase_r2),
        "phase_slope_over_expected_k": float(slope / k_wave),
        "mean_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].mean()),
        "median_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].median()),
        "interpretation": "Propagation succeeds if active mode tracks local h(x,t) and phase_slope_over_expected_k is close to 1 in magnitude."
    }

    pd.DataFrame([summary]).to_csv(outdir / "summary.csv", index=False)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Figures.
    plt.figure(figsize=(8, 5))
    for cx in args.ring_centers:
        sub = df[df["center_x"] == cx]
        plt.plot(sub["t"], sub[mode_col], label=f"x={cx:g}")
    plt.xlabel("t")
    plt.ylabel(mode_col)
    plt.title(f"Traveling BuP entanglement wave: {mode_col} by ring")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "fig_mode_timeseries_by_ring.png", dpi=170)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.scatter(ph["center_x"], ph["phase_fit_unwrapped"], label="fit")
    plt.plot(ph["center_x"], pred, label=f"linear fit slope={slope:.3f}")
    plt.xlabel("ring center x")
    plt.ylabel("unwrapped phase")
    plt.title("Phase propagation fit")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_phase_fit.png", dpi=170)
    plt.close()

    # Spacetime heatmap.
    pivot = df.pivot(index="frame", columns="center_x", values=mode_col)
    plt.figure(figsize=(7, 5))
    plt.imshow(
        pivot.values,
        aspect="auto",
        origin="lower",
        extent=[min(args.ring_centers), max(args.ring_centers), 0, args.period],
    )
    plt.colorbar(label=mode_col)
    plt.xlabel("x")
    plt.ylabel("t")
    plt.title(f"Spacetime pattern of {mode_col}")
    plt.tight_layout()
    plt.savefig(outdir / "fig_spacetime_mode.png", dpi=170)
    plt.close()

    md = [
        "# Paper 21 — BuP gravitational-wave propagation v3",
        "",
        "A traveling graph perturbation is imposed:",
        "",
        r"\[",
        r"h(x,t)=A\cos(kx-\omega t).",
        r"\]",
        "",
        "The response is extracted from the active quadrupole mode.",
        "",
        "## Summary",
        "",
    ]
    for k, v in summary.items():
        md.append(f"- `{k}`: `{v}`")
    (outdir / "summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Paper 21 — BuP gravitational-wave prototype v3 propagation")
    print("=" * 100)
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nFiles written:")
    for name in [
        "summary.csv",
        "summary.json",
        "summary.md",
        "ring_timeseries.csv",
        "phase_fit.csv",
        "fig_mode_timeseries_by_ring.png",
        "fig_phase_fit.png",
        "fig_spacetime_mode.png",
    ]:
        print(outdir / name)


if __name__ == "__main__":
    main()

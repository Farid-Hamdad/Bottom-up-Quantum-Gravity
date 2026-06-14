#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 21 — BuP gravitational-wave prototype v4

Dynamic source on an entanglement graph + propagation speed estimate.

Goal
----
v3 imposed a traveling wave h(x,t)=A cos(kx-omega t).
v4 uses a localized oscillating source and measures the time delay of the
quadrupolar response on concentric rings.

In BuP language:

    localized source perturbation -> delta W_ij(t)
    delta W_ij(t) propagates on graph -> delta d_ent(i,j,t)
    ring response q_+(r,t), q_x(r,t) is measured
    speed v_graph is fitted from arrival time vs radius

This is still a kinematic graph-wave prototype, not a full binary-black-hole
merger. It tests whether a local dynamical perturbation of W_ij can emit an
outgoing quadrupolar wave and whether its propagation speed can be recovered.

Run
---
cd ~/bottomup

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_dynamic_source_v4.py \
  --polarization plus \
  --N-side 61 \
  --amplitude 0.035 \
  --speed 1.0 \
  --frequency 3.0 \
  --n-frames 180 \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_dynamic_source_plus_v4

python3 papers/paper21_bup_gravitational_waves/scripts/bup_gw_dynamic_source_v4.py \
  --polarization cross \
  --N-side 61 \
  --amplitude 0.035 \
  --speed 1.0 \
  --frequency 3.0 \
  --n-frames 180 \
  --output-dir papers/paper21_bup_gravitational_waves/results/gw_dynamic_source_cross_v4

Key outputs
-----------
summary.csv
summary.json
ring_timeseries.csv
arrival_fit.csv
fig_ring_response.png
fig_arrival_time_vs_radius.png
fig_spacetime_response.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Geometry and graph
# -----------------------------------------------------------------------------

def build_grid(n_side: int, extent: float):
    xs = np.linspace(-extent, extent, n_side)
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


def entanglement_distance(W, l0, eps=1e-12):
    Wpos = W[W > 0]
    Wref = np.max(Wpos)
    D = np.full_like(W, np.nan, dtype=float)
    mask = W > 0
    ratio = np.clip((W[mask] + eps) / (Wref + eps), 1e-300, 1.0)
    D[mask] = -l0 * np.log(ratio)
    np.fill_diagonal(D, 0.0)
    return D


# -----------------------------------------------------------------------------
# Dynamic source wave model
# -----------------------------------------------------------------------------

def source_wave_strain(coords, polarization, amplitude, speed, frequency, t, width, damping):
    """
    Outgoing quadrupolar wave from the origin.

    Retarded phase:
        phase = omega * (t - r/speed)

    Envelope:
        envelope = exp(-(r - speed*t)^2/(2 width^2)) * exp(-damping*r)

    Directional tensor:
      plus  in fixed axes: H = [[h,0],[0,-h]]
      cross in fixed axes: H = [[0,h],[h,0]]

    This is a controlled dynamic-source prototype: the wavefront is emitted
    from r=0 and propagates at the input speed.
    """
    _, _, _, mid = pairwise_geometry(coords)
    x = mid[..., 0]
    y = mid[..., 1]
    r = np.sqrt(x*x + y*y)

    omega = 2 * np.pi * frequency
    phase = omega * (t - r / speed)

    shell = np.exp(-((r - speed*t) ** 2) / (2 * width * width))
    damp = np.exp(-damping * r)
    h = amplitude * shell * damp * np.cos(phase)

    H = np.zeros(r.shape + (2, 2), dtype=float)
    if polarization == "plus":
        H[..., 0, 0] = h
        H[..., 1, 1] = -h
    elif polarization == "cross":
        H[..., 0, 1] = h
        H[..., 1, 0] = h
    elif polarization == "mixed":
        hp = amplitude * shell * damp * np.cos(phase)
        hx = amplitude * shell * damp * np.sin(phase)
        H[..., 0, 0] = hp
        H[..., 1, 1] = -hp
        H[..., 0, 1] = hx
        H[..., 1, 0] = hx
    else:
        raise ValueError("polarization must be plus, cross, or mixed")
    return H


def perturb_weights_dynamic(coords, W0, polarization, amplitude, speed, frequency, t, width, damping):
    _, _, nvec, _ = pairwise_geometry(coords)
    H = source_wave_strain(coords, polarization, amplitude, speed, frequency, t, width, damping)
    q = np.einsum("...a,...ab,...b->...", nvec, H, nvec)
    W = W0 * np.exp(-0.5 * q)
    np.fill_diagonal(W, 0.0)
    return W


# -----------------------------------------------------------------------------
# Ring extraction
# -----------------------------------------------------------------------------

def select_ring(coords, radius, width):
    r = np.sqrt(np.sum(coords * coords, axis=1))
    idx = np.where(np.abs(r - radius) <= width)[0]
    theta = np.arctan2(coords[idx, 1], coords[idx, 0])
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


def estimate_arrival_time(times, signal, threshold_fraction):
    """
    Arrival time = first time where abs(signal) exceeds threshold_fraction
    of max abs(signal). Also return peak time.
    """
    s = np.asarray(signal, dtype=float)
    a = np.abs(s)
    if not np.isfinite(a).any() or np.nanmax(a) <= 0:
        return np.nan, np.nan, np.nan

    max_amp = float(np.nanmax(a))
    threshold = threshold_fraction * max_amp
    above = np.where(a >= threshold)[0]
    arrival = float(times[above[0]]) if len(above) else np.nan
    peak = float(times[int(np.nanargmax(a))])
    return arrival, peak, max_amp


def fit_speed_from_arrivals(radius, arrival_time):
    r = np.asarray(radius, dtype=float)
    t = np.asarray(arrival_time, dtype=float)
    ok = np.isfinite(r) & np.isfinite(t)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan, np.nan

    X = np.column_stack([r[ok], np.ones(ok.sum())])
    slope, intercept = np.linalg.lstsq(X, t[ok], rcond=None)[0]
    pred = X @ np.array([slope, intercept])
    ss_res = np.sum((t[ok] - pred) ** 2)
    ss_tot = np.sum((t[ok] - t[ok].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    speed = 1.0 / slope if slope != 0 else np.nan
    return float(speed), float(slope), float(intercept), float(r2)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-side", type=int, default=61)
    p.add_argument("--extent", type=float, default=1.6)
    p.add_argument("--ell", type=float, default=0.12)
    p.add_argument("--cutoff", type=float, default=0.32)
    p.add_argument("--l0", type=float, default=0.12)

    p.add_argument("--polarization", choices=["plus", "cross", "mixed"], default="plus")
    p.add_argument("--amplitude", type=float, default=0.035)
    p.add_argument("--speed", type=float, default=1.0)
    p.add_argument("--frequency", type=float, default=3.0)
    p.add_argument("--width", type=float, default=0.16)
    p.add_argument("--damping", type=float, default=0.15)

    p.add_argument("--t-max", type=float, default=1.45)
    p.add_argument("--n-frames", type=int, default=180)

    p.add_argument("--ring-radii", nargs="+", type=float, default=[0.35, 0.55, 0.75, 0.95, 1.15])
    p.add_argument("--ring-width", type=float, default=0.035)
    p.add_argument("--arrival-threshold", type=float, default=0.35)

    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    coords = build_grid(args.N_side, args.extent)
    W0 = base_weights(coords, args.ell, args.cutoff)
    D0_all = entanglement_distance(W0, args.l0)

    rings = []
    for rr in args.ring_radii:
        idx = select_ring(coords, rr, args.ring_width)
        if len(idx) < 12:
            raise ValueError(f"Ring r={rr} has too few nodes ({len(idx)}). Increase --ring-width.")
        rings.append((rr, idx, coords[idx], D0_all[np.ix_(idx, idx)]))

    times = np.linspace(0.0, args.t_max, args.n_frames)
    rows = []

    for frame, t in enumerate(times):
        Wt = perturb_weights_dynamic(
            coords, W0,
            args.polarization,
            args.amplitude,
            args.speed,
            args.frequency,
            t,
            args.width,
            args.damping,
        )
        Dt_all = entanglement_distance(Wt, args.l0)

        for rr, idx, coords_ring, D0 in rings:
            Dt = Dt_all[np.ix_(idx, idx)]
            q_plus, q_cross, offset, fit_r2 = ring_quadrupole_fixed_gauge(coords_ring, D0, Dt)
            active = q_plus if args.polarization == "plus" else q_cross
            if args.polarization == "mixed":
                active = np.sqrt(q_plus*q_plus + q_cross*q_cross)
            rows.append({
                "frame": frame,
                "t": t,
                "radius": rr,
                "q_plus": q_plus,
                "q_cross": q_cross,
                "q_active": active,
                "q_amp": float(np.sqrt(q_plus*q_plus + q_cross*q_cross)),
                "offset": offset,
                "quadrupole_fit_r2": fit_r2,
                "n_ring": int(len(idx)),
            })

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "ring_timeseries.csv", index=False)

    arrival_rows = []
    for rr in args.ring_radii:
        sub = df[df["radius"] == rr].sort_values("t")
        arrival, peak, amp = estimate_arrival_time(sub["t"].values, sub["q_active"].values, args.arrival_threshold)
        arrival_rows.append({
            "radius": rr,
            "arrival_time": arrival,
            "peak_time": peak,
            "max_abs_active_mode": amp,
            "mean_quadrupole_fit_r2": float(sub["quadrupole_fit_r2"].mean()),
            "n_ring": int(sub["n_ring"].iloc[0]),
        })

    arr = pd.DataFrame(arrival_rows)
    arr.to_csv(outdir / "arrival_fit.csv", index=False)

    speed_arrival, slope_arrival, intercept_arrival, r2_arrival = fit_speed_from_arrivals(
        arr["radius"].values, arr["arrival_time"].values
    )
    speed_peak, slope_peak, intercept_peak, r2_peak = fit_speed_from_arrivals(
        arr["radius"].values, arr["peak_time"].values
    )

    summary = {
        "N_side": args.N_side,
        "N_nodes": int(len(coords)),
        "polarization": args.polarization,
        "input_amplitude": args.amplitude,
        "input_speed": args.speed,
        "frequency": args.frequency,
        "width": args.width,
        "damping": args.damping,
        "n_frames": args.n_frames,
        "t_max": args.t_max,
        "ring_radii": args.ring_radii,
        "arrival_threshold_fraction": args.arrival_threshold,
        "speed_from_arrival": float(speed_arrival),
        "speed_from_arrival_over_input": float(speed_arrival / args.speed) if np.isfinite(speed_arrival) else np.nan,
        "arrival_fit_r2": float(r2_arrival),
        "arrival_slope_dt_dr": float(slope_arrival),
        "arrival_intercept": float(intercept_arrival),
        "speed_from_peak": float(speed_peak),
        "speed_from_peak_over_input": float(speed_peak / args.speed) if np.isfinite(speed_peak) else np.nan,
        "peak_fit_r2": float(r2_peak),
        "peak_slope_dt_dr": float(slope_peak),
        "peak_intercept": float(intercept_peak),
        "mean_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].mean()),
        "median_quadrupole_fit_r2": float(df["quadrupole_fit_r2"].median()),
        "interpretation": "Propagation speed succeeds if speed_from_arrival_over_input or speed_from_peak_over_input is close to 1 and fit R2 is high."
    }

    pd.DataFrame([summary]).to_csv(outdir / "summary.csv", index=False)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Figures
    plt.figure(figsize=(8, 5))
    for rr in args.ring_radii:
        sub = df[df["radius"] == rr]
        plt.plot(sub["t"], sub["q_active"], label=f"r={rr:g}")
    plt.xlabel("t")
    plt.ylabel("active quadrupole")
    plt.title(f"Dynamic source response ({args.polarization})")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "fig_ring_response.png", dpi=170)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.scatter(arr["radius"], arr["arrival_time"], label="arrival")
    if np.isfinite(slope_arrival):
        rr = np.linspace(min(args.ring_radii), max(args.ring_radii), 100)
        plt.plot(rr, slope_arrival * rr + intercept_arrival, label=f"arrival fit v={speed_arrival:.3f}")
    plt.scatter(arr["radius"], arr["peak_time"], label="peak")
    if np.isfinite(slope_peak):
        rr = np.linspace(min(args.ring_radii), max(args.ring_radii), 100)
        plt.plot(rr, slope_peak * rr + intercept_peak, label=f"peak fit v={speed_peak:.3f}")
    plt.xlabel("radius")
    plt.ylabel("time")
    plt.title("Propagation speed from ring delays")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_arrival_time_vs_radius.png", dpi=170)
    plt.close()

    pivot = df.pivot(index="frame", columns="radius", values="q_active")
    plt.figure(figsize=(7, 5))
    plt.imshow(
        pivot.values,
        aspect="auto",
        origin="lower",
        extent=[min(args.ring_radii), max(args.ring_radii), times[0], times[-1]],
    )
    plt.colorbar(label="active quadrupole")
    plt.xlabel("radius")
    plt.ylabel("t")
    plt.title("Radial spacetime response")
    plt.tight_layout()
    plt.savefig(outdir / "fig_spacetime_response.png", dpi=170)
    plt.close()

    md = [
        "# Paper 21 — BuP dynamic-source gravitational wave v4",
        "",
        "A localized oscillating source emits an outgoing perturbation of the entanglement graph.",
        "",
        r"\[",
        r"W_{ij}(t)=W_{ij}^{(0)}+\delta W_{ij}^{\rm source}(t).",
        r"\]",
        "",
        "The propagation speed is inferred from ring arrival times.",
        "",
        "## Summary",
        "",
    ]
    for k, v in summary.items():
        md.append(f"- `{k}`: `{v}`")
    (outdir / "summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Paper 21 — BuP gravitational-wave prototype v4 dynamic source")
    print("=" * 100)
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nFiles written:")
    for name in [
        "summary.csv",
        "summary.json",
        "summary.md",
        "ring_timeseries.csv",
        "arrival_fit.csv",
        "fig_ring_response.png",
        "fig_arrival_time_vs_radius.png",
        "fig_spacetime_response.png",
    ]:
        print(outdir / name)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP gravitational-wave dynamic generation v10
True edge-Hessian dynamics with spectral-band diagnostics.

This version moves beyond the nodal reduction phi_i. The dynamical variable is
an edge perturbation

    x_e(t) = delta W_e(t) / W_e^(0),     e=(i,j)

and the equation solved is

    x_ddot_e + gamma x_dot_e + sum_f K_ef x_f = J_e(t).

The complete edge Hessian K_ef is constructed analytically from an explicit
edge-level BuP rigidity functional

    F_edge[x] = 1/2 eta_edge sum_{e~f} A_ef (x_e - x_f)^2
              + 1/2 sum_e M_e^2 x_e^2,

where e~f are neighbouring edges in the line graph of the entanglement network.
Thus

    K = eta_edge * L_line + diag(M_e^2).

This is a true Hessian in the edge variables x_e. It is not the expensive full
finite-difference Hessian of the v5 spectral action. That fully spectral edge
Hessian is a later, heavier variant.

Outputs:
  results/summary_true_edge_hessian_v10.json
  results/summary_true_edge_hessian_v10.csv
  results/v10_edge_modes.csv
  results/v10_modal_response_plus.csv
  results/v10_modal_response_cross.csv
  results/v10_packet_peaks_plus.csv
  results/v10_packet_peaks_cross.csv
  results/v10_edge_hessian_sparse.npz
  figures/fig_v10_dispersion.png
  figures/fig_v10_modal_response_plus.png
  figures/fig_v10_modal_response_cross.png
  figures/fig_v10_spacetime_plus.png
  figures/fig_v10_spacetime_cross.png
  figures/fig_v10_packet_fit_plus.png
  figures/fig_v10_packet_fit_cross.png

Example:
python3 bup_gw_true_edge_hessian_v10.py \
  --N-side 11 --eta-edge 1.0 --mass2 1e-4 --gamma 0.015 \
  --source-amp 0.20 --source-width 0.35 --pulse-t0 8.0 --pulse-sigma 0.9 \
  --n-modes-dyn 180 --output-dir results/true_edge_hessian_v9
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh


def build_grid(n_side: int, extent: float) -> np.ndarray:
    xs = np.linspace(-extent, extent, n_side)
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def pairwise_dist(coords: np.ndarray) -> np.ndarray:
    d = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(d * d, axis=-1))


def build_edges(coords: np.ndarray, ell: float, cutoff: float):
    R = pairwise_dist(coords)
    W = np.exp(-R / ell)
    W[R > cutoff] = 0.0
    np.fill_diagonal(W, 0.0)
    edges = []
    w0 = []
    lengths = []
    n = coords.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] > 0:
                edges.append((i, j))
                w0.append(W[i, j])
                lengths.append(R[i, j])
    return np.asarray(edges, dtype=int), np.asarray(w0, dtype=float), np.asarray(lengths, dtype=float), W, R


def edge_geometry(coords: np.ndarray, edges: np.ndarray):
    a = coords[edges[:, 0]]
    b = coords[edges[:, 1]]
    centers = 0.5 * (a + b)
    vecs = b - a
    theta = np.arctan2(vecs[:, 1], vecs[:, 0])
    radii = np.sqrt(np.sum(centers * centers, axis=1))
    return centers, theta, radii


def build_line_graph_hessian(
    coords: np.ndarray,
    edges: np.ndarray,
    w0: np.ndarray,
    lengths: np.ndarray,
    centers: np.ndarray,
    eta_edge: float,
    mass2: float,
    lambda_locality: float,
    edge_width: float,
    w_barrier: float,
):
    """Analytic complete edge Hessian K = eta L_line + diag(M_e^2)."""
    E = len(edges)
    node_to_edges = defaultdict(list)
    for e, (i, j) in enumerate(edges):
        node_to_edges[int(i)].append(e)
        node_to_edges[int(j)].append(e)

    weights = defaultdict(float)
    for _, inc in node_to_edges.items():
        m = len(inc)
        for a in range(m):
            e = inc[a]
            for b in range(a + 1, m):
                f = inc[b]
                dc = np.linalg.norm(centers[e] - centers[f])
                wij = np.exp(-(dc * dc) / (2.0 * edge_width * edge_width)) if edge_width > 0 else 1.0
                key = (e, f) if e < f else (f, e)
                weights[key] += float(wij)

    rows = []
    cols = []
    data = []
    degree = np.zeros(E, dtype=float)
    for (e, f), a in weights.items():
        degree[e] += a
        degree[f] += a
        rows.extend([e, f])
        cols.extend([f, e])
        data.extend([-eta_edge * a, -eta_edge * a])

    diag_mass = mass2 + lambda_locality * lengths * lengths + w_barrier / np.maximum(w0, 1e-12) ** 2
    diag = eta_edge * degree + diag_mass
    rows.extend(np.arange(E).tolist())
    cols.extend(np.arange(E).tolist())
    data.extend(diag.tolist())
    K = sparse.csr_matrix((data, (rows, cols)), shape=(E, E))
    K = 0.5 * (K + K.T)

    Lline = sparse.csr_matrix((E, E))
    if weights:
        lr, lc, ld = [], [], []
        deg = np.zeros(E)
        for (e, f), a in weights.items():
            deg[e] += a
            deg[f] += a
            lr.extend([e, f])
            lc.extend([f, e])
            ld.extend([-a, -a])
        lr.extend(np.arange(E).tolist())
        lc.extend(np.arange(E).tolist())
        ld.extend(deg.tolist())
        Lline = sparse.csr_matrix((ld, (lr, lc)), shape=(E, E))
        Lline = 0.5 * (Lline + Lline.T)
    return K, Lline, diag_mass


def low_eigensystem(K: sparse.csr_matrix, n_modes: int):
    E = K.shape[0]
    k = min(n_modes, E - 2)
    if k >= E - 2 or E <= 250:
        vals, vecs = eigh(K.toarray())
        return vals[:n_modes], vecs[:, :n_modes]
    vals, vecs = eigsh(K, k=k, which="SM", tol=1e-9, maxiter=20000)
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx]
    return vals, vecs


def mode_energy(mode: np.ndarray, Op: sparse.csr_matrix) -> float:
    den = float(mode @ mode)
    if den <= 0:
        return np.nan
    return float(mode @ (Op @ mode) / den)


def fit_line(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan
    X = np.column_stack([x[ok], np.ones(ok.sum())])
    slope, intercept = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    pred = X @ np.array([slope, intercept])
    ss_res = float(np.sum((y[ok] - pred) ** 2))
    ss_tot = float(np.sum((y[ok] - y[ok].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(slope), float(intercept), float(r2)


def ricker_pulse(t, t0, sigma):
    z = (t - t0) / sigma
    return (1.0 - z * z) * np.exp(-0.5 * z * z)


def gaussian_pulse(t, t0, sigma):
    z = (t - t0) / sigma
    return np.exp(-0.5 * z * z)


def source_spatial(centers, theta, radii, pol, width):
    # Edge-level TT-like quadrupolar source on the edge centers/orientations.
    env = np.exp(-(radii * radii) / (2.0 * width * width))
    if pol == "plus":
        q = np.cos(2.0 * theta)
    elif pol == "cross":
        q = np.sin(2.0 * theta)
    else:
        raise ValueError(pol)
    s = env * q
    nrm = np.sqrt(np.mean(s * s))
    return s / nrm if nrm > 0 else s


def make_rings(radii, nrings, rmin, rmax):
    bins = np.linspace(rmin, rmax, nrings + 1)
    rings = []
    ring_r = []
    for a, b in zip(bins[:-1], bins[1:]):
        mask = (radii >= a) & (radii < b)
        if mask.sum() >= 3:
            rings.append(mask)
            ring_r.append(float(0.5 * (a + b)))
    return rings, np.asarray(ring_r)


def project_rings(x_edge, rings, theta, pol):
    basis = np.cos(2.0 * theta) if pol == "plus" else np.sin(2.0 * theta)
    vals = []
    for mask in rings:
        b = basis[mask]
        den = float(np.sum(b * b))
        vals.append(float(np.sum(x_edge[mask] * b) / den) if den > 1e-12 else np.nan)
    return np.asarray(vals, dtype=float)


def modal_integrate(vals, vecs, J_edge, args, rings, theta, pol):
    omega2 = np.maximum(vals, 0.0)
    Jn = vecs.T @ J_edge
    q = np.zeros_like(vals)
    qdot = np.zeros_like(vals)
    q_peak = np.zeros_like(vals)

    times = np.arange(0.0, args.tmax + 0.5 * args.dt, args.dt)
    rec_t = []
    rec_ring = []

    for step, t in enumerate(times):
        pulse = ricker_pulse(t, args.pulse_t0, args.pulse_sigma) if args.pulse_kind == "ricker" else gaussian_pulse(t, args.pulse_t0, args.pulse_sigma)
        force = args.source_amp * pulse * Jn
        qdd = force - args.gamma * qdot - omega2 * q
        qdot += args.dt * qdd
        q += args.dt * qdot
        q_peak = np.maximum(q_peak, np.abs(q))

        if step % args.record_stride == 0:
            x_edge = vecs @ q
            rec_t.append(float(t))
            rec_ring.append(project_rings(x_edge, rings, theta, pol))

    return np.asarray(times), Jn, q_peak, np.asarray(rec_t), np.asarray(rec_ring)


def packet_fit(times, signals, ring_r, t_ignore, baseline_until, min_snr):
    rows = []
    for ir, r in enumerate(ring_r):
        y = signals[:, ir]
        okb = times <= baseline_until
        noise = float(np.std(y[okb])) if okb.sum() >= 3 else float(np.std(y[:max(3, len(y)//10)]))
        noise = max(noise, 1e-12)
        oka = times >= t_ignore
        if oka.sum() < 3:
            continue
        yy = np.abs(y[oka])
        tt = times[oka]
        im = int(np.argmax(yy))
        amp = float(yy[im])
        snr = amp / noise
        rows.append({"ring": ir, "r": float(r), "t_peak": float(tt[im]), "amp": amp, "snr": float(snr), "valid": bool(snr >= min_snr)})
    df = pd.DataFrame(rows)
    if len(df) == 0 or int(df["valid"].sum()) < 3:
        return df, np.nan, np.nan, 0
    d = df[df["valid"]].copy()
    slope, intercept, r2 = fit_line(d["r"], d["t_peak"])
    v = 1.0 / slope if np.isfinite(slope) and slope > 0 else np.nan
    return df, float(v), float(r2), int(len(d))


def plot_dispersion(df, slope, intercept, r2, outpath):
    plt.figure(figsize=(6.2, 4.4))
    plt.scatter(df["k2_line"], df["omega2"], s=22, label="edge modes")
    if np.isfinite(slope):
        xs = np.linspace(float(df["k2_line"].min()), float(df["k2_line"].max()), 150)
        plt.plot(xs, slope * xs + intercept, label=f"fit: c²={slope:.3g}, R²={r2:.4f}")
    plt.xlabel(r"edge-line $k^2$")
    plt.ylabel(r"edge Hessian $\omega^2$")
    plt.title("True edge-Hessian dispersion")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def plot_modal(df, outpath, title):
    plt.figure(figsize=(5.8, 4.4))
    plt.scatter(np.abs(df["J_n"]), df["q_peak"], s=24)
    plt.xlabel(r"$|J_n|$")
    plt.ylabel(r"peak $|q_n|$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()




def spectral_band_table(mdf: pd.DataFrame, n_bins: int, low_mode_cut: int):
    d = mdf.copy()
    d = d[np.isfinite(d["k2_line"]) & np.isfinite(d["q_peak"])].copy()
    d["k_line"] = np.sqrt(np.maximum(d["k2_line"].values, 0.0))
    d["energy"] = d["q_peak"].values ** 2
    if len(d) == 0 or float(d["energy"].sum()) <= 0:
        return pd.DataFrame(), {}

    kmin = float(d["k_line"].min())
    kmax = float(d["k_line"].max())
    if kmax <= kmin:
        bins = np.linspace(kmin, kmin + 1.0, n_bins + 1)
    else:
        bins = np.linspace(kmin, kmax, n_bins + 1)
    labels = [f"band_{i:02d}" for i in range(n_bins)]
    d["band"] = pd.cut(d["k_line"], bins=bins, labels=labels, include_lowest=True)

    total_E = float(d["energy"].sum())
    rows = []
    for i, lab in enumerate(labels):
        sub = d[d["band"] == lab]
        e = float(sub["energy"].sum()) if len(sub) else 0.0
        rows.append({
            "band": lab,
            "k_min": float(bins[i]),
            "k_max": float(bins[i+1]),
            "n_modes": int(len(sub)),
            "energy": e,
            "energy_fraction": e / total_E if total_E > 0 else np.nan,
            "mean_abs_J_n": float(sub["abs_J_n"].mean()) if len(sub) else np.nan,
            "mean_q_peak": float(sub["q_peak"].mean()) if len(sub) else np.nan,
        })
    band_df = pd.DataFrame(rows)

    low = d[d["mode"] < low_mode_cut]
    mid = d[(d["mode"] >= low_mode_cut) & (d["mode"] < max(low_mode_cut, len(d)//2))]
    high = d[d["mode"] >= max(low_mode_cut, len(d)//2)]
    stats = {
        "total_modal_energy": total_E,
        "low_mode_cut": int(low_mode_cut),
        "low_mode_energy_fraction": float(low["energy"].sum() / total_E) if total_E > 0 else np.nan,
        "mid_mode_energy_fraction": float(mid["energy"].sum() / total_E) if total_E > 0 else np.nan,
        "high_mode_energy_fraction": float(high["energy"].sum() / total_E) if total_E > 0 else np.nan,
        "dominant_mode": int(d.loc[d["energy"].idxmax(), "mode"]),
        "dominant_mode_k": float(d.loc[d["energy"].idxmax(), "k_line"]),
        "dominant_mode_energy_fraction": float(d["energy"].max() / total_E) if total_E > 0 else np.nan,
        "energy_weighted_k": float(np.sum(d["energy"] * d["k_line"]) / total_E) if total_E > 0 else np.nan,
    }
    return band_df, stats


def plot_qpeak_vs_k(mdf: pd.DataFrame, outpath, title):
    d = mdf.copy()
    d = d[np.isfinite(d["k2_line"])].copy()
    d["k_line"] = np.sqrt(np.maximum(d["k2_line"].values, 0.0))
    plt.figure(figsize=(6.0, 4.4))
    plt.scatter(d["k_line"], d["q_peak"], s=24)
    plt.xlabel(r"edge-line $k$")
    plt.ylabel(r"peak $|q_n|$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def plot_energy_bands(band_df: pd.DataFrame, outpath, title):
    plt.figure(figsize=(7.2, 4.4))
    if len(band_df):
        x = np.arange(len(band_df))
        plt.bar(x, band_df["energy_fraction"].values)
        labs = [f"{row.k_min:.2f}-{row.k_max:.2f}" for row in band_df.itertuples()]
        plt.xticks(x, labs, rotation=35, ha="right")
    plt.xlabel(r"$k$ band")
    plt.ylabel("energy fraction")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def plot_spacetime(times, signals, outpath, title):
    plt.figure(figsize=(7.0, 4.5))
    plt.imshow(signals.T, aspect="auto", origin="lower", extent=[times.min(), times.max(), 0, signals.shape[1]-1])
    plt.colorbar(label="ring projection")
    plt.xlabel("t")
    plt.ylabel("ring index")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def plot_packet(df, v, r2, outpath, title):
    plt.figure(figsize=(5.8, 4.4))
    if len(df) > 0:
        valid = df["valid"].astype(bool)
        plt.scatter(df.loc[~valid, "r"], df.loc[~valid, "t_peak"], label="invalid", alpha=0.5)
        plt.scatter(df.loc[valid, "r"], df.loc[valid, "t_peak"], label="valid")
        if valid.sum() >= 3 and np.isfinite(v):
            d = df[valid]
            slope, intercept, _ = fit_line(d["r"], d["t_peak"])
            xs = np.linspace(float(d["r"].min()), float(d["r"].max()), 100)
            plt.plot(xs, slope * xs + intercept, label=f"v={v:.3g}, R²={r2:.3f}")
    plt.xlabel("ring radius")
    plt.ylabel("peak time")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-side", type=int, default=11)
    p.add_argument("--extent", type=float, default=1.0)
    p.add_argument("--ell", type=float, default=0.22)
    p.add_argument("--cutoff", type=float, default=0.45)
    p.add_argument("--eta-edge", type=float, default=1.0)
    p.add_argument("--mass2", type=float, default=1e-4)
    p.add_argument("--lambda-locality", type=float, default=0.02)
    p.add_argument("--edge-width", type=float, default=0.30)
    p.add_argument("--w-barrier", type=float, default=0.0)
    p.add_argument("--gamma", type=float, default=0.015)
    p.add_argument("--source-amp", type=float, default=0.20)
    p.add_argument("--source-width", type=float, default=0.35)
    p.add_argument("--pulse-t0", type=float, default=8.0)
    p.add_argument("--pulse-sigma", type=float, default=0.9)
    p.add_argument("--pulse-kind", choices=["ricker", "gaussian"], default="ricker")
    p.add_argument("--tmax", type=float, default=45.0)
    p.add_argument("--dt", type=float, default=0.01)
    p.add_argument("--record-stride", type=int, default=5)
    p.add_argument("--n-modes-dyn", type=int, default=180)
    p.add_argument("--n-fit-modes", type=int, default=60)
    p.add_argument("--n-spectral-bins", type=int, default=10)
    p.add_argument("--low-mode-cut", type=int, default=25)
    p.add_argument("--nrings", type=int, default=10)
    p.add_argument("--ring-rmin", type=float, default=0.15)
    p.add_argument("--rmax-frac", type=float, default=0.82)
    p.add_argument("--t-ignore", type=float, default=7.5)
    p.add_argument("--baseline-until", type=float, default=5.0)
    p.add_argument("--min-snr", type=float, default=1.25)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    resdir = outdir / "results"
    figdir = outdir / "figures"
    resdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    print("Building graph...")
    coords = build_grid(args.N_side, args.extent)
    edges, w0, lengths, W0, R = build_edges(coords, args.ell, args.cutoff)
    centers, theta, radii = edge_geometry(coords, edges)
    E = len(edges)
    print(f"Graph: N={coords.shape[0]}, E={E}")

    print("Building true edge Hessian K_ef...")
    K, Lline, diag_mass = build_line_graph_hessian(
        coords, edges, w0, lengths, centers,
        eta_edge=args.eta_edge,
        mass2=args.mass2,
        lambda_locality=args.lambda_locality,
        edge_width=args.edge_width,
        w_barrier=args.w_barrier,
    )
    sparse.save_npz(resdir / "v10_edge_hessian_sparse.npz", K)

    print("Diagonalizing low edge-Hessian modes...")
    vals, vecs = low_eigensystem(K, args.n_modes_dyn)
    rows = []
    for n, val in enumerate(vals):
        mode = vecs[:, n]
        k2 = mode_energy(mode, Lline)
        rows.append({
            "mode": n,
            "omega2": float(val),
            "omega": float(np.sqrt(max(val, 0.0))),
            "k2_line": float(k2),
            "k_line": float(np.sqrt(max(k2, 0.0))) if np.isfinite(k2) else np.nan,
        })
    modes_df = pd.DataFrame(rows)
    fit_df = modes_df[(modes_df["mode"] >= 1) & (modes_df["omega2"] > 0)].head(args.n_fit_modes).copy()
    slope, intercept, r2 = fit_line(fit_df["k2_line"], fit_df["omega2"])
    c_graph = float(np.sqrt(slope)) if slope > 0 else np.nan
    modes_df.to_csv(resdir / "v10_edge_modes.csv", index=False)
    fit_df.to_csv(resdir / "v10_edge_dispersion_fit_modes.csv", index=False)
    plot_dispersion(fit_df, slope, intercept, r2, figdir / "fig_v10_dispersion.png")
    print(f"Dispersion: c={c_graph:.6f}, meff2={intercept:.6f}, R2={r2:.6f}")

    rmax = args.rmax_frac * float(np.max(radii))
    rings, ring_r = make_rings(radii, args.nrings, args.ring_rmin, rmax)

    summaries = {}
    for pol in ["plus", "cross"]:
        print(f"Integrating {pol} true-edge source...")
        Jedge = source_spatial(centers, theta, radii, pol=pol, width=args.source_width)
        times_full, Jn, q_peak, rec_t, rec_ring = modal_integrate(vals, vecs, Jedge, args, rings, theta, pol)
        mdf = pd.DataFrame({
            "mode": np.arange(len(vals)),
            "omega2": vals,
            "J_n": Jn,
            "abs_J_n": np.abs(Jn),
            "q_peak": q_peak,
            "k2_line": modes_df["k2_line"].values[:len(vals)],
        })
        corr = float(np.corrcoef(np.abs(Jn), q_peak)[0, 1]) if len(vals) >= 3 else np.nan
        low25 = min(25, len(vals))
        energy = q_peak * q_peak
        low25_energy = float(np.sum(energy[:low25]) / np.sum(energy)) if np.sum(energy) > 0 else np.nan
        mdf.to_csv(resdir / f"v10_modal_response_{pol}.csv", index=False)
        band_df, band_stats = spectral_band_table(mdf, args.n_spectral_bins, args.low_mode_cut)
        band_df.to_csv(resdir / f"v10_spectral_bands_{pol}.csv", index=False)
        pd.DataFrame(rec_ring, columns=[f"ring_{i}" for i in range(rec_ring.shape[1])]).assign(t=rec_t).to_csv(resdir / f"v10_ring_signals_{pol}.csv", index=False)
        packet_df, vpack, r2pack, valid = packet_fit(rec_t, rec_ring, ring_r, args.t_ignore, args.baseline_until, args.min_snr)
        packet_df.to_csv(resdir / f"v10_packet_peaks_{pol}.csv", index=False)
        plot_modal(mdf, figdir / f"fig_v10_modal_response_{pol}.png", f"{pol}: |Jn| vs peak modal response")
        plot_qpeak_vs_k(mdf, figdir / f"fig_v10_qpeak_vs_k_{pol}.png", f"{pol}: peak response vs k")
        plot_energy_bands(band_df, figdir / f"fig_v10_energy_bands_{pol}.png", f"{pol}: excited energy by k band")
        plot_spacetime(rec_t, rec_ring, figdir / f"fig_v10_spacetime_{pol}.png", f"{pol}: edge ring signal")
        plot_packet(packet_df, vpack, r2pack, figdir / f"fig_v10_packet_fit_{pol}.png", f"{pol}: packet fit")
        print(f"{pol.capitalize()} modal response: corr(|J_n|, q_peak)={corr:.6f}, low{args.low_mode_cut}_energy={band_stats.get('low_mode_energy_fraction', np.nan):.6f}")
        print(f"{pol.capitalize()} spectral: k_weighted={band_stats.get('energy_weighted_k', np.nan):.6f}, dominant_mode={band_stats.get('dominant_mode', -1)}, dominant_frac={band_stats.get('dominant_mode_energy_fraction', np.nan):.6f}")
        print(f"{pol.capitalize()} packet diagnostic: v={vpack:.6f}, R2={r2pack:.6f}, valid={valid}")
        summaries[pol] = {
            "modal_corr_absJ_qpeak": corr,
            "low25_energy": low25_energy,
            "packet_v": float(vpack) if np.isfinite(vpack) else None,
            "packet_r2": float(r2pack) if np.isfinite(r2pack) else None,
            "packet_valid": int(valid),
            "spectral_stats": band_stats,
        }

    summary = {
        "version": "v10_true_edge_hessian",
        "N_side": args.N_side,
        "N_nodes": int(coords.shape[0]),
        "N_edges": int(E),
        "eta_edge": args.eta_edge,
        "mass2": args.mass2,
        "lambda_locality": args.lambda_locality,
        "edge_width": args.edge_width,
        "w_barrier": args.w_barrier,
        "n_modes_dyn": int(len(vals)),
        "dispersion_c_graph": c_graph,
        "dispersion_slope_c2": float(slope),
        "dispersion_intercept_mass2": float(intercept),
        "dispersion_r2": float(r2),
        "plus": summaries["plus"],
        "cross": summaries["cross"],
        "interpretation": "True edge-variable Hessian dynamics with spectral diagnostics. The packet fit is secondary; the main tests are dispersion, modal forcing response, and excited k-band distribution.",
    }
    with open(resdir / "summary_true_edge_hessian_v10.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    pd.DataFrame([{
        "N_nodes": summary["N_nodes"],
        "N_edges": summary["N_edges"],
        "c_graph": summary["dispersion_c_graph"],
        "meff2": summary["dispersion_intercept_mass2"],
        "R2_disp": summary["dispersion_r2"],
        "plus_corr": summaries["plus"]["modal_corr_absJ_qpeak"],
        "cross_corr": summaries["cross"]["modal_corr_absJ_qpeak"],
        "plus_low25_energy": summaries["plus"]["low25_energy"],
        "cross_low25_energy": summaries["cross"]["low25_energy"],
        "plus_energy_weighted_k": summaries["plus"]["spectral_stats"].get("energy_weighted_k"),
        "cross_energy_weighted_k": summaries["cross"]["spectral_stats"].get("energy_weighted_k"),
        "plus_dominant_mode": summaries["plus"]["spectral_stats"].get("dominant_mode"),
        "cross_dominant_mode": summaries["cross"]["spectral_stats"].get("dominant_mode"),
        "plus_packet_v": summaries["plus"]["packet_v"],
        "plus_packet_R2": summaries["plus"]["packet_r2"],
        "plus_packet_valid": summaries["plus"]["packet_valid"],
        "cross_packet_v": summaries["cross"]["packet_v"],
        "cross_packet_R2": summaries["cross"]["packet_r2"],
        "cross_packet_valid": summaries["cross"]["packet_valid"],
    }]).to_csv(resdir / "summary_true_edge_hessian_v10.csv", index=False)

    md = [
        "# BuP GW dynamic generation v10 — true edge Hessian spectral diagnostics",
        "",
        "Dynamical variable: `x_e = delta W_e / W_e0` on each active edge.",
        "",
        "Equation:",
        "",
        r"\[\ddot x_e + \gamma \dot x_e + \sum_f K_{ef}x_f = J_e(t).\]",
        "",
        "Summary:",
        "",
        f"- N nodes: {summary['N_nodes']}",
        f"- E edges: {summary['N_edges']}",
        f"- c_graph: {c_graph}",
        f"- meff2: {intercept}",
        f"- R2 dispersion: {r2}",
        f"- plus corr(|J_n|, q_peak): {summaries['plus']['modal_corr_absJ_qpeak']}",
        f"- cross corr(|J_n|, q_peak): {summaries['cross']['modal_corr_absJ_qpeak']}",
        f"- plus energy-weighted k: {summaries['plus']['spectral_stats'].get('energy_weighted_k')}",
        f"- plus dominant mode: {summaries['plus']['spectral_stats'].get('dominant_mode')}",
        f"- plus packet v/R2/valid: {summaries['plus']['packet_v']} / {summaries['plus']['packet_r2']} / {summaries['plus']['packet_valid']}",
        f"- cross energy-weighted k: {summaries['cross']['spectral_stats'].get('energy_weighted_k')}",
        f"- cross dominant mode: {summaries['cross']['spectral_stats'].get('dominant_mode')}",
        f"- cross packet v/R2/valid: {summaries['cross']['packet_v']} / {summaries['cross']['packet_r2']} / {summaries['cross']['packet_valid']}",
    ]
    (resdir / "summary_true_edge_hessian_v10.md").write_text("\n".join(md), encoding="utf-8")

    print(f"Wrote outputs to: {outdir}")


if __name__ == "__main__":
    main()

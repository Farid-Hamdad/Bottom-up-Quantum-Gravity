#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 16 — Flat torus spectrum convergence v1

Goal
----
Second concrete test for Paper 16:

    L_N = (D-W)/epsilon  ->  -c * Δ_{T^2}

on the flat 2D torus.

For the unit flat torus [0,1)^2 with periodic boundary conditions, the
positive Laplace-Beltrami spectrum is

    λ_{m,n} = 4π² (m²+n²),

with integer pairs (m,n), not both zero.

Graph construction
------------------
- m x m periodic grid, N_actual=m²
- intrinsic periodic distance on the flat torus
- local Gaussian weights on kNN graph
- rescaled graph Laplacian L_N=(D-W)/epsilon

As in the circle test, the raw graph Laplacian has an unknown multiplicative
normalization, so we fit one scalar c such that

    c * λ_graph ≈ λ_torus

on the first modes.

Recommended run
---------------
cd ~/bottomup

python3 papers/paper16_spectral_continuum_limit/scripts/paper16_flat_torus_spectrum_convergence_v1.py \
  --Ns 64 144 256 576 1024 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --n-modes 20 \
  --fit-modes 12 \
  --output-dir papers/paper16_spectral_continuum_limit/results/flat_torus_spectrum_convergence_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Paper 16 flat torus spectrum convergence v1.")
    p.add_argument("--Ns", nargs="+", type=int, default=[64, 144, 256, 576, 1024])
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=128)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--n-modes", type=int, default=20)
    p.add_argument("--fit-modes", type=int, default=12)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def choose_k(N, args):
    if args.k_mode == "fixed":
        k = args.k_fixed
    elif args.k_mode == "sqrt":
        k = int(round(np.sqrt(N)))
    elif args.k_mode == "log":
        k = int(round(4 * np.log(max(N, 3))))
    else:
        raise ValueError(args.k_mode)
    return int(max(args.k_min, min(args.k_max, k, N - 2)))


def make_flat_torus_points(N):
    m = int(round(np.sqrt(N)))
    u = np.linspace(0.0, 1.0, m, endpoint=False)
    v = np.linspace(0.0, 1.0, m, endpoint=False)
    uu, vv = np.meshgrid(u, v)
    pts = np.column_stack([uu.ravel(), vv.ravel()])
    return pts, m, len(pts)


def flat_torus_distance_matrix(pts):
    du = np.abs(pts[:, None, 0] - pts[None, :, 0])
    dv = np.abs(pts[:, None, 1] - pts[None, :, 1])
    du = np.minimum(du, 1.0 - du)
    dv = np.minimum(dv, 1.0 - dv)
    return np.sqrt(du*du + dv*dv)


def build_graph(pts, k, kernel_factor):
    D = flat_torus_distance_matrix(pts)
    n = D.shape[0]
    Dw = D.copy()
    np.fill_diagonal(Dw, np.inf)

    kk = min(k, n - 2)
    knn = np.partition(Dw, kk, axis=1)[:, kk]
    eps = kernel_factor * float(np.median(knn[np.isfinite(knn)])**2)
    eps = max(eps, 1e-14)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(Dw[i], kk)[:k]
        W[i, idx] = np.exp(-(Dw[i, idx]**2)/(4.0*eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)
    return W, eps


def graph_laplacian_rescaled(W, eps):
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return L / eps, deg


def analytic_flat_torus_spectrum(n_nonzero, mmax=20):
    vals = []
    for m in range(-mmax, mmax + 1):
        for n in range(-mmax, mmax + 1):
            if m == 0 and n == 0:
                continue
            vals.append(4.0 * np.pi**2 * (m*m + n*n))
    vals = np.array(sorted(vals), dtype=float)
    if len(vals) < n_nonzero:
        raise ValueError("mmax too small for requested number of modes")
    return vals[:n_nonzero]


def fit_scale(graph_nonzero, target_nonzero, fit_modes):
    x = graph_nonzero[:fit_modes]
    y = target_nonzero[:fit_modes]
    denom = float(np.dot(x, x))
    if denom <= 0:
        return np.nan
    return float(np.dot(x, y) / denom)


def run_case(N, args):
    pts, m_grid, N_actual = make_flat_torus_points(N)
    k = choose_k(N_actual, args)

    W, eps = build_graph(pts, k, args.kernel_factor)
    L, deg = graph_laplacian_rescaled(W, eps)
    evals = np.linalg.eigvalsh(L)
    evals = np.sort(np.maximum(evals, 0.0))

    nonzero = evals[evals > 1e-10]
    n_modes = min(args.n_modes, len(nonzero))
    graph_modes = nonzero[:n_modes]
    target_modes = analytic_flat_torus_spectrum(n_modes)

    fit_modes = min(args.fit_modes, n_modes)
    scale = fit_scale(graph_modes, target_modes, fit_modes)
    scaled = scale * graph_modes

    rel_err = np.abs(scaled - target_modes) / np.maximum(target_modes, 1e-12)
    abs_err = np.abs(scaled - target_modes)

    rows = []
    for idx in range(n_modes):
        rows.append({
            "N_input": N,
            "N_actual": N_actual,
            "m_grid": m_grid,
            "k": k,
            "epsilon": eps,
            "mode_index": idx + 1,
            "graph_lambda_raw": float(graph_modes[idx]),
            "scale_c": scale,
            "graph_lambda_scaled": float(scaled[idx]),
            "target_lambda": float(target_modes[idx]),
            "abs_error": float(abs_err[idx]),
            "rel_error": float(rel_err[idx]),
        })

    summary = {
        "N_input": N,
        "N_actual": N_actual,
        "m_grid": m_grid,
        "k": k,
        "epsilon": eps,
        "scale_c": scale,
        "n_modes": n_modes,
        "fit_modes": fit_modes,
        "mean_rel_error": float(np.mean(rel_err)),
        "median_rel_error": float(np.median(rel_err)),
        "max_rel_error": float(np.max(rel_err)),
        "mean_abs_error": float(np.mean(abs_err)),
        "median_abs_error": float(np.median(abs_err)),
        "degree_mean": float(deg.mean()),
        "degree_std": float(deg.std()),
        "lambda1_raw": float(graph_modes[0]) if n_modes else np.nan,
        "lambda1_scaled": float(scaled[0]) if n_modes else np.nan,
        "lambda1_target": float(target_modes[0]) if n_modes else np.nan,
    }
    return pd.DataFrame(rows), summary


def make_figures(rows_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    plt.figure(figsize=(7, 5))
    plt.plot(summary_df["N_actual"], summary_df["mean_rel_error"], marker="o", label="mean rel error")
    plt.plot(summary_df["N_actual"], summary_df["median_rel_error"], marker="o", label="median rel error")
    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xlabel("N actual")
    plt.ylabel("relative error")
    plt.title("Flat torus spectrum convergence — relative error")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_flat_torus_spectrum_error_vs_N.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(summary_df["N_actual"], summary_df["scale_c"], marker="o")
    plt.xscale("log", base=2)
    plt.xlabel("N actual")
    plt.ylabel("fitted scale c")
    plt.title("Fitted spectral normalization c vs N")
    plt.tight_layout()
    plt.savefig(figdir / "fig_flat_torus_spectrum_scale_vs_N.png", dpi=180)
    plt.close()

    Nmax = int(summary_df["N_actual"].max())
    sub = rows_df[rows_df["N_actual"] == Nmax].sort_values("mode_index")
    plt.figure(figsize=(7, 5))
    plt.plot(sub["mode_index"], sub["target_lambda"], marker="o", label="analytic T2")
    plt.plot(sub["mode_index"], sub["graph_lambda_scaled"], marker="x", label="scaled graph")
    plt.xlabel("nonzero mode index")
    plt.ylabel("eigenvalue")
    plt.title(f"Flat torus spectrum comparison, N={Nmax}")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_flat_torus_spectrum_comparison_bestN.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.bar(sub["mode_index"], sub["rel_error"])
    plt.xlabel("nonzero mode index")
    plt.ylabel("relative error")
    plt.title(f"Per-mode relative error, N={Nmax}")
    plt.tight_layout()
    plt.savefig(figdir / "fig_flat_torus_spectrum_per_mode_error_bestN.png", dpi=180)
    plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summaries = []

    print("=" * 100)
    print("Paper 16 — Flat torus spectrum convergence v1")
    print("=" * 100)

    for N in args.Ns:
        print(f"[run] N_input={N}")
        rows, summary = run_case(N, args)
        all_rows.append(rows)
        summaries.append(summary)
        print(
            f"      N_actual={summary['N_actual']} grid={summary['m_grid']}x{summary['m_grid']} "
            f"k={summary['k']} eps={summary['epsilon']:.6g} "
            f"scale={summary['scale_c']:.6g} "
            f"mean_rel={summary['mean_rel_error']:.6g} "
            f"median_rel={summary['median_rel_error']:.6g} "
            f"max_rel={summary['max_rel_error']:.6g}"
        )

    rows_df = pd.concat(all_rows, ignore_index=True)
    summary_df = pd.DataFrame(summaries)

    rows_df.to_csv(outdir / "flat_torus_spectrum_rows.csv", index=False)
    summary_df.to_csv(outdir / "flat_torus_spectrum_summary.csv", index=False)

    final = summary_df.sort_values("N_actual").tail(1).iloc[0].to_dict()

    with open(outdir / "flat_torus_spectrum_summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 16 flat torus spectrum convergence v1",
            "Ns": args.Ns,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "n_modes": args.n_modes,
            "fit_modes": args.fit_modes,
            "final": final,
        }, f, indent=2)

    make_figures(rows_df, summary_df, outdir)

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(summary_df.to_string(index=False))

    print("\n" + "=" * 100)
    print("FINAL N")
    print("=" * 100)
    print(pd.DataFrame([final]).to_string(index=False))

    print("\nFiles written:")
    print(outdir / "flat_torus_spectrum_rows.csv")
    print(outdir / "flat_torus_spectrum_summary.csv")
    print(outdir / "flat_torus_spectrum_summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

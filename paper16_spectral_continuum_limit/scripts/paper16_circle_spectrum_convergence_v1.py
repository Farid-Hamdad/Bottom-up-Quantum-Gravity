#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 16 — Circle spectrum convergence v1

Goal
----
First concrete test for Paper 16:

    L_N = (D - W) / epsilon  ->  -c * Δ_{S^1}

on the circle.

For a unit circle parameterized by theta in [0, 2π), the positive
Laplace-Beltrami spectrum is:

    λ_m = m^2,   m = 0, 1, 1, 2, 2, 3, 3, ...

Graph construction:
- N equally spaced points on S^1
- intrinsic geodesic distance
- local Gaussian weights on kNN graph
- rescaled graph Laplacian L_N=(D-W)/epsilon

Because the raw graph Laplacian has an unknown multiplicative normalization,
we fit one scalar scale c so that c*λ_graph ≈ λ_circle on the first modes.

Recommended run
---------------
cd ~/bottomup

python3 papers/paper16_spectral_continuum_limit/scripts/paper16_circle_spectrum_convergence_v1.py \
  --Ns 64 128 256 512 1024 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --n-modes 12 \
  --fit-modes 8 \
  --output-dir papers/paper16_spectral_continuum_limit/results/circle_spectrum_convergence_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Paper 16 circle spectrum convergence v1.")
    p.add_argument("--Ns", nargs="+", type=int, default=[64, 128, 256, 512, 1024])
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=96)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--n-modes", type=int, default=12)
    p.add_argument("--fit-modes", type=int, default=8)
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


def circle_distance_matrix(N):
    theta = np.linspace(0.0, 2.0*np.pi, N, endpoint=False)
    dtheta = np.abs(theta[:, None] - theta[None, :])
    D = np.minimum(dtheta, 2.0*np.pi - dtheta)
    return theta, D


def build_circle_graph(N, k, kernel_factor):
    theta, D = circle_distance_matrix(N)
    Dw = D.copy()
    np.fill_diagonal(Dw, np.inf)

    kk = min(k, N - 2)
    knn = np.partition(Dw, kk, axis=1)[:, kk]
    eps = kernel_factor * float(np.median(knn[np.isfinite(knn)]) ** 2)
    eps = max(eps, 1e-14)

    W = np.zeros((N, N), dtype=float)
    for i in range(N):
        idx = np.argpartition(Dw[i], kk)[:k]
        W[i, idx] = np.exp(-(Dw[i, idx] ** 2) / (4.0 * eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)
    return W, eps


def graph_laplacian_rescaled(W, eps):
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return L / eps, deg


def analytic_circle_spectrum(n_nonzero):
    vals = []
    m = 1
    while len(vals) < n_nonzero:
        vals.extend([m*m, m*m])
        m += 1
    return np.array(vals[:n_nonzero], dtype=float)


def fit_scale(graph_nonzero, target_nonzero, fit_modes):
    x = graph_nonzero[:fit_modes]
    y = target_nonzero[:fit_modes]
    denom = float(np.dot(x, x))
    if denom <= 0:
        return np.nan
    return float(np.dot(x, y) / denom)


def run_case(N, args):
    k = choose_k(N, args)
    W, eps = build_circle_graph(N, k, args.kernel_factor)
    L, deg = graph_laplacian_rescaled(W, eps)
    evals = np.linalg.eigvalsh(L)
    evals = np.sort(np.maximum(evals, 0.0))

    nonzero = evals[evals > 1e-10]
    n_modes = min(args.n_modes, len(nonzero))
    graph_modes = nonzero[:n_modes]
    target_modes = analytic_circle_spectrum(n_modes)

    fit_modes = min(args.fit_modes, n_modes)
    scale = fit_scale(graph_modes, target_modes, fit_modes)
    scaled = scale * graph_modes

    rel_err = np.abs(scaled - target_modes) / np.maximum(target_modes, 1e-12)
    abs_err = np.abs(scaled - target_modes)

    rows = []
    for idx in range(n_modes):
        rows.append({
            "N": N,
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
        "N": N,
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
    }
    return pd.DataFrame(rows), summary


def make_figures(rows_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    plt.figure(figsize=(7, 5))
    plt.plot(summary_df["N"], summary_df["mean_rel_error"], marker="o", label="mean rel error")
    plt.plot(summary_df["N"], summary_df["median_rel_error"], marker="o", label="median rel error")
    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xlabel("N")
    plt.ylabel("relative error")
    plt.title("Circle spectrum convergence — relative error")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_circle_spectrum_error_vs_N.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(summary_df["N"], summary_df["scale_c"], marker="o")
    plt.xscale("log", base=2)
    plt.xlabel("N")
    plt.ylabel("fitted scale c")
    plt.title("Fitted spectral normalization c vs N")
    plt.tight_layout()
    plt.savefig(figdir / "fig_circle_spectrum_scale_vs_N.png", dpi=180)
    plt.close()

    Nmax = int(summary_df["N"].max())
    sub = rows_df[rows_df["N"] == Nmax].sort_values("mode_index")
    plt.figure(figsize=(7, 5))
    plt.plot(sub["mode_index"], sub["target_lambda"], marker="o", label="analytic S1")
    plt.plot(sub["mode_index"], sub["graph_lambda_scaled"], marker="x", label="scaled graph")
    plt.xlabel("nonzero mode index")
    plt.ylabel("eigenvalue")
    plt.title(f"Circle spectrum comparison, N={Nmax}")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_circle_spectrum_comparison_bestN.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.bar(sub["mode_index"], sub["rel_error"])
    plt.xlabel("nonzero mode index")
    plt.ylabel("relative error")
    plt.title(f"Per-mode relative error, N={Nmax}")
    plt.tight_layout()
    plt.savefig(figdir / "fig_circle_spectrum_per_mode_error_bestN.png", dpi=180)
    plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summaries = []

    print("=" * 100)
    print("Paper 16 — Circle spectrum convergence v1")
    print("=" * 100)

    for N in args.Ns:
        print(f"[run] N={N}")
        rows, summary = run_case(N, args)
        all_rows.append(rows)
        summaries.append(summary)
        print(
            f"      k={summary['k']} eps={summary['epsilon']:.6g} "
            f"scale={summary['scale_c']:.6g} "
            f"mean_rel={summary['mean_rel_error']:.6g} "
            f"median_rel={summary['median_rel_error']:.6g} "
            f"max_rel={summary['max_rel_error']:.6g}"
        )

    rows_df = pd.concat(all_rows, ignore_index=True)
    summary_df = pd.DataFrame(summaries)

    rows_df.to_csv(outdir / "circle_spectrum_rows.csv", index=False)
    summary_df.to_csv(outdir / "circle_spectrum_summary.csv", index=False)

    final = summary_df.sort_values("N").tail(1).iloc[0].to_dict()

    with open(outdir / "circle_spectrum_summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 16 circle spectrum convergence v1",
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
    print(outdir / "circle_spectrum_rows.csv")
    print(outdir / "circle_spectrum_summary.csv")
    print(outdir / "circle_spectrum_summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

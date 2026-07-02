#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Spectral convergence test v1

First numerical test for Paper 15:

    W_ij -> L_ent -> Δ_g

We sample known geometries, build a local weighted graph W_ij, compute the
graph Laplacian, and test whether the heat trace and spectral dimension
stabilize as N increases.

Recommended run:
cd ~/bottomup
python3 papers/paper15_einstein_derivation/scripts/paper15_spectral_convergence_v1.py \
  --geometries interval circle grid2d sphere \
  --Ns 64 128 256 512 \
  --k 12 \
  --kernel-scale 1.5 \
  --output-dir papers/paper15_einstein_derivation/results/spectral_convergence_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Paper 15 spectral convergence test v1.")
    p.add_argument("--geometries", nargs="+", default=["interval", "circle", "grid2d", "sphere"])
    p.add_argument("--Ns", nargs="+", type=int, default=[64, 128, 256, 512])
    p.add_argument("--k", type=int, default=12)
    p.add_argument("--kernel-scale", type=float, default=1.5)
    p.add_argument("--t-min", type=float, default=1e-2)
    p.add_argument("--t-max", type=float, default=1e2)
    p.add_argument("--n-t", type=int, default=160)
    p.add_argument("--plateau-qlo", type=float, default=0.25)
    p.add_argument("--plateau-qhi", type=float, default=0.65)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def make_points(geometry, N, seed):
    if geometry == "interval":
        return np.linspace(0.0, 1.0, N)[:, None], 1

    if geometry == "circle":
        theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
        return np.column_stack([np.cos(theta), np.sin(theta)]), 1

    if geometry == "grid2d":
        m = int(np.round(np.sqrt(N)))
        xs = np.linspace(0.0, 1.0, m)
        xx, yy = np.meshgrid(xs, xs)
        return np.column_stack([xx.ravel(), yy.ravel()]), 2

    if geometry == "sphere":
        i = np.arange(N)
        phi = np.arccos(1.0 - 2.0 * (i + 0.5) / N)
        golden = np.pi * (3.0 - np.sqrt(5.0))
        theta = golden * i
        return np.column_stack([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi),
        ]), 2

    raise ValueError(f"Unknown geometry: {geometry}")


def pairwise_distances(pts):
    diff = pts[:, None, :] - pts[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def build_local_W(pts, k, kernel_scale):
    D = pairwise_distances(pts)
    n = D.shape[0]
    np.fill_diagonal(D, np.inf)

    kk = min(k, n - 2)
    kth = np.partition(D, kk, axis=1)[:, kk]
    xi = kernel_scale * np.median(kth[np.isfinite(kth)])
    xi = max(float(xi), 1e-12)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(D[i], kk)[:k]
        for j in idx:
            if np.isfinite(D[i, j]):
                W[i, j] = np.exp(-(D[i, j] ** 2) / (xi ** 2))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)
    density = np.count_nonzero(W) / (n * (n - 1))
    return W, xi, density


def normalized_laplacian(W):
    deg = W.sum(axis=1)
    invsqrt = np.zeros_like(deg)
    mask = deg > 0
    invsqrt[mask] = 1.0 / np.sqrt(deg[mask])
    L = np.eye(W.shape[0]) - invsqrt[:, None] * W * invsqrt[None, :]
    return L, deg


def heat_trace_ds(evals, times):
    evals = np.maximum(evals, 0.0)
    Z = np.array([np.sum(np.exp(-t * evals)) for t in times])
    logt = np.log(times)
    logZ = np.log(np.maximum(Z, 1e-300))
    ds = -2.0 * np.gradient(logZ, logt)
    return Z, ds


def plateau(times, ds, qlo, qhi, target_dim):
    n = len(times)
    ilo = int(np.floor(qlo * n))
    ihi = int(np.ceil(qhi * n))
    ilo = max(0, min(ilo, n - 1))
    ihi = max(ilo + 1, min(ihi, n))
    w = ds[ilo:ihi]
    return {
        "plateau_t_min": float(times[ilo]),
        "plateau_t_max": float(times[ihi - 1]),
        "ds_plateau_median": float(np.nanmedian(w)),
        "ds_plateau_mean": float(np.nanmean(w)),
        "ds_plateau_std": float(np.nanstd(w)),
        "ds_abs_error": float(abs(np.nanmedian(w) - target_dim)),
    }


def lambda2(evals):
    ev = np.sort(np.maximum(evals, 0.0))
    nz = ev[ev > 1e-10]
    return float(nz[0]) if len(nz) else np.nan


def run_case(geom, N, args):
    pts, target_dim = make_points(geom, N, args.seed + N)
    W, xi, density = build_local_W(pts, args.k, args.kernel_scale)
    L, deg = normalized_laplacian(W)
    evals = np.linalg.eigvalsh(L)

    times = np.logspace(np.log10(args.t_min), np.log10(args.t_max), args.n_t)
    Z, ds = heat_trace_ds(evals, times)

    row = {
        "geometry": geom,
        "N_input": N,
        "N_actual": len(pts),
        "target_dim": target_dim,
        "k": args.k,
        "kernel_scale": args.kernel_scale,
        "xi": xi,
        "density": density,
        "mean_degree_weighted": float(deg.mean()),
        "std_degree_weighted": float(deg.std()),
        "lambda2": lambda2(evals),
        "eval_max": float(evals.max()),
    }
    row.update(plateau(times, ds, args.plateau_qlo, args.plateau_qhi, target_dim))

    curve = pd.DataFrame({
        "geometry": geom,
        "N_actual": len(pts),
        "target_dim": target_dim,
        "t": times,
        "Z": Z,
        "ds": ds,
    })
    return row, curve


def make_figures(rows_df, curves_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    for geom in sorted(curves_df.geometry.unique()):
        sub = curves_df[curves_df.geometry == geom]

        plt.figure(figsize=(7, 5))
        for N in sorted(sub.N_actual.unique()):
            s = sub[sub.N_actual == N]
            plt.plot(s.t, s.ds, label=f"N={N}")
        target = sub.target_dim.iloc[0]
        plt.axhline(target, linestyle="--", label=f"target D={target}")
        plt.xscale("log")
        plt.xlabel("t")
        plt.ylabel("d_s(t)")
        plt.title(f"Spectral dimension — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_ds_curves_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(7, 5))
        for N in sorted(sub.N_actual.unique()):
            s = sub[sub.N_actual == N]
            plt.plot(s.t, s.Z, label=f"N={N}")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("t")
        plt.ylabel("Z(t)")
        plt.title(f"Heat trace — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_heat_trace_{geom}.png", dpi=180)
        plt.close()

    plt.figure(figsize=(7, 5))
    for geom in sorted(rows_df.geometry.unique()):
        s = rows_df[rows_df.geometry == geom].sort_values("N_actual")
        plt.plot(s.N_actual, s.ds_plateau_median, marker="o", label=geom)
    plt.xscale("log", base=2)
    plt.xlabel("N")
    plt.ylabel("median plateau d_s")
    plt.title("Spectral-dimension plateau vs N")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_plateau_vs_N.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    for geom in sorted(rows_df.geometry.unique()):
        s = rows_df[rows_df.geometry == geom].sort_values("N_actual")
        plt.plot(s.N_actual, s.ds_abs_error, marker="o", label=geom)
    plt.xscale("log", base=2)
    plt.xlabel("N")
    plt.ylabel("|plateau d_s - target D|")
    plt.title("Spectral-dimension error vs N")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_plateau_error_vs_N.png", dpi=180)
    plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    curves = []

    print("=" * 100)
    print("Paper 15 — Spectral convergence test v1")
    print("=" * 100)

    for geom in args.geometries:
        for N in args.Ns:
            print(f"[run] geometry={geom} N={N}")
            row, curve = run_case(geom, N, args)
            rows.append(row)
            curves.append(curve)
            print(f"      N_actual={row['N_actual']} targetD={row['target_dim']} "
                  f"ds_med={row['ds_plateau_median']:.4f} "
                  f"err={row['ds_abs_error']:.4f} lambda2={row['lambda2']:.4g}")

    rows_df = pd.DataFrame(rows)
    curves_df = pd.concat(curves, ignore_index=True)

    rows_df.to_csv(outdir / "spectral_convergence_rows.csv", index=False)
    curves_df.to_csv(outdir / "spectral_dimension_curves.csv", index=False)

    by_geom = rows_df.groupby("geometry").agg(
        n_cases=("geometry", "count"),
        target_dim=("target_dim", "first"),
        best_N=("N_actual", "max"),
        final_ds_plateau=("ds_plateau_median", "last"),
        final_ds_abs_error=("ds_abs_error", "last"),
        mean_abs_error=("ds_abs_error", "mean"),
        final_lambda2=("lambda2", "last"),
        final_density=("density", "last"),
    ).reset_index()
    by_geom.to_csv(outdir / "summary_by_geometry.csv", index=False)

    summary = {
        "experiment": "Paper 15 spectral convergence test v1",
        "geometries": args.geometries,
        "Ns": args.Ns,
        "k": args.k,
        "kernel_scale": args.kernel_scale,
        "n_rows": int(len(rows_df)),
        "mean_abs_error_all": float(rows_df.ds_abs_error.mean()),
        "median_abs_error_all": float(rows_df.ds_abs_error.median()),
        "by_geometry": by_geom.to_dict(orient="records"),
    }
    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    make_figures(rows_df, curves_df, outdir)

    print("\n" + "=" * 100)
    print("SUMMARY BY GEOMETRY")
    print("=" * 100)
    print(by_geom.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "spectral_convergence_rows.csv")
    print(outdir / "spectral_dimension_curves.csv")
    print(outdir / "summary_by_geometry.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

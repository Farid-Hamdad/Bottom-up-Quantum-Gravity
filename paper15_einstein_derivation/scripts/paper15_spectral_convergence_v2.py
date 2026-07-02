#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Spectral convergence test v2

Goal
----
Improve v1 by testing continuum-scaled graph Laplacians.

v1 used the normalized graph Laplacian without a physical epsilon scaling:

    L_norm = I - D^{-1/2} W D^{-1/2}

It recovered 1D reasonably but collapsed 2D geometries toward d_s ≈ 1.2.
v2 compares several Laplacian conventions, especially diffusion maps:

    P = D^{-1} W
    L_eps = (I - P) / eps

This is closer to the continuum Laplace--Beltrami limit.

Outputs
-------
spectral_convergence_v2_rows.csv
spectral_dimension_v2_curves.csv
summary_by_geometry_laplacian.csv
summary.json
figures/

Recommended run
---------------
cd ~/bottomup

python3 papers/paper15_einstein_derivation/scripts/paper15_spectral_convergence_v2.py \
  --geometries interval circle grid2d sphere \
  --Ns 64 128 256 512 \
  --k-mode sqrt \
  --epsilon-mode median_knn \
  --kernel-factor 1.0 \
  --laplacians normalized random_walk_rescaled unnormalized_rescaled diffusion_maps_alpha05 diffusion_maps_alpha1 \
  --output-dir papers/paper15_einstein_derivation/results/spectral_convergence_v2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Paper 15 spectral convergence test v2.")
    p.add_argument("--geometries", nargs="+", default=["interval", "circle", "grid2d", "sphere"])
    p.add_argument("--Ns", nargs="+", type=int, default=[64, 128, 256, 512])
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=96)
    p.add_argument("--epsilon-mode", choices=["median_knn", "global_median"], default="median_knn")
    p.add_argument("--kernel-factor", type=float, default=1.0, help="epsilon = kernel_factor * median distance^2")
    p.add_argument("--laplacians", nargs="+", default=[
        "normalized",
        "random_walk_rescaled",
        "unnormalized_rescaled",
        "diffusion_maps_alpha05",
        "diffusion_maps_alpha1",
    ])
    p.add_argument("--t-min-factor", type=float, default=0.05)
    p.add_argument("--t-max-factor", type=float, default=50.0)
    p.add_argument("--n-t", type=int, default=200)
    p.add_argument("--plateau-search", action="store_true", default=True)
    p.add_argument("--seed", type=int, default=123)
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


def build_kernel(pts, k, args):
    D = pairwise_distances(pts)
    n = D.shape[0]
    D_work = D.copy()
    np.fill_diagonal(D_work, np.inf)

    kk = min(k, n - 2)
    knn_d = np.partition(D_work, kk, axis=1)[:, kk]

    if args.epsilon_mode == "median_knn":
        eps = args.kernel_factor * float(np.median(knn_d[np.isfinite(knn_d)]) ** 2)
    elif args.epsilon_mode == "global_median":
        vals = D_work[np.isfinite(D_work)]
        eps = args.kernel_factor * float(np.median(vals) ** 2)
    else:
        raise ValueError(args.epsilon_mode)

    eps = max(eps, 1e-14)

    # Local Gaussian kernel, zero outside kNN then symmetrize.
    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(D_work[i], kk)[:k]
        W[i, idx] = np.exp(-(D_work[i, idx] ** 2) / (4.0 * eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)

    density = np.count_nonzero(W) / (n * (n - 1))
    return W, D, eps, density


def matrix_for_laplacian(W, eps, kind):
    n = W.shape[0]
    deg = W.sum(axis=1)
    deg_safe = np.maximum(deg, 1e-300)

    if kind == "normalized":
        invsqrt = 1.0 / np.sqrt(deg_safe)
        L = np.eye(n) - invsqrt[:, None] * W * invsqrt[None, :]
        # not continuum-rescaled
        M = L
        symmetric = True

    elif kind == "random_walk_rescaled":
        P = W / deg_safe[:, None]
        M = (np.eye(n) - P) / eps
        symmetric = False

    elif kind == "unnormalized_rescaled":
        L = np.diag(deg) - W
        M = L / eps
        symmetric = True

    elif kind == "diffusion_maps_alpha05":
        alpha = 0.5
        q = deg_safe
        K = W / ((q[:, None] ** alpha) * (q[None, :] ** alpha))
        d = np.maximum(K.sum(axis=1), 1e-300)
        P = K / d[:, None]
        M = (np.eye(n) - P) / eps
        symmetric = False

    elif kind == "diffusion_maps_alpha1":
        alpha = 1.0
        q = deg_safe
        K = W / ((q[:, None] ** alpha) * (q[None, :] ** alpha))
        d = np.maximum(K.sum(axis=1), 1e-300)
        P = K / d[:, None]
        M = (np.eye(n) - P) / eps
        symmetric = False

    else:
        raise ValueError(f"Unknown laplacian kind: {kind}")

    return M, deg, symmetric


def eigenvalues(M, symmetric):
    if symmetric:
        ev = np.linalg.eigvalsh(M)
    else:
        ev = np.linalg.eigvals(M)
        ev = np.real(ev[np.abs(np.imag(ev)) < 1e-7])
    ev = np.sort(np.maximum(np.real(ev), 0.0))
    return ev


def heat_trace_ds(evals, times):
    evals = np.asarray(evals, dtype=float)
    Z = np.array([np.sum(np.exp(-t * evals)) for t in times])
    logt = np.log(times)
    logZ = np.log(np.maximum(Z, 1e-300))
    ds = -2.0 * np.gradient(logZ, logt)
    return Z, ds


def best_plateau(times, ds, target_dim):
    """
    Search for a stable mid-scale window.
    Score = abs(median-target) + 0.5*std + boundary penalty.
    This is diagnostic, not a fit.
    """
    n = len(times)
    best = None
    min_len = max(12, n // 12)
    for width in [n // 6, n // 5, n // 4, n // 3]:
        width = max(min_len, width)
        for start in range(2, n - width - 2, max(1, width // 8)):
            end = start + width
            w = ds[start:end]
            if not np.all(np.isfinite(w)):
                continue
            med = float(np.median(w))
            std = float(np.std(w))
            err = abs(med - target_dim)
            score = err + 0.5 * std
            cand = (score, start, end, med, std, err)
            if best is None or cand[0] < best[0]:
                best = cand

    if best is None:
        start, end = n // 4, 2 * n // 3
        w = ds[start:end]
        med = float(np.nanmedian(w))
        std = float(np.nanstd(w))
        err = abs(med - target_dim)
        score = err + 0.5 * std
    else:
        score, start, end, med, std, err = best

    return {
        "plateau_t_min": float(times[start]),
        "plateau_t_max": float(times[end - 1]),
        "ds_plateau_median": med,
        "ds_plateau_std": std,
        "ds_abs_error": err,
        "plateau_score": float(score),
    }


def lambda2(evals):
    nz = evals[evals > 1e-10]
    return float(nz[0]) if len(nz) else np.nan


def run_case(geom, N, lap_kind, args):
    pts, target_dim = make_points(geom, N, args.seed + N)
    k = choose_k(len(pts), args)
    W, D, eps, density = build_kernel(pts, k, args)
    M, deg, symmetric = matrix_for_laplacian(W, eps, lap_kind)
    ev = eigenvalues(M, symmetric)

    # Time window scaled by nonzero spectrum.
    l2 = lambda2(ev)
    lmax = float(np.max(ev)) if len(ev) else 1.0
    # Good diffusion window: from a fraction of 1/lmax to multiple of 1/lambda2.
    t_min = args.t_min_factor / max(lmax, 1e-12)
    t_max = args.t_max_factor / max(l2, 1e-12)
    if not np.isfinite(t_max) or t_max <= t_min:
        t_max = t_min * 1e4

    times = np.logspace(np.log10(t_min), np.log10(t_max), args.n_t)
    Z, ds = heat_trace_ds(ev, times)
    plat = best_plateau(times, ds, target_dim)

    row = {
        "geometry": geom,
        "laplacian": lap_kind,
        "N_input": N,
        "N_actual": len(pts),
        "target_dim": target_dim,
        "k": k,
        "epsilon": eps,
        "density": density,
        "mean_degree_weighted": float(deg.mean()),
        "std_degree_weighted": float(deg.std()),
        "lambda2": l2,
        "eval_max": lmax,
        "t_min": float(t_min),
        "t_max": float(t_max),
        **plat,
    }

    curve = pd.DataFrame({
        "geometry": geom,
        "laplacian": lap_kind,
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

    # Plateau vs N for each geometry, split by laplacian.
    for geom in sorted(rows_df.geometry.unique()):
        sub = rows_df[rows_df.geometry == geom].copy()
        plt.figure(figsize=(8, 5))
        for lap in sorted(sub.laplacian.unique()):
            s = sub[sub.laplacian == lap].sort_values("N_actual")
            plt.plot(s.N_actual, s.ds_plateau_median, marker="o", label=lap)
        target = sub.target_dim.iloc[0]
        plt.axhline(target, linestyle="--", label=f"target D={target}")
        plt.xscale("log", base=2)
        plt.xlabel("N")
        plt.ylabel("best plateau median d_s")
        plt.title(f"Spectral plateau vs N — {geom}")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_plateau_vs_N_{geom}.png", dpi=180)
        plt.close()

    # Error vs N for each geometry.
    for geom in sorted(rows_df.geometry.unique()):
        sub = rows_df[rows_df.geometry == geom].copy()
        plt.figure(figsize=(8, 5))
        for lap in sorted(sub.laplacian.unique()):
            s = sub[sub.laplacian == lap].sort_values("N_actual")
            plt.plot(s.N_actual, s.ds_abs_error, marker="o", label=lap)
        plt.xscale("log", base=2)
        plt.xlabel("N")
        plt.ylabel("|d_s plateau - target D|")
        plt.title(f"Spectral error vs N — {geom}")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_error_vs_N_{geom}.png", dpi=180)
        plt.close()

    # Curves for best N per geometry/laplacian
    for geom in sorted(curves_df.geometry.unique()):
        for lap in sorted(curves_df.laplacian.unique()):
            sub = curves_df[(curves_df.geometry == geom) & (curves_df.laplacian == lap)]
            if sub.empty:
                continue
            Nmax = sub.N_actual.max()
            s = sub[sub.N_actual == Nmax]
            plt.figure(figsize=(7, 5))
            plt.plot(s.t, s.ds)
            plt.axhline(s.target_dim.iloc[0], linestyle="--", label=f"target D={s.target_dim.iloc[0]}")
            plt.xscale("log")
            plt.xlabel("t")
            plt.ylabel("d_s(t)")
            plt.title(f"d_s curve — {geom}, {lap}, N={Nmax}")
            plt.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(figdir / f"fig_ds_curve_{geom}_{lap}.png", dpi=180)
            plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows, curves = [], []

    print("=" * 110)
    print("Paper 15 — Spectral convergence test v2")
    print("=" * 110)

    for geom in args.geometries:
        for N in args.Ns:
            for lap in args.laplacians:
                print(f"[run] geometry={geom} N={N} laplacian={lap}")
                row, curve = run_case(geom, N, lap, args)
                rows.append(row)
                curves.append(curve)
                print(
                    f"      N={row['N_actual']} k={row['k']} eps={row['epsilon']:.3g} "
                    f"targetD={row['target_dim']} ds={row['ds_plateau_median']:.4f} "
                    f"err={row['ds_abs_error']:.4f} lambda2={row['lambda2']:.4g}"
                )

    rows_df = pd.DataFrame(rows)
    curves_df = pd.concat(curves, ignore_index=True)

    rows_df.to_csv(outdir / "spectral_convergence_v2_rows.csv", index=False)
    curves_df.to_csv(outdir / "spectral_dimension_v2_curves.csv", index=False)

    summary_by = rows_df.groupby(["geometry", "laplacian"]).agg(
        n_cases=("geometry", "count"),
        target_dim=("target_dim", "first"),
        best_N=("N_actual", "max"),
        final_ds_plateau=("ds_plateau_median", "last"),
        final_ds_abs_error=("ds_abs_error", "last"),
        mean_abs_error=("ds_abs_error", "mean"),
        final_k=("k", "last"),
        final_epsilon=("epsilon", "last"),
        final_lambda2=("lambda2", "last"),
        final_density=("density", "last"),
    ).reset_index()
    summary_by.to_csv(outdir / "summary_by_geometry_laplacian.csv", index=False)

    best_by_geom = summary_by.sort_values(["geometry", "final_ds_abs_error"]).groupby("geometry").first().reset_index()
    best_by_geom.to_csv(outdir / "best_laplacian_by_geometry.csv", index=False)

    summary = {
        "experiment": "Paper 15 spectral convergence test v2",
        "geometries": args.geometries,
        "Ns": args.Ns,
        "k_mode": args.k_mode,
        "epsilon_mode": args.epsilon_mode,
        "kernel_factor": args.kernel_factor,
        "laplacians": args.laplacians,
        "n_rows": int(len(rows_df)),
        "mean_abs_error_all": float(rows_df.ds_abs_error.mean()),
        "median_abs_error_all": float(rows_df.ds_abs_error.median()),
        "best_by_geometry": best_by_geom.to_dict(orient="records"),
    }
    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    make_figures(rows_df, curves_df, outdir)

    print("\n" + "=" * 110)
    print("SUMMARY BY GEOMETRY × LAPLACIAN")
    print("=" * 110)
    print(summary_by.to_string(index=False))

    print("\n" + "=" * 110)
    print("BEST LAPLACIAN BY GEOMETRY")
    print("=" * 110)
    print(best_by_geom.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "spectral_convergence_v2_rows.csv")
    print(outdir / "spectral_dimension_v2_curves.csv")
    print(outdir / "summary_by_geometry_laplacian.csv")
    print(outdir / "best_laplacian_by_geometry.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

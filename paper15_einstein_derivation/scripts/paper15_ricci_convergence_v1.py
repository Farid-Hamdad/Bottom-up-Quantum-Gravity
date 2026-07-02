#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Ricci curvature convergence test v1

Goal
----
Second numerical test for Paper 15:

    kappa_ij^OR -> Ricci curvature signal

We build local weighted graphs from controlled geometries and compute
Ollivier--Ricci curvature on sampled edges using optimal transport.

The target is not yet a rigorous continuum proof. The goal is a controlled
numerical separation:

    grid2d / torus  -> near-flat average curvature
    sphere          -> positive curvature signal
    circle          -> positive compact 1D signal / finite-size effects

Outputs
-------
ricci_curvature_rows.csv
summary_by_geometry.csv
summary.json
figures/

Recommended run
---------------
cd ~/bottomup

python3 papers/paper15_einstein_derivation/scripts/paper15_ricci_convergence_v1.py \
  --geometries circle grid2d sphere torus \
  --Ns 128 256 512 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --max-edges 500 \
  --alpha-idleness 0.5 \
  --output-dir papers/paper15_einstein_derivation/results/ricci_convergence_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.sparse.csgraph import shortest_path


def parse_args():
    p = argparse.ArgumentParser(description="Paper 15 Ricci curvature convergence test v1.")
    p.add_argument("--geometries", nargs="+", default=["circle", "grid2d", "sphere", "torus"])
    p.add_argument("--Ns", nargs="+", type=int, default=[128, 256, 512])
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=64)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--alpha-idleness", type=float, default=0.5, help="Ollivier idleness mass at source node")
    p.add_argument("--max-edges", type=int, default=500, help="max sampled edges per geometry/N")
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


def make_points(geometry, N):
    if geometry == "circle":
        theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
        return np.column_stack([np.cos(theta), np.sin(theta)]), 1, "positive_compact_1d"

    if geometry == "grid2d":
        m = int(np.round(np.sqrt(N)))
        xs = np.linspace(0.0, 1.0, m)
        xx, yy = np.meshgrid(xs, xs)
        return np.column_stack([xx.ravel(), yy.ravel()]), 2, "flat_with_boundary"

    if geometry == "sphere":
        i = np.arange(N)
        phi = np.arccos(1.0 - 2.0 * (i + 0.5) / N)
        golden = np.pi * (3.0 - np.sqrt(5.0))
        theta = golden * i
        pts = np.column_stack([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi),
        ])
        return pts, 2, "positive_curvature"

    if geometry == "torus":
        # roughly sqrt(N) x sqrt(N) torus embedded in R3
        m = int(np.round(np.sqrt(N)))
        u = np.linspace(0, 2*np.pi, m, endpoint=False)
        v = np.linspace(0, 2*np.pi, m, endpoint=False)
        uu, vv = np.meshgrid(u, v)
        R, r = 2.0, 0.7
        x = (R + r*np.cos(vv)) * np.cos(uu)
        y = (R + r*np.cos(vv)) * np.sin(uu)
        z = r * np.sin(vv)
        pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        return pts, 2, "mixed_mean_near_zero"

    raise ValueError(f"Unknown geometry: {geometry}")


def pairwise_distances(pts):
    diff = pts[:, None, :] - pts[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def build_weight_graph(pts, k, kernel_factor):
    D = pairwise_distances(pts)
    n = len(pts)
    D_work = D.copy()
    np.fill_diagonal(D_work, np.inf)

    kk = min(k, n - 2)
    knn = np.partition(D_work, kk, axis=1)[:, kk]
    eps = kernel_factor * float(np.median(knn[np.isfinite(knn)]) ** 2)
    eps = max(eps, 1e-14)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(D_work[i], kk)[:k]
        W[i, idx] = np.exp(-(D_work[i, idx] ** 2) / (4.0 * eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)

    # edge length for transport ground distance: use Euclidean/chord distances on edges.
    # Non-edges set to inf for shortest_path.
    A_len = np.where(W > 0, D, np.inf)
    np.fill_diagonal(A_len, 0.0)

    density = np.count_nonzero(W) / (n * (n - 1))
    return W, D, A_len, eps, density


def neighbor_distribution(W, node, alpha):
    n = W.shape[0]
    mu = np.zeros(n, dtype=float)
    mu[node] = alpha
    neigh = np.where(W[node] > 0)[0]
    total = W[node, neigh].sum()
    if total > 0:
        mu[neigh] += (1.0 - alpha) * W[node, neigh] / total
    else:
        mu[node] = 1.0
    return mu


def wasserstein_1(mu, nu, dist_matrix):
    # Restrict to support for speed.
    support_i = np.where(mu > 1e-15)[0]
    support_j = np.where(nu > 1e-15)[0]
    a = mu[support_i]
    b = nu[support_j]
    C = dist_matrix[np.ix_(support_i, support_j)]

    m, n = C.shape
    c = C.ravel()

    # equality constraints: row sums = a, column sums = b
    A_eq = []
    b_eq = []

    for i in range(m):
        row = np.zeros(m * n)
        row[i*n:(i+1)*n] = 1.0
        A_eq.append(row)
        b_eq.append(a[i])

    for j in range(n):
        col = np.zeros(m * n)
        col[j::n] = 1.0
        A_eq.append(col)
        b_eq.append(b[j])

    res = linprog(
        c,
        A_eq=np.array(A_eq),
        b_eq=np.array(b_eq),
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        return np.nan
    return float(res.fun)


def ollivier_ricci_edge(W, spdist, i, j, alpha):
    dij = spdist[i, j]
    if not np.isfinite(dij) or dij <= 0:
        return np.nan
    mu_i = neighbor_distribution(W, i, alpha)
    mu_j = neighbor_distribution(W, j, alpha)
    w1 = wasserstein_1(mu_i, mu_j, spdist)
    if not np.isfinite(w1):
        return np.nan
    return float(1.0 - w1 / dij)


def sample_edges(W, max_edges, seed):
    rng = np.random.default_rng(seed)
    ii, jj = np.where(np.triu(W > 0, k=1))
    edges = np.column_stack([ii, jj])
    if len(edges) > max_edges:
        idx = rng.choice(len(edges), size=max_edges, replace=False)
        edges = edges[idx]
    return edges


def run_case(geometry, N, args):
    pts, target_dim, expected = make_points(geometry, N)
    k = choose_k(len(pts), args)
    W, D, A_len, eps, density = build_weight_graph(pts, k, args.kernel_factor)

    # Graph shortest-path distance on edge lengths.
    spdist = shortest_path(A_len, directed=False, unweighted=False)
    edges = sample_edges(W, args.max_edges, args.seed + N + len(geometry))

    rows = []
    for i, j in edges:
        kappa = ollivier_ricci_edge(W, spdist, int(i), int(j), args.alpha_idleness)
        rows.append({
            "geometry": geometry,
            "N_input": N,
            "N_actual": len(pts),
            "target_dim": target_dim,
            "expected_curvature": expected,
            "k": k,
            "kernel_factor": args.kernel_factor,
            "epsilon": eps,
            "density": density,
            "alpha_idleness": args.alpha_idleness,
            "i": int(i),
            "j": int(j),
            "edge_weight": float(W[i, j]),
            "edge_length": float(D[i, j]),
            "graph_distance": float(spdist[i, j]),
            "kappa_OR": kappa,
        })

    return pd.DataFrame(rows)


def make_figures(edge_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    # Histogram by geometry at largest N
    for geom in sorted(edge_df.geometry.unique()):
        subg = edge_df[edge_df.geometry == geom]
        Nmax = subg.N_actual.max()
        sub = subg[subg.N_actual == Nmax]
        plt.figure(figsize=(7, 5))
        plt.hist(sub.kappa_OR.dropna(), bins=40)
        plt.axvline(sub.kappa_OR.mean(), linestyle="--", label=f"mean={sub.kappa_OR.mean():.3g}")
        plt.xlabel("Ollivier-Ricci kappa")
        plt.ylabel("count")
        plt.title(f"OR curvature distribution — {geom}, N={Nmax}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_kappa_hist_{geom}.png", dpi=180)
        plt.close()

    # Mean curvature vs N
    plt.figure(figsize=(8, 5))
    for geom in sorted(summary_df.geometry.unique()):
        s = summary_df[summary_df.geometry == geom].sort_values("N_actual")
        plt.plot(s.N_actual, s.mean_kappa, marker="o", label=geom)
    plt.axhline(0.0, linestyle="--")
    plt.xscale("log", base=2)
    plt.xlabel("N")
    plt.ylabel("mean OR curvature")
    plt.title("Mean Ollivier-Ricci curvature vs N")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_mean_kappa_vs_N.png", dpi=180)
    plt.close()

    # Median curvature vs N
    plt.figure(figsize=(8, 5))
    for geom in sorted(summary_df.geometry.unique()):
        s = summary_df[summary_df.geometry == geom].sort_values("N_actual")
        plt.plot(s.N_actual, s.median_kappa, marker="o", label=geom)
    plt.axhline(0.0, linestyle="--")
    plt.xscale("log", base=2)
    plt.xlabel("N")
    plt.ylabel("median OR curvature")
    plt.title("Median Ollivier-Ricci curvature vs N")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_median_kappa_vs_N.png", dpi=180)
    plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_edges = []

    print("=" * 110)
    print("Paper 15 — Ricci curvature convergence test v1")
    print("=" * 110)

    for geom in args.geometries:
        for N in args.Ns:
            print(f"[run] geometry={geom} N={N}")
            df = run_case(geom, N, args)
            all_edges.append(df)
            print(
                f"      N_actual={df.N_actual.iloc[0]} k={df.k.iloc[0]} "
                f"edges={len(df)} mean_kappa={df.kappa_OR.mean():.5f} "
                f"median_kappa={df.kappa_OR.median():.5f}"
            )

    edge_df = pd.concat(all_edges, ignore_index=True)
    edge_df.to_csv(outdir / "ricci_curvature_rows.csv", index=False)

    summary = edge_df.groupby(["geometry", "N_actual"]).agg(
        n_edges=("kappa_OR", "count"),
        target_dim=("target_dim", "first"),
        expected_curvature=("expected_curvature", "first"),
        k=("k", "first"),
        epsilon=("epsilon", "first"),
        density=("density", "first"),
        mean_kappa=("kappa_OR", "mean"),
        median_kappa=("kappa_OR", "median"),
        std_kappa=("kappa_OR", "std"),
        q10_kappa=("kappa_OR", lambda x: x.quantile(0.10)),
        q90_kappa=("kappa_OR", lambda x: x.quantile(0.90)),
        positive_fraction=("kappa_OR", lambda x: float((x > 0).mean())),
        negative_fraction=("kappa_OR", lambda x: float((x < 0).mean())),
    ).reset_index()

    summary.to_csv(outdir / "summary_by_geometry.csv", index=False)

    final = summary.sort_values("N_actual").groupby("geometry").tail(1).reset_index(drop=True)
    final.to_csv(outdir / "final_by_geometry.csv", index=False)

    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 15 Ricci curvature convergence test v1",
            "geometries": args.geometries,
            "Ns": args.Ns,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "alpha_idleness": args.alpha_idleness,
            "max_edges": args.max_edges,
            "final_by_geometry": final.to_dict(orient="records"),
        }, f, indent=2)

    make_figures(edge_df, summary, outdir)

    print("\n" + "=" * 110)
    print("SUMMARY BY GEOMETRY")
    print("=" * 110)
    print(summary.to_string(index=False))

    print("\n" + "=" * 110)
    print("FINAL BY GEOMETRY")
    print("=" * 110)
    print(final.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "ricci_curvature_rows.csv")
    print(outdir / "summary_by_geometry.csv")
    print(outdir / "final_by_geometry.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

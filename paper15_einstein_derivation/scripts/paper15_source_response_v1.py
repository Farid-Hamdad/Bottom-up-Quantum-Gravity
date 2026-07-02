#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Source response test v1

Goal
----
Third numerical test for Paper 15, aligned with Paper 8:

    δW_loc  ->  δκ_OR  ->  effective source response

We start from controlled geometries and a local weighted graph W_ij.
Then we inject a localized entanglement defect by increasing or decreasing
weights inside a small source region. We compute Ollivier--Ricci curvature
before/after and measure:

    Δκ_edge = κ_after - κ_before

as a function of distance from the source.

This is not yet a full tensor T_mu_nu derivation. It is a controlled
source-response test for the chain:

    local entanglement perturbation -> curvature response

Recommended run
---------------
cd ~/bottomup

python3 papers/paper15_einstein_derivation/scripts/paper15_source_response_v1.py \
  --geometries flat_torus2d sphere \
  --N 256 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --source-radius-factor 2.0 \
  --perturbation-strengths -0.30 -0.15 0.15 0.30 \
  --max-edges 700 \
  --alpha-idleness 0.5 \
  --output-dir papers/paper15_einstein_derivation/results/source_response_v1
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
    p = argparse.ArgumentParser(description="Paper 15 source response test v1.")
    p.add_argument("--geometries", nargs="+", default=["flat_torus2d", "sphere"])
    p.add_argument("--N", type=int, default=256)
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=64)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--source-radius-factor", type=float, default=2.0)
    p.add_argument("--perturbation-strengths", nargs="+", type=float, default=[-0.30, -0.15, 0.15, 0.30])
    p.add_argument("--alpha-idleness", type=float, default=0.5)
    p.add_argument("--max-edges", type=int, default=700)
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


def make_geometry(geometry, N):
    if geometry == "flat_torus2d":
        m = int(np.round(np.sqrt(N)))
        u = np.linspace(0.0, 1.0, m, endpoint=False)
        v = np.linspace(0.0, 1.0, m, endpoint=False)
        uu, vv = np.meshgrid(u, v)
        pts = np.column_stack([uu.ravel(), vv.ravel()])
        center = np.array([0.5, 0.5])
        return pts, center, {"m": m}

    if geometry == "sphere":
        i = np.arange(N)
        phi = np.arccos(1.0 - 2.0 * (i + 0.5) / N)
        golden = np.pi * (3.0 - np.sqrt(5.0))
        theta = golden * i
        xyz = np.column_stack([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi),
        ])
        center = np.array([0.0, 0.0, 1.0])
        return xyz, center, {}

    raise ValueError(f"Unknown geometry: {geometry}")


def intrinsic_distance_matrix(geometry, pts):
    if geometry == "flat_torus2d":
        du = np.abs(pts[:, None, 0] - pts[None, :, 0])
        dv = np.abs(pts[:, None, 1] - pts[None, :, 1])
        du = np.minimum(du, 1.0 - du)
        dv = np.minimum(dv, 1.0 - dv)
        return np.sqrt(du*du + dv*dv)

    if geometry == "sphere":
        dots = np.clip(pts @ pts.T, -1.0, 1.0)
        return np.arccos(dots)

    diff = pts[:, None, :] - pts[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def distance_to_center(geometry, pts, center):
    if geometry == "flat_torus2d":
        du = np.abs(pts[:, 0] - center[0])
        dv = np.abs(pts[:, 1] - center[1])
        du = np.minimum(du, 1.0 - du)
        dv = np.minimum(dv, 1.0 - dv)
        return np.sqrt(du*du + dv*dv)

    if geometry == "sphere":
        dots = np.clip(pts @ center, -1.0, 1.0)
        return np.arccos(dots)

    return np.sqrt(np.sum((pts - center) ** 2, axis=1))


def build_graph(D, k, kernel_factor):
    n = D.shape[0]
    Dw = D.copy()
    np.fill_diagonal(Dw, np.inf)

    kk = min(k, n - 2)
    knn = np.partition(Dw, kk, axis=1)[:, kk]
    eps = kernel_factor * float(np.median(knn[np.isfinite(knn)]) ** 2)
    eps = max(eps, 1e-14)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(Dw[i], kk)[:k]
        W[i, idx] = np.exp(-(Dw[i, idx] ** 2) / (4.0 * eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)

    A_len = np.where(W > 0, D, np.inf)
    np.fill_diagonal(A_len, 0.0)
    density = np.count_nonzero(W) / (n * (n - 1))
    return W, A_len, eps, density


def perturb_W(W, source_nodes, strength):
    Wp = W.copy()
    mask = np.zeros_like(W, dtype=bool)
    src = np.array(source_nodes, dtype=int)
    mask[np.ix_(src, src)] = True
    np.fill_diagonal(mask, False)

    # Multiplicative perturbation inside source region.
    factor = max(0.0, 1.0 + strength)
    Wp[mask] *= factor
    Wp = np.maximum(Wp, Wp.T)
    np.fill_diagonal(Wp, 0.0)
    return Wp


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
    support_i = np.where(mu > 1e-15)[0]
    support_j = np.where(nu > 1e-15)[0]
    a = mu[support_i]
    b = nu[support_j]
    C = dist_matrix[np.ix_(support_i, support_j)]

    m, n = C.shape
    c = C.ravel()

    A_eq = []
    b_eq = []

    for i in range(m):
        row = np.zeros(m*n)
        row[i*n:(i+1)*n] = 1.0
        A_eq.append(row)
        b_eq.append(a[i])

    for j in range(n):
        col = np.zeros(m*n)
        col[j::n] = 1.0
        A_eq.append(col)
        b_eq.append(b[j])

    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=(0, None), method="highs")
    return float(res.fun) if res.success else np.nan


def kappa_edge(W, spdist, i, j, alpha):
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


def edge_distance_to_source(i, j, node_dist_source):
    return float(min(node_dist_source[i], node_dist_source[j]))


def run_geometry_strength(geometry, strength, args):
    pts, center, meta = make_geometry(geometry, args.N)
    N_actual = len(pts)
    k = choose_k(N_actual, args)

    D = intrinsic_distance_matrix(geometry, pts)
    W, A_len, eps, density = build_graph(D, k, args.kernel_factor)

    node_dist = distance_to_center(geometry, pts, center)
    source_radius = args.source_radius_factor * np.sqrt(eps)
    source_nodes = np.where(node_dist <= source_radius)[0]
    if len(source_nodes) < 3:
        # Ensure non-empty region.
        source_nodes = np.argsort(node_dist)[:max(3, k // 2)]

    Wp = perturb_W(W, source_nodes, strength)

    sp0 = shortest_path(A_len, directed=False, unweighted=False)

    # Keep same edge-length graph topology for after perturbation, but shortest path still through same geometric edges.
    # This isolates curvature response due to probability redistribution, not topology changes.
    edges = sample_edges(W, args.max_edges, args.seed + N_actual + int(1000*abs(strength)) + len(geometry))

    rows = []
    for i, j in edges:
        i, j = int(i), int(j)
        k0 = kappa_edge(W, sp0, i, j, args.alpha_idleness)
        k1 = kappa_edge(Wp, sp0, i, j, args.alpha_idleness)
        rows.append({
            "geometry": geometry,
            "N_actual": N_actual,
            "k": k,
            "epsilon": eps,
            "density": density,
            "source_radius": source_radius,
            "n_source_nodes": int(len(source_nodes)),
            "perturbation_strength": strength,
            "i": i,
            "j": j,
            "edge_weight_before": float(W[i, j]),
            "edge_weight_after": float(Wp[i, j]),
            "edge_distance_to_source": edge_distance_to_source(i, j, node_dist),
            "edge_is_inside_source": bool(i in set(source_nodes) and j in set(source_nodes)),
            "kappa_before": k0,
            "kappa_after": k1,
            "delta_kappa": k1 - k0 if np.isfinite(k0) and np.isfinite(k1) else np.nan,
            "abs_delta_kappa": abs(k1 - k0) if np.isfinite(k0) and np.isfinite(k1) else np.nan,
        })

    return pd.DataFrame(rows)


def summarize(edge_df):
    summary = edge_df.groupby(["geometry", "perturbation_strength"]).agg(
        n_edges=("delta_kappa", "count"),
        N_actual=("N_actual", "first"),
        k=("k", "first"),
        epsilon=("epsilon", "first"),
        source_radius=("source_radius", "first"),
        n_source_nodes=("n_source_nodes", "first"),
        mean_delta_kappa=("delta_kappa", "mean"),
        median_delta_kappa=("delta_kappa", "median"),
        mean_abs_delta_kappa=("abs_delta_kappa", "mean"),
        median_abs_delta_kappa=("abs_delta_kappa", "median"),
        std_delta_kappa=("delta_kappa", "std"),
        near_mean_delta=("delta_kappa", lambda x: np.nan),
    ).reset_index()

    # Add near/far splits manually based on source radius multiples.
    extras = []
    for (geom, strength), sub in edge_df.groupby(["geometry", "perturbation_strength"]):
        r0 = float(sub.source_radius.iloc[0])
        near = sub[sub.edge_distance_to_source <= 1.5*r0]
        mid = sub[(sub.edge_distance_to_source > 1.5*r0) & (sub.edge_distance_to_source <= 3.0*r0)]
        far = sub[sub.edge_distance_to_source > 3.0*r0]

        extras.append({
            "geometry": geom,
            "perturbation_strength": strength,
            "near_n": len(near),
            "mid_n": len(mid),
            "far_n": len(far),
            "near_mean_delta": float(near.delta_kappa.mean()) if len(near) else np.nan,
            "mid_mean_delta": float(mid.delta_kappa.mean()) if len(mid) else np.nan,
            "far_mean_delta": float(far.delta_kappa.mean()) if len(far) else np.nan,
            "near_abs_mean_delta": float(near.abs_delta_kappa.mean()) if len(near) else np.nan,
            "mid_abs_mean_delta": float(mid.abs_delta_kappa.mean()) if len(mid) else np.nan,
            "far_abs_mean_delta": float(far.abs_delta_kappa.mean()) if len(far) else np.nan,
            "localization_ratio_near_far": (
                float(near.abs_delta_kappa.mean() / far.abs_delta_kappa.mean())
                if len(near) and len(far) and far.abs_delta_kappa.mean() > 0 else np.nan
            ),
        })
    extra_df = pd.DataFrame(extras)
    summary = summary.drop(columns=["near_mean_delta"]).merge(extra_df, on=["geometry", "perturbation_strength"], how="left")
    return summary


def make_figures(edge_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    for geom in sorted(edge_df.geometry.unique()):
        subg = edge_df[edge_df.geometry == geom]

        plt.figure(figsize=(8, 5))
        for strength in sorted(subg.perturbation_strength.unique()):
            sub = subg[subg.perturbation_strength == strength].copy()
            # bin by distance
            bins = np.linspace(0, sub.edge_distance_to_source.max(), 12)
            mids, vals = [], []
            for a, b in zip(bins[:-1], bins[1:]):
                s = sub[(sub.edge_distance_to_source >= a) & (sub.edge_distance_to_source < b)]
                if len(s):
                    mids.append(0.5*(a+b))
                    vals.append(s.abs_delta_kappa.mean())
            plt.plot(mids, vals, marker="o", label=f"strength={strength}")
        plt.xlabel("edge distance to source")
        plt.ylabel("mean |delta kappa|")
        plt.title(f"Curvature response localization — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_response_profile_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(8, 5))
        ssum = summary_df[summary_df.geometry == geom].sort_values("perturbation_strength")
        plt.plot(ssum.perturbation_strength, ssum.near_abs_mean_delta, marker="o", label="near")
        plt.plot(ssum.perturbation_strength, ssum.far_abs_mean_delta, marker="o", label="far")
        plt.xlabel("perturbation strength")
        plt.ylabel("mean |delta kappa|")
        plt.title(f"Near/far source response — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_near_far_vs_strength_{geom}.png", dpi=180)
        plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    print("=" * 110)
    print("Paper 15 — Source response test v1")
    print("=" * 110)

    for geom in args.geometries:
        for strength in args.perturbation_strengths:
            print(f"[run] geometry={geom} strength={strength}")
            df = run_geometry_strength(geom, strength, args)
            all_rows.append(df)
            print(
                f"      N={df.N_actual.iloc[0]} k={df.k.iloc[0]} "
                f"source_nodes={df.n_source_nodes.iloc[0]} "
                f"mean_delta={df.delta_kappa.mean():.5g} "
                f"mean_abs_delta={df.abs_delta_kappa.mean():.5g}"
            )

    edge_df = pd.concat(all_rows, ignore_index=True)
    edge_df.to_csv(outdir / "source_response_edges.csv", index=False)

    summary_df = summarize(edge_df)
    summary_df.to_csv(outdir / "source_response_summary.csv", index=False)

    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 15 source response test v1",
            "geometries": args.geometries,
            "N": args.N,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "source_radius_factor": args.source_radius_factor,
            "perturbation_strengths": args.perturbation_strengths,
            "alpha_idleness": args.alpha_idleness,
            "max_edges": args.max_edges,
            "summary": summary_df.to_dict(orient="records"),
        }, f, indent=2)

    make_figures(edge_df, summary_df, outdir)

    print("\n" + "=" * 110)
    print("SOURCE RESPONSE SUMMARY")
    print("=" * 110)
    print(summary_df.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "source_response_edges.csv")
    print(outdir / "source_response_summary.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

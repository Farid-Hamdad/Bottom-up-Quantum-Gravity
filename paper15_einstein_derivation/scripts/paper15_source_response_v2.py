#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 15 — Source response test v2

Goal
----
Improve v1 by replacing the ultra-local source-source perturbation with a
smooth radial entanglement source:

    phi_i = exp(-d(i,source)^2 / (2 sigma^2))

and

    W'_ij = W_ij * [1 + s * (phi_i + phi_j)/2]

This produces a controlled radial perturbation of the information network
and allows us to measure a curvature-response profile:

    Δκ(r) = κ_after(r) - κ_before(r)

Expected diagnostic
-------------------
1. Near-source response > far response.
2. Response amplitude grows with |s|.
3. Signed response should approximately flip when s changes sign.
4. The radial profile should decay with distance.

Recommended run
---------------
cd ~/bottomup

python3 papers/paper15_einstein_derivation/scripts/paper15_source_response_v2.py \
  --geometries flat_torus2d sphere \
  --N 256 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --sigma-factor 2.0 \
  --perturbation-strengths -0.30 -0.15 0.15 0.30 \
  --max-edges 900 \
  --alpha-idleness 0.5 \
  --output-dir papers/paper15_einstein_derivation/results/source_response_v2
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
from scipy.stats import spearmanr, pearsonr


def parse_args():
    p = argparse.ArgumentParser(description="Paper 15 source response test v2.")
    p.add_argument("--geometries", nargs="+", default=["flat_torus2d", "sphere"])
    p.add_argument("--N", type=int, default=256)
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=64)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--sigma-factor", type=float, default=2.0, help="sigma = sigma_factor * sqrt(epsilon)")
    p.add_argument("--perturbation-strengths", nargs="+", type=float, default=[-0.30, -0.15, 0.15, 0.30])
    p.add_argument("--alpha-idleness", type=float, default=0.5)
    p.add_argument("--max-edges", type=int, default=900)
    p.add_argument("--n-bins", type=int, default=10)
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
    return np.sqrt(np.sum(diff*diff, axis=-1))


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

    return np.sqrt(np.sum((pts-center)**2, axis=1))


def build_graph(D, k, kernel_factor):
    n = D.shape[0]
    Dw = D.copy()
    np.fill_diagonal(Dw, np.inf)

    kk = min(k, n-2)
    knn = np.partition(Dw, kk, axis=1)[:, kk]
    eps = kernel_factor * float(np.median(knn[np.isfinite(knn)])**2)
    eps = max(eps, 1e-14)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        idx = np.argpartition(Dw[i], kk)[:k]
        W[i, idx] = np.exp(-(Dw[i, idx]**2)/(4.0*eps))

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0.0)

    A_len = np.where(W > 0, D, np.inf)
    np.fill_diagonal(A_len, 0.0)

    density = np.count_nonzero(W)/(n*(n-1))
    return W, A_len, eps, density


def perturb_W_radial(W, phi, strength):
    """
    Smooth multiplicative perturbation on all existing links.
    factor_ij = 1 + strength * (phi_i + phi_j)/2
    clipped to positive values.
    """
    link_phi = 0.5*(phi[:, None] + phi[None, :])
    factor = 1.0 + strength * link_phi
    factor = np.maximum(factor, 1e-6)
    Wp = W * factor
    np.fill_diagonal(Wp, 0.0)
    Wp = np.maximum(Wp, Wp.T)
    return Wp, factor


def neighbor_distribution(W, node, alpha):
    n = W.shape[0]
    mu = np.zeros(n, dtype=float)
    mu[node] = alpha
    neigh = np.where(W[node] > 0)[0]
    total = W[node, neigh].sum()
    if total > 0:
        mu[neigh] += (1.0-alpha)*W[node, neigh]/total
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
    return float(1.0 - w1/dij)


def sample_edges(W, max_edges, seed):
    rng = np.random.default_rng(seed)
    ii, jj = np.where(np.triu(W > 0, k=1))
    edges = np.column_stack([ii, jj])
    if len(edges) > max_edges:
        idx = rng.choice(len(edges), size=max_edges, replace=False)
        edges = edges[idx]
    return edges


def edge_source_features(i, j, node_dist, phi):
    r_edge = float(min(node_dist[i], node_dist[j]))
    phi_edge = float(0.5*(phi[i] + phi[j]))
    return r_edge, phi_edge


def run_geometry_strength(geometry, strength, args):
    pts, center, meta = make_geometry(geometry, args.N)
    N_actual = len(pts)
    k = choose_k(N_actual, args)

    D = intrinsic_distance_matrix(geometry, pts)
    W, A_len, eps, density = build_graph(D, k, args.kernel_factor)
    node_dist = distance_to_center(geometry, pts, center)

    sigma = args.sigma_factor * np.sqrt(eps)
    sigma = max(sigma, 1e-12)
    phi = np.exp(-(node_dist**2)/(2.0*sigma**2))

    Wp, factor = perturb_W_radial(W, phi, strength)

    # Geometry distances are kept fixed; the source changes probability distributions, not the base geometry.
    spdist = shortest_path(A_len, directed=False, unweighted=False)
    edges = sample_edges(W, args.max_edges, args.seed + N_actual + int(1000*abs(strength)) + len(geometry))

    rows = []
    for i, j in edges:
        i, j = int(i), int(j)
        k0 = kappa_edge(W, spdist, i, j, args.alpha_idleness)
        k1 = kappa_edge(Wp, spdist, i, j, args.alpha_idleness)
        r_edge, phi_edge = edge_source_features(i, j, node_dist, phi)
        rows.append({
            "geometry": geometry,
            "N_actual": N_actual,
            "k": k,
            "epsilon": eps,
            "density": density,
            "sigma": sigma,
            "perturbation_strength": strength,
            "i": i,
            "j": j,
            "edge_distance_to_source": r_edge,
            "phi_edge": phi_edge,
            "edge_weight_before": float(W[i, j]),
            "edge_weight_after": float(Wp[i, j]),
            "edge_factor": float(factor[i, j]),
            "kappa_before": k0,
            "kappa_after": k1,
            "delta_kappa": k1-k0 if np.isfinite(k0) and np.isfinite(k1) else np.nan,
            "abs_delta_kappa": abs(k1-k0) if np.isfinite(k0) and np.isfinite(k1) else np.nan,
        })

    return pd.DataFrame(rows)


def binned_profile(edge_df, n_bins):
    rows = []
    for (geom, strength), sub in edge_df.groupby(["geometry", "perturbation_strength"]):
        max_r = sub.edge_distance_to_source.max()
        bins = np.linspace(0, max_r, n_bins+1)
        for b0, b1 in zip(bins[:-1], bins[1:]):
            s = sub[(sub.edge_distance_to_source >= b0) & (sub.edge_distance_to_source < b1)]
            if len(s) == 0:
                continue
            rows.append({
                "geometry": geom,
                "perturbation_strength": strength,
                "r_min": float(b0),
                "r_max": float(b1),
                "r_mid": float(0.5*(b0+b1)),
                "n_edges": int(len(s)),
                "mean_delta_kappa": float(s.delta_kappa.mean()),
                "median_delta_kappa": float(s.delta_kappa.median()),
                "mean_abs_delta_kappa": float(s.abs_delta_kappa.mean()),
                "mean_phi_edge": float(s.phi_edge.mean()),
            })
    return pd.DataFrame(rows)


def summarize(edge_df):
    rows = []
    for (geom, strength), sub in edge_df.groupby(["geometry", "perturbation_strength"]):
        sigma = float(sub.sigma.iloc[0])
        near = sub[sub.edge_distance_to_source <= sigma]
        mid = sub[(sub.edge_distance_to_source > sigma) & (sub.edge_distance_to_source <= 2.5*sigma)]
        far = sub[sub.edge_distance_to_source > 2.5*sigma]

        valid = sub.dropna(subset=["delta_kappa", "phi_edge"])
        if len(valid) > 3 and valid.phi_edge.std() > 0 and valid.delta_kappa.std() > 0:
            rho_s, p_s = spearmanr(valid.phi_edge, valid.delta_kappa)
            rho_abs, p_abs = spearmanr(valid.phi_edge, valid.abs_delta_kappa)
            r_p, p_p = pearsonr(valid.phi_edge, valid.delta_kappa)
        else:
            rho_s = p_s = rho_abs = p_abs = r_p = p_p = np.nan

        far_abs = float(far.abs_delta_kappa.mean()) if len(far) else np.nan
        near_abs = float(near.abs_delta_kappa.mean()) if len(near) else np.nan
        ratio = near_abs / far_abs if np.isfinite(far_abs) and far_abs > 0 else np.inf if np.isfinite(near_abs) and near_abs > 0 else np.nan

        rows.append({
            "geometry": geom,
            "perturbation_strength": strength,
            "N_actual": int(sub.N_actual.iloc[0]),
            "k": int(sub.k.iloc[0]),
            "epsilon": float(sub.epsilon.iloc[0]),
            "sigma": sigma,
            "n_edges": int(len(sub)),
            "mean_delta_kappa": float(sub.delta_kappa.mean()),
            "median_delta_kappa": float(sub.delta_kappa.median()),
            "mean_abs_delta_kappa": float(sub.abs_delta_kappa.mean()),
            "median_abs_delta_kappa": float(sub.abs_delta_kappa.median()),
            "near_n": int(len(near)),
            "mid_n": int(len(mid)),
            "far_n": int(len(far)),
            "near_mean_delta": float(near.delta_kappa.mean()) if len(near) else np.nan,
            "mid_mean_delta": float(mid.delta_kappa.mean()) if len(mid) else np.nan,
            "far_mean_delta": float(far.delta_kappa.mean()) if len(far) else np.nan,
            "near_abs_mean_delta": near_abs,
            "mid_abs_mean_delta": float(mid.abs_delta_kappa.mean()) if len(mid) else np.nan,
            "far_abs_mean_delta": far_abs,
            "localization_ratio_near_far": ratio,
            "spearman_phi_delta": float(rho_s) if np.isfinite(rho_s) else np.nan,
            "spearman_phi_delta_p": float(p_s) if np.isfinite(p_s) else np.nan,
            "spearman_phi_absdelta": float(rho_abs) if np.isfinite(rho_abs) else np.nan,
            "spearman_phi_absdelta_p": float(p_abs) if np.isfinite(p_abs) else np.nan,
            "pearson_phi_delta": float(r_p) if np.isfinite(r_p) else np.nan,
            "pearson_phi_delta_p": float(p_p) if np.isfinite(p_p) else np.nan,
        })
    return pd.DataFrame(rows)


def make_figures(edge_df, profile_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    for geom in sorted(edge_df.geometry.unique()):
        prof_g = profile_df[profile_df.geometry == geom]
        sum_g = summary_df[summary_df.geometry == geom]

        plt.figure(figsize=(8, 5))
        for strength in sorted(prof_g.perturbation_strength.unique()):
            s = prof_g[prof_g.perturbation_strength == strength].sort_values("r_mid")
            plt.plot(s.r_mid, s.mean_abs_delta_kappa, marker="o", label=f"s={strength}")
        plt.xlabel("edge distance to source")
        plt.ylabel("mean |delta kappa|")
        plt.title(f"Radial curvature response — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_radial_abs_response_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(8, 5))
        for strength in sorted(prof_g.perturbation_strength.unique()):
            s = prof_g[prof_g.perturbation_strength == strength].sort_values("r_mid")
            plt.plot(s.r_mid, s.mean_delta_kappa, marker="o", label=f"s={strength}")
        plt.axhline(0.0, linestyle="--")
        plt.xlabel("edge distance to source")
        plt.ylabel("mean delta kappa")
        plt.title(f"Signed radial curvature response — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_radial_signed_response_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(8, 5))
        plt.plot(sum_g.perturbation_strength, sum_g.near_abs_mean_delta, marker="o", label="near")
        plt.plot(sum_g.perturbation_strength, sum_g.mid_abs_mean_delta, marker="o", label="mid")
        plt.plot(sum_g.perturbation_strength, sum_g.far_abs_mean_delta, marker="o", label="far")
        plt.xlabel("perturbation strength")
        plt.ylabel("mean |delta kappa|")
        plt.title(f"Near/mid/far response — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_near_mid_far_vs_strength_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(8, 5))
        plt.plot(sum_g.perturbation_strength, sum_g.spearman_phi_absdelta, marker="o", label="Spearman(phi, |delta|)")
        plt.plot(sum_g.perturbation_strength, sum_g.spearman_phi_delta, marker="o", label="Spearman(phi, delta)")
        plt.axhline(0.0, linestyle="--")
        plt.xlabel("perturbation strength")
        plt.ylabel("correlation")
        plt.title(f"Source profile / curvature response correlation — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_phi_response_correlation_{geom}.png", dpi=180)
        plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_edges = []

    print("="*110)
    print("Paper 15 — Source response test v2")
    print("="*110)

    for geom in args.geometries:
        for strength in args.perturbation_strengths:
            print(f"[run] geometry={geom} strength={strength}")
            df = run_geometry_strength(geom, strength, args)
            all_edges.append(df)
            print(
                f"      N={df.N_actual.iloc[0]} k={df.k.iloc[0]} sigma={df.sigma.iloc[0]:.4g} "
                f"mean_delta={df.delta_kappa.mean():.5g} "
                f"mean_abs_delta={df.abs_delta_kappa.mean():.5g}"
            )

    edge_df = pd.concat(all_edges, ignore_index=True)
    summary_df = summarize(edge_df)
    profile_df = binned_profile(edge_df, args.n_bins)

    edge_df.to_csv(outdir / "source_response_v2_edges.csv", index=False)
    summary_df.to_csv(outdir / "source_response_v2_summary.csv", index=False)
    profile_df.to_csv(outdir / "source_response_v2_radial_profiles.csv", index=False)

    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 15 source response test v2",
            "geometries": args.geometries,
            "N": args.N,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "sigma_factor": args.sigma_factor,
            "perturbation_strengths": args.perturbation_strengths,
            "summary": summary_df.to_dict(orient="records"),
        }, f, indent=2)

    make_figures(edge_df, profile_df, summary_df, outdir)

    print("\n" + "="*110)
    print("SOURCE RESPONSE V2 SUMMARY")
    print("="*110)
    print(summary_df.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "source_response_v2_edges.csv")
    print(outdir / "source_response_v2_summary.csv")
    print(outdir / "source_response_v2_radial_profiles.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 18 — Modular source to curvature response v2

Tests the chain:
    δW_loc -> δ<K_A> -> δκ(r)

This script keeps the graph modular first-law diagnostic of v1, then measures
whether the modular response predicts the curvature response induced by the
same local entanglement perturbation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.optimize import linprog
from scipy.sparse.csgraph import shortest_path


def parse_args():
    p = argparse.ArgumentParser(description="Paper 18 modular source to curvature response v2.")
    p.add_argument("--geometries", nargs="+", default=["flat_torus2d", "sphere"])
    p.add_argument("--N", type=int, default=256)
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=80)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--region-radius-factor", type=float, default=2.5)
    p.add_argument("--sigma-factor", type=float, default=2.0)
    p.add_argument("--mu-factor", type=float, default=1.0)
    p.add_argument("--perturbation-strengths", nargs="+", type=float, default=[-0.20, -0.10, -0.05, 0.05, 0.10, 0.20])
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
        m = int(round(np.sqrt(N)))
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
        return np.sqrt(du * du + dv * dv)

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
        return np.sqrt(du * du + dv * dv)

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
    density = np.count_nonzero(W) / (n * (n - 1))

    A_len = np.where(W > 0, D, np.inf)
    np.fill_diagonal(A_len, 0.0)
    return W, A_len, eps, density


def perturb_W_radial(W, phi, strength):
    link_phi = 0.5 * (phi[:, None] + phi[None, :])
    factor = np.maximum(1.0 + strength * link_phi, 1e-8)
    Wp = np.maximum(W * factor, (W * factor).T)
    np.fill_diagonal(Wp, 0.0)
    return Wp


def laplacian(W):
    return np.diag(W.sum(axis=1)) - W


def density_from_laplacian(LA, mu):
    n = LA.shape[0]
    H = LA + mu * np.eye(n)
    evals, evecs = eigh(H)
    evals = np.maximum(evals, 1e-12)
    inv_evals = 1.0 / evals
    X = (evecs * inv_evals) @ evecs.T
    X = 0.5 * (X + X.T)
    rho = X / np.trace(X)
    return 0.5 * (rho + rho.T)


def entropy_and_K(rho):
    evals, evecs = eigh(rho)
    evals = np.maximum(evals, 1e-15)
    evals = evals / evals.sum()
    S = -float(np.sum(evals * np.log(evals)))
    K = -(evecs * np.log(evals)) @ evecs.T
    return S, 0.5 * (K + K.T)


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
    A_eq, b_eq = [], []

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


def compute_kappas_for_edges(W, spdist, edges, alpha):
    return np.array([kappa_edge(W, spdist, int(i), int(j), alpha) for i, j in edges], dtype=float)


def safe_mean(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if len(x) else np.nan


def run_geometry(geometry, args):
    pts, center, _ = make_geometry(geometry, args.N)
    N_actual = len(pts)
    k = choose_k(N_actual, args)
    D = intrinsic_distance_matrix(geometry, pts)
    W, A_len, eps, density = build_graph(D, k, args.kernel_factor)
    spdist = shortest_path(A_len, directed=False, unweighted=False)

    node_dist = distance_to_center(geometry, pts, center)
    sigma = args.sigma_factor * np.sqrt(eps)
    region_radius = args.region_radius_factor * np.sqrt(eps)
    region_nodes = np.where(node_dist <= region_radius)[0]
    if len(region_nodes) < max(8, k):
        region_nodes = np.argsort(node_dist)[:max(8, k)]

    phi = np.exp(-(node_dist**2) / (2.0 * sigma**2))
    L = laplacian(W)
    LA = L[np.ix_(region_nodes, region_nodes)]
    local_deg = W.sum(axis=1)[region_nodes]
    mu = args.mu_factor * max(float(np.mean(local_deg)), 1e-8)
    rho = density_from_laplacian(LA, mu)
    S0, K0 = entropy_and_K(rho)

    edges = sample_edges(W, args.max_edges, args.seed + N_actual + len(geometry))
    kappa0 = compute_kappas_for_edges(W, spdist, edges, args.alpha_idleness)
    edge_dist = np.array([min(node_dist[int(i)], node_dist[int(j)]) for i, j in edges], dtype=float)
    near_mask = edge_dist <= sigma
    mid_mask = (edge_dist > sigma) & (edge_dist <= 2.5 * sigma)
    far_mask = edge_dist > 2.5 * sigma

    rows, edge_rows = [], []
    for s in args.perturbation_strengths:
        Wp = perturb_W_radial(W, phi, s)
        Lp = laplacian(Wp)
        LAp = Lp[np.ix_(region_nodes, region_nodes)]
        rhop = density_from_laplacian(LAp, mu)
        S1, _ = entropy_and_K(rhop)
        delta_rho = rhop - rho
        delta_S = S1 - S0
        delta_K = float(np.trace(delta_rho @ K0))
        first_law_error = delta_S - delta_K
        first_law_relative_error = abs(first_law_error) / max(abs(delta_S), abs(delta_K), 1e-15)

        kappa1 = compute_kappas_for_edges(Wp, spdist, edges, args.alpha_idleness)
        delta_kappa = kappa1 - kappa0
        abs_delta_kappa = np.abs(delta_kappa)

        near_abs = safe_mean(abs_delta_kappa[near_mask])
        far_abs = safe_mean(abs_delta_kappa[far_mask])
        localization_ratio = near_abs / far_abs if np.isfinite(far_abs) and far_abs > 0 else np.inf

        rows.append({
            "geometry": geometry, "N_actual": N_actual, "k": k, "epsilon": eps,
            "density": density, "sigma": sigma, "region_radius": region_radius,
            "region_size": int(len(region_nodes)), "mu": mu, "n_edges": int(len(edges)),
            "near_n": int(np.sum(near_mask)), "mid_n": int(np.sum(mid_mask)), "far_n": int(np.sum(far_mask)),
            "perturbation_strength": float(s), "S0": S0, "S1": S1,
            "delta_S": delta_S, "delta_K_expectation": delta_K, "abs_delta_K": abs(delta_K),
            "first_law_error": first_law_error, "first_law_relative_error": first_law_relative_error,
            "mean_delta_kappa": safe_mean(delta_kappa), "mean_abs_delta_kappa": safe_mean(abs_delta_kappa),
            "near_mean_delta_kappa": safe_mean(delta_kappa[near_mask]),
            "mid_mean_delta_kappa": safe_mean(delta_kappa[mid_mask]),
            "far_mean_delta_kappa": safe_mean(delta_kappa[far_mask]),
            "near_abs_mean_delta_kappa": near_abs,
            "mid_abs_mean_delta_kappa": safe_mean(abs_delta_kappa[mid_mask]),
            "far_abs_mean_delta_kappa": far_abs,
            "localization_ratio_near_far": localization_ratio,
            "source_strength_region_mean": float(np.mean(phi[region_nodes])),
            "source_strength_region_sum": float(np.sum(phi[region_nodes])),
        })

        for idx, (i, j) in enumerate(edges):
            edge_rows.append({
                "geometry": geometry, "perturbation_strength": float(s),
                "i": int(i), "j": int(j), "edge_distance_to_source": float(edge_dist[idx]),
                "zone": "near" if near_mask[idx] else ("mid" if mid_mask[idx] else "far"),
                "kappa_before": float(kappa0[idx]) if np.isfinite(kappa0[idx]) else np.nan,
                "kappa_after": float(kappa1[idx]) if np.isfinite(kappa1[idx]) else np.nan,
                "delta_kappa": float(delta_kappa[idx]) if np.isfinite(delta_kappa[idx]) else np.nan,
                "abs_delta_kappa": float(abs_delta_kappa[idx]) if np.isfinite(abs_delta_kappa[idx]) else np.nan,
            })

    return pd.DataFrame(rows), pd.DataFrame(edge_rows)


def linear_fit(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 2 or np.std(x) <= 0 or np.std(y) <= 0:
        return {"intercept": np.nan, "slope": np.nan, "r2": np.nan, "pearson": np.nan}
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {"intercept": float(coef[0]), "slope": float(coef[1]),
            "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
            "pearson": float(np.corrcoef(x, y)[0, 1])}


def summarize(rows_df):
    parts = []
    for geom, sub in rows_df.groupby("geometry"):
        fit_first_law = linear_fit(sub["delta_K_expectation"], sub["delta_S"])
        fit_absK_near_abs = linear_fit(sub["abs_delta_K"], sub["near_abs_mean_delta_kappa"])
        fit_K_near_signed = linear_fit(sub["delta_K_expectation"], sub["near_mean_delta_kappa"])
        loc = sub.localization_ratio_near_far.replace([np.inf, -np.inf], np.nan)
        parts.append({
            "geometry": geom, "n_tests": int(len(sub)), "N_actual": int(sub.N_actual.iloc[0]),
            "k": int(sub.k.iloc[0]), "epsilon": float(sub.epsilon.iloc[0]),
            "region_size": int(sub.region_size.iloc[0]), "n_edges": int(sub.n_edges.iloc[0]),
            "mean_first_law_relative_error": float(sub.first_law_relative_error.mean()),
            "median_first_law_relative_error": float(sub.first_law_relative_error.median()),
            "max_first_law_relative_error": float(sub.first_law_relative_error.max()),
            "first_law_slope": fit_first_law["slope"], "first_law_r2": fit_first_law["r2"],
            "first_law_pearson": fit_first_law["pearson"],
            "mean_localization_ratio": float(loc.mean()), "median_localization_ratio": float(loc.median()),
            "absK_to_near_abs_slope": fit_absK_near_abs["slope"],
            "absK_to_near_abs_r2": fit_absK_near_abs["r2"],
            "absK_to_near_abs_pearson": fit_absK_near_abs["pearson"],
            "K_to_near_signed_slope": fit_K_near_signed["slope"],
            "K_to_near_signed_r2": fit_K_near_signed["r2"],
            "K_to_near_signed_pearson": fit_K_near_signed["pearson"],
        })
    return pd.DataFrame(parts)


def make_figures(rows_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)
    for geom in sorted(rows_df.geometry.unique()):
        sub = rows_df[rows_df.geometry == geom].sort_values("perturbation_strength")
        for name, ycols, labels, ylabel, title in [
            ("modular_response", ["delta_S", "delta_K_expectation"], ["delta S", "delta <K>"], "modular variation", "Modular response"),
            ("curvature_zones", ["near_abs_mean_delta_kappa", "mid_abs_mean_delta_kappa", "far_abs_mean_delta_kappa"], ["near", "mid", "far"], "mean |delta kappa|", "Curvature response zones"),
        ]:
            plt.figure(figsize=(7, 5))
            for col, lab in zip(ycols, labels):
                plt.plot(sub.perturbation_strength, sub[col], marker="o", label=lab)
            plt.axhline(0.0, linestyle="--")
            plt.xlabel("perturbation strength")
            plt.ylabel(ylabel)
            plt.title(f"{title} — {geom}")
            plt.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(figdir / f"fig_{name}_vs_strength_{geom}.png", dpi=180)
            plt.close()

        plt.figure(figsize=(7, 5))
        plt.scatter(sub.abs_delta_K, sub.near_abs_mean_delta_kappa)
        plt.xlabel("|delta <K>|")
        plt.ylabel("near mean |delta kappa|")
        plt.title(f"Modular source predicts curvature amplitude — {geom}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_absK_vs_near_abs_delta_kappa_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.scatter(sub.delta_K_expectation, sub.near_mean_delta_kappa)
        plt.axhline(0.0, linestyle="--")
        plt.axvline(0.0, linestyle="--")
        plt.xlabel("delta <K>")
        plt.ylabel("near mean delta kappa")
        plt.title(f"Signed modular source vs curvature response — {geom}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_K_vs_near_signed_delta_kappa_{geom}.png", dpi=180)
        plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows, all_edges = [], []
    print("=" * 110)
    print("Paper 18 — Modular source to curvature response v2")
    print("=" * 110)
    for geom in args.geometries:
        print(f"[run] geometry={geom}")
        rows, edge_rows = run_geometry(geom, args)
        all_rows.append(rows)
        all_edges.append(edge_rows)
        print(f"      N={rows.N_actual.iloc[0]} k={rows.k.iloc[0]} first_law_mean_rel={rows.first_law_relative_error.mean():.6g} localization_mean={rows.localization_ratio_near_far.replace([np.inf, -np.inf], np.nan).mean():.6g}")

    rows_df = pd.concat(all_rows, ignore_index=True)
    edges_df = pd.concat(all_edges, ignore_index=True)
    summary_df = summarize(rows_df)

    rows_df.to_csv(outdir / "modular_source_curvature_rows.csv", index=False)
    edges_df.to_csv(outdir / "modular_source_curvature_edges.csv", index=False)
    summary_df.to_csv(outdir / "modular_source_curvature_summary.csv", index=False)
    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({"experiment": "Paper 18 modular source to curvature response v2",
                   "geometries": args.geometries, "N": args.N, "k_mode": args.k_mode,
                   "kernel_factor": args.kernel_factor, "region_radius_factor": args.region_radius_factor,
                   "sigma_factor": args.sigma_factor, "mu_factor": args.mu_factor,
                   "perturbation_strengths": args.perturbation_strengths,
                   "summary": summary_df.to_dict(orient="records")}, f, indent=2)
    make_figures(rows_df, outdir)

    print("\n" + "=" * 110)
    print("MODULAR SOURCE CURVATURE ROWS")
    print("=" * 110)
    print(rows_df.to_string(index=False))
    print("\n" + "=" * 110)
    print("MODULAR SOURCE CURVATURE SUMMARY")
    print("=" * 110)
    print(summary_df.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "modular_source_curvature_rows.csv")
    print(outdir / "modular_source_curvature_edges.csv")
    print(outdir / "modular_source_curvature_summary.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

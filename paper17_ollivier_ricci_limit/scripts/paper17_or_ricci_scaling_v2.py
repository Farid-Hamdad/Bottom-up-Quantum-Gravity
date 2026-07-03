#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 17 — Ollivier-Ricci to Ricci scaling v2

Goal
----
Improve v1 by:

1. comparing flat torus and sphere by N_input rather than N_actual;
2. computing a calibrated Ricci proxy

    Rhat_OR = A_N * (kappa/epsilon - B_N)

where B_N is the flat-torus baseline and A_N is chosen so that the mean
sphere signal is 1 on the unit sphere:

    B_N = <kappa/epsilon>_flat
    A_N = 1 / ( <kappa/epsilon>_sphere - <kappa/epsilon>_flat )

This does NOT prove the continuum theorem. It identifies the finite-N affine
renormalization needed for

    kappa_OR / epsilon = B_N + C_N Ric(u,u) + o(1)

with Ric(u,u)=0 on flat torus and Ric(u,u)=1 on unit S^2.

Recommended run
---------------
cd ~/bottomup

python3 papers/paper17_ollivier_ricci_limit/scripts/paper17_or_ricci_scaling_v2.py \
  --geometries flat_torus2d sphere \
  --Ns 128 256 512 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --alpha-idleness 0.5 \
  --max-edges 900 \
  --output-dir papers/paper17_ollivier_ricci_limit/results/or_ricci_scaling_v2
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
    p = argparse.ArgumentParser(description="Paper 17 OR-to-Ricci scaling test v2.")
    p.add_argument("--geometries", nargs="+", default=["flat_torus2d", "sphere"])
    p.add_argument("--Ns", nargs="+", type=int, default=[128, 256, 512])
    p.add_argument("--k-mode", choices=["fixed", "sqrt", "log"], default="sqrt")
    p.add_argument("--k-fixed", type=int, default=16)
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=80)
    p.add_argument("--kernel-factor", type=float, default=0.5)
    p.add_argument("--alpha-idleness", type=float, default=0.5)
    p.add_argument("--max-edges", type=int, default=900)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def choose_k(N_actual, args):
    if args.k_mode == "fixed":
        k = args.k_fixed
    elif args.k_mode == "sqrt":
        k = int(round(np.sqrt(N_actual)))
    elif args.k_mode == "log":
        k = int(round(4 * np.log(max(N_actual, 3))))
    else:
        raise ValueError(args.k_mode)
    return int(max(args.k_min, min(args.k_max, k, N_actual - 2)))


def make_geometry(geometry, N_input):
    if geometry == "flat_torus2d":
        m = int(round(np.sqrt(N_input)))
        u = np.linspace(0.0, 1.0, m, endpoint=False)
        v = np.linspace(0.0, 1.0, m, endpoint=False)
        uu, vv = np.meshgrid(u, v)
        pts = np.column_stack([uu.ravel(), vv.ravel()])
        return pts, 2, 0.0, {"m": m}

    if geometry == "sphere":
        N = int(N_input)
        i = np.arange(N)
        phi = np.arccos(1.0 - 2.0 * (i + 0.5) / N)
        golden = np.pi * (3.0 - np.sqrt(5.0))
        theta = golden * i
        xyz = np.column_stack([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi),
        ])
        return xyz, 2, 1.0, {}

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


def run_case(geometry, N_input, args):
    pts, dim, ricci_target, meta = make_geometry(geometry, N_input)
    N_actual = len(pts)
    k = choose_k(N_actual, args)

    D = intrinsic_distance_matrix(geometry, pts)
    W, A_len, eps, density = build_graph(D, k, args.kernel_factor)
    spdist = shortest_path(A_len, directed=False, unweighted=False)
    edges = sample_edges(W, args.max_edges, args.seed + N_actual + len(geometry))

    rows = []
    for i, j in edges:
        i, j = int(i), int(j)
        ell = float(D[i, j])
        kappa = kappa_edge(W, spdist, i, j, args.alpha_idleness)
        rows.append({
            "geometry": geometry,
            "N_input": int(N_input),
            "N_actual": int(N_actual),
            "dim": dim,
            "ricci_target": ricci_target,
            "k": k,
            "epsilon": eps,
            "density": density,
            "alpha_idleness": args.alpha_idleness,
            "i": i,
            "j": j,
            "edge_length": ell,
            "edge_length2": ell * ell,
            "edge_weight": float(W[i, j]),
            "kappa_OR": kappa,
            "kappa_over_l2": kappa / (ell * ell) if np.isfinite(kappa) and ell > 0 else np.nan,
            "kappa_over_epsilon": kappa / eps if np.isfinite(kappa) else np.nan,
        })

    return pd.DataFrame(rows)


def summarize(edge_df):
    summary = edge_df.groupby(["geometry", "N_input"]).agg(
        n_edges=("kappa_OR", "count"),
        N_actual=("N_actual", "first"),
        dim=("dim", "first"),
        ricci_target=("ricci_target", "first"),
        k=("k", "first"),
        epsilon=("epsilon", "first"),
        density=("density", "first"),
        mean_edge_length=("edge_length", "mean"),
        median_edge_length=("edge_length", "median"),
        mean_kappa=("kappa_OR", "mean"),
        median_kappa=("kappa_OR", "median"),
        std_kappa=("kappa_OR", "std"),
        mean_kappa_over_l2=("kappa_over_l2", "mean"),
        median_kappa_over_l2=("kappa_over_l2", "median"),
        std_kappa_over_l2=("kappa_over_l2", "std"),
        mean_kappa_over_epsilon=("kappa_over_epsilon", "mean"),
        median_kappa_over_epsilon=("kappa_over_epsilon", "median"),
        std_kappa_over_epsilon=("kappa_over_epsilon", "std"),
        positive_fraction=("kappa_OR", lambda x: float((x > 0).mean())),
    ).reset_index()
    return summary


def calibrate_by_N(summary_df, edge_df):
    cal_rows = []
    edge_parts = []

    for N_input, sub in summary_df.groupby("N_input"):
        flat = sub[sub.geometry == "flat_torus2d"]
        sph = sub[sub.geometry == "sphere"]
        if not (len(flat) and len(sph)):
            continue

        f = flat.iloc[0]
        s = sph.iloc[0]

        B_eps = float(f.mean_kappa_over_epsilon)
        delta_eps = float(s.mean_kappa_over_epsilon - f.mean_kappa_over_epsilon)
        A_eps = 1.0 / delta_eps if abs(delta_eps) > 1e-15 else np.nan

        B_l2 = float(f.mean_kappa_over_l2)
        delta_l2 = float(s.mean_kappa_over_l2 - f.mean_kappa_over_l2)
        A_l2 = 1.0 / delta_l2 if abs(delta_l2) > 1e-15 else np.nan

        cal_rows.append({
            "N_input": int(N_input),
            "flat_N_actual": int(f.N_actual),
            "sphere_N_actual": int(s.N_actual),
            "flat_mean_kappa": float(f.mean_kappa),
            "sphere_mean_kappa": float(s.mean_kappa),
            "delta_mean_kappa": float(s.mean_kappa - f.mean_kappa),
            "B_epsilon_flat_baseline": B_eps,
            "sphere_mean_kappa_over_epsilon": float(s.mean_kappa_over_epsilon),
            "delta_kappa_over_epsilon": delta_eps,
            "A_epsilon": A_eps,
            "B_l2_flat_baseline": B_l2,
            "sphere_mean_kappa_over_l2": float(s.mean_kappa_over_l2),
            "delta_kappa_over_l2": delta_l2,
            "A_l2": A_l2,
        })

        part = edge_df[edge_df.N_input == N_input].copy()
        part["Rhat_epsilon"] = A_eps * (part["kappa_over_epsilon"] - B_eps)
        part["Rhat_l2"] = A_l2 * (part["kappa_over_l2"] - B_l2)
        edge_parts.append(part)

    cal_df = pd.DataFrame(cal_rows)
    if edge_parts:
        calibrated_edges = pd.concat(edge_parts, ignore_index=True)
    else:
        calibrated_edges = edge_df.copy()
        calibrated_edges["Rhat_epsilon"] = np.nan
        calibrated_edges["Rhat_l2"] = np.nan

    calibrated_summary = calibrated_edges.groupby(["geometry", "N_input"]).agg(
        n_edges=("kappa_OR", "count"),
        N_actual=("N_actual", "first"),
        ricci_target=("ricci_target", "first"),
        mean_Rhat_epsilon=("Rhat_epsilon", "mean"),
        median_Rhat_epsilon=("Rhat_epsilon", "median"),
        std_Rhat_epsilon=("Rhat_epsilon", "std"),
        mean_Rhat_l2=("Rhat_l2", "mean"),
        median_Rhat_l2=("Rhat_l2", "median"),
        std_Rhat_l2=("Rhat_l2", "std"),
    ).reset_index()

    calibrated_summary["error_mean_Rhat_epsilon"] = calibrated_summary["mean_Rhat_epsilon"] - calibrated_summary["ricci_target"]
    calibrated_summary["abs_error_mean_Rhat_epsilon"] = calibrated_summary["error_mean_Rhat_epsilon"].abs()
    calibrated_summary["error_mean_Rhat_l2"] = calibrated_summary["mean_Rhat_l2"] - calibrated_summary["ricci_target"]
    calibrated_summary["abs_error_mean_Rhat_l2"] = calibrated_summary["error_mean_Rhat_l2"].abs()

    return cal_df, calibrated_edges, calibrated_summary


def make_figures(summary_df, cal_df, calibrated_summary, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    # Raw and normalized signals.
    for col, ylabel, fname in [
        ("mean_kappa", "mean raw OR curvature", "fig_mean_raw_kappa_vs_Ninput.png"),
        ("mean_kappa_over_l2", "mean kappa / edge_length^2", "fig_mean_kappa_over_l2_vs_Ninput.png"),
        ("mean_kappa_over_epsilon", "mean kappa / epsilon", "fig_mean_kappa_over_epsilon_vs_Ninput.png"),
    ]:
        plt.figure(figsize=(7, 5))
        for geom in sorted(summary_df.geometry.unique()):
            s = summary_df[summary_df.geometry == geom].sort_values("N_input")
            plt.plot(s.N_input, s[col], marker="o", label=geom)
        plt.axhline(0.0, linestyle="--")
        plt.xscale("log", base=2)
        plt.xlabel("N input")
        plt.ylabel(ylabel)
        plt.title(ylabel + " vs N input")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / fname, dpi=180)
        plt.close()

    # Calibration constants.
    if len(cal_df):
        plt.figure(figsize=(7, 5))
        plt.plot(cal_df.N_input, cal_df.B_epsilon_flat_baseline, marker="o", label="B_epsilon flat baseline")
        plt.plot(cal_df.N_input, cal_df.delta_kappa_over_epsilon, marker="o", label="delta epsilon")
        plt.axhline(0.0, linestyle="--")
        plt.xscale("log", base=2)
        plt.xlabel("N input")
        plt.ylabel("value")
        plt.title("Affine calibration components")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / "fig_affine_calibration_components_vs_Ninput.png", dpi=180)
        plt.close()

    # Calibrated Ricci proxy.
    plt.figure(figsize=(7, 5))
    for geom in sorted(calibrated_summary.geometry.unique()):
        s = calibrated_summary[calibrated_summary.geometry == geom].sort_values("N_input")
        plt.plot(s.N_input, s.mean_Rhat_epsilon, marker="o", label=f"{geom} Rhat_epsilon")
    plt.axhline(0.0, linestyle="--", label="target flat")
    plt.axhline(1.0, linestyle=":", label="target sphere")
    plt.xscale("log", base=2)
    plt.xlabel("N input")
    plt.ylabel("mean calibrated Ricci proxy")
    plt.title("Calibrated Ricci proxy from kappa/epsilon")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(figdir / "fig_calibrated_Rhat_epsilon_vs_Ninput.png", dpi=180)
    plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_edges = []

    print("=" * 110)
    print("Paper 17 — Ollivier-Ricci to Ricci scaling v2")
    print("=" * 110)

    for geom in args.geometries:
        for N in args.Ns:
            print(f"[run] geometry={geom} N_input={N}")
            df = run_case(geom, N, args)
            all_edges.append(df)
            print(
                f"      N_actual={df.N_actual.iloc[0]} k={df.k.iloc[0]} "
                f"eps={df.epsilon.iloc[0]:.6g} "
                f"mean_kappa={df.kappa_OR.mean():.6g} "
                f"mean_kappa/eps={df.kappa_over_epsilon.mean():.6g}"
            )

    edge_df = pd.concat(all_edges, ignore_index=True)
    summary_df = summarize(edge_df)
    cal_df, calibrated_edges, calibrated_summary = calibrate_by_N(summary_df, edge_df)

    edge_df.to_csv(outdir / "or_ricci_scaling_v2_edges.csv", index=False)
    summary_df.to_csv(outdir / "or_ricci_scaling_v2_summary.csv", index=False)
    cal_df.to_csv(outdir / "or_ricci_scaling_v2_calibration.csv", index=False)
    calibrated_edges.to_csv(outdir / "or_ricci_scaling_v2_calibrated_edges.csv", index=False)
    calibrated_summary.to_csv(outdir / "or_ricci_scaling_v2_calibrated_summary.csv", index=False)

    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 17 OR-to-Ricci scaling v2",
            "geometries": args.geometries,
            "Ns": args.Ns,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "alpha_idleness": args.alpha_idleness,
            "max_edges": args.max_edges,
            "final_calibration": cal_df.sort_values("N_input").tail(1).to_dict(orient="records"),
            "final_calibrated_summary": calibrated_summary.sort_values("N_input").groupby("geometry").tail(1).to_dict(orient="records"),
        }, f, indent=2)

    make_figures(summary_df, cal_df, calibrated_summary, outdir)

    print("\n" + "=" * 110)
    print("SUMMARY BY GEOMETRY AND N_INPUT")
    print("=" * 110)
    print(summary_df.to_string(index=False))

    print("\n" + "=" * 110)
    print("AFFINE CALIBRATION BY N_INPUT")
    print("=" * 110)
    print(cal_df.to_string(index=False))

    print("\n" + "=" * 110)
    print("CALIBRATED RICCI PROXY SUMMARY")
    print("=" * 110)
    print(calibrated_summary.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "or_ricci_scaling_v2_edges.csv")
    print(outdir / "or_ricci_scaling_v2_summary.csv")
    print(outdir / "or_ricci_scaling_v2_calibration.csv")
    print(outdir / "or_ricci_scaling_v2_calibrated_edges.csv")
    print(outdir / "or_ricci_scaling_v2_calibrated_summary.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

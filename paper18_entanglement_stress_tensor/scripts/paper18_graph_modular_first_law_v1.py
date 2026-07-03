#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 18 — Graph modular first law v1

Goal
----
First numerical test for Paper 18:

    δW_loc  ->  δS_A^graph  ≈  δ<K_A^graph>

This is a controlled graph analogue of the modular first law:

    δS_A = δ<K_A>

Construction
------------
We build a controlled weighted graph on either:

    flat_torus2d
    sphere

Then choose a region A around a source center. From the restricted graph
W_A we construct a positive normalized matrix

    rho_A = (L_A + mu I)^(-1) / Tr[(L_A + mu I)^(-1)]

as a graph density proxy.

Then:

    S_A = -Tr rho_A log rho_A
    K_A = -log rho_A

We perturb W locally with a smooth radial source:

    phi_i = exp(-d(i,source)^2/(2 sigma^2))
    W'_ij = W_ij [1 + s (phi_i+phi_j)/2]

For the baseline rho_A and perturbed rho'_A, we compute:

    δS_A = S(rho'_A)-S(rho_A)
    δ<K_A> = Tr[(rho'_A-rho_A) K_A]

The modular first law should hold to first order in perturbation strength:

    δS_A ≈ δ<K_A>

Therefore this script scans small perturbations and reports:

    first_law_error = δS_A - δ<K_A>
    relative_error = |δS_A-δK| / max(|δS_A|, |δK|)

Recommended run
---------------
cd ~/bottomup

python3 papers/paper18_entanglement_stress_tensor/scripts/paper18_graph_modular_first_law_v1.py \
  --geometries flat_torus2d sphere \
  --N 256 \
  --k-mode sqrt \
  --kernel-factor 0.5 \
  --region-radius-factor 2.5 \
  --sigma-factor 2.0 \
  --mu-factor 1.0 \
  --perturbation-strengths -0.20 -0.10 -0.05 0.05 0.10 0.20 \
  --output-dir papers/paper18_entanglement_stress_tensor/results/graph_modular_first_law_v1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh


def parse_args():
    p = argparse.ArgumentParser(description="Paper 18 graph modular first law v1.")
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

    return np.sqrt(np.sum((pts - center)**2, axis=1))


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
    return W, eps, density


def perturb_W_radial(W, phi, strength):
    link_phi = 0.5 * (phi[:, None] + phi[None, :])
    factor = 1.0 + strength * link_phi
    factor = np.maximum(factor, 1e-8)
    Wp = W * factor
    Wp = np.maximum(Wp, Wp.T)
    np.fill_diagonal(Wp, 0.0)
    return Wp


def laplacian(W):
    deg = W.sum(axis=1)
    return np.diag(deg) - W


def restrict_matrix(M, nodes):
    return M[np.ix_(nodes, nodes)]


def density_from_laplacian(LA, mu):
    """
    rho = (LA + mu I)^(-1) / Tr(...)
    computed spectrally for stability.
    """
    n = LA.shape[0]
    H = LA + mu * np.eye(n)
    evals, evecs = eigh(H)
    evals = np.maximum(evals, 1e-12)
    inv_evals = 1.0 / evals
    X = (evecs * inv_evals) @ evecs.T
    X = 0.5 * (X + X.T)
    tr = np.trace(X)
    rho = X / tr
    rho = 0.5 * (rho + rho.T)
    return rho


def entropy_and_K(rho):
    evals, evecs = eigh(rho)
    evals = np.maximum(evals, 1e-15)
    evals = evals / evals.sum()

    S = -float(np.sum(evals * np.log(evals)))
    log_evals = np.log(evals)
    K = -(evecs * log_evals) @ evecs.T
    K = 0.5 * (K + K.T)
    return S, K, evals


def run_geometry(geometry, args):
    pts, center, meta = make_geometry(geometry, args.N)
    N_actual = len(pts)
    k = choose_k(N_actual, args)

    D = intrinsic_distance_matrix(geometry, pts)
    W, eps, density = build_graph(D, k, args.kernel_factor)
    node_dist = distance_to_center(geometry, pts, center)

    sigma = args.sigma_factor * np.sqrt(eps)
    region_radius = args.region_radius_factor * np.sqrt(eps)

    region_nodes = np.where(node_dist <= region_radius)[0]
    if len(region_nodes) < max(8, k):
        region_nodes = np.argsort(node_dist)[:max(8, k)]

    phi = np.exp(-(node_dist**2) / (2.0 * sigma**2))

    L = laplacian(W)
    LA = restrict_matrix(L, region_nodes)

    # mu sets the IR regularization of the graph density proxy.
    # Use local mean degree times mu_factor for scale adaptation.
    local_deg = W.sum(axis=1)[region_nodes]
    mu = args.mu_factor * max(float(np.mean(local_deg)), 1e-8)

    rho = density_from_laplacian(LA, mu)
    S0, K0, evals0 = entropy_and_K(rho)

    rows = []
    for s in args.perturbation_strengths:
        Wp = perturb_W_radial(W, phi, s)
        Lp = laplacian(Wp)
        LAp = restrict_matrix(Lp, region_nodes)

        rhop = density_from_laplacian(LAp, mu)
        S1, K1_unused, evals1 = entropy_and_K(rhop)

        delta_rho = rhop - rho
        delta_S = S1 - S0

        # Modular first law uses the unperturbed modular Hamiltonian K0.
        delta_K_exp = float(np.trace(delta_rho @ K0))

        first_law_error = delta_S - delta_K_exp
        denom = max(abs(delta_S), abs(delta_K_exp), 1e-15)
        relative_error = abs(first_law_error) / denom

        source_strength_region = float(np.mean(phi[region_nodes]))
        source_total_region = float(np.sum(phi[region_nodes]))

        rows.append({
            "geometry": geometry,
            "N_actual": N_actual,
            "k": k,
            "epsilon": eps,
            "density": density,
            "sigma": sigma,
            "region_radius": region_radius,
            "region_size": int(len(region_nodes)),
            "mu": mu,
            "perturbation_strength": float(s),
            "source_strength_region_mean": source_strength_region,
            "source_strength_region_sum": source_total_region,
            "S0": S0,
            "S1": S1,
            "delta_S": delta_S,
            "delta_K_expectation": delta_K_exp,
            "first_law_error": first_law_error,
            "relative_error": relative_error,
            "abs_delta_S": abs(delta_S),
            "abs_delta_K": abs(delta_K_exp),
            "rho_min_eig": float(np.min(evals0)),
            "rho_max_eig": float(np.max(evals0)),
            "rho_entropy_effective_rank": float(np.exp(S0)),
        })

    return pd.DataFrame(rows)


def summarize(rows_df):
    summary = rows_df.groupby("geometry").agg(
        n_tests=("perturbation_strength", "count"),
        N_actual=("N_actual", "first"),
        k=("k", "first"),
        epsilon=("epsilon", "first"),
        region_size=("region_size", "first"),
        mu=("mu", "first"),
        mean_abs_delta_S=("abs_delta_S", "mean"),
        mean_abs_delta_K=("abs_delta_K", "mean"),
        mean_abs_first_law_error=("first_law_error", lambda x: float(np.mean(np.abs(x)))),
        median_abs_first_law_error=("first_law_error", lambda x: float(np.median(np.abs(x)))),
        mean_relative_error=("relative_error", "mean"),
        median_relative_error=("relative_error", "median"),
        max_relative_error=("relative_error", "max"),
    ).reset_index()

    # Linear fit delta_S = a + b deltaK for each geometry
    fits = []
    for geom, sub in rows_df.groupby("geometry"):
        x = sub["delta_K_expectation"].to_numpy()
        y = sub["delta_S"].to_numpy()
        if len(x) >= 2 and np.std(x) > 0:
            A = np.column_stack([np.ones_like(x), x])
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            yhat = A @ coef
            ss_res = float(np.sum((y - yhat) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            corr = float(np.corrcoef(x, y)[0, 1])
            fits.append({
                "geometry": geom,
                "fit_intercept_deltaS_vs_deltaK": float(coef[0]),
                "fit_slope_deltaS_vs_deltaK": float(coef[1]),
                "fit_r2_deltaS_vs_deltaK": r2,
                "pearson_deltaS_deltaK": corr,
            })
    fit_df = pd.DataFrame(fits)
    out = summary.merge(fit_df, on="geometry", how="left")
    return out


def make_figures(rows_df, summary_df, outdir):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    for geom in sorted(rows_df.geometry.unique()):
        sub = rows_df[rows_df.geometry == geom].sort_values("perturbation_strength")

        plt.figure(figsize=(7, 5))
        plt.plot(sub.perturbation_strength, sub.delta_S, marker="o", label="delta S")
        plt.plot(sub.perturbation_strength, sub.delta_K_expectation, marker="x", label="delta <K>")
        plt.axhline(0.0, linestyle="--")
        plt.xlabel("perturbation strength")
        plt.ylabel("variation")
        plt.title(f"Graph modular first law — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_deltaS_deltaK_vs_strength_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.scatter(sub.delta_K_expectation, sub.delta_S)
        lo = min(sub.delta_K_expectation.min(), sub.delta_S.min())
        hi = max(sub.delta_K_expectation.max(), sub.delta_S.max())
        plt.plot([lo, hi], [lo, hi], linestyle="--", label="y=x")
        plt.xlabel("delta <K>")
        plt.ylabel("delta S")
        plt.title(f"delta S vs delta <K> — {geom}")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figdir / f"fig_deltaS_vs_deltaK_{geom}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.plot(sub.perturbation_strength, sub.relative_error, marker="o")
        plt.xlabel("perturbation strength")
        plt.ylabel("relative first-law error")
        plt.title(f"First-law relative error — {geom}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_first_law_relative_error_{geom}.png", dpi=180)
        plt.close()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    print("=" * 110)
    print("Paper 18 — Graph modular first law v1")
    print("=" * 110)

    for geom in args.geometries:
        print(f"[run] geometry={geom}")
        df = run_geometry(geom, args)
        all_rows.append(df)
        print(
            f"      N={df.N_actual.iloc[0]} k={df.k.iloc[0]} region={df.region_size.iloc[0]} "
            f"mean_rel_error={df.relative_error.mean():.6g}"
        )

    rows_df = pd.concat(all_rows, ignore_index=True)
    summary_df = summarize(rows_df)

    rows_df.to_csv(outdir / "graph_modular_first_law_rows.csv", index=False)
    summary_df.to_csv(outdir / "graph_modular_first_law_summary.csv", index=False)

    with open(outdir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "Paper 18 graph modular first law v1",
            "geometries": args.geometries,
            "N": args.N,
            "k_mode": args.k_mode,
            "kernel_factor": args.kernel_factor,
            "region_radius_factor": args.region_radius_factor,
            "sigma_factor": args.sigma_factor,
            "mu_factor": args.mu_factor,
            "perturbation_strengths": args.perturbation_strengths,
            "summary": summary_df.to_dict(orient="records"),
        }, f, indent=2)

    make_figures(rows_df, summary_df, outdir)

    print("\n" + "=" * 110)
    print("MODULAR FIRST LAW ROWS")
    print("=" * 110)
    print(rows_df.to_string(index=False))

    print("\n" + "=" * 110)
    print("MODULAR FIRST LAW SUMMARY")
    print("=" * 110)
    print(summary_df.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "graph_modular_first_law_rows.csv")
    print(outdir / "graph_modular_first_law_summary.csv")
    print(outdir / "summary.json")
    print(outdir / "figures")


if __name__ == "__main__":
    main()

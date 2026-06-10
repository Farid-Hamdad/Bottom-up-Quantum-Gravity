#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 11 candidate — Scaling of walk dimension d_w(N).

Goal:
    Test whether d_w(N) tends to 2 or to a non-brownian value > 2.

For each N:
    - generate a geometric MI matrix in 3D
    - build a kNN entanglement graph
    - compute spectral dimension d_s from heat trace
    - compute walk dimension d_w from MSD:
          <r^2(tau)> ~ tau^(2/d_w)
    - compute candidate alpha laws:
          alpha_simple = d_s - 2
          alpha_standard = (d_w/2)(d_s - 2)
          alpha_bup = 2*d_s/d_w + d_w - 4

Outputs:
    - dw_scaling_summary.csv
    - summary.json
    - figures/
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.sparse.csgraph import shortest_path
from scipy.stats import linregress


# ============================================================
# Utilities
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def auto_k(N, factor=1.30):
    return max(4, int(round(factor * (N ** (1.0 / 3.0)))))


def sample_points(N, dim=3, seed=0, periodic=True):
    rng = np.random.default_rng(seed)
    return rng.random((N, dim))


def pairwise_periodic_distances(X):
    N, dim = X.shape
    D = np.empty((N, N), dtype=float)

    for i in range(N):
        diff = np.abs(X - X[i])
        diff = np.minimum(diff, 1.0 - diff)
        D[i] = np.sqrt(np.sum(diff * diff, axis=1))

    return D


def generate_mi_matrix(N, lam=0.30, power=1.0, seed=0):
    X = sample_points(N, dim=3, seed=seed, periodic=True)
    R = pairwise_periodic_distances(X)

    W = np.exp(-np.power(R / max(lam, 1e-12), power))
    np.fill_diagonal(W, 0.0)

    W = 0.5 * (W + W.T)
    W[W < 0] = 0.0

    wmax = np.max(W)
    if wmax > 0:
        W /= wmax

    return W


def knn_graph(W, k):
    N = W.shape[0]
    mask = np.zeros((N, N), dtype=bool)

    for i in range(N):
        idx = np.argsort(W[i])[::-1]
        idx = [j for j in idx if j != i and W[i, j] > 0]
        for j in idx[:k]:
            mask[i, j] = True

    mask = np.logical_or(mask, mask.T)
    np.fill_diagonal(mask, False)

    A = np.where(mask, W, 0.0)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)

    deg = A.sum(axis=1)
    L = sparse.diags(deg) - sparse.csr_matrix(A)

    return A, L.tocsr(), mask, deg


def entanglement_distance(A, mask, eps=1e-12):
    Wmax = A[mask].max()
    lengths = np.full_like(A, np.inf, dtype=float)

    lengths[mask] = -np.log((A[mask] + eps) / (Wmax + eps))
    lengths[mask] = np.maximum(lengths[mask], 1e-9)
    np.fill_diagonal(lengths, 0.0)

    D = shortest_path(lengths, method="FW", directed=False)
    D[np.isinf(D)] = np.nan
    return D


# ============================================================
# Spectral dimension
# ============================================================

def spectral_dimension(L, tau_min, tau_max, tau_points, eig_k, dense_threshold=500):
    N = L.shape[0]
    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_points)

    if N <= dense_threshold:
        evals = np.linalg.eigvalsh(L.toarray())
    else:
        k = min(eig_k, N - 2)
        evals = eigsh(L, k=k, which="SM", return_eigenvectors=False)
        evals = np.sort(np.maximum(np.real(evals), 0.0))

    Z = np.array([np.sum(np.exp(-t * evals)) for t in taus], dtype=float)

    logt = np.log(taus)
    logZ = np.log(np.maximum(Z, 1e-300))
    slope = np.gradient(logZ, logt)
    ds_curve = -2.0 * slope

    n = len(ds_curve)
    lo = max(1, n // 4)
    hi = min(n - 1, 3 * n // 4)
    ds_eff = float(np.nanmedian(ds_curve[lo:hi]))

    return taus, Z, ds_curve, ds_eff, evals


# ============================================================
# Walk dimension
# ============================================================

def build_center_matrix(N, centers):
    B = np.zeros((N, len(centers)), dtype=float)
    for col, c in enumerate(centers):
        B[c, col] = 1.0
    return B


def msd_for_tau(L, D, centers, tau):
    """
    Compute heat kernel columns K(:, centers) using expm_multiply.
    Then compute <r^2(tau)> for each center.
    """
    N = L.shape[0]
    B = build_center_matrix(N, centers)

    # Kcols[j, a] = K(tau, j, center_a)
    Kcols = expm_multiply((-tau) * L, B)

    vals = []

    for col, c in enumerate(centers):
        weights = np.asarray(Kcols[:, col]).ravel()
        r = D[c]

        m = np.isfinite(r) & np.isfinite(weights) & (r > 0) & (weights > 0)

        if m.sum() < 5:
            continue

        denom = weights[m].sum()
        if denom <= 0:
            continue

        m2 = float(np.sum(weights[m] * (r[m] ** 2)) / denom)
        vals.append(m2)

    if len(vals) == 0:
        return np.nan, np.nan, 0

    vals = np.asarray(vals, dtype=float)
    return float(np.mean(vals)), float(np.std(vals)), int(len(vals))


def fit_walk_dimension(walk_df):
    sub = walk_df.copy()
    sub = sub[np.isfinite(sub["m2_mean"]) & (sub["m2_mean"] > 0)]
    sub = sub[np.isfinite(sub["tau"]) & (sub["tau"] > 0)]

    if len(sub) < 4:
        return None

    x = np.log(sub["tau"].values)
    y = np.log(sub["m2_mean"].values)

    reg = linregress(x, y)
    beta = float(reg.slope)
    dw = float(2.0 / beta) if beta != 0 else np.nan

    return {
        "n_tau_walk": int(len(sub)),
        "beta_msd": beta,
        "d_w": dw,
        "walk_r2": float(reg.rvalue ** 2),
        "walk_pvalue": float(reg.pvalue),
    }


# ============================================================
# One N run
# ============================================================

def run_one_N(args, N, seed):
    k = args.k if args.k > 0 else auto_k(N, args.k_auto_factor)

    print("=" * 100)
    print(f"Running N={N}, seed={seed}, lambda={args.lam}, k={k}")
    print("=" * 100)

    W = generate_mi_matrix(
        N=N,
        lam=args.lam,
        power=args.power,
        seed=seed
    )

    A, L, mask, deg = knn_graph(W, k)
    D = entanglement_distance(A, mask)

    taus, Z, ds_curve, ds_eff, evals = spectral_dimension(
        L,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        tau_points=args.tau_points,
        eig_k=args.eig_k,
        dense_threshold=args.dense_threshold
    )

    rng = np.random.default_rng(args.seed_centers + 1000 * seed + N)
    n_centers = min(args.n_centers, N)
    centers = rng.choice(N, size=n_centers, replace=False)

    tau_walk = np.logspace(
        np.log10(args.walk_tau_min),
        np.log10(args.walk_tau_max),
        args.walk_tau_points
    )

    walk_rows = []
    for idx, tau in enumerate(tau_walk):
        print(f"  [{idx+1}/{len(tau_walk)}] tau={tau:.5g}")
        m2_mean, m2_std, n_used = msd_for_tau(L, D, centers, tau)
        walk_rows.append({
            "N": N,
            "seed": seed,
            "k": k,
            "lambda": args.lam,
            "tau": float(tau),
            "m2_mean": m2_mean,
            "m2_std": m2_std,
            "n_centers_used": n_used
        })

    walk_df = pd.DataFrame(walk_rows)
    walk_fit = fit_walk_dimension(walk_df)

    if walk_fit is None:
        walk_fit = {
            "n_tau_walk": 0,
            "beta_msd": np.nan,
            "d_w": np.nan,
            "walk_r2": np.nan,
            "walk_pvalue": np.nan,
        }

    dw = walk_fit["d_w"]

    alpha_simple = ds_eff - 2.0
    alpha_standard = 0.5 * dw * (ds_eff - 2.0) if np.isfinite(dw) else np.nan
    alpha_bup = (2.0 * ds_eff / dw) + dw - 4.0 if np.isfinite(dw) and dw != 0 else np.nan

    row = {
        "N": N,
        "seed": seed,
        "lambda": args.lam,
        "k": k,
        "edges": int(mask.sum() // 2),
        "mean_degree": float(np.mean(deg)),
        "min_degree": float(np.min(deg)),
        "max_degree": float(np.max(deg)),

        "d_s": ds_eff,
        "d_w": dw,
        "beta_msd": walk_fit["beta_msd"],
        "walk_r2": walk_fit["walk_r2"],
        "walk_pvalue": walk_fit["walk_pvalue"],
        "n_centers": n_centers,

        "alpha_simple_ds_minus_2": alpha_simple,
        "alpha_standard_dw_over_2": alpha_standard,
        "alpha_bup_effective": alpha_bup,

        "lambda_2": float(np.sort(evals)[1]) if len(evals) > 1 else np.nan,
        "lambda_max_used": float(np.max(evals)),
        "eig_count": int(len(evals)),
    }

    return row, walk_df, pd.DataFrame({
        "N": N,
        "seed": seed,
        "lambda": args.lam,
        "k": k,
        "tau": taus,
        "Z": Z,
        "d_s_tau": ds_curve
    })


# ============================================================
# Figures
# ============================================================

def make_figures(summary_df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    # d_w vs N
    plt.figure(figsize=(8, 5))
    for seed in sorted(summary_df["seed"].unique()):
        sub = summary_df[summary_df["seed"] == seed].sort_values("N")
        plt.plot(sub["N"], sub["d_w"], "o-", label=f"seed={seed}")

    plt.axhline(2.0, linestyle="--", linewidth=1.5, label=r"Brownian $d_w=2$")
    plt.xlabel("N")
    plt.ylabel(r"$d_w$")
    plt.title("BuP Paper 11 — Walk dimension scaling")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_dw_vs_N.png"), dpi=250)
    plt.close()

    # d_s vs N
    plt.figure(figsize=(8, 5))
    for seed in sorted(summary_df["seed"].unique()):
        sub = summary_df[summary_df["seed"] == seed].sort_values("N")
        plt.plot(sub["N"], sub["d_s"], "o-", label=f"seed={seed}")

    plt.axhline(3.0, linestyle="--", linewidth=1.5, label=r"$d_s=3$")
    plt.xlabel("N")
    plt.ylabel(r"$d_s$")
    plt.title("BuP Paper 11 — Spectral dimension scaling")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_ds_vs_N.png"), dpi=250)
    plt.close()

    # alpha predictions vs N
    plt.figure(figsize=(9, 5))
    for seed in sorted(summary_df["seed"].unique()):
        sub = summary_df[summary_df["seed"] == seed].sort_values("N")
        plt.plot(sub["N"], sub["alpha_simple_ds_minus_2"], "o--", label=f"simple seed={seed}")
        plt.plot(sub["N"], sub["alpha_bup_effective"], "s-", label=f"BuP eff seed={seed}")

    plt.axhline(1.0, linestyle=":", linewidth=1.5, label=r"Newton $\alpha=1$")
    plt.xlabel("N")
    plt.ylabel(r"predicted $\alpha$")
    plt.title("BuP Paper 11 — Alpha prediction scaling")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_alpha_predictions_vs_N.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--output-dir", required=True)

    parser.add_argument("--N-list", nargs="+", type=int, default=[300, 500, 800, 1000])
    parser.add_argument("--seed-list", nargs="+", type=int, default=[0])
    parser.add_argument("--lam", type=float, default=0.30)
    parser.add_argument("--power", type=float, default=1.0)

    parser.add_argument("--k", type=int, default=-1, help="Use fixed k if >0, otherwise auto k.")
    parser.add_argument("--k-auto-factor", type=float, default=1.30)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=100.0)
    parser.add_argument("--tau-points", type=int, default=50)
    parser.add_argument("--eig-k", type=int, default=600)
    parser.add_argument("--dense-threshold", type=int, default=500)

    parser.add_argument("--walk-tau-min", type=float, default=0.05)
    parser.add_argument("--walk-tau-max", type=float, default=20.0)
    parser.add_argument("--walk-tau-points", type=int, default=20)
    parser.add_argument("--n-centers", type=int, default=80)
    parser.add_argument("--seed-centers", type=int, default=0)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    all_summary = []
    all_walk = []
    all_ds = []

    for seed in args.seed_list:
        for N in args.N_list:
            row, walk_df, ds_df = run_one_N(args, N, seed)
            all_summary.append(row)
            all_walk.append(walk_df)
            all_ds.append(ds_df)

    summary_df = pd.DataFrame(all_summary)
    walk_all = pd.concat(all_walk, ignore_index=True)
    ds_all = pd.concat(all_ds, ignore_index=True)

    summary_path = os.path.join(args.output_dir, "dw_scaling_summary.csv")
    walk_path = os.path.join(args.output_dir, "walk_msd_all.csv")
    ds_path = os.path.join(args.output_dir, "spectral_dimension_all.csv")

    summary_df.to_csv(summary_path, index=False)
    walk_all.to_csv(walk_path, index=False)
    ds_all.to_csv(ds_path, index=False)

    make_figures(summary_df, args.output_dir)

    result = {
        "title": "BuP Paper 11 candidate — d_w(N) scaling",
        "description": "Tests whether walk dimension tends to 2 or remains >2.",
        "parameters": vars(args),
        "files": {
            "summary_csv": "dw_scaling_summary.csv",
            "walk_msd_all": "walk_msd_all.csv",
            "spectral_dimension_all": "spectral_dimension_all.csv",
            "figures_dir": "figures/"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(result, f, indent=2)

    print("=" * 100)
    print("BuP Paper 11 — d_w(N) scaling summary")
    print("=" * 100)
    print(summary_df.to_string(index=False))
    print("-" * 100)
    print("Output:", args.output_dir)
    print("DONE")


if __name__ == "__main__":
    main()

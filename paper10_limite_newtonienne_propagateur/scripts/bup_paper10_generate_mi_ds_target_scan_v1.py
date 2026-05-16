#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 10 — Generate MI-like matrices and scan for d_s ≈ 3.

Goal:
    Generate controlled BuP-like mutual-information matrices W_ij,
    construct kNN entanglement graphs, compute spectral dimension d_s,
    and search for regimes where:

        d_s ≈ 3.

This is intended for Paper 10 / Paper 11:
    - find MI graph regimes close to a 3D Newtonian limit
    - export MI_N*_lam*.csv matrices for later Green-propagator tests

Model:
    Points x_i are sampled in an ambient dimension d_amb.
    Pairwise MI weights are generated as:

        W_ij = exp(-(r_ij / lambda)^power)

    with optional weak noise and optional long-range background.

Outputs:
    output_dir/
        ds_scan_summary.csv
        summary.json
        mi_matrices/
            MI_N*_d*_lam*_k*_seed*.csv
        figures/
            fig_ds_vs_lambda_by_k.png
            fig_best_ds_target.png
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse.linalg import eigsh


# ============================================================
# Utilities
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def clean_matrix(W):
    W = np.asarray(W, dtype=float)
    W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    W[W < 0] = 0.0
    return W


# ============================================================
# Geometry and MI generation
# ============================================================

def sample_points(N, dim, seed=0, mode="cube"):
    rng = np.random.default_rng(seed)

    if mode == "cube":
        X = rng.random((N, dim))

    elif mode == "gaussian":
        X = rng.normal(size=(N, dim))
        X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-12)

    elif mode == "sphere":
        X = rng.normal(size=(N, dim))
        X /= np.linalg.norm(X, axis=1, keepdims=True) + 1e-12
        X = 0.5 + 0.5 * X

    else:
        raise ValueError(f"Unknown point sampling mode: {mode}")

    return X


def pairwise_distances(X, periodic=False):
    N, dim = X.shape
    D = np.zeros((N, N), dtype=float)

    for i in range(N):
        diff = np.abs(X - X[i])

        if periodic:
            diff = np.minimum(diff, 1.0 - diff)

        D[i] = np.sqrt(np.sum(diff ** 2, axis=1))

    return D


def generate_mi_matrix(
    N,
    dim,
    lam,
    power=1.0,
    seed=0,
    noise=0.0,
    background=0.0,
    periodic=False,
    point_mode="cube"
):
    rng = np.random.default_rng(seed)

    X = sample_points(N, dim, seed=seed, mode=point_mode)
    D = pairwise_distances(X, periodic=periodic)

    W = np.exp(-np.power(D / max(lam, 1e-12), power))
    np.fill_diagonal(W, 0.0)

    if background > 0:
        B = rng.random((N, N))
        B = 0.5 * (B + B.T)
        np.fill_diagonal(B, 0.0)
        W = (1.0 - background) * W + background * B

    if noise > 0:
        E = rng.normal(scale=noise, size=(N, N))
        E = 0.5 * (E + E.T)
        np.fill_diagonal(E, 0.0)
        W = W + E

    W = clean_matrix(W)

    # Normalize to max 1 for stable comparison.
    wmax = np.max(W)
    if wmax > 0:
        W = W / wmax

    return W, X, D


# ============================================================
# kNN graph and Laplacian
# ============================================================

def knn_mask(W, k):
    N = W.shape[0]
    k = min(k, N - 1)
    mask = np.zeros_like(W, dtype=bool)

    for i in range(N):
        idx = np.argsort(W[i])[::-1]
        idx = [j for j in idx if j != i and W[i, j] > 0]
        for j in idx[:k]:
            mask[i, j] = True

    mask = np.logical_or(mask, mask.T)
    np.fill_diagonal(mask, False)
    return mask


def adjacency_from_mask(W, mask):
    A = np.zeros_like(W)
    A[mask] = W[mask]
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def sparse_laplacian_from_adjacency(A, normalized=False, eps=1e-12):
    A_sp = sparse.csr_matrix(A)
    deg = np.asarray(A_sp.sum(axis=1)).ravel()

    if normalized:
        inv_sqrt = np.zeros_like(deg)
        good = deg > eps
        inv_sqrt[good] = 1.0 / np.sqrt(deg[good])

        D_inv = sparse.diags(inv_sqrt)
        L = sparse.eye(A.shape[0], format="csr") - D_inv @ A_sp @ D_inv
        return L.tocsr(), deg

    L = sparse.diags(deg) - A_sp
    return L.tocsr(), deg


# ============================================================
# Spectral dimension from eigenvalues
# ============================================================

def spectral_dimension_from_eigs(
    L,
    tau_min=0.01,
    tau_max=100.0,
    tau_points=80,
    eig_k=200,
    dense_threshold=500
):
    """
    Estimate heat trace:

        Z(tau)=Tr exp(-tau L)

    using all eigenvalues if N <= dense_threshold,
    else the smallest eig_k eigenvalues.

    For large graphs, this is an approximation biased toward intermediate/large tau.
    """
    N = L.shape[0]

    if N <= dense_threshold:
        evals = np.linalg.eigvalsh(L.toarray())
    else:
        k = min(eig_k, N - 2)
        evals = eigsh(L, k=k, which="SM", return_eigenvectors=False)
        evals = np.sort(np.real(evals))

    evals = np.maximum(evals, 0.0)

    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_points)

    Z = []
    for t in taus:
        Z.append(float(np.sum(np.exp(-t * evals))))

    Z = np.asarray(Z)
    logt = np.log(taus)
    logZ = np.log(np.maximum(Z, 1e-300))

    slope = np.gradient(logZ, logt)
    ds = -2.0 * slope

    # Use central tau window.
    n = len(ds)
    lo = max(1, n // 4)
    hi = min(n - 1, 3 * n // 4)
    ds_eff = float(np.nanmedian(ds[lo:hi]))

    return {
        "taus": taus,
        "Z": Z,
        "ds_curve": ds,
        "ds_eff": ds_eff,
        "evals_used": int(len(evals)),
        "lambda_0": float(np.min(evals)),
        "lambda_1": float(np.sort(evals)[1]) if len(evals) > 1 else np.nan,
        "lambda_max_used": float(np.max(evals)),
    }


# ============================================================
# One scan point
# ============================================================

def run_one(args, N, dim, lam, k, seed):
    W, X, D = generate_mi_matrix(
        N=N,
        dim=dim,
        lam=lam,
        power=args.power,
        seed=seed,
        noise=args.noise,
        background=args.background,
        periodic=args.periodic,
        point_mode=args.point_mode
    )

    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)

    L, deg = sparse_laplacian_from_adjacency(
        A,
        normalized=args.normalized_laplacian
    )

    spec = spectral_dimension_from_eigs(
        L,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        tau_points=args.tau_points,
        eig_k=args.eig_k,
        dense_threshold=args.dense_threshold
    )

    edges = int(mask.sum() // 2)

    row = {
        "N": N,
        "dim_ambient": dim,
        "lambda": lam,
        "k": k,
        "seed": seed,
        "power": args.power,
        "noise": args.noise,
        "background": args.background,
        "periodic": bool(args.periodic),
        "point_mode": args.point_mode,
        "normalized_laplacian": bool(args.normalized_laplacian),

        "d_s_eff": spec["ds_eff"],
        "target_ds": args.target_ds,
        "abs_error_to_target": abs(spec["ds_eff"] - args.target_ds),

        "edges": edges,
        "mean_degree": float(np.mean(deg)),
        "min_degree": float(np.min(deg)),
        "max_degree": float(np.max(deg)),

        "evals_used": spec["evals_used"],
        "lambda_0": spec["lambda_0"],
        "lambda_1": spec["lambda_1"],
        "lambda_max_used": spec["lambda_max_used"],
    }

    matrix_file = None

    if args.save_all_mi or row["abs_error_to_target"] <= args.save_if_error_below:
        mi_dir = os.path.join(args.output_dir, "mi_matrices")
        ensure_dir(mi_dir)

        lam_str = str(lam).replace(".", "p")
        fname = f"MI_N{N}_d{dim}_lam{lam_str}_k{k}_seed{seed}.csv"
        matrix_file = os.path.join(mi_dir, fname)
        np.savetxt(matrix_file, W, delimiter=",")

    row["mi_file"] = matrix_file if matrix_file is not None else ""

    return row


# ============================================================
# Figures
# ============================================================

def make_figures(df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    # d_s vs lambda by k and dimension
    for dim in sorted(df["dim_ambient"].unique()):
        subd = df[df["dim_ambient"] == dim]

        plt.figure(figsize=(8, 5))
        for k in sorted(subd["k"].unique()):
            sub = subd[subd["k"] == k]
            grouped = sub.groupby("lambda")["d_s_eff"].mean().reset_index()
            plt.plot(grouped["lambda"], grouped["d_s_eff"], "o-", label=f"k={k}")

        plt.axhline(3.0, linestyle="--", linewidth=1, label=r"$d_s=3$")
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$d_s^{\rm eff}$")
        plt.title(f"Spectral dimension vs lambda — ambient d={dim}")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"fig_ds_vs_lambda_dim{dim}.png"), dpi=250)
        plt.close()

    # Best error to target
    best = df.sort_values("abs_error_to_target").head(30)

    plt.figure(figsize=(10, 5))
    labels = [
        f"N{int(r['N'])} d{int(r['dim_ambient'])} k{int(r['k'])} lam{float(r['lambda']):g}"
        for _, r in best.iterrows()    
    ]
    plt.bar(range(len(best)), best["abs_error_to_target"])
    plt.xticks(range(len(best)), labels, rotation=90, fontsize=7)
    plt.ylabel(r"$|d_s-3|$")
    plt.title("Best candidates near target spectral dimension")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_best_ds_target.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def parse_float_list(values):
    return [float(x) for x in values]


def parse_int_list(values):
    return [int(x) for x in values]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--output-dir", required=True)

    parser.add_argument("--N-list", nargs="+", type=int, default=[100, 200, 300])
    parser.add_argument("--dim-list", nargs="+", type=int, default=[2, 3])
    parser.add_argument("--lambda-list", nargs="+", type=float, default=[0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5])
    parser.add_argument("--k-list", nargs="+", type=int, default=[6, 8, 10, 12, 16])
    parser.add_argument("--seed-list", nargs="+", type=int, default=[0])

    parser.add_argument("--target-ds", type=float, default=3.0)

    parser.add_argument("--power", type=float, default=1.0)
    parser.add_argument("--noise", type=float, default=0.0)
    parser.add_argument("--background", type=float, default=0.0)
    parser.add_argument("--periodic", action="store_true")
    parser.add_argument("--point-mode", choices=["cube", "gaussian", "sphere"], default="cube")

    parser.add_argument("--normalized-laplacian", action="store_true")

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=100.0)
    parser.add_argument("--tau-points", type=int, default=80)

    parser.add_argument("--eig-k", type=int, default=200)
    parser.add_argument("--dense-threshold", type=int, default=500)

    parser.add_argument("--save-all-mi", action="store_true")
    parser.add_argument("--save-if-error-below", type=float, default=0.15)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    rows = []
    total = (
        len(args.N_list)
        * len(args.dim_list)
        * len(args.lambda_list)
        * len(args.k_list)
        * len(args.seed_list)
    )

    counter = 0

    print("=" * 100)
    print("BuP Paper 10 — Generate MI-like matrices and scan for d_s ≈ 3")
    print("=" * 100)
    print(f"output: {args.output_dir}")
    print(f"total runs: {total}")
    print("-" * 100)

    for N in args.N_list:
        for dim in args.dim_list:
            for lam in args.lambda_list:
                for k in args.k_list:
                    for seed in args.seed_list:
                        counter += 1
                        print(
                            f"[{counter}/{total}] "
                            f"N={N} dim={dim} lambda={lam:g} k={k} seed={seed}"
                        )

                        try:
                            row = run_one(args, N, dim, lam, k, seed)
                            rows.append(row)
                            print(
                                f"    d_s={row['d_s_eff']:.4f} "
                                f"|err|={row['abs_error_to_target']:.4f} "
                                f"edges={row['edges']}"
                            )
                        except Exception as e:
                            print(f"    ERROR: {e}")
                            rows.append({
                                "N": N,
                                "dim_ambient": dim,
                                "lambda": lam,
                                "k": k,
                                "seed": seed,
                                "error": str(e)
                            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(args.output_dir, "ds_scan_summary.csv")
    df.to_csv(csv_path, index=False)

    ok_df = df[df.get("d_s_eff").notna()] if "d_s_eff" in df.columns else pd.DataFrame()
    if len(ok_df) > 0:
        make_figures(ok_df, args.output_dir)

    best = []
    if len(ok_df) > 0:
        best = ok_df.sort_values("abs_error_to_target").head(20).to_dict(orient="records")

    summary = {
        "title": "BuP Paper 10 — MI generator scan for target spectral dimension",
        "description": "Generates MI-like matrices and scans (N, dim, lambda, k) for d_s close to target.",
        "target_ds": args.target_ds,
        "parameters": vars(args),
        "best_candidates": best,
        "files": {
            "scan_csv": "ds_scan_summary.csv",
            "mi_matrices_dir": "mi_matrices/",
            "figures_dir": "figures/"
        }
    }

    json_path = os.path.join(args.output_dir, "summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 100)
    print("BEST CANDIDATES")
    print("=" * 100)
    if len(ok_df) > 0:
        print(ok_df.sort_values("abs_error_to_target").head(20)[[
            "N",
            "dim_ambient",
            "lambda",
            "k",
            "seed",
            "d_s_eff",
            "abs_error_to_target",
            "mi_file"
        ]].to_string(index=False))
    print("=" * 100)
    print("Files:")
    print(csv_path)
    print(json_path)
    print("DONE")


if __name__ == "__main__":
    main()

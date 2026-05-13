#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_grad_source_v1.py

Paper 8 — Scan du terme de contrainte interne chi.

On teste une source matière complète :

    S(i) = T00(i) + omega*Taa(i) + chi*Tgrad(i)

où :
    T00(i)   = 1/2 sum_j (deltaW_ij)^2
    Taa(i)   = 1/2 sum_j (deltaW_ij)^2 ||v_ij||^2
    Tgrad(i) = sum_{j voisin} W_ij (T00_i - T00_j)^2

Objectif :
    Tester si l'ajout d'un terme de gradient interne améliore la corrélation
    avec la réponse de courbure nodale :

        |delta R_i|
        delta R_i

Usage :
    python3 bup_paper8_Tmatter_grad_source_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.01 0.02 0.03 0.05 0.08 0.10 0.15 \
      --omega -0.5 \
      --chi-list -10 -5 -3 -1 -0.5 0 0.5 1 3 5 10 \
      --output-dir results_paper8_Tmatter_grad_source_N20_lam057
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.sparse.csgraph import shortest_path
from scipy.stats import pearsonr, spearmanr
from sklearn.manifold import MDS


# ============================================================
# Utils
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


def load_matrix_csv(path):
    try:
        W = np.loadtxt(path, delimiter=",")
        if W.ndim == 2 and W.shape[0] == W.shape[1]:
            return clean_matrix(W)
    except Exception:
        pass

    df = pd.read_csv(path, header=None)
    W = df.values

    if not np.issubdtype(W.dtype, np.number):
        df = pd.read_csv(path, index_col=0)
        W = df.values

    return clean_matrix(W.astype(float))


def parse_label(path):
    return os.path.splitext(os.path.basename(path))[0]


def parse_lambda(path):
    name = os.path.basename(path)
    m = re.search(r"lam([0-9]+(?:\.[0-9]+)?)", name)
    if m:
        return float(m.group(1))
    return None


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    if len(x) < 4:
        return np.nan, np.nan, np.nan, np.nan

    if np.std(x) < 1e-14 or np.std(y) < 1e-14:
        return np.nan, np.nan, np.nan, np.nan

    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    return float(pr), float(pp), float(sr), float(sp)


# ============================================================
# Graph
# ============================================================

def knn_mask(W, k):
    n = W.shape[0]
    k = min(k, n - 1)
    mask = np.zeros_like(W, dtype=bool)

    for i in range(n):
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


def entanglement_lengths(W, mask, eps=1e-12):
    lengths = np.full_like(W, np.inf, dtype=float)
    vals = W[mask]

    if vals.size == 0 or np.max(vals) <= 0:
        raise ValueError("Aucune arête positive dans le graphe kNN.")

    Wmax = np.max(vals)
    lengths[mask] = -np.log((W[mask] + eps) / (Wmax + eps))
    lengths[mask] = np.maximum(lengths[mask], 1e-9)
    np.fill_diagonal(lengths, 0.0)
    return lengths


def entanglement_distance(W, mask, eps=1e-12):
    lengths = entanglement_lengths(W, mask, eps=eps)
    D = shortest_path(lengths, method="FW", directed=False)
    D[np.isinf(D)] = 10.0 * W.shape[0]
    return D


def choose_center(W):
    return int(np.argmax(np.sum(W, axis=1)))


# ============================================================
# MDS coordinates
# ============================================================

def mds_coordinates(D, ndim=3, random_state=0):
    D = np.asarray(D, dtype=float)
    finite = np.isfinite(D)
    max_finite = np.max(D[finite]) if np.any(finite) else 1.0
    D = np.nan_to_num(D, nan=max_finite, posinf=max_finite)

    model = MDS(
        n_components=ndim,
        dissimilarity="precomputed",
        random_state=random_state,
        normalized_stress="auto",
        n_init=4
    )
    return model.fit_transform(D)


# ============================================================
# Ollivier-Ricci proxy
# ============================================================

def neighborhood_measure(A, i):
    idx = np.where(A[i] > 0)[0]

    if idx.size == 0:
        return np.array([i]), np.array([1.0])

    w = A[i, idx].astype(float)
    s = np.sum(w)

    if s <= 1e-15:
        return np.array([i]), np.array([1.0])

    return idx, w / s


def wasserstein_greedy(idx_a, mass_a, idx_b, mass_b, D):
    supply = mass_a.copy()
    demand = mass_b.copy()
    cost = 0.0

    pairs = []
    for ia, a in enumerate(idx_a):
        for ib, b in enumerate(idx_b):
            pairs.append((D[a, b], ia, ib))

    pairs.sort(key=lambda x: x[0])

    for c, ia, ib in pairs:
        amount = min(supply[ia], demand[ib])

        if amount > 0:
            cost += amount * c
            supply[ia] -= amount
            demand[ib] -= amount

        if np.all(supply <= 1e-14) and np.all(demand <= 1e-14):
            break

    return cost


def ollivier_kappa(A, D, mask):
    n = A.shape[0]
    measures = [neighborhood_measure(A, i) for i in range(n)]
    kappas = {}

    edges = np.argwhere(np.triu(mask, 1))

    for i, j in edges:
        dij = D[i, j]
        if not np.isfinite(dij) or dij <= 1e-12:
            continue

        idx_i, mass_i = measures[i]
        idx_j, mass_j = measures[j]

        w1 = wasserstein_greedy(idx_i, mass_i, idx_j, mass_j, D)
        kappas[(int(i), int(j))] = float(1.0 - w1 / dij)

    return kappas


def node_curvature_R(W, k, eps=1e-12):
    W = clean_matrix(W)
    n = W.shape[0]

    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    D = entanglement_distance(W, mask, eps=eps)

    kappas = ollivier_kappa(A, D, mask)

    node_vals = [[] for _ in range(n)]
    for (i, j), kij in kappas.items():
        node_vals[i].append(kij)
        node_vals[j].append(kij)

    R = np.zeros(n)
    for i in range(n):
        R[i] = np.mean(node_vals[i]) if node_vals[i] else np.nan

    return R, A, D, mask


# ============================================================
# Excitation
# ============================================================

def inject_excitation(W0, D0, amp=0.15, sigma=0.05, center=-1):
    W0 = clean_matrix(W0)
    n = W0.shape[0]

    if center is None or center < 0:
        center = choose_center(W0)

    d = D0[center].copy()
    finite = np.isfinite(d)

    if not np.any(finite):
        raise ValueError("Distances non finies.")

    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    f = f / max(np.max(f), 1e-15)

    vals = W0[W0 > 0]
    scale = float(np.mean(vals)) if vals.size else 1.0

    deltaW = amp * scale * np.outer(f, f)
    np.fill_diagonal(deltaW, 0.0)

    W_exc = clean_matrix(W0 + deltaW)
    deltaW_actual = W_exc - W0
    np.fill_diagonal(deltaW_actual, 0.0)

    return W_exc, deltaW_actual, center, f


# ============================================================
# Tensor components
# ============================================================

def compute_T00(deltaW):
    return 0.5 * np.sum(deltaW ** 2, axis=1)


def compute_Taa(deltaW, D, X, eps=1e-12):
    n = deltaW.shape[0]
    Taa = np.zeros(n, dtype=float)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            dx = X[j] - X[i]
            v2 = float(np.dot(dx, dx) / ((dij + eps) ** 2))
            Taa[i] += 0.5 * (deltaW[i, j] ** 2) * v2

    return Taa


def compute_Tgrad(T00, A):
    """
    Terme de contrainte interne nodal :

        Tgrad(i) = sum_j A_ij (T00_i - T00_j)^2

    où A_ij est l'adjacence pondérée du graphe kNN de référence.

    Interprétation :
        mesure la tension / variation interne de densité d'énergie autour du noeud.
    """
    n = len(T00)
    Tgrad = np.zeros(n, dtype=float)

    for i in range(n):
        for j in range(n):
            wij = A[i, j]
            if wij <= 0:
                continue
            Tgrad[i] += wij * (T00[i] - T00[j]) ** 2

    return Tgrad


def norm01(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-15:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


# ============================================================
# One configuration
# ============================================================

def prepare_components(W0, args, sigma):
    R0, A0, D0, mask0 = node_curvature_R(W0, args.k, eps=args.eps)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    R1, A1, D1, mask1 = node_curvature_R(W_exc, args.k, eps=args.eps)

    X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)

    T00 = compute_T00(deltaW)
    Taa = compute_Taa(deltaW, D0, X, eps=args.eps)
    Tgrad = compute_Tgrad(T00, A0)

    delta_R = R1 - R0
    abs_delta_R = np.abs(delta_R)

    return {
        "T00": T00,
        "Taa": Taa,
        "Tgrad": Tgrad,
        "delta_R": delta_R,
        "abs_delta_R": abs_delta_R,
        "center": center,
        "D0": D0,
        "A0": A0,
        "R0": R0,
        "R1": R1,
        "deltaW": deltaW
    }


def run_one(W0, args, sigma, chi):
    comp = prepare_components(W0, args, sigma)

    T00 = comp["T00"]
    Taa = comp["Taa"]
    Tgrad = comp["Tgrad"]
    delta_R = comp["delta_R"]
    abs_delta_R = comp["abs_delta_R"]

    # Option : normaliser les composantes pour que chi soit comparable
    if args.normalize_components:
        T00_eff = norm01(T00)
        Taa_eff = norm01(Taa)
        Tgrad_eff = norm01(Tgrad)
    else:
        T00_eff = T00
        Taa_eff = Taa
        Tgrad_eff = Tgrad

    S_fluid = T00_eff + args.omega * Taa_eff
    S_full = T00_eff + args.omega * Taa_eff + chi * Tgrad_eff

    p_full_abs, pp_full_abs, s_full_abs, sp_full_abs = safe_corr(S_full, abs_delta_R)
    p_full_signed, pp_full_signed, s_full_signed, sp_full_signed = safe_corr(S_full, delta_R)

    p_fluid_abs, pp_fluid_abs, s_fluid_abs, sp_fluid_abs = safe_corr(S_fluid, abs_delta_R)
    p_fluid_signed, pp_fluid_signed, s_fluid_signed, sp_fluid_signed = safe_corr(S_fluid, delta_R)

    p_grad_abs, pp_grad_abs, s_grad_abs, sp_grad_abs = safe_corr(Tgrad_eff, abs_delta_R)
    p_grad_signed, pp_grad_signed, s_grad_signed, sp_grad_signed = safe_corr(Tgrad_eff, delta_R)

    return {
        "sigma": sigma,
        "chi": chi,
        "omega": args.omega,
        "normalize_components": args.normalize_components,

        "T00_sum": float(np.sum(T00)),
        "Taa_sum": float(np.sum(Taa)),
        "Tgrad_sum": float(np.sum(Tgrad)),

        "Taa_over_T00": float(np.sum(Taa) / max(np.sum(T00), 1e-15)),
        "Tgrad_over_T00": float(np.sum(Tgrad) / max(np.sum(T00), 1e-15)),

        "center": int(comp["center"]),

        "pearson_full_abs_delta_R": p_full_abs,
        "pearson_full_abs_delta_R_p": pp_full_abs,
        "spearman_full_abs_delta_R": s_full_abs,
        "spearman_full_abs_delta_R_p": sp_full_abs,

        "pearson_full_delta_R": p_full_signed,
        "pearson_full_delta_R_p": pp_full_signed,
        "spearman_full_delta_R": s_full_signed,
        "spearman_full_delta_R_p": sp_full_signed,

        "spearman_fluid_abs_delta_R": s_fluid_abs,
        "spearman_fluid_abs_delta_R_p": sp_fluid_abs,
        "spearman_fluid_delta_R": s_fluid_signed,
        "spearman_fluid_delta_R_p": sp_fluid_signed,

        "spearman_grad_abs_delta_R": s_grad_abs,
        "spearman_grad_abs_delta_R_p": sp_grad_abs,
        "spearman_grad_delta_R": s_grad_signed,
        "spearman_grad_delta_R_p": sp_grad_signed,
    }


# ============================================================
# Plots
# ============================================================

def make_plots(df, output_dir):
    plt.figure(figsize=(9, 5))
    for sigma in sorted(df["sigma"].unique()):
        sub = df[df["sigma"] == sigma].sort_values("chi")
        plt.plot(sub["chi"], sub["spearman_full_abs_delta_R"], "o-", label=f"sigma={sigma:g}")

    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\chi$")
    plt.ylabel(r"Spearman$(S_{\rm full},|\delta R|)$")
    plt.title(r"Scan $\chi$ : réponse absolue")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_chi_abs_response.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(9, 5))
    for sigma in sorted(df["sigma"].unique()):
        sub = df[df["sigma"] == sigma].sort_values("chi")
        plt.plot(sub["chi"], sub["spearman_full_delta_R"], "o-", label=f"sigma={sigma:g}")

    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\chi$")
    plt.ylabel(r"Spearman$(S_{\rm full},\delta R)$")
    plt.title(r"Scan $\chi$ : réponse signée")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_chi_signed_response.png"), dpi=250)
    plt.close()

    best_abs = df.sort_values("spearman_full_abs_delta_R", ascending=False).head(20)
    plt.figure(figsize=(10, 5))
    labels = [f"s={r.sigma:g},c={r.chi:g}" for _, r in best_abs.iterrows()]
    vals = best_abs["spearman_full_abs_delta_R"].values
    plt.bar(range(len(vals)), vals)
    plt.xticks(range(len(vals)), labels, rotation=60, ha="right", fontsize=8)
    plt.ylabel(r"Spearman$(S_{\rm full},|\delta R|)$")
    plt.title("Top 20 — source complète")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_top20_full_abs_response.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_grad_source_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)

    parser.add_argument("--omega", type=float, default=-0.5)
    parser.add_argument("--chi-list", type=float, nargs="+", required=True)

    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    parser.add_argument("--ndim", type=int, default=3)
    parser.add_argument("--random-state", type=int, default=0)

    parser.add_argument("--normalize-components", action="store_true",
                        help="Normalise T00, Taa, Tgrad in [0,1] before combining.")

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 100)
    print("BuP Paper 8 — Scan chi pour terme de contrainte interne")
    print("=" * 100)
    print(f"MI file     : {args.mi_file}")
    print(f"label       : {label}")
    print(f"lambda      : {lam}")
    print(f"N           : {W0.shape[0]}")
    print(f"k           : {args.k}")
    print(f"amp         : {args.amp}")
    print(f"omega       : {args.omega}")
    print(f"normalize   : {args.normalize_components}")
    print(f"sigma list  : {args.sigma_list}")
    print(f"chi list    : {args.chi_list}")
    print(f"output      : {args.output_dir}")
    print("-" * 100)

    rows = []

    for sigma in args.sigma_list:
        for chi in args.chi_list:
            row = run_one(W0, args, sigma=sigma, chi=chi)
            rows.append(row)

            print(
                f"sigma={sigma:.3g} chi={chi:+.3g} "
                f"Sp_abs={row['spearman_full_abs_delta_R']:+.4f} "
                f"p_abs={row['spearman_full_abs_delta_R_p']:.3g} "
                f"Sp_signed={row['spearman_full_delta_R']:+.4f} "
                f"p_signed={row['spearman_full_delta_R_p']:.3g}"
            )

    df = pd.DataFrame(rows)
    csv_path = os.path.join(args.output_dir, "grad_source_chi_scan.csv")
    df.to_csv(csv_path, index=False)

    best_abs = df.sort_values("spearman_full_abs_delta_R", ascending=False).iloc[0].to_dict()
    best_signed = df.sort_values("spearman_full_delta_R", ascending=False).iloc[0].to_dict()

    best_abs_p = df.sort_values(
        ["spearman_full_abs_delta_R_p", "spearman_full_abs_delta_R"],
        ascending=[True, False]
    ).iloc[0].to_dict()

    best_signed_p = df.sort_values(
        ["spearman_full_delta_R_p", "spearman_full_delta_R"],
        ascending=[True, False]
    ).iloc[0].to_dict()

    summary = {
        "title": "BuP Paper 8 — Gradient/internal stress source scan",
        "description": "Tests S_full = T00 + omega*Taa + chi*Tgrad against curvature response.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "omega": args.omega,
        "normalize_components": args.normalize_components,
        "sigma_list": args.sigma_list,
        "chi_list": args.chi_list,
        "definition": {
            "T00": "1/2 sum_j deltaW_ij^2",
            "Taa": "1/2 sum_j deltaW_ij^2 ||v_ij||^2",
            "Tgrad": "sum_j A_ij (T00_i - T00_j)^2",
            "S_full": "T00 + omega*Taa + chi*Tgrad"
        },
        "best_abs_response_by_spearman": best_abs,
        "best_signed_response_by_spearman": best_signed,
        "best_abs_response_by_pvalue": best_abs_p,
        "best_signed_response_by_pvalue": best_signed_p,
        "files": {
            "scan_csv": "grad_source_chi_scan.csv",
            "fig_abs": "fig_chi_abs_response.png",
            "fig_signed": "fig_chi_signed_response.png",
            "fig_top20": "fig_top20_full_abs_response.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    make_plots(df, args.output_dir)

    print("-" * 100)
    print("MEILLEURE RÉPONSE ABSOLUE |delta R|")
    print("-" * 100)
    print(f"sigma     : {best_abs['sigma']}")
    print(f"chi       : {best_abs['chi']}")
    print(f"Spearman  : {best_abs['spearman_full_abs_delta_R']:+.6f}")
    print(f"p         : {best_abs['spearman_full_abs_delta_R_p']:.6g}")
    print(f"Pearson   : {best_abs['pearson_full_abs_delta_R']:+.6f}")
    print(f"p         : {best_abs['pearson_full_abs_delta_R_p']:.6g}")

    print("-" * 100)
    print("MEILLEURE RÉPONSE SIGNÉE delta R")
    print("-" * 100)
    print(f"sigma     : {best_signed['sigma']}")
    print(f"chi       : {best_signed['chi']}")
    print(f"Spearman  : {best_signed['spearman_full_delta_R']:+.6f}")
    print(f"p         : {best_signed['spearman_full_delta_R_p']:.6g}")
    print(f"Pearson   : {best_signed['pearson_full_delta_R']:+.6f}")
    print(f"p         : {best_signed['pearson_full_delta_R_p']:.6g}")

    print("-" * 100)
    print("Fichiers générés")
    print("-" * 100)
    print(csv_path)
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_chi_abs_response.png"))
    print(os.path.join(args.output_dir, "fig_chi_signed_response.png"))
    print("DONE")


if __name__ == "__main__":
    main()

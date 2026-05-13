#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_fluid_source_v1.py

Paper 8 — Test fluide effectif :
    T00 seul est positif/localisé mais insuffisant pour prédire la courbure.
    On teste donc une source de type fluide :

        S_fluid(i) = T00(i) + omega * Taa(i)

    avec :
        T00(i) = 1/2 sum_j (deltaW_ij)^2
        Taa(i) = 1/2 sum_j (deltaW_ij)^2 ||v_ij||^2

    On compare S_fluid à :
        |delta R_i|
        delta R_i signé

Usage :
    python3 bup_paper8_Tmatter_fluid_source_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.01 0.02 0.03 0.05 0.08 0.10 0.15 \
      --omega-list -3 -1 -0.5 0 0.5 1 3 \
      --output-dir results_paper8_Tmatter_fluid_source_N20_lam057
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
    """
    Coordonnées émergentes à partir de la distance d'intrication.
    Utilisées seulement pour construire ||v_ij||^2.
    """
    D = np.asarray(D, dtype=float)
    D = np.nan_to_num(D, nan=np.nanmax(D[np.isfinite(D)]), posinf=np.nanmax(D[np.isfinite(D)]))

    model = MDS(
        n_components=ndim,
        dissimilarity="precomputed",
        random_state=random_state,
        normalized_stress="auto"
    )
    X = model.fit_transform(D)
    return X


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
    edge_rows = []

    for (i, j), kij in kappas.items():
        node_vals[i].append(kij)
        node_vals[j].append(kij)
        edge_rows.append({
            "i": i,
            "j": j,
            "kappa": kij,
            "W": W[i, j],
            "d_ent": D[i, j]
        })

    R = np.zeros(n)
    for i in range(n):
        R[i] = np.mean(node_vals[i]) if node_vals[i] else np.nan

    return R, pd.DataFrame(edge_rows), D, mask


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
    """
    Taa(i) = 1/2 sum_j (deltaW_ij)^2 ||v_ij||^2

    v_ij^a = (x_j^a - x_i^a) / (d_ent_ij + eps)

    Dans un embedding MDS cohérent, ||v||^2 mesure la part spatiale isotrope.
    """
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


def localization_radius(source, D0, center):
    total = np.sum(source)
    if total <= 1e-15:
        return np.nan

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    return float(np.sqrt(np.sum(source * d ** 2) / total))


def participation_ratio(source):
    s1 = np.sum(source)
    s2 = np.sum(source ** 2)
    if s2 <= 1e-15:
        return np.nan
    return float((s1 ** 2) / s2)


# ============================================================
# Run one
# ============================================================

def run_one(W0, args, sigma, omega, X_mode="mds"):
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)
    n = W0.shape[0]

    R0, edge0, D0, mask0 = node_curvature_R(W0, args.k, eps=args.eps)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    R1, edge1, D1, mask1 = node_curvature_R(W_exc, args.k, eps=args.eps)

    if X_mode == "mds":
        X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)
    else:
        raise ValueError("Only X_mode='mds' is implemented.")

    T00 = compute_T00(deltaW)
    Taa = compute_Taa(deltaW, D0, X, eps=args.eps)

    S_fluid = T00 + omega * Taa
    S_trace_like = T00 + Taa

    delta_R = R1 - R0
    abs_delta_R = np.abs(delta_R)

    # corrélations principales
    p_T00_abs, pp_T00_abs, s_T00_abs, sp_T00_abs = safe_corr(T00, abs_delta_R)
    p_Taa_abs, pp_Taa_abs, s_Taa_abs, sp_Taa_abs = safe_corr(Taa, abs_delta_R)
    p_fluid_abs, pp_fluid_abs, s_fluid_abs, sp_fluid_abs = safe_corr(S_fluid, abs_delta_R)
    p_trace_abs, pp_trace_abs, s_trace_abs, sp_trace_abs = safe_corr(S_trace_like, abs_delta_R)

    p_T00_signed, pp_T00_signed, s_T00_signed, sp_T00_signed = safe_corr(T00, delta_R)
    p_Taa_signed, pp_Taa_signed, s_Taa_signed, sp_Taa_signed = safe_corr(Taa, delta_R)
    p_fluid_signed, pp_fluid_signed, s_fluid_signed, sp_fluid_signed = safe_corr(S_fluid, delta_R)
    p_trace_signed, pp_trace_signed, s_trace_signed, sp_trace_signed = safe_corr(S_trace_like, delta_R)

    out = {
        "label": label,
        "lambda": lam,
        "N": n,
        "k": args.k,
        "amp": args.amp,
        "sigma": sigma,
        "omega": omega,
        "center": center,

        "T00_sum": float(np.sum(T00)),
        "Taa_sum": float(np.sum(Taa)),
        "Taa_over_T00_sum": float(np.sum(Taa) / max(np.sum(T00), 1e-15)),
        "S_fluid_sum": float(np.sum(S_fluid)),
        "S_trace_like_sum": float(np.sum(S_trace_like)),

        "T00_Rloc": localization_radius(T00, D0, center),
        "Taa_Rloc": localization_radius(np.abs(Taa), D0, center),
        "S_fluid_Rloc": localization_radius(np.abs(S_fluid), D0, center),

        "T00_PR": participation_ratio(T00),
        "Taa_PR": participation_ratio(np.abs(Taa)),
        "S_fluid_PR": participation_ratio(np.abs(S_fluid)),

        "delta_R_min": float(np.nanmin(delta_R)),
        "delta_R_max": float(np.nanmax(delta_R)),
        "abs_delta_R_mean": float(np.nanmean(abs_delta_R)),

        # abs curvature response
        "pearson_T00_abs_delta_R": p_T00_abs,
        "pearson_T00_abs_delta_R_p": pp_T00_abs,
        "spearman_T00_abs_delta_R": s_T00_abs,
        "spearman_T00_abs_delta_R_p": sp_T00_abs,

        "pearson_Taa_abs_delta_R": p_Taa_abs,
        "pearson_Taa_abs_delta_R_p": pp_Taa_abs,
        "spearman_Taa_abs_delta_R": s_Taa_abs,
        "spearman_Taa_abs_delta_R_p": sp_Taa_abs,

        "pearson_fluid_abs_delta_R": p_fluid_abs,
        "pearson_fluid_abs_delta_R_p": pp_fluid_abs,
        "spearman_fluid_abs_delta_R": s_fluid_abs,
        "spearman_fluid_abs_delta_R_p": sp_fluid_abs,

        "pearson_trace_abs_delta_R": p_trace_abs,
        "pearson_trace_abs_delta_R_p": pp_trace_abs,
        "spearman_trace_abs_delta_R": s_trace_abs,
        "spearman_trace_abs_delta_R_p": sp_trace_abs,

        # signed curvature response
        "pearson_T00_delta_R": p_T00_signed,
        "pearson_T00_delta_R_p": pp_T00_signed,
        "spearman_T00_delta_R": s_T00_signed,
        "spearman_T00_delta_R_p": sp_T00_signed,

        "pearson_Taa_delta_R": p_Taa_signed,
        "pearson_Taa_delta_R_p": pp_Taa_signed,
        "spearman_Taa_delta_R": s_Taa_signed,
        "spearman_Taa_delta_R_p": sp_Taa_signed,

        "pearson_fluid_delta_R": p_fluid_signed,
        "pearson_fluid_delta_R_p": pp_fluid_signed,
        "spearman_fluid_delta_R": s_fluid_signed,
        "spearman_fluid_delta_R_p": sp_fluid_signed,

        "pearson_trace_delta_R": p_trace_signed,
        "pearson_trace_delta_R_p": pp_trace_signed,
        "spearman_trace_delta_R": s_trace_signed,
        "spearman_trace_delta_R_p": sp_trace_signed,
    }

    return out


# ============================================================
# Plots
# ============================================================

def make_plots(df, output_dir):
    # abs response
    plt.figure(figsize=(9, 5))
    for omega in sorted(df["omega"].unique()):
        sub = df[df["omega"] == omega].sort_values("sigma")
        plt.plot(sub["sigma"], sub["spearman_fluid_abs_delta_R"], "o-", label=f"omega={omega:g}")

    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"Spearman$(T_{00}+\omega T_{aa},|\delta R|)$")
    plt.title("Source fluide vs réponse de courbure absolue")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_fluid_abs_response_vs_sigma.png"), dpi=250)
    plt.close()

    # signed response
    plt.figure(figsize=(9, 5))
    for omega in sorted(df["omega"].unique()):
        sub = df[df["omega"] == omega].sort_values("sigma")
        plt.plot(sub["sigma"], sub["spearman_fluid_delta_R"], "o-", label=f"omega={omega:g}")

    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"Spearman$(T_{00}+\omega T_{aa},\delta R)$")
    plt.title("Source fluide vs réponse de courbure signée")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_fluid_signed_response_vs_sigma.png"), dpi=250)
    plt.close()

    # best by omega
    best_abs = df.sort_values("spearman_fluid_abs_delta_R", ascending=False).head(20)

    plt.figure(figsize=(10, 5))
    labels = [f"s={r.sigma:g},w={r.omega:g}" for _, r in best_abs.iterrows()]
    vals = best_abs["spearman_fluid_abs_delta_R"].values
    plt.bar(range(len(vals)), vals)
    plt.xticks(range(len(vals)), labels, rotation=60, ha="right", fontsize=8)
    plt.ylabel(r"Spearman$(S_{\rm fluid},|\delta R|)$")
    plt.title("Top 20 sources fluides — réponse absolue")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_top20_fluid_abs_response.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_fluid_source_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)

    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)
    parser.add_argument("--omega-list", type=float, nargs="+", default=[-3, -1, -0.5, 0, 0.5, 1, 3])

    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    parser.add_argument("--ndim", type=int, default=3)
    parser.add_argument("--random-state", type=int, default=0)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 100)
    print("BuP Paper 8 — Source fluide T00 + omega Taa")
    print("=" * 100)
    print(f"MI file    : {args.mi_file}")
    print(f"label      : {label}")
    print(f"lambda     : {lam}")
    print(f"N          : {W0.shape[0]}")
    print(f"k          : {args.k}")
    print(f"amp        : {args.amp}")
    print(f"sigma list : {args.sigma_list}")
    print(f"omega list : {args.omega_list}")
    print(f"output     : {args.output_dir}")
    print("-" * 100)

    rows = []

    for sigma in args.sigma_list:
        for omega in args.omega_list:
            row = run_one(W0, args, sigma=sigma, omega=omega)
            rows.append(row)

            print(
                f"sigma={sigma:.3g} omega={omega:+.3g} "
                f"Sp_abs={row['spearman_fluid_abs_delta_R']:+.4f} "
                f"p_abs={row['spearman_fluid_abs_delta_R_p']:.3g} "
                f"Sp_signed={row['spearman_fluid_delta_R']:+.4f} "
                f"p_signed={row['spearman_fluid_delta_R_p']:.3g}"
            )

    df = pd.DataFrame(rows)
    csv_path = os.path.join(args.output_dir, "fluid_source_scan.csv")
    df.to_csv(csv_path, index=False)

    best_abs = df.sort_values("spearman_fluid_abs_delta_R", ascending=False).iloc[0].to_dict()
    best_signed = df.sort_values("spearman_fluid_delta_R", ascending=False).iloc[0].to_dict()
    best_abs_p = df.sort_values(["spearman_fluid_abs_delta_R_p", "spearman_fluid_abs_delta_R"], ascending=[True, False]).iloc[0].to_dict()
    best_signed_p = df.sort_values(["spearman_fluid_delta_R_p", "spearman_fluid_delta_R"], ascending=[True, False]).iloc[0].to_dict()

    summary = {
        "title": "BuP Paper 8 — Fluid source scan T00 + omega Taa",
        "description": "Tests whether adding isotropic pressure Taa improves curvature response prediction.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "sigma_list": args.sigma_list,
        "omega_list": args.omega_list,
        "definition": {
            "T00": "T00(i)=1/2 sum_j deltaW_ij^2",
            "Taa": "Taa(i)=1/2 sum_j deltaW_ij^2 ||v_ij||^2",
            "S_fluid": "S_fluid(i)=T00(i)+omega*Taa(i)"
        },
        "best_abs_response_by_spearman": best_abs,
        "best_signed_response_by_spearman": best_signed,
        "best_abs_response_by_pvalue": best_abs_p,
        "best_signed_response_by_pvalue": best_signed_p,
        "files": {
            "scan_csv": "fluid_source_scan.csv",
            "fig_abs": "fig_fluid_abs_response_vs_sigma.png",
            "fig_signed": "fig_fluid_signed_response_vs_sigma.png",
            "fig_top20_abs": "fig_top20_fluid_abs_response.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    make_plots(df, args.output_dir)

    print("-" * 100)
    print("MEILLEURE RÉPONSE ABSOLUE |delta R|")
    print("-" * 100)
    print(f"sigma     : {best_abs['sigma']}")
    print(f"omega     : {best_abs['omega']}")
    print(f"Spearman  : {best_abs['spearman_fluid_abs_delta_R']:+.6f}")
    print(f"p         : {best_abs['spearman_fluid_abs_delta_R_p']:.6g}")
    print(f"Pearson   : {best_abs['pearson_fluid_abs_delta_R']:+.6f}")
    print(f"p         : {best_abs['pearson_fluid_abs_delta_R_p']:.6g}")

    print("-" * 100)
    print("MEILLEURE RÉPONSE SIGNÉE delta R")
    print("-" * 100)
    print(f"sigma     : {best_signed['sigma']}")
    print(f"omega     : {best_signed['omega']}")
    print(f"Spearman  : {best_signed['spearman_fluid_delta_R']:+.6f}")
    print(f"p         : {best_signed['spearman_fluid_delta_R_p']:.6g}")
    print(f"Pearson   : {best_signed['pearson_fluid_delta_R']:+.6f}")
    print(f"p         : {best_signed['pearson_fluid_delta_R_p']:.6g}")

    print("-" * 100)
    print("Fichiers générés")
    print("-" * 100)
    print(csv_path)
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_fluid_abs_response_vs_sigma.png"))
    print(os.path.join(args.output_dir, "fig_fluid_signed_response_vs_sigma.png"))
    print("DONE")


if __name__ == "__main__":
    main()

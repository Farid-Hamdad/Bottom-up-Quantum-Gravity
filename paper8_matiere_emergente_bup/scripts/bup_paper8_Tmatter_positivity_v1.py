#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_positivity_v1.py

Paper 8 — Test 1 : positivité de T_00^matter à partir de delta W_loc.

Objectif :
    Vérifier la première définition du tenseur de matière émergent :

        T_00^matter(i) = 1/2 * sum_j (delta W_ij)^2

    avec :

        M_matter = sum_i T_00^matter(i)
                 = sum_{i<j} (delta W_ij)^2

Ce test vérifie :
    1. T_00(i) >= 0 pour tous les noeuds.
    2. La normalisation de masse est correcte.
    3. L'excitation est bien localisée autour du centre choisi.
    4. Les fichiers de résultats sont exportés pour Paper 8.

Usage :
    python3 bup_paper8_Tmatter_positivity_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma 1.0 \
      --output-dir results_paper8_Tmatter_positivity_N20_lam057
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import shortest_path


# ============================================================
# Utilitaires
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


# ============================================================
# Graphe d'intrication
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
    strength = np.sum(W, axis=1)
    return int(np.argmax(strength))


# ============================================================
# Excitation locale
# ============================================================

def inject_excitation(W0, D0, amp=0.15, sigma=1.0, center=-1):
    W0 = clean_matrix(W0)
    n = W0.shape[0]

    if center is None or center < 0:
        center = choose_center(W0)

    d = D0[center].copy()
    finite = np.isfinite(d)

    if not np.any(finite):
        raise ValueError("Distances non finies : impossible d'injecter l'excitation.")

    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    f = f / max(np.max(f), 1e-15)

    vals = W0[W0 > 0]
    scale = float(np.mean(vals)) if vals.size else 1.0

    deltaW = amp * scale * np.outer(f, f)
    np.fill_diagonal(deltaW, 0.0)

    W_exc = clean_matrix(W0 + deltaW)

    # deltaW réel après nettoyage
    deltaW_actual = W_exc - W0
    np.fill_diagonal(deltaW_actual, 0.0)

    return W_exc, deltaW_actual, center, f


# ============================================================
# T_00 matter
# ============================================================

def compute_T00(deltaW):
    """
    Définition Paper 8 :

        T_00(i) = 1/2 * sum_j (deltaW_ij)^2

    Le facteur 1/2 évite le double comptage des arêtes.
    """
    deltaW = np.asarray(deltaW, dtype=float)
    T00 = 0.5 * np.sum(deltaW ** 2, axis=1)
    return T00


def mass_from_nodes(T00):
    return float(np.sum(T00))


def mass_from_edges(deltaW):
    n = deltaW.shape[0]
    return float(np.sum(np.triu(deltaW ** 2, 1)))


def localization_radius(T00, D0, center):
    total = np.sum(T00)

    if total <= 1e-15:
        return np.nan

    d = D0[center].copy()
    finite = np.isfinite(d)

    if not np.any(finite):
        return np.nan

    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    R2 = np.sum(T00 * d ** 2) / total
    return float(np.sqrt(max(R2, 0.0)))


def participation_ratio(T00):
    """
    Mesure grossière du nombre effectif de noeuds porteurs de l'énergie.
    """
    s1 = np.sum(T00)
    s2 = np.sum(T00 ** 2)

    if s2 <= 1e-15:
        return np.nan

    return float((s1 ** 2) / s2)


# ============================================================
# Figures
# ============================================================

def make_figures(node_df, output_dir, label):
    ensure_dir(output_dir)

    # T00 par noeud
    plt.figure(figsize=(9, 5))
    plt.bar(node_df["node"], node_df["T00"])
    plt.xlabel("Noeud")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title(f"Paper 8 — densité d'énergie émergente {label}")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_by_node.png"), dpi=250)
    plt.close()

    # T00 vs distance au centre
    plt.figure(figsize=(8, 5))
    plt.scatter(node_df["d_to_center"], node_df["T00"])
    plt.xlabel("Distance d'intrication au centre")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title("Localisation de la densité d'énergie")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_vs_distance.png"), dpi=250)
    plt.close()

    # Profil radial simple
    df = node_df.copy()
    df = df.sort_values("d_to_center")
    plt.figure(figsize=(8, 5))
    plt.plot(df["d_to_center"], df["T00"], "o-")
    plt.xlabel("Distance d'intrication au centre")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title("Profil radial discret de T00")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_radial_profile.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_positivity_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()

    ensure_dir(args.output_dir)

    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 92)
    print("BuP Paper 8 — Test de positivité de T_00^matter")
    print("=" * 92)
    print(f"MI file     : {args.mi_file}")
    print(f"label       : {label}")
    print(f"lambda      : {lam}")
    print(f"k           : {args.k}")
    print(f"amp         : {args.amp}")
    print(f"sigma       : {args.sigma}")
    print(f"output      : {args.output_dir}")
    print("-" * 92)

    # Charger W0
    W0 = load_matrix_csv(args.mi_file)
    n = W0.shape[0]

    # Géométrie de référence
    mask0 = knn_mask(W0, args.k)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)

    # Injection locale
    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=args.sigma,
        center=args.center
    )

    # T00
    T00 = compute_T00(deltaW)

    M_nodes = mass_from_nodes(T00)
    M_edges = mass_from_edges(deltaW)

    abs_err = abs(M_nodes - M_edges)
    rel_err = abs_err / max(abs(M_edges), 1e-15)

    rho_min = float(np.min(T00))
    rho_max = float(np.max(T00))
    rho_mean = float(np.mean(T00))
    rho_sum = float(np.sum(T00))

    R_loc = localization_radius(T00, D0, center)
    PR = participation_ratio(T00)

    positivity_pass = bool(rho_min >= -1e-14)
    normalization_pass = bool(rel_err < 1e-10)

    verdict = "PASS" if positivity_pass and normalization_pass else "FAIL"

    # Tables noeuds
    node_rows = []
    for i in range(n):
        node_rows.append({
            "node": i,
            "T00": float(T00[i]),
            "d_to_center": float(D0[center, i]),
            "profile_f": float(profile[i]),
            "strength_W0": float(np.sum(W0[i])),
            "strength_deltaW": float(np.sum(np.abs(deltaW[i]))),
        })

    node_df = pd.DataFrame(node_rows)
    node_df.to_csv(os.path.join(args.output_dir, "T00_node_table.csv"), index=False)

    # Sauvegardes matrices
    np.savetxt(os.path.join(args.output_dir, "W0.csv"), W0, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "W_excited.csv"), W_exc, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "deltaW_loc.csv"), deltaW, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "T00_vector.csv"), T00, delimiter=",")

    # Résumé
    summary = {
        "title": "BuP Paper 8 — Test de positivité de T_00^matter",
        "description": "Premier test de cohérence du tenseur énergie-impulsion émergent construit à partir de delta W_loc.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": n,
        "k": args.k,
        "amp": args.amp,
        "sigma": args.sigma,
        "center": center,
        "definition": {
            "T00": "T00(i) = 1/2 * sum_j (delta W_ij)^2",
            "mass_nodes": "M_nodes = sum_i T00(i)",
            "mass_edges": "M_edges = sum_{i<j} (delta W_ij)^2"
        },
        "results": {
            "T00_min": rho_min,
            "T00_max": rho_max,
            "T00_mean": rho_mean,
            "T00_sum": rho_sum,
            "M_from_nodes": M_nodes,
            "M_from_edges": M_edges,
            "normalization_abs_error": abs_err,
            "normalization_rel_error": rel_err,
            "localization_radius": R_loc,
            "participation_ratio": PR,
            "positivity_pass": positivity_pass,
            "normalization_pass": normalization_pass,
            "verdict": verdict
        },
        "interpretation": {
            "positivity": "T00(i) is non-negative because it is a sum of squared local entanglement perturbations.",
            "normalization": "The factor 1/2 avoids double counting and makes sum_i T00(i) equal to sum_{i<j} deltaW_ij^2.",
            "paper8_role": "This is the first consistency test for defining T_matter from delta W_loc."
        },
        "files": {
            "node_table": "T00_node_table.csv",
            "T00_vector": "T00_vector.csv",
            "deltaW": "deltaW_loc.csv",
            "figure_T00_by_node": "fig_T00_by_node.png",
            "figure_T00_vs_distance": "fig_T00_vs_distance.png",
            "figure_T00_radial_profile": "fig_T00_radial_profile.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    make_figures(node_df, args.output_dir, label)

    print("RÉSULTATS")
    print("-" * 92)
    print(f"center                         : {center}")
    print(f"T00 min                        : {rho_min:.12g}")
    print(f"T00 max                        : {rho_max:.12g}")
    print(f"T00 mean                       : {rho_mean:.12g}")
    print(f"Somme_i T00_i                  : {M_nodes:.12g}")
    print(f"Somme_i<j deltaW_ij^2           : {M_edges:.12g}")
    print(f"Erreur absolue normalisation    : {abs_err:.3e}")
    print(f"Erreur relative normalisation   : {rel_err:.3e}")
    print(f"Rayon localisation T00          : {R_loc:.12g}")
    print(f"Participation ratio             : {PR:.12g}")
    print("-" * 92)
    print(f"POSITIVITÉ                     : {'PASS' if positivity_pass else 'FAIL'}")
    print(f"NORMALISATION                  : {'PASS' if normalization_pass else 'FAIL'}")
    print(f"VERDICT                        : {verdict}")
    print("-" * 92)
    print("Fichiers générés :")
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "T00_node_table.csv"))
    print(os.path.join(args.output_dir, "fig_T00_by_node.png"))
    print(os.path.join(args.output_dir, "fig_T00_vs_distance.png"))
    print(os.path.join(args.output_dir, "fig_T00_radial_profile.png"))
    print("DONE")


if __name__ == "__main__":
    main()

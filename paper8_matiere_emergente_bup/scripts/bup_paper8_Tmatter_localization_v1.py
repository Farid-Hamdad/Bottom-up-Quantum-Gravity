#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_localization_v1.py

Paper 8 — Test 2 : localisation de T_00^matter.

Objectif :
    Vérifier que la densité d'énergie émergente

        T00(i) = 1/2 * sum_j (delta W_ij)^2

    est localisée autour du centre d'excitation.

Tests :
    1. T00 >= 0
    2. Normalisation correcte
    3. Corrélation négative entre distance au centre et T00
    4. Rayon de localisation
    5. Participation ratio
    6. Scan optionnel en sigma

Usage simple :
    python3 bup_paper8_Tmatter_localization_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma 0.5 \
      --output-dir results_paper8_Tmatter_localization_N20_lam057_sigma05

Usage scan sigma :
    python3 bup_paper8_Tmatter_localization_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.2 0.3 0.5 0.8 1.0 1.5 2.0 \
      --output-dir results_paper8_Tmatter_localization_scan
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.sparse.csgraph import shortest_path
from scipy.stats import spearmanr, pearsonr


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
# Excitation and T00
# ============================================================

def inject_excitation(W0, D0, amp=0.15, sigma=1.0, center=-1):
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


def compute_T00(deltaW):
    return 0.5 * np.sum(deltaW ** 2, axis=1)


def mass_from_edges(deltaW):
    return float(np.sum(np.triu(deltaW ** 2, 1)))


def localization_radius(T00, D0, center):
    total = np.sum(T00)
    if total <= 1e-15:
        return np.nan

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    return float(np.sqrt(np.sum(T00 * d ** 2) / total))


def participation_ratio(T00):
    s1 = np.sum(T00)
    s2 = np.sum(T00 ** 2)
    if s2 <= 1e-15:
        return np.nan
    return float((s1 ** 2) / s2)


def concentration_fraction(T00, D0, center, radius):
    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    total = np.sum(T00)
    if total <= 1e-15:
        return np.nan

    return float(np.sum(T00[d <= radius]) / total)


def run_one(W0, label, lam, args, sigma, output_dir=None):
    n = W0.shape[0]

    mask0 = knn_mask(W0, args.k)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    T00 = compute_T00(deltaW)
    d = D0[center].copy()

    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    M_nodes = float(np.sum(T00))
    M_edges = mass_from_edges(deltaW)
    rel_err = abs(M_nodes - M_edges) / max(abs(M_edges), 1e-15)

    pearson_r, pearson_p, spearman_r, spearman_p = safe_corr(d, T00)

    Rloc = localization_radius(T00, D0, center)
    PR = participation_ratio(T00)

    frac_Rloc = concentration_fraction(T00, D0, center, Rloc)
    frac_sigma = concentration_fraction(T00, D0, center, sigma)

    positivity_pass = bool(np.min(T00) >= -1e-14)
    normalization_pass = bool(rel_err < 1e-10)
    localization_pass = bool(np.isfinite(spearman_r) and spearman_r < 0)

    verdict = "PASS" if positivity_pass and normalization_pass and localization_pass else "CHECK"

    node_df = pd.DataFrame({
        "node": np.arange(n),
        "d_to_center": d,
        "T00": T00,
        "profile_f": profile,
        "strength_deltaW": np.sum(np.abs(deltaW), axis=1),
        "strength_W0": np.sum(W0, axis=1),
    }).sort_values("d_to_center")

    summary = {
        "label": label,
        "lambda": lam,
        "N": n,
        "k": args.k,
        "amp": args.amp,
        "sigma": sigma,
        "center": center,
        "T00_min": float(np.min(T00)),
        "T00_max": float(np.max(T00)),
        "T00_mean": float(np.mean(T00)),
        "M_from_nodes": M_nodes,
        "M_from_edges": M_edges,
        "normalization_rel_error": float(rel_err),
        "localization_radius": Rloc,
        "participation_ratio": PR,
        "fraction_inside_Rloc": frac_Rloc,
        "fraction_inside_sigma": frac_sigma,
        "pearson_d_vs_T00": pearson_r,
        "pearson_p": pearson_p,
        "spearman_d_vs_T00": spearman_r,
        "spearman_p": spearman_p,
        "positivity_pass": positivity_pass,
        "normalization_pass": normalization_pass,
        "localization_pass": localization_pass,
        "verdict": verdict,
    }

    if output_dir is not None:
        ensure_dir(output_dir)

        node_df.to_csv(os.path.join(output_dir, "T00_localization_node_table.csv"), index=False)
        np.savetxt(os.path.join(output_dir, "T00_vector.csv"), T00, delimiter=",")
        np.savetxt(os.path.join(output_dir, "deltaW_loc.csv"), deltaW, delimiter=",")

        with open(os.path.join(output_dir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=2)

        make_figures(node_df, output_dir, label, sigma, summary)

    return summary, node_df


# ============================================================
# Figures
# ============================================================

def make_figures(node_df, output_dir, label, sigma, summary):
    plt.figure(figsize=(8, 5))
    plt.scatter(node_df["d_to_center"], node_df["T00"])
    plt.xlabel("Distance d'intrication au centre")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title(f"T00 vs distance — {label}, sigma={sigma}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_vs_distance.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(node_df["d_to_center"], node_df["T00"], "o-")
    plt.xlabel("Distance d'intrication au centre")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title("Profil radial discret de T00")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_radial_profile.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.bar(node_df["node"], node_df["T00"])
    plt.xlabel("Noeud")
    plt.ylabel(r"$T_{00}^{matter}(i)$")
    plt.title("T00 par noeud")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_by_node.png"), dpi=250)
    plt.close()


def make_scan_figures(scan_df, output_dir):
    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["spearman_d_vs_T00"], "o-")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"Spearman$(d_{\rm ent},T_{00})$")
    plt.title("Localisation : corrélation distance / T00")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_scan_sigma_spearman.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["localization_radius"], "o-", label=r"$R_{\rm loc}$")
    plt.plot(scan_df["sigma"], scan_df["participation_ratio"], "s-", label="participation ratio")
    plt.xlabel(r"$\sigma$")
    plt.title("Rayon de localisation et participation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_scan_sigma_localization.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["M_from_nodes"], "o-")
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$M_{\rm matter}$")
    plt.title("Masse informationnelle vs sigma")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_scan_sigma_mass.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_localization_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--sigma-list", type=float, nargs="*", default=None)

    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    W0 = load_matrix_csv(args.mi_file)

    print("=" * 92)
    print("BuP Paper 8 — Test de localisation de T_00^matter")
    print("=" * 92)
    print(f"MI file : {args.mi_file}")
    print(f"label   : {label}")
    print(f"lambda  : {lam}")
    print(f"N       : {W0.shape[0]}")
    print(f"k       : {args.k}")
    print(f"amp     : {args.amp}")
    print(f"output  : {args.output_dir}")
    print("-" * 92)

    if args.sigma_list is not None and len(args.sigma_list) > 0:
        rows = []

        for sig in args.sigma_list:
            subdir = os.path.join(args.output_dir, f"sigma_{sig:.3g}".replace(".", "p"))
            summary, _ = run_one(W0, label, lam, args, sig, output_dir=subdir)
            rows.append(summary)

            print(
                f"sigma={sig:.3g} "
                f"Spearman={summary['spearman_d_vs_T00']:+.4f} "
                f"p={summary['spearman_p']:.3g} "
                f"Rloc={summary['localization_radius']:.4g} "
                f"PR={summary['participation_ratio']:.4g} "
                f"verdict={summary['verdict']}"
            )

        scan_df = pd.DataFrame(rows)
        scan_df.to_csv(os.path.join(args.output_dir, "localization_sigma_scan.csv"), index=False)

        global_summary = {
            "title": "BuP Paper 8 — Sigma scan for T00 localization",
            "description": "Checks whether T00 decreases with entanglement distance from excitation center.",
            "mi_file": args.mi_file,
            "label": label,
            "lambda": lam,
            "N": int(W0.shape[0]),
            "k": args.k,
            "amp": args.amp,
            "sigma_list": args.sigma_list,
            "best_localized": scan_df.sort_values("spearman_d_vs_T00").iloc[0].to_dict(),
            "files": {
                "scan_csv": "localization_sigma_scan.csv",
                "fig_spearman": "fig_scan_sigma_spearman.png",
                "fig_localization": "fig_scan_sigma_localization.png",
                "fig_mass": "fig_scan_sigma_mass.png"
            }
        }

        with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
            json.dump(global_summary, f, indent=2)

        make_scan_figures(scan_df, args.output_dir)

        print("-" * 92)
        print("RÉSUMÉ SCAN SIGMA")
        print(scan_df[[
            "sigma",
            "spearman_d_vs_T00",
            "spearman_p",
            "localization_radius",
            "participation_ratio",
            "M_from_nodes",
            "verdict"
        ]].to_string(index=False))
        print("-" * 92)
        print(f"Fichiers : {args.output_dir}/localization_sigma_scan.csv")
        print("DONE")

    else:
        summary, node_df = run_one(W0, label, lam, args, args.sigma, output_dir=args.output_dir)

        print("RÉSULTATS")
        print("-" * 92)
        print(f"sigma                         : {summary['sigma']}")
        print(f"center                        : {summary['center']}")
        print(f"T00 min                       : {summary['T00_min']:.12g}")
        print(f"T00 max                       : {summary['T00_max']:.12g}")
        print(f"M nodes                       : {summary['M_from_nodes']:.12g}")
        print(f"M edges                       : {summary['M_from_edges']:.12g}")
        print(f"Erreur rel normalisation      : {summary['normalization_rel_error']:.3e}")
        print(f"Rayon localisation            : {summary['localization_radius']:.12g}")
        print(f"Participation ratio           : {summary['participation_ratio']:.12g}")
        print(f"Pearson d vs T00              : {summary['pearson_d_vs_T00']:+.6f}  p={summary['pearson_p']:.3g}")
        print(f"Spearman d vs T00             : {summary['spearman_d_vs_T00']:+.6f}  p={summary['spearman_p']:.3g}")
        print("-" * 92)
        print(f"POSITIVITÉ                    : {'PASS' if summary['positivity_pass'] else 'FAIL'}")
        print(f"NORMALISATION                 : {'PASS' if summary['normalization_pass'] else 'FAIL'}")
        print(f"LOCALISATION                  : {'PASS' if summary['localization_pass'] else 'CHECK'}")
        print(f"VERDICT                       : {summary['verdict']}")
        print("-" * 92)
        print(f"Fichiers : {args.output_dir}/summary.json")
        print("DONE")


if __name__ == "__main__":
    main()

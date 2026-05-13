#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_curvature_response_v1.py

Paper 8 — Test 3 : T00^matter vs réponse de courbure nodale.

Objectif :
    Construire

        T00(i) = 1/2 * sum_j (delta W_ij)^2

    puis mesurer si cette densité d'énergie émergente est corrélée à la réponse
    de courbure locale :

        |delta R_i| = |R_i[W_excited] - R_i[W0]|

    où R_i est un scalaire de courbure nodal construit comme moyenne des
    courbures d'Ollivier-Ricci des arêtes incidentes.

Usage :
    python3 bup_paper8_Tmatter_curvature_response_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma 0.03 \
      --output-dir results_paper8_Tmatter_curvature_N20_lam057_sigma003
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
# Graph geometry
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
    """
    R_i = moyenne des kappa_ij sur les arêtes incidentes.
    """
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
        edge_rows.append({"i": i, "j": j, "kappa": kij, "W": W[i, j], "d_ent": D[i, j]})

    R = np.zeros(n)
    degree_curv = np.zeros(n, dtype=int)

    for i in range(n):
        degree_curv[i] = len(node_vals[i])
        R[i] = np.mean(node_vals[i]) if node_vals[i] else np.nan

    return R, pd.DataFrame(edge_rows), D, mask


# ============================================================
# Excitation and T00
# ============================================================

def inject_excitation(W0, D0, amp=0.15, sigma=0.03, center=-1):
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


# ============================================================
# Figures
# ============================================================

def make_figures(node_df, output_dir, label, sigma):
    plt.figure(figsize=(8, 5))
    plt.scatter(node_df["T00"], node_df["abs_delta_R"])
    plt.xlabel(r"$T_{00}^{matter}(i)$")
    plt.ylabel(r"$|\delta R_i|$")
    plt.title(f"T00 vs réponse de courbure — {label}, sigma={sigma}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_vs_abs_delta_R.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(node_df["T00"], node_df["delta_R"])
    plt.xlabel(r"$T_{00}^{matter}(i)$")
    plt.ylabel(r"$\delta R_i$")
    plt.title(f"T00 vs delta R signé — {label}, sigma={sigma}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_T00_vs_delta_R_signed.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.plot(node_df["node"], node_df["T00"], "o-", label=r"$T_{00}$")
    plt.plot(node_df["node"], node_df["abs_delta_R"], "s-", label=r"$|\delta R|$")
    plt.xlabel("Noeud")
    plt.title("Profil nodal : T00 et réponse de courbure")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_node_profiles_T00_deltaR.png"), dpi=250)
    plt.close()


def make_scan_figures(scan_df, output_dir):
    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["spearman_T00_abs_delta_R"], "o-", label=r"$T_{00}$ vs $|\delta R|$")
    plt.plot(scan_df["sigma"], scan_df["spearman_T00_delta_R"], "s-", label=r"$T_{00}$ vs $\delta R$")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel("Spearman")
    plt.title("Réponse de courbure vs largeur de l'excitation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_scan_sigma_curvature_response.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["localization_radius"], "o-", label=r"$R_{\rm loc}$")
    plt.plot(scan_df["sigma"], scan_df["participation_ratio"], "s-", label="PR")
    plt.xlabel(r"$\sigma$")
    plt.title("Localisation de T00")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_scan_sigma_localization.png"), dpi=250)
    plt.close()


# ============================================================
# Run one / scan
# ============================================================

def run_one(W0, args, sigma, output_dir=None):
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

    T00 = compute_T00(deltaW)
    delta_R = R1 - R0
    abs_delta_R = np.abs(delta_R)

    pr_abs, pp_abs, sr_abs, sp_abs = safe_corr(T00, abs_delta_R)
    pr_signed, pp_signed, sr_signed, sp_signed = safe_corr(T00, delta_R)

    Rloc = localization_radius(T00, D0, center)
    PR = participation_ratio(T00)

    positivity_pass = bool(np.min(T00) >= -1e-14)
    curvature_abs_pass = bool(np.isfinite(sr_abs) and sr_abs > 0)

    verdict = "PASS" if positivity_pass and curvature_abs_pass else "CHECK"

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    node_df = pd.DataFrame({
        "node": np.arange(n),
        "d_to_center": d,
        "T00": T00,
        "R0": R0,
        "R1": R1,
        "delta_R": delta_R,
        "abs_delta_R": abs_delta_R,
        "profile_f": profile,
        "strength_deltaW": np.sum(np.abs(deltaW), axis=1),
    })

    summary = {
        "label": label,
        "lambda": lam,
        "N": n,
        "k": args.k,
        "amp": args.amp,
        "sigma": sigma,
        "center": center,
        "definition": {
            "T00": "T00(i) = 1/2 * sum_j (delta W_ij)^2",
            "R_node": "R_i = mean Ollivier-Ricci curvature of incident edges",
            "delta_R": "R_i[W_excited] - R_i[W0]"
        },
        "T00_min": float(np.min(T00)),
        "T00_max": float(np.max(T00)),
        "T00_sum": float(np.sum(T00)),
        "delta_R_min": float(np.nanmin(delta_R)),
        "delta_R_max": float(np.nanmax(delta_R)),
        "abs_delta_R_mean": float(np.nanmean(abs_delta_R)),
        "localization_radius": Rloc,
        "participation_ratio": PR,
        "pearson_T00_abs_delta_R": pr_abs,
        "pearson_T00_abs_delta_R_p": pp_abs,
        "spearman_T00_abs_delta_R": sr_abs,
        "spearman_T00_abs_delta_R_p": sp_abs,
        "pearson_T00_delta_R": pr_signed,
        "pearson_T00_delta_R_p": pp_signed,
        "spearman_T00_delta_R": sr_signed,
        "spearman_T00_delta_R_p": sp_signed,
        "positivity_pass": positivity_pass,
        "curvature_abs_pass": curvature_abs_pass,
        "verdict": verdict
    }

    if output_dir is not None:
        ensure_dir(output_dir)

        node_df.to_csv(os.path.join(output_dir, "T00_curvature_node_table.csv"), index=False)
        edge0.to_csv(os.path.join(output_dir, "edge_curvature_W0.csv"), index=False)
        edge1.to_csv(os.path.join(output_dir, "edge_curvature_W_excited.csv"), index=False)

        np.savetxt(os.path.join(output_dir, "T00_vector.csv"), T00, delimiter=",")
        np.savetxt(os.path.join(output_dir, "deltaW_loc.csv"), deltaW, delimiter=",")
        np.savetxt(os.path.join(output_dir, "R0_vector.csv"), R0, delimiter=",")
        np.savetxt(os.path.join(output_dir, "R1_vector.csv"), R1, delimiter=",")
        np.savetxt(os.path.join(output_dir, "deltaR_vector.csv"), delta_R, delimiter=",")

        with open(os.path.join(output_dir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=2)

        make_figures(node_df, output_dir, label, sigma)

    return summary, node_df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_curvature_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=0.03)
    parser.add_argument("--sigma-list", type=float, nargs="*", default=None)

    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 94)
    print("BuP Paper 8 — T00^matter vs réponse de courbure nodale")
    print("=" * 94)
    print(f"MI file : {args.mi_file}")
    print(f"label   : {label}")
    print(f"lambda  : {lam}")
    print(f"N       : {W0.shape[0]}")
    print(f"k       : {args.k}")
    print(f"amp     : {args.amp}")
    print(f"output  : {args.output_dir}")
    print("-" * 94)

    if args.sigma_list is not None and len(args.sigma_list) > 0:
        rows = []

        for sig in args.sigma_list:
            subdir = os.path.join(args.output_dir, f"sigma_{sig:.3g}".replace(".", "p"))
            summary, _ = run_one(W0, args, sig, output_dir=subdir)
            rows.append(summary)

            print(
                f"sigma={sig:.3g} "
                f"Spearman(T00,|dR|)={summary['spearman_T00_abs_delta_R']:+.4f} "
                f"p={summary['spearman_T00_abs_delta_R_p']:.3g} "
                f"Spearman(T00,dR)={summary['spearman_T00_delta_R']:+.4f} "
                f"Rloc={summary['localization_radius']:.4g} "
                f"PR={summary['participation_ratio']:.4g} "
                f"verdict={summary['verdict']}"
            )

        scan_df = pd.DataFrame(rows)
        scan_df.to_csv(os.path.join(args.output_dir, "curvature_response_sigma_scan.csv"), index=False)

        global_summary = {
            "title": "BuP Paper 8 — T00 matter vs nodal curvature response",
            "description": "Tests whether emergent matter density predicts local curvature response.",
            "mi_file": args.mi_file,
            "label": label,
            "lambda": lam,
            "N": int(W0.shape[0]),
            "k": args.k,
            "amp": args.amp,
            "sigma_list": args.sigma_list,
            "best_abs_response": scan_df.sort_values("spearman_T00_abs_delta_R", ascending=False).iloc[0].to_dict(),
            "files": {
                "scan_csv": "curvature_response_sigma_scan.csv",
                "fig_response": "fig_scan_sigma_curvature_response.png",
                "fig_localization": "fig_scan_sigma_localization.png"
            }
        }

        with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
            json.dump(global_summary, f, indent=2)

        make_scan_figures(scan_df, args.output_dir)

        print("-" * 94)
        print("RÉSUMÉ SCAN SIGMA")
        print(scan_df[[
            "sigma",
            "spearman_T00_abs_delta_R",
            "spearman_T00_abs_delta_R_p",
            "spearman_T00_delta_R",
            "spearman_T00_delta_R_p",
            "localization_radius",
            "participation_ratio",
            "verdict"
        ]].to_string(index=False))
        print("-" * 94)
        print(f"Fichiers : {args.output_dir}/curvature_response_sigma_scan.csv")
        print("DONE")

    else:
        summary, _ = run_one(W0, args, args.sigma, output_dir=args.output_dir)

        print("RÉSULTATS")
        print("-" * 94)
        print(f"sigma                         : {summary['sigma']}")
        print(f"center                        : {summary['center']}")
        print(f"T00 min                       : {summary['T00_min']:.12g}")
        print(f"T00 max                       : {summary['T00_max']:.12g}")
        print(f"T00 sum                       : {summary['T00_sum']:.12g}")
        print(f"delta_R min                   : {summary['delta_R_min']:.12g}")
        print(f"delta_R max                   : {summary['delta_R_max']:.12g}")
        print(f"mean |delta_R|                : {summary['abs_delta_R_mean']:.12g}")
        print(f"Rayon localisation            : {summary['localization_radius']:.12g}")
        print(f"Participation ratio           : {summary['participation_ratio']:.12g}")
        print(f"Pearson T00 vs |dR|            : {summary['pearson_T00_abs_delta_R']:+.6f} p={summary['pearson_T00_abs_delta_R_p']:.3g}")
        print(f"Spearman T00 vs |dR|           : {summary['spearman_T00_abs_delta_R']:+.6f} p={summary['spearman_T00_abs_delta_R_p']:.3g}")
        print(f"Pearson T00 vs dR signé        : {summary['pearson_T00_delta_R']:+.6f} p={summary['pearson_T00_delta_R_p']:.3g}")
        print(f"Spearman T00 vs dR signé       : {summary['spearman_T00_delta_R']:+.6f} p={summary['spearman_T00_delta_R_p']:.3g}")
        print("-" * 94)
        print(f"POSITIVITÉ                    : {'PASS' if summary['positivity_pass'] else 'FAIL'}")
        print(f"RÉPONSE COURBURE |dR|          : {'PASS' if summary['curvature_abs_pass'] else 'CHECK'}")
        print(f"VERDICT                       : {summary['verdict']}")
        print("-" * 94)
        print(f"Fichiers : {args.output_dir}/summary.json")
        print("DONE")


if __name__ == "__main__":
    main()

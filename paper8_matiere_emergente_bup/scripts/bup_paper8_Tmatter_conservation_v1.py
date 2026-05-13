#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_Tmatter_conservation_v1.py

Paper 8 — Test de conservation discrète de la source tensorielle candidate.

But :
    Tester une condition de conservation spatiale effective associée au flux T_0a.

La conservation covariante complète serait :

    ∇_μ T^{μν} = 0

Mais nos tests actuels sont statiques : on ne dispose pas encore de ∂_t T^{0ν}.
On teste donc un proxy spatial :

    div f(i) ≈ 0

où :

    f_vec(i) = sum_j (deltaW_ij)^2 v_ij

et :

    v_ij = (x_j - x_i) / (d_ent(i,j) + eps)

Le script calcule :
    - T00
    - Taa
    - Tgrad
    - flux vectoriel F_vec
    - norme du flux ||F||
    - divergence discrète du flux
    - résidus relatifs de conservation

Usage :
    python3 bup_paper8_Tmatter_conservation_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.02 0.05 0.08 0.10 0.15 \
      --output-dir results_paper8_Tmatter_conservation_N20_lam057
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.sparse.csgraph import shortest_path
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


def norm01(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)

    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-15:
        return np.zeros_like(x)

    return (x - xmin) / (xmax - xmin)


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
# Excitation
# ============================================================

def inject_excitation(W0, D0, amp=0.15, sigma=0.05, center=-1):
    W0 = clean_matrix(W0)

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
    n = len(T00)
    Tgrad = np.zeros(n, dtype=float)

    for i in range(n):
        for j in range(n):
            wij = A[i, j]
            if wij <= 0:
                continue

            Tgrad[i] += wij * (T00[i] - T00[j]) ** 2

    return Tgrad


def compute_flux_vector(deltaW, D, X, eps=1e-12):
    """
    f_vec(i) = sum_j (deltaW_ij)^2 v_ij

    v_ij = (x_j - x_i)/(d_ij + eps)
    """
    n, ndim = X.shape
    flux_vec = np.zeros((n, ndim), dtype=float)

    for i in range(n):
        f = np.zeros(ndim, dtype=float)

        for j in range(n):
            if i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            f += (deltaW[i, j] ** 2) * v

        flux_vec[i] = f

    flux_norm = np.linalg.norm(flux_vec, axis=1)
    return flux_vec, flux_norm


# ============================================================
# Discrete divergence tests
# ============================================================

def divergence_flux_node(flux_vec, A, D, X, eps=1e-12):
    """
    Divergence discrète nodale du flux vectoriel :

        div f(i) = sum_j A_ij [f(j)-f(i)] · v_ij

    avec :

        v_ij = (x_j - x_i)/(d_ij + eps)

    Interprétation :
        Si div f ≈ 0, le flux est localement équilibré.
    """
    n = A.shape[0]
    div = np.zeros(n, dtype=float)

    for i in range(n):
        acc = 0.0

        for j in range(n):
            wij = A[i, j]
            if wij <= 0 or i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            acc += wij * np.dot(flux_vec[j] - flux_vec[i], v)

        div[i] = acc

    return div


def edge_antisymmetric_current(deltaW, D, X, A, eps=1e-12):
    """
    Courant orienté sur les arêtes :

        J_ij = (deltaW_ij)^2 * ||v_ij||

    puis divergence scalaire orientée approximative :

        divJ(i) = sum_j A_ij (J_ij - J_ji)

    Pour une deltaW symétrique purement statique, J_ij = J_ji,
    donc ce test peut être trivialement nul.
    Il est inclus comme diagnostic.
    """
    n = A.shape[0]
    J = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):
            wij = A[i, j]
            if wij <= 0 or i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            J[i, j] = (deltaW[i, j] ** 2) * np.linalg.norm(v)

    divJ = np.sum(A * (J - J.T), axis=1)
    return J, divJ


def conservation_metrics(div, flux_norm, eps=1e-15):
    l1_div = float(np.sum(np.abs(div)))
    l2_div = float(np.sqrt(np.sum(div ** 2)))
    linf_div = float(np.max(np.abs(div)))

    l1_flux = float(np.sum(np.abs(flux_norm)))
    l2_flux = float(np.sqrt(np.sum(flux_norm ** 2)))
    linf_flux = float(np.max(np.abs(flux_norm)))

    return {
        "div_l1": l1_div,
        "div_l2": l2_div,
        "div_linf": linf_div,
        "flux_l1": l1_flux,
        "flux_l2": l2_flux,
        "flux_linf": linf_flux,
        "rel_l1": l1_div / max(l1_flux, eps),
        "rel_l2": l2_div / max(l2_flux, eps),
        "rel_linf": linf_div / max(linf_flux, eps)
    }


# ============================================================
# One run
# ============================================================

def run_one(W0, args, sigma):
    mask0 = knn_mask(W0, args.k)
    A0 = adjacency_from_mask(W0, mask0)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)

    X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    T00 = compute_T00(deltaW)
    Taa = compute_Taa(deltaW, D0, X, eps=args.eps)
    Tgrad = compute_Tgrad(T00, A0)
    flux_vec, flux_norm = compute_flux_vector(deltaW, D0, X, eps=args.eps)

    div_flux = divergence_flux_node(flux_vec, A0, D0, X, eps=args.eps)
    Jedge, divJ = edge_antisymmetric_current(deltaW, D0, X, A0, eps=args.eps)

    metrics_flux = conservation_metrics(div_flux, flux_norm)
    metrics_edge = conservation_metrics(divJ, np.sum(np.abs(Jedge), axis=1))

    source_flux = T00 + args.omega * Taa + args.chi * Tgrad + args.psi * flux_norm

    if args.normalize_components:
        source_flux_normed = (
            norm01(T00)
            + args.omega * norm01(Taa)
            + args.chi * norm01(Tgrad)
            + args.psi * norm01(flux_norm)
        )
    else:
        source_flux_normed = source_flux

    node_df = pd.DataFrame({
        "node": np.arange(W0.shape[0]),
        "T00": T00,
        "Taa": Taa,
        "Tgrad": Tgrad,
        "flux_norm": flux_norm,
        "div_flux": div_flux,
        "abs_div_flux": np.abs(div_flux),
        "divJ_edge": divJ,
        "source_flux": source_flux,
        "source_flux_normed": source_flux_normed,
        "profile_f": profile,
    })

    summary = {
        "sigma": sigma,
        "center": int(center),
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "omega": args.omega,
        "chi": args.chi,
        "psi": args.psi,
        "normalize_components": args.normalize_components,

        "T00_sum": float(np.sum(T00)),
        "Taa_sum": float(np.sum(Taa)),
        "Tgrad_sum": float(np.sum(Tgrad)),
        "flux_norm_sum": float(np.sum(flux_norm)),

        "Taa_over_T00": float(np.sum(Taa) / max(np.sum(T00), 1e-15)),
        "Tgrad_over_T00": float(np.sum(Tgrad) / max(np.sum(T00), 1e-15)),
        "flux_over_T00": float(np.sum(flux_norm) / max(np.sum(T00), 1e-15)),

        "div_flux_l1": metrics_flux["div_l1"],
        "div_flux_l2": metrics_flux["div_l2"],
        "div_flux_linf": metrics_flux["div_linf"],
        "flux_l1": metrics_flux["flux_l1"],
        "flux_l2": metrics_flux["flux_l2"],
        "flux_linf": metrics_flux["flux_linf"],
        "conservation_rel_l1": metrics_flux["rel_l1"],
        "conservation_rel_l2": metrics_flux["rel_l2"],
        "conservation_rel_linf": metrics_flux["rel_linf"],

        "edge_div_l1": metrics_edge["div_l1"],
        "edge_div_l2": metrics_edge["div_l2"],
        "edge_conservation_rel_l1": metrics_edge["rel_l1"],
        "edge_conservation_rel_l2": metrics_edge["rel_l2"],
    }

    # Verdict qualitatif
    rel = summary["conservation_rel_l2"]

    if rel < 1e-2:
        verdict = "APPROX_CONSERVED"
    elif rel < 1e-1:
        verdict = "WEAKLY_VIOLATED"
    else:
        verdict = "NOT_CONSERVED"

    summary["verdict_spatial_flux_conservation"] = verdict

    return summary, node_df


# ============================================================
# Figures
# ============================================================

def make_figures(scan_df, output_dir):
    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["conservation_rel_l2"], "o-", label="rel L2")
    plt.plot(scan_df["sigma"], scan_df["conservation_rel_l1"], "s-", label="rel L1")
    plt.axhline(1e-2, linestyle="--", linewidth=1, label="1e-2")
    plt.axhline(1e-1, linestyle="--", linewidth=1, label="1e-1")
    plt.xlabel(r"$\sigma$")
    plt.ylabel("résidu relatif de divergence")
    plt.title("Test de conservation spatiale du flux")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_conservation_residual_vs_sigma.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(scan_df["sigma"], scan_df["flux_norm_sum"], "o-", label=r"$\sum_i ||f_i||$")
    plt.plot(scan_df["sigma"], scan_df["div_flux_l1"], "s-", label=r"$\sum_i |\mathrm{div} f_i|$")
    plt.xlabel(r"$\sigma$")
    plt.title("Flux total et divergence totale")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_flux_and_divergence_vs_sigma.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_Tmatter_conservation_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)

    # meilleurs poids trouvés précédemment
    parser.add_argument("--omega", type=float, default=-0.5)
    parser.add_argument("--chi", type=float, default=0.5)
    parser.add_argument("--psi", type=float, default=1.0)

    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--eps", type=float, default=1e-12)

    parser.add_argument("--ndim", type=int, default=3)
    parser.add_argument("--random-state", type=int, default=0)

    parser.add_argument("--normalize-components", action="store_true")

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 100)
    print("BuP Paper 8 — Test de conservation discrète")
    print("=" * 100)
    print(f"MI file     : {args.mi_file}")
    print(f"label       : {label}")
    print(f"lambda      : {lam}")
    print(f"N           : {W0.shape[0]}")
    print(f"k           : {args.k}")
    print(f"amp         : {args.amp}")
    print(f"omega       : {args.omega}")
    print(f"chi         : {args.chi}")
    print(f"psi         : {args.psi}")
    print(f"normalize   : {args.normalize_components}")
    print(f"sigma list  : {args.sigma_list}")
    print(f"output      : {args.output_dir}")
    print("-" * 100)

    rows = []

    for sigma in args.sigma_list:
        summary, node_df = run_one(W0, args, sigma)
        rows.append(summary)

        subdir = os.path.join(args.output_dir, f"sigma_{sigma:.3g}".replace(".", "p"))
        ensure_dir(subdir)
        node_df.to_csv(os.path.join(subdir, "conservation_node_table.csv"), index=False)

        print(
            f"sigma={sigma:.3g} "
            f"rel_L2={summary['conservation_rel_l2']:.4g} "
            f"rel_L1={summary['conservation_rel_l1']:.4g} "
            f"rel_Linf={summary['conservation_rel_linf']:.4g} "
            f"verdict={summary['verdict_spatial_flux_conservation']}"
        )

    scan_df = pd.DataFrame(rows)
    csv_path = os.path.join(args.output_dir, "conservation_sigma_scan.csv")
    scan_df.to_csv(csv_path, index=False)

    best = scan_df.sort_values("conservation_rel_l2", ascending=True).iloc[0].to_dict()
    worst = scan_df.sort_values("conservation_rel_l2", ascending=False).iloc[0].to_dict()

    summary_global = {
        "title": "BuP Paper 8 — Discrete conservation test",
        "description": "Tests spatial divergence of the flux component T0a. This is not a full covariant conservation test.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "omega": args.omega,
        "chi": args.chi,
        "psi": args.psi,
        "normalize_components": args.normalize_components,
        "sigma_list": args.sigma_list,
        "definition": {
            "flux": "f_vec(i)=sum_j deltaW_ij^2 v_ij",
            "divergence": "div f(i)=sum_j A_ij [f(j)-f(i)] dot v_ij",
            "relative_residual_L2": "||div f||_2 / ||f||_2",
            "status": "spatial static proxy only; full covariant conservation requires time evolution"
        },
        "best_conserved_by_L2": best,
        "worst_conserved_by_L2": worst,
        "files": {
            "scan_csv": "conservation_sigma_scan.csv",
            "fig_residual": "fig_conservation_residual_vs_sigma.png",
            "fig_flux_divergence": "fig_flux_and_divergence_vs_sigma.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary_global, f, indent=2)

    make_figures(scan_df, args.output_dir)

    print("-" * 100)
    print("MEILLEUR CAS — conservation spatiale du flux")
    print("-" * 100)
    print(f"sigma      : {best['sigma']}")
    print(f"rel_L2     : {best['conservation_rel_l2']:.6g}")
    print(f"rel_L1     : {best['conservation_rel_l1']:.6g}")
    print(f"verdict    : {best['verdict_spatial_flux_conservation']}")

    print("-" * 100)
    print("Fichiers générés")
    print("-" * 100)
    print(csv_path)
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_conservation_residual_vs_sigma.png"))
    print("DONE")


if __name__ == "__main__":
    main()

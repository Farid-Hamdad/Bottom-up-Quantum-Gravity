#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_total_conservation_test_v1.py

Paper 8 — Test de conservation totale matière + intrication.

Idée :
    La matière seule n'est pas forcément conservée dans BuP :

        div T_matter != 0

    car elle peut échanger avec le fond d'intrication.

    La conservation attendue est plutôt :

        div( T_matter + T_ent ) ≈ 0

Proxy testé :
    - Flux matière :
        f_matter(i) = sum_j (deltaW_ij)^2 v_ij

    - Champ entropique local :
        phi_ent(i) = delta d_s(i)
                  = d_s[W_exc](i) - d_s[W0](i)

    - Flux d'intrication :
        f_ent(i) = - sum_j A_ij (phi_j - phi_i) v_ij

    - Flux total :
        f_total(i) = f_matter(i) + alpha f_ent(i)

    Le script cherche alpha_opt qui minimise :

        || div f_total ||_2

Interprétation :
    Si le résidu total devient beaucoup plus petit que le résidu matière seul,
    alors la non-conservation de la matière peut être compensée par un flux
    d'intrication. Cela soutient l'idée d'une conservation totale BuP.

Usage :
    python3 bup_paper8_total_conservation_test_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.02 0.05 0.08 0.10 0.15 \
      --tau-min 0.01 \
      --tau-max 50 \
      --tau-points 20 \
      --output-dir results_paper8_total_conservation_N20_lam057
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import eigh
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
# Spectral dimension locale
# ============================================================

def normalized_laplacian(A, eps=1e-12):
    A = clean_matrix(A)
    deg = np.sum(A, axis=1)
    inv_sqrt = 1.0 / np.sqrt(deg + eps)
    L = np.eye(A.shape[0]) - (inv_sqrt[:, None] * A * inv_sqrt[None, :])
    L = 0.5 * (L + L.T)
    return L


def local_spectral_dimension(W, k, tau_min=0.01, tau_max=50.0, tau_points=20, eps=1e-12):
    """
    Calcule un proxy local de dimension spectrale.

    K_i(tau) = [exp(-tau L)]_ii

    d_s(i) = -2 * pente de log K_i(tau) vs log tau
    """
    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    L = normalized_laplacian(A, eps=eps)

    evals, evecs = eigh(L)
    evals = np.maximum(evals, 0.0)

    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_points)
    log_tau = np.log(taus)

    n = W.shape[0]
    Kdiag = np.zeros((tau_points, n), dtype=float)

    # diag(exp(-tau L)) = sum_m exp(-tau lambda_m) * v_im^2
    V2 = evecs ** 2

    for a, tau in enumerate(taus):
        weights = np.exp(-tau * evals)
        Kdiag[a, :] = V2 @ weights

    Kdiag = np.maximum(Kdiag, eps)

    ds = np.zeros(n, dtype=float)

    for i in range(n):
        y = np.log(Kdiag[:, i])
        slope, intercept = np.polyfit(log_tau, y, 1)
        ds[i] = -2.0 * slope

    return ds


# ============================================================
# Fluxes and divergence
# ============================================================

def compute_matter_flux(deltaW, D, X, eps=1e-12):
    """
    f_matter(i) = sum_j (deltaW_ij)^2 v_ij
    """
    n, ndim = X.shape
    flux = np.zeros((n, ndim), dtype=float)

    for i in range(n):
        f = np.zeros(ndim)

        for j in range(n):
            if i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            f += (deltaW[i, j] ** 2) * v

        flux[i] = f

    return flux


def compute_ent_flux(phi, A, D, X, eps=1e-12):
    """
    f_ent(i) = - sum_j A_ij (phi_j - phi_i) v_ij

    C'est un flux de compensation porté par le gradient du champ entropique phi.
    Ici phi = delta d_s.
    """
    n, ndim = X.shape
    flux = np.zeros((n, ndim), dtype=float)

    for i in range(n):
        f = np.zeros(ndim)

        for j in range(n):
            wij = A[i, j]
            if wij <= 0 or i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            f += -wij * (phi[j] - phi[i]) * v

        flux[i] = f

    return flux


def divergence_flux(flux_vec, A, D, X, eps=1e-12):
    """
    div f(i) = sum_j A_ij [f(j)-f(i)] · v_ij
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


def norm_l2(x):
    return float(np.sqrt(np.sum(np.asarray(x, dtype=float) ** 2)))


def norm_l1(x):
    return float(np.sum(np.abs(np.asarray(x, dtype=float))))


def conservation_residual(div_total, div_ref, eps=1e-15):
    return {
        "total_div_l1": norm_l1(div_total),
        "total_div_l2": norm_l2(div_total),
        "matter_div_l1": norm_l1(div_ref),
        "matter_div_l2": norm_l2(div_ref),
        "rel_l1": norm_l1(div_total) / max(norm_l1(div_ref), eps),
        "rel_l2": norm_l2(div_total) / max(norm_l2(div_ref), eps),
    }


def optimal_alpha(div_matter, div_ent, eps=1e-15):
    """
    Minimise || div_matter + alpha div_ent ||_2.
    """
    denom = float(np.dot(div_ent, div_ent))

    if denom < eps:
        return 0.0

    return float(-np.dot(div_matter, div_ent) / denom)


# ============================================================
# One sigma run
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

    # matière
    f_matter = compute_matter_flux(deltaW, D0, X, eps=args.eps)
    div_matter = divergence_flux(f_matter, A0, D0, X, eps=args.eps)

    # champ entropique phi = delta d_s local
    ds0 = local_spectral_dimension(
        W0,
        args.k,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        tau_points=args.tau_points,
        eps=args.eps
    )

    ds1 = local_spectral_dimension(
        W_exc,
        args.k,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        tau_points=args.tau_points,
        eps=args.eps
    )

    delta_ds = ds1 - ds0

    # flux entropique
    f_ent = compute_ent_flux(delta_ds, A0, D0, X, eps=args.eps)
    div_ent = divergence_flux(f_ent, A0, D0, X, eps=args.eps)

    # meilleur alpha analytique
    alpha_opt = optimal_alpha(div_matter, div_ent)
    div_total_opt = div_matter + alpha_opt * div_ent
    f_total_opt = f_matter + alpha_opt * f_ent

    metrics_opt = conservation_residual(div_total_opt, div_matter)

    # scan alpha optionnel
    scan_rows = []
    for alpha in args.alpha_list:
        div_total = div_matter + alpha * div_ent
        f_total = f_matter + alpha * f_ent
        m = conservation_residual(div_total, div_matter)

        scan_rows.append({
            "sigma": sigma,
            "alpha": alpha,
            "total_div_l1": m["total_div_l1"],
            "total_div_l2": m["total_div_l2"],
            "rel_l1": m["rel_l1"],
            "rel_l2": m["rel_l2"],
            "matter_div_l1": m["matter_div_l1"],
            "matter_div_l2": m["matter_div_l2"],
            "flux_total_l2": norm_l2(np.linalg.norm(f_total, axis=1)),
        })

    scan_df = pd.DataFrame(scan_rows)

    best_scan = scan_df.sort_values("rel_l2").iloc[0].to_dict()

    node_df = pd.DataFrame({
        "node": np.arange(W0.shape[0]),
        "profile_f": profile,
        "ds0": ds0,
        "ds1": ds1,
        "delta_ds": delta_ds,
        "div_matter": div_matter,
        "div_ent": div_ent,
        "div_total_opt": div_total_opt,
        "matter_flux_norm": np.linalg.norm(f_matter, axis=1),
        "ent_flux_norm": np.linalg.norm(f_ent, axis=1),
        "total_flux_norm_opt": np.linalg.norm(f_total_opt, axis=1),
    })

    summary = {
        "sigma": sigma,
        "center": int(center),

        "ds0_mean": float(np.mean(ds0)),
        "ds1_mean": float(np.mean(ds1)),
        "delta_ds_mean": float(np.mean(delta_ds)),
        "delta_ds_min": float(np.min(delta_ds)),
        "delta_ds_max": float(np.max(delta_ds)),

        "matter_div_l1": norm_l1(div_matter),
        "matter_div_l2": norm_l2(div_matter),
        "ent_div_l1": norm_l1(div_ent),
        "ent_div_l2": norm_l2(div_ent),

        "alpha_opt": alpha_opt,

        "total_opt_div_l1": metrics_opt["total_div_l1"],
        "total_opt_div_l2": metrics_opt["total_div_l2"],
        "rel_l1_opt": metrics_opt["rel_l1"],
        "rel_l2_opt": metrics_opt["rel_l2"],

        "best_scan_alpha": best_scan["alpha"],
        "best_scan_rel_l2": best_scan["rel_l2"],
        "best_scan_total_div_l2": best_scan["total_div_l2"],

        "improvement_factor_l2": norm_l2(div_matter) / max(metrics_opt["total_div_l2"], 1e-15),
    }

    if summary["rel_l2_opt"] < 0.1:
        verdict = "TOTAL_APPROX_CONSERVED"
    elif summary["rel_l2_opt"] < 0.5:
        verdict = "TOTAL_PARTIALLY_COMPENSATED"
    else:
        verdict = "TOTAL_NOT_COMPENSATED"

    summary["verdict_total_conservation"] = verdict

    return summary, node_df, scan_df


# ============================================================
# Figures
# ============================================================

def make_figures(global_df, output_dir):
    plt.figure(figsize=(8, 5))
    plt.plot(global_df["sigma"], global_df["matter_div_l2"], "o-", label="matter div L2")
    plt.plot(global_df["sigma"], global_df["total_opt_div_l2"], "s-", label="total div L2 opt")
    plt.xlabel(r"$\sigma$")
    plt.ylabel("norme L2 divergence")
    plt.title("Conservation totale : divergence matière vs totale")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_total_conservation_div_l2.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(global_df["sigma"], global_df["rel_l2_opt"], "o-")
    plt.axhline(1.0, linewidth=1)
    plt.axhline(0.5, linestyle="--", linewidth=1)
    plt.axhline(0.1, linestyle="--", linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$||\mathrm{div}_{tot}||_2 / ||\mathrm{div}_{matter}||_2$")
    plt.title("Résidu relatif après compensation entropique")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_total_conservation_relative_residual.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(global_df["sigma"], global_df["alpha_opt"], "o-")
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$\alpha_{\rm opt}$")
    plt.title("Coefficient optimal du flux d'intrication")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_alpha_opt_vs_sigma.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_total_conservation_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--tau-points", type=int, default=20)

    parser.add_argument("--alpha-list", type=float, nargs="*", default=[-10, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 10])

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
    print("BuP Paper 8 — Test conservation totale matière + intrication")
    print("=" * 100)
    print(f"MI file     : {args.mi_file}")
    print(f"label       : {label}")
    print(f"lambda      : {lam}")
    print(f"N           : {W0.shape[0]}")
    print(f"k           : {args.k}")
    print(f"amp         : {args.amp}")
    print(f"sigma list  : {args.sigma_list}")
    print(f"tau window  : [{args.tau_min}, {args.tau_max}] with {args.tau_points} points")
    print(f"output      : {args.output_dir}")
    print("-" * 100)

    global_rows = []

    for sigma in args.sigma_list:
        summary, node_df, alpha_scan_df = run_one(W0, args, sigma)
        global_rows.append(summary)

        subdir = os.path.join(args.output_dir, f"sigma_{sigma:.3g}".replace(".", "p"))
        ensure_dir(subdir)

        node_df.to_csv(os.path.join(subdir, "total_conservation_node_table.csv"), index=False)
        alpha_scan_df.to_csv(os.path.join(subdir, "alpha_scan.csv"), index=False)

        print(
            f"sigma={sigma:.3g} "
            f"matter_L2={summary['matter_div_l2']:.4g} "
            f"total_L2={summary['total_opt_div_l2']:.4g} "
            f"rel={summary['rel_l2_opt']:.4g} "
            f"alpha_opt={summary['alpha_opt']:+.4g} "
            f"improv={summary['improvement_factor_l2']:.3g} "
            f"verdict={summary['verdict_total_conservation']}"
        )

    global_df = pd.DataFrame(global_rows)
    csv_path = os.path.join(args.output_dir, "total_conservation_sigma_scan.csv")
    global_df.to_csv(csv_path, index=False)

    best = global_df.sort_values("rel_l2_opt").iloc[0].to_dict()

    summary_global = {
        "title": "BuP Paper 8 — Total conservation test",
        "description": "Tests whether matter non-conservation can be compensated by an entanglement flux built from delta local spectral dimension.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "sigma_list": args.sigma_list,
        "tau_min": args.tau_min,
        "tau_max": args.tau_max,
        "tau_points": args.tau_points,
        "definition": {
            "matter_flux": "f_matter(i)=sum_j deltaW_ij^2 v_ij",
            "phi_ent": "delta d_s(i)=d_s[W_exc](i)-d_s[W0](i)",
            "ent_flux": "f_ent(i)=-sum_j A_ij (phi_j-phi_i) v_ij",
            "total_flux": "f_total=f_matter+alpha f_ent",
            "alpha_opt": "least-squares minimizer of ||div(f_matter+alpha f_ent)||_2"
        },
        "best_total_conservation_by_rel_l2": best,
        "files": {
            "scan_csv": "total_conservation_sigma_scan.csv",
            "fig_div_l2": "fig_total_conservation_div_l2.png",
            "fig_relative": "fig_total_conservation_relative_residual.png",
            "fig_alpha": "fig_alpha_opt_vs_sigma.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary_global, f, indent=2)

    make_figures(global_df, args.output_dir)

    print("-" * 100)
    print("MEILLEUR CAS — conservation totale")
    print("-" * 100)
    print(f"sigma       : {best['sigma']}")
    print(f"alpha_opt   : {best['alpha_opt']:+.6g}")
    print(f"matter L2   : {best['matter_div_l2']:.6g}")
    print(f"total L2    : {best['total_opt_div_l2']:.6g}")
    print(f"rel L2      : {best['rel_l2_opt']:.6g}")
    print(f"improvement : {best['improvement_factor_l2']:.6g}")
    print(f"verdict     : {best['verdict_total_conservation']}")

    print("-" * 100)
    print("Fichiers générés")
    print("-" * 100)
    print(csv_path)
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_total_conservation_relative_residual.png"))
    print("DONE")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper8_total_conservation_test_v2.py

Paper 8 — Test de conservation totale v2.

But :
    Tester si la non-conservation de la source matière candidate peut être
    compensée par un flux d'intrication plus riche.

Conservation cherchée :
    div( f_matter + alpha f_ent ) ≈ 0

Flux matière :
    f_matter(i) = sum_j (deltaW_ij)^2 v_ij

Flux d'intrication testés :
    1. phi = delta d_s
    2. phi = delta R
    3. phi = delta d_s + beta * delta R
    4. phi = delta d_s + beta * delta R + gamma * Lap(delta R)

Flux associé :
    f_ent(i) = - sum_j A_ij (phi_j - phi_i) v_ij

Pour chaque proxy, le script trouve alpha_opt qui minimise :
    || div(f_matter + alpha f_ent) ||_2

Résultat à battre :
    v1 avec phi=delta d_s :
        rel_L2 ≈ 0.668

Usage :
    python3 bup_paper8_total_conservation_test_v2.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.02 0.05 0.08 0.10 0.15 \
      --tau-min 0.01 \
      --tau-max 50 \
      --tau-points 20 \
      --beta-list -3 -1 -0.5 0 0.5 1 3 \
      --gamma-list -3 -1 -0.5 0 0.5 1 3 \
      --output-dir results_paper8_total_conservation_v2_N20_lam057
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


def norm_l2(x):
    return float(np.sqrt(np.sum(np.asarray(x, dtype=float) ** 2)))


def norm_l1(x):
    return float(np.sum(np.abs(np.asarray(x, dtype=float))))


def standardize(x, eps=1e-15):
    """
    Standardisation robuste pour combiner delta d_s, delta R, Lap(delta R).
    """
    x = np.asarray(x, dtype=float)
    mu = np.nanmean(x)
    sig = np.nanstd(x)
    if sig < eps:
        return np.zeros_like(x)
    return (x - mu) / sig


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
# Local spectral dimension
# ============================================================

def normalized_laplacian(A, eps=1e-12):
    A = clean_matrix(A)
    deg = np.sum(A, axis=1)
    inv_sqrt = 1.0 / np.sqrt(deg + eps)
    L = np.eye(A.shape[0]) - (inv_sqrt[:, None] * A * inv_sqrt[None, :])
    L = 0.5 * (L + L.T)
    return L


def local_spectral_dimension(W, k, tau_min=0.01, tau_max=50.0, tau_points=20, eps=1e-12):
    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    L = normalized_laplacian(A, eps=eps)

    evals, evecs = eigh(L)
    evals = np.maximum(evals, 0.0)

    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_points)
    log_tau = np.log(taus)

    n = W.shape[0]
    Kdiag = np.zeros((tau_points, n), dtype=float)
    V2 = evecs ** 2

    for a, tau in enumerate(taus):
        weights = np.exp(-tau * evals)
        Kdiag[a, :] = V2 @ weights

    Kdiag = np.maximum(Kdiag, eps)

    ds = np.zeros(n, dtype=float)
    for i in range(n):
        y = np.log(Kdiag[:, i])
        slope, _ = np.polyfit(log_tau, y, 1)
        ds[i] = -2.0 * slope

    return ds


# ============================================================
# Ollivier-Ricci curvature proxy
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
    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    D = entanglement_distance(W, mask, eps=eps)
    kappas = ollivier_kappa(A, D, mask)

    n = W.shape[0]
    node_vals = [[] for _ in range(n)]
    for (i, j), kij in kappas.items():
        node_vals[i].append(kij)
        node_vals[j].append(kij)

    R = np.zeros(n)
    for i in range(n):
        R[i] = np.mean(node_vals[i]) if node_vals[i] else np.nan

    R = np.nan_to_num(R, nan=0.0, posinf=0.0, neginf=0.0)
    return R


# ============================================================
# Laplacian of scalar
# ============================================================

def graph_laplacian_scalar(phi, A, eps=1e-12):
    """
    Lap(phi)(i) = sum_j A_ij (phi_j - phi_i)
    """
    phi = np.asarray(phi, dtype=float)
    n = len(phi)
    out = np.zeros(n, dtype=float)

    for i in range(n):
        s = 0.0
        for j in range(n):
            wij = A[i, j]
            if wij <= 0 or i == j:
                continue
            s += wij * (phi[j] - phi[i])
        out[i] = s

    return out


# ============================================================
# Flux and divergence
# ============================================================

def compute_matter_flux(deltaW, D, X, eps=1e-12):
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


def compute_ent_flux_from_phi(phi, A, D, X, eps=1e-12):
    """
    f_ent(i) = - sum_j A_ij (phi_j - phi_i) v_ij
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


def optimal_alpha(div_matter, div_ent, eps=1e-15):
    denom = float(np.dot(div_ent, div_ent))
    if denom < eps:
        return 0.0
    return float(-np.dot(div_matter, div_ent) / denom)


def compensation_metrics(div_matter, div_ent, eps=1e-15):
    alpha = optimal_alpha(div_matter, div_ent, eps=eps)
    div_total = div_matter + alpha * div_ent

    matter_l1 = norm_l1(div_matter)
    matter_l2 = norm_l2(div_matter)
    total_l1 = norm_l1(div_total)
    total_l2 = norm_l2(div_total)

    return {
        "alpha_opt": alpha,
        "matter_div_l1": matter_l1,
        "matter_div_l2": matter_l2,
        "ent_div_l1": norm_l1(div_ent),
        "ent_div_l2": norm_l2(div_ent),
        "total_div_l1": total_l1,
        "total_div_l2": total_l2,
        "rel_l1": total_l1 / max(matter_l1, eps),
        "rel_l2": total_l2 / max(matter_l2, eps),
        "improvement_factor_l2": matter_l2 / max(total_l2, eps)
    }


# ============================================================
# One sigma run
# ============================================================

def run_one(W0, args, sigma):
    # geometry from background
    mask0 = knn_mask(W0, args.k)
    A0 = adjacency_from_mask(W0, mask0)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)
    X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)

    # excitation
    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    # matter divergence
    f_matter = compute_matter_flux(deltaW, D0, X, eps=args.eps)
    div_matter = divergence_flux(f_matter, A0, D0, X, eps=args.eps)

    # geometric fields
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

    R0 = node_curvature_R(W0, args.k, eps=args.eps)
    R1 = node_curvature_R(W_exc, args.k, eps=args.eps)
    delta_R = R1 - R0

    lap_delta_R = graph_laplacian_scalar(delta_R, A0, eps=args.eps)

    # Standardized fields for combinations
    z_ds = standardize(delta_ds)
    z_R = standardize(delta_R)
    z_lapR = standardize(lap_delta_R)

    proxy_rows = []

    def eval_proxy(proxy_name, phi, beta=None, gamma=None):
        f_ent = compute_ent_flux_from_phi(phi, A0, D0, X, eps=args.eps)
        div_ent = divergence_flux(f_ent, A0, D0, X, eps=args.eps)
        m = compensation_metrics(div_matter, div_ent, eps=args.eps)

        row = {
            "sigma": sigma,
            "proxy": proxy_name,
            "beta": beta,
            "gamma": gamma,
            "center": int(center),
            "delta_ds_mean": float(np.mean(delta_ds)),
            "delta_ds_min": float(np.min(delta_ds)),
            "delta_ds_max": float(np.max(delta_ds)),
            "delta_R_mean": float(np.mean(delta_R)),
            "delta_R_min": float(np.min(delta_R)),
            "delta_R_max": float(np.max(delta_R)),
            "lap_delta_R_mean": float(np.mean(lap_delta_R)),
            "lap_delta_R_min": float(np.min(lap_delta_R)),
            "lap_delta_R_max": float(np.max(lap_delta_R)),
        }
        row.update(m)

        if row["rel_l2"] < 0.1:
            verdict = "TOTAL_APPROX_CONSERVED"
        elif row["rel_l2"] < 0.5:
            verdict = "TOTAL_PARTIALLY_COMPENSATED"
        else:
            verdict = "TOTAL_NOT_COMPENSATED"

        row["verdict"] = verdict
        proxy_rows.append(row)

    # 1. delta ds
    eval_proxy("delta_ds", z_ds)

    # 2. delta R
    eval_proxy("delta_R", z_R)

    # 3. lap delta R
    eval_proxy("lap_delta_R", z_lapR)

    # 4. delta ds + beta delta R
    for beta in args.beta_list:
        phi = z_ds + beta * z_R
        eval_proxy("delta_ds_plus_beta_delta_R", phi, beta=beta)

    # 5. delta R + gamma lapR
    for gamma in args.gamma_list:
        phi = z_R + gamma * z_lapR
        eval_proxy("delta_R_plus_gamma_lapR", phi, gamma=gamma)

    # 6. delta ds + beta delta R + gamma lapR
    for beta in args.beta_list:
        for gamma in args.gamma_list:
            phi = z_ds + beta * z_R + gamma * z_lapR
            eval_proxy("delta_ds_plus_beta_delta_R_plus_gamma_lapR", phi, beta=beta, gamma=gamma)

    node_df = pd.DataFrame({
        "node": np.arange(W0.shape[0]),
        "profile_f": profile,
        "delta_ds": delta_ds,
        "delta_R": delta_R,
        "lap_delta_R": lap_delta_R,
        "div_matter": div_matter,
        "matter_flux_norm": np.linalg.norm(f_matter, axis=1),
    })

    proxy_df = pd.DataFrame(proxy_rows)
    best = proxy_df.sort_values("rel_l2").iloc[0].to_dict()

    summary = {
        "sigma": sigma,
        "center": int(center),
        "best_proxy": best["proxy"],
        "best_beta": best.get("beta"),
        "best_gamma": best.get("gamma"),
        "best_alpha_opt": best["alpha_opt"],
        "best_rel_l2": best["rel_l2"],
        "best_rel_l1": best["rel_l1"],
        "best_improvement_factor_l2": best["improvement_factor_l2"],
        "best_verdict": best["verdict"],
        "matter_div_l2": best["matter_div_l2"],
        "best_total_div_l2": best["total_div_l2"],
    }

    return summary, proxy_df, node_df


# ============================================================
# Figures
# ============================================================

def make_figures(global_df, all_df, output_dir):
    plt.figure(figsize=(8, 5))
    plt.plot(global_df["sigma"], global_df["best_rel_l2"], "o-")
    plt.axhline(1.0, linewidth=1)
    plt.axhline(0.5, linestyle="--", linewidth=1)
    plt.axhline(0.1, linestyle="--", linewidth=1)
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"meilleur résidu relatif \(L^2\)")
    plt.title("Conservation totale v2 — meilleur proxy par sigma")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_best_rel_l2_vs_sigma.png"), dpi=250)
    plt.close()

    # Top 20 proxies overall
    top = all_df.sort_values("rel_l2").head(20)
    labels = []
    vals = []

    for _, r in top.iterrows():
        b = "" if pd.isna(r.get("beta")) else f",b={r.get('beta'):g}"
        g = "" if pd.isna(r.get("gamma")) else f",g={r.get('gamma'):g}"
        labels.append(f"s={r['sigma']:g},{r['proxy']}{b}{g}")
        vals.append(r["rel_l2"])

    plt.figure(figsize=(12, 5))
    plt.bar(range(len(vals)), vals)
    plt.xticks(range(len(vals)), labels, rotation=65, ha="right", fontsize=7)
    plt.ylabel(r"rel \(L^2\)")
    plt.title("Top 20 compensations totales")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_top20_total_conservation_v2.png"), dpi=250)
    plt.close()

    # Proxy comparison median/best
    best_by_proxy = all_df.groupby("proxy")["rel_l2"].min().sort_values()
    plt.figure(figsize=(9, 5))
    best_by_proxy.plot(kind="bar")
    plt.ylabel(r"meilleur rel \(L^2\)")
    plt.title("Meilleure compensation par type de proxy")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_best_by_proxy.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper8_total_conservation_v2")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--tau-points", type=int, default=20)

    parser.add_argument("--beta-list", type=float, nargs="+", default=[-3, -1, -0.5, 0, 0.5, 1, 3])
    parser.add_argument("--gamma-list", type=float, nargs="+", default=[-3, -1, -0.5, 0, 0.5, 1, 3])

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
    print("BuP Paper 8 — Test conservation totale v2")
    print("=" * 100)
    print(f"MI file     : {args.mi_file}")
    print(f"label       : {label}")
    print(f"lambda      : {lam}")
    print(f"N           : {W0.shape[0]}")
    print(f"k           : {args.k}")
    print(f"amp         : {args.amp}")
    print(f"sigma list  : {args.sigma_list}")
    print(f"beta list   : {args.beta_list}")
    print(f"gamma list  : {args.gamma_list}")
    print(f"tau window  : [{args.tau_min}, {args.tau_max}] with {args.tau_points} points")
    print(f"output      : {args.output_dir}")
    print("-" * 100)

    global_rows = []
    all_proxy_rows = []

    for sigma in args.sigma_list:
        summary, proxy_df, node_df = run_one(W0, args, sigma)
        global_rows.append(summary)

        proxy_df["sigma"] = sigma
        all_proxy_rows.append(proxy_df)

        subdir = os.path.join(args.output_dir, f"sigma_{sigma:.3g}".replace(".", "p"))
        ensure_dir(subdir)
        proxy_df.to_csv(os.path.join(subdir, "proxy_scan.csv"), index=False)
        node_df.to_csv(os.path.join(subdir, "node_fields.csv"), index=False)

        print(
            f"sigma={sigma:.3g} "
            f"best={summary['best_proxy']} "
            f"beta={summary['best_beta']} "
            f"gamma={summary['best_gamma']} "
            f"alpha={summary['best_alpha_opt']:+.4g} "
            f"rel_L2={summary['best_rel_l2']:.4g} "
            f"improv={summary['best_improvement_factor_l2']:.3g} "
            f"verdict={summary['best_verdict']}"
        )

    global_df = pd.DataFrame(global_rows)
    all_df = pd.concat(all_proxy_rows, ignore_index=True)

    global_csv = os.path.join(args.output_dir, "total_conservation_v2_sigma_summary.csv")
    all_csv = os.path.join(args.output_dir, "total_conservation_v2_all_proxy_scan.csv")

    global_df.to_csv(global_csv, index=False)
    all_df.to_csv(all_csv, index=False)

    best_overall = all_df.sort_values("rel_l2").iloc[0].to_dict()

    summary_global = {
        "title": "BuP Paper 8 — Total conservation test v2",
        "description": "Tests richer entanglement flux proxies built from delta d_s, delta R, and Lap(delta R).",
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
        "beta_list": args.beta_list,
        "gamma_list": args.gamma_list,
        "definition": {
            "matter_flux": "f_matter(i)=sum_j deltaW_ij^2 v_ij",
            "ent_flux": "f_ent(i)=-sum_j A_ij(phi_j-phi_i)v_ij",
            "proxies": [
                "phi=delta d_s",
                "phi=delta R",
                "phi=Lap(delta R)",
                "phi=delta d_s + beta delta R",
                "phi=delta R + gamma Lap(delta R)",
                "phi=delta d_s + beta delta R + gamma Lap(delta R)"
            ],
            "objective": "minimize ||div(f_matter + alpha f_ent)||_2"
        },
        "best_overall": best_overall,
        "files": {
            "sigma_summary": "total_conservation_v2_sigma_summary.csv",
            "all_proxy_scan": "total_conservation_v2_all_proxy_scan.csv",
            "fig_best_rel": "fig_best_rel_l2_vs_sigma.png",
            "fig_top20": "fig_top20_total_conservation_v2.png",
            "fig_best_by_proxy": "fig_best_by_proxy.png"
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary_global, f, indent=2)

    make_figures(global_df, all_df, args.output_dir)

    print("-" * 100)
    print("MEILLEUR CAS GLOBAL — conservation totale v2")
    print("-" * 100)
    print(f"sigma       : {best_overall['sigma']}")
    print(f"proxy       : {best_overall['proxy']}")
    print(f"beta        : {best_overall.get('beta')}")
    print(f"gamma       : {best_overall.get('gamma')}")
    print(f"alpha_opt   : {best_overall['alpha_opt']:+.6g}")
    print(f"rel_L2      : {best_overall['rel_l2']:.6g}")
    print(f"rel_L1      : {best_overall['rel_l1']:.6g}")
    print(f"improvement : {best_overall['improvement_factor_l2']:.6g}")
    print(f"verdict     : {best_overall['verdict']}")

    print("-" * 100)
    print("Fichiers générés")
    print("-" * 100)
    print(global_csv)
    print(all_csv)
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_best_rel_l2_vs_sigma.png"))
    print("DONE")


if __name__ == "__main__":
    main()

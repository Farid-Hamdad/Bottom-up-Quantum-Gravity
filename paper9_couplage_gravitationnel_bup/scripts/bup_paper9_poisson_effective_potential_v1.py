#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 9 — Poisson effective potential test.

Goal:
    Couple the Paper 8 emergent matter source

        S_flux = T00 - 1/2 Taa + 1/2 Tgrad + ||T0a||

    to a discrete Poisson equation on the entanglement graph:

        L_ent Phi = S_flux

    Then test whether Phi and |grad Phi| correlate with the geometric response |delta R|.

Main outputs:
    - paper9_poisson_summary.csv
    - node_table.csv
    - summary.json
    - figures

Usage:
    python3 papers/paper9_couplage_gravitationnel_bup/scripts/bup_paper9_poisson_effective_potential_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma 0.15 \
      --output-dir papers/paper9_couplage_gravitationnel_bup/results/run_N20_lam057_sigma015
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
from sklearn.manifold import MDS


# ============================================================
# Basic utilities
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


def norm01(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-15:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


# ============================================================
# Graph construction
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
        raise ValueError("No positive edges in kNN graph.")

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


def combinatorial_laplacian(A):
    deg = np.sum(A, axis=1)
    return np.diag(deg) - A


# ============================================================
# MDS embedding
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

def inject_excitation(W0, D0, amp=0.15, sigma=0.15, center=-1):
    W0 = clean_matrix(W0)

    if center is None or center < 0:
        center = choose_center(W0)

    d = D0[center].copy()
    finite = np.isfinite(d)

    if not np.any(finite):
        raise ValueError("Non-finite distances.")

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
# Paper 8 source components
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


def compute_flux_norm(deltaW, D, X, eps=1e-12):
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
    return flux_norm, flux_vec


# ============================================================
# Poisson solve
# ============================================================

def solve_poisson(L, source, ridge=1e-8):
    """
    Solve L Phi = source on a graph Laplacian.

    Since L has a zero mode, we remove the mean source and solve with
    a small ridge. Then we remove the mean of Phi.
    """
    source = np.asarray(source, dtype=float)
    source0 = source - np.mean(source)

    n = L.shape[0]
    A = L + ridge * np.eye(n)

    Phi = np.linalg.solve(A, source0)
    Phi = Phi - np.mean(Phi)

    return Phi


def gradient_norm_scalar(phi, A, D, X, eps=1e-12):
    n, ndim = X.shape
    grad_norm = np.zeros(n, dtype=float)

    for i in range(n):
        g = np.zeros(ndim, dtype=float)

        for j in range(n):
            wij = A[i, j]
            if wij <= 0 or i == j:
                continue

            dij = D[i, j]
            if not np.isfinite(dij) or dij <= eps:
                continue

            v = (X[j] - X[i]) / (dij + eps)
            g += wij * (phi[j] - phi[i]) * v

        grad_norm[i] = np.linalg.norm(g)

    return grad_norm


# ============================================================
# Radial profile
# ============================================================

def radial_profile(values, distances, n_bins=8):
    values = np.asarray(values, dtype=float)
    distances = np.asarray(distances, dtype=float)

    finite = np.isfinite(values) & np.isfinite(distances)
    values = values[finite]
    distances = distances[finite]

    if len(values) == 0:
        return pd.DataFrame(columns=["r_min", "r_max", "r_mid", "mean", "std", "count"])

    rmin, rmax = np.min(distances), np.max(distances)
    if abs(rmax - rmin) < 1e-15:
        bins = np.linspace(rmin, rmin + 1e-9, n_bins + 1)
    else:
        bins = np.linspace(rmin, rmax, n_bins + 1)

    rows = []
    for a, b in zip(bins[:-1], bins[1:]):
        mask = (distances >= a) & (distances < b)
        if b == bins[-1]:
            mask = (distances >= a) & (distances <= b)

        vals = values[mask]
        rows.append({
            "r_min": float(a),
            "r_max": float(b),
            "r_mid": float(0.5 * (a + b)),
            "mean": float(np.mean(vals)) if vals.size else np.nan,
            "std": float(np.std(vals)) if vals.size else np.nan,
            "count": int(vals.size)
        })

    return pd.DataFrame(rows)


# ============================================================
# Plotting
# ============================================================

def make_figures(node_df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    # Source vs |deltaR|
    plt.figure(figsize=(7, 5))
    plt.scatter(node_df["S_flux"], node_df["abs_delta_R"])
    plt.xlabel(r"$S_{\rm flux}$")
    plt.ylabel(r"$|\delta R|$")
    plt.title("Paper 9 — Source effective vs réponse de courbure")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_source_vs_curvature.png"), dpi=250)
    plt.close()

    # Potential vs |deltaR|
    plt.figure(figsize=(7, 5))
    plt.scatter(node_df["Phi"], node_df["abs_delta_R"])
    plt.xlabel(r"$\Phi_{\rm BuP}$")
    plt.ylabel(r"$|\delta R|$")
    plt.title("Paper 9 — Potentiel effectif vs réponse de courbure")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_potential_vs_curvature.png"), dpi=250)
    plt.close()

    # Acceleration vs |deltaR|
    plt.figure(figsize=(7, 5))
    plt.scatter(node_df["a_norm"], node_df["abs_delta_R"])
    plt.xlabel(r"$|\nabla \Phi_{\rm BuP}|$")
    plt.ylabel(r"$|\delta R|$")
    plt.title("Paper 9 — Accélération effective vs réponse de courbure")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_acceleration_vs_curvature.png"), dpi=250)
    plt.close()

    # Radial profiles
    for col, fname, ylabel in [
        ("S_flux", "fig_source_profile.png", r"$S_{\rm flux}(r)$"),
        ("Phi", "fig_phi_profile.png", r"$\Phi_{\rm BuP}(r)$"),
        ("a_norm", "fig_acceleration_profile.png", r"$|\nabla \Phi_{\rm BuP}|(r)$"),
        ("abs_delta_R", "fig_curvature_response_profile.png", r"$|\delta R|(r)$"),
    ]:
        prof = radial_profile(node_df[col].values, node_df["r_ent"].values, n_bins=8)
        plt.figure(figsize=(7, 5))
        plt.plot(prof["r_mid"], prof["mean"], "o-")
        plt.xlabel(r"$r_{\rm ent}$")
        plt.ylabel(ylabel)
        plt.title(f"Paper 9 — Profil radial : {col}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, fname), dpi=250)
        plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", required=True)
    parser.add_argument("--output-dir", default="results_paper9_poisson_effective_potential_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=0.15)
    parser.add_argument("--center", type=int, default=-1)

    parser.add_argument("--omega", type=float, default=-0.5)
    parser.add_argument("--chi", type=float, default=0.5)
    parser.add_argument("--psi", type=float, default=1.0)

    parser.add_argument("--normalize-components", action="store_true")
    parser.add_argument("--ridge", type=float, default=1e-8)

    parser.add_argument("--ndim", type=int, default=3)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    mask0 = knn_mask(W0, args.k)
    A0 = adjacency_from_mask(W0, mask0)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)
    X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)
    L = combinatorial_laplacian(A0)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=args.sigma,
        center=args.center
    )

    R0 = node_curvature_R(W0, args.k, eps=args.eps)
    R1 = node_curvature_R(W_exc, args.k, eps=args.eps)
    delta_R = R1 - R0
    abs_delta_R = np.abs(delta_R)

    T00 = compute_T00(deltaW)
    Taa = compute_Taa(deltaW, D0, X, eps=args.eps)
    Tgrad = compute_Tgrad(T00, A0)
    Fflux, flux_vec = compute_flux_norm(deltaW, D0, X, eps=args.eps)

    if args.normalize_components:
        T00_eff = norm01(T00)
        Taa_eff = norm01(Taa)
        Tgrad_eff = norm01(Tgrad)
        Fflux_eff = norm01(Fflux)
    else:
        T00_eff = T00
        Taa_eff = Taa
        Tgrad_eff = Tgrad
        Fflux_eff = Fflux

    S_flux = (
        T00_eff
        + args.omega * Taa_eff
        + args.chi * Tgrad_eff
        + args.psi * Fflux_eff
    )

    # Positive mass convention for Poisson.
    # We shift only if the source has negative values due to pressure subtraction.
    S_poisson = S_flux - np.min(S_flux)
    if np.max(S_poisson) > 1e-15:
        S_poisson = S_poisson / np.sum(S_poisson)
    else:
        S_poisson = S_flux - np.mean(S_flux)

    Phi = solve_poisson(L, S_poisson, ridge=args.ridge)
    a_norm = gradient_norm_scalar(Phi, A0, D0, X, eps=args.eps)

    r_ent = D0[center].copy()

    node_df = pd.DataFrame({
        "node": np.arange(W0.shape[0]),
        "r_ent": r_ent,
        "profile_f": profile,
        "T00": T00,
        "Taa": Taa,
        "Tgrad": Tgrad,
        "Fflux": Fflux,
        "S_flux": S_flux,
        "S_poisson": S_poisson,
        "Phi": Phi,
        "a_norm": a_norm,
        "R0": R0,
        "R1": R1,
        "delta_R": delta_R,
        "abs_delta_R": abs_delta_R,
    })

    node_path = os.path.join(args.output_dir, "node_table.csv")
    node_df.to_csv(node_path, index=False)

    # Correlations
    metrics = {}

    pairs = [
        ("S_flux", S_flux, abs_delta_R),
        ("S_poisson", S_poisson, abs_delta_R),
        ("Phi", Phi, abs_delta_R),
        ("abs_Phi", np.abs(Phi), abs_delta_R),
        ("a_norm", a_norm, abs_delta_R),
        ("r_ent_vs_Phi", r_ent, Phi),
        ("r_ent_vs_a_norm", r_ent, a_norm),
    ]

    for name, x, y in pairs:
        pr, pp, sr, sp = safe_corr(x, y)
        metrics[f"{name}_pearson"] = pr
        metrics[f"{name}_pearson_p"] = pp
        metrics[f"{name}_spearman"] = sr
        metrics[f"{name}_spearman_p"] = sp

    summary = {
        "title": "BuP Paper 9 — Poisson effective potential test",
        "description": "Couples the Paper 8 S_flux source to a discrete Poisson equation L_ent Phi = S_flux.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "sigma": args.sigma,
        "center": int(center),
        "omega": args.omega,
        "chi": args.chi,
        "psi": args.psi,
        "normalize_components": args.normalize_components,
        "ridge": args.ridge,
        "source": {
            "definition": "S_flux = T00 + omega*Taa + chi*Tgrad + psi*||T0a||",
            "omega": args.omega,
            "chi": args.chi,
            "psi": args.psi,
            "S_flux_min": float(np.min(S_flux)),
            "S_flux_max": float(np.max(S_flux)),
            "S_flux_sum": float(np.sum(S_flux)),
            "S_poisson_sum": float(np.sum(S_poisson)),
        },
        "poisson": {
            "equation": "L_ent Phi = S_poisson",
            "Phi_min": float(np.min(Phi)),
            "Phi_max": float(np.max(Phi)),
            "Phi_mean": float(np.mean(Phi)),
            "a_norm_min": float(np.min(a_norm)),
            "a_norm_max": float(np.max(a_norm)),
            "a_norm_mean": float(np.mean(a_norm)),
        },
        "correlations": metrics,
        "files": {
            "node_table": "node_table.csv",
            "summary_csv": "paper9_poisson_summary.csv",
            "figures_dir": "figures/"
        }
    }

    summary_path = os.path.join(args.output_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    row = {
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "sigma": args.sigma,
        "center": int(center),
        "omega": args.omega,
        "chi": args.chi,
        "psi": args.psi,
        "normalize_components": args.normalize_components,
        "S_flux_spearman_abs_deltaR": metrics["S_flux_spearman"],
        "S_flux_spearman_p": metrics["S_flux_spearman_p"],
        "Phi_spearman_abs_deltaR": metrics["Phi_spearman"],
        "Phi_spearman_p": metrics["Phi_spearman_p"],
        "absPhi_spearman_abs_deltaR": metrics["abs_Phi_spearman"],
        "absPhi_spearman_p": metrics["abs_Phi_spearman_p"],
        "a_norm_spearman_abs_deltaR": metrics["a_norm_spearman"],
        "a_norm_spearman_p": metrics["a_norm_spearman_p"],
    }

    summary_csv = os.path.join(args.output_dir, "paper9_poisson_summary.csv")
    pd.DataFrame([row]).to_csv(summary_csv, index=False)

    make_figures(node_df, args.output_dir)

    print("=" * 100)
    print("BuP Paper 9 — Poisson effective potential test")
    print("=" * 100)
    print(f"MI file   : {args.mi_file}")
    print(f"label     : {label}")
    print(f"lambda    : {lam}")
    print(f"N         : {W0.shape[0]}")
    print(f"k         : {args.k}")
    print(f"sigma     : {args.sigma}")
    print(f"center    : {center}")
    print(f"output    : {args.output_dir}")
    print("-" * 100)
    print("Corrélations avec |delta R|")
    print("-" * 100)
    print(f"S_flux     Spearman = {metrics['S_flux_spearman']:+.6f}, p={metrics['S_flux_spearman_p']:.6g}")
    print(f"Phi        Spearman = {metrics['Phi_spearman']:+.6f}, p={metrics['Phi_spearman_p']:.6g}")
    print(f"|Phi|      Spearman = {metrics['abs_Phi_spearman']:+.6f}, p={metrics['abs_Phi_spearman_p']:.6g}")
    print(f"|grad Phi| Spearman = {metrics['a_norm_spearman']:+.6f}, p={metrics['a_norm_spearman_p']:.6g}")
    print("-" * 100)
    print("Fichiers")
    print(node_path)
    print(summary_csv)
    print(summary_path)
    print("DONE")


if __name__ == "__main__":
    main()

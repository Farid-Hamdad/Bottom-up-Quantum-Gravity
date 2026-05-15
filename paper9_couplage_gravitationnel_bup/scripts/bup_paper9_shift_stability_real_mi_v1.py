#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 9 — Shift stability test on real MI matrix + Ollivier-Ricci curvature.

Goal:
    Test whether shifting the Paper 8 source

        S_shift = S_flux - min(S_flux)

    changes the Poisson potential:

        L_ent Phi = S

    when the input is a real MI matrix and the geometric response is computed
    with Ollivier-Ricci curvature.

Tested variants:
    S_raw   = S_flux
    S_shift = S_flux - min(S_flux)
    S_norm  = S_shift / sum(S_shift)

For each sigma:
    - build S_flux
    - compute delta R with Ollivier-Ricci curvature
    - solve Poisson for raw/shift/norm
    - compare rho(Phi, |delta R|)
    - measure Delta rho

Expected:
    Since L_ent 1 = 0, adding a constant to the source should not change
    the physical potential after projection/recentering.

Usage:
    python3 papers/paper9_couplage_gravitationnel_bup/scripts/bup_paper9_shift_stability_real_mi_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma-list 0.01 0.02 0.05 0.08 0.10 0.15 \
      --omega -0.5 \
      --chi 0.5 \
      --psi 1.0 \
      --normalize-components \
      --output-dir papers/paper9_couplage_gravitationnel_bup/results/shift_stability_real_mi_N20_lam057
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


def norm01(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-15:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    if len(x) < 4:
        return np.nan, np.nan, np.nan, np.nan

    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return np.nan, np.nan, np.nan, np.nan

    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    return float(pr), float(pp), float(sr), float(sp)


def centered(x):
    x = np.asarray(x, dtype=float)
    return x - np.mean(x)


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
# Poisson solve by pseudo-inverse
# ============================================================

def solve_poisson_pinv(L, source, eig_tol=1e-10):
    """
    Solve L Phi = source using spectral pseudo-inverse.

    The constant mode is projected out:
        source -> source - mean(source)

    Therefore adding c*1 to source is physically invisible.
    """
    source = np.asarray(source, dtype=float)
    b = source - np.mean(source)

    evals, evecs = np.linalg.eigh(L)
    inv = np.zeros_like(evals)

    good = evals > eig_tol
    inv[good] = 1.0 / evals[good]

    Phi = evecs @ (inv * (evecs.T @ b))
    Phi = Phi - np.mean(Phi)

    return Phi


# ============================================================
# One sigma test
# ============================================================

def run_one_sigma(W0, args, sigma):
    mask0 = knn_mask(W0, args.k)
    A0 = adjacency_from_mask(W0, mask0)
    D0 = entanglement_distance(W0, mask0, eps=args.eps)
    X = mds_coordinates(D0, ndim=args.ndim, random_state=args.random_state)
    L = combinatorial_laplacian(A0)

    W_exc, deltaW, center, profile = inject_excitation(
        W0,
        D0,
        amp=args.amp,
        sigma=sigma,
        center=args.center
    )

    # Ollivier-Ricci response
    R0 = node_curvature_R(W0, args.k, eps=args.eps)
    R1 = node_curvature_R(W_exc, args.k, eps=args.eps)
    delta_R = R1 - R0
    abs_delta_R = np.abs(delta_R)

    # Paper 8 source
    T00 = compute_T00(deltaW)
    Taa = compute_Taa(deltaW, D0, X, eps=args.eps)
    Tgrad = compute_Tgrad(T00, A0)
    Fflux, _ = compute_flux_norm(deltaW, D0, X, eps=args.eps)

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

    S_raw = (
        T00_eff
        + args.omega * Taa_eff
        + args.chi * Tgrad_eff
        + args.psi * Fflux_eff
    )

    S_shift = S_raw - np.min(S_raw)

    if np.sum(S_shift) > args.eps:
        S_norm = S_shift / np.sum(S_shift)
    else:
        S_norm = S_shift.copy()

    # Poisson potentials
    Phi_raw = solve_poisson_pinv(L, S_raw, eig_tol=args.eig_tol)
    Phi_shift = solve_poisson_pinv(L, S_shift, eig_tol=args.eig_tol)
    Phi_norm = solve_poisson_pinv(L, S_norm, eig_tol=args.eig_tol)

    # Correlations source vs curvature
    _, _, rho_Sraw, p_Sraw = safe_corr(S_raw, abs_delta_R)
    _, _, rho_Sshift, p_Sshift = safe_corr(S_shift, abs_delta_R)
    _, _, rho_Snorm, p_Snorm = safe_corr(S_norm, abs_delta_R)

    # Correlations potential vs curvature
    _, _, rho_Phi_raw, p_Phi_raw = safe_corr(Phi_raw, abs_delta_R)
    _, _, rho_Phi_shift, p_Phi_shift = safe_corr(Phi_shift, abs_delta_R)
    _, _, rho_Phi_norm, p_Phi_norm = safe_corr(Phi_norm, abs_delta_R)

    # Potential equivalence tests
    # Shift should be identical up to numerical precision.
    max_abs_phi_raw_shift = float(np.max(np.abs(Phi_raw - Phi_shift)))

    # Norm is a positive rescaling of shift, so Spearman should be identical.
    # For direct vector comparison, rescale Phi_norm to best match Phi_raw.
    denom = float(np.dot(Phi_norm, Phi_norm))
    if denom > args.eps:
        scale_norm_to_raw = float(np.dot(Phi_raw, Phi_norm) / denom)
        Phi_norm_rescaled = scale_norm_to_raw * Phi_norm
        max_abs_phi_raw_norm_rescaled = float(np.max(np.abs(Phi_raw - Phi_norm_rescaled)))
    else:
        scale_norm_to_raw = np.nan
        max_abs_phi_raw_norm_rescaled = np.nan

    row = {
        "sigma": sigma,
        "center": int(center),

        "S_raw_min": float(np.min(S_raw)),
        "S_raw_max": float(np.max(S_raw)),
        "S_raw_sum": float(np.sum(S_raw)),
        "S_shift_sum": float(np.sum(S_shift)),
        "S_norm_sum": float(np.sum(S_norm)),

        "rho_Sraw_abs_deltaR": rho_Sraw,
        "p_Sraw_abs_deltaR": p_Sraw,
        "rho_Sshift_abs_deltaR": rho_Sshift,
        "p_Sshift_abs_deltaR": p_Sshift,
        "rho_Snorm_abs_deltaR": rho_Snorm,
        "p_Snorm_abs_deltaR": p_Snorm,

        "rho_Phi_raw_abs_deltaR": rho_Phi_raw,
        "p_Phi_raw_abs_deltaR": p_Phi_raw,
        "rho_Phi_shift_abs_deltaR": rho_Phi_shift,
        "p_Phi_shift_abs_deltaR": p_Phi_shift,
        "rho_Phi_norm_abs_deltaR": rho_Phi_norm,
        "p_Phi_norm_abs_deltaR": p_Phi_norm,

        "delta_rho_shift": abs(rho_Phi_shift - rho_Phi_raw),
        "delta_rho_norm": abs(rho_Phi_norm - rho_Phi_raw),

        "max_abs_phi_raw_minus_shift": max_abs_phi_raw_shift,
        "scale_norm_to_raw": scale_norm_to_raw,
        "max_abs_phi_raw_minus_norm_rescaled": max_abs_phi_raw_norm_rescaled,

        "shift_stable_005": bool(abs(rho_Phi_shift - rho_Phi_raw) < 0.05),
        "shift_stable_010": bool(abs(rho_Phi_shift - rho_Phi_raw) < 0.10),
    }

    node_df = pd.DataFrame({
        "node": np.arange(W0.shape[0]),
        "profile_f": profile,
        "R0": R0,
        "R1": R1,
        "delta_R": delta_R,
        "abs_delta_R": abs_delta_R,
        "T00": T00,
        "Taa": Taa,
        "Tgrad": Tgrad,
        "Fflux": Fflux,
        "S_raw": S_raw,
        "S_shift": S_shift,
        "S_norm": S_norm,
        "Phi_raw": Phi_raw,
        "Phi_shift": Phi_shift,
        "Phi_norm": Phi_norm,
    })

    return row, node_df


# ============================================================
# Figures
# ============================================================

def make_figures(df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    # Figure 1: correlations vs sigma
    plt.figure(figsize=(9, 5))
    plt.plot(df["sigma"], df["rho_Sraw_abs_deltaR"], "o-", label=r"$\rho(S_{\rm raw},|\delta R|)$")
    plt.plot(df["sigma"], df["rho_Phi_raw_abs_deltaR"], "s-", label=r"$\rho(\Phi_{\rm raw},|\delta R|)$")
    plt.plot(df["sigma"], df["rho_Phi_shift_abs_deltaR"], "D-", label=r"$\rho(\Phi_{\rm shift},|\delta R|)$")
    plt.plot(df["sigma"], df["rho_Phi_norm_abs_deltaR"], "^-", label=r"$\rho(\Phi_{\rm norm},|\delta R|)$")
    plt.axhline(0.0, linewidth=1)
    plt.axvspan(0.05, 0.10, alpha=0.15, label="zone morte")
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$\rho$ Spearman")
    plt.title("BuP Paper 9 — Stabilité du shift sur MI réelle + Ollivier-Ricci")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_shift_stability_real_mi_correlations.png"), dpi=250)
    plt.close()

    # Figure 2: impact of shift
    plt.figure(figsize=(9, 5))
    plt.plot(df["sigma"], df["delta_rho_shift"], "o-", label=r"$|\Delta\rho_{\rm shift}|$")
    plt.plot(df["sigma"], df["delta_rho_norm"], "s-", label=r"$|\Delta\rho_{\rm norm}|$")
    plt.axhline(0.05, linestyle="--", linewidth=1, label="seuil 0.05")
    plt.axhline(0.10, linestyle="--", linewidth=1, label="seuil 0.10")
    plt.axvspan(0.05, 0.10, alpha=0.15, label="zone morte")
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$|\Delta \rho|$")
    plt.title("Impact du shift / normalisation sur le potentiel")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_shift_stability_real_mi_delta_rho.png"), dpi=250)
    plt.close()

    # Figure 3: raw vs shifted potential numerical difference
    plt.figure(figsize=(9, 5))
    plt.semilogy(df["sigma"], df["max_abs_phi_raw_minus_shift"], "o-", label=r"$\max|\Phi_{\rm raw}-\Phi_{\rm shift}|$")
    plt.semilogy(df["sigma"], df["max_abs_phi_raw_minus_norm_rescaled"], "s-", label=r"$\max|\Phi_{\rm raw}-c\Phi_{\rm norm}|$")
    plt.xlabel(r"$\sigma$")
    plt.ylabel("écart potentiel")
    plt.title("Équivalence numérique des potentiels")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_shift_stability_real_mi_phi_difference.png"), dpi=250)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", required=True)
    parser.add_argument("--output-dir", default="results_paper9_shift_stability_real_mi_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma-list", type=float, nargs="+", required=True)
    parser.add_argument("--center", type=int, default=-1)

    parser.add_argument("--omega", type=float, default=-0.5)
    parser.add_argument("--chi", type=float, default=0.5)
    parser.add_argument("--psi", type=float, default=1.0)

    parser.add_argument("--normalize-components", action="store_true")

    parser.add_argument("--ndim", type=int, default=3)
    parser.add_argument("--random-state", type=int, default=0)

    parser.add_argument("--eps", type=float, default=1e-12)
    parser.add_argument("--eig-tol", type=float, default=1e-10)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    W0 = load_matrix_csv(args.mi_file)
    label = parse_label(args.mi_file)
    lam = parse_lambda(args.mi_file)

    print("=" * 100)
    print("BuP Paper 9 — Shift stability test on real MI + Ollivier-Ricci")
    print("=" * 100)
    print(f"MI file   : {args.mi_file}")
    print(f"label     : {label}")
    print(f"lambda    : {lam}")
    print(f"N         : {W0.shape[0]}")
    print(f"k         : {args.k}")
    print(f"sigmas    : {args.sigma_list}")
    print(f"output    : {args.output_dir}")
    print("-" * 100)

    rows = []

    for sigma in args.sigma_list:
        row, node_df = run_one_sigma(W0, args, sigma)
        rows.append(row)

        subdir = os.path.join(args.output_dir, f"sigma_{sigma:.3g}".replace(".", "p"))
        ensure_dir(subdir)
        node_df.to_csv(os.path.join(subdir, "node_table_shift_stability.csv"), index=False)

        print(
            f"sigma={sigma:g} "
            f"rho_Sraw={row['rho_Sraw_abs_deltaR']:+.6f} "
            f"rho_Phi_raw={row['rho_Phi_raw_abs_deltaR']:+.6f} "
            f"rho_Phi_shift={row['rho_Phi_shift_abs_deltaR']:+.6f} "
            f"rho_Phi_norm={row['rho_Phi_norm_abs_deltaR']:+.6f} "
            f"d_shift={row['delta_rho_shift']:.3g} "
            f"d_norm={row['delta_rho_norm']:.3g} "
            f"stable={row['shift_stable_005']}"
        )

    df = pd.DataFrame(rows)
    csv_path = os.path.join(args.output_dir, "shift_stability_real_mi_summary.csv")
    df.to_csv(csv_path, index=False)

    make_figures(df, args.output_dir)

    best = df.sort_values("delta_rho_shift", ascending=False).iloc[0].to_dict()
    worst_delta_shift = float(best["delta_rho_shift"])
    worst_delta_norm = float(df["delta_rho_norm"].max())

    summary = {
        "title": "BuP Paper 9 — Shift stability test on real MI + Ollivier-Ricci",
        "description": "Tests whether S_flux shift and normalization affect the Poisson potential correlations.",
        "mi_file": args.mi_file,
        "label": label,
        "lambda": lam,
        "N": int(W0.shape[0]),
        "k": args.k,
        "amp": args.amp,
        "sigma_list": args.sigma_list,
        "omega": args.omega,
        "chi": args.chi,
        "psi": args.psi,
        "normalize_components": args.normalize_components,
        "definition": {
            "S_raw": "S_flux",
            "S_shift": "S_flux - min(S_flux)",
            "S_norm": "S_shift / sum(S_shift)",
            "poisson": "Phi = L_ent^+ S",
            "curvature": "|delta R| from Ollivier-Ricci node curvature"
        },
        "results": {
            "max_delta_rho_shift": worst_delta_shift,
            "max_delta_rho_norm": worst_delta_norm,
            "shift_stable_all_005": bool((df["delta_rho_shift"] < 0.05).all()),
            "shift_stable_all_010": bool((df["delta_rho_shift"] < 0.10).all()),
            "norm_stable_all_005": bool((df["delta_rho_norm"] < 0.05).all()),
            "norm_stable_all_010": bool((df["delta_rho_norm"] < 0.10).all())
        },
        "files": {
            "summary_csv": "shift_stability_real_mi_summary.csv",
            "fig_correlations": "figures/fig_shift_stability_real_mi_correlations.png",
            "fig_delta_rho": "figures/fig_shift_stability_real_mi_delta_rho.png",
            "fig_phi_difference": "figures/fig_shift_stability_real_mi_phi_difference.png"
        }
    }

    summary_path = os.path.join(args.output_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("-" * 100)
    print("RÉSUMÉ")
    print("-" * 100)
    print(f"max Delta rho shift : {worst_delta_shift:.12g}")
    print(f"max Delta rho norm  : {worst_delta_norm:.12g}")
    print(f"shift stable <0.05  : {summary['results']['shift_stable_all_005']}")
    print(f"norm stable <0.05   : {summary['results']['norm_stable_all_005']}")
    print("-" * 100)
    print("Fichiers")
    print(csv_path)
    print(summary_path)
    print(os.path.join(args.output_dir, "figures"))
    print("DONE")


if __name__ == "__main__":
    main()

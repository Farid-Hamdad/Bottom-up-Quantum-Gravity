#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_pure_matter_excitation_v1.py

BuP Paper 6 — Premier test numérique :
matière comme excitation localisée du réseau d'intrication.

Objectif :
    Partir d'une matrice MI de fond W0.
    Injecter une excitation localisée deltaW.
    Mesurer la réponse géométrique :
        delta d_s
        delta kappa
        delta G_proxy
    Comparer cette réponse avec un proxy de matière :
        T_matter_proxy = (deltaW)^2

Test principal :
    Corr( T_ij^matter_proxy , |delta G_ij^proxy| )
    Corr( rho_i^exc , |delta d_s(i)| )

Usage :
    python3 bup_pure_matter_excitation_v1.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --amp 0.15 \
      --sigma 2.0 \
      --center 0 \
      --tau-min 0.01 \
      --tau-max 50 \
      --output-dir results_pure_matter_excitation_v1

Batch simple :
    for f in results_mi_N20_full/*.csv; do
      python3 bup_pure_matter_excitation_v1.py \
        --mi-file "$f" \
        --k 5 \
        --amp 0.15 \
        --sigma 2.0 \
        --output-dir results_pure_matter_excitation_v1/$(basename "$f" .csv)
    done
"""

import os
import re
import json
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import expm
from scipy.stats import pearsonr, spearmanr


# ============================================================
# UTILITAIRES
# ============================================================

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def clean_matrix(W: np.ndarray) -> np.ndarray:
    W = np.asarray(W, dtype=float)

    if W.shape[0] != W.shape[1]:
        raise ValueError(f"Matrice non carrée : shape={W.shape}")

    W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    W[W < 0] = 0.0
    return W


def load_matrix_csv(path: str) -> np.ndarray:
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


def parse_label(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def infer_N_lambda(label: str):
    N = None
    lam = None

    m = re.search(r"N0?(\d+)", label)
    if m:
        N = int(m.group(1))

    m = re.search(r"lam(?:bda)?([0-9]+(?:\.[0-9]+)?)", label)
    if m:
        lam = float(m.group(1))

    return N, lam


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    out = {
        "n": int(x.size),
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "spearman_r": np.nan,
        "spearman_p": np.nan,
    }

    if x.size < 4:
        return out

    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return out

    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    out.update({
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "spearman_r": float(sr),
        "spearman_p": float(sp),
    })
    return out


# ============================================================
# GRAPHE
# ============================================================

def knn_mask(W: np.ndarray, k: int) -> np.ndarray:
    n = W.shape[0]
    mask = np.zeros_like(W, dtype=bool)

    for i in range(n):
        idx = np.argsort(W[i])[::-1]
        idx = [j for j in idx if j != i and W[i, j] > 0]
        for j in idx[:k]:
            mask[i, j] = True

    mask = np.logical_or(mask, mask.T)
    np.fill_diagonal(mask, False)
    return mask


def weighted_adjacency_from_mask(W: np.ndarray, mask: np.ndarray) -> np.ndarray:
    A = np.zeros_like(W)
    A[mask] = W[mask]
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def build_edge_lengths(W: np.ndarray, mask: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    Lell = np.full_like(W, np.inf, dtype=float)
    vals = W[mask]

    if vals.size == 0 or np.max(vals) <= 0:
        raise ValueError("Aucune arête positive dans le graphe.")

    Wmax = np.max(vals)
    Lell[mask] = -np.log((W[mask] + eps) / (Wmax + eps))

    good = np.isfinite(Lell)
    Lell[good] = np.maximum(Lell[good], 1e-9)

    np.fill_diagonal(Lell, 0.0)
    return Lell


def shortest_paths_floyd(lengths: np.ndarray) -> np.ndarray:
    D = lengths.copy()
    n = D.shape[0]

    for i in range(n):
        D[i, i] = 0.0

    for k in range(n):
        D = np.minimum(D, D[:, [k]] + D[[k], :])

    return D


# ============================================================
# DIMENSION SPECTRALE
# ============================================================

def graph_laplacian(A: np.ndarray, normalized: bool = True) -> np.ndarray:
    deg = np.sum(A, axis=1)

    if normalized:
        invsqrt = np.zeros_like(deg)
        good = deg > 1e-15
        invsqrt[good] = 1.0 / np.sqrt(deg[good])
        Dm = np.diag(invsqrt)
        L = np.eye(A.shape[0]) - Dm @ A @ Dm

        for i, g in enumerate(good):
            if not g:
                L[i, i] = 0.0

        return L

    return np.diag(deg) - A


def local_spectral_dimension(A: np.ndarray, tau_min: float, tau_max: float, n_tau: int = 24):
    n = A.shape[0]
    L = graph_laplacian(A, normalized=True)

    taus = np.geomspace(tau_min, tau_max, n_tau)
    logt = np.log(taus)

    P = np.zeros((n_tau, n))

    for a, tau in enumerate(taus):
        K = expm(-tau * L)
        P[a, :] = np.clip(np.diag(K), 1e-300, None)

    ds = np.zeros(n)
    r2 = np.zeros(n)

    for i in range(n):
        y = np.log(P[:, i])
        coeff = np.polyfit(logt, y, 1)
        slope = coeff[0]
        fit = coeff[0] * logt + coeff[1]

        ss_res = np.sum((y - fit) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)

        r2[i] = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else np.nan
        ds[i] = -2.0 * slope

    return ds, r2


# ============================================================
# COURBURE OLLIVIER-RICCI PROXY
# ============================================================

def neighborhood_measure(A: np.ndarray, i: int):
    idx = np.where(A[i] > 0)[0]

    if idx.size == 0:
        return np.array([i], dtype=int), np.array([1.0])

    weights = A[i, idx].astype(float)
    s = np.sum(weights)

    if s <= 0:
        return np.array([i], dtype=int), np.array([1.0])

    return idx, weights / s


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


def ollivier_ricci_proxy(A: np.ndarray, lengths: np.ndarray, mask: np.ndarray):
    Dsp = shortest_paths_floyd(lengths)
    n = A.shape[0]
    kappa = {}

    measures = [neighborhood_measure(A, i) for i in range(n)]

    edges = np.argwhere(np.triu(mask, 1))
    for i, j in edges:
        dij = Dsp[i, j]

        if not np.isfinite(dij) or dij <= 1e-12:
            continue

        idx_i, mass_i = measures[i]
        idx_j, mass_j = measures[j]

        w1 = wasserstein_greedy(idx_i, mass_i, idx_j, mass_j, Dsp)
        k = 1.0 - w1 / dij
        kappa[(int(i), int(j))] = float(k)

    return kappa, Dsp


def node_scalar_curvature(n: int, kappa: dict):
    vals = [[] for _ in range(n)]

    for (i, j), k in kappa.items():
        vals[i].append(k)
        vals[j].append(k)

    R = np.zeros(n)
    for i in range(n):
        R[i] = float(np.mean(vals[i])) if vals[i] else np.nan

    return R


def edge_geometry_table(label, W, A, lengths, mask, ds):
    """
    Calcule kappa, R_edge, G_proxy sur les arêtes.

    Convention Paper 5 stable :
        G_proxy = -kappa + 1/2 R_edge
    """
    n = W.shape[0]
    kappa, Dsp = ollivier_ricci_proxy(A, lengths, mask)
    R = node_scalar_curvature(n, kappa)

    rows = []

    for (i, j), kij in kappa.items():
        Redge = np.nanmean([R[i], R[j]])
        Gproxy = -kij + 0.5 * Redge

        rows.append({
            "label": label,
            "i": i,
            "j": j,
            "W_ij": W[i, j],
            "ell_ij": lengths[i, j],
            "kappa_ij": kij,
            "R_i": R[i],
            "R_j": R[j],
            "R_edge": Redge,
            "G_proxy": Gproxy,
            "ds_i": ds[i],
            "ds_j": ds[j],
        })

    return pd.DataFrame(rows), R


# ============================================================
# EXCITATION LOCALISÉE
# ============================================================

def choose_center(W0, center_arg):
    n = W0.shape[0]
    strength = np.sum(W0, axis=1)

    if center_arg is not None and center_arg >= 0:
        if center_arg >= n:
            raise ValueError(f"center={center_arg} hors range N={n}")
        return int(center_arg)

    # par défaut : nœud le plus central au sens strength
    return int(np.argmax(strength))


def localized_profile_from_distances(Dsp, center, sigma):
    d = Dsp[center].copy()

    finite = np.isfinite(d)
    if not np.any(finite):
        raise RuntimeError("Distances infinies : graphe non connecté ?")

    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2.0 * sigma ** 2))

    # normalisation max=1
    if np.max(f) > 0:
        f = f / np.max(f)

    return f


def inject_local_excitation(W0, f, amp, mask_mode="all", base_mask=None):
    """
    Perturbation :
        deltaW_ij = amp * f_i f_j * scale

    scale = mean positive W0 pour rester dans l'ordre de grandeur.
    """
    W0 = clean_matrix(W0)
    n = W0.shape[0]

    positive = W0[W0 > 0]
    scale = float(np.mean(positive)) if positive.size else 1.0

    delta = amp * scale * np.outer(f, f)
    np.fill_diagonal(delta, 0.0)

    if mask_mode == "edges":
        if base_mask is None:
            raise ValueError("base_mask requis si mask_mode='edges'")
        delta = delta * base_mask

    W1 = W0 + delta
    W1 = clean_matrix(W1)

    return W1, delta


def node_excitation_density(deltaW):
    return np.sum(np.abs(deltaW), axis=1)


def edge_excitation_table(edge_df0, edge_df1, deltaW):
    """
    Fusionne géométrie avant/après sur les arêtes communes.
    """
    key = ["i", "j"]

    a = edge_df0.copy()
    b = edge_df1.copy()

    a = a.rename(columns={
        "W_ij": "W0_ij",
        "kappa_ij": "kappa0_ij",
        "R_edge": "R0_edge",
        "G_proxy": "G0_proxy",
        "ds_i": "ds0_i",
        "ds_j": "ds0_j",
    })

    b = b.rename(columns={
        "W_ij": "W1_ij",
        "kappa_ij": "kappa1_ij",
        "R_edge": "R1_edge",
        "G_proxy": "G1_proxy",
        "ds_i": "ds1_i",
        "ds_j": "ds1_j",
    })

    merged = pd.merge(
        a[["i", "j", "W0_ij", "kappa0_ij", "R0_edge", "G0_proxy", "ds0_i", "ds0_j"]],
        b[["i", "j", "W1_ij", "kappa1_ij", "R1_edge", "G1_proxy", "ds1_i", "ds1_j"]],
        on=key,
        how="inner"
    )

    vals = []
    for _, row in merged.iterrows():
        i = int(row["i"])
        j = int(row["j"])

        dWij = deltaW[i, j]
        Tmatter = dWij ** 2

        deltaG = row["G1_proxy"] - row["G0_proxy"]
        deltaKappa = row["kappa1_ij"] - row["kappa0_ij"]
        deltaR = row["R1_edge"] - row["R0_edge"]

        vals.append({
            **row.to_dict(),
            "deltaW_ij": dWij,
            "T_matter_edge": Tmatter,
            "delta_G_proxy": deltaG,
            "abs_delta_G_proxy": abs(deltaG),
            "delta_kappa": deltaKappa,
            "abs_delta_kappa": abs(deltaKappa),
            "delta_R_edge": deltaR,
            "abs_delta_R_edge": abs(deltaR),
        })

    return pd.DataFrame(vals)


# ============================================================
# FIGURES
# ============================================================

def make_figures(node_df, edge_df, output_dir):
    ensure_dir(output_dir)

    # Figure 1 : densité d'excitation vs |delta ds|
    x = node_df["rho_exc"].values
    y = node_df["abs_delta_ds"].values
    good = np.isfinite(x) & np.isfinite(y)

    plt.figure(figsize=(7, 5))
    plt.scatter(x[good], y[good], s=35, alpha=0.75)
    plt.xlabel(r"$\rho_i^{\rm exc}$")
    plt.ylabel(r"$|\delta d_s(i)|$")
    plt.title("Excitation locale vs réponse de dimension spectrale")
    plt.grid(True, alpha=0.3)

    if np.sum(good) >= 4 and np.std(x[good]) > 1e-15:
        coeff = np.polyfit(x[good], y[good], 1)
        xx = np.linspace(np.min(x[good]), np.max(x[good]), 200)
        yy = coeff[0] * xx + coeff[1]
        plt.plot(xx, yy, linewidth=2)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_node_rho_vs_delta_ds.png"), dpi=220)
    plt.close()

    # Figure 2 : T_matter_edge vs |delta G|
    x = edge_df["T_matter_edge"].values
    y = edge_df["abs_delta_G_proxy"].values
    good = np.isfinite(x) & np.isfinite(y)

    plt.figure(figsize=(7, 5))
    plt.scatter(x[good], y[good], s=22, alpha=0.75)
    plt.xlabel(r"$T_{ij}^{\rm matter,proxy}=(\delta W_{ij})^2$")
    plt.ylabel(r"$|\delta G_{ij}^{\rm proxy}|$")
    plt.title("Proxy matière vs réponse de courbure")
    plt.grid(True, alpha=0.3)

    if np.sum(good) >= 4 and np.std(x[good]) > 1e-15:
        coeff = np.polyfit(x[good], y[good], 1)
        xx = np.linspace(np.min(x[good]), np.max(x[good]), 200)
        yy = coeff[0] * xx + coeff[1]
        plt.plot(xx, yy, linewidth=2)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_edge_Tmatter_vs_deltaG.png"), dpi=220)
    plt.close()

    # Figure 3 : profil par noeud
    plt.figure(figsize=(8, 5))
    idx = np.arange(len(node_df))
    plt.plot(idx, node_df["rho_exc"].values, marker="o", label=r"$\rho_i^{exc}$")
    plt.plot(idx, node_df["abs_delta_ds"].values, marker="s", label=r"$|\delta d_s|$")
    plt.xlabel("node")
    plt.ylabel("value")
    plt.title("Profil nodal excitation / réponse")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_node_profiles.png"), dpi=220)
    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_pure_matter_excitation_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=2.0)
    parser.add_argument("--center", type=int, default=-1,
                        help="Nœud centre de l'excitation. -1 = nœud de strength max.")

    parser.add_argument("--mask-mode", choices=["all", "edges"], default="all",
                        help="all = deltaW sur toutes paires ; edges = seulement arêtes kNN de fond.")

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--n-tau", type=int, default=24)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    label = parse_label(args.mi_file)
    N_guess, lam_guess = infer_N_lambda(label)

    print("=" * 88)
    print("BuP Paper 6 — v1 matière comme excitation localisée du réseau")
    print("=" * 88)
    print(f"MI file    : {args.mi_file}")
    print(f"label      : {label}")
    print(f"k          : {args.k}")
    print(f"amp        : {args.amp}")
    print(f"sigma      : {args.sigma}")
    print(f"center     : {args.center}")
    print(f"mask mode  : {args.mask_mode}")
    print(f"tau        : [{args.tau_min}, {args.tau_max}]")
    print(f"output     : {args.output_dir}")
    print("-" * 88)

    # ------------------------------------------------------------
    # Fond W0
    # ------------------------------------------------------------
    W0 = load_matrix_csv(args.mi_file)
    n = W0.shape[0]
    k_eff = min(args.k, n - 1)

    mask0 = knn_mask(W0, k_eff)
    A0 = weighted_adjacency_from_mask(W0, mask0)
    lengths0 = build_edge_lengths(W0, mask0, eps=args.eps)
    Dsp0 = shortest_paths_floyd(lengths0)

    center = choose_center(W0, args.center)
    f = localized_profile_from_distances(Dsp0, center=center, sigma=args.sigma)

    # ------------------------------------------------------------
    # Excitation
    # ------------------------------------------------------------
    W1, deltaW = inject_local_excitation(
        W0,
        f=f,
        amp=args.amp,
        mask_mode=args.mask_mode,
        base_mask=mask0
    )

    # Graphe après excitation : on reconstruit avec le même k
    mask1 = knn_mask(W1, k_eff)
    A1 = weighted_adjacency_from_mask(W1, mask1)
    lengths1 = build_edge_lengths(W1, mask1, eps=args.eps)

    # ------------------------------------------------------------
    # Observables fond / excité
    # ------------------------------------------------------------
    ds0, ds0_r2 = local_spectral_dimension(A0, args.tau_min, args.tau_max, args.n_tau)
    ds1, ds1_r2 = local_spectral_dimension(A1, args.tau_min, args.tau_max, args.n_tau)

    edge0, R0 = edge_geometry_table(label + "_W0", W0, A0, lengths0, mask0, ds0)
    edge1, R1 = edge_geometry_table(label + "_W1", W1, A1, lengths1, mask1, ds1)

    rho_exc = node_excitation_density(deltaW)
    delta_ds = ds1 - ds0

    node_rows = []
    for i in range(n):
        node_rows.append({
            "label": label,
            "node": i,
            "center": center,
            "profile_f": f[i],
            "rho_exc": rho_exc[i],
            "T_matter_node": rho_exc[i] ** 2,
            "ds0": ds0[i],
            "ds1": ds1[i],
            "delta_ds": delta_ds[i],
            "abs_delta_ds": abs(delta_ds[i]),
            "ds0_r2": ds0_r2[i],
            "ds1_r2": ds1_r2[i],
            "R0_node": R0[i],
            "R1_node": R1[i] if i < len(R1) else np.nan,
            "delta_R_node": (R1[i] - R0[i]) if i < len(R1) and np.isfinite(R1[i]) and np.isfinite(R0[i]) else np.nan,
        })

    node_df = pd.DataFrame(node_rows)

    edge_df = edge_excitation_table(edge0, edge1, deltaW)

    # ------------------------------------------------------------
    # Corrélations principales
    # ------------------------------------------------------------
    corr_node_ds = safe_corr(node_df["rho_exc"], node_df["abs_delta_ds"])
    corr_node_R = safe_corr(node_df["T_matter_node"], np.abs(node_df["delta_R_node"]))

    corr_edge_G = safe_corr(edge_df["T_matter_edge"], edge_df["abs_delta_G_proxy"])
    corr_edge_kappa = safe_corr(edge_df["T_matter_edge"], edge_df["abs_delta_kappa"])
    corr_edge_R = safe_corr(edge_df["T_matter_edge"], edge_df["abs_delta_R_edge"])

    # Corr avec signe, utile pour diagnostic
    corr_edge_G_signed = safe_corr(edge_df["T_matter_edge"], edge_df["delta_G_proxy"])

    # Localisation
    finite_d = np.isfinite(Dsp0[center])
    d_center = Dsp0[center].copy()
    d_center[~finite_d] = np.nanmax(d_center[finite_d])

    if np.sum(rho_exc) > 1e-15:
        R_exc2 = np.sum(rho_exc * d_center**2) / np.sum(rho_exc)
        R_exc = float(np.sqrt(R_exc2))
        M_exc = float(np.sum(rho_exc))
    else:
        R_exc = np.nan
        M_exc = 0.0

    graph_radius = float(np.nanmax(d_center[finite_d]))
    localization_ratio = R_exc / graph_radius if graph_radius > 0 else np.nan

    # ------------------------------------------------------------
    # Sauvegardes
    # ------------------------------------------------------------
    node_csv = os.path.join(args.output_dir, "node_excitation_table.csv")
    edge_csv = os.path.join(args.output_dir, "edge_response_table.csv")
    W0_csv = os.path.join(args.output_dir, "W0.csv")
    W1_csv = os.path.join(args.output_dir, "W1_excited.csv")
    delta_csv = os.path.join(args.output_dir, "deltaW_excitation.csv")
    summary_json = os.path.join(args.output_dir, "summary.json")

    node_df.to_csv(node_csv, index=False)
    edge_df.to_csv(edge_csv, index=False)
    np.savetxt(W0_csv, W0, delimiter=",")
    np.savetxt(W1_csv, W1, delimiter=",")
    np.savetxt(delta_csv, deltaW, delimiter=",")

    summary = {
        "config": {
            "mi_file": args.mi_file,
            "label": label,
            "N": n,
            "lambda": lam_guess,
            "k": args.k,
            "amp": args.amp,
            "sigma": args.sigma,
            "center": center,
            "mask_mode": args.mask_mode,
            "tau_min": args.tau_min,
            "tau_max": args.tau_max,
            "n_tau": args.n_tau,
        },
        "localization": {
            "M_exc": M_exc,
            "R_exc": R_exc,
            "graph_radius": graph_radius,
            "localization_ratio": localization_ratio,
        },
        "correlations": {
            "node_rho_vs_abs_delta_ds": corr_node_ds,
            "node_Tmatter_vs_abs_delta_R": corr_node_R,
            "edge_Tmatter_vs_abs_delta_G": corr_edge_G,
            "edge_Tmatter_vs_delta_G_signed": corr_edge_G_signed,
            "edge_Tmatter_vs_abs_delta_kappa": corr_edge_kappa,
            "edge_Tmatter_vs_abs_delta_R": corr_edge_R,
        },
        "interpretation": {
            "positive_node_corr": "rho_exc corrélé à |delta d_s| signifie que l'excitation modifie localement la dimension spectrale.",
            "positive_edge_corr": "T_matter_edge corrélé à |delta G_proxy| signifie que l'excitation produit une réponse de courbure.",
            "localization_ratio": "R_exc / R_graph << 1 indique une excitation localisée."
        }
    }

    with open(summary_json, "w") as fjson:
        json.dump(summary, fjson, indent=2)

    make_figures(node_df, edge_df, args.output_dir)

    # ------------------------------------------------------------
    # Impression résultat
    # ------------------------------------------------------------
    print("RÉSULTATS — LOCALISATION")
    print("-" * 88)
    print(f"center              : {center}")
    print(f"M_exc               : {M_exc:.6g}")
    print(f"R_exc               : {R_exc:.6g}")
    print(f"graph_radius        : {graph_radius:.6g}")
    print(f"localization ratio  : {localization_ratio:.6g}")

    print("-" * 88)
    print("RÉSULTATS — CORRÉLATIONS NODALES")
    print("-" * 88)
    print(
        "rho_exc vs |delta d_s|     : "
        f"Pearson={corr_node_ds['pearson_r']:+.6f} "
        f"p={corr_node_ds['pearson_p']:.3g} "
        f"Spearman={corr_node_ds['spearman_r']:+.6f} "
        f"p={corr_node_ds['spearman_p']:.3g}"
    )
    print(
        "T_node vs |delta R_node|   : "
        f"Pearson={corr_node_R['pearson_r']:+.6f} "
        f"p={corr_node_R['pearson_p']:.3g} "
        f"Spearman={corr_node_R['spearman_r']:+.6f} "
        f"p={corr_node_R['spearman_p']:.3g}"
    )

    print("-" * 88)
    print("RÉSULTATS — CORRÉLATIONS ARÊTES")
    print("-" * 88)
    print(
        "T_edge vs |delta G_proxy|  : "
        f"Pearson={corr_edge_G['pearson_r']:+.6f} "
        f"p={corr_edge_G['pearson_p']:.3g} "
        f"Spearman={corr_edge_G['spearman_r']:+.6f} "
        f"p={corr_edge_G['spearman_p']:.3g}"
    )
    print(
        "T_edge vs delta G signed   : "
        f"Pearson={corr_edge_G_signed['pearson_r']:+.6f} "
        f"p={corr_edge_G_signed['pearson_p']:.3g} "
        f"Spearman={corr_edge_G_signed['spearman_r']:+.6f} "
        f"p={corr_edge_G_signed['spearman_p']:.3g}"
    )
    print(
        "T_edge vs |delta kappa|    : "
        f"Pearson={corr_edge_kappa['pearson_r']:+.6f} "
        f"p={corr_edge_kappa['pearson_p']:.3g} "
        f"Spearman={corr_edge_kappa['spearman_r']:+.6f} "
        f"p={corr_edge_kappa['spearman_p']:.3g}"
    )
    print(
        "T_edge vs |delta R_edge|   : "
        f"Pearson={corr_edge_R['pearson_r']:+.6f} "
        f"p={corr_edge_R['pearson_p']:.3g} "
        f"Spearman={corr_edge_R['spearman_r']:+.6f} "
        f"p={corr_edge_R['spearman_p']:.3g}"
    )

    if (
        np.isfinite(corr_node_ds["spearman_r"])
        and np.isfinite(corr_edge_G["spearman_r"])
        and corr_node_ds["spearman_r"] > 0
        and corr_edge_G["spearman_r"] > 0
    ):
        verdict = "signal positif : excitation localisée produit une réponse géométrique"
    else:
        verdict = "signal faible ou instable : tester amp/sigma/k ou proxy"

    print("-" * 88)
    print(f"VERDICT : {verdict}")
    print("-" * 88)
    print("Fichiers générés")
    print("-" * 88)
    print(node_csv)
    print(edge_csv)
    print(summary_json)
    print(os.path.join(args.output_dir, "fig_node_rho_vs_delta_ds.png"))
    print(os.path.join(args.output_dir, "fig_edge_Tmatter_vs_deltaG.png"))
    print(os.path.join(args.output_dir, "fig_node_profiles.png"))
    print("DONE")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()

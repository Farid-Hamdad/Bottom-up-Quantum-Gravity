#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_lambda_critical_scan_v1.py

BuP Paper 6 — Scan critique de lambda_c

Objectif :
    Tester si le paramètre lambda agit comme un paramètre d'ordre
    contrôlant l'émergence de la matière localisée dans BuP.

Observable principale :
    O_matter(lambda) =
        Spearman( T_ij^matter_proxy , |delta G_ij^proxy| )

où :
    T_ij^matter_proxy = (delta W_ij)^2
    G_ij^proxy = -kappa_ij + 1/2 R_edge

Interprétation :
    O_matter(lambda) > 0  : phase locale supportant une source matérielle effective
    O_matter(lambda) < 0  : phase non locale / réponse anti-corrélée
    O_matter(lambda)=0    : seuil critique candidat lambda_c

Usage :
    python3 bup_lambda_critical_scan_v1.py \
      --mi-files "results_mi_N20_full/*.csv" \
      --k 5 \
      --amp 0.15 \
      --sigma 1.0 \
      --output-dir results_lambda_critical_scan_v1

Avec des matrices plus fines autour du seuil :
    python3 bup_lambda_critical_scan_v1.py \
      --mi-files "results_mi_N20_critical/*.csv" \
      --k 5 \
      --amp 0.15 \
      --sigma 1.0 \
      --output-dir results_lambda_critical_scan_fine
"""

import os
import re
import glob
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
    """
    Exemples acceptés :
        MI_N20_lam0.57
        MI_N20_lambda0.57
        N20_lam0.64_seed1
    """
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

    return int(np.argmax(strength))


def localized_profile_from_distances(Dsp, center, sigma):
    d = Dsp[center].copy()

    finite = np.isfinite(d)
    if not np.any(finite):
        raise RuntimeError("Distances infinies : graphe non connecté ?")

    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2.0 * sigma ** 2))

    if np.max(f) > 0:
        f = f / np.max(f)

    return f


def inject_local_excitation(W0, f, amp, mode="positive", mask_mode="all", base_mask=None):
    W0 = clean_matrix(W0)
    sign = +1.0 if mode == "positive" else -1.0

    positive_vals = W0[W0 > 0]
    scale = float(np.mean(positive_vals)) if positive_vals.size else 1.0

    raw_delta = sign * abs(amp) * scale * np.outer(f, f)
    np.fill_diagonal(raw_delta, 0.0)

    if mask_mode == "edges":
        if base_mask is None:
            raise ValueError("base_mask requis si mask_mode='edges'")
        raw_delta = raw_delta * base_mask

    W_raw = W0 + raw_delta

    clipped_mask = W_raw < 0
    n_clipped = int(np.sum(clipped_mask))
    total_entries = int(W_raw.size - W_raw.shape[0])

    W1 = W_raw.copy()
    W1[W1 < 0] = 0.0
    W1 = 0.5 * (W1 + W1.T)
    np.fill_diagonal(W1, 0.0)

    actual_delta = W1 - W0

    clipping = {
        "n_clipped_entries": n_clipped,
        "total_offdiag_entries": total_entries,
        "clip_fraction": float(n_clipped / max(total_entries, 1)),
        "actual_delta_sum": float(np.sum(actual_delta)),
        "actual_delta_abs_sum": float(np.sum(np.abs(actual_delta))),
        "actual_delta_min": float(np.min(actual_delta)),
        "actual_delta_max": float(np.max(actual_delta)),
    }

    return W1, actual_delta, clipping


def node_excitation_density(deltaW):
    return np.sum(np.abs(deltaW), axis=1)


def edge_response_table(edge_df0, edge_df1, deltaW):
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
        on=["i", "j"],
        how="inner"
    )

    rows = []

    for _, row in merged.iterrows():
        i = int(row["i"])
        j = int(row["j"])

        dWij = deltaW[i, j]
        Tmatter = dWij ** 2

        deltaG = row["G1_proxy"] - row["G0_proxy"]
        deltaR = row["R1_edge"] - row["R0_edge"]
        deltaKappa = row["kappa1_ij"] - row["kappa0_ij"]

        rows.append({
            **row.to_dict(),
            "deltaW_ij": dWij,
            "abs_deltaW_ij": abs(dWij),
            "T_matter_edge": Tmatter,
            "delta_G_proxy": deltaG,
            "abs_delta_G_proxy": abs(deltaG),
            "delta_R_edge": deltaR,
            "abs_delta_R_edge": abs(deltaR),
            "delta_kappa": deltaKappa,
            "abs_delta_kappa": abs(deltaKappa),
        })

    return pd.DataFrame(rows)


# ============================================================
# RUN UNIQUE
# ============================================================

def run_single(path, args):
    label = parse_label(path)
    N_guess, lam = infer_N_lambda(label)

    W0 = load_matrix_csv(path)
    n = W0.shape[0]
    k_eff = min(args.k, n - 1)

    mask0 = knn_mask(W0, k_eff)
    A0 = weighted_adjacency_from_mask(W0, mask0)
    lengths0 = build_edge_lengths(W0, mask0, eps=args.eps)
    Dsp0 = shortest_paths_floyd(lengths0)

    center = choose_center(W0, args.center)
    f = localized_profile_from_distances(Dsp0, center=center, sigma=args.sigma)

    W1, deltaW, clipping = inject_local_excitation(
        W0,
        f=f,
        amp=args.amp,
        mode=args.mode,
        mask_mode=args.mask_mode,
        base_mask=mask0
    )

    mask1 = knn_mask(W1, k_eff)
    A1 = weighted_adjacency_from_mask(W1, mask1)
    lengths1 = build_edge_lengths(W1, mask1, eps=args.eps)

    ds0, ds0_r2 = local_spectral_dimension(A0, args.tau_min, args.tau_max, args.n_tau)
    ds1, ds1_r2 = local_spectral_dimension(A1, args.tau_min, args.tau_max, args.n_tau)

    edge0, R0 = edge_geometry_table(label + "_W0", W0, A0, lengths0, mask0, ds0)
    edge1, R1 = edge_geometry_table(label + "_W1", W1, A1, lengths1, mask1, ds1)

    edge_df = edge_response_table(edge0, edge1, deltaW)

    corr_edge_absG = safe_corr(edge_df["T_matter_edge"], edge_df["abs_delta_G_proxy"])
    corr_edge_signedG = safe_corr(edge_df["T_matter_edge"], edge_df["delta_G_proxy"])
    corr_edge_absR = safe_corr(edge_df["T_matter_edge"], edge_df["abs_delta_R_edge"])

    rho_exc = node_excitation_density(deltaW)

    finite_d = np.isfinite(Dsp0[center])
    d_center = Dsp0[center].copy()
    d_center[~finite_d] = np.nanmax(d_center[finite_d])

    M_exc = float(np.sum(rho_exc))
    if M_exc > 1e-15:
        R_exc = float(np.sqrt(np.sum(rho_exc * d_center ** 2) / M_exc))
    else:
        R_exc = np.nan

    graph_radius = float(np.nanmax(d_center[finite_d]))
    localization_ratio = R_exc / graph_radius if graph_radius > 0 else np.nan

    row = {
        "label": label,
        "path": path,
        "N": n,
        "lambda": lam,
        "k": args.k,
        "mode": args.mode,
        "amp": args.amp,
        "sigma": args.sigma,
        "center": center,
        "mask_mode": args.mask_mode,
        "M_exc": M_exc,
        "R_exc": R_exc,
        "graph_radius": graph_radius,
        "localization_ratio": localization_ratio,
        "clip_fraction": clipping["clip_fraction"],
        "delta_sum": clipping["actual_delta_sum"],

        "O_matter": corr_edge_absG["spearman_r"],
        "O_matter_p": corr_edge_absG["spearman_p"],
        "O_matter_n": corr_edge_absG["n"],

        "pearson_absG": corr_edge_absG["pearson_r"],
        "pearson_absG_p": corr_edge_absG["pearson_p"],

        "spearman_signedG": corr_edge_signedG["spearman_r"],
        "spearman_signedG_p": corr_edge_signedG["spearman_p"],
        "pearson_signedG": corr_edge_signedG["pearson_r"],
        "pearson_signedG_p": corr_edge_signedG["pearson_p"],

        "spearman_absR": corr_edge_absR["spearman_r"],
        "spearman_absR_p": corr_edge_absR["spearman_p"],
    }

    return row


# ============================================================
# ESTIMATION LAMBDA_C
# ============================================================

def estimate_lambda_c(df):
    """
    Estime lambda_c par interpolation linéaire entre deux points
    dont O_matter change de signe.

    Si plusieurs crossings, prend celui le plus proche de O=0.
    """
    sub = df[["lambda", "O_matter"]].dropna().sort_values("lambda")

    if len(sub) < 2:
        return None

    vals = sub.values
    candidates = []

    for a, b in zip(vals[:-1], vals[1:]):
        lam1, o1 = a
        lam2, o2 = b

        if o1 == 0:
            candidates.append((lam1, lam1, lam2, o1, o2))
        elif o1 * o2 < 0:
            # interpolation linéaire
            lc = lam1 + (0.0 - o1) * (lam2 - lam1) / (o2 - o1)
            candidates.append((lc, lam1, lam2, o1, o2))

    if not candidates:
        return None

    # prend crossing dont les deux valeurs sont les plus proches de zéro en moyenne
    candidates.sort(key=lambda x: abs(x[3]) + abs(x[4]))
    lc, lam1, lam2, o1, o2 = candidates[0]

    return {
        "lambda_c_estimate": float(lc),
        "bracket": [float(lam1), float(lam2)],
        "O_bracket": [float(o1), float(o2)]
    }


# ============================================================
# FIGURES
# ============================================================

def make_figures(df, output_dir):
    ensure_dir(output_dir)

    plot_df = df.dropna(subset=["lambda", "O_matter"]).sort_values("lambda")

    # Figure 1 : paramètre d'ordre
    plt.figure(figsize=(7.5, 5.2))
    plt.plot(
        plot_df["lambda"],
        plot_df["O_matter"],
        marker="o",
        linewidth=2,
        label=r"$\mathcal{O}_{matter}(\lambda)$"
    )
    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\mathcal{O}_{matter}=\rho_{\rm Spearman}(T_{edge},|\delta G|)$")
    plt.title("BuP Paper 6 — candidat paramètre d'ordre matière/localité")
    plt.grid(True, alpha=0.3)

    lc_info = estimate_lambda_c(plot_df)
    if lc_info is not None:
        lc = lc_info["lambda_c_estimate"]
        plt.axvline(lc, linestyle="--", linewidth=1.5, label=rf"$\lambda_c\simeq {lc:.3f}$")
        plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_order_parameter_vs_lambda.png"), dpi=230)
    plt.close()

    # Figure 2 : localisation
    plt.figure(figsize=(7.5, 5.2))
    plt.plot(
        plot_df["lambda"],
        plot_df["localization_ratio"],
        marker="s",
        linewidth=2
    )
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$R_{exc}/R_{graph}$")
    plt.title("Localisation de l'excitation vs lambda")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_localization_vs_lambda.png"), dpi=230)
    plt.close()

    # Figure 3 : signed response
    plt.figure(figsize=(7.5, 5.2))
    plt.plot(
        plot_df["lambda"],
        plot_df["spearman_signedG"],
        marker="o",
        linewidth=2,
        label=r"$\rho(T_{edge},\delta G)$"
    )
    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\lambda$")
    plt.ylabel("Spearman signed response")
    plt.title("Réponse signée de courbure vs lambda")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_signed_response_vs_lambda.png"), dpi=230)
    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-files", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_lambda_critical_scan_v1")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--mode", choices=["positive", "negative"], default="positive")
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--mask-mode", choices=["all", "edges"], default="all")

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--n-tau", type=int, default=24)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    files = sorted(glob.glob(args.mi_files))

    print("=" * 92)
    print("BuP Paper 6 — Lambda critical scan v1")
    print("=" * 92)
    print(f"MI files   : {args.mi_files}")
    print(f"count      : {len(files)}")
    print(f"k          : {args.k}")
    print(f"mode       : {args.mode}")
    print(f"amp        : {args.amp}")
    print(f"sigma      : {args.sigma}")
    print(f"output     : {args.output_dir}")
    print("-" * 92)

    if not files:
        raise FileNotFoundError(f"Aucun fichier trouvé : {args.mi_files}")

    rows = []

    for path in files:
        try:
            row = run_single(path, args)
            rows.append(row)

            print(
                f"OK {row['label']:25s} "
                f"lambda={row['lambda']} "
                f"O={row['O_matter']:+.6f} "
                f"p={row['O_matter_p']:.3g} "
                f"loc={row['localization_ratio']:.3f} "
                f"signed={row['spearman_signedG']:+.4f}"
            )

        except Exception as e:
            print(f"FAIL {path}: {e}")

    if not rows:
        raise RuntimeError("Aucun run exploitable.")

    df = pd.DataFrame(rows)

    # Si plusieurs fichiers ont même lambda, résumé par lambda
    agg = df.groupby("lambda", dropna=False).agg(
        n_runs=("O_matter", "count"),
        O_mean=("O_matter", "mean"),
        O_median=("O_matter", "median"),
        O_std=("O_matter", "std"),
        O_frac_pos=("O_matter", lambda x: float(np.mean(np.asarray(x) > 0))),
        p_median=("O_matter_p", "median"),
        loc_mean=("localization_ratio", "mean"),
        loc_median=("localization_ratio", "median"),
        signed_mean=("spearman_signedG", "mean"),
        signed_median=("spearman_signedG", "median"),
    ).reset_index()

    # Pour estimation lambda_c, utiliser la médiane par lambda
    lc_info = estimate_lambda_c(
        agg.rename(columns={"O_median": "O_matter"})[["lambda", "O_matter"]]
    )

    # Sauvegardes
    scan_csv = os.path.join(args.output_dir, "lambda_critical_scan.csv")
    agg_csv = os.path.join(args.output_dir, "lambda_critical_summary_by_lambda.csv")
    summary_json = os.path.join(args.output_dir, "summary.json")

    df.to_csv(scan_csv, index=False)
    agg.to_csv(agg_csv, index=False)

    make_figures(
        agg.rename(columns={
            "O_median": "O_matter",
            "loc_median": "localization_ratio",
            "signed_median": "spearman_signedG",
        }),
        args.output_dir
    )

    summary = {
        "config": {
            "mi_files": args.mi_files,
            "k": args.k,
            "mode": args.mode,
            "amp": args.amp,
            "sigma": args.sigma,
            "center": args.center,
            "mask_mode": args.mask_mode,
            "tau_min": args.tau_min,
            "tau_max": args.tau_max,
            "n_tau": args.n_tau,
        },
        "lambda_c": lc_info,
        "definition": {
            "O_matter": "Spearman(T_edge_matter_proxy, abs(delta G_proxy))",
            "phase_positive": "O_matter > 0 : localized excitation behaves as effective matter source",
            "phase_negative": "O_matter < 0 : non-local / anti-correlated response"
        },
        "n_runs": int(len(df)),
    }

    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print("-" * 92)
    print("RÉSULTAT GLOBAL")
    print("-" * 92)
    if lc_info is None:
        print("lambda_c : non détecté dans l'intervalle scanné")
    else:
        print(f"lambda_c estimate : {lc_info['lambda_c_estimate']:.6f}")
        print(f"bracket           : {lc_info['bracket']}")
        print(f"O bracket         : {lc_info['O_bracket']}")

    print("-" * 92)
    print("Fichiers générés")
    print("-" * 92)
    print(scan_csv)
    print(agg_csv)
    print(summary_json)
    print(os.path.join(args.output_dir, "fig_order_parameter_vs_lambda.png"))
    print(os.path.join(args.output_dir, "fig_localization_vs_lambda.png"))
    print(os.path.join(args.output_dir, "fig_signed_response_vs_lambda.png"))
    print("DONE")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_einstein_tensor_correlation_v2.py

BuP Paper 5/6 — Scan de conventions pour la corrélation :

    G_ij^ent        : proxy discret de courbure / Einstein
    T_ij^ent[d_s]  : proxy discret construit depuis d_s

Nouveauté v2 :
    - teste plusieurs conventions de G_ij
    - teste plusieurs formes de T_ij
    - produit un tableau de classement des corrélations
    - cherche les combinaisons stables en Pearson et Spearman

Usage :
    python3 bup_einstein_tensor_correlation_v2.py \
      --mi-files "results_mi_matrices/*.csv" \
      --k 5 \
      --tau-min 0.01 \
      --tau-max 50 \
      --output-dir results_einstein_corr_v2

    python3 bup_einstein_tensor_correlation_v2.py \
      --mi-files "results_mi_N20_full/*.csv" \
      --k 5 \
      --tau-min 0.01 \
      --tau-max 50 \
      --output-dir results_einstein_corr_N20_v2
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


def graph_laplacian_field(A: np.ndarray, field: np.ndarray):
    deg = np.sum(A, axis=1)
    out = np.zeros_like(field, dtype=float)

    for i in range(A.shape[0]):
        if deg[i] > 1e-15:
            out[i] = np.sum(A[i] * (field - field[i])) / deg[i]
        else:
            out[i] = 0.0
    return out


# ============================================================
# OLLIVIER-RICCI PROXY
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


# ============================================================
# EDGE OBSERVABLES
# ============================================================

def build_base_edge_table(label, W, A, lengths, mask, ds, ds_r2):
    n = W.shape[0]

    kappa, Dsp = ollivier_ricci_proxy(A, lengths, mask)
    R = node_scalar_curvature(n, kappa)
    lap_ds = graph_laplacian_field(A, ds)

    rows = []

    for (i, j), kij in kappa.items():
        ell = lengths[i, j]
        if not np.isfinite(ell) or ell <= 1e-12:
            continue

        Ri = R[i]
        Rj = R[j]
        Redge = np.nanmean([Ri, Rj])

        grad_ds = (ds[j] - ds[i]) / ell
        grad2 = grad_ds ** 2
        lap_edge = 0.5 * (lap_ds[i] + lap_ds[j])
        abs_lap_edge = abs(lap_edge)
        abs_grad = abs(grad_ds)

        rows.append({
            "label": label,
            "N": n,
            "i": i,
            "j": j,
            "W_ij": W[i, j],
            "ell_ij": ell,
            "kappa": kij,
            "minus_kappa": -kij,
            "R_i": Ri,
            "R_j": Rj,
            "R_edge": Redge,
            "minus_R_edge": -Redge,
            "G_kappa_minus_halfR": kij - 0.5 * Redge,
            "G_minus_kappa_plus_halfR": -kij + 0.5 * Redge,
            "G_kappa_plus_halfR": kij + 0.5 * Redge,
            "G_minus_kappa_minus_halfR": -kij - 0.5 * Redge,
            "ds_i": ds[i],
            "ds_j": ds[j],
            "ds_edge_mean": 0.5 * (ds[i] + ds[j]),
            "grad_ds": grad_ds,
            "abs_grad_ds": abs_grad,
            "grad_ds_sq": grad2,
            "lap_ds_i": lap_ds[i],
            "lap_ds_j": lap_ds[j],
            "lap_ds_edge": lap_edge,
            "abs_lap_ds_edge": abs_lap_edge,
            "T_grad2_plus_lap": grad2 + lap_edge,
            "T_grad2_minus_lap": grad2 - lap_edge,
            "T_grad2_plus_abslap": grad2 + abs_lap_edge,
            "T_grad2_only": grad2,
            "T_lap_only": lap_edge,
            "T_minus_lap_only": -lap_edge,
            "T_abs_lap_only": abs_lap_edge,
            "T_abs_grad_plus_abslap": abs_grad + abs_lap_edge,
            "ds_r2_i": ds_r2[i],
            "ds_r2_j": ds_r2[j],
        })

    return pd.DataFrame(rows)


def build_node_table(label, A, ds, ds_r2, edge_df):
    n = A.shape[0]
    node_rows = []

    for i in range(n):
        sub = edge_df[(edge_df["i"] == i) | (edge_df["j"] == i)]
        R_node = np.nanmean(
            np.concatenate([
                sub.loc[sub["i"] == i, "R_i"].values,
                sub.loc[sub["j"] == i, "R_j"].values,
            ])
        ) if len(sub) else np.nan

        node_rows.append({
            "label": label,
            "N": n,
            "node": i,
            "strength": float(np.sum(A[i])),
            "ds": float(ds[i]),
            "ds_r2": float(ds_r2[i]),
            "R_node": float(R_node) if np.isfinite(R_node) else np.nan,
        })

    return pd.DataFrame(node_rows)


# ============================================================
# SCAN DES COMBINAISONS
# ============================================================

G_COLUMNS = [
    "kappa",
    "minus_kappa",
    "R_edge",
    "minus_R_edge",
    "G_kappa_minus_halfR",
    "G_minus_kappa_plus_halfR",
    "G_kappa_plus_halfR",
    "G_minus_kappa_minus_halfR",
]

T_COLUMNS = [
    "T_grad2_plus_lap",
    "T_grad2_minus_lap",
    "T_grad2_plus_abslap",
    "T_grad2_only",
    "T_lap_only",
    "T_minus_lap_only",
    "T_abs_lap_only",
    "T_abs_grad_plus_abslap",
]


def scan_correlations(edge_df: pd.DataFrame, group_cols=None):
    if group_cols is None:
        group_cols = []

    rows = []

    if group_cols:
        grouped = edge_df.groupby(group_cols, dropna=False)
    else:
        grouped = [((), edge_df)]

    for key, sub in grouped:
        if not isinstance(key, tuple):
            key = (key,)

        group_info = {}
        for c, v in zip(group_cols, key):
            group_info[c] = v

        for gcol in G_COLUMNS:
            for tcol in T_COLUMNS:
                corr = safe_corr(sub[tcol], sub[gcol])

                score = np.nan
                if np.isfinite(corr["pearson_r"]) and np.isfinite(corr["spearman_r"]):
                    # Score favorise corr positive et cohérence Pearson/Spearman
                    score = 0.5 * corr["pearson_r"] + 0.5 * corr["spearman_r"]

                rows.append({
                    **group_info,
                    "G_proxy": gcol,
                    "T_proxy": tcol,
                    **corr,
                    "score_mean_r": score,
                    "abs_pearson": abs(corr["pearson_r"]) if np.isfinite(corr["pearson_r"]) else np.nan,
                    "abs_spearman": abs(corr["spearman_r"]) if np.isfinite(corr["spearman_r"]) else np.nan,
                })

    return pd.DataFrame(rows)


def stability_summary(scan_by_file: pd.DataFrame):
    """
    Résume la stabilité d'une combinaison G/T sur tous les fichiers.
    """
    rows = []

    for (g, t), sub in scan_by_file.groupby(["G_proxy", "T_proxy"]):
        pr = sub["pearson_r"].values.astype(float)
        sr = sub["spearman_r"].values.astype(float)

        pr_good = pr[np.isfinite(pr)]
        sr_good = sr[np.isfinite(sr)]

        if len(pr_good) == 0 or len(sr_good) == 0:
            continue

        rows.append({
            "G_proxy": g,
            "T_proxy": t,
            "n_files": int(len(sub)),
            "pearson_mean": float(np.nanmean(pr)),
            "pearson_median": float(np.nanmedian(pr)),
            "pearson_frac_pos": float(np.mean(pr_good > 0)),
            "spearman_mean": float(np.nanmean(sr)),
            "spearman_median": float(np.nanmedian(sr)),
            "spearman_frac_pos": float(np.mean(sr_good > 0)),
            "score": float(
                0.25 * np.nanmean(pr)
                + 0.25 * np.nanmedian(pr)
                + 0.25 * np.nanmean(sr)
                + 0.25 * np.nanmedian(sr)
            ),
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("score", ascending=False)
    return df


# ============================================================
# FIGURES
# ============================================================

def make_figures(edge_all, global_scan, stability_df, output_dir):
    ensure_dir(output_dir)

    # Heatmap score global
    if not global_scan.empty:
        pivot = global_scan.pivot(index="G_proxy", columns="T_proxy", values="score_mean_r")

        plt.figure(figsize=(12, 6))
        plt.imshow(pivot.values, aspect="auto")
        plt.colorbar(label="score = 0.5 Pearson + 0.5 Spearman")
        plt.xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=75, ha="right", fontsize=8)
        plt.yticks(np.arange(len(pivot.index)), pivot.index, fontsize=8)
        plt.title("Scan global des conventions G/T")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fig_scan_heatmap.png"), dpi=220)
        plt.close()

    # Top combo scatter
    if stability_df is not None and not stability_df.empty:
        best = stability_df.iloc[0]
        gcol = best["G_proxy"]
        tcol = best["T_proxy"]

        x = edge_all[tcol].values
        y = edge_all[gcol].values
        good = np.isfinite(x) & np.isfinite(y)

        plt.figure(figsize=(7, 5))
        plt.scatter(x[good], y[good], s=18, alpha=0.7)
        plt.xlabel(tcol)
        plt.ylabel(gcol)
        plt.title(f"Meilleure combinaison stable : {gcol} vs {tcol}")
        plt.grid(True, alpha=0.3)

        if np.sum(good) >= 4 and np.std(x[good]) > 0:
            coeff = np.polyfit(x[good], y[good], 1)
            xx = np.linspace(np.min(x[good]), np.max(x[good]), 200)
            yy = coeff[0] * xx + coeff[1]
            plt.plot(xx, yy, linewidth=2)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fig_best_G_vs_T.png"), dpi=220)
        plt.close()

    # Stability top 10
    if stability_df is not None and not stability_df.empty:
        top = stability_df.head(10).copy()
        labels = [f"{r.G_proxy}\nvs\n{r.T_proxy}" for _, r in top.iterrows()]
        vals = top["score"].values

        plt.figure(figsize=(11, 6))
        plt.bar(np.arange(len(vals)), vals)
        plt.axhline(0, linewidth=1)
        plt.xticks(np.arange(len(vals)), labels, rotation=70, ha="right", fontsize=8)
        plt.ylabel("stability score")
        plt.title("Top 10 combinaisons G/T")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fig_top10_stability.png"), dpi=220)
        plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-files", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_einstein_corr_v2")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--n-tau", type=int, default=24)
    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    files = sorted(glob.glob(args.mi_files))

    print("=" * 86)
    print("BuP Paper 5/6 — V2 scan conventions G_ent_ij vs T_ent_ij[d_s]")
    print("=" * 86)
    print(f"MI files   : {args.mi_files}")
    print(f"Nombre     : {len(files)}")
    print(f"k          : {args.k}")
    print(f"tau        : [{args.tau_min}, {args.tau_max}]")
    print(f"output     : {args.output_dir}")
    print("-" * 86)

    if not files:
        raise FileNotFoundError(f"Aucun fichier trouvé : {args.mi_files}")

    all_edges = []
    all_nodes = []
    file_info = []

    for path in files:
        label = parse_label(path)
        N_guess, lam_guess = infer_N_lambda(label)

        try:
            W = load_matrix_csv(path)
            n = W.shape[0]
            k_eff = min(args.k, n - 1)

            mask = knn_mask(W, k_eff)
            A = weighted_adjacency_from_mask(W, mask)
            lengths = build_edge_lengths(W, mask, eps=args.eps)

            ds, ds_r2 = local_spectral_dimension(
                A,
                tau_min=args.tau_min,
                tau_max=args.tau_max,
                n_tau=args.n_tau
            )

            edge_df = build_base_edge_table(
                label=label,
                W=W,
                A=A,
                lengths=lengths,
                mask=mask,
                ds=ds,
                ds_r2=ds_r2
            )

            if edge_df.empty:
                print(f"SKIP {label}: aucune arête exploitable")
                continue

            edge_df["lambda"] = lam_guess
            edge_df["file_N"] = n

            node_df = build_node_table(label, A, ds, ds_r2, edge_df)
            node_df["lambda"] = lam_guess
            node_df["file_N"] = n

            all_edges.append(edge_df)
            all_nodes.append(node_df)

            file_info.append({
                "label": label,
                "path": path,
                "N": n,
                "lambda": lam_guess,
                "edges": int(len(edge_df)),
                "ds_mean": float(np.nanmean(ds)),
                "ds_std": float(np.nanstd(ds)),
                "ds_r2_mean": float(np.nanmean(ds_r2)),
            })

            print(
                f"OK {label:30s} "
                f"N={n:2d} edges={len(edge_df):3d} "
                f"ds={np.nanmean(ds):.4f}±{np.nanstd(ds):.4f}"
            )

        except Exception as e:
            print(f"FAIL {label}: {e}")

    if not all_edges:
        raise RuntimeError("Aucune matrice exploitable.")

    edge_all = pd.concat(all_edges, ignore_index=True)
    node_all = pd.concat(all_nodes, ignore_index=True)
    file_df = pd.DataFrame(file_info)

    # Scan global
    global_scan = scan_correlations(edge_all)
    global_scan = global_scan.sort_values("score_mean_r", ascending=False)

    # Scan par fichier
    scan_by_file = scan_correlations(edge_all, group_cols=["label", "file_N", "lambda"])

    # Stabilité par combinaison
    stability = stability_summary(scan_by_file)

    # Scan par N
    scan_by_N = scan_correlations(edge_all, group_cols=["file_N"])
    scan_by_N = scan_by_N.sort_values(["file_N", "score_mean_r"], ascending=[True, False])

    # Sauvegardes
    edge_csv = os.path.join(args.output_dir, "edge_table_v2.csv")
    node_csv = os.path.join(args.output_dir, "node_table_v2.csv")
    file_csv = os.path.join(args.output_dir, "files_summary.csv")
    global_csv = os.path.join(args.output_dir, "scan_global.csv")
    byfile_csv = os.path.join(args.output_dir, "scan_by_file.csv")
    byN_csv = os.path.join(args.output_dir, "scan_by_N.csv")
    stability_csv = os.path.join(args.output_dir, "stability_summary.csv")
    json_path = os.path.join(args.output_dir, "summary.json")

    edge_all.to_csv(edge_csv, index=False)
    node_all.to_csv(node_csv, index=False)
    file_df.to_csv(file_csv, index=False)
    global_scan.to_csv(global_csv, index=False)
    scan_by_file.to_csv(byfile_csv, index=False)
    scan_by_N.to_csv(byN_csv, index=False)
    stability.to_csv(stability_csv, index=False)

    make_figures(edge_all, global_scan, stability, args.output_dir)

    best_global = global_scan.iloc[0].to_dict() if not global_scan.empty else {}
    best_stable = stability.iloc[0].to_dict() if not stability.empty else {}

    summary = {
        "config": {
            "mi_files": args.mi_files,
            "k": args.k,
            "tau_min": args.tau_min,
            "tau_max": args.tau_max,
            "n_tau": args.n_tau,
            "eps": args.eps,
        },
        "n_files_ok": int(len(file_df)),
        "n_edges_total": int(len(edge_all)),
        "best_global": best_global,
        "best_stable": best_stable,
        "interpretation": {
            "best_global": "Meilleure corrélation sur toutes les arêtes empilées.",
            "best_stable": "Meilleure combinaison stable en moyenne/médiane sur les fichiers.",
            "positive_stability": "Un score stable positif suggère une convention G/T cohérente.",
            "negative_or_mixed": "Un signe instable indique proxy discret ou convention à raffiner."
        }
    }

    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("-" * 86)
    print("MEILLEURE COMBINAISON GLOBALE")
    print("-" * 86)
    if best_global:
        print(f"G proxy      : {best_global['G_proxy']}")
        print(f"T proxy      : {best_global['T_proxy']}")
        print(f"n            : {best_global['n']}")
        print(f"Pearson r    : {best_global['pearson_r']:+.6f}")
        print(f"Pearson p    : {best_global['pearson_p']:.6g}")
        print(f"Spearman rho : {best_global['spearman_r']:+.6f}")
        print(f"Spearman p   : {best_global['spearman_p']:.6g}")
        print(f"score        : {best_global['score_mean_r']:+.6f}")

    print("-" * 86)
    print("MEILLEURE COMBINAISON STABLE")
    print("-" * 86)
    if best_stable:
        print(f"G proxy             : {best_stable['G_proxy']}")
        print(f"T proxy             : {best_stable['T_proxy']}")
        print(f"n_files             : {best_stable['n_files']}")
        print(f"Pearson mean/median : {best_stable['pearson_mean']:+.6f} / {best_stable['pearson_median']:+.6f}")
        print(f"Pearson frac pos    : {best_stable['pearson_frac_pos']:.3f}")
        print(f"Spearman mean/med   : {best_stable['spearman_mean']:+.6f} / {best_stable['spearman_median']:+.6f}")
        print(f"Spearman frac pos   : {best_stable['spearman_frac_pos']:.3f}")
        print(f"score               : {best_stable['score']:+.6f}")

    print("-" * 86)
    print("Fichiers générés")
    print("-" * 86)
    print(edge_csv)
    print(node_csv)
    print(file_csv)
    print(global_csv)
    print(byfile_csv)
    print(byN_csv)
    print(stability_csv)
    print(json_path)
    print(os.path.join(args.output_dir, "fig_scan_heatmap.png"))
    print(os.path.join(args.output_dir, "fig_best_G_vs_T.png"))
    print(os.path.join(args.output_dir, "fig_top10_stability.png"))
    print("DONE")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()

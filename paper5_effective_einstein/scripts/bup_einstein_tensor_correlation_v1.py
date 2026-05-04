#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_einstein_tensor_correlation_v1.py

BuP Paper 5 / Paper 6 — Test numérique de corrélation entre :

    G_ij^ent        : proxy discret du tenseur d'Einstein émergent
    T_ij^ent[d_s]  : proxy discret du tenseur d'intrication construit depuis d_s

Objectif :
    Tester si la courbure émergente du graphe d'intrication est corrélée
    avec les gradients / variations de dimension spectrale locale.

Entrée :
    Matrices MI W_ij au format CSV.

Sorties :
    results_einstein_corr/
        edge_table.csv
        node_table.csv
        summary.json
        fig_G_vs_T.png
        fig_corr_by_file.png

Usage :
    python3 bup_einstein_tensor_correlation_v1.py \
        --mi-files "results_mi_matrices/*.csv" \
        --output-dir results_einstein_corr

    python3 bup_einstein_tensor_correlation_v1.py \
        --mi-files "results_mi_N20_full/*.csv" \
        --k 5 \
        --tau-min 0.01 \
        --tau-max 50 \
        --output-dir results_einstein_corr_N20
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


def load_matrix_csv(path: str) -> np.ndarray:
    """
    Charge une matrice CSV.
    Accepte :
      - CSV brut numérique
      - CSV avec index/colonnes
    """
    try:
        W = np.loadtxt(path, delimiter=",")
        if W.ndim == 2 and W.shape[0] == W.shape[1]:
            return clean_matrix(W)
    except Exception:
        pass

    df = pd.read_csv(path, header=None)
    W = df.values

    # Si première ligne/colonne non numérique, retry pandas avec header/index
    if not np.issubdtype(W.dtype, np.number):
        df = pd.read_csv(path, index_col=0)
        W = df.values

    W = W.astype(float)
    return clean_matrix(W)


def clean_matrix(W: np.ndarray) -> np.ndarray:
    W = np.asarray(W, dtype=float)

    if W.shape[0] != W.shape[1]:
        raise ValueError(f"Matrice non carrée : shape={W.shape}")

    W = np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)

    # symétrisation
    W = 0.5 * (W + W.T)

    # diagonale nulle
    np.fill_diagonal(W, 0.0)

    # poids négatifs interdits
    W[W < 0] = 0.0

    return W


def parse_label(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def infer_N_lambda(label: str):
    """
    Extrait N et lambda si le nom ressemble à :
      MI_N20_lam0.57.csv
      MI_N16_lambda0.40_rep1.csv
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


# ============================================================
# GRAPHE LOCAL
# ============================================================

def knn_mask(W: np.ndarray, k: int) -> np.ndarray:
    """
    Construit un masque kNN symétrisé.
    Garde les k plus grands voisins par nœud.
    """
    n = W.shape[0]
    mask = np.zeros_like(W, dtype=bool)

    for i in range(n):
        idx = np.argsort(W[i])[::-1]
        idx = [j for j in idx if j != i and W[i, j] > 0]
        for j in idx[:k]:
            mask[i, j] = True

    # symétrisation : si i choisit j OU j choisit i
    mask = np.logical_or(mask, mask.T)
    np.fill_diagonal(mask, False)
    return mask


def build_edge_lengths(W: np.ndarray, mask: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Longueur d'intrication :
        ell_ij = -log((W_ij + eps)/(Wmax + eps))
    uniquement sur les arêtes du masque.
    """
    Lell = np.full_like(W, np.inf, dtype=float)

    vals = W[mask]
    if vals.size == 0 or np.max(vals) <= 0:
        raise ValueError("Aucune arête positive dans le graphe.")

    Wmax = np.max(vals)

    Lell[mask] = -np.log((W[mask] + eps) / (Wmax + eps))

    # éviter longueur exactement nulle pour l'arête max
    positive = np.isfinite(Lell)
    Lell[positive] = np.maximum(Lell[positive], 1e-9)

    np.fill_diagonal(Lell, 0.0)
    return Lell


def weighted_adjacency_from_mask(W: np.ndarray, mask: np.ndarray) -> np.ndarray:
    A = np.zeros_like(W)
    A[mask] = W[mask]
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


# ============================================================
# DIMENSION SPECTRALE LOCALE
# ============================================================

def graph_laplacian(A: np.ndarray, normalized: bool = True) -> np.ndarray:
    deg = np.sum(A, axis=1)

    if normalized:
        invsqrt = np.zeros_like(deg)
        good = deg > 1e-15
        invsqrt[good] = 1.0 / np.sqrt(deg[good])
        Dm = np.diag(invsqrt)
        L = np.eye(A.shape[0]) - Dm @ A @ Dm
        # pour nœuds isolés, remettre ligne cohérente
        for i, g in enumerate(good):
            if not g:
                L[i, i] = 0.0
        return L

    return np.diag(deg) - A


def local_spectral_dimension(A: np.ndarray, tau_min: float, tau_max: float, n_tau: int = 24):
    """
    Calcule d_s(i) depuis :
        P_i(tau) = [exp(-tau L)]_ii
        d_s = -2 d log P / d log tau

    Fit linéaire log(P_i) vs log(tau) pour chaque nœud.
    """
    n = A.shape[0]
    L = graph_laplacian(A, normalized=True)

    taus = np.geomspace(tau_min, tau_max, n_tau)
    logt = np.log(taus)

    P = np.zeros((n_tau, n))

    for a, tau in enumerate(taus):
        K = expm(-tau * L)
        diag = np.clip(np.diag(K), 1e-300, None)
        P[a, :] = diag

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
# COURBURE DISCRÈTE : PROXY OLLIVIER SIMPLIFIÉ
# ============================================================

def shortest_paths_floyd(lengths: np.ndarray) -> np.ndarray:
    """
    Floyd-Warshall simple pour N petit.
    """
    D = lengths.copy()
    n = D.shape[0]

    for i in range(n):
        D[i, i] = 0.0

    for k in range(n):
        D = np.minimum(D, D[:, [k]] + D[[k], :])

    return D


def neighborhood_measure(A: np.ndarray, i: int):
    """
    Mesure m_i(j) = A_ij / sum_k A_ik.
    Retourne indices et masses.
    """
    idx = np.where(A[i] > 0)[0]
    if idx.size == 0:
        return np.array([i], dtype=int), np.array([1.0])
    weights = A[i, idx].astype(float)
    s = np.sum(weights)
    if s <= 0:
        return np.array([i], dtype=int), np.array([1.0])
    return idx, weights / s


def wasserstein_greedy(idx_a, mass_a, idx_b, mass_b, D):
    """
    Approximation greedy de W1 entre deux mesures discrètes.

    Ce n'est pas un solveur optimal exact, mais suffisant comme proxy rapide.
    Pour N petit, on pourrait remplacer par scipy.optimize.linprog.
    """
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
    """
    Courbure d'Ollivier-Ricci simplifiée :
        kappa_ij = 1 - W1(m_i,m_j)/d_ij

    Retourne dict {(i,j): kappa}
    """
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


# ============================================================
# PROXIES G_ij ET T_ij
# ============================================================

def node_scalar_curvature(n: int, kappa: dict):
    """
    R_i = moyenne des kappa_ij incidents.
    """
    vals = [[] for _ in range(n)]
    for (i, j), k in kappa.items():
        vals[i].append(k)
        vals[j].append(k)

    R = np.zeros(n)
    for i in range(n):
        if vals[i]:
            R[i] = float(np.mean(vals[i]))
        else:
            R[i] = np.nan
    return R


def graph_laplacian_field(A: np.ndarray, field: np.ndarray):
    """
    Laplacien discret pondéré :
        Δ f_i = sum_j A_ij (f_j - f_i) / sum_j A_ij
    normalisé par degré pour éviter l'effet taille.
    """
    deg = np.sum(A, axis=1)
    out = np.zeros_like(field, dtype=float)

    for i in range(A.shape[0]):
        if deg[i] > 1e-15:
            out[i] = np.sum(A[i] * (field - field[i])) / deg[i]
        else:
            out[i] = 0.0
    return out


def build_edge_observables(label, W, A, lengths, mask, ds, ds_r2, alpha=1.0, beta=1.0):
    """
    Construit edge_table avec :
      G_ent_ij : proxy Einstein discret
      T_ent_ij : proxy tenseur spectral discret
    """
    n = W.shape[0]

    kappa, Dsp = ollivier_ricci_proxy(A, lengths, mask)
    R = node_scalar_curvature(n, kappa)
    lap_ds = graph_laplacian_field(A, ds)

    rows = []

    for (i, j), kij in kappa.items():
        ell = lengths[i, j]
        if not np.isfinite(ell) or ell <= 1e-12:
            continue

        Rij = np.nanmean([R[i], R[j]])

        # proxy Einstein discret :
        # G_ij = kappa_ij - 1/2 R_edge
        # Le facteur métrique exact est inconnu ; on garde la forme centrée.
        Gij = kij - 0.5 * Rij

        grad_ds = (ds[j] - ds[i]) / ell
        lap_edge = 0.5 * (lap_ds[i] + lap_ds[j])

        # Proxy T_ent :
        # α (∇d)^2 + β Δd
        Tij = alpha * grad_ds**2 + beta * lap_edge

        rows.append({
            "label": label,
            "N": n,
            "i": i,
            "j": j,
            "W_ij": W[i, j],
            "ell_ij": ell,
            "kappa_ij": kij,
            "R_i": R[i],
            "R_j": R[j],
            "R_edge": Rij,
            "G_ent_ij": Gij,
            "ds_i": ds[i],
            "ds_j": ds[j],
            "grad_ds_ij": grad_ds,
            "lap_ds_i": lap_ds[i],
            "lap_ds_j": lap_ds[j],
            "T_ent_ij": Tij,
            "ds_r2_i": ds_r2[i],
            "ds_r2_j": ds_r2[j],
        })

    node_rows = []
    for i in range(n):
        node_rows.append({
            "label": label,
            "N": n,
            "node": i,
            "strength": float(np.sum(A[i])),
            "ds": float(ds[i]),
            "ds_r2": float(ds_r2[i]),
            "R_node": float(R[i]) if np.isfinite(R[i]) else np.nan,
            "lap_ds": float(lap_ds[i]),
        })

    return pd.DataFrame(rows), pd.DataFrame(node_rows)


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)

    x = x[good]
    y = y[good]

    if x.size < 4:
        return {
            "n": int(x.size),
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
        }

    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return {
            "n": int(x.size),
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
        }

    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)

    return {
        "n": int(x.size),
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "spearman_r": float(sr),
        "spearman_p": float(sp),
    }


# ============================================================
# FIGURES
# ============================================================

def make_figures(edge_df: pd.DataFrame, summary_by_file: pd.DataFrame, output_dir: str):
    if edge_df.empty:
        return

    # Figure 1 : G vs T
    plt.figure(figsize=(7, 5))
    x = edge_df["T_ent_ij"].values
    y = edge_df["G_ent_ij"].values
    good = np.isfinite(x) & np.isfinite(y)

    plt.scatter(x[good], y[good], s=18, alpha=0.7)
    plt.xlabel(r"$T_{ij}^{\rm ent}[d_s]$ proxy")
    plt.ylabel(r"$G_{ij}^{\rm ent}$ proxy")
    plt.title("BuP — corrélation Einstein discret / tenseur spectral")
    plt.grid(True, alpha=0.3)

    if np.sum(good) >= 4:
        coeff = np.polyfit(x[good], y[good], 1)
        xx = np.linspace(np.min(x[good]), np.max(x[good]), 200)
        yy = coeff[0] * xx + coeff[1]
        plt.plot(xx, yy, linewidth=2)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_G_vs_T.png"), dpi=200)
    plt.close()

    # Figure 2 : corr par fichier
    if summary_by_file is not None and not summary_by_file.empty:
        plt.figure(figsize=(9, 5))
        labels = summary_by_file["label"].astype(str).tolist()
        vals = summary_by_file["pearson_r"].values
        xidx = np.arange(len(vals))

        plt.bar(xidx, vals)
        plt.axhline(0, linewidth=1)
        plt.xticks(xidx, labels, rotation=75, ha="right", fontsize=8)
        plt.ylabel("Pearson r")
        plt.title("Corrélation par matrice MI")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fig_corr_by_file.png"), dpi=200)
        plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-files", type=str, required=True,
                        help='Glob des matrices MI, ex: "results_mi_matrices/*.csv"')

    parser.add_argument("--output-dir", type=str, default="results_einstein_corr")

    parser.add_argument("--k", type=int, default=5,
                        help="Nombre de voisins locaux kNN.")

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--n-tau", type=int, default=24)

    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=1.0)

    parser.add_argument("--eps", type=float, default=1e-12)

    args = parser.parse_args()

    ensure_dir(args.output_dir)

    files = sorted(glob.glob(args.mi_files))

    print("=" * 78)
    print("BuP Paper 5/6 — Corrélation G_ent_ij vs T_ent_ij[d_s]")
    print("=" * 78)
    print(f"MI files     : {args.mi_files}")
    print(f"Nombre       : {len(files)}")
    print(f"k local      : {args.k}")
    print(f"tau window   : [{args.tau_min}, {args.tau_max}]")
    print(f"alpha,beta   : {args.alpha}, {args.beta}")
    print(f"output       : {args.output_dir}")
    print("-" * 78)

    if not files:
        raise FileNotFoundError(f"Aucun fichier trouvé pour : {args.mi_files}")

    all_edges = []
    all_nodes = []
    file_summaries = []

    for path in files:
        label = parse_label(path)
        N_guess, lam_guess = infer_N_lambda(label)

        try:
            W = load_matrix_csv(path)
            n = W.shape[0]

            if n < 5:
                print(f"SKIP {path} : N trop petit ({n})")
                continue

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

            edge_df, node_df = build_edge_observables(
                label=label,
                W=W,
                A=A,
                lengths=lengths,
                mask=mask,
                ds=ds,
                ds_r2=ds_r2,
                alpha=args.alpha,
                beta=args.beta
            )

            if edge_df.empty:
                print(f"SKIP {path} : aucune arête exploitable")
                continue

            edge_df["lambda"] = lam_guess
            node_df["lambda"] = lam_guess

            corr = safe_corr(edge_df["T_ent_ij"], edge_df["G_ent_ij"])

            file_summaries.append({
                "label": label,
                "path": path,
                "N": n,
                "lambda": lam_guess,
                "n_edges": int(len(edge_df)),
                "ds_mean": float(np.nanmean(ds)),
                "ds_std": float(np.nanstd(ds)),
                "R_node_mean": float(np.nanmean(node_df["R_node"])),
                "R_node_std": float(np.nanstd(node_df["R_node"])),
                **corr
            })

            all_edges.append(edge_df)
            all_nodes.append(node_df)

            print(
                f"OK {label:30s} "
                f"N={n:2d} edges={len(edge_df):3d} "
                f"r={corr['pearson_r']:+.4f} "
                f"p={corr['pearson_p']:.3g} "
                f"rho={corr['spearman_r']:+.4f}"
            )

        except Exception as e:
            print(f"FAIL {path} : {e}")

    if not all_edges:
        raise RuntimeError("Aucune matrice exploitable.")

    edge_all = pd.concat(all_edges, ignore_index=True)
    node_all = pd.concat(all_nodes, ignore_index=True)
    summary_df = pd.DataFrame(file_summaries)

    global_corr = safe_corr(edge_all["T_ent_ij"], edge_all["G_ent_ij"])

    # Sauvegardes
    edge_csv = os.path.join(args.output_dir, "edge_table.csv")
    node_csv = os.path.join(args.output_dir, "node_table.csv")
    summary_csv = os.path.join(args.output_dir, "summary_by_file.csv")
    summary_json = os.path.join(args.output_dir, "summary.json")

    edge_all.to_csv(edge_csv, index=False)
    node_all.to_csv(node_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)

    summary = {
        "config": {
            "mi_files": args.mi_files,
            "k": args.k,
            "tau_min": args.tau_min,
            "tau_max": args.tau_max,
            "n_tau": args.n_tau,
            "alpha": args.alpha,
            "beta": args.beta,
            "eps": args.eps,
        },
        "global": global_corr,
        "n_files_ok": int(len(summary_df)),
        "n_edges_total": int(len(edge_all)),
        "files": file_summaries,
        "interpretation": {
            "positive_corr": "Un r positif stable suggère que la courbure émergente est corrélée au tenseur spectral T_ent[d_s].",
            "weak_corr": "Une corrélation faible indique que le proxy G_ij ou T_ij doit être amélioré, ou que N est trop petit.",
            "negative_corr": "Une corrélation négative stable suggère une convention de signe à revoir ou une forme de T_ent incomplète."
        }
    }

    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    make_figures(edge_all, summary_df, args.output_dir)

    print("-" * 78)
    print("RÉSULTAT GLOBAL")
    print("-" * 78)
    print(f"N edges       : {global_corr['n']}")
    print(f"Pearson r     : {global_corr['pearson_r']:+.6f}")
    print(f"Pearson p     : {global_corr['pearson_p']:.6g}")
    print(f"Spearman rho  : {global_corr['spearman_r']:+.6f}")
    print(f"Spearman p    : {global_corr['spearman_p']:.6g}")

    if np.isfinite(global_corr["pearson_r"]):
        r = global_corr["pearson_r"]
        if r > 0.25:
            verdict = "corrélation positive non triviale"
        elif r > 0.05:
            verdict = "corrélation positive faible"
        elif r > -0.05:
            verdict = "corrélation compatible avec zéro"
        else:
            verdict = "corrélation négative ou convention de signe à vérifier"
    else:
        verdict = "corrélation non estimable"

    print(f"VERDICT       : {verdict}")
    print("-" * 78)
    print("Fichiers générés")
    print("-" * 78)
    print(f"Edges   : {edge_csv}")
    print(f"Nodes   : {node_csv}")
    print(f"CSV     : {summary_csv}")
    print(f"JSON    : {summary_json}")
    print(f"Figure  : {os.path.join(args.output_dir, 'fig_G_vs_T.png')}")
    print(f"Figure  : {os.path.join(args.output_dir, 'fig_corr_by_file.png')}")
    print("DONE")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_action_sim_v2.py

BuP Paper 6/7 — Action pure BuP v2 avec stabilisation d'excitation.

Différence avec v1 :
    v1 relaxait la géométrie mais ne stabilisait pas l'excitation.
    v2 ajoute deux contraintes :

        S_mass = chi_mass * (M_exc - M0)^2
        S_size = chi_size * (R_exc - R0)^2

    avec :
        M_exc = sum_ij (W_ij - W0_ij)^2
        R_exc = rayon effectif de l'excitation autour du centre initial

Objectif :
    Tester si une excitation localisée peut rester cohérente pendant que
    S_BuP diminue.

Usage conseillé :
    python3 bup_action_sim_v2.py \
      --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
      --k 5 \
      --inject-excitation \
      --amp 0.15 \
      --sigma 1.0 \
      --steps 20 \
      --eta 0.001 \
      --d-star 0.9 \
      --alpha-spec 1.0 \
      --beta-geom 0.1 \
      --gamma-loc 0.01 \
      --chi-mass 1.0 \
      --chi-size 1.0 \
      --output-dir results_bup_action_v2_lam057
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
from scipy.stats import spearmanr


# ============================================================
# Utilitaires
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def clean_matrix(W):
    W = np.asarray(W, dtype=float)
    if W.shape[0] != W.shape[1]:
        raise ValueError(f"Matrice non carrée : {W.shape}")

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


def infer_lambda(label):
    m = re.search(r"lam(?:bda)?([0-9]+(?:\.[0-9]+)?)", label)
    if m:
        return float(m.group(1))
    return None


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]

    if len(x) < 4:
        return np.nan, np.nan
    if np.std(x) < 1e-14 or np.std(y) < 1e-14:
        return np.nan, np.nan

    r, p = spearmanr(x, y)
    return float(r), float(p)


# ============================================================
# Graphe
# ============================================================

def knn_mask(W, k):
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


def adjacency_from_mask(W, mask):
    A = np.zeros_like(W)
    A[mask] = W[mask]
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def edge_lengths(W, mask, eps=1e-12):
    Lell = np.full_like(W, np.inf, dtype=float)
    vals = W[mask]

    if vals.size == 0 or np.max(vals) <= 0:
        raise ValueError("Aucune arête positive.")

    Wmax = np.max(vals)
    Lell[mask] = -np.log((W[mask] + eps) / (Wmax + eps))
    good = np.isfinite(Lell)
    Lell[good] = np.maximum(Lell[good], 1e-9)
    np.fill_diagonal(Lell, 0.0)
    return Lell


def floyd_shortest_paths(lengths):
    D = lengths.copy()
    n = D.shape[0]

    for i in range(n):
        D[i, i] = 0.0

    for kk in range(n):
        D = np.minimum(D, D[:, [kk]] + D[[kk], :])

    return D


def graph_laplacian(A, normalized=True):
    deg = np.sum(A, axis=1)

    if normalized:
        invsqrt = np.zeros_like(deg)
        good = deg > 1e-15
        invsqrt[good] = 1.0 / np.sqrt(deg[good])
        Dm = np.diag(invsqrt)
        L = np.eye(A.shape[0]) - Dm @ A @ Dm

        for i, ok in enumerate(good):
            if not ok:
                L[i, i] = 0.0
        return L

    return np.diag(deg) - A


# ============================================================
# Dimension spectrale
# ============================================================

def spectral_dimension_mean(A, tau_min, tau_max, n_tau):
    L = graph_laplacian(A, normalized=True)
    taus = np.geomspace(tau_min, tau_max, n_tau)

    Z = []
    for tau in taus:
        K = expm(-tau * L)
        Z.append(np.trace(K))

    Z = np.clip(np.asarray(Z), 1e-300, None)
    logt = np.log(taus)
    logz = np.log(Z)

    coeff = np.polyfit(logt, logz, 1)
    slope = coeff[0]
    d_mean = -2.0 * slope

    fit = coeff[0] * logt + coeff[1]
    ss_res = np.sum((logz - fit) ** 2)
    ss_tot = np.sum((logz - np.mean(logz)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else np.nan

    return float(d_mean), float(r2)


# ============================================================
# Courbure Ollivier-Ricci approx
# ============================================================

def neighborhood_measure(A, i):
    idx = np.where(A[i] > 0)[0]

    if idx.size == 0:
        return np.array([i], dtype=int), np.array([1.0])

    w = A[i, idx].astype(float)
    s = np.sum(w)

    if s <= 0:
        return np.array([i], dtype=int), np.array([1.0])

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


def ollivier_ricci_proxy(A, lengths, mask):
    Dsp = floyd_shortest_paths(lengths)
    n = A.shape[0]
    measures = [neighborhood_measure(A, i) for i in range(n)]

    kappa = {}
    edges = np.argwhere(np.triu(mask, 1))

    for i, j in edges:
        dij = Dsp[i, j]

        if not np.isfinite(dij) or dij <= 1e-12:
            continue

        idx_i, mass_i = measures[i]
        idx_j, mass_j = measures[j]

        w1 = wasserstein_greedy(idx_i, mass_i, idx_j, mass_j, Dsp)
        kij = 1.0 - w1 / dij
        kappa[(int(i), int(j))] = float(kij)

    return kappa, Dsp


def node_curvature(n, kappa):
    vals = [[] for _ in range(n)]

    for (i, j), k in kappa.items():
        vals[i].append(k)
        vals[j].append(k)

    R = np.zeros(n)
    for i in range(n):
        R[i] = np.mean(vals[i]) if vals[i] else np.nan

    return R


def edge_geometric_quantities(W, k, eps=1e-12):
    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    lengths = edge_lengths(W, mask, eps=eps)
    kappa, Dsp = ollivier_ricci_proxy(A, lengths, mask)
    R = node_curvature(W.shape[0], kappa)

    rows = []
    for (i, j), kij in kappa.items():
        Redge = np.nanmean([R[i], R[j]])
        Gproxy = -kij + 0.5 * Redge
        rows.append({
            "i": i,
            "j": j,
            "kappa": kij,
            "R_edge": Redge,
            "G_proxy": Gproxy,
            "W": W[i, j],
            "d_ent": Dsp[i, j],
        })

    return pd.DataFrame(rows), A, mask, Dsp


# ============================================================
# Excitation
# ============================================================

def choose_center(W):
    strength = np.sum(W, axis=1)
    return int(np.argmax(strength))


def inject_excitation(W0, Dsp, amp, sigma, center=None):
    W0 = clean_matrix(W0)
    n = W0.shape[0]

    if center is None or center < 0:
        center = choose_center(W0)

    d = Dsp[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2 * sigma ** 2))
    f = f / max(np.max(f), 1e-15)

    positive_vals = W0[W0 > 0]
    scale = float(np.mean(positive_vals)) if positive_vals.size else 1.0

    deltaW = amp * scale * np.outer(f, f)
    np.fill_diagonal(deltaW, 0.0)

    W1 = clean_matrix(W0 + deltaW)
    return W1, deltaW, center


def excitation_density(deltaW):
    return np.sum(np.abs(deltaW), axis=1)


def excitation_mass(deltaW):
    return float(np.sum(deltaW ** 2))


def excitation_radius(deltaW, Dsp, center):
    rho = excitation_density(deltaW)
    total = np.sum(rho)
    if total <= 1e-15:
        return np.nan

    d = Dsp[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite])
    d[~finite] = dmax

    R2 = np.sum(rho * d ** 2) / total
    return float(np.sqrt(max(R2, 0.0)))


# ============================================================
# Matter response
# ============================================================

def matter_response(W_ref, W, args):
    k_eff = min(args.k, W.shape[0] - 1)

    edge0, _, _, _ = edge_geometric_quantities(W_ref, k_eff, eps=args.eps)
    edge1, _, _, _ = edge_geometric_quantities(W, k_eff, eps=args.eps)

    a = edge0.rename(columns={"G_proxy": "G0"})
    b = edge1.rename(columns={"G_proxy": "G1"})

    merged = pd.merge(
        a[["i", "j", "G0"]],
        b[["i", "j", "G1"]],
        on=["i", "j"],
        how="inner"
    )

    deltaW = W - W_ref

    vals_T = []
    vals_dG_abs = []
    vals_dG_signed = []

    for _, row in merged.iterrows():
        i = int(row["i"])
        j = int(row["j"])

        T = deltaW[i, j] ** 2
        dG = row["G1"] - row["G0"]

        vals_T.append(T)
        vals_dG_abs.append(abs(dG))
        vals_dG_signed.append(dG)

    rho_abs, p_abs = safe_spearman(vals_T, vals_dG_abs)
    rho_signed, p_signed = safe_spearman(vals_T, vals_dG_signed)

    return {
        "O_matter_absG": rho_abs,
        "O_matter_absG_p": p_abs,
        "O_matter_signedG": rho_signed,
        "O_matter_signedG_p": p_signed,
        "n_edges_common": int(len(vals_T)),
    }


# ============================================================
# Action BuP v2
# ============================================================

def action_components(W, W_ref, Dsp_ref, center, M0, R0, args):
    W = clean_matrix(W)
    n = W.shape[0]
    k_eff = min(args.k, n - 1)

    edge_df, A, mask, Dsp = edge_geometric_quantities(W, k_eff, eps=args.eps)

    d_mean, ds_r2 = spectral_dimension_mean(
        A,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        n_tau=args.n_tau
    )

    # 1. Spectral
    S_spec = (d_mean - args.d_star) ** 2

    # 2. Geometry
    if len(edge_df) > 0:
        S_geom = float(np.mean(edge_df["G_proxy"].values ** 2))
    else:
        S_geom = 1e6

    # 3. Locality
    if len(edge_df) > 0:
        S_loc = float(np.mean(edge_df["W"].values * edge_df["d_ent"].values ** 2))
    else:
        S_loc = 1e6

    # 4. Excitation mass and size relative to W_ref
    deltaW = W - W_ref
    M_exc = excitation_mass(deltaW)
    R_exc = excitation_radius(deltaW, Dsp_ref, center)

    S_mass = (M_exc - M0) ** 2

    if np.isfinite(R_exc) and np.isfinite(R0):
        S_size = (R_exc - R0) ** 2
    else:
        S_size = 0.0

    # 5. Optional smoothness of excitation over nodes
    # rho_i = nodal excitation density; penalize sharp noisy distributions
    rho = excitation_density(deltaW)
    S_smooth = 0.0
    if args.chi_smooth != 0:
        # use adjacency of reference kNN graph from W_ref
        mask_ref = knn_mask(W_ref, k_eff)
        Aref = adjacency_from_mask(W_ref, mask_ref)
        deg = np.sum(Aref, axis=1)
        for i in range(n):
            for j in range(i + 1, n):
                if Aref[i, j] > 0:
                    S_smooth += Aref[i, j] * (rho[i] - rho[j]) ** 2
        S_smooth = float(S_smooth / max(np.sum(Aref > 0), 1))

    S_total = (
        args.alpha_spec * S_spec
        + args.beta_geom * S_geom
        + args.gamma_loc * S_loc
        + args.chi_mass * S_mass
        + args.chi_size * S_size
        + args.chi_smooth * S_smooth
    )

    return {
        "S_total": float(S_total),
        "S_spec": float(S_spec),
        "S_geom": float(S_geom),
        "S_loc": float(S_loc),
        "S_mass": float(S_mass),
        "S_size": float(S_size),
        "S_smooth": float(S_smooth),
        "d_mean": float(d_mean),
        "ds_r2": float(ds_r2),
        "M_exc": float(M_exc),
        "R_exc": float(R_exc) if np.isfinite(R_exc) else np.nan,
        "M0": float(M0),
        "R0": float(R0),
        "n_edges": int(len(edge_df)),
    }


def finite_difference_gradient(W, W_ref, Dsp_ref, center, M0, R0, args):
    W = clean_matrix(W)
    n = W.shape[0]
    grad = np.zeros_like(W)

    base = action_components(W, W_ref, Dsp_ref, center, M0, R0, args)
    h = args.fd_eps

    for i in range(n):
        for j in range(i + 1, n):
            if args.grad_positive_edges_only and W[i, j] <= 0:
                continue

            Wp = W.copy()
            Wm = W.copy()

            Wp[i, j] += h
            Wp[j, i] += h

            Wm[i, j] -= h
            Wm[j, i] -= h
            Wm[Wm < 0] = 0.0
            np.fill_diagonal(Wm, 0.0)

            Sp = action_components(Wp, W_ref, Dsp_ref, center, M0, R0, args)["S_total"]
            Sm = action_components(Wm, W_ref, Dsp_ref, center, M0, R0, args)["S_total"]

            gij = (Sp - Sm) / (2.0 * h)
            grad[i, j] = gij
            grad[j, i] = gij

    return grad, base


def propose_step(W, grad, eta, clip_max):
    gnorm = np.linalg.norm(grad)
    if gnorm > 1e-15:
        grad_eff = grad / gnorm
    else:
        grad_eff = grad

    Wnew = W - eta * grad_eff
    Wnew = clean_matrix(Wnew)

    if clip_max is not None and clip_max > 0:
        Wnew = np.minimum(Wnew, clip_max)
        Wnew = clean_matrix(Wnew)

    return Wnew


# ============================================================
# Figures
# ============================================================

def make_figures(history_df, output_dir):
    ensure_dir(output_dir)

    plt.figure(figsize=(8, 5))
    plt.plot(history_df["step"], history_df["S_total"], marker="o")
    plt.xlabel("step")
    plt.ylabel(r"$S_{\rm BuP}^{(v2)}$")
    plt.title("Relaxation de l'action BuP v2")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_action_total.png"), dpi=220)
    plt.close()

    plt.figure(figsize=(9, 5))
    for col in ["S_spec", "S_geom", "S_loc", "S_mass", "S_size", "S_smooth"]:
        plt.plot(history_df["step"], history_df[col], marker="o", label=col)
    plt.xlabel("step")
    plt.ylabel("component")
    plt.title("Composantes de l'action BuP v2")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_action_components.png"), dpi=220)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(history_df["step"], history_df["M_exc"], marker="o", label=r"$M_{exc}$")
    plt.axhline(history_df["M0"].iloc[0], linestyle="--", label=r"$M_0$")
    plt.xlabel("step")
    plt.ylabel("mass")
    plt.title("Masse de l'excitation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_excitation_mass.png"), dpi=220)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(history_df["step"], history_df["R_exc"], marker="o", label=r"$R_{exc}$")
    plt.axhline(history_df["R0"].iloc[0], linestyle="--", label=r"$R_0$")
    plt.xlabel("step")
    plt.ylabel("radius")
    plt.title("Rayon de l'excitation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_excitation_radius.png"), dpi=220)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(history_df["step"], history_df["O_matter_absG"], marker="o", label=r"$\rho(T,|\delta G|)$")
    plt.plot(history_df["step"], history_df["O_matter_signedG"], marker="s", label=r"$\rho(T,\delta G)$")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("step")
    plt.ylabel("Spearman")
    plt.title("Réponse matière / courbure")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_matter_response.png"), dpi=220)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(history_df["step"], history_df["d_mean"], marker="o")
    plt.axhline(history_df["d_star"].iloc[0], linestyle="--")
    plt.xlabel("step")
    plt.ylabel(r"$d_s$")
    plt.title("Dimension spectrale moyenne")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_spectral_dimension.png"), dpi=220)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_bup_action_sim_v2")

    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--eta", type=float, default=0.001)
    parser.add_argument("--fd-eps", type=float, default=1e-4)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--n-tau", type=int, default=20)

    parser.add_argument("--d-star", type=float, default=0.9)

    parser.add_argument("--alpha-spec", type=float, default=1.0)
    parser.add_argument("--beta-geom", type=float, default=0.1)
    parser.add_argument("--gamma-loc", type=float, default=0.01)
    parser.add_argument("--chi-mass", type=float, default=1.0)
    parser.add_argument("--chi-size", type=float, default=1.0)
    parser.add_argument("--chi-smooth", type=float, default=0.0)

    parser.add_argument("--inject-excitation", action="store_true")
    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--center", type=int, default=-1)

    parser.add_argument("--clip-max", type=float, default=1.0)
    parser.add_argument("--eps", type=float, default=1e-12)

    parser.add_argument("--grad-positive-edges-only", action="store_true")
    parser.add_argument("--line-search", action="store_true",
                        help="Accepte uniquement les pas qui diminuent S_total.")

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    label = parse_label(args.mi_file)
    lam = infer_lambda(label)

    print("=" * 96)
    print("BuP Action Simulation v2 — excitation-stabilized")
    print("=" * 96)
    print(f"MI file    : {args.mi_file}")
    print(f"label      : {label}")
    print(f"lambda     : {lam}")
    print(f"k          : {args.k}")
    print(f"steps      : {args.steps}")
    print(f"eta        : {args.eta}")
    print(f"d_star     : {args.d_star}")
    print(f"chi_mass   : {args.chi_mass}")
    print(f"chi_size   : {args.chi_size}")
    print(f"chi_smooth : {args.chi_smooth}")
    print(f"line search: {args.line_search}")
    print(f"output     : {args.output_dir}")
    print("-" * 96)

    W0 = load_matrix_csv(args.mi_file)
    W_ref = W0.copy()

    k_eff = min(args.k, W0.shape[0] - 1)
    _, _, _, Dsp_ref = edge_geometric_quantities(W_ref, k_eff, eps=args.eps)

    W = W0.copy()

    if args.inject_excitation:
        W, deltaW0, center = inject_excitation(
            W0,
            Dsp_ref,
            amp=args.amp,
            sigma=args.sigma,
            center=args.center
        )
    else:
        deltaW0 = W - W_ref
        center = choose_center(W0)

    M0 = excitation_mass(deltaW0)
    R0 = excitation_radius(deltaW0, Dsp_ref, center)

    if M0 <= 1e-18:
        print("ATTENTION : M0≈0. La contrainte de masse n'aura pas d'effet sans excitation.")

    np.savetxt(os.path.join(args.output_dir, "W_ref.csv"), W_ref, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "W_initial.csv"), W, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "deltaW_initial.csv"), deltaW0, delimiter=",")

    print(f"center     : {center}")
    print(f"M0         : {M0:.6g}")
    print(f"R0         : {R0:.6g}")
    print("-" * 96)

    history = []

    for step in range(args.steps + 1):
        comps = action_components(W, W_ref, Dsp_ref, center, M0, R0, args)
        response = matter_response(W_ref, W, args)

        row = {
            "step": step,
            "lambda": lam,
            "d_star": args.d_star,
            **comps,
            **response,
        }
        history.append(row)

        print(
            f"step={step:03d} "
            f"S={comps['S_total']:.6g} "
            f"Sg={comps['S_geom']:.3g} "
            f"Sloc={comps['S_loc']:.3g} "
            f"Smass={comps['S_mass']:.3g} "
            f"Ssize={comps['S_size']:.3g} "
            f"M={comps['M_exc']:.3g}/{M0:.3g} "
            f"R={comps['R_exc']:.3g}/{R0:.3g} "
            f"ds={comps['d_mean']:.4f} "
            f"O={response['O_matter_absG']:+.4f}"
        )

        if step == args.steps:
            break

        grad, base = finite_difference_gradient(W, W_ref, Dsp_ref, center, M0, R0, args)

        if args.line_search:
            eta_try = args.eta
            S_old = comps["S_total"]
            accepted = False

            for _ in range(12):
                W_try = propose_step(W, grad, eta_try, args.clip_max)
                S_try = action_components(W_try, W_ref, Dsp_ref, center, M0, R0, args)["S_total"]

                if S_try <= S_old:
                    W = W_try
                    accepted = True
                    break

                eta_try *= 0.5

            if not accepted:
                print("  line-search: aucun pas accepté, arrêt anticipé.")
                break

        else:
            W = propose_step(W, grad, args.eta, args.clip_max)

    history_df = pd.DataFrame(history)
    history_df.to_csv(os.path.join(args.output_dir, "action_history.csv"), index=False)

    np.savetxt(os.path.join(args.output_dir, "W_final.csv"), W, delimiter=",")
    np.savetxt(os.path.join(args.output_dir, "deltaW_final.csv"), W - W_ref, delimiter=",")

    summary = {
        "config": {
            "mi_file": args.mi_file,
            "label": label,
            "lambda": lam,
            "k": args.k,
            "steps": args.steps,
            "eta": args.eta,
            "fd_eps": args.fd_eps,
            "d_star": args.d_star,
            "alpha_spec": args.alpha_spec,
            "beta_geom": args.beta_geom,
            "gamma_loc": args.gamma_loc,
            "chi_mass": args.chi_mass,
            "chi_size": args.chi_size,
            "chi_smooth": args.chi_smooth,
            "inject_excitation": bool(args.inject_excitation),
            "amp": args.amp,
            "sigma": args.sigma,
            "center": center,
            "M0": M0,
            "R0": R0,
            "line_search": bool(args.line_search),
        },
        "initial": history[0],
        "final": history[-1],
        "delta": {
            "S_total": history[-1]["S_total"] - history[0]["S_total"],
            "d_mean": history[-1]["d_mean"] - history[0]["d_mean"],
            "M_exc": history[-1]["M_exc"] - history[0]["M_exc"],
            "R_exc": history[-1]["R_exc"] - history[0]["R_exc"],
            "O_matter_absG": history[-1]["O_matter_absG"] - history[0]["O_matter_absG"],
        },
        "interpretation": {
            "S_mass": "Constrains excitation mass M_exc close to its initial value M0.",
            "S_size": "Constrains excitation radius R_exc close to its initial value R0.",
            "success_condition": "Good run if S decreases, d_s remains stable, M_exc and R_exc remain close to initial values, and O_matter stays positive."
        }
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    make_figures(history_df, args.output_dir)

    print("-" * 96)
    print("RÉSUMÉ")
    print("-" * 96)
    print(f"S initial : {history[0]['S_total']:.6g}")
    print(f"S final   : {history[-1]['S_total']:.6g}")
    print(f"Delta S   : {history[-1]['S_total'] - history[0]['S_total']:.6g}")
    print(f"d_s init  : {history[0]['d_mean']:.6g}")
    print(f"d_s final : {history[-1]['d_mean']:.6g}")
    print(f"M init    : {history[0]['M_exc']:.6g}")
    print(f"M final   : {history[-1]['M_exc']:.6g}")
    print(f"R init    : {history[0]['R_exc']:.6g}")
    print(f"R final   : {history[-1]['R_exc']:.6g}")
    print(f"O init    : {history[0]['O_matter_absG']:+.6f}")
    print(f"O final   : {history[-1]['O_matter_absG']:+.6f}")
    print("-" * 96)
    print("Fichiers générés")
    print(os.path.join(args.output_dir, "action_history.csv"))
    print(os.path.join(args.output_dir, "summary.json"))
    print(os.path.join(args.output_dir, "fig_action_total.png"))
    print(os.path.join(args.output_dir, "fig_excitation_mass.png"))
    print(os.path.join(args.output_dir, "fig_excitation_radius.png"))
    print(os.path.join(args.output_dir, "fig_matter_response.png"))
    print("DONE")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()

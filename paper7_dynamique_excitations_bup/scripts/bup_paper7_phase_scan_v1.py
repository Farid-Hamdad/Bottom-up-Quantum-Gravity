#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper7_phase_scan_v1.py

Paper 7 — Dynamique variationnelle et phases internes des excitations BuP.

Objectif :
    Explorer comment une excitation localisée δW évolue sous relaxation
    d'une action candidate S_BuP[W].

Idée centrale :
    La chute de O_matter n'est pas forcément un bug :
    elle peut indiquer une transformation de phase interne de l'excitation.

Ce script mesure :
    O_init(lambda)
    O_final(lambda)
    Delta_O(lambda)
    M_exc, R_exc
    d_s
    phase label : (++), (+-), (-+), (--)

Contrairement à l'ancien script hybride :
    - pas de fit gaussien de lambda_c ;
    - extraction par crossing si changement de signe ;
    - observable corrigée avec G_proxy = -kappa + 1/2 R_edge ;
    - d_star peut être auto ;
    - line search optionnelle ;
    - interprétation en diagramme de phases, pas en validation directe de Paper 6.

Usage :
    python3 bup_paper7_phase_scan_v1.py \
      --mi-dir results_mi_N20_full/ \
      --output-dir results_paper7_phase_scan_v1 \
      --steps 20 \
      --k 5 \
      --d-star auto \
      --line-search
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
from scipy.sparse.csgraph import shortest_path
from scipy.stats import spearmanr


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


def parse_lambda(path):
    name = os.path.basename(path)
    m = re.search(r"lam(?:bda)?([0-9]+(?:\.[0-9]+)?)", name)
    return float(m.group(1)) if m else None


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
# Graph geometry
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
        raise ValueError("Aucune arête positive.")

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


def graph_laplacian(A):
    deg = np.sum(A, axis=1)
    invsqrt = np.zeros_like(deg)
    good = deg > 1e-15
    invsqrt[good] = 1.0 / np.sqrt(deg[good])
    Dm = np.diag(invsqrt)
    L = np.eye(A.shape[0]) - Dm @ A @ Dm

    for i, ok in enumerate(good):
        if not ok:
            L[i, i] = 0.0

    return L


def spectral_dimension(A, tau_min, tau_max, tau_points):
    L = graph_laplacian(A)
    taus = np.geomspace(tau_min, tau_max, tau_points)

    Z = []
    for tau in taus:
        K = expm(-tau * L)
        Z.append(max(np.trace(K), 1e-300))

    logt = np.log(taus)
    logz = np.log(Z)

    coeff = np.polyfit(logt, logz, 1)
    ds = -2.0 * coeff[0]

    return float(ds)


# ============================================================
# Ollivier-Ricci proxy with greedy transport
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


def edge_geometry(W, k, tau_min, tau_max, tau_points, eps=1e-12):
    mask = knn_mask(W, k)
    A = adjacency_from_mask(W, mask)
    D = entanglement_distance(W, mask, eps=eps)
    kappas = ollivier_kappa(A, D, mask)

    n = W.shape[0]
    node_vals = [[] for _ in range(n)]

    for (i, j), kij in kappas.items():
        node_vals[i].append(kij)
        node_vals[j].append(kij)

    R_node = np.zeros(n)
    for i in range(n):
        R_node[i] = np.mean(node_vals[i]) if node_vals[i] else np.nan

    rows = []
    for (i, j), kij in kappas.items():
        R_edge = np.nanmean([R_node[i], R_node[j]])
        G_proxy = -kij + 0.5 * R_edge

        rows.append({
            "i": i,
            "j": j,
            "kappa": kij,
            "R_edge": R_edge,
            "G_proxy": G_proxy,
            "W": W[i, j],
            "d_ent": D[i, j]
        })

    ds = spectral_dimension(A, tau_min, tau_max, tau_points)

    return pd.DataFrame(rows), A, mask, D, ds


# ============================================================
# Excitation
# ============================================================

def choose_center(W):
    return int(np.argmax(np.sum(W, axis=1)))


def inject_excitation(W0, D0, amp=0.15, sigma=1.0, center=-1):
    W0 = clean_matrix(W0)
    n = W0.shape[0]

    if center is None or center < 0:
        center = choose_center(W0)

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    f = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    f = f / max(np.max(f), 1e-15)

    vals = W0[W0 > 0]
    scale = float(np.mean(vals)) if vals.size else 1.0

    deltaW = amp * scale * np.outer(f, f)
    np.fill_diagonal(deltaW, 0.0)

    W1 = clean_matrix(W0 + deltaW)
    return W1, deltaW, center


def excitation_density(deltaW):
    return np.sum(np.abs(deltaW), axis=1)


def excitation_mass(deltaW):
    return float(np.sum(deltaW ** 2))


def excitation_radius(deltaW, D0, center):
    rho = excitation_density(deltaW)
    total = np.sum(rho)

    if total <= 1e-15:
        return np.nan

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    return float(np.sqrt(np.sum(rho * d ** 2) / total))


def excitation_kurtosis(deltaW, D0, center):
    rho = excitation_density(deltaW)
    total = np.sum(rho)

    if total <= 1e-15:
        return np.nan

    d = D0[center].copy()
    finite = np.isfinite(d)
    dmax = np.nanmax(d[finite]) if np.any(finite) else 1.0
    d[~finite] = dmax

    m2 = np.sum(rho * d ** 2) / total
    m4 = np.sum(rho * d ** 4) / total

    if m2 <= 1e-15:
        return np.nan

    return float(m4 / (m2 ** 2))


# ============================================================
# Matter observable
# ============================================================

def matter_observable(W_ref, W, k, tau_min, tau_max, tau_points):
    geom0, _, _, _, _ = edge_geometry(W_ref, k, tau_min, tau_max, tau_points)
    geom1, _, _, _, _ = edge_geometry(W, k, tau_min, tau_max, tau_points)

    g0 = geom0[["i", "j", "G_proxy"]].rename(columns={"G_proxy": "G0"})
    g1 = geom1[["i", "j", "G_proxy"]].rename(columns={"G_proxy": "G1"})

    merged = pd.merge(g0, g1, on=["i", "j"], how="inner")

    deltaW = W - W_ref

    T_vals = []
    dG_vals = []
    dG_signed = []

    for _, row in merged.iterrows():
        i = int(row["i"])
        j = int(row["j"])

        T = deltaW[i, j] ** 2
        dG = row["G1"] - row["G0"]

        T_vals.append(T)
        dG_vals.append(abs(dG))
        dG_signed.append(dG)

    rho_abs, p_abs = safe_spearman(T_vals, dG_vals)
    rho_signed, p_signed = safe_spearman(T_vals, dG_signed)

    return {
        "O_absG": rho_abs,
        "p_absG": p_abs,
        "O_signedG": rho_signed,
        "p_signedG": p_signed,
        "n_common_edges": int(len(T_vals))
    }


# ============================================================
# Action candidate
# ============================================================

def action_components(W, W_ref, args, d_star_value):
    geom, A, mask, D, ds = edge_geometry(
        W,
        args.k,
        args.tau_min,
        args.tau_max,
        args.tau_points
    )

    S_spec = (ds - d_star_value) ** 2

    if len(geom) > 0:
        S_geom = float(np.mean(geom["G_proxy"].values ** 2))
        S_loc = float(np.mean(geom["W"].values * geom["d_ent"].values ** 2))
    else:
        S_geom = 1e6
        S_loc = 1e6

    W_sum = float(np.sum(np.triu(W, 1)))
    W_ref_sum = float(np.sum(np.triu(W_ref, 1)))
    S_norm = ((W_sum - W_ref_sum) / max(W_ref_sum, 1e-15)) ** 2

    deltaW = W - W_ref
    S_exc_mass = float(np.mean(deltaW ** 2))

    # smoothness on excitation density
    rho = excitation_density(deltaW)
    S_smooth = 0.0
    n_terms = 0

    for _, row in geom.iterrows():
        i = int(row["i"])
        j = int(row["j"])
        wij = row["W"]
        S_smooth += wij * (rho[i] - rho[j]) ** 2
        n_terms += 1

    S_smooth = float(S_smooth / max(n_terms, 1))

    S_exc = S_exc_mass + args.eta_grad * S_smooth

    S_total = (
        args.alpha * S_spec
        + args.beta * S_geom
        + args.gamma * S_loc
        + args.mu * S_norm
        + args.m_exc * S_exc
    )

    return {
        "S_total": float(S_total),
        "S_spec": float(S_spec),
        "S_geom": float(S_geom),
        "S_loc": float(S_loc),
        "S_norm": float(S_norm),
        "S_exc": float(S_exc),
        "d_s": float(ds),
        "W_sum": float(W_sum),
        "n_edges": int(len(geom))
    }


def finite_difference_gradient(W, W_ref, args, d_star_value):
    W = clean_matrix(W)
    n = W.shape[0]
    h = args.fd_eps
    grad = np.zeros_like(W)

    S0 = action_components(W, W_ref, args, d_star_value)["S_total"]

    for i in range(n):
        for j in range(i + 1, n):
            if args.grad_positive_only and W[i, j] <= 0:
                continue

            Wp = W.copy()
            Wp[i, j] += h
            Wp[j, i] += h
            Wp = clean_matrix(Wp)

            Sp = action_components(Wp, W_ref, args, d_star_value)["S_total"]
            gij = (Sp - S0) / h

            grad[i, j] = gij
            grad[j, i] = gij

    return grad


def relax(W_init, W_ref, args, d_star_value, center, D_ref):
    W = clean_matrix(W_init.copy())
    rows = []

    for step in range(args.steps + 1):
        comps = action_components(W, W_ref, args, d_star_value)
        obs = matter_observable(W_ref, W, args.k, args.tau_min, args.tau_max, args.tau_points)

        deltaW = W - W_ref

        row = {
            "step": step,
            **comps,
            **obs,
            "M_exc": excitation_mass(deltaW),
            "R_exc": excitation_radius(deltaW, D_ref, center),
            "K_exc": excitation_kurtosis(deltaW, D_ref, center)
        }

        rows.append(row)

        if step == args.steps:
            break

        grad = finite_difference_gradient(W, W_ref, args, d_star_value)

        gnorm = np.linalg.norm(grad)
        if args.grad_normalize and gnorm > 1e-15:
            grad = grad / gnorm

        eta_try = args.eta_relax
        S_old = comps["S_total"]
        accepted = False

        for _ in range(12 if args.line_search else 1):
            W_try = clean_matrix(W - eta_try * grad)

            if args.clip_max > 0:
                W_try = np.minimum(W_try, args.clip_max)
                W_try = clean_matrix(W_try)

            S_try = action_components(W_try, W_ref, args, d_star_value)["S_total"]

            if (not args.line_search) or (S_try <= S_old):
                W = W_try
                accepted = True
                break

            eta_try *= 0.5

        if not accepted:
            break

    return W, pd.DataFrame(rows)


# ============================================================
# Crossings and phases
# ============================================================

def estimate_crossing(df, xcol, ycol):
    sub = df[[xcol, ycol]].dropna().sort_values(xcol)

    vals = sub.values
    candidates = []

    for a, b in zip(vals[:-1], vals[1:]):
        x1, y1 = a
        x2, y2 = b

        if y1 == 0:
            candidates.append(x1)
        elif y1 * y2 < 0:
            xc = x1 + (0.0 - y1) * (x2 - x1) / (y2 - y1)
            candidates.append(xc)

    if not candidates:
        return None

    return float(candidates[0])


def phase_label(o_init, o_final):
    si = "+" if o_init > 0 else "-"
    sf = "+" if o_final > 0 else "-"

    mapping = {
        "++": "source_stable",
        "+-": "source_decohered",
        "-+": "source_emergent_after_relaxation",
        "--": "non_source"
    }

    return mapping[si + sf]


# ============================================================
# Plot
# ============================================================

def make_plots(df, output_dir):
    ensure_dir(output_dir)

    plt.figure(figsize=(9, 5.5))
    plt.plot(df["lambda"], df["O_init"], "o-", linewidth=2, label=r"$O_{\rm init}$")
    plt.plot(df["lambda"], df["O_final"], "s-", linewidth=2, label=r"$O_{\rm final}$")
    plt.axhline(0.0, linewidth=1)

    c_init = estimate_crossing(df, "lambda", "O_init")
    c_final = estimate_crossing(df, "lambda", "O_final")

    if c_init is not None:
        plt.axvline(c_init, linestyle="--", linewidth=1.5, label=rf"cross init $\lambda\simeq{c_init:.3f}$")

    if c_final is not None:
        plt.axvline(c_final, linestyle=":", linewidth=1.8, label=rf"cross final $\lambda\simeq{c_final:.3f}$")

    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\mathcal{O}_{\rm matter}$")
    plt.title("Paper 7 — Transformation de phase sous relaxation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_phase_scan_lambda.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["lambda"], df["Delta_O"], "o-", linewidth=2)
    plt.axhline(0.0, linewidth=1)
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$\Delta O = O_{\rm final}-O_{\rm init}$")
    plt.title("Variation de cohérence géométrique")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_delta_O_lambda.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["lambda"], df["d_s_final"], "o-", label=r"$d_s^{final}$")
    plt.plot(df["lambda"], df["M_final"], "s-", label=r"$M_{\rm exc}^{final}$")
    plt.plot(df["lambda"], df["R_final"], "^-", label=r"$R_{\rm exc}^{final}$")
    plt.xlabel(r"$\lambda$")
    plt.title("Observables internes après relaxation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fig_internal_observables_lambda.png"), dpi=250)
    plt.close()


# ============================================================
# Main scan lambda
# ============================================================

def run_scan_lambda(args):
    ensure_dir(args.output_dir)

    files = sorted(glob.glob(os.path.join(args.mi_dir, "*.csv")))

    items = []
    for f in files:
        lam = parse_lambda(f)
        if lam is not None:
            items.append((lam, f))

    items = sorted(items, key=lambda x: x[0])

    if not items:
        raise FileNotFoundError(f"Aucun fichier MI avec lambda trouvé dans {args.mi_dir}")

    rows = []

    print("=" * 90)
    print("Paper 7 — Phase scan sous relaxation de S_BuP")
    print("=" * 90)
    print(f"MI dir     : {args.mi_dir}")
    print(f"count      : {len(items)}")
    print(f"steps      : {args.steps}")
    print(f"k          : {args.k}")
    print(f"d_star     : {args.d_star}")
    print(f"output     : {args.output_dir}")
    print("-" * 90)

    for lam, path in items:
        W0 = load_matrix_csv(path)

        geom0, A0, mask0, D0, ds0 = edge_geometry(
            W0,
            args.k,
            args.tau_min,
            args.tau_max,
            args.tau_points
        )

        if args.d_star == "auto":
            d_star_value = ds0
        else:
            d_star_value = float(args.d_star)

        W_init, deltaW0, center = inject_excitation(
            W0,
            D0,
            amp=args.amp,
            sigma=args.sigma,
            center=args.center
        )

        obs_init = matter_observable(
            W0,
            W_init,
            args.k,
            args.tau_min,
            args.tau_max,
            args.tau_points
        )

        W_final, hist = relax(
            W_init,
            W0,
            args,
            d_star_value=d_star_value,
            center=center,
            D_ref=D0
        )

        obs_final = matter_observable(
            W0,
            W_final,
            args.k,
            args.tau_min,
            args.tau_max,
            args.tau_points
        )

        deltaW_final = W_final - W0

        out_sub = os.path.join(args.output_dir, f"lambda_{lam:.4f}")
        ensure_dir(out_sub)
        hist.to_csv(os.path.join(out_sub, "history.csv"), index=False)
        np.savetxt(os.path.join(out_sub, "W_final.csv"), W_final, delimiter=",")

        o_init = obs_init["O_absG"]
        o_final = obs_final["O_absG"]

        row = {
            "lambda": lam,
            "file": os.path.basename(path),
            "O_init": o_init,
            "O_final": o_final,
            "Delta_O": o_final - o_init,
            "O_signed_init": obs_init["O_signedG"],
            "O_signed_final": obs_final["O_signedG"],
            "phase": phase_label(o_init, o_final),
            "d_s_initial": ds0,
            "d_s_final": hist["d_s"].iloc[-1],
            "S_initial": hist["S_total"].iloc[0],
            "S_final": hist["S_total"].iloc[-1],
            "Delta_S": hist["S_total"].iloc[-1] - hist["S_total"].iloc[0],
            "M_initial": excitation_mass(deltaW0),
            "M_final": excitation_mass(deltaW_final),
            "R_initial": excitation_radius(deltaW0, D0, center),
            "R_final": excitation_radius(deltaW_final, D0, center),
            "K_initial": excitation_kurtosis(deltaW0, D0, center),
            "K_final": excitation_kurtosis(deltaW_final, D0, center),
            "d_star_used": d_star_value
        }

        rows.append(row)

        print(
            f"lam={lam:.4f} "
            f"O_init={o_init:+.4f} "
            f"O_final={o_final:+.4f} "
            f"Delta={row['Delta_O']:+.4f} "
            f"phase={row['phase']} "
            f"dS={row['Delta_S']:+.3g}"
        )

    df = pd.DataFrame(rows).sort_values("lambda")
    df.to_csv(os.path.join(args.output_dir, "paper7_phase_scan_summary.csv"), index=False)

    cross_init = estimate_crossing(df, "lambda", "O_init")
    cross_final = estimate_crossing(df, "lambda", "O_final")

    summary = {
        "title": "Paper 7 phase scan under S_BuP relaxation",
        "description": "This scan tracks how localized BuP excitations change phase under relaxation of a candidate action.",
        "observable": "O_matter = Spearman(T_matter_proxy, |delta G_proxy|)",
        "interpretation": {
            "++": "source remains source after relaxation",
            "+-": "source decoheres or changes internal phase",
            "-+": "source-like behavior emerges after relaxation",
            "--": "non-source regime"
        },
        "config": vars(args),
        "crossing_init": cross_init,
        "crossing_final": cross_final,
        "n_points": int(len(df))
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    make_plots(df, args.output_dir)

    print("-" * 90)
    print("RÉSUMÉ")
    print("-" * 90)
    print(df[["lambda", "O_init", "O_final", "Delta_O", "phase", "d_s_final", "Delta_S"]].to_string(index=False))
    print("-" * 90)
    print(f"crossing init  : {cross_init}")
    print(f"crossing final : {cross_final}")
    print(f"output         : {args.output_dir}")
    print("DONE")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mi-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results_paper7_phase_scan_v1")

    parser.add_argument("--k", type=int, default=5)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=50.0)
    parser.add_argument("--tau-points", type=int, default=20)

    parser.add_argument("--d-star", type=str, default="auto",
                        help="'auto' or a float value, e.g. 0.9")

    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=0.1)
    parser.add_argument("--gamma", type=float, default=0.01)
    parser.add_argument("--mu", type=float, default=0.0)
    parser.add_argument("--m-exc", type=float, default=0.05)
    parser.add_argument("--eta-grad", type=float, default=0.1)

    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--eta-relax", type=float, default=0.001)
    parser.add_argument("--fd-eps", type=float, default=1e-4)
    parser.add_argument("--grad-normalize", action="store_true", default=True)
    parser.add_argument("--grad-positive-only", action="store_true")

    parser.add_argument("--amp", type=float, default=0.15)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--center", type=int, default=-1)

    parser.add_argument("--line-search", action="store_true")
    parser.add_argument("--clip-max", type=float, default=0.0)

    args = parser.parse_args()
    run_scan_lambda(args)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()

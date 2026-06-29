#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP Paper 14 — graph invariants + d_s/d_w + beta regression v2

Purpose
-------
Compute exact graph-level invariants from the mutual-information graph,
including:
    - rho_graph
    - mean degree
    - Laplacian algebraic connectivity lambda_2
    - spectral gap normalized
    - clustering coefficient
    - effective diameter proxy
    - spectral dimension d_s from heat-kernel scaling
    - walk dimension proxy d_w from diffusion MSD scaling on the graph

Then test whether:

    beta_mod = F(rho_graph, <k>, lambda_2, d_s, d_w, Delta_3, ...)

improves leave-one-topology-out generalization.

This script recomputes the quantum states, MI matrices, candidate L_A fits,
SFF features, graph invariants, and regression in one pipeline.

Recommended N=9 run
-------------------
python3 papers/paper14_modular_spectral_action/scripts/bup_paper14_ds_dw_beta_regression_v2.py \
    --topologies chain grid er \
    --n 9 \
    --subsystem-size 4 \
    --seeds 0 1 2 3 4 \
    --eta-zz-list 0.5 0.7 1.0 \
    --hx-list 0.4 0.6 0.8 \
    --hz-random-list 0.2 0.4 0.6 \
    --candidate schur \
    --min-r2-fit 0.90 \
    --rmt-min 0.45 \
    --rmt-max 0.60 \
    --j-disorder 0.2 \
    --output-dir papers/paper14_modular_spectral_action/results/ds_dw_beta_regression_v2_N9_schur_RMT

Also test induced:
------------------
python3 papers/paper14_modular_spectral_action/scripts/bup_paper14_ds_dw_beta_regression_v2.py \
    --topologies chain grid er \
    --n 9 \
    --subsystem-size 4 \
    --seeds 0 1 2 3 4 \
    --eta-zz-list 0.5 0.7 1.0 \
    --hx-list 0.4 0.6 0.8 \
    --hz-random-list 0.2 0.4 0.6 \
    --candidate induced \
    --min-r2-fit 0.90 \
    --rmt-min 0.45 \
    --rmt-max 0.60 \
    --j-disorder 0.2 \
    --output-dir papers/paper14_modular_spectral_action/results/ds_dw_beta_regression_v2_N9_induced_RMT
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import eigsh
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import KFold, LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.inspection import permutation_importance


EPS = 1e-12


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="BuP Paper 14 ds/dw graph invariant beta regression v2")

    p.add_argument("--topologies", nargs="+", default=["chain", "grid", "er"],
                   choices=["chain", "grid", "er"])
    p.add_argument("--n", type=int, default=9)
    p.add_argument("--subsystem-size", type=int, default=4)
    p.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])

    p.add_argument("--j", type=float, default=1.0)
    p.add_argument("--eta-zz-list", nargs="+", type=float, default=[0.7])
    p.add_argument("--hx-list", nargs="+", type=float, default=[0.4])
    p.add_argument("--hz", type=float, default=0.0)
    p.add_argument("--hz-random-list", nargs="+", type=float, default=[0.6])
    p.add_argument("--j-disorder", type=float, default=0.2)
    p.add_argument("--er-p", type=float, default=0.35)

    p.add_argument("--candidate", default="schur",
                   choices=["induced", "schur", "mi_normalized", "input_induced", "all"])
    p.add_argument("--min-r2-fit", type=float, default=0.90)
    p.add_argument("--rmt-min", type=float, default=None)
    p.add_argument("--rmt-max", type=float, default=None)

    p.add_argument("--max-qubits", type=int, default=16)

    p.add_argument("--sff-tmax", type=float, default=120.0)
    p.add_argument("--sff-nt", type=int, default=1024)
    p.add_argument("--plateau-frac", type=float, default=0.66)

    p.add_argument("--heat-t-min", type=float, default=0.05)
    p.add_argument("--heat-t-max", type=float, default=20.0)
    p.add_argument("--heat-n", type=int, default=50)
    p.add_argument("--walk-t-min", type=float, default=0.05)
    p.add_argument("--walk-t-max", type=float, default=20.0)
    p.add_argument("--walk-n", type=int, default=50)

    p.add_argument("--output-dir", required=True)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ---------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------

def make_edges(topology, n, seed, er_p):
    rng = np.random.default_rng(seed)
    if topology == "chain":
        return [(i, i + 1) for i in range(n - 1)]
    if topology == "grid":
        side = int(round(math.sqrt(n)))
        if side * side != n:
            return [(i, i + 1) for i in range(n - 1)]
        edges = []
        for r in range(side):
            for c in range(side):
                i = r * side + c
                if c + 1 < side:
                    edges.append((i, r * side + c + 1))
                if r + 1 < side:
                    edges.append((i, (r + 1) * side + c))
        return edges
    if topology == "er":
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < er_p:
                    edges.append((i, j))
        s = {tuple(sorted(e)) for e in edges}
        for i in range(n - 1):
            s.add((i, i + 1))
        return sorted(s)
    raise ValueError(topology)


def weighted_adjacency(n, edges, rng, disorder):
    A = np.zeros((n, n), dtype=float)
    for i, j in edges:
        w = 1.0 + disorder * rng.normal()
        w = max(w, 0.05)
        A[i, j] = A[j, i] = w
    return A


def laplacian(W):
    W = np.array(W, dtype=float, copy=True)
    np.fill_diagonal(W, 0.0)
    return np.diag(W.sum(axis=1)) - W


def normalized_laplacian(W):
    W = np.array(W, dtype=float, copy=True)
    np.fill_diagonal(W, 0.0)
    deg = W.sum(axis=1)
    inv = np.zeros_like(deg)
    m = deg > EPS
    inv[m] = 1.0 / np.sqrt(deg[m])
    return np.eye(W.shape[0]) - inv[:, None] * W * inv[None, :]


# ---------------------------------------------------------------------
# Sparse quantum operators
# ---------------------------------------------------------------------

def paulis():
    X = csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float))
    Z = csr_matrix(np.array([[1.0, 0.0], [0.0, -1.0]], dtype=float))
    I = identity(2, format="csr", dtype=float)
    return X, Z, I


def one_site_op(op, site, n):
    _, _, I = paulis()
    out = None
    for k in range(n):
        f = op if k == site else I
        out = f if out is None else kron(out, f, format="csr")
    return out


def two_site_op(op_a, op_b, i, j, n):
    _, _, I = paulis()
    out = None
    for k in range(n):
        if k == i:
            f = op_a
        elif k == j:
            f = op_b
        else:
            f = I
        out = f if out is None else kron(out, f, format="csr")
    return out


def precompute_ops(n):
    X, Z, _ = paulis()
    x_ops = [one_site_op(X, i, n) for i in range(n)]
    z_ops = [one_site_op(Z, i, n) for i in range(n)]
    xx_ops, zz_ops = {}, {}
    for i in range(n):
        for j in range(i + 1, n):
            xx_ops[(i, j)] = two_site_op(X, X, i, j, n)
            zz_ops[(i, j)] = two_site_op(Z, Z, i, j, n)
    return x_ops, z_ops, xx_ops, zz_ops


def build_hamiltonian(n, A_weight, j, eta_zz, hx, hz, hz_random_amp, rng, x_ops, z_ops, xx_ops, zz_ops):
    H = csr_matrix((2 ** n, 2 ** n), dtype=float)
    hz_random_values = hz_random_amp * rng.normal(size=n)

    for i in range(n):
        H += -hx * x_ops[i]
        H += -(hz + hz_random_values[i]) * z_ops[i]

    for i in range(n):
        for k in range(i + 1, n):
            if A_weight[i, k] > 0:
                Jij = j * A_weight[i, k]
                H += -Jij * xx_ops[(i, k)]
                H += -eta_zz * Jij * zz_ops[(i, k)]
    return H, hz_random_values


def ground_state(H):
    vals, vecs = eigsh(H, k=1, which="SA")
    psi = np.asarray(vecs[:, 0], dtype=complex)
    psi /= np.linalg.norm(psi)
    return float(vals[0]), psi


# ---------------------------------------------------------------------
# Quantum reductions
# ---------------------------------------------------------------------

def rdm(psi, keep, n):
    keep = list(keep)
    trace = [i for i in range(n) if i not in keep]
    T = psi.reshape([2] * n)
    P = np.transpose(T, axes=keep + trace)
    M = P.reshape(2 ** len(keep), 2 ** (n - len(keep)))
    rho = M @ M.conj().T
    return 0.5 * (rho + rho.conj().T)


def entropy(rho):
    vals = np.linalg.eigvalsh(rho).real
    vals = vals[vals > EPS]
    return float(-np.sum(vals * np.log(vals))) if len(vals) else 0.0


def mutual_information(psi, n):
    S1 = np.array([entropy(rdm(psi, [i], n)) for i in range(n)])
    MI = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            sij = entropy(rdm(psi, [i, j], n))
            mij = S1[i] + S1[j] - sij
            if mij < 0 and abs(mij) < 1e-10:
                mij = 0.0
            MI[i, j] = MI[j, i] = mij
    return MI, S1


def modular_spectrum(psi, A_sites, n):
    rho_A = rdm(psi, A_sites, n)
    vals = np.linalg.eigvalsh(rho_A).real
    vals = np.maximum(vals, EPS)
    vals /= vals.sum()
    return np.sort(-np.log(vals))


# ---------------------------------------------------------------------
# Spectral/SFF features
# ---------------------------------------------------------------------

def gap_ratio(vals):
    vals = np.sort(np.asarray(vals, dtype=float))
    gaps = np.diff(vals)
    gaps = gaps[gaps > 1e-10]
    if len(gaps) < 2:
        return np.nan
    r = np.minimum(gaps[:-1], gaps[1:]) / np.maximum(gaps[:-1], gaps[1:])
    return float(np.mean(r))


def sff_curve(kappa, tmax=120.0, nt=1024):
    k = np.asarray(kappa, dtype=float)
    dA = len(k)
    x = k - k.mean()
    if x.std() > EPS:
        x = x / x.std()
    t = np.linspace(0.0, tmax, nt)
    g2 = np.abs(np.exp(-1j * np.outer(t, x)).sum(axis=1)) ** 2 / dA ** 2
    return t, g2


def delta3_proxy(kappa):
    x = np.sort(np.asarray(kappa, dtype=float))
    if len(x) < 5 or np.std(x) < EPS:
        return np.nan
    u = (x - x.min()) / (x.max() - x.min() + EPS)
    N = np.arange(1, len(u) + 1, dtype=float) / len(u)
    X = np.vstack([u, np.ones_like(u)]).T
    coef, _, _, _ = np.linalg.lstsq(X, N, rcond=None)
    fit = X @ coef
    return float(np.mean((N - fit) ** 2))


def sff_features(kappa, tmax, nt, plateau_frac):
    t, g2 = sff_curve(kappa, tmax=tmax, nt=nt)
    dA = len(kappa)
    start = int(plateau_frac * nt)
    plateau = float(np.mean(g2[start:]))
    C = float(dA * plateau)

    valid_start = max(1, int(0.02 * nt))
    dip_idx = valid_start + int(np.argmin(g2[valid_start:]))
    t_dip = float(t[dip_idx])
    g2_dip = float(g2[dip_idx])

    win = max(5, nt // 80)
    smooth = np.convolve(g2, np.ones(win) / win, mode="same")
    ramp_idx = None
    for idx in range(dip_idx + 1, nt):
        if smooth[idx] >= plateau:
            ramp_idx = idx
            break
    if ramp_idx is None:
        ramp_idx = nt - 1
    if ramp_idx <= dip_idx:
        ramp_idx = min(nt - 1, dip_idx + max(1, nt // 100))

    t_ramp = float(t[ramp_idx])
    g2_ramp = float(g2[ramp_idx])
    if t_ramp <= t_dip:
        t_ramp = float(t[-1])
        g2_ramp = float(g2[-1])

    slope = float((g2_ramp - g2_dip) / max(t_ramp - t_dip, EPS))

    return dict(
        g2_plateau=plateau,
        C_modular=C,
        t_dip=t_dip,
        g2_dip=g2_dip,
        t_ramp=t_ramp,
        g2_ramp=g2_ramp,
        slope_ramp=slope,
        log_t_dip=float(np.log(t_dip + EPS)),
        log_t_ramp=float(np.log(t_ramp + EPS)),
        delta3_proxy=delta3_proxy(kappa),
    )


# ---------------------------------------------------------------------
# Graph invariants, d_s, d_w
# ---------------------------------------------------------------------

def graph_shortest_paths_from_weights(W):
    n = W.shape[0]
    # Convert stronger weight to shorter distance.
    dist = np.full((n, n), np.inf, dtype=float)
    np.fill_diagonal(dist, 0.0)
    for i in range(n):
        for j in range(n):
            if i != j and W[i, j] > EPS:
                dist[i, j] = 1.0 / W[i, j]
    # Floyd-Warshall, fine for n <= 16.
    for k in range(n):
        dist = np.minimum(dist, dist[:, [k]] + dist[[k], :])
    return dist


def clustering_weighted_binary(W):
    A = (W > EPS).astype(float)
    np.fill_diagonal(A, 0.0)
    n = A.shape[0]
    vals = []
    for i in range(n):
        neigh = np.where(A[i] > 0)[0]
        k = len(neigh)
        if k < 2:
            vals.append(0.0)
            continue
        sub = A[np.ix_(neigh, neigh)]
        e = sub.sum() / 2.0
        vals.append(2.0 * e / (k * (k - 1)))
    return float(np.mean(vals))


def fit_slope_loglog(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > EPS) & (y > EPS)
    x, y = x[m], y[m]
    if len(x) < 5:
        return np.nan, np.nan
    lx, ly = np.log(x), np.log(y)
    # Use middle 60% to avoid UV/IR artifacts.
    order = np.argsort(lx)
    lx, ly = lx[order], ly[order]
    a = int(0.2 * len(lx))
    b = int(0.8 * len(lx))
    if b <= a + 2:
        a, b = 0, len(lx)
    coef = np.polyfit(lx[a:b], ly[a:b], 1)
    pred = np.polyval(coef, lx[a:b])
    ss_res = np.sum((ly[a:b] - pred) ** 2)
    ss_tot = np.sum((ly[a:b] - ly[a:b].mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > EPS else np.nan
    return float(coef[0]), float(r2)


def spectral_dimension_from_L(L, t_min, t_max, n_t):
    eig = np.linalg.eigvalsh(L)
    eig = np.maximum(eig, 0.0)
    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)
    Z = np.array([np.sum(np.exp(-tt * eig)) for tt in t]) / len(eig)
    slope, r2 = fit_slope_loglog(t, Z)
    ds = -2.0 * slope if np.isfinite(slope) else np.nan
    return float(ds), float(r2)


def walk_dimension_proxy_from_L_and_dist(L, dist, t_min, t_max, n_t):
    eig, U = np.linalg.eigh(L)
    eig = np.maximum(eig, 0.0)
    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)

    D2 = dist ** 2
    finite = np.isfinite(D2)
    if not np.all(finite):
        max_f = np.nanmax(D2[finite]) if finite.any() else 1.0
        D2 = np.where(finite, D2, max_f)

    msd = []
    for tt in t:
        K = (U * np.exp(-tt * eig)) @ U.T
        K = np.maximum(K, 0.0)
        row_sums = K.sum(axis=1, keepdims=True)
        P = K / np.maximum(row_sums, EPS)
        m = np.mean(np.sum(P * D2, axis=1))
        msd.append(m)
    msd = np.asarray(msd)

    slope, r2 = fit_slope_loglog(t, msd)
    # MSD ~ t^(2/dw), so d_w = 2/slope.
    dw = 2.0 / slope if np.isfinite(slope) and abs(slope) > EPS else np.nan
    return float(dw), float(r2)


def graph_invariants(W, args):
    W = np.array(W, dtype=float, copy=True)
    np.fill_diagonal(W, 0.0)
    n = W.shape[0]
    E = float(np.count_nonzero(np.triu(W > EPS, k=1)))
    rho = 2.0 * E / max(n * (n - 1), 1)
    deg = W.sum(axis=1)
    bin_deg = (W > EPS).sum(axis=1)

    L = laplacian(W)
    Ln = normalized_laplacian(W)

    eigL = np.linalg.eigvalsh(L)
    eigLn = np.linalg.eigvalsh(Ln)
    eigL_sorted = np.sort(np.maximum(eigL, 0.0))
    eigLn_sorted = np.sort(np.maximum(eigLn, 0.0))

    lambda2 = float(eigL_sorted[1]) if len(eigL_sorted) > 1 else np.nan
    norm_lambda2 = float(eigLn_sorted[1]) if len(eigLn_sorted) > 1 else np.nan

    dist = graph_shortest_paths_from_weights(W)
    finite = dist[np.isfinite(dist) & (dist > 0)]
    eff_diam = float(np.percentile(finite, 90)) if finite.size else np.nan
    mean_sp = float(np.mean(finite)) if finite.size else np.nan

    clustering = clustering_weighted_binary(W)

    ds, ds_r2 = spectral_dimension_from_L(Ln, args.heat_t_min, args.heat_t_max, args.heat_n)
    dw, dw_r2 = walk_dimension_proxy_from_L_and_dist(Ln, dist, args.walk_t_min, args.walk_t_max, args.walk_n)

    return dict(
        graph_edges=E,
        rho_graph=rho,
        mean_degree_weighted=float(np.mean(deg)),
        std_degree_weighted=float(np.std(deg)),
        mean_degree_binary=float(np.mean(bin_deg)),
        std_degree_binary=float(np.std(bin_deg)),
        lambda2=float(lambda2),
        norm_lambda2=float(norm_lambda2),
        laplacian_trace=float(np.trace(L)),
        clustering_binary=clustering,
        effective_diameter_90=eff_diam,
        mean_shortest_path=mean_sp,
        ds_graph=ds,
        ds_fit_r2=ds_r2,
        dw_graph=dw,
        dw_fit_r2=dw_r2,
        ds_over_dw=float(ds / dw) if np.isfinite(ds) and np.isfinite(dw) and abs(dw) > EPS else np.nan,
        alpha_eff_graph=float(2.0 * ds / dw + dw - 4.0) if np.isfinite(ds) and np.isfinite(dw) and abs(dw) > EPS else np.nan,
    )


# ---------------------------------------------------------------------
# Candidate Laplacians
# ---------------------------------------------------------------------

def schur_laplacian(L, A_sites):
    n = L.shape[0]
    A = list(A_sites)
    B = [i for i in range(n) if i not in A]
    LAA = L[np.ix_(A, A)]
    if len(B) == 0:
        return LAA
    LAB = L[np.ix_(A, B)]
    LBA = L[np.ix_(B, A)]
    LBB = L[np.ix_(B, B)]
    out = LAA - LAB @ np.linalg.pinv(LBB, rcond=1e-10) @ LBA
    return 0.5 * (out + out.T)


def candidate_spectra(MI, A_weight, A_sites):
    A = list(A_sites)
    L_global_mi = laplacian(MI)
    W_A_mi = MI[np.ix_(A, A)]
    W_A_input = A_weight[np.ix_(A, A)]
    return {
        "induced": np.linalg.eigvalsh(laplacian(W_A_mi)),
        "schur": np.linalg.eigvalsh(schur_laplacian(L_global_mi, A)),
        "input_induced": np.linalg.eigvalsh(laplacian(W_A_input)),
        "mi_normalized": np.linalg.eigvalsh(normalized_laplacian(W_A_mi)),
    }


def quantile_match(lam, kap):
    lam = np.sort(np.real(np.asarray(lam)[np.isfinite(lam)]))
    kap = np.sort(np.real(np.asarray(kap)[np.isfinite(kap)]))
    if len(lam) < 2 or len(kap) < 4:
        return np.array([]), np.array([])
    qk = np.linspace(0.0, 1.0, len(kap))
    ql = np.linspace(0.0, 1.0, len(lam))
    x = np.interp(qk, ql, lam)
    x = x - x.min() + 1e-6
    y = kap - kap.min() + 1e-6
    return x, y


def f_pos(x, A, beta, B):
    return A * np.power(x + 1e-9, beta) + B


def f_neg(x, A, beta, B):
    return A * np.power(x + 1e-9, -beta) + B


def r2_manual(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > EPS else np.nan


def fit_family(x, y, family):
    func = f_pos if family == "pos" else f_neg
    try:
        p0 = [max(y.max() - y.min(), 1e-3), 1.0, y.min()]
        popt, _ = curve_fit(
            func, x, y,
            p0=p0,
            bounds=([0.0, 0.0, -np.inf], [np.inf, 10.0, np.inf]),
            maxfev=20000,
        )
        yhat = func(x, *popt)
        return dict(fit_ok=1, r2=float(r2_manual(y, yhat)),
                    A_fit=float(popt[0]), beta=float(popt[1]), B_fit=float(popt[2]))
    except Exception:
        return dict(fit_ok=0, r2=np.nan, A_fit=np.nan, beta=np.nan, B_fit=np.nan)


def evaluate_candidate(lam, kap):
    x, y = quantile_match(lam, kap)
    if len(x) == 0:
        return []
    sp = spearmanr(x, y).statistic
    sp = float(sp) if np.isfinite(sp) else np.nan
    out = []
    for fam in ["pos", "neg"]:
        fit = fit_family(x, y, fam)
        out.append(dict(family=fam, spearman_lambda_kappa=sp, **fit))
    return out


# ---------------------------------------------------------------------
# Regression
# ---------------------------------------------------------------------

def model_dict(random_state):
    return {
        "linear": LinearRegression(),
        "ridge_1": Ridge(alpha=1.0),
        "ridge_10": Ridge(alpha=10.0),
        "lasso_001": Lasso(alpha=0.01, max_iter=20000),
        "lasso_01": Lasso(alpha=0.1, max_iter=20000),
        "random_forest": RandomForestRegressor(n_estimators=500, min_samples_leaf=3, random_state=random_state),
        "gradient_boosting": GradientBoostingRegressor(n_estimators=300, learning_rate=0.035, max_depth=3, random_state=random_state),
    }


def make_preprocessor(numeric, categorical):
    try:
        enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        enc = OneHotEncoder(handle_unknown="ignore", sparse=False)
    return ColumnTransformer(
        transformers=[
            ("num", Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())]), numeric),
            ("cat", Pipeline([("imputer", SimpleImputer(strategy="most_frequent")), ("onehot", enc)]), categorical),
        ],
        remainder="drop",
    )


def metrics(y, pred):
    return dict(
        r2=float(r2_score(y, pred)),
        rmse=float(np.sqrt(mean_squared_error(y, pred))),
        mae=float(mean_absolute_error(y, pred)),
    )


def leave_group_eval(model, pre, X, y, groups):
    logo = LeaveOneGroupOut()
    pred = np.full_like(y, np.nan, dtype=float)
    rows = []
    for fold, (tr, te) in enumerate(logo.split(X, y, groups=groups)):
        pipe = Pipeline([("pre", pre), ("model", model)])
        pipe.fit(X.iloc[tr], y[tr])
        p = pipe.predict(X.iloc[te])
        pred[te] = p
        rows.append(dict(fold=fold, heldout_group=str(np.asarray(groups)[te][0]),
                         n_train=int(len(tr)), n_test=int(len(te)), **metrics(y[te], p)))
    return pred, pd.DataFrame(rows)


def evaluate_feature_set(df, numeric, categorical, set_name, out, random_state):
    X = df[numeric + categorical].copy()
    y = df["beta"].astype(float).values
    pre = make_preprocessor(numeric, categorical)
    rows = []
    kfold = KFold(n_splits=min(5, len(df)), shuffle=True, random_state=random_state)
    preds = {}
    fold_dir = out / f"folds_{set_name}"
    fold_dir.mkdir(exist_ok=True)

    for name, model in model_dict(random_state).items():
        pipe = Pipeline([("pre", pre), ("model", model)])
        pred = cross_val_predict(pipe, X, y, cv=kfold)
        rows.append(dict(feature_set=set_name, model=name, cv="kfold", n=len(df), **metrics(y, pred)))
        preds[(name, "kfold")] = pred

    if "topology" in df.columns and df["topology"].nunique() >= 2:
        groups = df["topology"].astype(str).values
        for name, model in model_dict(random_state).items():
            pred, tab = leave_group_eval(model, pre, X, y, groups)
            mask = np.isfinite(pred)
            rows.append(dict(feature_set=set_name, model=name, cv="leave_one_topology", n=int(mask.sum()), **metrics(y[mask], pred[mask])))
            tab.to_csv(fold_dir / f"leave_one_topology_{name}.csv", index=False)
            preds[(name, "leave_one_topology")] = pred

    if "regime" in df.columns and df["regime"].nunique() >= 2:
        groups = df["regime"].astype(str).values
        for name, model in model_dict(random_state).items():
            pred, tab = leave_group_eval(model, pre, X, y, groups)
            mask = np.isfinite(pred)
            rows.append(dict(feature_set=set_name, model=name, cv="leave_one_regime", n=int(mask.sum()), **metrics(y[mask], pred[mask])))
            tab.to_csv(fold_dir / f"leave_one_regime_{name}.csv", index=False)
            preds[(name, "leave_one_regime")] = pred

    res = pd.DataFrame(rows)

    # Best kfold feature importance.
    best = res[res["cv"] == "kfold"].sort_values("r2", ascending=False).iloc[0]
    best_model = model_dict(random_state)[best["model"]]
    pipe = Pipeline([("pre", pre), ("model", best_model)])
    pipe.fit(X, y)

    try:
        imp = permutation_importance(pipe, X, y, n_repeats=30, random_state=random_state, scoring="r2")
        imp_df = pd.DataFrame({
            "feature": numeric + categorical,
            "importance_mean": imp.importances_mean,
            "importance_std": imp.importances_std,
            "feature_set": set_name,
            "best_model": best["model"],
        }).sort_values("importance_mean", ascending=False)
    except Exception as exc:
        imp_df = pd.DataFrame({"error": [str(exc)], "feature_set": [set_name]})

    imp_df.to_csv(out / f"permutation_importance_{set_name}.csv", index=False)

    # Plots for best per cv.
    for cv in ["kfold", "leave_one_topology", "leave_one_regime"]:
        sub = res[res["cv"] == cv].sort_values("r2", ascending=False)
        if len(sub) == 0:
            continue
        row = sub.iloc[0]
        key = (row["model"], cv)
        if key not in preds:
            continue
        pred = preds[key]
        mask = np.isfinite(pred)
        if mask.sum() < 3:
            continue
        plt.figure(figsize=(6, 6))
        plt.scatter(y[mask], pred[mask], s=35, alpha=0.75)
        mn = min(y[mask].min(), pred[mask].min())
        mx = max(y[mask].max(), pred[mask].max())
        plt.plot([mn, mx], [mn, mx], "--")
        plt.xlabel("beta true")
        plt.ylabel("beta predicted")
        plt.title(f"{set_name} | {cv} | {row['model']} | R2={row['r2']:.3f}")
        plt.tight_layout()
        plt.savefig(out / f"fig_pred_{set_name}_{cv}_{row['model']}.png", dpi=220)
        plt.close()

    return res, imp_df


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()
    if args.n > args.max_qubits:
        raise ValueError(f"Requested n={args.n}, max-qubits={args.max_qubits}")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("\n=== BuP Paper 14 — d_s/d_w beta regression v2 ===")
    print(f"n={args.n}, dim={2**args.n}, |A|={args.subsystem_size}")
    print(f"candidate={args.candidate}, min_fit_r2={args.min_r2_fit}")
    print(f"RMT filter={args.rmt_min}, {args.rmt_max}")
    print(f"output={out}")

    print("\nPrecomputing sparse operators...")
    x_ops, z_ops, xx_ops, zz_ops = precompute_ops(args.n)

    rows = []
    t0 = time.time()
    param_grid = list(itertools.product(args.eta_zz_list, args.hx_list, args.hz_random_list))
    total = len(param_grid) * len(args.topologies) * len(args.seeds)
    counter = 0

    for eta_zz, hx, hz_random in param_grid:
        for topology in args.topologies:
            for seed in args.seeds:
                counter += 1
                rng = np.random.default_rng(seed)
                print(f"\n[{counter}/{total}] eta={eta_zz}, hx={hx}, hzr={hz_random}, topology={topology}, seed={seed}")

                try:
                    edges = make_edges(topology, args.n, seed, args.er_p)
                    A_weight = weighted_adjacency(args.n, edges, rng, args.j_disorder)
                    A_sites = sorted(rng.choice(np.arange(args.n), size=args.subsystem_size, replace=False).tolist())

                    H, hz_vals = build_hamiltonian(args.n, A_weight, args.j, eta_zz, hx, args.hz, hz_random, rng, x_ops, z_ops, xx_ops, zz_ops)
                    E0, psi = ground_state(H)
                    MI, S1 = mutual_information(psi, args.n)
                    kappa = modular_spectrum(psi, A_sites, args.n)
                    gr = gap_ratio(kappa)
                    sff = sff_features(kappa, args.sff_tmax, args.sff_nt, args.plateau_frac)

                    # Graph invariants on global MI graph.
                    inv_mi = graph_invariants(MI, args)
                    inv_input = graph_invariants(A_weight, args)

                    spectra = candidate_spectra(MI, A_weight, A_sites)
                    for cand, lam in spectra.items():
                        if args.candidate != "all" and cand != args.candidate:
                            continue
                        for erow in evaluate_candidate(lam, kappa):
                            row = dict(
                                n=args.n,
                                hilbert_dim=2**args.n,
                                subsystem_size=args.subsystem_size,
                                d_A=2**args.subsystem_size,
                                topology=topology,
                                seed=seed,
                                eta_zz=eta_zz,
                                hx=hx,
                                hz=args.hz,
                                hz_random=hz_random,
                                j_disorder=args.j_disorder,
                                j=args.j,
                                regime=f"{eta_zz:.6g}_{hx:.6g}_{hz_random:.6g}",
                                A_sites=",".join(map(str, A_sites)),
                                n_edges_input=len(edges),
                                ground_energy=float(E0),
                                mean_hz_random=float(np.mean(hz_vals)),
                                std_hz_random=float(np.std(hz_vals)),
                                mean_single_entropy=float(np.mean(S1)),
                                mean_MI=float(np.mean(MI[np.triu_indices(args.n, k=1)])),
                                max_MI=float(np.max(MI)),
                                K_gap_ratio=float(gr) if np.isfinite(gr) else np.nan,
                                candidate_L=cand,
                            )
                            row.update(sff)
                            # prefix invariants.
                            row.update({f"mi_{k}": v for k, v in inv_mi.items()})
                            row.update({f"input_{k}": v for k, v in inv_input.items()})
                            row.update(erow)
                            rows.append(row)

                    print(f"  E0={E0:.5f} meanMI={row['mean_MI']:.4f} gap={gr:.4f} ds={inv_mi['ds_graph']:.3f} dw={inv_mi['dw_graph']:.3f}")

                except Exception as exc:
                    print("  FAILED:", exc)
                    rows.append(dict(topology=topology, seed=seed, eta_zz=eta_zz, hx=hx, hz_random=hz_random,
                                     candidate_L=args.candidate, family="failed", fit_ok=0, error=str(exc)))

    df = pd.DataFrame(rows)
    df.to_csv(out / "paper14_v2_ds_dw_full_results.csv", index=False)

    # Filter for regression.
    reg = df[(df["fit_ok"] == 1) & (df["family"] == "pos")].copy()
    reg = reg[np.isfinite(reg["beta"]) & np.isfinite(reg["r2"])].copy()
    reg = reg[reg["r2"] >= args.min_r2_fit].copy()
    if args.rmt_min is not None:
        reg = reg[reg["K_gap_ratio"] >= args.rmt_min].copy()
    if args.rmt_max is not None:
        reg = reg[reg["K_gap_ratio"] <= args.rmt_max].copy()

    reg.to_csv(out / "paper14_v2_ds_dw_regression_dataset.csv", index=False)

    if len(reg) < 10:
        print(f"\nNot enough regression rows after filtering: {len(reg)}")
        return

    # Feature sets.
    sff_feature_cols = [
        "C_modular", "K_gap_ratio", "log_t_dip", "log_t_ramp",
        "slope_ramp", "delta3_proxy",
    ]

    graph_basic = [
        "mi_rho_graph", "mi_mean_degree_weighted", "mi_mean_degree_binary",
        "mi_lambda2", "mi_norm_lambda2", "mi_clustering_binary",
        "mi_effective_diameter_90", "mi_mean_shortest_path",
        "mi_laplacian_trace", "mean_MI", "max_MI", "mean_single_entropy",
        "n_edges_input",
    ]

    graph_ds_dw = [
        "mi_ds_graph", "mi_ds_fit_r2", "mi_dw_graph", "mi_dw_fit_r2",
        "mi_ds_over_dw", "mi_alpha_eff_graph",
    ]

    input_graph = [
        "input_rho_graph", "input_mean_degree_weighted", "input_mean_degree_binary",
        "input_lambda2", "input_norm_lambda2", "input_clustering_binary",
        "input_effective_diameter_90", "input_mean_shortest_path",
    ]

    params = ["eta_zz", "hx", "hz_random", "j_disorder", "n", "subsystem_size"]

    feature_sets = {
        "basic_no_topology": (sff_feature_cols + graph_basic + params, []),
        "dsdw_no_topology": (sff_feature_cols + graph_basic + graph_ds_dw + params, []),
        "dsdw_input_no_topology": (sff_feature_cols + graph_basic + graph_ds_dw + input_graph + params, []),
        "dsdw_with_topology": (sff_feature_cols + graph_basic + graph_ds_dw + params, ["topology"]),
    }

    all_res, all_imp = [], []
    for name, (num, cat) in feature_sets.items():
        num = [c for c in num if c in reg.columns]
        cat = [c for c in cat if c in reg.columns]
        print(f"\n--- Regression feature set: {name} ---")
        res, imp = evaluate_feature_set(reg, num, cat, name, out, args.random_state)
        print(res.sort_values(["cv", "r2"], ascending=[True, False]).to_string(index=False))
        all_res.append(res)
        all_imp.append(imp)

    results = pd.concat(all_res, ignore_index=True)
    importance = pd.concat(all_imp, ignore_index=True)
    results.to_csv(out / "ds_dw_beta_regression_results.csv", index=False)
    importance.to_csv(out / "ds_dw_permutation_importance_all.csv", index=False)

    best = results.sort_values("r2", ascending=False).groupby(["feature_set", "cv"]).head(1)
    best = best.sort_values(["cv", "r2"], ascending=[True, False])
    best.to_csv(out / "best_by_feature_set_and_cv.csv", index=False)

    summary = dict(
        experiment="BuP Paper 14 graph invariants d_s/d_w beta regression v2",
        n=args.n,
        candidate=args.candidate,
        min_r2_fit=args.min_r2_fit,
        rmt_min=args.rmt_min,
        rmt_max=args.rmt_max,
        n_full_rows=int(len(df)),
        n_regression_rows=int(len(reg)),
        runtime_seconds=float(time.time() - t0),
        best_by_cv=best.to_dict(orient="records"),
    )
    with open(out / "ds_dw_beta_regression_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 90)
    print("BEST BY FEATURE SET AND CV")
    print("=" * 90)
    print(best.to_string(index=False))

    print("\nFiles written:")
    for p in sorted(out.iterdir()):
        if p.is_file():
            print("  -", p)

    print("\nDone.")


if __name__ == "__main__":
    main()

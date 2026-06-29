#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP Paper 14 — modular/spectral + SFF features v3

Goal
----
Move beyond the failed simple hypothesis C -> beta.

V3 tests the richer hypothesis:

    beta = F(C, gap_ratio, t_dip, t_ramp, slope_ramp, delta3, topology, ...)

Pipeline
--------
For each quantum realization:
1. Build a chaotic-ish XX+ZZ Hamiltonian with disorder.
2. Compute the ground state.
3. Compute K_A = -log rho_A and its spectrum {kappa_n}.
4. Compute the MI graph and candidate restricted Laplacians L_A.
5. Fit Spec(K_A) ≃ A Spec(L_A)^beta + B.
6. Extract SFF features from {kappa_n}:
   - C = d_A * g2_plateau
   - t_dip
   - t_ramp
   - slope_ramp
   - delta3_proxy
   - gap_ratio
7. Write one CSV suitable for multivariate regression.

Recommended N=9 scan
--------------------
python3 papers/paper14_modular_spectral_action/scripts/bup_paper14_sff_features_v3.py \
    --topologies chain grid er \
    --n 9 \
    --subsystem-size 4 \
    --seeds 0 1 2 3 4 \
    --eta-zz-list 0.5 0.7 1.0 \
    --hx-list 0.4 0.6 0.8 \
    --hz-random-list 0.2 0.4 0.6 \
    --j-disorder 0.2 \
    --output-dir papers/paper14_modular_spectral_action/results/sff_features_v3_N9

Recommended targeted N=16
-------------------------
python3 papers/paper14_modular_spectral_action/scripts/bup_paper14_sff_features_v3.py \
    --topologies chain grid er \
    --n 16 \
    --subsystem-size 4 \
    --seeds 0 1 2 \
    --eta-zz-list 0.7 \
    --hx-list 0.4 \
    --hz-random-list 0.6 \
    --j-disorder 0.2 \
    --max-qubits 16 \
    --output-dir papers/paper14_modular_spectral_action/results/sff_features_v3_N16_optimal
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


EPS = 1e-12


# ---------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="BuP Paper 14 SFF features v3")

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

    p.add_argument("--max-qubits", type=int, default=16)
    p.add_argument("--sff-tmax", type=float, default=120.0)
    p.add_argument("--sff-nt", type=int, default=1024)
    p.add_argument("--plateau-frac", type=float, default=0.66)
    p.add_argument("--output-dir", type=str, required=True)
    p.add_argument("--no-plots", action="store_true")

    return p.parse_args()


# ---------------------------------------------------------------------
# Graphs
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
        edge_set = {tuple(sorted(e)) for e in edges}
        for i in range(n - 1):
            edge_set.add((i, i + 1))
        return sorted(edge_set)

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
# Sparse operators
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
    xx_ops = {}
    zz_ops = {}
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


def gap_ratio(vals):
    vals = np.sort(np.asarray(vals, dtype=float))
    gaps = np.diff(vals)
    gaps = gaps[gaps > 1e-10]
    if len(gaps) < 2:
        return np.nan
    r = np.minimum(gaps[:-1], gaps[1:]) / np.maximum(gaps[:-1], gaps[1:])
    return float(np.mean(r))


# ---------------------------------------------------------------------
# SFF features
# ---------------------------------------------------------------------

def sff_curve(kappa, tmax=120.0, nt=1024):
    k = np.asarray(kappa, dtype=float)
    dA = len(k)
    x = k - k.mean()
    if x.std() > EPS:
        x = x / x.std()
    t = np.linspace(0.0, tmax, nt)
    g2 = np.abs(np.exp(-1j * np.outer(t, x)).sum(axis=1)) ** 2 / dA ** 2
    return t, g2


def spectral_rigidity_delta3_proxy(kappa):
    """
    Lightweight proxy for Delta_3 rigidity.

    True Delta_3 requires unfolding and window averaging.
    Here we use the mean squared deviation of the unfolded staircase from
    the best linear fit. It is a stable diagnostic feature, not a final
    analytic Delta_3 measurement.
    """
    x = np.sort(np.asarray(kappa, dtype=float))
    if len(x) < 5 or np.std(x) < EPS:
        return np.nan

    # crude unfolding: normalize spectrum to [0, 1]
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

    # Ignore t=0 for dip.
    valid_start = max(1, int(0.02 * nt))
    dip_idx = valid_start + int(np.argmin(g2[valid_start:]))
    t_dip = float(t[dip_idx])
    g2_dip = float(g2[dip_idx])

    # ramp reaches plateau: first time after dip where smoothed g2 >= plateau
    # use simple moving average
    win = max(5, nt // 80)
    kernel = np.ones(win) / win
    smooth = np.convolve(g2, kernel, mode="same")

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

    dt = max(t_ramp - t_dip, EPS)
    slope_ramp = float((g2_ramp - g2_dip) / dt)

    delta3 = spectral_rigidity_delta3_proxy(kappa)

    return {
        "g2_plateau": plateau,
        "C_modular": C,
        "t_dip": t_dip,
        "g2_dip": g2_dip,
        "t_ramp": t_ramp,
        "g2_ramp": g2_ramp,
        "slope_ramp": slope_ramp,
        "log_t_dip": float(np.log(t_dip + EPS)),
        "log_t_ramp": float(np.log(t_ramp + EPS)),
        "delta3_proxy": delta3,
    }


# ---------------------------------------------------------------------
# Candidate L_A spectra
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


# ---------------------------------------------------------------------
# Fitting K_A spectrum to L_A spectrum
# ---------------------------------------------------------------------

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


def r2_score_manual(y, yhat):
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > EPS else np.nan


def fit_family(x, y, family):
    func = f_pos if family == "pos" else f_neg
    try:
        p0 = [max(y.max() - y.min(), 1e-3), 1.0, y.min()]
        popt, _ = curve_fit(
            func, x, y, p0=p0,
            bounds=([0.0, 0.0, -np.inf], [np.inf, 10.0, np.inf]),
            maxfev=20000,
        )
        yhat = func(x, *popt)
        return {
            "fit_ok": 1,
            "r2": float(r2_score_manual(y, yhat)),
            "A_fit": float(popt[0]),
            "beta": float(popt[1]),
            "B_fit": float(popt[2]),
        }
    except Exception:
        return {"fit_ok": 0, "r2": np.nan, "A_fit": np.nan, "beta": np.nan, "B_fit": np.nan}


def evaluate_candidate(lam, kap):
    x, y = quantile_match(lam, kap)
    if len(x) == 0:
        return []

    sp = spearmanr(x, y).statistic
    sp = float(sp) if np.isfinite(sp) else np.nan

    rows = []
    for family in ["pos", "neg"]:
        fit = fit_family(x, y, family)
        rows.append({
            "family": family,
            "spearman_lambda_kappa": sp,
            **fit,
        })
    return rows


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

def plot_beta_features(df, out):
    valid = df[(df["fit_ok"] == 1) & (df["family"] == "pos") & (df["r2"] >= 0.9)].copy()
    if len(valid) == 0:
        return

    features = ["C_modular", "K_gap_ratio", "t_dip", "t_ramp", "slope_ramp", "delta3_proxy"]
    for feat in features:
        if feat not in valid.columns:
            continue
        plt.figure(figsize=(7, 5))
        for cand in sorted(valid["candidate_L"].unique()):
            sub = valid[valid["candidate_L"] == cand]
            plt.scatter(sub[feat], sub["beta"], s=40, alpha=0.75, label=cand)
        plt.xlabel(feat)
        plt.ylabel("beta")
        plt.title(f"beta versus {feat} — high quality rows")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out / f"fig_beta_vs_{feat}.png", dpi=220)
        plt.close()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()
    if args.n > args.max_qubits:
        raise ValueError(f"Requested n={args.n}, max-qubits={args.max_qubits}")
    if args.subsystem_size >= args.n:
        raise ValueError("--subsystem-size must be smaller than n")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("\n=== BuP Paper 14 — SFF features v3 ===")
    print(f"n qubits          : {args.n}")
    print(f"Hilbert dimension : {2 ** args.n}")
    print(f"|A|, d_A          : {args.subsystem_size}, {2 ** args.subsystem_size}")
    print(f"topologies        : {args.topologies}")
    print(f"seeds             : {args.seeds}")
    print(f"eta_zz_list       : {args.eta_zz_list}")
    print(f"hx_list           : {args.hx_list}")
    print(f"hz_random_list    : {args.hz_random_list}")
    print(f"output            : {out}")

    print("\nPrecomputing sparse operators...")
    x_ops, z_ops, xx_ops, zz_ops = precompute_ops(args.n)

    rows = []
    t0 = time.time()

    param_grid = list(itertools.product(args.eta_zz_list, args.hx_list, args.hz_random_list))
    total = len(args.topologies) * len(args.seeds) * len(param_grid)
    counter = 0

    for eta_zz, hx, hz_random in param_grid:
        for topology in args.topologies:
            for seed in args.seeds:
                counter += 1
                rng = np.random.default_rng(seed)
                print(f"\n[{counter}/{total}] eta={eta_zz}, hx={hx}, hzr={hz_random}, topology={topology}, seed={seed}")

                edges = make_edges(topology, args.n, seed, args.er_p)
                A_weight = weighted_adjacency(args.n, edges, rng, args.j_disorder)
                A_sites = sorted(
                    rng.choice(np.arange(args.n), size=args.subsystem_size, replace=False).tolist()
                )

                try:
                    H, hz_random_values = build_hamiltonian(
                        args.n, A_weight, args.j, eta_zz, hx, args.hz,
                        hz_random, rng, x_ops, z_ops, xx_ops, zz_ops
                    )
                    E0, psi0 = ground_state(H)
                    MI, S1 = mutual_information(psi0, args.n)
                    kappa = modular_spectrum(psi0, A_sites, args.n)

                    sff = sff_features(
                        kappa,
                        tmax=args.sff_tmax,
                        nt=args.sff_nt,
                        plateau_frac=args.plateau_frac,
                    )
                    gr = gap_ratio(kappa)

                    spectra = candidate_spectra(MI, A_weight, A_sites)

                    for cand_name, lam in spectra.items():
                        for erow in evaluate_candidate(lam, kappa):
                            row = {
                                "n": args.n,
                                "hilbert_dim": 2 ** args.n,
                                "subsystem_size": args.subsystem_size,
                                "d_A": 2 ** args.subsystem_size,
                                "topology": topology,
                                "seed": seed,
                                "eta_zz": eta_zz,
                                "hx": hx,
                                "hz": args.hz,
                                "hz_random": hz_random,
                                "j_disorder": args.j_disorder,
                                "j": args.j,
                                "A_sites": ",".join(map(str, A_sites)),
                                "n_edges_input": len(edges),
                                "ground_energy": float(E0),
                                "mean_hz_random": float(np.mean(hz_random_values)),
                                "std_hz_random": float(np.std(hz_random_values)),
                                "mean_single_entropy": float(np.mean(S1)),
                                "mean_MI": float(np.mean(MI[np.triu_indices(args.n, k=1)])),
                                "max_MI": float(np.max(MI)),
                                "K_gap_ratio": float(gr) if np.isfinite(gr) else np.nan,
                                "candidate_L": cand_name,
                            }
                            row.update(sff)
                            row.update(erow)
                            rows.append(row)

                    print(f"  E0={E0:.6f} C={sff['C_modular']:.4f} gap={gr:.4f} meanMI={row['mean_MI']:.4f}")

                except Exception as exc:
                    rows.append({
                        "n": args.n,
                        "topology": topology,
                        "seed": seed,
                        "eta_zz": eta_zz,
                        "hx": hx,
                        "hz_random": hz_random,
                        "candidate_L": "failed",
                        "family": "failed",
                        "fit_ok": 0,
                        "error": str(exc),
                    })
                    print("  FAILED:", exc)

    df = pd.DataFrame(rows)
    df.to_csv(out / "paper14_v3_sff_features.csv", index=False)

    valid = df[(df["fit_ok"] == 1) & (df["family"] == "pos")].copy() if "fit_ok" in df else pd.DataFrame()
    high = valid[valid["r2"] >= 0.9].copy() if len(valid) else pd.DataFrame()

    if len(valid):
        summary_by_candidate = (
            valid.groupby(["candidate_L"])
            .agg(
                n_rows=("r2", "count"),
                mean_r2=("r2", "mean"),
                median_r2=("r2", "median"),
                mean_beta=("beta", "mean"),
                std_beta=("beta", "std"),
                mean_C=("C_modular", "mean"),
                std_C=("C_modular", "std"),
                mean_gap=("K_gap_ratio", "mean"),
                max_gap=("K_gap_ratio", "max"),
                mean_t_dip=("t_dip", "mean"),
                mean_t_ramp=("t_ramp", "mean"),
                mean_slope_ramp=("slope_ramp", "mean"),
                mean_delta3=("delta3_proxy", "mean"),
            )
            .reset_index()
        )
    else:
        summary_by_candidate = pd.DataFrame()

    summary_by_candidate.to_csv(out / "paper14_v3_summary_by_candidate.csv", index=False)

    summary = {
        "experiment": "BuP Paper 14 SFF features v3",
        "n": args.n,
        "hilbert_dim": 2 ** args.n,
        "subsystem_size": args.subsystem_size,
        "d_A": 2 ** args.subsystem_size,
        "topologies": args.topologies,
        "seeds": args.seeds,
        "eta_zz_list": args.eta_zz_list,
        "hx_list": args.hx_list,
        "hz_random_list": args.hz_random_list,
        "n_rows": int(len(df)),
        "n_valid_pos": int(len(valid)),
        "n_high_quality_pos": int(len(high)),
        "runtime_seconds": float(time.time() - t0),
    }

    if len(high):
        best = high.sort_values("r2", ascending=False).iloc[0].to_dict()
        safe = {}
        for k, v in best.items():
            if isinstance(v, (int, np.integer)):
                safe[k] = int(v)
            elif isinstance(v, (float, np.floating)):
                safe[k] = float(v) if np.isfinite(v) else None
            else:
                safe[k] = str(v)
        summary["best_high_quality"] = safe

    with open(out / "paper14_v3_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    if not args.no_plots:
        plot_beta_features(df, out)

    print("\n=== Summary ===")
    print(f"rows                 : {len(df)}")
    print(f"valid positive fits  : {len(valid)}")
    print(f"high quality R2>=0.9 : {len(high)}")
    if len(summary_by_candidate):
        print("\nBy candidate:")
        print(summary_by_candidate.to_string(index=False))

    if "best_high_quality" in summary:
        b = summary["best_high_quality"]
        print("\nBest high-quality row:")
        for key in [
            "topology", "seed", "eta_zz", "hx", "hz_random", "candidate_L",
            "r2", "beta", "C_modular", "K_gap_ratio", "t_dip", "t_ramp",
            "slope_ramp", "delta3_proxy",
        ]:
            print(f"  {key}: {b.get(key)}")

    print("\nFiles written:")
    for p in sorted(out.iterdir()):
        if p.is_file():
            print("  -", p)

    print("\nDone.")


if __name__ == "__main__":
    main()

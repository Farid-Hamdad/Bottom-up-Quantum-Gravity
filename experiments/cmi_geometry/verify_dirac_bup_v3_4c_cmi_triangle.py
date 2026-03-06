#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================================================
BUp Geometry Test v3.4c — CMI (triplets) vs MI-triangle + clean nulls
========================================================================
Target:
  CMI(i:j | k) = S(ik) + S(jk) - S(k) - S(ijk)

Leak-aware design:
  - Split PAIRS into TRAIN/TEST, evaluate only on TEST TRIPLETS whose base
    pair (i,j) is NOT in TRAIN PAIRS.

Score modes (--score-mode):
  - sum          : I(i:k) + I(j:k)
  - triangle     : I(i:k) + I(j:k) - I(i:j)
  - triangle_inv : I(i:j) - I(i:k) - I(j:k)   (= -triangle)

Null models (--null-mode):
  - shuffle_score (recommended): permute scores among selected samples (pos+neg)
  - shuffle_cmi              : permute CMI among selected samples (pos+neg)
  - random_k                 : replace k by random k' != i,j (previous v3.4/v3.4b)

Metrics per seed:
  - rho_obs, auc_obs
  - z_rho: (rho_obs - mean(null_rho))/std(null_rho)
  - p_rho: P(null_rho >= rho_obs) empirical (1-sided)
  - z_auc, p_auc similarly

Outputs:
  - results_v3_4c/batch_v3_4c_cmi_triangle_summary.json
  - plots in results_v3_4c/ if --plot

Example:
  python verify_dirac_bup_v3_4c_cmi_triangle.py \
    --n 16 --layers 15 --gates 15 \
    --lams 0.0,0.2,0.4,0.5,0.7 \
    --density-list 0.333 \
    --qpos 0.20 --qneg 0.20 \
    --pairs-sample 4000 --triplets-sample 4000 \
    --train-frac 0.30 \
    --n-seeds 30 --seed0 42 --null 200 \
    --score-mode triangle_inv \
    --null-mode shuffle_score \
    --plot --save
"""

import os
import json
import math
import argparse
from dataclasses import dataclass, asdict
from typing import Tuple, List, Dict

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -------------------------
# Utils
# -------------------------

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def parse_float_list(s: str) -> List[float]:
    if s is None or s.strip() == "":
        return []
    return [float(x.strip()) for x in s.split(",") if x.strip() != ""]

def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if np.sum(m) < 5:
        return float("nan")
    x = x[m]; y = y[m]
    rx = x.argsort().argsort().astype(float)
    ry = y.argsort().argsort().astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    den = float(np.sqrt(np.sum(rx*rx) * np.sum(ry*ry)))
    if den <= 0:
        return float("nan")
    return float(np.sum(rx*ry) / den)

def auc_from_scores_posneg(scores_pos: np.ndarray, scores_neg: np.ndarray) -> float:
    """Probability(score_pos > score_neg) with tie=0.5."""
    sp = np.asarray(scores_pos, float)
    sn = np.asarray(scores_neg, float)
    sp = sp[np.isfinite(sp)]
    sn = sn[np.isfinite(sn)]
    if sp.size < 5 or sn.size < 5:
        return float("nan")
    scores = np.concatenate([sp, sn])
    labels = np.concatenate([np.ones_like(sp), np.zeros_like(sn)])
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(order.size, dtype=float)

    vals = scores[order]
    i = 0
    while i < len(vals):
        j = i + 1
        while j < len(vals) and vals[j] == vals[i]:
            j += 1
        if j - i > 1:
            r = ranks[order[i:j]].mean()
            ranks[order[i:j]] = r
        i = j

    r_pos = ranks[labels == 1]
    n_pos = r_pos.size
    n_neg = sn.size
    auc = (r_pos.sum() - n_pos*(n_pos-1)/2.0) / (n_pos*n_neg)
    return float(auc)

def zscore(x: float, null_values: np.ndarray) -> float:
    v = np.asarray(null_values, float)
    v = v[np.isfinite(v)]
    if not np.isfinite(x) or v.size < 10:
        return float("nan")
    mu = float(np.mean(v))
    sd = float(np.std(v))
    if sd <= 1e-12:
        return float("nan")
    return (x - mu) / sd

def pvalue_one_sided_ge(x: float, null_values: np.ndarray) -> float:
    """Empirical p = P(null >= x), with +1 smoothing."""
    v = np.asarray(null_values, float)
    v = v[np.isfinite(v)]
    if not np.isfinite(x) or v.size < 10:
        return float("nan")
    ge = int(np.sum(v >= x))
    return float((ge + 1) / (v.size + 1))


# -------------------------
# Quantum state generation
# -------------------------

def haar_unitary(dim: int, rng: np.random.Generator) -> np.ndarray:
    z = (rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))) / math.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diag(r)
    ph = d / np.abs(d)
    return q * ph

def apply_two_qubit_gate(psi: np.ndarray, n: int, q1: int, q2: int, U: np.ndarray) -> np.ndarray:
    if q1 == q2:
        return psi
    if q1 > q2:
        q1, q2 = q2, q1
        SW = np.array([[1,0,0,0],
                       [0,0,1,0],
                       [0,1,0,0],
                       [0,0,0,1]], dtype=complex)
        U = SW @ U @ SW

    psi_t = psi.reshape([2]*n)
    axes = list(range(n))
    axes.remove(q1); axes.remove(q2)
    perm = [q1, q2] + axes
    inv = np.argsort(perm)

    x = np.transpose(psi_t, perm).reshape(4, -1)
    x = (U @ x)
    psi_out = np.transpose(x.reshape([2,2] + [2]*(n-2)), inv).reshape(-1)
    return psi_out

def generate_state(n: int, layers: int, gates_per_layer: int, lam: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dim = 2**n
    psi = np.ones(dim, dtype=complex) / math.sqrt(dim)

    neigh = [(i, (i+1) % n) for i in range(n)]
    for _ in range(layers):
        for _g in range(gates_per_layer):
            if rng.random() < (1.0 - lam):
                q1, q2 = neigh[int(rng.integers(0, len(neigh)))]
            else:
                q1, q2 = rng.choice(n, size=2, replace=False).tolist()
            U = haar_unitary(4, rng)
            psi = apply_two_qubit_gate(psi, n, q1, q2, U)

    return psi / np.linalg.norm(psi)


# -------------------------
# Entropy + MI + CMI
# -------------------------

def reduced_rho(psi: np.ndarray, n: int, subset: List[int]) -> np.ndarray:
    subset = sorted(subset)
    k = len(subset)
    if k == 0:
        return np.array([[1.0]], dtype=complex)

    psi_t = psi.reshape([2]*n)
    axes = subset + [i for i in range(n) if i not in subset]
    psi_perm = np.transpose(psi_t, axes).reshape(2**k, 2**(n-k))
    rho = psi_perm @ psi_perm.conj().T
    return (rho + rho.conj().T) * 0.5

def von_neumann_entropy_bits(rho: np.ndarray) -> float:
    if rho.shape == (1, 1):
        return 0.0
    w = np.linalg.eigvalsh(rho)
    w = np.real(w)
    w = np.clip(w, 0.0, 1.0)
    w = w[w > 1e-15]
    if w.size == 0:
        return 0.0
    return float(-np.sum(w * (np.log(w) / np.log(2.0))))

def entropy_subset(psi: np.ndarray, n: int, subset: List[int]) -> float:
    return von_neumann_entropy_bits(reduced_rho(psi, n, subset))

def mutual_information_ij(
    psi: np.ndarray, n: int, i: int, j: int, S1: np.ndarray,
    S2_cache: Dict[Tuple[int, int], float],
) -> float:
    if i == j:
        return 0.0
    a, b = (i, j) if i < j else (j, i)
    key = (a, b)
    if key not in S2_cache:
        S2_cache[key] = entropy_subset(psi, n, [a, b])
    return float(S1[a] + S1[b] - S2_cache[key])

def cmi_i_j_given_k(
    psi: np.ndarray, n: int, i: int, j: int, k: int, S1: np.ndarray,
    S2_cache: Dict[Tuple[int, int], float],
    S3_cache: Dict[Tuple[int, int, int], float],
) -> float:
    if len({i, j, k}) < 3:
        return float("nan")
    a, b = (i, k) if i < k else (k, i)
    c, d = (j, k) if j < k else (k, j)

    if (a, b) not in S2_cache:
        S2_cache[(a, b)] = entropy_subset(psi, n, [a, b])
    if (c, d) not in S2_cache:
        S2_cache[(c, d)] = entropy_subset(psi, n, [c, d])

    tri = tuple(sorted([i, j, k]))
    if tri not in S3_cache:
        S3_cache[tri] = entropy_subset(psi, n, list(tri))

    return float(S2_cache[(a, b)] + S2_cache[(c, d)] - S1[k] - S3_cache[tri])


# -------------------------
# Core evaluation
# -------------------------

@dataclass
class SeedMetrics:
    seed: int
    lam: float
    density: float
    score_mode: str
    null_mode: str
    n_eval: int
    auc: float
    rho: float
    z_rho: float
    p_rho: float
    z_auc: float
    p_auc: float
    null_ok: int

def score_triplet(
    psi: np.ndarray, n: int, i: int, j: int, k: int,
    S1: np.ndarray, S2_cache: Dict[Tuple[int, int], float],
    score_mode: str,
) -> float:
    Iik = mutual_information_ij(psi, n, i, k, S1, S2_cache)
    Ijk = mutual_information_ij(psi, n, j, k, S1, S2_cache)
    if score_mode == "sum":
        return float(Iik + Ijk)

    Iij = mutual_information_ij(psi, n, i, j, S1, S2_cache)

    if score_mode == "triangle":
        return float(Iik + Ijk - Iij)
    if score_mode == "triangle_inv":
        return float(Iij - Iik - Ijk)
    raise ValueError(f"Unknown score_mode: {score_mode}")

def random_k_not_ij(rng: np.random.Generator, n: int, i: int, j: int) -> int:
    while True:
        k = int(rng.integers(0, n))
        if k != i and k != j:
            return k

def evaluate_seed_lam_density(
    n: int, layers: int, gates: int, lam: float, density: float,
    qpos: float, qneg: float,
    pairs_sample: int, triplets_sample: int, train_frac: float,
    n_null: int, seed: int,
    score_mode: str, null_mode: str
) -> SeedMetrics:
    rng = np.random.default_rng(seed)
    psi = generate_state(n=n, layers=layers, gates_per_layer=gates, lam=lam, seed=seed)

    # S1
    S1 = np.zeros(n, float)
    for i in range(n):
        S1[i] = entropy_subset(psi, n, [i])

    S2_cache: Dict[Tuple[int, int], float] = {}
    S3_cache: Dict[Tuple[int, int, int], float] = {}

    # leak-aware pair split
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    rng.shuffle(pairs)
    if pairs_sample and pairs_sample < len(pairs):
        pairs = pairs[:pairs_sample]

    n_train = int(max(1, round(train_frac * len(pairs))))
    train_set = set(pairs[:n_train])

    # sample triplets
    triplets = []
    for _ in range(triplets_sample):
        i, j, k = rng.choice(n, size=3, replace=False).tolist()
        a, b = (i, j) if i < j else (j, i)
        triplets.append((a, b, k))

    test_triplets = [(a, b, k) for (a, b, k) in triplets if (a, b) not in train_set]
    if len(test_triplets) < 40:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           score_mode=score_mode, null_mode=null_mode,
                           n_eval=len(test_triplets),
                           auc=float("nan"), rho=float("nan"),
                           z_rho=float("nan"), p_rho=float("nan"),
                           z_auc=float("nan"), p_auc=float("nan"),
                           null_ok=0)

    # compute cmi + score
    cmi = np.zeros(len(test_triplets), float)
    score = np.zeros(len(test_triplets), float)
    for idx, (a, b, k) in enumerate(test_triplets):
        cmi[idx] = cmi_i_j_given_k(psi, n, a, b, k, S1, S2_cache, S3_cache)
        score[idx] = score_triplet(psi, n, a, b, k, S1, S2_cache, score_mode)

    m = np.isfinite(cmi) & np.isfinite(score)
    cmi = cmi[m]; score = score[m]
    trip_f = [t for ok, t in zip(m.tolist(), test_triplets) if ok]

    if cmi.size < 80:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           score_mode=score_mode, null_mode=null_mode,
                           n_eval=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho=float("nan"), p_rho=float("nan"),
                           z_auc=float("nan"), p_auc=float("nan"),
                           null_ok=0)

    # select extremes for evaluation (as before)
    hi = float(np.quantile(cmi, 1.0 - qpos))
    lo = float(np.quantile(cmi, qneg))
    pos = np.where(cmi >= hi)[0]
    neg = np.where(cmi <= lo)[0]
    if pos.size < 10 or neg.size < 10:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           score_mode=score_mode, null_mode=null_mode,
                           n_eval=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho=float("nan"), p_rho=float("nan"),
                           z_auc=float("nan"), p_auc=float("nan"),
                           null_ok=0)

    sel = np.concatenate([pos, neg])
    cmi_sel = cmi[sel]
    score_sel = score[sel]

    auc_obs = auc_from_scores_posneg(score[pos], score[neg])
    rho_obs = safe_spearman(score_sel, cmi_sel)

    # null distributions
    null_rhos = []
    null_aucs = []
    ok = 0

    for kidx in range(n_null):
        rng0 = np.random.default_rng(seed * 100000 + kidx)

        if null_mode == "shuffle_score":
            s0 = score_sel.copy()
            rng0.shuffle(s0)
            rho0 = safe_spearman(s0, cmi_sel)
            # AUC under permuted scores: recompute with same pos/neg indices
            auc0 = auc_from_scores_posneg(s0[:pos.size], s0[pos.size:])  # WRONG ordering
            # Fix: build s0_pos/s0_neg based on original pos/neg membership in sel
            s0_pos = s0[:pos.size]  # placeholder (will be overwritten below)

        elif null_mode == "shuffle_cmi":
            y0 = cmi_sel.copy()
            rng0.shuffle(y0)
            rho0 = safe_spearman(score_sel, y0)
            auc0 = auc_from_scores_posneg(score_sel[:pos.size], score_sel[pos.size:])  # placeholder

        elif null_mode == "random_k":
            s0 = np.zeros_like(score_sel)
            for ii, idx_sel in enumerate(sel):
                a, b, _k = trip_f[int(idx_sel)]
                k2 = random_k_not_ij(rng0, n, a, b)
                s0[ii] = score_triplet(psi, n, a, b, k2, S1, S2_cache, score_mode)
            rho0 = safe_spearman(s0, cmi_sel)
            # AUC with same split (pos vs neg) on sel
            s0_pos = s0[:pos.size]  # placeholder
            auc0 = auc_from_scores_posneg(s0_pos, s0[pos.size:])  # placeholder

        else:
            raise ValueError(f"Unknown null_mode: {null_mode}")

        # Correct AUC computation for all null modes:
        # sel = [pos, neg] concatenation in that order
        # so first len(pos) entries correspond to pos, remaining to neg
        if null_mode == "shuffle_score":
            s0_pos = s0[:pos.size]
            s0_neg = s0[pos.size:]
            auc0 = auc_from_scores_posneg(s0_pos, s0_neg)
        elif null_mode == "shuffle_cmi":
            # labels are fixed; score fixed; AUC doesn't change under shuffling cmi.
            # Instead, define AUC-null for shuffle_cmi by recomputing labels from shuffled cmi_sel:
            y0 = y0  # already shuffled
            hi0 = float(np.quantile(y0, 1.0 - qpos))
            lo0 = float(np.quantile(y0, qneg))
            p0 = np.where(y0 >= hi0)[0]
            n0 = np.where(y0 <= lo0)[0]
            if p0.size >= 10 and n0.size >= 10:
                auc0 = auc_from_scores_posneg(score_sel[p0], score_sel[n0])
            else:
                auc0 = float("nan")
        elif null_mode == "random_k":
            s0_pos = s0[:pos.size]
            s0_neg = s0[pos.size:]
            auc0 = auc_from_scores_posneg(s0_pos, s0_neg)

        if np.isfinite(rho0) and np.isfinite(auc0):
            ok += 1
        null_rhos.append(rho0)
        null_aucs.append(auc0)

    null_rhos = np.array(null_rhos, float)
    null_aucs = np.array(null_aucs, float)

    z_rho = zscore(rho_obs, null_rhos)
    p_rho = pvalue_one_sided_ge(rho_obs, null_rhos)

    z_auc = zscore(auc_obs, null_aucs)
    p_auc = pvalue_one_sided_ge(auc_obs, null_aucs)

    return SeedMetrics(seed=seed, lam=lam, density=density,
                       score_mode=score_mode, null_mode=null_mode,
                       n_eval=int(cmi.size),
                       auc=float(auc_obs), rho=float(rho_obs),
                       z_rho=float(z_rho), p_rho=float(p_rho),
                       z_auc=float(z_auc), p_auc=float(p_auc),
                       null_ok=int(ok))


def aggregate(metrics: List[SeedMetrics]) -> Dict[str, Dict[str, float]]:
    def agg(x):
        x = np.asarray(x, float)
        msk = np.isfinite(x)
        if np.sum(msk) == 0:
            return {"n": 0, "mean": float("nan"), "std": float("nan")}
        return {"n": int(np.sum(msk)), "mean": float(np.mean(x[msk])), "std": float(np.std(x[msk]))}

    return {
        "auc": agg([m.auc for m in metrics]),
        "rho": agg([m.rho for m in metrics]),
        "z_rho": agg([m.z_rho for m in metrics]),
        "p_rho": agg([m.p_rho for m in metrics]),
        "z_auc": agg([m.z_auc for m in metrics]),
        "p_auc": agg([m.p_auc for m in metrics]),
        "n_eval": agg([m.n_eval for m in metrics]),
        "null_ok": agg([m.null_ok for m in metrics]),
    }


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--layers", type=int, default=15)
    ap.add_argument("--gates", type=int, default=15)

    ap.add_argument("--lams", type=str, default="0.0,0.2,0.4,0.5,0.7")
    ap.add_argument("--density-list", type=str, default="0.333")

    ap.add_argument("--qpos", type=float, default=0.20)
    ap.add_argument("--qneg", type=float, default=0.20)

    ap.add_argument("--pairs-sample", type=int, default=4000)
    ap.add_argument("--triplets-sample", type=int, default=4000)
    ap.add_argument("--train-frac", type=float, default=0.30)

    ap.add_argument("--n-seeds", type=int, default=30)
    ap.add_argument("--seed0", type=int, default=42)
    ap.add_argument("--null", type=int, default=200)

    ap.add_argument("--score-mode", type=str, default="triangle_inv",
                    choices=["sum", "triangle", "triangle_inv"])
    ap.add_argument("--null-mode", type=str, default="shuffle_score",
                    choices=["shuffle_score", "shuffle_cmi", "random_k"])

    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--outdir", type=str, default="results_v3_4c")
    args = ap.parse_args()

    lams = parse_float_list(args.lams)
    densities = parse_float_list(args.density_list) or [0.333]

    print("========================================================================")
    print("BUp Geometry Test v3.4c — CMI(i:j|k) vs MI-triangle + clean nulls")
    print("========================================================================")
    print(f"n={args.n} layers={args.layers} gates={args.gates} lams={lams}")
    print(f"densities={densities} qpos={args.qpos} qneg={args.qneg}")
    print(f"pairs_sample={args.pairs_sample} triplets_sample={args.triplets_sample} train_frac={args.train_frac}")
    print(f"n_seeds={args.n_seeds} seed0={args.seed0} null={args.null}")
    print(f"score_mode={args.score_mode} null_mode={args.null_mode} plot={args.plot} save={args.save}")
    print("------------------------------------------------------------------------")

    ensure_dir(args.outdir)

    all_results = {
        "meta": {
            "version": "v3.4c",
            "n": args.n,
            "layers": args.layers,
            "gates": args.gates,
            "lams": lams,
            "densities": densities,
            "qpos": args.qpos,
            "qneg": args.qneg,
            "pairs_sample": args.pairs_sample,
            "triplets_sample": args.triplets_sample,
            "train_frac": args.train_frac,
            "n_seeds": args.n_seeds,
            "seed0": args.seed0,
            "n_null": args.null,
            "score_mode": args.score_mode,
            "null_mode": args.null_mode,
            "target": "CMI(i:j|k)",
            "leak_aware": True
        },
        "by_lam_density": {}
    }

    plot_rows = []

    for den in densities:
        for lam in lams:
            key = f"lam={lam:.6g}|density={den:.6g}|score={args.score_mode}|null={args.null_mode}"
            print(f"\n=== λ = {lam:.3f} | density={den:.4f} | score={args.score_mode} | null={args.null_mode} ===")

            metrics = []
            for sidx in range(args.n_seeds):
                seed = args.seed0 + sidx
                m = evaluate_seed_lam_density(
                    n=args.n, layers=args.layers, gates=args.gates,
                    lam=lam, density=den,
                    qpos=args.qpos, qneg=args.qneg,
                    pairs_sample=args.pairs_sample,
                    triplets_sample=args.triplets_sample,
                    train_frac=args.train_frac,
                    n_null=args.null,
                    seed=seed,
                    score_mode=args.score_mode,
                    null_mode=args.null_mode
                )
                metrics.append(m)
                if (sidx < 3) or (sidx == args.n_seeds - 1):
                    print(f"[{sidx+1:02d}/{args.n_seeds}] seed={seed}  "
                          f"AUC={m.auc:.3f} rho={m.rho:.3f}  "
                          f"p_rho={m.p_rho:.3g} p_auc={m.p_auc:.3g} "
                          f"z_rho={m.z_rho:.2f} z_auc={m.z_auc:.2f} "
                          f"(n_eval={m.n_eval}, null_ok={m.null_ok}/{args.null})")

            agg = aggregate(metrics)
            all_results["by_lam_density"][key] = {
                "lam": lam,
                "density": den,
                "score_mode": args.score_mode,
                "null_mode": args.null_mode,
                "aggregate": agg,
                "per_seed": [asdict(x) for x in metrics],
            }

            print(f"--- SUMMARY λ={lam:.3f} density={den:.4f} ---")
            print(f"AUC:   {agg['auc']['mean']:.3f} ± {agg['auc']['std']:.3f} (n={agg['auc']['n']})")
            print(f"rho:   {agg['rho']['mean']:.3f} ± {agg['rho']['std']:.3f} (n={agg['rho']['n']})")
            print(f"p_rho: {agg['p_rho']['mean']:.3g} ± {agg['p_rho']['std']:.3g} (n={agg['p_rho']['n']})")
            print(f"p_auc: {agg['p_auc']['mean']:.3g} ± {agg['p_auc']['std']:.3g} (n={agg['p_auc']['n']})")

            plot_rows.append((den, lam, agg["auc"]["mean"], agg["rho"]["mean"],
                              agg["p_rho"]["mean"], agg["p_auc"]["mean"]))

    if args.save:
        out_json = os.path.join(args.outdir, "batch_v3_4c_cmi_triangle_summary.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nSaved: {out_json}")

    if args.plot:
        for den in densities:
            xs   = [lam for (d, lam, aucm, rhom, pr, pa) in plot_rows if abs(d - den) < 1e-12]
            aucs = [aucm for (d, lam, aucm, rhom, pr, pa) in plot_rows if abs(d - den) < 1e-12]
            rhos = [rhom for (d, lam, aucm, rhom, pr, pa) in plot_rows if abs(d - den) < 1e-12]
            prs  = [pr   for (d, lam, aucm, rhom, pr, pa) in plot_rows if abs(d - den) < 1e-12]
            pas  = [pa   for (d, lam, aucm, rhom, pr, pa) in plot_rows if abs(d - den) < 1e-12]

            plt.figure()
            plt.plot(xs, aucs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("AUC (CMI high vs low)")
            plt.title(f"v3.4c AUC vs λ (density={den})")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4c_auc_density_{den:.3f}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(xs, rhos, marker="o")
            plt.xlabel("λ")
            plt.ylabel("Spearman rho(score, CMI)")
            plt.title(f"v3.4c rho vs λ (density={den})")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4c_rho_density_{den:.3f}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(xs, prs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("p-value (rho)")
            plt.title(f"v3.4c p_rho vs λ (density={den})")
            plt.yscale("log")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4c_p_rho_density_{den:.3f}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(xs, pas, marker="o")
            plt.xlabel("λ")
            plt.ylabel("p-value (AUC)")
            plt.title(f"v3.4c p_auc vs λ (density={den})")
            plt.yscale("log")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4c_p_auc_density_{den:.3f}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

        print(f"Saved plots to: {args.outdir}/")


if __name__ == "__main__":
    main()

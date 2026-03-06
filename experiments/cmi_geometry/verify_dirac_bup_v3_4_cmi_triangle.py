#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================================================
BUp Geometry Test v3.4b — CMI (triplets) vs MI-triangle (incl. inverted)
========================================================================
Goal:
  Predict multipartite conditional mutual information (CMI) using 2-body MI
  "entanglement-of-entanglement" scores with a leak-aware split.

Target:
  CMI(i:j | k) = S(ik) + S(jk) - S(k) - S(ijk)

Leak-aware design:
  - Split PAIRS into TRAIN/TEST, and evaluate only on TEST TRIPLETS whose
    base pair (i,j) is NOT in TRAIN PAIRS (same as v3.3/v3.4).
  - Score depends on (i,j,k) via MI terms; no diffusion geometry here.

Score modes (--score-mode):
  - sum          : I(i:k) + I(j:k)
  - triangle     : I(i:k) + I(j:k) - I(i:j)
  - triangle_inv : I(i:j) - I(i:k) - I(j:k)   = -triangle

Metrics per (seed, lam, density):
  - rho: Spearman correlation between score and CMI on labeled eval set
  - AUC: classify top-qpos CMI vs bottom-qneg CMI using score
  - z_rho: z-score of rho vs null model (k-permutation)

Null model:
  - For each evaluated (i,j,k), replace k by random k' != i,j (repeated n_null)
    This preserves the (i,j) base pairs and destroys conditional dependence via k.

Outputs:
  - results_v3_4b/batch_v3_4b_cmi_triangle_summary.json
  - plots in results_v3_4b/ if --plot

Dependencies:
  numpy, matplotlib

Example:
  python verify_dirac_bup_v3_4b_cmi_triangle.py \
    --n 16 --layers 15 --gates 15 \
    --lams 0.0,0.2,0.4,0.5,0.7 \
    --density-list 0.333 \
    --qpos 0.20 --qneg 0.20 \
    --pairs-sample 4000 --triplets-sample 4000 \
    --train-frac 0.30 \
    --n-seeds 30 --seed0 42 --null 100 \
    --score-mode triangle_inv \
    --plot --save
"""

import os
import json
import math
import argparse
from dataclasses import dataclass, asdict
from typing import Tuple, List, Dict

import numpy as np

# matplotlib only if --plot
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

    # average ranks for ties (simple)
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


# -------------------------
# Quantum state generation
# -------------------------

def haar_unitary(dim: int, rng: np.random.Generator) -> np.ndarray:
    """Haar random unitary via QR (Mezzadri)."""
    z = (rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))) / math.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diag(r)
    ph = d / np.abs(d)
    return q * ph

def apply_two_qubit_gate(psi: np.ndarray, n: int, q1: int, q2: int, U: np.ndarray) -> np.ndarray:
    """Apply 4x4 unitary U on qubits (q1,q2) of statevector psi."""
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
    """
    Bottom-up circuit:
      - Start |+>^n
      - Each 2-qubit gate: local neighbor w.p. (1-lam), random pair w.p. lam
      - Haar random SU(4)
    """
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
    psi: np.ndarray,
    n: int,
    i: int,
    j: int,
    S1: np.ndarray,
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
    psi: np.ndarray,
    n: int,
    i: int,
    j: int,
    k: int,
    S1: np.ndarray,
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
    n_triplets_test: int
    auc: float
    rho: float
    z_rho_nullk: float
    null_ok: int

def score_triplet(
    psi: np.ndarray,
    n: int,
    i: int,
    j: int,
    k: int,
    S1: np.ndarray,
    S2_cache: Dict[Tuple[int, int], float],
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
    n: int,
    layers: int,
    gates: int,
    lam: float,
    density: float,
    qpos: float,
    qneg: float,
    pairs_sample: int,
    triplets_sample: int,
    train_frac: float,
    n_null: int,
    seed: int,
    score_mode: str,
) -> SeedMetrics:
    rng = np.random.default_rng(seed)
    psi = generate_state(n=n, layers=layers, gates_per_layer=gates, lam=lam, seed=seed)

    # single-qubit entropies
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
    train_pairs = pairs[:n_train]
    train_set = set(train_pairs)

    # sample triplets, keep those whose base pair not in train
    triplets = []
    for _ in range(triplets_sample):
        i, j, k = rng.choice(n, size=3, replace=False).tolist()
        a, b = (i, j) if i < j else (j, i)
        triplets.append((a, b, k))

    test_triplets = [(a, b, k) for (a, b, k) in triplets if (a, b) not in train_set]
    if len(test_triplets) < 40:
        return SeedMetrics(seed=seed, lam=lam, density=density, score_mode=score_mode,
                           n_triplets_test=len(test_triplets),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_nullk=float("nan"), null_ok=0)

    cmi = np.zeros(len(test_triplets), float)
    score = np.zeros(len(test_triplets), float)
    for idx, (a, b, k) in enumerate(test_triplets):
        cmi[idx] = cmi_i_j_given_k(psi, n, a, b, k, S1, S2_cache, S3_cache)
        score[idx] = score_triplet(psi, n, a, b, k, S1, S2_cache, score_mode)

    m = np.isfinite(cmi) & np.isfinite(score)
    cmi = cmi[m]; score = score[m]
    test_triplets_f = [t for ok, t in zip(m.tolist(), test_triplets) if ok]

    if cmi.size < 40:
        return SeedMetrics(seed=seed, lam=lam, density=density, score_mode=score_mode,
                           n_triplets_test=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_nullk=float("nan"), null_ok=0)

    # labels by quantiles
    hi = float(np.quantile(cmi, 1.0 - qpos))
    lo = float(np.quantile(cmi, qneg))
    pos = np.where(cmi >= hi)[0]
    neg = np.where(cmi <= lo)[0]
    if pos.size < 10 or neg.size < 10:
        return SeedMetrics(seed=seed, lam=lam, density=density, score_mode=score_mode,
                           n_triplets_test=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_nullk=float("nan"), null_ok=0)

    auc = auc_from_scores_posneg(score[pos], score[neg])
    sel = np.concatenate([pos, neg])
    rho = safe_spearman(score[sel], cmi[sel])

    # null: randomize k for each selected (i,j)
    null_rhos = []
    ok = 0
    for kidx in range(n_null):
        rng0 = np.random.default_rng(seed * 100000 + kidx)
        score0 = np.zeros(len(sel), float)
        for ii, idx in enumerate(sel):
            a, b, _k = test_triplets_f[int(idx)]
            k2 = random_k_not_ij(rng0, n, a, b)
            score0[ii] = score_triplet(psi, n, a, b, k2, S1, S2_cache, score_mode)
        rho0 = safe_spearman(score0, cmi[sel])
        if np.isfinite(rho0):
            ok += 1
        null_rhos.append(rho0)

    z_rho_nullk = zscore(rho, np.array(null_rhos, float))

    return SeedMetrics(seed=seed, lam=lam, density=density, score_mode=score_mode,
                       n_triplets_test=int(cmi.size),
                       auc=float(auc), rho=float(rho),
                       z_rho_nullk=float(z_rho_nullk), null_ok=int(ok))


def aggregate(metrics: List[SeedMetrics]) -> Dict[str, float]:
    auc = np.array([m.auc for m in metrics], float)
    rho = np.array([m.rho for m in metrics], float)
    z  = np.array([m.z_rho_nullk for m in metrics], float)
    nt = np.array([m.n_triplets_test for m in metrics], float)

    def agg(x):
        msk = np.isfinite(x)
        if np.sum(msk) == 0:
            return {"n": 0, "mean": float("nan"), "std": float("nan")}
        return {"n": int(np.sum(msk)), "mean": float(np.mean(x[msk])), "std": float(np.std(x[msk]))}

    return {
        "auc": agg(auc),
        "rho": agg(rho),
        "z_rho_nullk": agg(z),
        "n_triplets_test": agg(nt),
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
    ap.add_argument("--density-list", type=str, default="0.333")  # compatibility/logging

    ap.add_argument("--qpos", type=float, default=0.20)
    ap.add_argument("--qneg", type=float, default=0.20)

    ap.add_argument("--pairs-sample", type=int, default=4000)
    ap.add_argument("--triplets-sample", type=int, default=4000)
    ap.add_argument("--train-frac", type=float, default=0.30)

    ap.add_argument("--n-seeds", type=int, default=30)
    ap.add_argument("--seed0", type=int, default=42)
    ap.add_argument("--null", type=int, default=100)

    ap.add_argument("--score-mode", type=str, default="triangle_inv",
                    choices=["sum", "triangle", "triangle_inv"])

    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--outdir", type=str, default="results_v3_4b")
    args = ap.parse_args()

    lams = parse_float_list(args.lams)
    densities = parse_float_list(args.density_list) or [0.333]

    print("========================================================================")
    print("BUp Geometry Test v3.4b — CMI(i:j|k) vs MI-triangle (incl. inverted)")
    print("========================================================================")
    print(f"n={args.n} layers={args.layers} gates={args.gates} lams={lams}")
    print(f"densities={densities} qpos={args.qpos} qneg={args.qneg}")
    print(f"pairs_sample={args.pairs_sample} triplets_sample={args.triplets_sample} train_frac={args.train_frac}")
    print(f"n_seeds={args.n_seeds} seed0={args.seed0} null={args.null}")
    print(f"score_mode={args.score_mode} plot={args.plot} save={args.save}")
    print("------------------------------------------------------------------------")

    ensure_dir(args.outdir)

    all_results = {
        "meta": {
            "version": "v3.4b",
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
            "target": "CMI(i:j|k)",
            "score_mode": args.score_mode,
            "score": {
                "sum": "I(i:k)+I(j:k)",
                "triangle": "I(i:k)+I(j:k)-I(i:j)",
                "triangle_inv": "I(i:j)-I(i:k)-I(j:k) = -triangle",
            }[args.score_mode],
            "null": "randomize k (k != i,j)",
            "leak_aware": True
        },
        "by_lam_density": {}
    }

    plot_rows = []

    for den in densities:
        for lam in lams:
            key = f"lam={lam:.6g}|density={den:.6g}|score={args.score_mode}"
            print(f"\n=== λ = {lam:.3f} | density={den:.4f} | score={args.score_mode} ===")

            metrics = []
            for sidx in range(args.n_seeds):
                seed = args.seed0 + sidx
                m = evaluate_seed_lam_density(
                    n=args.n,
                    layers=args.layers,
                    gates=args.gates,
                    lam=lam,
                    density=den,
                    qpos=args.qpos,
                    qneg=args.qneg,
                    pairs_sample=args.pairs_sample,
                    triplets_sample=args.triplets_sample,
                    train_frac=args.train_frac,
                    n_null=args.null,
                    seed=seed,
                    score_mode=args.score_mode
                )
                metrics.append(m)
                if (sidx < 3) or (sidx == args.n_seeds - 1):
                    print(f"[{sidx+1:02d}/{args.n_seeds}] seed={seed}  "
                          f"AUC={m.auc if np.isfinite(m.auc) else float('nan'):.3f} "
                          f"rho={m.rho if np.isfinite(m.rho) else float('nan'):.3f}  "
                          f"z_rho(null-k)={m.z_rho_nullk if np.isfinite(m.z_rho_nullk) else float('nan'):.2f} "
                          f"(n_trip_test={m.n_triplets_test}, null_ok={m.null_ok}/{args.null})")

            agg = aggregate(metrics)
            all_results["by_lam_density"][key] = {
                "lam": lam,
                "density": den,
                "score_mode": args.score_mode,
                "aggregate": agg,
                "per_seed": [asdict(x) for x in metrics],
            }

            print(f"--- SUMMARY λ={lam:.3f} density={den:.4f} score={args.score_mode} ---")
            print(f"AUC:   {agg['auc']['mean']:.3f} ± {agg['auc']['std']:.3f} (n={agg['auc']['n']})")
            print(f"rho:   {agg['rho']['mean']:.3f} ± {agg['rho']['std']:.3f} (n={agg['rho']['n']})")
            print(f"z_rho: {agg['z_rho_nullk']['mean']:.3f} ± {agg['z_rho_nullk']['std']:.3f} (n={agg['z_rho_nullk']['n']})")

            plot_rows.append((den, lam, agg["auc"]["mean"], agg["rho"]["mean"], agg["z_rho_nullk"]["mean"]))

    if args.save:
        out_json = os.path.join(args.outdir, "batch_v3_4b_cmi_triangle_summary.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nSaved: {out_json}")

    if args.plot:
        for den in densities:
            xs = [lam for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            aucs = [aucm for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            rhos = [rhom for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            zs   = [zm   for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]            
            

            plt.figure()
            plt.plot(xs, aucs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("AUC (CMI high vs low)")
            plt.title(f"v3.4b AUC (density={den}, score={args.score_mode})")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4b_auc_density_{den:.3f}_{args.score_mode}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(xs, rhos, marker="o")
            plt.xlabel("λ")
            plt.ylabel("Spearman rho(score, CMI)")
            plt.title(f"v3.4b rho (density={den}, score={args.score_mode})")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4b_rho_density_{den:.3f}_{args.score_mode}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(xs, zs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("z_rho vs null-k")
            plt.title(f"v3.4b z-score vs null-k (density={den}, score={args.score_mode})")
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(args.outdir, f"v3_4b_zrho_density_{den:.3f}_{args.score_mode}.png"),
                        dpi=160, bbox_inches="tight")
            plt.close()

        print(f"Saved plots to: {args.outdir}/")

if __name__ == "__main__":
    main()

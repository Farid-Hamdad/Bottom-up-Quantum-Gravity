#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
========================================================================
BUp Geometry Test v3.3 — CMI (triplets) vs Diffusion Geometry (leak-free)
========================================================================
Goal:
  Test whether geometry built from pairwise MI (TRAIN pairs only)
  predicts multipartite conditional mutual information:
      CMI(i:j | k) = S(ik) + S(jk) - S(k) - S(ijk)

Leak-free design (important):
  - Build feature-graph G from MI on TRAIN PAIRS only (a subset of pairs).
  - Evaluate on TEST TRIPLETS whose base pair (i,j) is NOT in TRAIN PAIRS.
  - Geometry feature: diffusion distance d_diff(i,j) on G.
  - Target: CMI(i:j|k) computed from the quantum state (triplet entropies).
  - Score: -d_diff(i,j) should correlate with CMI(i:j|k).

Metrics per (seed, lam, density):
  - rho: Spearman correlation between score and CMI on labeled eval set
  - AUC: classify top-qpos CMI vs bottom-qneg CMI using score
  - z_rho: z-score of rho vs degree-preserving rewired null graphs

Outputs:
  - results_v3_3/batch_v3_3_cmi_geometry_summary.json
  - plots in results_v3_3/ if --plot

Dependencies:
  numpy, matplotlib, networkx

Example:
  python verify_dirac_bup_v3_3_cmi_geometry.py \
    --n 16 --layers 15 --gates 15 \
    --lams 0.0,0.2,0.4,0.5,0.7 \
    --density-list 0.333 \
    --diff-t 2 --qpos 0.20 --qneg 0.20 \
    --pairs-sample 4000 --triplets-sample 4000 \
    --train-frac 0.30 \
    --n-seeds 30 --seed0 42 --null 100 \
    --plot --save
"""

import os
import json
import math
import argparse
from dataclasses import dataclass, asdict
from typing import Tuple, List, Dict, Optional

import numpy as np
import networkx as nx

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
    # efficient rank-based AUC
    scores = np.concatenate([sp, sn])
    labels = np.concatenate([np.ones_like(sp), np.zeros_like(sn)])
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(order.size, dtype=float)
    # average ranks for ties
    # (small n; simple tie handling)
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
        # swap ordering inside U: implement by conjugation with swap matrix
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
    Simple bottom-up circuit:
      - Start |+>^n (uniform)
      - Each 2-qubit gate: choose local neighbor w.p. (1-lam), random pair w.p. lam
      - Apply Haar random SU(4)
    """
    rng = np.random.default_rng(seed)
    dim = 2**n
    psi = np.ones(dim, dtype=complex) / math.sqrt(dim)

    # 1D ring neighbors (fast + stable)
    neigh = [(i, (i+1) % n) for i in range(n)]

    for _ in range(layers):
        for _g in range(gates_per_layer):
            if rng.random() < (1.0 - lam):
                q1, q2 = neigh[int(rng.integers(0, len(neigh)))]
            else:
                q1, q2 = rng.choice(n, size=2, replace=False).tolist()
            U = haar_unitary(4, rng)
            psi = apply_two_qubit_gate(psi, n, q1, q2, U)

    # normalize
    psi = psi / np.linalg.norm(psi)
    return psi


# -------------------------
# Entropy + MI + CMI
# -------------------------

def reduced_rho(psi: np.ndarray, n: int, subset: List[int]) -> np.ndarray:
    """Reduced density matrix rho_subset from pure state |psi>."""
    subset = sorted(subset)
    k = len(subset)
    if k == 0:
        return np.array([[1.0]], dtype=complex)

    psi_t = psi.reshape([2]*n)
    axes = subset + [i for i in range(n) if i not in subset]
    psi_perm = np.transpose(psi_t, axes).reshape(2**k, 2**(n-k))
    rho = psi_perm @ psi_perm.conj().T
    rho = (rho + rho.conj().T) * 0.5
    return rho

def von_neumann_entropy_bits(rho: np.ndarray) -> float:
    if rho.shape == (1,1):
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

def mutual_information_ij(psi: np.ndarray, n: int, i: int, j: int, S1: np.ndarray, S2_cache: Dict[Tuple[int,int], float]) -> float:
    if i == j:
        return 0.0
    a, b = (i, j) if i < j else (j, i)
    key = (a, b)
    if key not in S2_cache:
        S2_cache[key] = entropy_subset(psi, n, [a, b])
    return float(S1[a] + S1[b] - S2_cache[key])

def cmi_i_j_given_k(psi: np.ndarray, n: int, i: int, j: int, k: int,
                    S1: np.ndarray,
                    S2_cache: Dict[Tuple[int,int], float],
                    S3_cache: Dict[Tuple[int,int,int], float]) -> float:
    # CMI(i:j|k) = S(ik)+S(jk)-S(k)-S(ijk)
    if len({i, j, k}) < 3:
        return float("nan")
    a, b = (i, k) if i < k else (k, i)
    c, d = (j, k) if j < k else (k, j)

    key_ik = (a, b)
    key_jk = (c, d)
    if key_ik not in S2_cache:
        S2_cache[key_ik] = entropy_subset(psi, n, [a, b])
    if key_jk not in S2_cache:
        S2_cache[key_jk] = entropy_subset(psi, n, [c, d])

    tri = tuple(sorted([i, j, k]))
    if tri not in S3_cache:
        S3_cache[tri] = entropy_subset(psi, n, list(tri))

    return float(S2_cache[key_ik] + S2_cache[key_jk] - S1[k] - S3_cache[tri])


# -------------------------
# Graph + diffusion distance
# -------------------------

def build_feature_graph_from_train_pairs(n: int,
                                        train_pairs: List[Tuple[int,int]],
                                        mi_values: np.ndarray,
                                        density: float,
                                        seed: int) -> nx.Graph:
    """
    Build weighted graph from train pair MI values only, keeping top edges by density.
    Also enforce connectivity by adding a maximum spanning tree on available edges.
    """
    rng = np.random.default_rng(seed)
    G = nx.Graph()
    G.add_nodes_from(range(n))

    # how many edges to keep (relative to complete graph)
    m_target = int(max(1, round(density * (n*(n-1)//2))))

    # sort train edges by MI weight desc
    order = np.argsort(-mi_values)
    chosen = []
    for idx in order:
        if len(chosen) >= m_target:
            break
        u, v = train_pairs[int(idx)]
        w = float(mi_values[int(idx)])
        if not np.isfinite(w):
            continue
        chosen.append((u, v, w))

    # if too few (small train set), just take what we have
    if len(chosen) == 0:
        # fallback: connect line
        for i in range(n-1):
            G.add_edge(i, i+1, weight=1.0)
        return G

    G.add_weighted_edges_from(chosen, weight="weight")

    # enforce connectivity using maximum spanning tree over ALL available chosen edges
    if not nx.is_connected(G):
        H = nx.Graph()
        H.add_nodes_from(range(n))
        H.add_weighted_edges_from(chosen, weight="weight")
        try:
            mst = nx.maximum_spanning_tree(H, weight="weight")
            for (u, v, d) in mst.edges(data=True):
                if not G.has_edge(u, v):
                    G.add_edge(u, v, weight=float(d.get("weight", 1.0)))
        except Exception:
            pass

    # final fallback: connect components arbitrarily
    if not nx.is_connected(G):
        comps = [list(c) for c in nx.connected_components(G)]
        comps.sort(key=len, reverse=True)
        for i in range(len(comps)-1):
            a = int(rng.choice(comps[i]))
            b = int(rng.choice(comps[i+1]))
            G.add_edge(a, b, weight=1e-6)

    return G

def diffusion_distance_matrix(G: nx.Graph, t: float) -> np.ndarray:
    """
    Compute full diffusion distance matrix using combinatorial Laplacian eigendecomp:
      d_t^2(i,j) = sum_k exp(-2 t lambda_k) (v_k(i)-v_k(j))^2
    """
    n = G.number_of_nodes()
    A = nx.to_numpy_array(G, nodelist=range(n), weight="weight", dtype=float)
    A = (A + A.T) * 0.5
    deg = A.sum(axis=1)
    L = np.diag(deg) - A
    L = (L + L.T) * 0.5
    w, V = np.linalg.eigh(L)
    w = np.clip(np.real(w), 0.0, None)
    ew = np.exp(-t * w)
    ew2 = ew * ew

    # compute pairwise distances via matrix trick:
    # d^2(i,j) = sum_k ew2_k (v_i - v_j)^2
    # = sum_k ew2_k v_i^2 + sum_k ew2_k v_j^2 - 2 sum_k ew2_k v_i v_j
    W = V * ew2.reshape(1, -1)  # scale columns by ew2
    s = np.sum(W * V, axis=1)   # sum_k ew2_k v_i^2
    M = V @ (V.T * ew2.reshape(-1, 1))  # V diag(ew2) V^T
    d2 = (s.reshape(-1, 1) + s.reshape(1, -1) - 2.0 * M)
    d2 = np.clip(d2, 0.0, None)
    D = np.sqrt(d2)
    return D

def null_rewire_degree_preserving(G: nx.Graph, seed: int, nswap_mult: float = 10.0) -> nx.Graph:
    """Degree-preserving swaps, then shuffle weights across edges."""
    rng = np.random.default_rng(seed)
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    H.add_edges_from(G.edges())
    m = H.number_of_edges()
    if m <= 1:
        return G.copy()

    nswap = int(max(1, round(nswap_mult * m)))
    try:
        nx.double_edge_swap(H, nswap=nswap, max_tries=nswap * 20,
                            seed=int(rng.integers(0, 2**31 - 1)))
    except Exception:
        pass

    w = [float(d.get("weight", 1.0)) for _, _, d in G.edges(data=True)]
    rng.shuffle(w)
    for (u, v), ww in zip(H.edges(), w):
        H[u][v]["weight"] = float(ww)
    return H


# -------------------------
# Core evaluation
# -------------------------

@dataclass
class SeedMetrics:
    seed: int
    lam: float
    density: float
    n_triplets_test: int
    auc: float
    rho: float
    z_rho_ws: float
    null_ok: int

def evaluate_seed_lam_density(n: int,
                              layers: int,
                              gates: int,
                              lam: float,
                              density: float,
                              diff_t: float,
                              qpos: float,
                              qneg: float,
                              pairs_sample: int,
                              triplets_sample: int,
                              train_frac: float,
                              n_null: int,
                              seed: int) -> SeedMetrics:
    rng = np.random.default_rng(seed)

    psi = generate_state(n=n, layers=layers, gates_per_layer=gates, lam=lam, seed=seed)

    # cache entropies
    S1 = np.zeros(n, float)
    for i in range(n):
        S1[i] = entropy_subset(psi, n, [i])

    S2_cache: Dict[Tuple[int,int], float] = {}
    S3_cache: Dict[Tuple[int,int,int], float] = {}

    # build pair list
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    rng.shuffle(pairs)
    if pairs_sample and pairs_sample < len(pairs):
        pairs = pairs[:pairs_sample]

    # split train pairs
    n_train = int(max(1, round(train_frac * len(pairs))))
    train_pairs = pairs[:n_train]
    train_set = set(train_pairs)

    # compute MI for train pairs only (for graph weights)
    mi_vals = np.zeros(len(train_pairs), float)
    for idx, (i, j) in enumerate(train_pairs):
        mi_vals[idx] = mutual_information_ij(psi, n, i, j, S1, S2_cache)

    G = build_feature_graph_from_train_pairs(n=n, train_pairs=train_pairs, mi_values=mi_vals,
                                             density=density, seed=seed)

    # diffusion distances from G
    D = diffusion_distance_matrix(G, t=float(diff_t))

    # sample triplets (i,j,k) with distinct indices
    triplets = []
    for _ in range(triplets_sample):
        i, j, k = rng.choice(n, size=3, replace=False).tolist()
        # canonical base pair is (min(i,j), max(i,j))
        a, b = (i, j) if i < j else (j, i)
        triplets.append((a, b, k))
    # keep only test triplets whose base pair NOT in train set
    test_triplets = [(a, b, k) for (a, b, k) in triplets if (a, b) not in train_set]

    # compute targets + scores
    if len(test_triplets) < 40:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           n_triplets_test=len(test_triplets),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_ws=float("nan"), null_ok=0)

    cmi = np.zeros(len(test_triplets), float)
    score = np.zeros(len(test_triplets), float)
    for idx, (a, b, k) in enumerate(test_triplets):
        cmi[idx] = cmi_i_j_given_k(psi, n, a, b, k, S1, S2_cache, S3_cache)
        score[idx] = -float(D[a, b])  # smaller distance => higher score

    m = np.isfinite(cmi) & np.isfinite(score)
    cmi = cmi[m]; score = score[m]
    if cmi.size < 40:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           n_triplets_test=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_ws=float("nan"), null_ok=0)

    # labels by quantiles (robust)
    hi = float(np.quantile(cmi, 1.0 - qpos))
    lo = float(np.quantile(cmi, qneg))
    pos = np.where(cmi >= hi)[0]
    neg = np.where(cmi <= lo)[0]
    if pos.size < 10 or neg.size < 10:
        return SeedMetrics(seed=seed, lam=lam, density=density,
                           n_triplets_test=int(cmi.size),
                           auc=float("nan"), rho=float("nan"),
                           z_rho_ws=float("nan"), null_ok=0)

    auc = auc_from_scores_posneg(score[pos], score[neg])
    sel = np.concatenate([pos, neg])
    rho = safe_spearman(score[sel], cmi[sel])

    # nulls (degree-preserving rewire + weight shuffle)
    null_rhos = []
    ok = 0
    for kidx in range(n_null):
        H = null_rewire_degree_preserving(G, seed=seed * 100000 + kidx)
        D0 = diffusion_distance_matrix(H, t=float(diff_t))
        score0 = np.zeros(len(sel), float)
        for ii, idx in enumerate(sel):
            a, b, _kk = test_triplets[int(np.where(m)[0][idx])]  # map back to triplet indices
            score0[ii] = -float(D0[a, b])
        rho0 = safe_spearman(score0, cmi[sel])
        if np.isfinite(rho0):
            ok += 1
        null_rhos.append(rho0)

    z_rho_ws = zscore(rho, np.array(null_rhos, float))

    return SeedMetrics(seed=seed, lam=lam, density=density,
                       n_triplets_test=int(cmi.size),
                       auc=float(auc), rho=float(rho),
                       z_rho_ws=float(z_rho_ws), null_ok=int(ok))


def aggregate(metrics: List[SeedMetrics]) -> Dict[str, float]:
    auc = np.array([m.auc for m in metrics], float)
    rho = np.array([m.rho for m in metrics], float)
    z  = np.array([m.z_rho_ws for m in metrics], float)
    nt = np.array([m.n_triplets_test for m in metrics], float)

    def agg(x):
        msk = np.isfinite(x)
        if np.sum(msk) == 0:
            return {"n": 0, "mean": float("nan"), "std": float("nan")}
        return {"n": int(np.sum(msk)), "mean": float(np.mean(x[msk])), "std": float(np.std(x[msk]))}

    return {
        "auc": agg(auc),
        "rho": agg(rho),
        "z_rho_ws": agg(z),
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
    ap.add_argument("--density-list", type=str, default="0.333")

    ap.add_argument("--diff-t", type=float, default=2.0)
    ap.add_argument("--qpos", type=float, default=0.20)
    ap.add_argument("--qneg", type=float, default=0.20)

    ap.add_argument("--pairs-sample", type=int, default=4000)
    ap.add_argument("--triplets-sample", type=int, default=4000)
    ap.add_argument("--train-frac", type=float, default=0.30)

    ap.add_argument("--n-seeds", type=int, default=30)
    ap.add_argument("--seed0", type=int, default=42)
    ap.add_argument("--null", type=int, default=100)

    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--outdir", type=str, default="results_v3_3")
    args = ap.parse_args()

    lams = parse_float_list(args.lams)
    densities = parse_float_list(args.density_list)

    print("========================================================================")
    print("BUp Geometry Test v3.3 — CMI(i:j|k) vs Diffusion Geometry (leak-free)")
    print("========================================================================")
    print(f"n={args.n} layers={args.layers} gates={args.gates} lams={lams}")
    print(f"densities={densities} diff_t={args.diff_t} qpos={args.qpos} qneg={args.qneg}")
    print(f"pairs_sample={args.pairs_sample} triplets_sample={args.triplets_sample} train_frac={args.train_frac}")
    print(f"n_seeds={args.n_seeds} seed0={args.seed0} null={args.null}")
    print(f"plot={args.plot} save={args.save}")
    print("------------------------------------------------------------------------")

    ensure_dir(args.outdir)

    all_results = {
        "meta": {
            "version": "v3.3",
            "n": args.n,
            "layers": args.layers,
            "gates": args.gates,
            "lams": lams,
            "densities": densities,
            "diff_t": args.diff_t,
            "qpos": args.qpos,
            "qneg": args.qneg,
            "pairs_sample": args.pairs_sample,
            "triplets_sample": args.triplets_sample,
            "train_frac": args.train_frac,
            "n_seeds": args.n_seeds,
            "seed0": args.seed0,
            "n_null": args.null,
            "leak_free": True,
            "target": "CMI(i:j|k)",
            "score": "-diffusion_distance(i,j)",
        },
        "by_lam_density": {}
    }

    # store for plotting
    plot_rows = []

    for den in densities:
        for lam in lams:
            key = f"lam={lam:.6g}|density={den:.6g}"
            print(f"\n=== λ = {lam:.3f} | density={den:.4f} ===")
            metrics = []
            for sidx in range(args.n_seeds):
                seed = args.seed0 + sidx
                m = evaluate_seed_lam_density(
                    n=args.n,
                    layers=args.layers,
                    gates=args.gates,
                    lam=lam,
                    density=den,
                    diff_t=args.diff_t,
                    qpos=args.qpos,
                    qneg=args.qneg,
                    pairs_sample=args.pairs_sample,
                    triplets_sample=args.triplets_sample,
                    train_frac=args.train_frac,
                    n_null=args.null,
                    seed=seed
                )
                metrics.append(m)
                if (sidx < 3) or (sidx == args.n_seeds - 1):
                    print(f"[{sidx+1:02d}/{args.n_seeds}] seed={seed}  "
                          f"AUC={m.auc if np.isfinite(m.auc) else float('nan'):.3f} "
                          f"rho={m.rho if np.isfinite(m.rho) else float('nan'):.3f}  "
                          f"z_rho(ws)={m.z_rho_ws if np.isfinite(m.z_rho_ws) else float('nan'):.2f} "
                          f"(n_trip_test={m.n_triplets_test}, null_ok={m.null_ok}/{args.null})")

            agg = aggregate(metrics)
            all_results["by_lam_density"][key] = {
                "lam": lam,
                "density": den,
                "aggregate": agg,
                "per_seed": [asdict(x) for x in metrics],
            }

            print(f"--- SUMMARY λ={lam:.3f} density={den:.4f} ---")
            print(f"AUC:   {agg['auc']['mean']:.3f} ± {agg['auc']['std']:.3f} (n={agg['auc']['n']})")
            print(f"rho:   {agg['rho']['mean']:.3f} ± {agg['rho']['std']:.3f} (n={agg['rho']['n']})")
            print(f"z_rho: {agg['z_rho_ws']['mean']:.3f} ± {agg['z_rho_ws']['std']:.3f} (n={agg['z_rho_ws']['n']})")

            plot_rows.append((den, lam, agg["auc"]["mean"], agg["rho"]["mean"], agg["z_rho_ws"]["mean"]))

    # save json
    if args.save:
        out_json = os.path.join(args.outdir, "batch_v3_3_cmi_geometry_summary.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(all_results, f, indent=2)
        print(f"\nSaved: {out_json}")

    # plots
    if args.plot:
        # group by density
        for den in densities:
            xs = [lam for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            aucs = [aucm for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            rhos = [rhom for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]
            zs   = [zm   for (d, lam, aucm, rhom, zm) in plot_rows if abs(d - den) < 1e-12]

            # AUC plot
            plt.figure()
            plt.plot(xs, aucs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("AUC (CMI high vs low)")
            plt.title(f"v3.3 CMI-Geometry AUC vs λ (density={den})")
            plt.grid(True, alpha=0.3)
            p = os.path.join(args.outdir, f"v3_3_auc_density_{den:.3f}.png")
            plt.savefig(p, dpi=160, bbox_inches="tight")
            plt.close()

            # rho plot
            plt.figure()
            plt.plot(xs, rhos, marker="o")
            plt.xlabel("λ")
            plt.ylabel("Spearman rho(score, CMI)")
            plt.title(f"v3.3 CMI-Geometry rho vs λ (density={den})")
            plt.grid(True, alpha=0.3)
            p = os.path.join(args.outdir, f"v3_3_rho_density_{den:.3f}.png")
            plt.savefig(p, dpi=160, bbox_inches="tight")
            plt.close()

            # z plot
            plt.figure()
            plt.plot(xs, zs, marker="o")
            plt.xlabel("λ")
            plt.ylabel("z_rho vs rewired null")
            plt.title(f"v3.3 CMI-Geometry z-score vs λ (density={den})")
            plt.grid(True, alpha=0.3)
            p = os.path.join(args.outdir, f"v3_3_zrho_density_{den:.3f}.png")
            plt.savefig(p, dpi=160, bbox_inches="tight")
            plt.close()

        print(f"Saved plots to: {args.outdir}/")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bottom-Up BH Benchmark (N=16)
+ HIE fixed-size enrichi avec seeds par communautés (Louvain) + fallback networkx
Script complet prêt à lancer.

Dépendances optionnelles (recommandées) :
  pip install networkx python-louvain

Si python-louvain n'est pas dispo, on utilise un fallback networkx (greedy_modularity_communities).
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable, Optional, Sequence, Tuple, List, Dict
from numpy.linalg import eigvalsh

# =========================
# Community detection imports (optional)
# =========================
HAS_NX = False
HAS_LOUVAIN = False
try:
    import networkx as nx
    HAS_NX = True
except Exception:
    HAS_NX = False

community_louvain = None
if HAS_NX:
    # python-louvain has different import styles depending on install
    try:
        import community as community_louvain  # type: ignore
        # community might be something else; ensure best_partition exists
        if hasattr(community_louvain, "best_partition"):
            HAS_LOUVAIN = True
        else:
            HAS_LOUVAIN = False
    except Exception:
        HAS_LOUVAIN = False

    if not HAS_LOUVAIN:
        try:
            import community.community_louvain as community_louvain  # type: ignore
            HAS_LOUVAIN = hasattr(community_louvain, "best_partition")
        except Exception:
            HAS_LOUVAIN = False


# =========================
# PARAMÈTRES
# =========================
N = 16
SEED = 20260219
rng = np.random.default_rng(SEED)

GLOBAL_SCRAMBLE_LAYERS = 8
GLOBAL_GATES_PER_LAYER = 10
RHO_TARGET = 1 / 3
EPS = 1e-12
A_SIZE = N // 2  # k=8

POOL_FOR_NORM = 600
N_CANDIDATES = 9000
HILL_STEPS = 6000
TEMPERATURE = 0.10

W_S = 1.0
W_PHI = 1.3
W_INT = 0.8

N_RANDOM_REGIONS = 3000
MAX_DIM_FOR_MODULAR = 256

# HIE classic
HIE_MAX_DEPTH = 3
HIE_MIN_SIZE = 3
HIE_TOP_K = 10

# HIE fixed-size (enrichi)
HIE8_SEEDS = 140          # augmente nettement l'exploration
HIE8_REFINE_STEPS = 12000
HIE8_REFINE_T = 0.10
HIE8_SEED_MODE = "mixed"  # "communities"|"spectral"|"random"|"mixed"
HIE8_COMM_METHOD = "louvain"  # "louvain" ou "greedy"

# NOUVEAU: Paramètres HIE horizon-aware (B)
W_PHI_HORIZON = 3.0       # poids fort sur phi pour privilégier les bottlenecks
W_S_HORIZON = 1.0         # bonus entropie élevée
W_CUT_HORIZON = 1.0       # bonus cut faible (négatif sur z(cut))

OUTDIR = "results_bh_benchmark_N16_horizon_aware"
os.makedirs(OUTDIR, exist_ok=True)


# =========================
# HIE / METRICS
# =========================
@dataclass
class HIECandidate:
    nodes: Tuple[int, ...]
    size: int
    I_int: float
    A_cut: float
    A_per_node: float
    volA: float
    volB: float
    conductance: float
    cohesion: float
    S_A: Optional[float] = None
    area_law_ratio: Optional[float] = None
    score: float = 0.0


def _validate_W(W: np.ndarray) -> np.ndarray:
    W = np.asarray(W, dtype=float)
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError("Matrix must be square (N,N).")
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, 0.0)
    return W


def _connected_components_from_W(W: np.ndarray) -> List[List[int]]:
    W = _validate_W(W)
    n = W.shape[0]
    support = (W > 0)
    visited = np.zeros(n, dtype=bool)
    comps: List[List[int]] = []
    adj = [np.flatnonzero(support[i]).tolist() for i in range(n)]
    for s in range(n):
        if visited[s]:
            continue
        stack = [s]
        visited[s] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
        comps.append(comp)
    return comps


def _greedy_split_by_spectral(W: np.ndarray, nodes: Sequence[int]) -> List[List[int]]:
    idx = np.array(list(nodes), dtype=int)
    if idx.size < 4:
        return [idx.tolist()]
    Wsub = W[np.ix_(idx, idx)]
    deg = Wsub.sum(axis=1)
    if np.all(deg < 1e-12):
        return [idx.tolist()]
    d_inv_sqrt = np.zeros_like(deg)
    nz = deg > 1e-12
    d_inv_sqrt[nz] = 1.0 / np.sqrt(deg[nz])
    Dm = np.diag(d_inv_sqrt)
    L = np.eye(len(idx)) - Dm @ Wsub @ Dm
    vals, vecs = np.linalg.eigh(L)
    if len(vals) < 2:
        return [idx.tolist()]
    fiedler = vecs[:, 1]
    A_mask = fiedler >= 0
    if A_mask.sum() in (0, len(idx)):
        med = np.median(fiedler)
        A_mask = fiedler >= med
    A = idx[A_mask].tolist()
    B = idx[~A_mask].tolist()
    if len(A) == 0 or len(B) == 0:
        return [idx.tolist()]
    return [A, B]


def propose_regions(W: np.ndarray, max_depth: int = 2, min_size: int = 3) -> List[Tuple[int, ...]]:
    W = _validate_W(W)
    regions: List[List[int]] = []
    comps = _connected_components_from_W(W)

    def rec(nodes: List[int], depth: int):
        if len(nodes) < 2 * min_size or depth >= max_depth:
            regions.append(nodes)
            return
        parts = _greedy_split_by_spectral(W, nodes)
        if len(parts) == 1:
            regions.append(nodes)
            return
        for p in parts:
            rec(p, depth + 1)

    for c in comps:
        if len(c) >= min_size:
            rec(c, 0)

    uniq = {tuple(sorted(r)) for r in regions if len(r) >= min_size}
    return sorted(list(uniq), key=lambda x: (len(x), x))


def compute_metrics(W: np.ndarray, region: Sequence[int],
                    entropy_fn: Optional[Callable[[Sequence[int]], float]] = None) -> HIECandidate:
    W = _validate_W(W)
    n = W.shape[0]
    A = np.array(sorted(region), dtype=int)
    if A.size == 0 or A.size >= n:
        raise ValueError("Region must be non-empty and not equal to all nodes.")
    B_mask = np.ones(n, dtype=bool)
    B_mask[A] = False
    B = np.flatnonzero(B_mask)

    W_AA = W[np.ix_(A, A)]
    I_int = float(np.triu(W_AA, 1).sum())  # (i<j) internal weight once
    W_AB = W[np.ix_(A, B)]
    A_cut = float(W_AB.sum())

    deg = W.sum(axis=1)
    volA = float(deg[A].sum())
    volB = float(deg[B].sum())
    denom = min(volA, volB) if min(volA, volB) > 1e-12 else 1e-12
    conductance = float(A_cut / denom)

    size = int(A.size)
    cohesion = float(I_int / max(size * (size - 1) / 2, 1e-12))  # density-like
    A_per_node = float(A_cut / size)

    S_A = None
    area_law_ratio = None
    if entropy_fn is not None:
        S_A = float(entropy_fn(A.tolist()))
        area_law_ratio = float(S_A / (A_cut + 1e-12))

    return HIECandidate(
        nodes=tuple(int(x) for x in A.tolist()),
        size=size,
        I_int=I_int,
        A_cut=A_cut,
        A_per_node=A_per_node,
        volA=volA,
        volB=volB,
        conductance=conductance,
        cohesion=cohesion,
        S_A=S_A,
        area_law_ratio=area_law_ratio,
        score=0.0
    )


def rank_hie_candidates(
    W: np.ndarray,
    regions: List[Tuple[int, ...]],
    entropy_fn: Optional[Callable[[Sequence[int]], float]] = None,
    *,
    mad_floor: float = 1e-6,
    std_floor: float = 1e-6,
    clip_z: float = 20.0,
    w_phi: float = 1.0,
    w_coh: float = 1.0,
    w_apn: float = 0.5,
    w_entropy_bonus: float = 0.75,
    # NOUVEAU: paramètres horizon-aware (C)
    w_S: float = 0.0,
    w_cut: float = 0.0,
    mode: str = "classic",  # "classic" ou "horizon"
) -> List[HIECandidate]:
    W = _validate_W(W)
    if not regions:
        return []

    cand = [compute_metrics(W, r, entropy_fn=entropy_fn) for r in regions]

    cond = np.array([c.conductance for c in cand], dtype=float)
    coh = np.array([c.cohesion for c in cand], dtype=float)
    apn = np.array([c.A_per_node for c in cand], dtype=float)
    
    # NOUVEAU: extraire S_A et A_cut pour mode horizon (C)
    S_vals = np.array([c.S_A if c.S_A is not None else np.nan for c in cand], dtype=float)
    cut_vals = np.array([c.A_cut for c in cand], dtype=float)

    def z_robust(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        med = float(np.median(x))
        mad = float(np.median(np.abs(x - med)))
        if not np.isfinite(mad) or mad < mad_floor:
            mu = float(np.mean(x))
            sig = float(np.std(x))
            sig = sig if (np.isfinite(sig) and sig > std_floor) else std_floor
            z = (x - mu) / sig
        else:
            z = (x - med) / max(mad, mad_floor)
        z = np.nan_to_num(z, nan=0.0, posinf=clip_z, neginf=-clip_z)
        z = np.clip(z, -clip_z, clip_z)
        return z

    z_phi = z_robust(cond)
    z_coh = z_robust(coh)
    z_apn = z_robust(apn)
    
    # NOUVEAU: z-scores pour S et cut (C)
    z_S = z_robust(S_vals)
    z_cut = z_robust(cut_vals)

    if mode == "horizon":
        # Score horizon: pénalise phi fort, récompense S élevé, pénalise cut élevé
        s = (-w_phi * z_phi) + (w_coh * z_coh) + (-w_apn * z_apn) + (w_S * z_S) + (-w_cut * z_cut)
    else:
        # Score classic original
        s = (-w_phi * z_phi) + (w_coh * z_coh) + (-w_apn * z_apn)

    if entropy_fn is not None and mode != "horizon":
        # Bonus area law seulement en mode classic
        ratios = np.array(
            [c.area_law_ratio if (c.area_law_ratio is not None and np.isfinite(c.area_law_ratio)) else np.nan
             for c in cand],
            dtype=float
        )
        ok = np.isfinite(ratios)
        if ok.sum() >= 2:
            target = float(np.nanmedian(ratios))
            dist = np.abs(ratios - target)
            z_dist = z_robust(dist)
            s += (-w_entropy_bonus * z_dist)

    for i, c in enumerate(cand):
        c.score = float(s[i])

    cand.sort(key=lambda x: x.score, reverse=True)
    return cand


def detect_hie(
    W: np.ndarray,
    max_depth: int = 2,
    min_size: int = 3,
    entropy_fn: Optional[Callable[[Sequence[int]], float]] = None,
    top_k: int = 10
) -> Dict[str, object]:
    W = _validate_W(W)
    regions = propose_regions(W, max_depth=max_depth, min_size=min_size)
    ranked = rank_hie_candidates(W, regions, entropy_fn=entropy_fn)
    return {
        "N": W.shape[0],
        "n_regions": len(regions),
        "regions": regions,
        "top": ranked[:top_k],
        "all": ranked
    }


# ============================================================
# FIXED-SIZE HIE WITH LOUVAIN SEEDS
# ============================================================

def _spectral_bipartition(W: np.ndarray) -> Tuple[List[int], List[int]]:
    W = _validate_W(W)
    n = W.shape[0]
    deg = W.sum(axis=1)
    if np.all(deg < 1e-12):
        return list(range(n // 2)), list(range(n // 2, n))

    d_inv_sqrt = np.zeros_like(deg)
    nz = deg > 1e-12
    d_inv_sqrt[nz] = 1.0 / np.sqrt(deg[nz])
    Dm = np.diag(d_inv_sqrt)
    L = np.eye(n) - Dm @ W @ Dm
    vals, vecs = np.linalg.eigh(L)
    fiedler = vecs[:, 1] if n >= 2 else np.ones(n)

    A_mask = fiedler >= 0
    if A_mask.sum() in (0, n):
        med = np.median(fiedler)
        A_mask = fiedler >= med

    A = np.flatnonzero(A_mask).tolist()
    B = np.flatnonzero(~A_mask).tolist()
    if len(A) == 0 or len(B) == 0:
        order = np.argsort(fiedler)
        A = order[:n // 2].tolist()
        B = order[n // 2:].tolist()
    return A, B


def _ensure_size_k_by_greedy(W: np.ndarray, A: Sequence[int], k: int) -> List[int]:
    W = _validate_W(W)
    n = W.shape[0]
    Aset = set(int(x) for x in A)
    if k <= 0 or k >= n:
        raise ValueError("k must be in [1, N-1].")

    deg = W.sum(axis=1)

    def phi_of(S: set) -> float:
        S_list = sorted(S)
        T_list = [i for i in range(n) if i not in S]
        cut = float(np.sum(W[np.ix_(S_list, T_list)])) if S_list and T_list else 0.0
        volA = float(deg[S_list].sum()) if S_list else 0.0
        volB = float(deg[T_list].sum()) if T_list else 0.0
        return float(cut / max(min(volA, volB), 1e-12))

    # Too large -> remove nodes
    while len(Aset) > k:
        current_phi = phi_of(Aset)
        best_node = None
        best_phi = current_phi
        for u in list(Aset):
            S2 = set(Aset)
            S2.remove(u)
            if len(S2) == 0:
                continue
            p2 = phi_of(S2)
            if p2 < best_phi:
                best_phi = p2
                best_node = u
        if best_node is None:
            best_node = min(Aset, key=lambda x: deg[x])
        Aset.remove(best_node)

    # Too small -> add nodes
    while len(Aset) < k:
        current_phi = phi_of(Aset)
        best_node = None
        best_phi = current_phi
        outside = [i for i in range(n) if i not in Aset]
        for u in outside:
            S2 = set(Aset)
            S2.add(u)
            if len(S2) == n:
                continue
            p2 = phi_of(S2)
            if p2 < best_phi:
                best_phi = p2
                best_node = u
        if best_node is None:
            best_node = max(outside, key=lambda x: deg[x])
        Aset.add(best_node)

    return sorted(Aset)


def _local_refine_swaps(
    W: np.ndarray,
    A: Sequence[int],
    k: int,
    steps: int = 8000,
    T: float = 0.10,
    w_phi: float = 1.0,
    w_int: float = 0.8,
    seed: Optional[int] = None
) -> List[int]:
    W = _validate_W(W)
    n = W.shape[0]
    rr = np.random.default_rng(seed)

    Aset = set(int(x) for x in A)
    if len(Aset) != k:
        Aset = set(_ensure_size_k_by_greedy(W, sorted(Aset), k))

    deg = W.sum(axis=1)

    def internal(S: List[int]) -> float:
        W_AA = W[np.ix_(S, S)]
        return float(np.triu(W_AA, 1).sum())

    def phi(S: List[int]) -> float:
        Sset = set(S)
        B = [i for i in range(n) if i not in Sset]
        cut = float(np.sum(W[np.ix_(S, B)])) if S and B else 0.0
        volA = float(deg[S].sum()) if S else 0.0
        volB = float(deg[B].sum()) if B else 0.0
        return float(cut / max(min(volA, volB), 1e-12))

    # pool for z-normalization
    pool_phi = []
    pool_int = []
    for _ in range(250):
        S = sorted(rr.choice(n, size=k, replace=False).tolist())
        pool_phi.append(phi(S))
        pool_int.append(internal(S))
    muP, sigP = float(np.mean(pool_phi)), float(np.std(pool_phi) + 1e-12)
    muI, sigI = float(np.mean(pool_int)), float(np.std(pool_int) + 1e-12)

    def z(x, mu, sig):
        return (x - mu) / sig

    def score(S: List[int]) -> float:
        return (-w_phi * z(phi(S), muP, sigP) + w_int * z(internal(S), muI, sigI))

    cur = sorted(Aset)
    cur_sc = score(cur)

    for _ in range(steps):
        out = int(rr.choice(cur))
        outside = [i for i in range(n) if i not in set(cur)]
        inn = int(rr.choice(outside))
        cand = sorted([x for x in cur if x != out] + [inn])

        sc = score(cand)
        if sc > cur_sc:
            cur, cur_sc = cand, sc
        else:
            if rr.random() < np.exp((sc - cur_sc) / max(T, 1e-6)):
                cur, cur_sc = cand, sc

    return cur


def _community_partition(W: np.ndarray, method: str = "louvain", seed: int = 0) -> List[List[int]]:
    if not HAS_NX:
        raise ImportError("networkx is required for community detection.")

    W = _validate_W(W)
    n = W.shape[0]

    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            w = float(W[i, j])
            if w > 0:
                G.add_edge(i, j, weight=w)

    method = method.lower().strip()
    comms: List[List[int]] = []

    if method == "louvain" and HAS_LOUVAIN:
        try:
            part = community_louvain.best_partition(G, weight="weight", random_state=seed)
        except TypeError:
            part = community_louvain.best_partition(G, weight="weight")
        tmp: Dict[int, List[int]] = {}
        for node, cid in part.items():
            tmp.setdefault(int(cid), []).append(int(node))
        comms = [sorted(v) for v in tmp.values()]
        return comms

    # fallback: greedy modularity communities
    from networkx.algorithms.community import greedy_modularity_communities
    sets = greedy_modularity_communities(G, weight="weight")
    comms = [sorted(list(s)) for s in sets]
    return comms


def _fixed_size_from_communities(
    W: np.ndarray,
    k: int,
    comms: List[List[int]],
    *,
    max_unions: int = 250,
    seed: int = 0
) -> List[Tuple[int, ...]]:
    W = _validate_W(W)
    n = W.shape[0]
    rr = np.random.default_rng(seed)

    base_sets = [set(c) for c in comms if len(c) >= 1]
    unions: List[set] = []
    unions.extend(base_sets)

    if len(base_sets) >= 2:
        pairs = []
        for i in range(len(base_sets)):
            for j in range(i + 1, len(base_sets)):
                pairs.append((i, j))
        rr.shuffle(pairs)
        for (i, j) in pairs[:max_unions // 2]:
            unions.append(set(base_sets[i] | base_sets[j]))

    if len(base_sets) >= 3:
        for _ in range(min(max_unions // 2, 200)):
            i, j, l = rr.choice(len(base_sets), size=3, replace=False).tolist()
            unions.append(set(base_sets[i] | base_sets[j] | base_sets[l]))

    cand = []
    for S in unions:
        S_list = sorted(S)
        if len(S_list) == 0:
            continue
        if len(S_list) != k:
            S_list = _ensure_size_k_by_greedy(W, S_list, k)
        cand.append(tuple(sorted(S_list)))

    cand = sorted(list(set(cand)), key=lambda x: (len(x), x))
    return cand


def propose_regions_fixed_size(
    W: np.ndarray,
    k: int,
    *,
    n_seeds: int = 80,
    refine_steps: int = 8000,
    refine_T: float = 0.10,
    seed_mode: str = "mixed",          # "random" | "spectral" | "communities" | "mixed"
    community_method: str = "louvain",
    forced_seeds: Optional[List[Sequence[int]]] = None,
    seed: int = 0
) -> List[Tuple[int, ...]]:
    W = _validate_W(W)
    n = W.shape[0]
    if not (1 <= k <= n - 1):
        raise ValueError("k must be in [1, N-1].")

    rr = np.random.default_rng(seed)
    candidates: List[Tuple[int, ...]] = []

    # --- KEEP FORCED SEEDS (exact)
    forced_exact: List[Tuple[int, ...]] = []
    if forced_seeds:
        for fs in forced_seeds:
            fs = sorted({int(x) for x in fs})
            if len(fs) == 0:
                continue
            if len(fs) != k:
                fs = _ensure_size_k_by_greedy(W, fs, k)
            t = tuple(sorted(fs))
            forced_exact.append(t)
            candidates.append(t)

    seed_mode = seed_mode.lower().strip()

    # spectral seeds
    if seed_mode in ("spectral", "mixed"):
        A0, B0 = _spectral_bipartition(W)
        candidates.append(tuple(_ensure_size_k_by_greedy(W, A0, k)))
        candidates.append(tuple(_ensure_size_k_by_greedy(W, B0, k)))

    # community seeds
    if seed_mode in ("communities", "mixed"):
        if HAS_NX:
            try:
                comms = _community_partition(W, method=community_method, seed=seed + 123)
                comm_cand = _fixed_size_from_communities(W, k, comms, max_unions=250, seed=seed + 456)
                candidates.extend(comm_cand)
            except Exception as e:
                print(f"Warning: community seeds failed ({e}). Continuing without.")
        else:
            print("Warning: networkx not available. Continuing without community seeds.")

    # random seeds
    if seed_mode in ("random", "mixed"):
        for _ in range(max(0, n_seeds)):
            reg = sorted(rr.choice(n, size=k, replace=False).tolist())
            candidates.append(tuple(reg))

    # deduplicate before refinement
    candidates = list(set(tuple(sorted(c)) for c in candidates if len(c) == k))
    rr.shuffle(candidates)

    refined: List[Tuple[int, ...]] = []

    # refine ONLY non-forced seeds (so forced remain exact)
    forced_set = set(forced_exact)
    for i, cand in enumerate(candidates):
        if cand in forced_set:
            refined.append(cand)  # keep exact, no refine
            continue
        ref = _local_refine_swaps(
            W, list(cand), k,
            steps=refine_steps,
            T=refine_T,
            w_phi=1.0,
            w_int=0.8,
            seed=seed + 1000 + i
        )
        refined.append(tuple(sorted(ref)))

    # ensure forced exact always present
    refined.extend(forced_exact)

    uniq = sorted(list(set(refined)), key=lambda x: (len(x), x))
    return uniq


def detect_hie_fixed_size(
    W: np.ndarray,
    k: int,
    entropy_fn: Optional[Callable[[Sequence[int]], float]] = None,
    top_k: int = 10,
    *,
    n_seeds: int = 80,
    refine_steps: int = 8000,
    refine_T: float = 0.10,
    seed_mode: str = "mixed",
    community_method: str = "louvain",
    forced_seeds: Optional[List[Sequence[int]]] = None,
    seed: int = 0,
    # NOUVEAU: paramètres pour mode horizon
    mode: str = "classic",
    w_phi: float = 1.0,
    w_S: float = 0.0,
    w_cut: float = 0.0,
) -> Dict[str, object]:
    W = _validate_W(W)
    regions = propose_regions_fixed_size(
        W, k,
        n_seeds=n_seeds,
        refine_steps=refine_steps,
        refine_T=refine_T,
        seed_mode=seed_mode,
        community_method=community_method,
        forced_seeds=forced_seeds,
        seed=seed
    )
    ranked = rank_hie_candidates(W, regions, entropy_fn=entropy_fn, 
                                  mode=mode, w_phi=w_phi, w_S=w_S, w_cut=w_cut)
    return {
        "N": W.shape[0],
        "fixed_size": k,
        "n_regions": len(regions),
        "regions": regions,
        "top": ranked[:top_k],
        "all": ranked
    }


# =========================
# OUTILS QUANTIQUES
# =========================
def random_two_qubit_unitary(rng_):
    X = (rng_.normal(size=(4, 4)) + 1j * rng_.normal(size=(4, 4)))
    Q, _ = np.linalg.qr(X)
    return Q


def apply_two_qubit_gate(psi, U, q1, q2, n):
    if q1 == q2:
        return psi
    if q1 > q2:
        q1, q2 = q2, q1
    tensor = psi.reshape([2] * n)
    axes = list(range(n))
    axes[0], axes[q1] = axes[q1], axes[0]
    axes[1], axes[q2] = axes[q2], axes[1]
    tensor = np.transpose(tensor, axes).reshape(4, -1)
    tensor = (U @ tensor).reshape([2, 2] + [2] * (n - 2))
    inv = np.argsort(axes)
    return np.transpose(tensor, inv).reshape(-1)


def reduced_density_matrix(psi, keep, n):
    keep = sorted(list(keep))
    rest = [i for i in range(n) if i not in keep]
    psi_t = psi.reshape([2] * n)
    perm = keep + rest
    psi_p = np.transpose(psi_t, perm)
    dk = 2 ** len(keep)
    dr = 2 ** len(rest)
    psi_m = psi_p.reshape(dk, dr)
    rho = psi_m @ np.conjugate(psi_m.T)
    rho = 0.5 * (rho + rho.conjugate().T)
    tr = np.trace(rho)
    if np.abs(tr) > 1e-12:
        rho = rho / tr
    return rho


def von_neumann_entropy(rho, eps=1e-15):
    vals = eigvalsh(rho)
    vals = np.clip(vals, 0.0, 1.0)
    vals = vals[vals > eps]
    return 0.0 if len(vals) == 0 else float(-np.sum(vals * np.log2(vals)))


def r_stat_from_hermitian(H, eps=1e-12):
    e = np.sort(eigvalsh(H))
    s = np.diff(e)
    if len(s) < 3:
        return np.nan
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:] + eps)
    return float(np.mean(r))


def modular_r_stat_from_rho(rho):
    dim = rho.shape[0]
    if dim > MAX_DIM_FOR_MODULAR:
        return np.nan
    rho_reg = rho + EPS * np.eye(dim)
    vals, vecs = np.linalg.eigh(rho_reg)
    vals = np.clip(vals, EPS, None)
    log_rho = vecs @ np.diag(np.log(vals)) @ vecs.conjugate().T
    K = -log_rho
    K = 0.5 * (K + K.conjugate().T)
    return r_stat_from_hermitian(K)


def scramble_global(psi):
    for _ in range(GLOBAL_SCRAMBLE_LAYERS):
        for _ in range(GLOBAL_GATES_PER_LAYER):
            q1, q2 = rng.choice(N, size=2, replace=False)
            psi = apply_two_qubit_gate(psi, random_two_qubit_unitary(rng), int(q1), int(q2), N)
    return psi


def compute_MI_matrix(psi):
    S1 = np.zeros(N, dtype=float)
    for i in range(N):
        rho_i = reduced_density_matrix(psi, [i], N)
        S1[i] = von_neumann_entropy(rho_i)

    W = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            rho_ij = reduced_density_matrix(psi, [i, j], N)
            Sij = von_neumann_entropy(rho_ij)
            Iij = S1[i] + S1[j] - Sij
            W[i, j] = W[j, i] = max(0.0, float(Iij))

    np.fill_diagonal(W, 0.0)
    upper = W[np.triu_indices(N, 1)]
    q99 = float(np.quantile(upper, 0.99)) if upper.size else 1.0
    scale = max(q99, 1e-12)
    Wn = np.clip(W / scale, 0.0, 1.0)
    return Wn


def adjacency_by_density(W, rho_target):
    W = _validate_W(W)
    n = W.shape[0]
    upper = W[np.triu_indices(n, 1)]
    if upper.size == 0:
        return np.zeros_like(W)
    thr = float(np.quantile(upper, 1.0 - rho_target))
    A = np.where(W >= thr, W, 0.0)
    np.fill_diagonal(A, 0.0)
    return np.maximum(A, A.T)


def graph_stats(A):
    A = _validate_W(A)
    n = A.shape[0]
    deg = (A > 0).sum(axis=1)
    dens = np.count_nonzero(np.triu(A, 1)) * 2 / (n * (n - 1))
    return {
        "density": float(dens),
        "min_degree": int(deg.min()) if deg.size else 0,
        "mean_degree": float(deg.mean()) if deg.size else 0.0,
        "max_degree": int(deg.max()) if deg.size else 0,
    }


def region_cut_internal(A, region):
    A = _validate_W(A)
    n = A.shape[0]
    region = sorted(region)
    outside = [i for i in range(n) if i not in region]
    cut = float(np.sum(A[np.ix_(region, outside)])) if outside else 0.0
    internal = float(np.triu(A[np.ix_(region, region)], 1).sum())
    return cut, internal


def region_conductance(A, region):
    A = _validate_W(A)
    n = A.shape[0]
    region = sorted(region)
    outside = [i for i in range(n) if i not in region]
    cut, _ = region_cut_internal(A, region)
    deg = A.sum(axis=1)
    volA = float(deg[region].sum()) if region else 0.0
    volB = float(deg[outside].sum()) if outside else 0.0
    denom = max(min(volA, volB), 1e-12)
    return float(cut / denom)


def region_metrics(psi, A, region):
    rho = reduced_density_matrix(psi, region, N)
    S = von_neumann_entropy(rho)
    cut, internal = region_cut_internal(A, region)
    phi = region_conductance(A, region)
    r = modular_r_stat_from_rho(rho)
    return S, cut, internal, phi, r


def normalize(x, mu, sig):
    return 0.0 if sig < 1e-12 else (x - mu) / sig


def find_bh_region(psi, A):
    S_pool, phi_pool, I_pool = [], [], []
    for _ in range(POOL_FOR_NORM):
        reg = sorted(rng.choice(N, size=A_SIZE, replace=False).tolist())
        S, _, internal, phi, _ = region_metrics(psi, A, reg)
        S_pool.append(S)
        phi_pool.append(phi)
        I_pool.append(internal)

    muS, sigS = float(np.mean(S_pool)), float(np.std(S_pool))
    muP, sigP = float(np.mean(phi_pool)), float(np.std(phi_pool))
    muI, sigI = float(np.mean(I_pool)), float(np.std(I_pool))

    def score(reg):
        S, _, internal, phi, _ = region_metrics(psi, A, reg)
        return (W_S * normalize(S, muS, sigS)
                - W_PHI * normalize(phi, muP, sigP)
                + W_INT * normalize(internal, muI, sigI))

    best = None
    best_sc = -1e18
    for _ in range(N_CANDIDATES):
        reg = sorted(rng.choice(N, size=A_SIZE, replace=False).tolist())
        sc = score(reg)
        if sc > best_sc:
            best_sc, best = sc, reg

    cur = best[:]
    cur_sc = best_sc
    for _ in range(HILL_STEPS):
        out = int(rng.choice(cur))
        outside = [i for i in range(N) if i not in cur]
        inn = int(rng.choice(outside))
        cand = sorted([x for x in cur if x != out] + [inn])
        sc = score(cand)
        if sc > cur_sc:
            cur, cur_sc = cand, sc
        else:
            if rng.random() < np.exp((sc - cur_sc) / max(TEMPERATURE, 1e-6)):
                cur, cur_sc = cand, sc

    return cur


def percentile(samples, value):
    s = np.array(samples, float)
    s = s[np.isfinite(s)]
    return np.nan if len(s) == 0 or not np.isfinite(value) else 100.0 * float(np.mean(s <= value))


def hist_with_marker(ax, data, value, title, xlabel, label_value):
    x = [t for t in data if np.isfinite(t)]
    ax.hist(x, bins=35, density=True, alpha=0.75)
    ax.axvline(value, linewidth=2.5, label=label_value)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend()
    ax.grid(True, alpha=0.3)


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def neighbors_of(region, n, flips=1, n_samples=80, seed=0):
    rr = np.random.default_rng(seed)
    region = sorted(region)
    S = set(region)
    out = []
    for _ in range(n_samples):
        cur = set(S)
        for __ in range(flips):
            out_node = rr.choice(list(cur))
            cur.remove(int(out_node))
            outside = [i for i in range(n) if i not in cur]
            in_node = rr.choice(outside)
            cur.add(int(in_node))
        out.append(sorted(cur))
    return out


# =========================
# MAIN
# =========================
def main():
    print("=== Community availability ===")
    print(f"networkx: {'OK' if HAS_NX else 'NO'}")
    print(f"python-louvain: {'OK' if HAS_LOUVAIN else 'NO'}")
    if HAS_NX and not HAS_LOUVAIN:
        print("Note: Louvain non dispo -> fallback greedy_modularity_communities (networkx).")
    if not HAS_NX:
        print("Note: pas de networkx -> pas de seeds communautés (mode mixed gardera spectral+random).")

    # 1) état initial |0..0>
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0

    # 2) scrambling global
    psi = scramble_global(psi)

    # 3) MI matrix W (normalisée robustement) puis adjacency A (densité fixée)
    W = compute_MI_matrix(psi)
    A = adjacency_by_density(W, RHO_TARGET)

    print("\n=== Graph stats (A, density-fixed) ===")
    print(graph_stats(A))

    # 4) BH-like region
    bh = find_bh_region(psi, A)
    S_bh, cut_bh, int_bh, phi_bh, r_bh = region_metrics(psi, A, bh)

    print("\n=== BH-like candidate (N=16) ===")
    print("Region:", bh)
    print(f"S={S_bh:.4f} | cut={cut_bh:.4f} | internal={int_bh:.4f} | phi={phi_bh:.4f} | r={r_bh:.4f}")

    # 5) HIE classic on SAME A
    print("\n=== HIE detection (classic spectral hierarchy, on A) ===")

    def entropy_fn(region):
        if len(region) == 0:
            return 0.0
        rho = reduced_density_matrix(psi, region, N)
        return von_neumann_entropy(rho)

    out_hie = detect_hie(
        A,
        max_depth=HIE_MAX_DEPTH,
        min_size=HIE_MIN_SIZE,
        entropy_fn=entropy_fn,
        top_k=HIE_TOP_K
    )

    print(f"Nombre de régions candidates: {out_hie['n_regions']}")
    print(f"\nTop {HIE_TOP_K} candidats HIE:")
    print(f"{'Rank':<6}{'Nodes':<44}{'Size':<6}{'Phi':<8}{'Cohesion':<12}{'S(A)':<10}{'Score':<8}")
    print("-" * 96)
    for i, c in enumerate(out_hie['top']):
        nodes_str = str(list(c.nodes))
        if len(nodes_str) > 42:
            nodes_str = nodes_str[:39] + "...]"
        s_a = f"{c.S_A:.2f}" if c.S_A is not None else "N/A"
        print(f"{i+1:<6}{nodes_str:<44}{c.size:<6}{c.conductance:<8.4f}{c.cohesion:<12.4f}{s_a:<10}{c.score:<8.2f}")

    bh_tuple = tuple(sorted(bh))
    bh_in_top = any(tuple(c.nodes) == bh_tuple for c in out_hie["top"])
    print(f"\nBH-like region in HIE top-{HIE_TOP_K}: {'✓ OUI' if bh_in_top else '✗ Non'}")

    # 6) HIE fixed-size (k=8) with Louvain seeds + neighbors
    forced = [bh] + neighbors_of(bh, N, flips=1, n_samples=120, seed=SEED+777)
    
    print("\n=== HIE fixed-size detection (k=8) with Louvain/community seeds ===")
    print(f"Mode: CLASSIC (w_phi=1.0)")
    
    out_hie8_classic = detect_hie_fixed_size(
        A,
        k=A_SIZE,
        entropy_fn=entropy_fn,
        top_k=10,
        n_seeds=HIE8_SEEDS,
        refine_steps=HIE8_REFINE_STEPS,
        refine_T=HIE8_REFINE_T,
        seed_mode=HIE8_SEED_MODE,
        community_method=HIE8_COMM_METHOD,
        forced_seeds=forced,
        seed=SEED,
        mode="classic",
        w_phi=1.0,
    )

    print(f"Nombre de régions candidates (k=8): {out_hie8_classic['n_regions']}")
    
    # Distribution phi pour analyse
    all_phi = [c.conductance for c in out_hie8_classic["all"]]
    print(f"Distribution φ: min={min(all_phi):.4f}, max={max(all_phi):.4f}, mean={np.mean(all_phi):.4f}")
    print(f"BH φ = {phi_bh:.4f} (doit être proche du min pour être un bottleneck)")

    print(f"\nTop 10 candidats HIE (k=8) - Mode CLASSIC:")
    print(f"{'Rank':<6}{'Nodes':<44}{'Size':<6}{'Phi':<8}{'Cohesion':<12}{'S(A)':<10}{'Score':<8}")
    print("-" * 96)

    for i, c in enumerate(out_hie8_classic['top']):
        nodes_str = str(list(c.nodes))
        if len(nodes_str) > 42:
            nodes_str = nodes_str[:39] + "...]"
        s_a = f"{c.S_A:.2f}" if c.S_A is not None else "N/A"
        marker = " <-- BH" if tuple(c.nodes) == bh_tuple else ""
        print(f"{i+1:<6}{nodes_str:<44}{c.size:<6}{c.conductance:<8.4f}{c.cohesion:<12.4f}{s_a:<10}{c.score:<8.2f}{marker}")

    # (A) RANG EXACT DE BH DANS LA LISTE
    rank_bh_classic = None
    for i, c in enumerate(out_hie8_classic["all"], start=1):
        if tuple(c.nodes) == bh_tuple:
            rank_bh_classic = i
            break
    print(f"\n>>> Rang BH dans HIE(k=8) CLASSIC: {rank_bh_classic} / {len(out_hie8_classic['all'])}")

    bh_in_hie8_classic = bh_tuple in set(out_hie8_classic["regions"])
    print(f"BH-like region present in HIE_fixed_size candidates: {'✓ OUI' if bh_in_hie8_classic else '✗ Non'}")

    # 6b) MODE HORIZON-AWARE avec poids fort sur phi
    print(f"\n=== HIE fixed-size detection (k=8) - Mode HORIZON-AWARE ===")
    print(f"Paramètres: w_phi={W_PHI_HORIZON}, w_S={W_S_HORIZON}, w_cut={W_CUT_HORIZON}")
    
    out_hie8_horizon = detect_hie_fixed_size(
        A,
        k=A_SIZE,
        entropy_fn=entropy_fn,
        top_k=10,
        n_seeds=HIE8_SEEDS,
        refine_steps=HIE8_REFINE_STEPS,
        refine_T=HIE8_REFINE_T,
        seed_mode=HIE8_SEED_MODE,
        community_method=HIE8_COMM_METHOD,
        forced_seeds=forced,
        seed=SEED,
        mode="horizon",
        w_phi=W_PHI_HORIZON,
        w_S=W_S_HORIZON,
        w_cut=W_CUT_HORIZON,
    )

    print(f"\nTop 10 candidats HIE (k=8) - Mode HORIZON:")
    print(f"{'Rank':<6}{'Nodes':<44}{'Size':<6}{'Phi':<8}{'Cohesion':<12}{'S(A)':<10}{'Score':<8}")
    print("-" * 96)

    for i, c in enumerate(out_hie8_horizon['top']):
        nodes_str = str(list(c.nodes))
        if len(nodes_str) > 42:
            nodes_str = nodes_str[:39] + "...]"
        s_a = f"{c.S_A:.2f}" if c.S_A is not None else "N/A"
        marker = " <-- BH" if tuple(c.nodes) == bh_tuple else ""
        print(f"{i+1:<6}{nodes_str:<44}{c.size:<6}{c.conductance:<8.4f}{c.cohesion:<12.4f}{s_a:<10}{c.score:<8.2f}{marker}")

    # (A) RANG EXACT DE BH EN MODE HORIZON
    rank_bh_horizon = None
    for i, c in enumerate(out_hie8_horizon["all"], start=1):
        if tuple(c.nodes) == bh_tuple:
            rank_bh_horizon = i
            break
    print(f"\n>>> Rang BH dans HIE(k=8) HORIZON: {rank_bh_horizon} / {len(out_hie8_horizon['all'])}")

    # Comparaison des deux modes
    print(f"\n=== COMPARAISON DES MODES ===")
    print(f"Mode CLASSIC:  BH rang {rank_bh_classic}/{len(out_hie8_classic['all'])}")
    print(f"Mode HORIZON:  BH rang {rank_bh_horizon}/{len(out_hie8_horizon['all'])}")
    improvement = rank_bh_classic - rank_bh_horizon if rank_bh_horizon and rank_bh_classic else 0
    print(f"Amélioration: {improvement} places" if improvement > 0 else f"Dégradation: {-improvement} places" if improvement < 0 else "Identique")

    if len(out_hie8_horizon["top"]) > 0:
        print("\nChevauchements Jaccard top-5 HIE HORIZON vs BH:")
        for i, c in enumerate(out_hie8_horizon["top"][:5]):
            jac = jaccard(bh, c.nodes)
            print(f"  Rank {i+1}: {jac:.3f} - {list(c.nodes)}")

    # 7) Benchmark random (same size)
    S_rand, cut_rand, int_rand, phi_rand, r_rand = [], [], [], [], []
    for _ in range(N_RANDOM_REGIONS):
        reg = sorted(rng.choice(N, size=A_SIZE, replace=False).tolist())
        S, cut, internal, phi, r = region_metrics(psi, A, reg)
        S_rand.append(S)
        cut_rand.append(cut)
        int_rand.append(internal)
        phi_rand.append(phi)
        r_rand.append(r)

    pS = percentile(S_rand, S_bh)
    pC = percentile(cut_rand, cut_bh)
    pI = percentile(int_rand, int_bh)
    pP = percentile(phi_rand, phi_bh)
    pr = percentile(r_rand, r_bh)

    print("\n=== Percentiles BH-like vs random regions (same size) ===")
    print(f"S: {pS:.1f}% | cut: {pC:.1f}% | internal: {pI:.1f}% | phi: {pP:.1f}% | r: {pr:.1f}%")
    print("Note: horizon-like bottleneck => phi should be LOW (so percentile should be SMALL).")

    # 8) Figures
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))

    # Panel 1: Matrix W avec overlays
    ax = axes[0, 0]
    im = ax.imshow(W, aspect='auto')
    ax.set_title("Matrice MI normalisée $W_{ij}$ (visualisation)", fontsize=11, fontweight='bold')

    # Overlay BH (cyan)
    for i in bh:
        for j in bh:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor='cyan', linewidth=2, alpha=0.7))

    # Overlay HIE8 HORIZON top1 (red)
    if len(out_hie8_horizon["top"]) > 0:
        hie_top_nodes = list(out_hie8_horizon["top"][0].nodes)
        for i in hie_top_nodes:
            for j in hie_top_nodes:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor='red', linewidth=2, alpha=0.7))

    ax.set_xlabel("Qubit j")
    ax.set_ylabel("Qubit i")
    plt.colorbar(im, ax=ax)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor='none', edgecolor='cyan', linewidth=2, label='BH-like'),
        Patch(facecolor='none', edgecolor='red', linewidth=2, label='HIE8 HORIZON Top-1'),
    ], loc='upper right')

    hist_with_marker(axes[0, 1], S_rand, S_bh, "Benchmark S(A)", "S(A) bits", f"BH-like (p={pS:.0f}%)")
    hist_with_marker(axes[0, 2], phi_rand, phi_bh, "Benchmark conductance φ(A)", "φ(A)", f"BH-like (p={pP:.0f}%)")
    hist_with_marker(axes[1, 0], int_rand, int_bh, "Benchmark internal(A)", "internal(A)", f"BH-like (p={pI:.0f}%)")
    hist_with_marker(axes[1, 1], cut_rand, cut_bh, "Benchmark cut(A)", "cut(A)", f"BH-like (p={pC:.0f}%)")

    # Panel 6: Comparaison CLASSIC vs HORIZON pour BH
    ax = axes[1, 2]
    categories = ['φ (low=good)', 'S(A)', 'internal', 'cut (low=good)']
    bh_metrics = np.array([phi_bh, S_bh, int_bh, cut_bh])
    
    # Normaliser pour visualisation (inverser phi et cut pour que "high" = "good")
    display_vals = np.array([1.0/phi_bh, S_bh, int_bh, 1.0/cut_bh])
    
    x = np.arange(len(categories))
    width = 0.6
    ax.bar(x, display_vals / max(display_vals), width, color=['red', 'blue', 'green', 'orange'], alpha=0.7)
    ax.set_ylabel('Normalized (BH profile)')
    ax.set_title('BH-like region profile\n(high = desirable for horizon)')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=15, ha='right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "comparison_bh_hie8.png"), dpi=170)
    plt.close()

    # Distribution des rangs
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Rang de BH dans les deux modes
    modes = ['CLASSIC', 'HORIZON']
    ranks = [rank_bh_classic if rank_bh_classic else len(out_hie8_classic['all']), 
             rank_bh_horizon if rank_bh_horizon else len(out_hie8_horizon['all'])]
    
    axes[0].bar(modes, ranks, color=['gray', 'green'], alpha=0.7)
    axes[0].set_ylabel('Rang de BH (plus petit = mieux)')
    axes[0].set_title('Comparaison du rang BH: CLASSIC vs HORIZON')
    axes[0].axhline(y=1, color='r', linestyle='--', label='Top-1')
    axes[0].legend()
    
    # Distribution des phi pour tous les candidats
    axes[1].hist([c.conductance for c in out_hie8_horizon["all"]], bins=20, alpha=0.5, label='HORIZON candidates')
    axes[1].axvline(phi_bh, color='cyan', linewidth=3, label=f'BH φ={phi_bh:.3f}')
    axes[1].set_xlabel('Conductance φ')
    axes[1].set_ylabel('Nombre de régions')
    axes[1].set_title('Distribution de φ (HORIZON mode)')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "rank_comparison.png"), dpi=170)
    plt.close()

    # Percentiles bar
    plt.figure(figsize=(9, 4))
    plt.bar(["S", "phi", "internal", "cut", "r"], [pS, pP, pI, pC, pr])
    plt.axhline(95, linestyle="--", label="95%")
    plt.axhline(5, linestyle="--", color='red', label="5% (bottleneck zone)")
    plt.ylim(0, 100)
    plt.title("Percentiles BH-like vs random regions (N=16)")
    plt.ylabel("Percentile (%)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "percentiles.png"), dpi=170)
    plt.close()

    # 9) JSON report
    report = {
        "N": N,
        "seed": SEED,
        "A_SIZE": A_SIZE,
        "community": {
            "has_networkx": bool(HAS_NX),
            "has_louvain": bool(HAS_LOUVAIN),
            "seed_mode": HIE8_SEED_MODE,
            "community_method": HIE8_COMM_METHOD,
            "n_seeds": HIE8_SEEDS,
            "refine_steps": HIE8_REFINE_STEPS,
            "refine_T": HIE8_REFINE_T,
        },
        "graph_stats_A": graph_stats(A),
        "bh_region": bh,
        "bh_metrics": {"S": S_bh, "cut": cut_bh, "internal": int_bh, "phi": phi_bh, "r": r_bh},
        "bh_percentiles": {"S": pS, "cut": pC, "internal": pI, "phi": pP, "r": pr},
        "hie_classic": {
            "n_candidates": out_hie["n_regions"],
            "top_candidates": [
                {
                    "nodes": list(c.nodes),
                    "size": c.size,
                    "phi": c.conductance,
                    "cohesion": c.cohesion,
                    "I_int": c.I_int,
                    "cut": c.A_cut,
                    "S_A": c.S_A,
                    "score": c.score,
                } for c in out_hie["top"]
            ],
            "bh_in_top": bool(bh_in_top),
        },
        "hie_fixed_size_classic": {
            "n_candidates": out_hie8_classic["n_regions"],
            "bh_rank": rank_bh_classic,
            "top_candidates": [
                {
                    "nodes": list(c.nodes),
                    "size": c.size,
                    "phi": c.conductance,
                    "cohesion": c.cohesion,
                    "I_int": c.I_int,
                    "cut": c.A_cut,
                    "S_A": c.S_A,
                    "score": c.score,
                } for c in out_hie8_classic["top"]
            ],
            "bh_in_candidates": bool(bh_in_hie8_classic),
        },
        "hie_fixed_size_horizon": {
            "n_candidates": out_hie8_horizon["n_regions"],
            "bh_rank": rank_bh_horizon,
            "params": {"w_phi": W_PHI_HORIZON, "w_S": W_S_HORIZON, "w_cut": W_CUT_HORIZON},
            "top_candidates": [
                {
                    "nodes": list(c.nodes),
                    "size": c.size,
                    "phi": c.conductance,
                    "cohesion": c.cohesion,
                    "I_int": c.I_int,
                    "cut": c.A_cut,
                    "S_A": c.S_A,
                    "score": c.score,
                } for c in out_hie8_horizon["top"]
            ],
            "bh_in_candidates": bh_tuple in set(out_hie8_horizon["regions"]),
        },
    }

    with open(os.path.join(OUTDIR, "report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print("\nSorties:", OUTDIR)
    print("- report.json")
    print("- comparison_bh_hie8.png")
    print("- rank_comparison.png")
    print("- percentiles.png")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
BUp Ollivier-Ricci Local Response Batch v1.3
-------------------------------------------------------------------------------
Objectif
--------
Valider statistiquement, sur plusieurs seeds, si une perturbation
informationnelle locale modifie la courbure d'Ollivier-Ricci du graphe
corrélation/intrication, avec séparation entre :

1) réponse totale (graphe libre),
2) réponse à connectivité fixée (graphes gelés),
3) contrôle par arêtes communes.

Nouveautés v1.3
---------------
- Batch multi-seeds (--n-seeds, --seed0)
- Agrégation statistique par mode
- CSV par seed + JSON résumé
- Boxplots / scatter agrégés pour juger la robustesse du signal
- Conservation des mêmes modes que v1.2 :
  free, frozen_baseline, frozen_perturbed, frozen_intersection, frozen_union

Sorties principales
-------------------
- report_batch.json
- per_seed_results.csv
- summary_by_mode.csv
- boxplot_delta_edge_mean.png
- boxplot_delta_near_minus_far.png
- boxplot_common_edge_delta_mean.png
- scatter_near_vs_far_<mode>.png

Dépendances minimales
---------------------
pip install numpy scipy networkx matplotlib

Optionnel pour backend graphricci
---------------------------------
pip install GraphRicciCurvature pot
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

try:
    from scipy.optimize import linprog
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except Exception:  # pragma: no cover
    HAS_SCIPY = False
    linprog = None
    scipy_stats = None

try:
    from GraphRicciCurvature.OllivierRicci import OllivierRicci

    HAS_GRAPH_RICCI = True
except Exception:  # pragma: no cover
    HAS_GRAPH_RICCI = False
    OllivierRicci = None


# -----------------------------------------------------------------------------
# Utilitaires
# -----------------------------------------------------------------------------


def parse_region(text: str) -> List[int]:
    text = text.strip()
    if not text:
        return []
    return [int(x) for x in text.split(",") if x.strip()]



def save_json(obj: Dict, path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)



def ring_distance(a: int, b: int, n: int) -> int:
    d = abs(a - b)
    return min(d, n - d)



def distance_to_region_on_ring(node: int, region: Sequence[int], n: int) -> int:
    return min(ring_distance(node, r, n) for r in region)



def _safe_mean(x: List[float]) -> float:
    return float(np.mean(x)) if x else float("nan")



def _safe_stat(x: List[float], fn) -> float:
    return float(fn(x)) if x else float("nan")


# -----------------------------------------------------------------------------
# Portes quantiques et simulateur d'état pur
# -----------------------------------------------------------------------------


def random_su2(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(z)
    d = np.diag(r)
    phase = d / np.where(np.abs(d) > 0, np.abs(d), 1.0)
    q = q @ np.diag(np.conjugate(phase))
    det = np.linalg.det(q)
    q = q / np.sqrt(det)
    return q.astype(np.complex128)



def partial_iswap(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, c, 1j * s, 0.0],
            [0.0, 1j * s, c, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=np.complex128,
    )



def apply_one_qubit_gate(psi: np.ndarray, gate: np.ndarray, q: int, n: int) -> np.ndarray:
    tensor = psi.reshape([2] * n)
    axes = [q] + [i for i in range(n) if i != q]
    inv_axes = np.argsort(axes)
    temp = np.transpose(tensor, axes).reshape(2, -1)
    temp = gate @ temp
    temp = temp.reshape([2] + [2] * (n - 1))
    tensor = np.transpose(temp, inv_axes)
    return tensor.reshape(-1)



def apply_two_qubit_gate(psi: np.ndarray, gate: np.ndarray, q1: int, q2: int, n: int) -> np.ndarray:
    if q1 == q2:
        raise ValueError("q1 et q2 doivent être distincts")
    if q1 > q2:
        q1, q2 = q2, q1
    tensor = psi.reshape([2] * n)
    axes = [q1, q2] + [i for i in range(n) if i not in (q1, q2)]
    inv_axes = np.argsort(axes)
    temp = np.transpose(tensor, axes).reshape(4, -1)
    temp = gate @ temp
    temp = temp.reshape([2, 2] + [2] * (n - 2))
    tensor = np.transpose(temp, inv_axes)
    return tensor.reshape(-1)


@dataclass
class LayerSpec:
    one_qubit_gates: List[np.ndarray]
    pairs: List[Tuple[int, int, float]]



def sample_disjoint_pairs(
    n: int,
    gates_per_layer: int,
    lam: float,
    rng: np.random.Generator,
) -> List[Tuple[int, int]]:
    chosen: List[Tuple[int, int]] = []
    used: set[int] = set()
    attempts = 0
    max_attempts = 50 * max(1, gates_per_layer)

    while len(chosen) < gates_per_layer and attempts < max_attempts:
        attempts += 1
        if rng.random() < lam:
            q1, q2 = sorted(rng.choice(n, size=2, replace=False).tolist())
        else:
            q1 = int(rng.integers(0, n))
            q2 = (q1 + 1) % n
            q1, q2 = sorted((q1, q2))

        if q1 in used or q2 in used:
            continue
        used.add(q1)
        used.add(q2)
        chosen.append((q1, q2))

    return chosen



def make_schedule(
    n: int,
    layers: int,
    gates_per_layer: int,
    lam: float,
    theta_base: float,
    seed: int,
) -> List[LayerSpec]:
    rng = np.random.default_rng(seed)
    schedule: List[LayerSpec] = []
    for _ in range(layers):
        one_q = [random_su2(rng) for _ in range(n)]
        pairs = sample_disjoint_pairs(n, gates_per_layer, lam, rng)
        pair_specs: List[Tuple[int, int, float]] = []
        for q1, q2 in pairs:
            theta = float(np.clip(theta_base + 0.15 * rng.normal(), 0.02, np.pi / 2 - 0.02))
            pair_specs.append((q1, q2, theta))
        schedule.append(LayerSpec(one_qubit_gates=one_q, pairs=pair_specs))
    return schedule



def run_schedule(
    n: int,
    schedule: Sequence[LayerSpec],
    source_region: Sequence[int],
    defect_radius: int,
    defect_strength: float,
    defect_mode: str,
) -> np.ndarray:
    psi = np.zeros(2**n, dtype=np.complex128)
    psi[0] = 1.0

    for layer in schedule:
        for q, gate in enumerate(layer.one_qubit_gates):
            psi = apply_one_qubit_gate(psi, gate, q, n)

        for q1, q2, theta in layer.pairs:
            pair_dist = min(
                distance_to_region_on_ring(q1, source_region, n),
                distance_to_region_on_ring(q2, source_region, n),
            )
            theta_eff = theta
            if pair_dist <= defect_radius:
                if defect_mode == "enhance":
                    theta_eff = theta * (1.0 + defect_strength)
                elif defect_mode == "suppress":
                    theta_eff = theta * max(0.0, 1.0 - defect_strength)
                else:
                    raise ValueError(f"defect_mode inconnu: {defect_mode}")
            theta_eff = float(np.clip(theta_eff, 0.0, np.pi / 2))
            psi = apply_two_qubit_gate(psi, partial_iswap(theta_eff), q1, q2, n)

    norm = np.linalg.norm(psi)
    if norm == 0:
        raise RuntimeError("État final nul, ce qui ne devrait pas arriver.")
    return psi / norm


# -----------------------------------------------------------------------------
# Entropies et matrice d'information mutuelle
# -----------------------------------------------------------------------------


def subsystem_entropy(psi: np.ndarray, subset: Sequence[int], n: int, eps: float = 1e-15) -> float:
    subset = tuple(sorted(subset))
    rest = tuple(i for i in range(n) if i not in subset)
    perm = subset + rest
    tensor = psi.reshape([2] * n)
    mat = np.transpose(tensor, perm).reshape(2 ** len(subset), -1)
    svals = np.linalg.svd(mat, compute_uv=False)
    probs = np.maximum((svals**2).real, 0.0)
    probs = probs / max(probs.sum(), eps)
    probs = probs[probs > eps]
    return float(-(probs * np.log2(probs)).sum())



def mutual_information_matrix(psi: np.ndarray, n: int) -> Tuple[np.ndarray, List[float]]:
    s1 = [subsystem_entropy(psi, [i], n) for i in range(n)]
    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            s2 = subsystem_entropy(psi, [i, j], n)
            mi = max(0.0, s1[i] + s1[j] - s2)
            W[i, j] = W[j, i] = mi
    return W, s1


# -----------------------------------------------------------------------------
# Construction des graphes
# -----------------------------------------------------------------------------


def build_graph_from_edges_and_weights(
    edges: Iterable[Tuple[int, int]],
    W: np.ndarray,
    eps_length: float = 1e-9,
) -> nx.Graph:
    n = W.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))

    selected_pairs = []
    for i, j in sorted({tuple(sorted((int(a), int(b)))) for a, b in edges}):
        sim = float(max(W[i, j], 0.0))
        selected_pairs.append((i, j, sim))

    smax = max((s for _, _, s in selected_pairs), default=1.0)
    smax = max(smax, eps_length)

    for i, j, sim in selected_pairs:
        sim_norm = sim / smax
        length = 1.0 / max(sim_norm, eps_length)
        G.add_edge(i, j, similarity=sim, similarity_norm=sim_norm, length=length)

    return G



def build_connected_threshold_graph(W: np.ndarray, density: float, eps_length: float = 1e-9) -> nx.Graph:
    n = W.shape[0]
    pairs: List[Tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, float(max(W[i, j], 0.0))))

    target_edges = int(round(density * n * (n - 1) / 2.0))
    target_edges = min(max(target_edges, n - 1), n * (n - 1) // 2)

    K = nx.Graph()
    K.add_nodes_from(range(n))
    for i, j, sim in pairs:
        K.add_edge(i, j, similarity=sim)

    mst = nx.maximum_spanning_tree(K, weight="similarity")
    selected = {tuple(sorted((u, v))) for u, v in mst.edges()}

    remaining = sorted(pairs, key=lambda x: (x[2], -abs(x[0] - x[1])), reverse=True)
    for i, j, _ in remaining:
        if len(selected) >= target_edges:
            break
        e = tuple(sorted((i, j)))
        if e not in selected:
            selected.add(e)

    return build_graph_from_edges_and_weights(selected, W, eps_length=eps_length)



def edge_set(G: nx.Graph) -> set[Tuple[int, int]]:
    return {tuple(sorted((int(u), int(v)))) for u, v in G.edges()}


# -----------------------------------------------------------------------------
# Courbure Ollivier-Ricci
# -----------------------------------------------------------------------------


def local_measure(G: nx.Graph, node: int, alpha: float) -> Tuple[List[int], np.ndarray]:
    nbrs = list(G.neighbors(node))
    support = [node] + nbrs
    mass = np.zeros(len(support), dtype=float)
    mass[0] = alpha

    if nbrs:
        weights = []
        for nbr in nbrs:
            data = G[node][nbr]
            sim = float(data.get("similarity_norm", 0.0))
            if sim <= 0.0:
                length = float(data.get("length", 1.0))
                sim = 1.0 / max(length, 1e-12)
            weights.append(sim)
        weights = np.asarray(weights, dtype=float)
        if weights.sum() <= 0:
            weights = np.ones_like(weights)
        weights = weights / weights.sum()
        mass[1:] = (1.0 - alpha) * weights
    else:
        mass[0] = 1.0

    return support, mass



def wasserstein_distance_lp(mu: np.ndarray, nu: np.ndarray, cost: np.ndarray) -> float:
    if not HAS_SCIPY:
        raise RuntimeError("scipy est requis pour le backend internal (linprog).")

    m, n = cost.shape
    c = cost.reshape(-1)

    A_eq = []
    b_eq = []

    for i in range(m):
        row = np.zeros(m * n, dtype=float)
        row[i * n : (i + 1) * n] = 1.0
        A_eq.append(row)
        b_eq.append(mu[i])

    for j in range(n):
        col = np.zeros(m * n, dtype=float)
        col[j::n] = 1.0
        A_eq.append(col)
        b_eq.append(nu[j])

    res = linprog(
        c=c,
        A_eq=np.asarray(A_eq),
        b_eq=np.asarray(b_eq),
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        raise RuntimeError(f"Échec du transport optimal: {res.message}")
    return float(res.fun)



def compute_ollivier_ricci_internal(G: nx.Graph, alpha: float) -> nx.Graph:
    G = G.copy()
    if not nx.is_connected(G):
        raise RuntimeError("Le backend internal requiert un graphe connexe.")
    apsp = dict(nx.all_pairs_dijkstra_path_length(G, weight="length"))

    for u, v in G.edges():
        support_u, mu = local_measure(G, u, alpha=alpha)
        support_v, nu = local_measure(G, v, alpha=alpha)

        cost = np.zeros((len(support_u), len(support_v)), dtype=float)
        for i, a in enumerate(support_u):
            for j, b in enumerate(support_v):
                cost[i, j] = float(apsp[a][b])

        w1 = wasserstein_distance_lp(mu, nu, cost)
        d_uv = float(G[u][v]["length"])
        kappa = 1.0 - (w1 / max(d_uv, 1e-12))
        G[u][v]["ricciCurvature"] = float(kappa)

    for n in G.nodes():
        vals = [float(G[n][nbr]["ricciCurvature"]) for nbr in G.neighbors(n)]
        G.nodes[n]["ricciCurvature"] = float(np.mean(vals)) if vals else float("nan")

    return G



def compute_ollivier_ricci_graphricci(
    G: nx.Graph,
    alpha: float,
    method: str,
    proc: int,
    verbose: str,
) -> nx.Graph:
    if not HAS_GRAPH_RICCI:
        raise RuntimeError(
            "GraphRicciCurvature indisponible. Installe: pip install GraphRicciCurvature pot"
        )
    orc = OllivierRicci(
        G,
        weight="length",
        alpha=alpha,
        method=method,
        proc=proc,
        verbose=verbose,
        shortest_path="all_pairs",
    )
    return orc.compute_ricci_curvature()



def compute_ollivier_ricci(
    G: nx.Graph,
    alpha: float,
    backend: str,
    method: str,
    proc: int,
    verbose: str,
) -> Tuple[nx.Graph, str]:
    if backend == "auto":
        backend = "graphricci" if HAS_GRAPH_RICCI else "internal"

    if backend == "graphricci":
        return compute_ollivier_ricci_graphricci(G, alpha, method, proc, verbose), backend
    if backend == "internal":
        return compute_ollivier_ricci_internal(G, alpha), backend
    raise ValueError(f"backend inconnu: {backend}")


# -----------------------------------------------------------------------------
# Analyse comparative baseline vs perturbé
# -----------------------------------------------------------------------------


def node_source_distances(G: nx.Graph, source_region: Sequence[int]) -> Dict[int, float]:
    dist = {int(n): math.inf for n in G.nodes()}
    for s in source_region:
        if s not in G:
            continue
        d = nx.single_source_shortest_path_length(G, s)
        for k, v in d.items():
            dist[int(k)] = min(dist[int(k)], float(v))
    return dist



def edge_source_distance(u: int, v: int, node_dist: Dict[int, float]) -> float:
    return min(float(node_dist.get(u, math.inf)), float(node_dist.get(v, math.inf)))



def summarize_graph_curvature(G: nx.Graph, node_dist: Dict[int, float], near_radius: int) -> Dict[str, float]:
    edge_vals = []
    edge_near = []
    edge_far = []
    for u, v, data in G.edges(data=True):
        rc = float(data.get("ricciCurvature", np.nan))
        if np.isnan(rc):
            continue
        edge_vals.append(rc)
        d = edge_source_distance(int(u), int(v), node_dist)
        if (not math.isinf(d)) and d <= near_radius:
            edge_near.append(rc)
        elif not math.isinf(d):
            edge_far.append(rc)

    node_vals = [float(G.nodes[n].get("ricciCurvature", np.nan)) for n in G.nodes()]
    node_vals = [x for x in node_vals if not np.isnan(x)]

    near_mean = _safe_mean(edge_near)
    far_mean = _safe_mean(edge_far)

    return {
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "edge_mean": _safe_stat(edge_vals, np.mean),
        "node_mean": _safe_stat(node_vals, np.mean),
        "near_edge_mean": near_mean,
        "far_edge_mean": far_mean,
        "near_minus_far": float(near_mean - far_mean) if (not np.isnan(near_mean) and not np.isnan(far_mean)) else float("nan"),
    }



def edge_dict(G: nx.Graph) -> Dict[Tuple[int, int], Dict[str, float]]:
    out: Dict[Tuple[int, int], Dict[str, float]] = {}
    for u, v, data in G.edges(data=True):
        e = tuple(sorted((int(u), int(v))))
        out[e] = {
            "ricciCurvature": float(data.get("ricciCurvature", np.nan)),
            "length": float(data.get("length", np.nan)),
            "similarity": float(data.get("similarity", np.nan)),
        }
    return out



def compare_graphs(
    G0: nx.Graph,
    G1: nx.Graph,
    source_region: Sequence[int],
    near_radius: int,
) -> Dict:
    node_dist = node_source_distances(G0, source_region)
    sum0 = summarize_graph_curvature(G0, node_dist, near_radius)
    sum1 = summarize_graph_curvature(G1, node_dist, near_radius)

    e0 = edge_dict(G0)
    e1 = edge_dict(G1)
    common_edges = sorted(set(e0) & set(e1))

    common_deltas = []
    near_common = []
    far_common = []
    for e in common_edges:
        rc0 = float(e0[e]["ricciCurvature"])
        rc1 = float(e1[e]["ricciCurvature"])
        if np.isnan(rc0) or np.isnan(rc1):
            continue
        delta = rc1 - rc0
        common_deltas.append(delta)
        d = edge_source_distance(e[0], e[1], node_dist)
        if (not math.isinf(d)) and d <= near_radius:
            near_common.append(delta)
        elif not math.isinf(d):
            far_common.append(delta)

    delta_edge_mean = float(sum1["edge_mean"] - sum0["edge_mean"]) if (not np.isnan(sum0["edge_mean"]) and not np.isnan(sum1["edge_mean"])) else float("nan")
    delta_node_mean = float(sum1["node_mean"] - sum0["node_mean"]) if (not np.isnan(sum0["node_mean"]) and not np.isnan(sum1["node_mean"])) else float("nan")
    delta_near_edge_mean = float(sum1["near_edge_mean"] - sum0["near_edge_mean"]) if (not np.isnan(sum0["near_edge_mean"]) and not np.isnan(sum1["near_edge_mean"])) else float("nan")
    delta_far_edge_mean = float(sum1["far_edge_mean"] - sum0["far_edge_mean"]) if (not np.isnan(sum0["far_edge_mean"]) and not np.isnan(sum1["far_edge_mean"])) else float("nan")
    delta_near_minus_far = float(sum1["near_minus_far"] - sum0["near_minus_far"]) if (not np.isnan(sum0["near_minus_far"]) and not np.isnan(sum1["near_minus_far"])) else float("nan")

    return {
        "baseline": sum0,
        "perturbed": sum1,
        "delta_summary": {
            "delta_edge_mean": delta_edge_mean,
            "delta_node_mean": delta_node_mean,
            "delta_near_edge_mean": delta_near_edge_mean,
            "delta_far_edge_mean": delta_far_edge_mean,
            "delta_near_minus_far": delta_near_minus_far,
            "common_edge_delta_mean": _safe_mean(common_deltas),
            "common_edge_delta_mean_near": _safe_mean(near_common),
            "common_edge_delta_mean_far": _safe_mean(far_common),
            "edge_jaccard": float(len(common_edges) / max(1, len(set(e0) | set(e1)))),
            "edges_only_baseline": int(len(set(e0) - set(e1))),
            "edges_only_perturbed": int(len(set(e1) - set(e0))),
            "n_common_edges": int(len(common_edges)),
        },
    }


# -----------------------------------------------------------------------------
# Modes de support de graphe
# -----------------------------------------------------------------------------


def make_mode_graphs(mode: str, G0_free: nx.Graph, G1_free: nx.Graph, W0: np.ndarray, W1: np.ndarray) -> Tuple[nx.Graph, nx.Graph, Dict[str, float]]:
    e0 = edge_set(G0_free)
    e1 = edge_set(G1_free)

    if mode == "free":
        return G0_free.copy(), G1_free.copy(), {
            "edge_count_baseline": float(len(e0)),
            "edge_count_perturbed": float(len(e1)),
            "edge_count_support": float(len(e0 | e1)),
            "connected_baseline": float(nx.is_connected(G0_free)),
            "connected_perturbed": float(nx.is_connected(G1_free)),
            "jaccard_free_graphs": float(len(e0 & e1) / max(1, len(e0 | e1))),
        }

    if mode == "frozen_baseline":
        support = e0
    elif mode == "frozen_perturbed":
        support = e1
    elif mode == "frozen_union":
        support = e0 | e1
    elif mode == "frozen_intersection":
        support = e0 & e1
    else:
        raise ValueError(f"Mode de graphe inconnu: {mode}")

    G0 = build_graph_from_edges_and_weights(support, W0)
    G1 = build_graph_from_edges_and_weights(support, W1)
    meta = {
        "edge_count_baseline": float(len(e0)),
        "edge_count_perturbed": float(len(e1)),
        "edge_count_support": float(len(support)),
        "jaccard_free_graphs": float(len(e0 & e1) / max(1, len(e0 | e1))),
        "connected_baseline": float(nx.is_connected(G0)) if G0.number_of_edges() > 0 else 0.0,
        "connected_perturbed": float(nx.is_connected(G1)) if G1.number_of_edges() > 0 else 0.0,
    }
    return G0, G1, meta



def requested_modes(graph_mode: str) -> List[str]:
    if graph_mode == "all":
        return ["free", "frozen_baseline", "frozen_perturbed", "frozen_intersection", "frozen_union"]
    return [graph_mode]


# -----------------------------------------------------------------------------
# Batch et statistiques
# -----------------------------------------------------------------------------


SUMMARY_METRICS = [
    "delta_edge_mean",
    "delta_node_mean",
    "delta_near_edge_mean",
    "delta_far_edge_mean",
    "delta_near_minus_far",
    "common_edge_delta_mean",
    "common_edge_delta_mean_near",
    "common_edge_delta_mean_far",
    "edge_jaccard",
    "edges_only_baseline",
    "edges_only_perturbed",
    "n_common_edges",
    "support_jaccard",
]



def build_states_for_seed(args: argparse.Namespace, seed: int, source_region: Sequence[int]) -> Tuple[np.ndarray, np.ndarray]:
    schedule = make_schedule(
        n=args.n,
        layers=args.layers,
        gates_per_layer=args.gates,
        lam=args.lam,
        theta_base=args.theta,
        seed=seed,
    )
    psi0 = run_schedule(
        n=args.n,
        schedule=schedule,
        source_region=source_region,
        defect_radius=-1,
        defect_strength=0.0,
        defect_mode=args.defect_mode,
    )
    psi1 = run_schedule(
        n=args.n,
        schedule=schedule,
        source_region=source_region,
        defect_radius=args.defect_radius,
        defect_strength=args.defect_strength,
        defect_mode=args.defect_mode,
    )
    return psi0, psi1



def compute_mode_analysis(
    mode: str,
    G0_free: nx.Graph,
    G1_free: nx.Graph,
    W0: np.ndarray,
    W1: np.ndarray,
    source_region: Sequence[int],
    args: argparse.Namespace,
) -> Tuple[Dict, str]:
    G0_raw, G1_raw, mode_meta = make_mode_graphs(mode, G0_free, G1_free, W0, W1)
    mode_report: Dict[str, object] = {
        "graph_mode": mode,
        "graph_meta": mode_meta,
        "status": "ok",
    }

    if G0_raw.number_of_edges() == 0 or G1_raw.number_of_edges() == 0:
        mode_report["status"] = "skipped_empty_graph"
        return mode_report, args.backend

    if (not nx.is_connected(G0_raw)) or (not nx.is_connected(G1_raw)):
        mode_report["status"] = "skipped_disconnected_graph"
        return mode_report, args.backend

    G0, backend_used = compute_ollivier_ricci(
        G0_raw,
        alpha=args.alpha,
        backend=args.backend,
        method=args.method,
        proc=args.proc,
        verbose=args.verbose,
    )
    G1, _ = compute_ollivier_ricci(
        G1_raw,
        alpha=args.alpha,
        backend=args.backend,
        method=args.method,
        proc=args.proc,
        verbose=args.verbose,
    )

    comparison = compare_graphs(G0, G1, source_region=source_region, near_radius=args.near_radius)
    mode_report["comparison"] = comparison
    return mode_report, backend_used



def row_from_mode_report(seed: int, mode: str, mode_report: Dict) -> Dict[str, float | int | str]:
    row: Dict[str, float | int | str] = {"seed": int(seed), "mode": mode, "status": str(mode_report.get("status", "?"))}
    graph_meta = mode_report.get("graph_meta", {}) or {}
    row["support_jaccard"] = float(graph_meta.get("jaccard_free_graphs", np.nan))

    if row["status"] != "ok":
        for metric in SUMMARY_METRICS:
            row.setdefault(metric, np.nan)
        return row

    delta = mode_report["comparison"]["delta_summary"]
    for metric in SUMMARY_METRICS:
        if metric == "support_jaccard":
            continue
        row[metric] = float(delta.get(metric, np.nan))
    return row



def one_sample_stats(values: List[float]) -> Dict[str, float | int | None]:
    arr = np.asarray([v for v in values if not np.isnan(v)], dtype=float)
    n = int(arr.size)
    out: Dict[str, float | int | None] = {
        "n": n,
        "mean": float(np.mean(arr)) if n else float("nan"),
        "std": float(np.std(arr, ddof=1)) if n >= 2 else float("nan"),
        "median": float(np.median(arr)) if n else float("nan"),
        "q25": float(np.quantile(arr, 0.25)) if n else float("nan"),
        "q75": float(np.quantile(arr, 0.75)) if n else float("nan"),
        "min": float(np.min(arr)) if n else float("nan"),
        "max": float(np.max(arr)) if n else float("nan"),
        "positive_frac": float(np.mean(arr > 0.0)) if n else float("nan"),
        "negative_frac": float(np.mean(arr < 0.0)) if n else float("nan"),
        "nonzero_frac": float(np.mean(arr != 0.0)) if n else float("nan"),
        "t_stat": None,
        "t_pvalue": None,
        "sign_pvalue": None,
    }
    if n >= 2 and HAS_SCIPY and scipy_stats is not None:
        t_res = scipy_stats.ttest_1samp(arr, popmean=0.0, nan_policy="omit")
        out["t_stat"] = float(t_res.statistic)
        out["t_pvalue"] = float(t_res.pvalue)
    nz = arr[arr != 0.0]
    if nz.size and HAS_SCIPY and scipy_stats is not None:
        npos = int(np.sum(nz > 0.0))
        try:
            sign_res = scipy_stats.binomtest(npos, n=int(nz.size), p=0.5, alternative="two-sided")
            out["sign_pvalue"] = float(sign_res.pvalue)
        except Exception:
            out["sign_pvalue"] = None
    return out



def aggregate_rows(rows: List[Dict[str, float | int | str]], modes: Sequence[str]) -> Dict[str, Dict[str, Dict[str, float | int | None]]]:
    out: Dict[str, Dict[str, Dict[str, float | int | None]]] = {}
    for mode in modes:
        mode_rows = [r for r in rows if r["mode"] == mode and r["status"] == "ok"]
        out[mode] = {}
        for metric in SUMMARY_METRICS:
            vals = [float(r[metric]) for r in mode_rows if metric in r and not np.isnan(float(r[metric]))]
            out[mode][metric] = one_sample_stats(vals)
    return out


# -----------------------------------------------------------------------------
# Visualisation batch
# -----------------------------------------------------------------------------


def boxplot_metric(rows: List[Dict[str, float | int | str]], modes: Sequence[str], metric: str, outpath: Path) -> None:
    data = []
    labels = []
    for mode in modes:
        vals = [float(r[metric]) for r in rows if r["mode"] == mode and r["status"] == "ok" and not np.isnan(float(r[metric]))]
        if vals:
            data.append(vals)
            labels.append(mode)
    if not data:
        return
    plt.figure(figsize=(8.4, 4.8))
    plt.boxplot(data, tick_labels=labels, showmeans=True)
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.ylabel(metric)
    plt.title(f"Distribution par mode : {metric}")
    plt.grid(True, axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()



def scatter_near_vs_far(rows: List[Dict[str, float | int | str]], mode: str, outpath: Path) -> None:
    mode_rows = [r for r in rows if r["mode"] == mode and r["status"] == "ok"]
    if not mode_rows:
        return
    x = np.array([float(r["delta_far_edge_mean"]) for r in mode_rows], dtype=float)
    y = np.array([float(r["delta_near_edge_mean"]) for r in mode_rows], dtype=float)
    seeds = [int(r["seed"]) for r in mode_rows]

    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        return
    x = x[finite]
    y = y[finite]
    seeds = [s for s, keep in zip(seeds, finite) if keep]

    lim = max(1e-6, float(np.max(np.abs(np.concatenate([x, y])))))
    lim *= 1.08

    plt.figure(figsize=(5.4, 5.0))
    plt.scatter(x, y, alpha=0.8)
    for sx, sy, sd in zip(x, y, seeds):
        plt.text(sx, sy, str(sd), fontsize=7, alpha=0.7)
    plt.plot([-lim, lim], [-lim, lim], linestyle="--", linewidth=1)
    plt.axhline(0.0, linestyle=":", linewidth=1)
    plt.axvline(0.0, linestyle=":", linewidth=1)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.xlabel("Δ far edge mean")
    plt.ylabel("Δ near edge mean")
    plt.title(f"Near vs far — {mode}")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()


# -----------------------------------------------------------------------------
# IO tabulaire
# -----------------------------------------------------------------------------


def write_rows_csv(rows: List[Dict[str, float | int | str]], outpath: Path) -> None:
    fieldnames = ["seed", "mode", "status"] + SUMMARY_METRICS
    with outpath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})



def write_summary_csv(summary: Dict[str, Dict[str, Dict[str, float | int | None]]], outpath: Path) -> None:
    fieldnames = [
        "mode",
        "metric",
        "n",
        "mean",
        "std",
        "median",
        "q25",
        "q75",
        "min",
        "max",
        "positive_frac",
        "negative_frac",
        "nonzero_frac",
        "t_stat",
        "t_pvalue",
        "sign_pvalue",
    ]
    with outpath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for mode, mode_metrics in summary.items():
            for metric, stats_dict in mode_metrics.items():
                row = {"mode": mode, "metric": metric}
                row.update(stats_dict)
                writer.writerow(row)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch multi-seeds: réponse de la courbure Ollivier-Ricci à une perturbation informationnelle locale."
    )
    parser.add_argument("--n", type=int, default=16, help="Nombre de qubits")
    parser.add_argument("--layers", type=int, default=15, help="Nombre de couches du circuit")
    parser.add_argument("--gates", type=int, default=15, help="Nombre de portes à deux qubits par couche")
    parser.add_argument("--lam", type=float, default=0.5, help="Probabilité de couplage non local")
    parser.add_argument("--theta", type=float, default=0.45, help="Angle moyen du partial-iSWAP")
    parser.add_argument("--seed0", type=int, default=42, help="Première seed")
    parser.add_argument("--n-seeds", type=int, default=20, help="Nombre de seeds")

    parser.add_argument("--source-region", type=str, default="5,6", help="Région source, ex: 5,6")
    parser.add_argument("--defect-radius", type=int, default=1, help="Rayon local du défaut sur l'anneau")
    parser.add_argument("--defect-strength", type=float, default=0.60, help="Intensité du défaut local")
    parser.add_argument(
        "--defect-mode",
        type=str,
        default="enhance",
        choices=["enhance", "suppress"],
        help="Renforcer ou supprimer les couplages près de la source",
    )

    parser.add_argument("--density", type=float, default=0.333, help="Densité cible du graphe")
    parser.add_argument("--near-radius", type=int, default=1, help="Rayon local pour les stats near/far")
    parser.add_argument("--alpha", type=float, default=0.5, help="Paramètre alpha d'Ollivier-Ricci")
    parser.add_argument(
        "--backend",
        type=str,
        default="auto",
        choices=["auto", "internal", "graphricci"],
        help="Backend pour calculer la courbure",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="OTD",
        choices=["OTD", "ATD", "Sinkhorn", "OTDSinkhornMix"],
        help="Méthode de transport pour backend graphricci",
    )
    parser.add_argument("--proc", type=int, default=1, help="Nombre de processus pour graphricci")
    parser.add_argument(
        "--verbose",
        type=str,
        default="ERROR",
        choices=["INFO", "TRACE", "DEBUG", "ERROR"],
        help="Niveau de verbosité de GraphRicciCurvature",
    )
    parser.add_argument(
        "--graph-mode",
        type=str,
        default="all",
        choices=["all", "free", "frozen_baseline", "frozen_perturbed", "frozen_intersection", "frozen_union"],
        help="Mode d'analyse du support du graphe",
    )
    parser.add_argument("--out", type=str, default="results_orc_local_batch_v1_3", help="Dossier de sortie")
    args = parser.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    source_region = parse_region(args.source_region)
    if not source_region:
        raise ValueError("La région source ne peut pas être vide.")
    if any(q < 0 or q >= args.n for q in source_region):
        raise ValueError(f"Région source invalide pour n={args.n}: {source_region}")
    if args.n_seeds <= 0:
        raise ValueError("--n-seeds doit être > 0")

    modes = requested_modes(args.graph_mode)

    print("=" * 72)
    print("BUp Ollivier-Ricci Local Response Batch v1.3")
    print("=" * 72)
    print(f"n={args.n} layers={args.layers} gates={args.gates} lam={args.lam} theta={args.theta}")
    print(
        f"source_region={source_region} defect_radius={args.defect_radius} "
        f"defect_strength={args.defect_strength} defect_mode={args.defect_mode}"
    )
    print(
        f"density={args.density} alpha={args.alpha} backend={args.backend} "
        f"graph_mode={args.graph_mode} seed0={args.seed0} n_seeds={args.n_seeds}"
    )
    print()

    rows: List[Dict[str, float | int | str]] = []
    batch_logs: List[Dict[str, object]] = []
    backend_used_final = args.backend

    for k in range(args.n_seeds):
        seed = args.seed0 + k
        print(f"[{k+1}/{args.n_seeds}] Seed {seed}...")
        try:
            psi0, psi1 = build_states_for_seed(args, seed, source_region)
            W0, _ = mutual_information_matrix(psi0, args.n)
            W1, _ = mutual_information_matrix(psi1, args.n)
            G0_free = build_connected_threshold_graph(W0, density=args.density)
            G1_free = build_connected_threshold_graph(W1, density=args.density)

            seed_report = {"seed": seed, "modes": {}}
            for mode in modes:
                mode_report, backend_used = compute_mode_analysis(mode, G0_free, G1_free, W0, W1, source_region, args)
                seed_report["modes"][mode] = mode_report
                rows.append(row_from_mode_report(seed, mode, mode_report))
                backend_used_final = backend_used
            batch_logs.append(seed_report)

            short_parts = []
            for mode in modes:
                mr = seed_report["modes"][mode]
                if mr.get("status") == "ok":
                    d = mr["comparison"]["delta_summary"]
                    short_parts.append(
                        f"{mode}: Δ={d['delta_edge_mean']:+.3f}, near-far={d['delta_near_minus_far']:+.3f}"
                    )
                else:
                    short_parts.append(f"{mode}: {mr.get('status')}")
            print("   " + " | ".join(short_parts))

        except Exception as exc:
            print(f"   ERREUR seed {seed}: {exc}")
            for mode in modes:
                rows.append({"seed": seed, "mode": mode, "status": f"error: {exc}"})
            batch_logs.append({"seed": seed, "error": str(exc)})

    summary = aggregate_rows(rows, modes)

    write_rows_csv(rows, outdir / "per_seed_results.csv")
    write_summary_csv(summary, outdir / "summary_by_mode.csv")

    boxplot_metric(rows, modes, "delta_edge_mean", outdir / "boxplot_delta_edge_mean.png")
    boxplot_metric(rows, modes, "delta_near_minus_far", outdir / "boxplot_delta_near_minus_far.png")
    boxplot_metric(rows, modes, "common_edge_delta_mean", outdir / "boxplot_common_edge_delta_mean.png")
    for mode in modes:
        scatter_near_vs_far(rows, mode, outdir / f"scatter_near_vs_far_{mode}.png")

    report = {
        "config": {
            "n": args.n,
            "layers": args.layers,
            "gates": args.gates,
            "lam": args.lam,
            "theta": args.theta,
            "seed0": args.seed0,
            "n_seeds": args.n_seeds,
            "source_region": source_region,
            "defect_radius": args.defect_radius,
            "defect_strength": args.defect_strength,
            "defect_mode": args.defect_mode,
            "density": args.density,
            "near_radius": args.near_radius,
            "alpha": args.alpha,
            "backend_requested": args.backend,
            "backend_used": backend_used_final,
            "graph_mode": args.graph_mode,
            "graphricci_available": HAS_GRAPH_RICCI,
            "scipy_available": HAS_SCIPY,
        },
        "modes": summary,
        "logs": batch_logs,
        "artifacts": {
            "per_seed_results": str(outdir / "per_seed_results.csv"),
            "summary_by_mode": str(outdir / "summary_by_mode.csv"),
            "boxplot_delta_edge_mean": str(outdir / "boxplot_delta_edge_mean.png"),
            "boxplot_delta_near_minus_far": str(outdir / "boxplot_delta_near_minus_far.png"),
            "boxplot_common_edge_delta_mean": str(outdir / "boxplot_common_edge_delta_mean.png"),
        },
    }
    save_json(report, outdir / "report_batch.json")

    print("\nRésumé agrégé:")
    for mode in modes:
        stats1 = summary[mode]["delta_edge_mean"]
        stats2 = summary[mode]["delta_near_minus_far"]
        print(
            f"- {mode}: "
            f"Δedge mean={stats1['mean']:+.4f} ± {stats1['std']:.4f} "
            f"(pos_frac={stats1['positive_frac']:.2f}, p_sign={stats1['sign_pvalue']}) | "
            f"Δ(near-far)={stats2['mean']:+.4f} ± {stats2['std']:.4f} "
            f"(pos_frac={stats2['positive_frac']:.2f}, p_sign={stats2['sign_pvalue']})"
        )

    print(f"\nSorties sauvegardées dans: {outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
BUp Ollivier-Ricci Local Response v1.2
-------------------------------------------------------------------------------
Objectif
--------
Tester si une perturbation informationnelle localisée modifie la courbure
Ollivier-Ricci du graphe d'intrication / corrélation, en séparant :

1) la réponse géométrique totale (graphe reconstruit librement),
2) la réponse de courbure à connectivité fixée (graphe gelé),
3) un contrôle par arêtes communes.

Nouveautés v1.2
---------------
- Mode "free" : reconstruction libre baseline vs perturbé.
- Modes "frozen_*" : même adjacency, poids/similarités propres à chaque état.
- Mode "intersection" : contrôle par arêtes communes uniquement.
- Option --graph-mode all pour calculer tous les modes d'un coup.
- Rapport JSON séparant clairement courbure pure vs rewiring.

Deux modes d'usage
------------------
1) Génération interne (par défaut)
   - Construit deux états purs à partir du même circuit aléatoire.
   - L'état "perturbé" reçoit un défaut local d'intrication autour d'une région
     source (renforcement ou suppression locale des couplages).

2) Chargement d'états externes
   - Fournir --psi0 et --psi1 (fichiers .npy de vecteurs d'état complexes).

Backends ORC
------------
- internal : implémentation auto-contenue de la courbure d'Ollivier-Ricci
  via Wasserstein-1 exact (programmation linéaire).
- graphricci : backend GraphRicciCurvature si installé.
- auto : utilise graphricci si disponible, sinon internal.

Sorties principales
-------------------
- report.json
- mi_baseline.npy / mi_perturbed.npy
- mode_<nom>/curvature_node_distance.png
- mode_<nom>/curvature_edge_hist.png
- mode_<nom>/delta_edge_vs_distance.png
- mode_<nom>/top_changed_edges.csv

Dépendances minimales
---------------------
pip install numpy scipy networkx matplotlib

Optionnel pour backend graphricci
---------------------------------
pip install GraphRicciCurvature pot

Remarque physique importante
---------------------------
Une simple rotation locale unitaire *après coup* ne modifie pas l'information
mutuelle bipartite. Ici, la perturbation locale est modélisée comme un défaut
local dans le *processus de génération* de l'état (renforcement/suppression des
couplages touchant la région source), ce qui change réellement les corrélations.
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

    HAS_SCIPY = True
except Exception:  # pragma: no cover
    HAS_SCIPY = False
    linprog = None

try:
    from GraphRicciCurvature.OllivierRicci import OllivierRicci

    HAS_GRAPH_RICCI = True
except Exception:  # pragma: no cover
    HAS_GRAPH_RICCI = False
    OllivierRicci = None


# -----------------------------------------------------------------------------
# Utilitaires généraux
# -----------------------------------------------------------------------------


def parse_region(text: str) -> List[int]:
    text = text.strip()
    if not text:
        return []
    return [int(x) for x in text.split(",") if x.strip()]



def complex_npy_load(path: str) -> np.ndarray:
    arr = np.load(path)
    arr = np.asarray(arr, dtype=np.complex128).reshape(-1)
    norm = np.linalg.norm(arr)
    if norm == 0:
        raise ValueError(f"Vecteur d'état nul dans {path}")
    return arr / norm



def save_json(obj: Dict, path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)



def ring_distance(a: int, b: int, n: int) -> int:
    d = abs(a - b)
    return min(d, n - d)



def distance_to_region_on_ring(node: int, region: Sequence[int], n: int) -> int:
    return min(ring_distance(node, r, n) for r in region)


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
# Construction du graphe de corrélation / intrication
# -----------------------------------------------------------------------------


def build_connected_threshold_graph(W: np.ndarray, density: float, eps_length: float = 1e-9) -> nx.Graph:
    n = W.shape[0]
    assert W.shape == (n, n)

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



def _safe_mean(x: List[float]) -> float:
    return float(np.mean(x)) if x else float("nan")



def _safe_stat(x: List[float], fn) -> float:
    return float(fn(x)) if x else float("nan")



def summarize_graph_curvature(G: nx.Graph, node_dist: Dict[int, float], near_radius: int) -> Dict[str, float]:
    edge_vals = []
    edge_near = []
    edge_far = []
    edge_unreachable = []
    for u, v, data in G.edges(data=True):
        rc = float(data.get("ricciCurvature", np.nan))
        if np.isnan(rc):
            continue
        edge_vals.append(rc)
        d = edge_source_distance(int(u), int(v), node_dist)
        if math.isinf(d):
            edge_unreachable.append(rc)
        elif d <= near_radius:
            edge_near.append(rc)
        else:
            edge_far.append(rc)

    node_vals = [float(G.nodes[n].get("ricciCurvature", np.nan)) for n in G.nodes()]
    node_vals = [x for x in node_vals if not np.isnan(x)]

    near_mean = _safe_mean(edge_near)
    far_mean = _safe_mean(edge_far)

    return {
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "connected": bool(nx.is_connected(G)) if G.number_of_nodes() > 0 else False,
        "edge_mean": _safe_stat(edge_vals, np.mean),
        "edge_median": _safe_stat(edge_vals, np.median),
        "edge_min": _safe_stat(edge_vals, np.min),
        "edge_max": _safe_stat(edge_vals, np.max),
        "edge_std": _safe_stat(edge_vals, np.std),
        "node_mean": _safe_stat(node_vals, np.mean),
        "node_median": _safe_stat(node_vals, np.median),
        "node_std": _safe_stat(node_vals, np.std),
        "near_edge_mean": near_mean,
        "far_edge_mean": far_mean,
        "near_edge_std": _safe_stat(edge_near, np.std),
        "far_edge_std": _safe_stat(edge_far, np.std),
        "unreachable_edge_mean": _safe_stat(edge_unreachable, np.mean),
        "n_edge_near": int(len(edge_near)),
        "n_edge_far": int(len(edge_far)),
        "n_edge_unreachable": int(len(edge_unreachable)),
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
) -> Tuple[Dict, List[Dict[str, float]], Dict[int, Dict[str, float]]]:
    node_dist = node_source_distances(G0, source_region)
    sum0 = summarize_graph_curvature(G0, node_dist, near_radius)
    sum1 = summarize_graph_curvature(G1, node_dist, near_radius)

    e0 = edge_dict(G0)
    e1 = edge_dict(G1)
    all_edges = sorted(set(e0) | set(e1))
    common_edges = sorted(set(e0) & set(e1))

    edge_rows: List[Dict[str, float]] = []
    for e in all_edges:
        r0 = e0.get(e)
        r1 = e1.get(e)
        rc0 = float(r0["ricciCurvature"]) if r0 else np.nan
        rc1 = float(r1["ricciCurvature"]) if r1 else np.nan
        d = edge_source_distance(e[0], e[1], node_dist)
        edge_rows.append(
            {
                "u": int(e[0]),
                "v": int(e[1]),
                "present_baseline": int(r0 is not None),
                "present_perturbed": int(r1 is not None),
                "edge_distance": float(d),
                "rc_baseline": rc0,
                "rc_perturbed": rc1,
                "delta_rc": float(rc1 - rc0) if (not np.isnan(rc0) and not np.isnan(rc1)) else np.nan,
                "sim_baseline": float(r0["similarity"]) if r0 else np.nan,
                "sim_perturbed": float(r1["similarity"]) if r1 else np.nan,
            }
        )

    all_nodes = sorted(set(G0.nodes()) | set(G1.nodes()))
    node_rows: Dict[int, Dict[str, float]] = {}
    for n in all_nodes:
        rc0 = float(G0.nodes[n].get("ricciCurvature", np.nan)) if n in G0 else np.nan
        rc1 = float(G1.nodes[n].get("ricciCurvature", np.nan)) if n in G1 else np.nan
        d = float(node_dist.get(int(n), math.inf))
        node_rows[int(n)] = {
            "distance": d,
            "rc_baseline": rc0,
            "rc_perturbed": rc1,
            "delta_rc": float(rc1 - rc0) if (not np.isnan(rc0) and not np.isnan(rc1)) else np.nan,
        }

    common_deltas = [r["delta_rc"] for r in edge_rows if not np.isnan(r["delta_rc"])]
    near_common = [
        r["delta_rc"]
        for r in edge_rows
        if (not np.isnan(r["delta_rc"])) and (not math.isinf(r["edge_distance"])) and r["edge_distance"] <= near_radius
    ]
    far_common = [
        r["delta_rc"]
        for r in edge_rows
        if (not np.isnan(r["delta_rc"])) and (not math.isinf(r["edge_distance"])) and r["edge_distance"] > near_radius
    ]

    delta_edge_mean = float(sum1["edge_mean"] - sum0["edge_mean"]) if (not np.isnan(sum0["edge_mean"]) and not np.isnan(sum1["edge_mean"])) else float("nan")
    delta_node_mean = float(sum1["node_mean"] - sum0["node_mean"]) if (not np.isnan(sum0["node_mean"]) and not np.isnan(sum1["node_mean"])) else float("nan")
    delta_near_edge_mean = float(sum1["near_edge_mean"] - sum0["near_edge_mean"]) if (not np.isnan(sum0["near_edge_mean"]) and not np.isnan(sum1["near_edge_mean"])) else float("nan")
    delta_far_edge_mean = float(sum1["far_edge_mean"] - sum0["far_edge_mean"]) if (not np.isnan(sum0["far_edge_mean"]) and not np.isnan(sum1["far_edge_mean"])) else float("nan")
    delta_near_minus_far = float(sum1["near_minus_far"] - sum0["near_minus_far"]) if (not np.isnan(sum0["near_minus_far"]) and not np.isnan(sum1["near_minus_far"])) else float("nan")

    comparison = {
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
    return comparison, edge_rows, node_rows


# -----------------------------------------------------------------------------
# Contrôles de graphe gelé
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
# Visualisation
# -----------------------------------------------------------------------------


def plot_node_curvature_vs_distance(node_rows: Dict[int, Dict[str, float]], outpath: Path) -> None:
    nodes = sorted(node_rows)
    x = np.array([node_rows[n]["distance"] for n in nodes], dtype=float)
    y0 = np.array([node_rows[n]["rc_baseline"] for n in nodes], dtype=float)
    y1 = np.array([node_rows[n]["rc_perturbed"] for n in nodes], dtype=float)

    finite_x = np.where(np.isfinite(x), x, np.nan)

    plt.figure(figsize=(7, 4.5))
    plt.scatter(finite_x - 0.04, y0, label="baseline", alpha=0.8)
    plt.scatter(finite_x + 0.04, y1, label="perturbed", alpha=0.8)
    plt.xlabel("distance en sauts depuis la région source")
    plt.ylabel("courbure nodale Ollivier-Ricci")
    plt.title("Courbure nodale vs distance à la source")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()



def plot_edge_curvature_hist(edge_rows: List[Dict[str, float]], outpath: Path) -> None:
    rc0 = np.array([r["rc_baseline"] for r in edge_rows if not np.isnan(r["rc_baseline"])] , dtype=float)
    rc1 = np.array([r["rc_perturbed"] for r in edge_rows if not np.isnan(r["rc_perturbed"])] , dtype=float)

    plt.figure(figsize=(7, 4.5))
    if rc0.size:
        plt.hist(rc0, bins=24, alpha=0.55, label="baseline")
    if rc1.size:
        plt.hist(rc1, bins=24, alpha=0.55, label="perturbed")
    plt.xlabel("courbure d'arête Ollivier-Ricci")
    plt.ylabel("nombre d'arêtes")
    plt.title("Distribution des courbures d'arête")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()



def plot_delta_edge_vs_distance(edge_rows: List[Dict[str, float]], outpath: Path) -> None:
    rows = [r for r in edge_rows if not np.isnan(r["delta_rc"]) and not math.isinf(r["edge_distance"])]
    if not rows:
        return
    x = np.array([r["edge_distance"] for r in rows], dtype=float)
    y = np.array([r["delta_rc"] for r in rows], dtype=float)

    by_d: Dict[int, List[float]] = {}
    for r in rows:
        by_d.setdefault(int(r["edge_distance"]), []).append(float(r["delta_rc"]))

    xs = sorted(by_d)
    ys = [float(np.mean(by_d[d])) for d in xs]

    plt.figure(figsize=(7, 4.5))
    plt.scatter(x, y, alpha=0.55, label="arêtes communes")
    plt.plot(xs, ys, marker="o", linewidth=2, label="moyenne par distance")
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("distance d'arête à la région source")
    plt.ylabel("Δ courbure = RC_perturbed - RC_baseline")
    plt.title("Réponse locale de la courbure")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=160)
    plt.close()



def save_top_changed_edges(edge_rows: List[Dict[str, float]], outpath: Path, topk: int = 30) -> None:
    rows = [r for r in edge_rows if not np.isnan(r["delta_rc"])]
    rows.sort(key=lambda r: abs(r["delta_rc"]), reverse=True)
    rows = rows[:topk]
    with outpath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "u",
                "v",
                "edge_distance",
                "rc_baseline",
                "rc_perturbed",
                "delta_rc",
                "sim_baseline",
                "sim_perturbed",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in writer.fieldnames})


# -----------------------------------------------------------------------------
# Pipeline principal
# -----------------------------------------------------------------------------


def build_states(args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray]:
    if args.psi0 and args.psi1:
        psi0 = complex_npy_load(args.psi0)
        psi1 = complex_npy_load(args.psi1)
        expected_dim = 2**args.n
        if psi0.size != expected_dim or psi1.size != expected_dim:
            raise ValueError(
                f"Dimension incohérente: attendu 2**n = {expected_dim}, "
                f"reçu {psi0.size} et {psi1.size}"
            )
        return psi0, psi1

    if bool(args.psi0) ^ bool(args.psi1):
        raise ValueError("Fournir soit les deux fichiers --psi0 et --psi1, soit aucun.")

    source_region = parse_region(args.source_region)
    if not source_region:
        raise ValueError("source_region vide")

    schedule = make_schedule(
        n=args.n,
        layers=args.layers,
        gates_per_layer=args.gates,
        lam=args.lam,
        theta_base=args.theta,
        seed=args.seed,
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
        mode_report["baseline_connected"] = bool(nx.is_connected(G0_raw))
        mode_report["perturbed_connected"] = bool(nx.is_connected(G1_raw))
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

    comparison, edge_rows, node_rows = compare_graphs(
        G0,
        G1,
        source_region=source_region,
        near_radius=args.near_radius,
    )

    mode_dir = Path(args.out) / f"mode_{mode}"
    mode_dir.mkdir(parents=True, exist_ok=True)
    plot_node_curvature_vs_distance(node_rows, mode_dir / "curvature_node_distance.png")
    plot_edge_curvature_hist(edge_rows, mode_dir / "curvature_edge_hist.png")
    plot_delta_edge_vs_distance(edge_rows, mode_dir / "delta_edge_vs_distance.png")
    save_top_changed_edges(edge_rows, mode_dir / "top_changed_edges.csv", topk=30)

    mode_report["comparison"] = comparison
    mode_report["artifacts"] = {
        "dir": str(mode_dir),
        "curvature_node_distance": str(mode_dir / "curvature_node_distance.png"),
        "curvature_edge_hist": str(mode_dir / "curvature_edge_hist.png"),
        "delta_edge_vs_distance": str(mode_dir / "delta_edge_vs_distance.png"),
        "top_changed_edges": str(mode_dir / "top_changed_edges.csv"),
    }

    return mode_report, backend_used



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tester la réponse de la courbure Ollivier-Ricci à une perturbation informationnelle locale."
    )
    parser.add_argument("--n", type=int, default=12, help="Nombre de qubits")
    parser.add_argument("--layers", type=int, default=12, help="Nombre de couches du circuit")
    parser.add_argument("--gates", type=int, default=6, help="Nombre de portes à deux qubits par couche")
    parser.add_argument("--lam", type=float, default=0.35, help="Probabilité de couplage non local")
    parser.add_argument("--theta", type=float, default=0.45, help="Angle moyen du partial-iSWAP")
    parser.add_argument("--seed", type=int, default=42, help="Seed RNG")

    parser.add_argument("--source-region", type=str, default="5,6", help="Région source, ex: 5,6")
    parser.add_argument("--defect-radius", type=int, default=1, help="Rayon local du défaut sur l'anneau")
    parser.add_argument(
        "--defect-strength",
        type=float,
        default=0.60,
        help="Intensité du défaut local (facteur relatif)",
    )
    parser.add_argument(
        "--defect-mode",
        type=str,
        default="enhance",
        choices=["enhance", "suppress"],
        help="Renforcer ou supprimer les couplages près de la source",
    )

    parser.add_argument("--density", type=float, default=0.33, help="Densité cible du graphe")
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

    parser.add_argument("--psi0", type=str, default="", help=".npy état baseline externe")
    parser.add_argument("--psi1", type=str, default="", help=".npy état perturbé externe")
    parser.add_argument("--out", type=str, default="results_ollivier_ricci_local_response_v1_2", help="Dossier de sortie")
    args = parser.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    source_region = parse_region(args.source_region)
    if not source_region:
        raise ValueError("La région source ne peut pas être vide.")
    if any(q < 0 or q >= args.n for q in source_region):
        raise ValueError(f"Région source invalide pour n={args.n}: {source_region}")

    print("=" * 72)
    print("BUp Ollivier-Ricci Local Response v1.2")
    print("=" * 72)
    print(f"n={args.n} layers={args.layers} gates={args.gates} lam={args.lam} theta={args.theta}")
    print(
        f"source_region={source_region} defect_radius={args.defect_radius} "
        f"defect_strength={args.defect_strength} defect_mode={args.defect_mode}"
    )
    print(
        f"density={args.density} alpha={args.alpha} backend={args.backend} "
        f"graph_mode={args.graph_mode} seed={args.seed}"
    )
    print()

    psi0, psi1 = build_states(args)
    np.save(outdir / "psi_baseline.npy", psi0)
    np.save(outdir / "psi_perturbed.npy", psi1)

    print("[1/5] Calcul des matrices d'information mutuelle...")
    W0, s10 = mutual_information_matrix(psi0, args.n)
    W1, s11 = mutual_information_matrix(psi1, args.n)
    np.save(outdir / "mi_baseline.npy", W0)
    np.save(outdir / "mi_perturbed.npy", W1)

    print("[2/5] Construction des graphes libres connectés à densité fixée...")
    G0_free = build_connected_threshold_graph(W0, density=args.density)
    G1_free = build_connected_threshold_graph(W1, density=args.density)

    print("[3/5] Calcul de la courbure Ollivier-Ricci selon les modes demandés...")
    mode_reports: Dict[str, Dict] = {}
    backend_used_final = args.backend
    for mode in requested_modes(args.graph_mode):
        print(f"  - mode {mode}...")
        mode_report, backend_used = compute_mode_analysis(mode, G0_free, G1_free, W0, W1, source_region, args)
        mode_reports[mode] = mode_report
        backend_used_final = backend_used

    print("[4/5] Agrégation du rapport...")
    free_edges0 = edge_set(G0_free)
    free_edges1 = edge_set(G1_free)
    free_graph_overlap = {
        "edge_jaccard": float(len(free_edges0 & free_edges1) / max(1, len(free_edges0 | free_edges1))),
        "edges_only_baseline": int(len(free_edges0 - free_edges1)),
        "edges_only_perturbed": int(len(free_edges1 - free_edges0)),
        "n_common_edges": int(len(free_edges0 & free_edges1)),
        "n_union_edges": int(len(free_edges0 | free_edges1)),
    }

    print("[5/5] Sauvegarde des sorties...")
    report = {
        "config": {
            "n": args.n,
            "layers": args.layers,
            "gates": args.gates,
            "lam": args.lam,
            "theta": args.theta,
            "seed": args.seed,
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
        "single_qubit_entropy": {
            "baseline": [float(x) for x in s10],
            "perturbed": [float(x) for x in s11],
        },
        "free_graph_overlap": free_graph_overlap,
        "modes": mode_reports,
    }
    save_json(report, outdir / "report.json")

    print("\nRésumé synthétique:")
    for mode in requested_modes(args.graph_mode):
        mr = mode_reports.get(mode, {})
        status = mr.get("status", "?")
        print(f"- {mode}: status={status}")
        if status == "ok":
            delta = mr["comparison"]["delta_summary"]
            print(
                "  "
                + json.dumps(
                    {
                        "delta_edge_mean": delta["delta_edge_mean"],
                        "delta_near_edge_mean": delta["delta_near_edge_mean"],
                        "delta_far_edge_mean": delta["delta_far_edge_mean"],
                        "delta_near_minus_far": delta["delta_near_minus_far"],
                        "common_edge_delta_mean": delta["common_edge_delta_mean"],
                        "common_edge_delta_mean_near": delta["common_edge_delta_mean_near"],
                        "common_edge_delta_mean_far": delta["common_edge_delta_mean_far"],
                    },
                    ensure_ascii=False,
                )
            )

    print(f"\nSorties sauvegardées dans: {outdir}")


if __name__ == "__main__":
    main()

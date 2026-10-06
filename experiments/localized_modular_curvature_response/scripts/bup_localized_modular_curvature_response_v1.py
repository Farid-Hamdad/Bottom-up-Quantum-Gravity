#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP — Localized Modular Response and Emergent Ollivier–Ricci Curvature

Implementation target:
  experiments/localized_modular_curvature_response/PREREGISTRATION_FROZEN.md

Safety:
- Default mode is baseline-only.
- The 36-condition fresh N=12 confirmatory grid requires BOTH:
    --mode confirmatory
    --confirm-unblind I_UNDERSTAND_THIS_EXECUTES_THE_FRESH_N12_GRID
- No network access, no Git operations, no adaptive parameter selection.
- The new localized curvature observable is never evaluated on the prior N=10 raw dataset.

Frozen provenance:
- preregistration commit:
  00aceb88318f42fbc3c2a59f99d483a74b8bfabb
- previous frozen confirmatory-results commit:
  2aab1eb557010e4d6edc4305ffd992f40d41159a
- previous completed-experiment synthesis commit:
  3dccfbc7c342c1b55428a1e3f4cb4139f8a51e2c
- Paper 28 source blob:
  5fc9e845115a68135199e31f65324290ae0f210e
- Ollivier–Ricci source blob:
  3fd6eceb9c107075c80abf9eefd8fa73db1b69cd
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import networkx as nx
import numpy as np

try:
    from scipy.optimize import linprog
    from scipy.stats import rankdata
except Exception as exc:  # pragma: no cover
    raise RuntimeError("scipy is required for this preregistered benchmark") from exc


# =============================================================================
# Frozen scientific constants
# =============================================================================

N = 12
J = 1.0
H_FIELD = 1.0
P = 6
DT = 0.35
NOMINAL_T = P * DT

A = (0, 1, 2, 3)
A_PAIRS = tuple((i, j) for i in A for j in A if i < j)

EPSILONS = (0.075, 0.15, 0.30)
SITES = tuple(range(N))

GRAPH_DENSITY = 0.333
OR_ALPHA = 0.5
EPS_LENGTH = 1e-9
STATE_NORM_TOL = 1e-12
HERMITICITY_TOL = 1e-12
NEG_EIG_TOL = 1e-12

PERMUTATIONS = 100_000
PERMUTATION_SEED = 20261006
SIGNIFICANCE_ALPHA = 0.05

PREREG_COMMIT = "00aceb88318f42fbc3c2a59f99d483a74b8bfabb"
PREREG_SHA256 = "a96703fee17d0bdd53bbe910aca687672c8144a3cbaf836e87be10547031ae8d"
PAPER28_BLOB = "5fc9e845115a68135199e31f65324290ae0f210e"
OR_BLOB = "3fd6eceb9c107075c80abf9eefd8fa73db1b69cd"
PREVIOUS_CONFIRMATORY_COMMIT = "2aab1eb557010e4d6edc4305ffd992f40d41159a"
PREVIOUS_SYNTHESIS_COMMIT = "3dccfbc7c342c1b55428a1e3f4cb4139f8a51e2c"

CONFIRM_TOKEN = "I_UNDERSTAND_THIS_EXECUTES_THE_FRESH_N12_GRID"

I2 = np.eye(2, dtype=np.complex128)
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
PAULI = {"X": X, "Y": Y, "Z": Z}


# =============================================================================
# IO helpers
# =============================================================================

def write_json(obj: object, path: Path) -> None:
    path.write_text(
        json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_csv(rows: Sequence[Dict[str, object]], path: Path) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# =============================================================================
# Quantum dynamics — audited Paper 28 conventions
# =============================================================================

def plus_state(n: int) -> np.ndarray:
    state = np.array([1.0 + 0.0j])
    q = np.array([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
    for _ in range(n):
        state = np.kron(state, q)
    return state


def one_qubit_gate(
    state: np.ndarray,
    gate: np.ndarray,
    q: int,
    n: int,
) -> np.ndarray:
    tensor = np.moveaxis(state.reshape([2] * n), q, 0)
    out = np.tensordot(gate, tensor, axes=([1], [0]))
    return np.moveaxis(out, 0, q).reshape(-1)


def ry_gate(epsilon: float) -> np.ndarray:
    return np.cos(epsilon / 2.0) * I2 - 1j * np.sin(epsilon / 2.0) * Y


def initial_state(q: int | None = None, epsilon: float = 0.0) -> np.ndarray:
    state = plus_state(N)
    if q is not None and epsilon != 0.0:
        state = one_qubit_gate(state, ry_gate(epsilon), q, N)
    return state


def xhalf(state: np.ndarray, dt: float, n: int) -> np.ndarray:
    gate = np.cos(H_FIELD * dt / 2.0) * I2 + 1j * np.sin(H_FIELD * dt / 2.0) * X
    out = state
    for q in range(n):
        out = one_qubit_gate(out, gate, q, n)
    return out


def zz(state: np.ndarray, dt: float, n: int) -> np.ndarray:
    idx = np.arange(2**n, dtype=np.uint64)
    energy = np.zeros(2**n, dtype=float)
    for q in range(n - 1):
        bq = (idx >> np.uint64(n - 1 - q)) & np.uint64(1)
        br = (idx >> np.uint64(n - 2 - q)) & np.uint64(1)
        zq = 1.0 - 2.0 * bq.astype(float)
        zr = 1.0 - 2.0 * br.astype(float)
        energy += zq * zr
    return state * np.exp(1j * J * dt * energy)


def evolve_from_initial(initial: np.ndarray, dt: float) -> np.ndarray:
    state = np.asarray(initial, dtype=np.complex128).copy()
    for _ in range(P):
        state = xhalf(state, dt, N)
        state = zz(state, dt, N)
        state = xhalf(state, dt, N)
    norm = float(np.real(np.vdot(state, state)))
    if abs(norm - 1.0) > STATE_NORM_TOL:
        raise RuntimeError(f"STATE_NORM_DRIFT norm={norm:.17g}")
    return state


def pure_reduced_density(
    state: np.ndarray,
    keep: Sequence[int],
    n: int,
) -> np.ndarray:
    keep = tuple(int(q) for q in keep)
    rest = tuple(q for q in range(n) if q not in keep)
    perm = keep + rest
    tensor = state.reshape([2] * n)
    mat = np.transpose(tensor, perm).reshape(2 ** len(keep), -1)
    return mat @ mat.conj().T


def mixed_reduced_density_from_branches(
    psi_plus: np.ndarray,
    psi_minus: np.ndarray,
    keep: Sequence[int],
) -> np.ndarray:
    rp = pure_reduced_density(psi_plus, keep, N)
    rm = pure_reduced_density(psi_minus, keep, N)
    return 0.5 * (rp + rm)


@dataclass(frozen=True)
class AnalysisBranches:
    psi_plus: np.ndarray
    psi_minus: np.ndarray


def build_analysis_branches(q: int | None, epsilon: float) -> AnalysisBranches:
    init = initial_state(q=q, epsilon=epsilon)
    return AnalysisBranches(
        psi_plus=evolve_from_initial(init, +DT),
        psi_minus=evolve_from_initial(init, -DT),
    )


# =============================================================================
# Density-matrix diagnostics and mutual information in nats
# =============================================================================

def hermitian_part(rho: np.ndarray) -> np.ndarray:
    return 0.5 * (rho + rho.conj().T)


def density_diagnostics(rho: np.ndarray) -> Dict[str, float]:
    rh = hermitian_part(rho)
    vals = np.real(np.linalg.eigvalsh(rh))
    herm_resid = float(np.linalg.norm(rho - rho.conj().T))
    trace = float(np.real(np.trace(rho)))
    lam_min = float(vals.min())
    lam_max = float(vals.max())
    positive = vals[vals > 0.0]
    cond = float(lam_max / positive.min()) if positive.size else float("inf")
    return {
        "trace": trace,
        "hermiticity_residual": herm_resid,
        "lambda_min": lam_min,
        "lambda_max": lam_max,
        "condition_number_positive_spectrum": cond,
    }


def entropy_nats(rho: np.ndarray) -> float:
    rh = hermitian_part(rho)
    vals = np.real(np.linalg.eigvalsh(rh))
    if float(vals.min()) < -NEG_EIG_TOL:
        raise RuntimeError(f"DENSITY_NEGATIVE_EIGENVALUE {float(vals.min()):.17g}")
    vals = np.clip(vals, 0.0, None)
    vals = vals[vals > 1e-15]
    return float(-np.sum(vals * np.log(vals)))


def mutual_information_matrix_mixed(
    branches: AnalysisBranches,
) -> Tuple[np.ndarray, Dict[Tuple[int, int], float], List[Dict[str, object]]]:
    single_entropy: Dict[int, float] = {}
    density_rows: List[Dict[str, object]] = []

    for i in range(N):
        rho_i = mixed_reduced_density_from_branches(
            branches.psi_plus, branches.psi_minus, (i,)
        )
        single_entropy[i] = entropy_nats(rho_i)
        density_rows.append({
            "subset_type": "single",
            "subset": str(i),
            **density_diagnostics(rho_i),
        })

    W = np.zeros((N, N), dtype=float)
    pairs: Dict[Tuple[int, int], float] = {}
    for i in range(N):
        for j in range(i + 1, N):
            rho_ij = mixed_reduced_density_from_branches(
                branches.psi_plus, branches.psi_minus, (i, j)
            )
            value = single_entropy[i] + single_entropy[j] - entropy_nats(rho_ij)
            density_rows.append({
                "subset_type": "pair",
                "subset": f"{i},{j}",
                **density_diagnostics(rho_ij),
            })
            # Numerical negative MI is not physically meaningful; values below
            # tiny roundoff are set to zero, matching historical graph handling.
            if value < -1e-12:
                raise RuntimeError(f"NEGATIVE_MUTUAL_INFORMATION pair={(i, j)} value={value}")
            value = float(max(0.0, value))
            W[i, j] = W[j, i] = value
            pairs[(i, j)] = value
    return W, pairs, density_rows


# =============================================================================
# Modular observable — audited Paper 28 conventions
# =============================================================================

def subsystem_pauli_ops(a_size: int) -> Dict[Tuple[int, str], np.ndarray]:
    out: Dict[Tuple[int, str], np.ndarray] = {}
    for q in range(a_size):
        for label, sigma in PAULI.items():
            op = np.array([[1.0 + 0.0j]])
            for k in range(a_size):
                op = np.kron(op, sigma if k == q else I2)
            out[(q, label)] = op
    return out


def comm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def normalized_modular_hamiltonian(
    rho_a: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, float]]:
    rh = hermitian_part(rho_a)
    vals, vecs = np.linalg.eigh(rh)
    vals = np.real(vals)

    lam_min = float(vals.min())
    lam_max = float(vals.max())

    # Frozen Paper 28 positivity gate: strict positivity, no floor/pseudolog.
    if lam_min <= 0.0:
        raise RuntimeError(f"POSITIVITY_GATE_FAIL lambda_min={lam_min:.17g}")

    kappas = -np.log(vals)
    kappa_mean = float(np.mean(kappas))
    kappa_std = float(np.std(kappas, ddof=0))
    if not np.isfinite(kappa_std) or kappa_std <= 0.0:
        raise RuntimeError("MODULAR_STD_FAIL")

    normalized_eigs = (kappas - kappa_mean) / kappa_std
    k_tilde = vecs @ np.diag(normalized_eigs) @ vecs.conj().T

    diag = {
        "lambda_min": lam_min,
        "lambda_max": lam_max,
        "condition_number": float(lam_max / lam_min),
        "kappa_mean": kappa_mean,
        "kappa_std": kappa_std,
    }
    return k_tilde, diag


def modular_pair_coefficients(
    rho_a: np.ndarray,
    k_tilde: np.ndarray,
) -> Dict[Tuple[int, int], float]:
    a_size = len(A)
    ops = subsystem_pauli_ops(a_size)
    coeffs: Dict[Tuple[int, int], float] = {}

    # A is the prefix (0,1,2,3), so local labels equal global labels here.
    for i, j in A_PAIRS:
        total = 0.0
        for pa in "XYZ":
            inner = comm(k_tilde, ops[(i, pa)])
            for pb in "XYZ":
                qop = comm(inner, ops[(j, pb)])
                total += float(np.real(np.trace(rho_a @ qop.conj().T @ qop)))
        coeff = total / 9.0
        if coeff < -1e-12:
            raise RuntimeError(f"NEGATIVE_MODULAR_COEFFICIENT pair={(i, j)} coeff={coeff}")
        coeffs[(i, j)] = float(max(0.0, coeff))
    return coeffs


def modular_observable(
    branches: AnalysisBranches,
) -> Tuple[float, Dict[Tuple[int, int], float], Dict[str, float], np.ndarray]:
    rho_a = mixed_reduced_density_from_branches(
        branches.psi_plus, branches.psi_minus, A
    )
    dd = density_diagnostics(rho_a)
    if abs(dd["trace"] - 1.0) > 1e-12:
        raise RuntimeError(f"RHO_A_TRACE_FAIL trace={dd['trace']:.17g}")
    if dd["hermiticity_residual"] > HERMITICITY_TOL:
        raise RuntimeError(
            f"RHO_A_HERMITICITY_FAIL residual={dd['hermiticity_residual']:.17g}"
        )

    k_tilde, mod_diag = normalized_modular_hamiltonian(rho_a)
    coeffs = modular_pair_coefficients(rho_a, k_tilde)
    v = {pair: float(np.sqrt(value)) for pair, value in coeffs.items()}
    m_k = float(np.mean(list(v.values())))

    diagnostics = dict(dd)
    diagnostics.update(mod_diag)
    return m_k, v, diagnostics, rho_a


# =============================================================================
# Graph construction — audited historical OR conventions
# =============================================================================

def build_graph_from_edges_and_weights(
    edges: Iterable[Tuple[int, int]],
    W: np.ndarray,
    eps_length: float = EPS_LENGTH,
) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(W.shape[0]))

    selected_pairs: List[Tuple[int, int, float]] = []
    for i, j in sorted({tuple(sorted((int(a), int(b)))) for a, b in edges}):
        sim = float(max(W[i, j], 0.0))
        selected_pairs.append((i, j, sim))

    smax = max((s for _, _, s in selected_pairs), default=1.0)
    smax = max(smax, eps_length)

    for i, j, sim in selected_pairs:
        sim_norm = sim / smax
        length = 1.0 / max(sim_norm, eps_length)
        graph.add_edge(
            i,
            j,
            similarity=sim,
            similarity_norm=sim_norm,
            length=length,
        )
    return graph


def build_connected_threshold_graph(
    W: np.ndarray,
    density: float = GRAPH_DENSITY,
) -> nx.Graph:
    n = W.shape[0]
    pairs: List[Tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, float(max(W[i, j], 0.0))))

    target_edges = int(round(density * n * (n - 1) / 2.0))
    target_edges = min(max(target_edges, n - 1), n * (n - 1) // 2)

    complete = nx.Graph()
    complete.add_nodes_from(range(n))
    for i, j, sim in pairs:
        complete.add_edge(i, j, similarity=sim)

    mst = nx.maximum_spanning_tree(complete, weight="similarity")
    selected = {tuple(sorted((int(u), int(v)))) for u, v in mst.edges()}

    remaining = sorted(
        pairs,
        key=lambda x: (x[2], -abs(x[0] - x[1])),
        reverse=True,
    )
    for i, j, _ in remaining:
        if len(selected) >= target_edges:
            break
        edge = tuple(sorted((i, j)))
        if edge not in selected:
            selected.add(edge)

    return build_graph_from_edges_and_weights(selected, W)


def edge_set(graph: nx.Graph) -> Tuple[Tuple[int, int], ...]:
    return tuple(
        sorted(tuple(sorted((int(u), int(v)))) for u, v in graph.edges())
    )


def local_measure(
    graph: nx.Graph,
    node: int,
    alpha: float,
) -> Tuple[List[int], np.ndarray]:
    neighbors = list(graph.neighbors(node))
    support = [node] + neighbors
    mass = np.zeros(len(support), dtype=float)
    mass[0] = alpha

    if neighbors:
        weights = []
        for nbr in neighbors:
            data = graph[node][nbr]
            sim = float(data.get("similarity_norm", 0.0))
            if sim <= 0.0:
                length = float(data.get("length", 1.0))
                sim = 1.0 / max(length, 1e-12)
            weights.append(sim)
        weights_arr = np.asarray(weights, dtype=float)
        if weights_arr.sum() <= 0.0:
            weights_arr = np.ones_like(weights_arr)
        weights_arr = weights_arr / weights_arr.sum()
        mass[1:] = (1.0 - alpha) * weights_arr
    else:
        mass[0] = 1.0

    return support, mass


def wasserstein_distance_lp(
    mu: np.ndarray,
    nu: np.ndarray,
    cost: np.ndarray,
) -> float:
    m, n = cost.shape
    c = cost.reshape(-1)

    A_eq: List[np.ndarray] = []
    b_eq: List[float] = []

    for i in range(m):
        row = np.zeros(m * n, dtype=float)
        row[i * n : (i + 1) * n] = 1.0
        A_eq.append(row)
        b_eq.append(float(mu[i]))

    for j in range(n):
        col = np.zeros(m * n, dtype=float)
        col[j::n] = 1.0
        A_eq.append(col)
        b_eq.append(float(nu[j]))

    res = linprog(
        c=c,
        A_eq=np.asarray(A_eq),
        b_eq=np.asarray(b_eq),
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        raise RuntimeError(f"OPTIMAL_TRANSPORT_FAIL {res.message}")
    return float(res.fun)


def compute_ollivier_ricci_internal(
    graph: nx.Graph,
    alpha: float = OR_ALPHA,
) -> nx.Graph:
    graph = graph.copy()
    if not nx.is_connected(graph):
        raise RuntimeError("OR_GRAPH_DISCONNECTED")

    apsp = dict(nx.all_pairs_dijkstra_path_length(graph, weight="length"))

    for u, v in graph.edges():
        support_u, mu = local_measure(graph, int(u), alpha)
        support_v, nu = local_measure(graph, int(v), alpha)

        cost = np.zeros((len(support_u), len(support_v)), dtype=float)
        for i, a in enumerate(support_u):
            for j, b in enumerate(support_v):
                cost[i, j] = float(apsp[a][b])

        w1 = wasserstein_distance_lp(mu, nu, cost)
        d_uv = float(graph[u][v]["length"])
        kappa = 1.0 - w1 / max(d_uv, 1e-12)
        graph[u][v]["ricciCurvature"] = float(kappa)

    for node in graph.nodes():
        vals = [
            float(graph[node][nbr]["ricciCurvature"])
            for nbr in graph.neighbors(node)
        ]
        graph.nodes[node]["ricciCurvature"] = (
            float(np.mean(vals)) if vals else float("nan")
        )
    return graph


def curvature_edge_mean(graph: nx.Graph) -> float:
    vals = [
        float(data["ricciCurvature"])
        for _, _, data in graph.edges(data=True)
    ]
    return float(np.mean(vals))


def graph_diagnostics(graph: nx.Graph) -> Dict[str, object]:
    sims = [
        float(data["similarity"])
        for _, _, data in graph.edges(data=True)
    ]
    sims_norm = [
        float(data["similarity_norm"])
        for _, _, data in graph.edges(data=True)
    ]
    lengths = [
        float(data["length"])
        for _, _, data in graph.edges(data=True)
    ]
    return {
        "n_nodes": int(graph.number_of_nodes()),
        "n_edges": int(graph.number_of_edges()),
        "connected": bool(nx.is_connected(graph)),
        "edges": [list(e) for e in edge_set(graph)],
        "similarity_min": float(min(sims)) if sims else None,
        "similarity_max": float(max(sims)) if sims else None,
        "similarity_norm_min": float(min(sims_norm)) if sims_norm else None,
        "similarity_norm_max": float(max(sims_norm)) if sims_norm else None,
        "length_min": float(min(lengths)) if lengths else None,
        "length_max": float(max(lengths)) if lengths else None,
    }


# =============================================================================
# Condition evaluation
# =============================================================================

@dataclass
class Baseline:
    branches: AnalysisBranches
    W: np.ndarray
    W_pairs: Dict[Tuple[int, int], float]
    support: Tuple[Tuple[int, int], ...]
    local_support: Tuple[Tuple[int, int], ...]
    m_w: float
    m_k: float
    v_pairs: Dict[Tuple[int, int], float]
    modular_diag: Dict[str, float]
    kappa_global: float
    kappa_a: float
    curvature_graph: nx.Graph
    density_rows: List[Dict[str, object]]


def m_w_from_pairs(pairs: Dict[Tuple[int, int], float]) -> float:
    return float(np.mean([pairs[p] for p in A_PAIRS]))


def localized_edge_set(
    support: Sequence[Tuple[int, int]],
) -> Tuple[Tuple[int, int], ...]:
    a_nodes = set(A)
    return tuple(
        sorted(
            e
            for e in (tuple(sorted((int(u), int(v)))) for u, v in support)
            if e[0] in a_nodes or e[1] in a_nodes
        )
    )


def curvature_mean_on_edges(
    graph: nx.Graph,
    edges: Sequence[Tuple[int, int]],
) -> float:
    values = [
        float(graph[u][v]["ricciCurvature"])
        for u, v in edges
    ]
    if not values:
        raise RuntimeError("EMPTY_CURVATURE_EDGE_SET")
    return float(np.mean(values))


def distance_to_A(q: int) -> int:
    if q in A:
        return 0
    return min(abs(int(q) - int(a)) for a in A)


def evaluate_baseline() -> Baseline:
    branches = build_analysis_branches(q=None, epsilon=0.0)
    W, W_pairs, density_rows = mutual_information_matrix_mixed(branches)
    m_w = m_w_from_pairs(W_pairs)

    m_k, v_pairs, modular_diag, rho_a = modular_observable(branches)
    density_rows.append({
        "subset_type": "modular_A",
        "subset": ",".join(str(x) for x in A),
        **density_diagnostics(rho_a),
    })
    for row in density_rows:
        row["condition_kind"] = "baseline"
        row["q"] = ""
        row["epsilon"] = 0.0

    free_graph = build_connected_threshold_graph(W, GRAPH_DENSITY)
    support = edge_set(free_graph)
    if len(support) != 22:
        raise RuntimeError(f"BASELINE_EDGE_COUNT_FAIL got={len(support)} expected=22")

    local_support = localized_edge_set(support)
    if len(local_support) <= 0:
        raise RuntimeError("BASELINE_LOCAL_SUPPORT_EMPTY")
    if any(e not in support for e in local_support):
        raise RuntimeError("BASELINE_LOCAL_SUPPORT_NOT_SUBSET")

    frozen_graph = build_graph_from_edges_and_weights(support, W)
    if not nx.is_connected(frozen_graph):
        raise RuntimeError("BASELINE_FROZEN_GRAPH_DISCONNECTED")

    curvature_graph = compute_ollivier_ricci_internal(frozen_graph, OR_ALPHA)
    kappa_global = curvature_edge_mean(curvature_graph)
    kappa_a = curvature_mean_on_edges(curvature_graph, local_support)

    return Baseline(
        branches=branches,
        W=W,
        W_pairs=W_pairs,
        support=support,
        local_support=local_support,
        m_w=m_w,
        m_k=m_k,
        v_pairs=v_pairs,
        modular_diag=modular_diag,
        kappa_global=kappa_global,
        kappa_a=kappa_a,
        curvature_graph=curvature_graph,
        density_rows=density_rows,
    )


def evaluate_condition(
    q: int,
    epsilon: float,
    baseline: Baseline,
) -> Tuple[
    Dict[str, object],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
]:
    branches = build_analysis_branches(q=q, epsilon=epsilon)
    W, W_pairs, density_rows = mutual_information_matrix_mixed(branches)
    for drow in density_rows:
        drow["condition_kind"] = "perturbed"
        drow["q"] = int(q)
        drow["epsilon"] = float(epsilon)
    m_w = m_w_from_pairs(W_pairs)

    status = "ok"
    m_k = float("nan")
    v_pairs: Dict[Tuple[int, int], float] = {}
    mod_diag: Dict[str, float] = {
        "lambda_min": float("nan"),
        "lambda_max": float("nan"),
        "condition_number": float("nan"),
        "kappa_mean": float("nan"),
        "kappa_std": float("nan"),
    }

    try:
        m_k, v_pairs, mod_diag, rho_a = modular_observable(branches)
        density_rows.append({
            "condition_kind": "perturbed",
            "q": int(q),
            "epsilon": float(epsilon),
            "subset_type": "modular_A",
            "subset": ",".join(str(x) for x in A),
            **density_diagnostics(rho_a),
        })
    except RuntimeError as exc:
        if str(exc).startswith("POSITIVITY_GATE_FAIL"):
            status = "modular_rank_undefined"
        else:
            raise

    graph_raw = build_graph_from_edges_and_weights(baseline.support, W)
    support_now = edge_set(graph_raw)
    if support_now != baseline.support:
        raise RuntimeError("FROZEN_SUPPORT_IDENTITY_FAIL")
    local_support_now = localized_edge_set(support_now)
    if local_support_now != baseline.local_support:
        raise RuntimeError("FROZEN_LOCAL_SUPPORT_IDENTITY_FAIL")
    if not nx.is_connected(graph_raw):
        raise RuntimeError("PERTURBED_FROZEN_GRAPH_DISCONNECTED")

    graph_curv = compute_ollivier_ricci_internal(graph_raw, OR_ALPHA)
    kappa_global = curvature_edge_mean(graph_curv)
    kappa_a = curvature_mean_on_edges(graph_curv, baseline.local_support)
    gdiag = graph_diagnostics(graph_raw)

    row: Dict[str, object] = {
        "q": int(q),
        "epsilon": float(epsilon),
        "distance_to_A": int(distance_to_A(q)),
        "status": status,
        "lambda_min": float(mod_diag["lambda_min"]),
        "lambda_max": float(mod_diag["lambda_max"]),
        "condition_number": float(mod_diag["condition_number"]),
        "M_K": float(m_k),
        "delta_M_K": float(m_k - baseline.m_k) if np.isfinite(m_k) else float("nan"),
        "M_W": float(m_w),
        "delta_M_W": float(m_w - baseline.m_w),
        "kappa_A": float(kappa_a),
        "delta_kappa_A": float(kappa_a - baseline.kappa_a),
        "kappa_global": float(kappa_global),
        "delta_kappa_global": float(kappa_global - baseline.kappa_global),
        "global_support_count": len(baseline.support),
        "local_support_count": len(baseline.local_support),
        "graph_connected": bool(nx.is_connected(graph_raw)),
        "global_support_matches_baseline": bool(support_now == baseline.support),
        "local_support_matches_baseline": bool(local_support_now == baseline.local_support),
        "similarity_min": gdiag["similarity_min"],
        "similarity_max": gdiag["similarity_max"],
        "similarity_norm_min": gdiag["similarity_norm_min"],
        "similarity_norm_max": gdiag["similarity_norm_max"],
        "length_min": gdiag["length_min"],
        "length_max": gdiag["length_max"],
        "psi_plus_norm": float(np.real(np.vdot(branches.psi_plus, branches.psi_plus))),
        "psi_minus_norm": float(np.real(np.vdot(branches.psi_minus, branches.psi_minus))),
    }

    modular_rows: List[Dict[str, object]] = []
    for pair in A_PAIRS:
        value = float(v_pairs[pair]) if pair in v_pairs else float("nan")
        base_value = float(baseline.v_pairs[pair])
        modular_rows.append({
            "q": int(q),
            "epsilon": float(epsilon),
            "i": pair[0],
            "j": pair[1],
            "v_ij": value,
            "baseline_v_ij": base_value,
            "delta_v_ij": value - base_value if np.isfinite(value) else float("nan"),
            "status": status,
        })

    information_rows: List[Dict[str, object]] = []
    for pair in sorted(W_pairs):
        information_rows.append({
            "q": int(q),
            "epsilon": float(epsilon),
            "i": pair[0],
            "j": pair[1],
            "in_A_pair": bool(pair in A_PAIRS),
            "incident_to_A": bool(pair[0] in A or pair[1] in A),
            "W_ij": float(W_pairs[pair]),
            "baseline_W_ij": float(baseline.W_pairs[pair]),
            "delta_W_ij": float(W_pairs[pair] - baseline.W_pairs[pair]),
        })

    curvature_rows: List[Dict[str, object]] = []
    base_edge_data = {
        tuple(sorted((int(u), int(v)))): data
        for u, v, data in baseline.curvature_graph.edges(data=True)
    }
    local_support_set = set(baseline.local_support)
    for u, v, data in graph_curv.edges(data=True):
        e = tuple(sorted((int(u), int(v))))
        base = base_edge_data[e]
        curvature_rows.append({
            "q": int(q),
            "epsilon": float(epsilon),
            "u": e[0],
            "v": e[1],
            "in_local_support": bool(e in local_support_set),
            "similarity": float(data["similarity"]),
            "similarity_norm": float(data["similarity_norm"]),
            "length": float(data["length"]),
            "ricci_curvature": float(data["ricciCurvature"]),
            "baseline_ricci_curvature": float(base["ricciCurvature"]),
            "delta_ricci_curvature": float(data["ricciCurvature"] - base["ricciCurvature"]),
        })

    return row, modular_rows, information_rows, curvature_rows, density_rows


# =============================================================================
# Frozen confirmatory statistics
# =============================================================================

def _pearson_on_ranks(rx: np.ndarray, ry: np.ndarray) -> float:
    ax = rx - np.mean(rx)
    ay = ry - np.mean(ry)
    denom = float(np.linalg.norm(ax) * np.linalg.norm(ay))
    if denom <= 0.0:
        return float("nan")
    return float(np.dot(ax, ay) / denom)


def spearman_stat(x: np.ndarray, y: np.ndarray) -> float:
    return _pearson_on_ranks(
        np.asarray(rankdata(x, method="average"), dtype=float),
        np.asarray(rankdata(y, method="average"), dtype=float),
    )


def stratified_permutation_test(
    rows: Sequence[Dict[str, object]],
    target_key: str,
    *,
    n_perm: int = PERMUTATIONS,
    seed: int = PERMUTATION_SEED,
) -> Dict[str, object]:
    x = np.asarray([float(r["delta_M_K"]) for r in rows], dtype=float)
    y = np.asarray([float(r[target_key]) for r in rows], dtype=float)
    eps = np.asarray([float(r["epsilon"]) for r in rows], dtype=float)

    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return {
            "status": "undefined_nonfinite",
            "target": target_key,
            "observed_spearman": None,
            "p_perm_two_sided": None,
            "n_permutations": 0,
            "seed": None,
            "confirmatory_test_performed": False,
        }

    rx = np.asarray(rankdata(x, method="average"), dtype=float)
    ry = np.asarray(rankdata(y, method="average"), dtype=float)
    observed = _pearson_on_ranks(rx, ry)
    rng = np.random.default_rng(seed)
    extreme = 0

    strata = [np.where(eps == e)[0] for e in EPSILONS]

    for _ in range(n_perm):
        rxp = rx.copy()
        for idx in strata:
            rxp[idx] = rx[rng.permutation(idx)]
        stat = _pearson_on_ranks(rxp, ry)
        if abs(stat) >= abs(observed):
            extreme += 1

    p_value = (1 + extreme) / (n_perm + 1)
    return {
        "status": "ok",
        "target": target_key,
        "observed_spearman": float(observed),
        "p_perm_two_sided": float(p_value),
        "n_permutations": int(n_perm),
        "seed": int(seed),
        "extreme_count": int(extreme),
        "confirmatory_test_performed": True,
    }


def leave_one_epsilon_out(
    rows: Sequence[Dict[str, object]],
    target_key: str,
) -> Dict[str, object]:
    full_x = np.asarray([float(r["delta_M_K"]) for r in rows], dtype=float)
    full_y = np.asarray([float(r[target_key]) for r in rows], dtype=float)
    full = spearman_stat(full_x, full_y)

    vals: List[Dict[str, object]] = []
    for excluded in EPSILONS:
        subset = [r for r in rows if float(r["epsilon"]) != float(excluded)]
        x = np.asarray([float(r["delta_M_K"]) for r in subset], dtype=float)
        y = np.asarray([float(r[target_key]) for r in subset], dtype=float)
        rho = spearman_stat(x, y)
        vals.append({
            "excluded_epsilon": float(excluded),
            "spearman": float(rho),
            "same_sign_as_full": bool(np.sign(rho) == np.sign(full)),
        })

    arr = np.asarray([float(v["spearman"]) for v in vals], dtype=float)
    return {
        "target": target_key,
        "full_spearman": float(full),
        "leave_one_level_out": vals,
        "all_same_sign_as_full": bool(all(bool(v["same_sign_as_full"]) for v in vals)),
        "min_spearman": float(arr.min()),
        "max_spearman": float(arr.max()),
    }


def classify_confirmatory(
    rows: Sequence[Dict[str, object]],
    h1: Dict[str, object],
    h2: Dict[str, object],
) -> str:
    if any(str(r["status"]) != "ok" for r in rows):
        return "LOCALIZED_CURVATURE_CONFIRMATORY_TEST_UNDEFINED_DUE_TO_MODULAR_RANK"

    p1_raw = h1.get("p_perm_two_sided")
    p1 = float(p1_raw) if p1_raw is not None else float("nan")
    if not np.isfinite(p1) or p1 > SIGNIFICANCE_ALPHA:
        return "CROSS_SIZE_MODULAR_INFORMATION_ALIGNMENT_NOT_SUPPORTED"

    p2_raw = h2.get("p_perm_two_sided")
    p2 = float(p2_raw) if p2_raw is not None else float("nan")
    if not np.isfinite(p2) or p2 > SIGNIFICANCE_ALPHA:
        return "MODULAR_INFORMATION_ALIGNMENT_SUPPORTED_LOCAL_CURVATURE_ALIGNMENT_NOT_SUPPORTED"

    return "LOCALIZED_MODULAR_CURVATURE_RESPONSE_SUPPORTED"


# =============================================================================
# Baseline and full output writers
# =============================================================================

def frozen_config() -> Dict[str, object]:
    return {
        "N": N,
        "J": J,
        "h": H_FIELD,
        "P": P,
        "dt": DT,
        "nominal_t": NOMINAL_T,
        "subsystem_A": list(A),
        "sites": list(SITES),
        "epsilons": list(EPSILONS),
        "perturbation": "Ry(q, epsilon) applied to |+>^N before both +/-t evolutions",
        "analysis_state": "0.5*(|psi(+t)><psi(+t)| + |psi(-t)><psi(-t)|)",
        "entropy_log_base": "natural",
        "graph_density": GRAPH_DENSITY,
        "graph_target_edges": 22,
        "graph_mode": "frozen_baseline",
        "localized_edge_rule": "baseline support edges with at least one endpoint in A",
        "or_alpha": OR_ALPHA,
        "or_backend": "internal",
        "eps_length": EPS_LENGTH,
        "permutations": PERMUTATIONS,
        "permutation_seed": PERMUTATION_SEED,
        "significance_alpha": SIGNIFICANCE_ALPHA,
    }


def provenance(script_path: Path) -> Dict[str, object]:
    import scipy
    import platform

    return {
        "preregistration_commit": PREREG_COMMIT,
        "preregistration_sha256": PREREG_SHA256,
        "previous_confirmatory_results_commit": PREVIOUS_CONFIRMATORY_COMMIT,
        "previous_experiment_synthesis_commit": PREVIOUS_SYNTHESIS_COMMIT,
        "paper28_source_blob": PAPER28_BLOB,
        "ollivier_ricci_source_blob": OR_BLOB,
        "implementation_script": str(script_path),
        "implementation_script_sha256": sha256_file(script_path),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "networkx_version": nx.__version__,
    }


def baseline_payload(baseline: Baseline) -> Dict[str, object]:
    return {
        "M_W": baseline.m_w,
        "M_K": baseline.m_k,
        "kappa_A": baseline.kappa_a,
        "kappa_global": baseline.kappa_global,
        "psi_plus_norm": float(np.real(np.vdot(baseline.branches.psi_plus, baseline.branches.psi_plus))),
        "psi_minus_norm": float(np.real(np.vdot(baseline.branches.psi_minus, baseline.branches.psi_minus))),
        "modular_diagnostics": baseline.modular_diag,
        "global_support_edges": [list(e) for e in baseline.support],
        "global_support_edge_count": len(baseline.support),
        "local_support_edges": [list(e) for e in baseline.local_support],
        "local_support_edge_count": len(baseline.local_support),
        "graph_diagnostics": graph_diagnostics(baseline.curvature_graph),
        "A_pair_W": {f"{i}-{j}": baseline.W_pairs[(i, j)] for i, j in A_PAIRS},
        "A_pair_v": {f"{i}-{j}": baseline.v_pairs[(i, j)] for i, j in A_PAIRS},
    }


def run_baseline_only(outdir: Path, script_path: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    baseline = evaluate_baseline()

    write_json(frozen_config(), outdir / "config.json")
    write_json(provenance(script_path), outdir / "provenance.json")
    write_json(baseline_payload(baseline), outdir / "baseline.json")
    write_csv(baseline.density_rows, outdir / "density_diagnostics.csv")

    print("BASELINE_ONLY_PASS")
    print(f"lambda_min={baseline.modular_diag['lambda_min']:.17g}")
    print(f"condition_number={baseline.modular_diag['condition_number']:.17g}")
    print(f"M_W={baseline.m_w:.17g}")
    print(f"M_K={baseline.m_k:.17g}")
    print(f"kappa_A={baseline.kappa_a:.17g}")
    print(f"kappa_global={baseline.kappa_global:.17g}")
    print(f"global_support_edge_count={len(baseline.support)}")
    print(f"local_support_edge_count={len(baseline.local_support)}")


def run_confirmatory(outdir: Path, script_path: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    baseline = evaluate_baseline()
    rows: List[Dict[str, object]] = []
    modular_rows: List[Dict[str, object]] = []
    information_rows: List[Dict[str, object]] = []
    curvature_rows: List[Dict[str, object]] = []
    density_rows: List[Dict[str, object]] = list(baseline.density_rows)

    for epsilon in EPSILONS:
        for q in SITES:
            row, mod, info, curv, dens = evaluate_condition(q, epsilon, baseline)
            rows.append(row)
            modular_rows.extend(mod)
            information_rows.extend(info)
            curvature_rows.extend(curv)
            density_rows.extend(dens)
            print(
                f"condition q={q} epsilon={epsilon:.3f} "
                f"status={row['status']} "
                f"dMK={float(row['delta_M_K']):+.8e} "
                f"dMW={float(row['delta_M_W']):+.8e} "
                f"dKappaA={float(row['delta_kappa_A']):+.8e} "
                f"dKappaGlobal={float(row['delta_kappa_global']):+.8e}"
            )

    write_json(frozen_config(), outdir / "config.json")
    write_json(provenance(script_path), outdir / "provenance.json")
    write_json(baseline_payload(baseline), outdir / "baseline.json")
    write_csv(rows, outdir / "per_condition.csv")
    write_csv(modular_rows, outdir / "modular_pairwise.csv")
    write_csv(information_rows, outdir / "information_pairwise.csv")
    write_csv(curvature_rows, outdir / "curvature_edges.csv")
    write_csv(density_rows, outdir / "density_diagnostics.csv")

    rank_undefined = any(str(r["status"]) != "ok" for r in rows)

    if rank_undefined:
        h1 = {
            "status": "undefined_due_to_modular_rank",
            "target": "delta_M_W",
            "observed_spearman": None,
            "p_perm_two_sided": None,
            "n_permutations": 0,
            "seed": None,
            "confirmatory_test_performed": False,
        }
        h2 = {
            "status": "undefined_due_to_modular_rank",
            "target": "delta_kappa_A",
            "observed_spearman": None,
            "p_perm_two_sided": None,
            "n_permutations": 0,
            "seed": None,
            "confirmatory_test_performed": False,
        }
        robustness: Dict[str, object] = {}
    else:
        h1 = stratified_permutation_test(rows, "delta_M_W")
        robustness = {}
        p1 = float(h1["p_perm_two_sided"])

        if np.isfinite(p1) and p1 <= SIGNIFICANCE_ALPHA:
            h2 = stratified_permutation_test(rows, "delta_kappa_A")
            robustness["H1"] = leave_one_epsilon_out(rows, "delta_M_W")

            p2 = float(h2["p_perm_two_sided"])
            if np.isfinite(p2) and p2 <= SIGNIFICANCE_ALPHA:
                robustness["H2"] = leave_one_epsilon_out(rows, "delta_kappa_A")
        else:
            h2 = {
                "status": "not_performed_h1_failed",
                "target": "delta_kappa_A",
                "observed_spearman": None,
                "p_perm_two_sided": None,
                "n_permutations": 0,
                "seed": None,
                "confirmatory_test_performed": False,
            }

    classification = classify_confirmatory(rows, h1, h2)

    permutation_payload = {
        "H1_modular_information": h1,
        "H2_modular_local_curvature": h2,
        "hierarchical_testing": True,
    }
    write_json(permutation_payload, outdir / "permutation_results.json")

    summary = {
        "classification": classification,
        "n_conditions": len(rows),
        "n_ok": sum(str(r["status"]) == "ok" for r in rows),
        "n_modular_rank_undefined": sum(str(r["status"]) != "ok" for r in rows),
        "H1": h1,
        "H2": h2,
        "robustness": robustness,
        "global_curvature_control_is_confirmatory": False,
        "claim_boundary": {
            "W_to_OR_is_pipeline_construction": True,
            "localized_curvature_definition_frozen_pre_N12": True,
            "independent_curvature_source_claim": False,
            "Einstein_equation_claim": False,
            "continuum_claim": False,
            "universality_claim": False,
        },
    }
    write_json(summary, outdir / "summary.json")

    print(classification)
    print(json.dumps({"H1": h1, "H2": h2}, indent=2, sort_keys=True, allow_nan=False))


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Preregistered fresh N=12 TFIM modular / MI / localized Ollivier–Ricci response benchmark."
    )
    parser.add_argument(
        "--mode",
        choices=("baseline", "confirmatory"),
        default="baseline",
        help="Default is baseline-only. Confirmatory mode executes the frozen 36-condition fresh N=12 grid.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output directory. Defaults are mode-specific under this experiment directory.",
    )
    parser.add_argument(
        "--confirm-unblind",
        default="",
        help="Required exact safety token for confirmatory execution.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_path = Path(__file__).resolve()
    experiment_dir = script_path.parents[1]

    if args.mode == "baseline":
        outdir = (
            Path(args.out)
            if args.out is not None
            else experiment_dir / "results" / "baseline_preunblind"
        )
        run_baseline_only(outdir, script_path)
        return

    if args.confirm_unblind != CONFIRM_TOKEN:
        raise RuntimeError(
            "CONFIRMATORY_EXECUTION_BLOCKED: pass "
            f"--confirm-unblind {CONFIRM_TOKEN}"
        )

    outdir = (
        Path(args.out)
        if args.out is not None
        else experiment_dir / "results" / "confirmatory_localized_v1"
    )
    run_confirmatory(outdir, script_path)


if __name__ == "__main__":
    main()

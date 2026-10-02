#!/usr/bin/env python3
"""
Paper 29 — Information Loss and the Arrow of Modular Time
Prospectively frozen implementation.

IMPORTANT:
- This script implements the preregistered protocol.
- Do not modify primary metrics, grid, pair ordering, controls, or exclusion rules after unblinding.
- Preserve null, negative, failed, and undefined results.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------
# Frozen protocol constants
# ---------------------------------------------------------------------

J = 1.0
H_FIELD = 1.0
P = 6
DT = 0.35
N_PRIMARY = 8
A_SIZE = 4
A_KEEP = (0, 1, 2, 3)
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
S_GRID = np.arange(0.0, 3.0 + 1e-12, 0.1)

STATE_NORM_TOL = 1e-12
ENTROPY_EIG_CUTOFF = 1e-15
STEP_EPS = 1e-12

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"X": X, "Y": Y, "Z": Z}

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "paper29_arrow_of_modular_time_v1"


# ---------------------------------------------------------------------
# Paper 28 state preparation
# ---------------------------------------------------------------------

def plus_state(n: int) -> np.ndarray:
    state = np.array([1.0 + 0.0j])
    q = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)
    for _ in range(n):
        state = np.kron(state, q)
    return state


def one_gate(state: np.ndarray, gate: np.ndarray, q: int, n: int) -> np.ndarray:
    tensor = np.moveaxis(state.reshape([2] * n), q, 0)
    updated = np.tensordot(gate, tensor, axes=([1], [0]))
    return np.moveaxis(updated, 0, q).reshape(-1)


def x_half(state: np.ndarray, dt: float, n: int) -> np.ndarray:
    gate = np.cos(H_FIELD * dt / 2.0) * I2 + 1j * np.sin(H_FIELD * dt / 2.0) * X
    for q in range(n):
        state = one_gate(state, gate, q, n)
    return state


def zz_layer(state: np.ndarray, dt: float, n: int) -> np.ndarray:
    out = state.copy()
    idx = np.arange(2 ** n, dtype=np.uint64)
    energy = np.zeros(2 ** n, dtype=float)
    for q in range(n - 1):
        bq = (idx >> np.uint64(n - 1 - q)) & np.uint64(1)
        br = (idx >> np.uint64(n - 2 - q)) & np.uint64(1)
        zq = 1.0 - 2.0 * bq.astype(float)
        zr = 1.0 - 2.0 * br.astype(float)
        energy += zq * zr
    return out * np.exp(1j * J * dt * energy)


def evolve_tfim(n: int, dt: float) -> np.ndarray:
    state = plus_state(n)
    for _ in range(P):
        state = x_half(state, dt, n)
        state = zz_layer(state, dt, n)
        state = x_half(state, dt, n)
    norm = float(np.real(np.vdot(state, state)))
    if abs(norm - 1.0) > STATE_NORM_TOL:
        raise RuntimeError("STATE_NORM_DRIFT")
    return state


def initial_density(n: int) -> np.ndarray:
    psi_plus = evolve_tfim(n, DT)
    psi_minus = evolve_tfim(n, -DT)
    rho = 0.5 * (
        np.outer(psi_plus, psi_plus.conj())
        + np.outer(psi_minus, psi_minus.conj())
    )
    return hermitize(rho)


# ---------------------------------------------------------------------
# Generic density-matrix utilities
# ---------------------------------------------------------------------

def hermitize(rho: np.ndarray) -> np.ndarray:
    return 0.5 * (rho + rho.conj().T)


def reduced_density(rho: np.ndarray, keep: tuple[int, ...], n: int) -> np.ndarray:
    keep = tuple(keep)
    rest = tuple(q for q in range(n) if q not in keep)
    perm = keep + rest + tuple(n + q for q in keep) + tuple(n + q for q in rest)
    tensor = np.transpose(rho.reshape([2] * (2 * n)), perm)
    dk = 2 ** len(keep)
    dr = 2 ** len(rest)
    out = np.einsum("abcb->ac", tensor.reshape(dk, dr, dk, dr), optimize=True)
    return hermitize(out)


def eigvals_hermitian(rho: np.ndarray) -> np.ndarray:
    return np.real(np.linalg.eigvalsh(hermitize(rho)))


def entropy(rho: np.ndarray) -> float:
    vals = np.clip(eigvals_hermitian(rho), 0.0, None)
    vals = vals[vals > ENTROPY_EIG_CUTOFF]
    return float(-np.sum(vals * np.log(vals)))


def trace_real(rho: np.ndarray) -> float:
    return float(np.real(np.trace(rho)))


def positivity_min(rho: np.ndarray) -> float:
    return float(eigvals_hermitian(rho).min())


def validate_density(rho: np.ndarray, trace_tol: float = 1e-12, positivity_tol: float = 1e-12) -> None:
    tr = trace_real(rho)
    lam_min = positivity_min(rho)
    if abs(tr - 1.0) > trace_tol:
        raise RuntimeError(f"TRACE_PRESERVATION_FAIL:{tr:.17e}")
    if lam_min < -positivity_tol:
        raise RuntimeError(f"DENSITY_POSITIVITY_FAIL:{lam_min:.17e}")


# ---------------------------------------------------------------------
# Irreversible dephasing channel
# ---------------------------------------------------------------------

def z_on_qubit(q: int, n: int) -> np.ndarray:
    op = np.array([[1.0 + 0.0j]])
    for k in range(n):
        op = np.kron(op, Z if k == q else I2)
    return op


def dephase_one_qubit(rho: np.ndarray, s: float, q: int, n: int) -> np.ndarray:
    zq = z_on_qubit(q, n)
    a = 0.5 * (1.0 + np.exp(-s))
    b = 0.5 * (1.0 - np.exp(-s))
    return a * rho + b * (zq @ rho @ zq)


_HAMMING_CACHE = {}

def hamming_distance_matrix(n: int) -> np.ndarray:
    if n in _HAMMING_CACHE:
        return _HAMMING_CACHE[n]
    dim = 2 ** n
    idx = np.arange(dim, dtype=np.uint16)
    xor = np.bitwise_xor(idx[:, None], idx[None, :])
    lut = np.array([int(i).bit_count() for i in range(256)], dtype=np.uint8)
    h = lut[(xor & 0x00FF).astype(np.uint8)]
    if n > 8:
        h = h + lut[((xor >> 8) & 0x00FF).astype(np.uint8)]
    h = h.astype(float)
    _HAMMING_CACHE[n] = h
    return h


def dephase_all(rho0: np.ndarray, s: float, n: int) -> np.ndarray:
    # Independent local dephasing multiplies rho_xy by exp[-s * Hamming(x,y)].
    # This is exactly equivalent to applying the frozen single-qubit channel
    # to every qubit, while avoiding repeated dense Z rho Z multiplications.
    h = hamming_distance_matrix(n)
    rho = rho0 * np.exp(-float(s) * h)
    rho = hermitize(rho)
    validate_density(rho)
    return rho


# ---------------------------------------------------------------------
# Reference-information coordinate
# ---------------------------------------------------------------------

def relative_entropy_to_maximally_mixed(rho: np.ndarray, n: int) -> float:
    # D(rho || I/d) = log(d) - S(rho) = n log 2 - S(rho)
    return float(n * np.log(2.0) - entropy(rho))


# ---------------------------------------------------------------------
# Informational geometry
# ---------------------------------------------------------------------

def pair_mutual_information(rho_a: np.ndarray, a_size: int, pairs=PAIRS) -> np.ndarray:
    values = []
    for i, j in pairs:
        rho_i = reduced_density(rho_a, (i,), a_size)
        rho_j = reduced_density(rho_a, (j,), a_size)
        rho_ij = reduced_density(rho_a, (i, j), a_size)
        values.append(entropy(rho_i) + entropy(rho_j) - entropy(rho_ij))
    return np.asarray(values, dtype=float)


# ---------------------------------------------------------------------
# Modular sector
# ---------------------------------------------------------------------

def local_ops(a_size: int) -> dict[tuple[int, str], np.ndarray]:
    ops = {}
    for q in range(a_size):
        for label, sigma in PAULI.items():
            op = np.array([[1.0 + 0.0j]])
            for k in range(a_size):
                op = np.kron(op, sigma if k == q else I2)
            ops[(q, label)] = op
    return ops


def modular_hamiltonian_normalized(rho_a: np.ndarray):
    vals, vecs = np.linalg.eigh(hermitize(rho_a))
    vals = np.real(vals)
    if vals.min() <= 0.0:
        raise RuntimeError("POSITIVITY_GATE_FAIL")
    kappa = -np.log(vals)
    kappa_mean = float(kappa.mean())
    kappa_std = float(kappa.std(ddof=0))
    if kappa_std <= 0.0:
        raise RuntimeError("MODULAR_STD_ZERO")
    k_raw = vecs @ np.diag(kappa) @ vecs.conj().T
    kappa_norm = (kappa - kappa_mean) / kappa_std
    k_tilde = vecs @ np.diag(kappa_norm) @ vecs.conj().T
    diagnostics = {
        "lambda_min": float(vals.min()),
        "lambda_max": float(vals.max()),
        "condition_number": float(vals.max() / vals.min()),
        "kappa_mean": kappa_mean,
        "kappa_std": kappa_std,
    }
    return hermitize(k_raw), hermitize(k_tilde), diagnostics


def comm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def modular_coefficients(
    rho_a: np.ndarray,
    k_tilde: np.ndarray,
    a_size: int,
    pairs=PAIRS,
) -> np.ndarray:
    ops = local_ops(a_size)
    coeffs = []
    for i, j in pairs:
        total = 0.0
        for lab_a in "XYZ":
            inner = comm(k_tilde, ops[(i, lab_a)])
            for lab_b in "XYZ":
                qop = comm(inner, ops[(j, lab_b)])
                total += float(np.real(np.trace(rho_a @ qop.conj().T @ qop)))
        coeffs.append(total / 9.0)
    coeffs = np.asarray(coeffs, dtype=float)
    if np.any(coeffs < -1e-12):
        raise RuntimeError("NEGATIVE_MODULAR_COEFFICIENT")
    return np.clip(coeffs, 0.0, None)


# ---------------------------------------------------------------------
# Normalization and directional comparison
# ---------------------------------------------------------------------

def normalize_vector(vec: np.ndarray, name: str):
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0:
        raise RuntimeError(f"ZERO_VECTOR_NORM:{name}")
    return vec / norm, norm


def directional_step(prev_w_hat, next_w_hat, prev_v_hat, next_v_hat):
    dw = next_w_hat - prev_w_hat
    dv = next_v_hat - prev_v_hat
    nw = float(np.linalg.norm(dw))
    nv = float(np.linalg.norm(dv))
    if nw <= STEP_EPS or nv <= STEP_EPS:
        return {
            "delta_w": dw,
            "delta_v": dv,
            "d_W": nw,
            "d_mod": nv,
            "chi": np.nan,
            "status": "DIRECTION_UNDEFINED_ZERO_STEP",
        }
    chi = float(np.dot(dw, dv) / (nw * nv))
    chi = float(np.clip(chi, -1.0, 1.0))
    return {
        "delta_w": dw,
        "delta_v": dv,
        "d_W": nw,
        "d_mod": nv,
        "chi": chi,
        "status": "VALID",
    }


# ---------------------------------------------------------------------
# Reversible TFIM control
# ---------------------------------------------------------------------

def full_tfim_hamiltonian(n: int) -> np.ndarray:
    dim = 2 ** n
    h = np.zeros((dim, dim), dtype=complex)

    # -J sum Z_i Z_{i+1}
    for q in range(n - 1):
        op = np.array([[1.0 + 0.0j]])
        for k in range(n):
            if k == q or k == q + 1:
                op = np.kron(op, Z)
            else:
                op = np.kron(op, I2)
        h += -J * op

    # -h sum X_i
    for q in range(n):
        op = np.array([[1.0 + 0.0j]])
        for k in range(n):
            op = np.kron(op, X if k == q else I2)
        h += -H_FIELD * op

    return hermitize(h)


def unitary_control_state(rho0: np.ndarray, h: np.ndarray, s: float) -> np.ndarray:
    vals, vecs = np.linalg.eigh(h)
    phase = np.exp(-1j * vals * s)
    u = vecs @ np.diag(phase) @ vecs.conj().T
    rho = u @ rho0 @ u.conj().T
    rho = hermitize(rho)
    validate_density(rho)
    return rho


# ---------------------------------------------------------------------
# Static protocol manifest
# ---------------------------------------------------------------------

def frozen_manifest() -> dict:
    return {
        "status": "IMPLEMENTATION_NOT_YET_UNBLINDED",
        "paper": 29,
        "title": "Information Loss and the Arrow of Modular Time",
        "primary_system": {
            "N": N_PRIMARY,
            "A_keep": list(A_KEEP),
            "pair_order": [list(p) for p in PAIRS],
        },
        "paper28_state_preparation": {
            "J": J,
            "h": H_FIELD,
            "p": P,
            "dt": DT,
            "nominal_t": P * DT,
            "state": "time-reversal mixture of Strang-evolved |+>^N",
        },
        "irreversible_channel": {
            "type": "independent local dephasing",
            "single_qubit_formula": "(1+exp(-s))/2 * rho + (1-exp(-s))/2 * Z rho Z",
            "s_grid": [float(x) for x in S_GRID],
        },
        "reference": "I / 2^N",
        "information_coordinate": "Q(s)=D(rho(s)||I/2^N)=N log 2 - S(rho(s))",
        "sigma_step": "Sigma_n=Q(s_n)-Q(s_{n+1})",
        "geometry": "w(s)=[I01,I02,I03,I12,I13,I23]",
        "modular": {
            "K": "-log(rho_A)",
            "normalization": "(K-mean(kappa)I)/std(kappa)",
            "primary_profile": "v_ij=sqrt(A_ij) from exact double commutator coefficient",
        },
        "primary_endpoint": "mean chi over valid consecutive steps",
        "secondary": [
            "fraction chi > 0",
            "Spearman(d_W,d_mod)",
            "Spearman(Sigma,d_W)",
            "Spearman(Sigma,d_mod)",
            "cumulative normalized distances",
            "raw magnitudes",
        ],
        "zero_step_rule": {
            "epsilon": STEP_EPS,
            "status": "DIRECTION_UNDEFINED_ZERO_STEP",
            "chi": "NaN",
        },
        "numerical": {
            "state_norm_tol": STATE_NORM_TOL,
            "entropy_eig_cutoff": ENTROPY_EIG_CUTOFF,
            "strict_modular_positivity": True,
        },
        "controls": [
            "s=0 Paper 28 recovery",
            "reversible TFIM unitary trajectory",
            "raw versus normalized quantities",
            "N10/A4 replication",
        ],
    }


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------
# Result-producing functions
# ---------------------------------------------------------------------

def analyze_state(rho_global: np.ndarray, n: int):
    rho_a = reduced_density(rho_global, A_KEEP, n)
    validate_density(rho_a)

    q_info = relative_entropy_to_maximally_mixed(rho_global, n)

    w = pair_mutual_information(rho_a, A_SIZE, PAIRS)
    w_hat, m_w = normalize_vector(w, "w")

    k_raw, k_tilde, modular_diag = modular_hamiltonian_normalized(rho_a)
    a_coeff = modular_coefficients(rho_a, k_tilde, A_SIZE, PAIRS)
    v = np.sqrt(a_coeff)
    v_hat, m_mod = normalize_vector(v, "v")

    return {
        "Q": q_info,
        "rho_A": rho_a,
        "w": w,
        "w_hat": w_hat,
        "M_W": m_w,
        "K_raw": k_raw,
        "K_tilde": k_tilde,
        "A_coeff": a_coeff,
        "v": v,
        "v_hat": v_hat,
        "M_mod": m_mod,
        "modular_diag": modular_diag,
    }


PAPER28_EXPECTED = {
    8: {
        "Delta_R_primary": 0.12292356524542974,
        "cosine_W": 0.9660284140131971,
        "residual_W": 0.258435878544654,
        "cosine_chain": 0.9244268357377549,
        "residual_chain": 0.38135944379008374,
        "lambda_min": 0.00031808556235593744,
        "lambda_max": 0.38457742531689254,
        "condition_number": 1209.0376641695884,
        "kappa_mean": 4.587517093532144,
        "kappa_std": 2.163383760368966,
    },
    10: {
        "Delta_R_primary": 0.12419919426286236,
        "cosine_W": 0.9668328093185072,
        "residual_W": 0.2554100993016585,
        "cosine_chain": 0.9251468987352471,
        "residual_chain": 0.37960929356452083,
        "lambda_min": 0.0003142460094062414,
        "lambda_max": 0.3829563020160241,
        "condition_number": 1218.6512813308552,
        "kappa_mean": 4.583362007862227,
        "kappa_std": 2.167375790925045,
    },
}
PAPER28_RECOVERY_ATOL = 1e-12


def paper28_metric(v: np.ndarray, g: np.ndarray):
    nv = np.linalg.norm(v)
    dot = float(v @ g)
    scale = max(0.0, dot / float(g @ g))
    cosine = dot / (nv * np.linalg.norm(g))
    residual = float(np.linalg.norm(v - scale * g) / nv)
    return float(cosine), float(residual)


def paper28_recovery(result: dict, n: int) -> dict:
    g = np.array([1.0 if j == i + 1 else 0.0 for i, j in PAIRS], dtype=float)
    cosine_w, residual_w = paper28_metric(result["v"], result["w"])
    cosine_chain, residual_chain = paper28_metric(result["v"], g)
    observed = {
        "Delta_R_primary": residual_chain - residual_w,
        "cosine_W": cosine_w,
        "residual_W": residual_w,
        "cosine_chain": cosine_chain,
        "residual_chain": residual_chain,
        **result["modular_diag"],
    }
    expected = PAPER28_EXPECTED[n]
    deltas = {k: float(observed[k] - expected[k]) for k in expected}
    passed = all(abs(deltas[k]) <= PAPER28_RECOVERY_ATOL for k in expected)
    return {
        "passed": bool(passed),
        "atol": PAPER28_RECOVERY_ATOL,
        "observed": observed,
        "expected": expected,
        "delta": deltas,
    }


def spearman_no_scipy(x: np.ndarray, y: np.ndarray) -> float:
    # Frozen deterministic rank implementation with average ranks for ties.
    def ranks(a):
        a = np.asarray(a, dtype=float)
        order = np.argsort(a, kind="mergesort")
        r = np.empty(len(a), dtype=float)
        i = 0
        while i < len(a):
            j = i + 1
            while j < len(a) and a[order[j]] == a[order[i]]:
                j += 1
            rank = 0.5 * ((i + 1) + j)
            r[order[i:j]] = rank
            i = j
        return r

    rx = ranks(x)
    ry = ranks(y)
    sx = rx.std(ddof=0)
    sy = ry.std(ddof=0)
    if sx == 0.0 or sy == 0.0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def run_primary(n: int):
    rho0 = initial_density(n)
    validate_density(rho0)

    state_rows = []
    state_payloads = []

    for s in S_GRID:
        rho_s = dephase_all(rho0, float(s), n)
        result = analyze_state(rho_s, n)
        state_payloads.append(result)

        row = {
            "s": float(s),
            "Q": result["Q"],
            "M_W": result["M_W"],
            "M_mod": result["M_mod"],
            **result["modular_diag"],
        }

        for idx, pair in enumerate(PAIRS):
            label = f"{pair[0]}{pair[1]}"
            row[f"W_{label}"] = float(result["w"][idx])
            row[f"W_hat_{label}"] = float(result["w_hat"][idx])
            row[f"A_{label}"] = float(result["A_coeff"][idx])
            row[f"v_{label}"] = float(result["v"][idx])
            row[f"v_hat_{label}"] = float(result["v_hat"][idx])

        state_rows.append(row)

    w0_hat = state_payloads[0]["w_hat"]
    v0_hat = state_payloads[0]["v_hat"]
    for row, payload in zip(state_rows, state_payloads):
        row["D_W_from_s0"] = float(np.linalg.norm(payload["w_hat"] - w0_hat))
        row["D_mod_from_s0"] = float(np.linalg.norm(payload["v_hat"] - v0_hat))

    recovery = paper28_recovery(state_payloads[0], n)
    if not recovery["passed"]:
        raise RuntimeError(f"PAPER28_S0_RECOVERY_FAIL_N{n}")

    transition_rows = []
    valid_chi = []

    for k in range(len(S_GRID) - 1):
        cur = state_payloads[k]
        nxt = state_payloads[k + 1]

        sigma_n = float(cur["Q"] - nxt["Q"])
        step = directional_step(cur["w_hat"], nxt["w_hat"], cur["v_hat"], nxt["v_hat"])

        row = {
            "n": k,
            "s_from": float(S_GRID[k]),
            "s_to": float(S_GRID[k + 1]),
            "Sigma": sigma_n,
            "d_W": step["d_W"],
            "d_mod": step["d_mod"],
            "chi": step["chi"],
            "status": step["status"],
        }

        for idx, pair in enumerate(PAIRS):
            label = f"{pair[0]}{pair[1]}"
            row[f"delta_W_hat_{label}"] = float(step["delta_w"][idx])
            row[f"delta_v_hat_{label}"] = float(step["delta_v"][idx])

        if step["status"] == "VALID":
            valid_chi.append(step["chi"])

        transition_rows.append(row)

    valid_chi = np.asarray(valid_chi, dtype=float)
    d_w = np.asarray([r["d_W"] for r in transition_rows], dtype=float)
    d_mod = np.asarray([r["d_mod"] for r in transition_rows], dtype=float)
    sigma = np.asarray([r["Sigma"] for r in transition_rows], dtype=float)

    summary = {
        "N": n,
        "A_size": A_SIZE,
        "n_states": len(state_rows),
        "n_transitions": len(transition_rows),
        "N_valid": int(len(valid_chi)),
        "N_undefined": int(len(transition_rows) - len(valid_chi)),
        "mean_chi": float(np.mean(valid_chi)) if len(valid_chi) else float("nan"),
        "fraction_chi_positive": float(np.mean(valid_chi > 0.0)) if len(valid_chi) else float("nan"),
        "spearman_dW_dmod": spearman_no_scipy(d_w, d_mod),
        "spearman_Sigma_dW": spearman_no_scipy(sigma, d_w),
        "spearman_Sigma_dmod": spearman_no_scipy(sigma, d_mod),
        "Sigma_min": float(np.min(sigma)),
        "Sigma_max": float(np.max(sigma)),
        "Sigma_negative_count": int(np.sum(sigma < -1e-12)),
        "paper28_s0_recovery": recovery,
    }

    return state_rows, transition_rows, summary, state_payloads


def run_unitary_control(n: int):
    rho0 = initial_density(n)
    h = full_tfim_hamiltonian(n)
    rows = []

    for s in S_GRID:
        rho_u = unitary_control_state(rho0, h, float(s))
        q = relative_entropy_to_maximally_mixed(rho_u, n)
        rows.append({"s": float(s), "Q": q})

    for i in range(len(rows) - 1):
        rows[i]["Sigma_to_next"] = float(rows[i]["Q"] - rows[i + 1]["Q"])
    rows[-1]["Sigma_to_next"] = float("nan")

    sigma = np.asarray([r["Sigma_to_next"] for r in rows[:-1]], dtype=float)
    summary = {
        "N": n,
        "Q_span": float(max(r["Q"] for r in rows) - min(r["Q"] for r in rows)),
        "max_abs_Sigma": float(np.max(np.abs(sigma))),
    }
    return rows, summary


def write_state_npz(path: Path, payloads: list[dict]) -> None:
    np.savez_compressed(
        path,
        s=np.asarray(S_GRID, dtype=float),
        rho_A=np.stack([p["rho_A"] for p in payloads], axis=0),
        K_raw=np.stack([p["K_raw"] for p in payloads], axis=0),
        K_tilde=np.stack([p["K_tilde"] for p in payloads], axis=0),
        w=np.stack([p["w"] for p in payloads], axis=0),
        w_hat=np.stack([p["w_hat"] for p in payloads], axis=0),
        A_coeff=np.stack([p["A_coeff"] for p in payloads], axis=0),
        v=np.stack([p["v"] for p in payloads], axis=0),
        v_hat=np.stack([p["v_hat"] for p in payloads], axis=0),
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Write the frozen manifest first.
    write_json(OUT_DIR / "protocol_manifest.json", frozen_manifest())

    # Result-producing execution begins below this line.
    primary_states, primary_transitions, primary_summary, primary_payloads = run_primary(N_PRIMARY)
    write_csv(OUT_DIR / "primary_states_N8_A4.csv", primary_states)
    write_csv(OUT_DIR / "primary_transitions_N8_A4.csv", primary_transitions)
    write_state_npz(OUT_DIR / "primary_matrices_N8_A4.npz", primary_payloads)
    write_json(OUT_DIR / "primary_summary_N8_A4.json", primary_summary)

    control_rows, control_summary = run_unitary_control(N_PRIMARY)
    write_csv(OUT_DIR / "unitary_control_N8.csv", control_rows)
    write_json(OUT_DIR / "unitary_control_summary_N8.json", control_summary)

    replication_states, replication_transitions, replication_summary, replication_payloads = run_primary(10)
    write_csv(OUT_DIR / "replication_states_N10_A4.csv", replication_states)
    write_csv(OUT_DIR / "replication_transitions_N10_A4.csv", replication_transitions)
    write_state_npz(OUT_DIR / "replication_matrices_N10_A4.npz", replication_payloads)
    write_json(OUT_DIR / "replication_summary_N10_A4.json", replication_summary)

    combined = {
        "protocol": frozen_manifest(),
        "primary_N8_A4": primary_summary,
        "unitary_control_N8": control_summary,
        "replication_N10_A4": replication_summary,
        "claim_boundary": {
            "irreversibility_is_externally_imposed": True,
            "modular_flow_arrow_claim_allowed": False,
            "bup_explains_thermodynamic_arrow_claim_allowed": False,
        },
    }
    write_json(OUT_DIR / "summary.json", combined)

    print("PAPER29_EXECUTION_COMPLETED")


if __name__ == "__main__":
    main()

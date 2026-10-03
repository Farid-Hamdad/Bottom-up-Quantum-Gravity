#!/usr/bin/env python3
"""
Paper 30 — Robustness of Information-Sector Compensation under Amplitude Damping

Prospectively frozen implementation candidate.

IMPORTANT
---------
- Implements the frozen Paper 30 preregistration.
- Do not modify primary metrics, gamma grid, pair ordering, bandwidths,
  thresholds, controls, or success rules after scientific unblinding.
- Preserve null, negative, failed, and undefined outcomes.
- Paper 29 is the discovery study.
- N=12 is the primary prospective Paper 30 validation system.
- N=8 and N=10 are secondary same-size channel replications.
"""

from __future__ import annotations

import csv
import hashlib
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

N_PRIMARY = 12
N_REPLICATIONS = (8, 10)

A_SIZE = 4
A_KEEP = (0, 1, 2, 3)
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

S_GRID = np.arange(0.0, 3.0 + 1e-12, 0.1)
GAMMA_GRID = 1.0 - np.exp(-S_GRID)

BANDWIDTH_FRACTIONS = (0.05, 0.10, 0.20, 0.30, 0.40)
WEIGHT_FLOOR = 1e-14

STATE_NORM_TOL = 1e-12
DENSITY_TRACE_TOL = 1e-12
DENSITY_POSITIVITY_TOL = 1e-12
HERMITICITY_TOL = 1e-12
ENTROPY_EIG_CUTOFF = 1e-15
STEP_EPS = 1e-12

REL_ENT_EPS = 1e-6

ALPHA_MIN = 0.8
ALPHA_MAX = 1.2
RC_MAX = 0.20

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"X": X, "Y": Y, "Z": Z}

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "paper30_amplitude_damping_compensation_v1"
PREREG_PATH = ROOT / "notes" / "PREREGISTRATION_FROZEN.md"


# ---------------------------------------------------------------------
# Paper 28 / Paper 29 state preparation — copied without scientific change
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


def entropy_from_eigvals(vals: np.ndarray) -> float:
    vals = np.clip(np.asarray(vals, dtype=float), 0.0, None)
    vals = vals[vals > ENTROPY_EIG_CUTOFF]
    return float(-np.sum(vals * np.log(vals)))


def entropy(rho: np.ndarray) -> float:
    return entropy_from_eigvals(eigvals_hermitian(rho))


def trace_real(rho: np.ndarray) -> float:
    return float(np.real(np.trace(rho)))


def hermiticity_error(rho: np.ndarray) -> float:
    return float(np.max(np.abs(rho - rho.conj().T)))


def global_density_diagnostics(rho: np.ndarray) -> dict:
    """
    One global eigendecomposition per state supplies both entropy and
    positivity control. This avoids duplicate N=12 diagonalizations.
    """
    tr = trace_real(rho)
    herm_err = hermiticity_error(rho)

    if abs(tr - 1.0) > DENSITY_TRACE_TOL:
        raise RuntimeError(f"TRACE_PRESERVATION_FAIL:{tr:.17e}")
    if herm_err > HERMITICITY_TOL:
        raise RuntimeError(f"HERMITICITY_FAIL:{herm_err:.17e}")

    vals = eigvals_hermitian(rho)
    lam_min = float(vals.min())
    if lam_min < -DENSITY_POSITIVITY_TOL:
        raise RuntimeError(f"DENSITY_POSITIVITY_FAIL:{lam_min:.17e}")

    return {
        "trace": tr,
        "hermiticity_error": herm_err,
        "lambda_min_global": lam_min,
        "entropy_global": entropy_from_eigvals(vals),
    }


def validate_small_density(
    rho: np.ndarray,
    trace_tol: float = DENSITY_TRACE_TOL,
    positivity_tol: float = DENSITY_POSITIVITY_TOL,
) -> None:
    tr = trace_real(rho)
    lam_min = float(eigvals_hermitian(rho).min())
    if abs(tr - 1.0) > trace_tol:
        raise RuntimeError(f"TRACE_PRESERVATION_FAIL:{tr:.17e}")
    if lam_min < -positivity_tol:
        raise RuntimeError(f"DENSITY_POSITIVITY_FAIL:{lam_min:.17e}")


# ---------------------------------------------------------------------
# Local amplitude-damping channel
# ---------------------------------------------------------------------

def amplitude_damp_one_inplace(rho: np.ndarray, gamma: float, q: int, n: int) -> None:
    """
    Apply one-qubit amplitude damping in place.

    For the local ket/bra indices of qubit q:
      rho_00 <- rho_00 + gamma rho_11
      rho_01 <- sqrt(1-gamma) rho_01
      rho_10 <- sqrt(1-gamma) rho_10
      rho_11 <- (1-gamma) rho_11

    rho_11 is read before it is rescaled, so the update is Kraus-equivalent.
    """
    if gamma < -1e-15 or gamma > 1.0 + 1e-15:
        raise ValueError(f"GAMMA_OUT_OF_RANGE:{gamma:.17e}")

    gamma = float(np.clip(gamma, 0.0, 1.0))
    damp = float(np.sqrt(1.0 - gamma))
    survive = 1.0 - gamma

    tensor = rho.reshape([2] * (2 * n))

    s00 = [slice(None)] * (2 * n)
    s01 = [slice(None)] * (2 * n)
    s10 = [slice(None)] * (2 * n)
    s11 = [slice(None)] * (2 * n)

    s00[q], s00[n + q] = 0, 0
    s01[q], s01[n + q] = 0, 1
    s10[q], s10[n + q] = 1, 0
    s11[q], s11[n + q] = 1, 1

    t00 = tuple(s00)
    t01 = tuple(s01)
    t10 = tuple(s10)
    t11 = tuple(s11)

    tensor[t00] += gamma * tensor[t11]
    tensor[t01] *= damp
    tensor[t10] *= damp
    tensor[t11] *= survive


def amplitude_damp_all(rho0: np.ndarray, gamma: float, n: int) -> np.ndarray:
    rho = rho0.copy()
    for q in range(n):
        amplitude_damp_one_inplace(rho, gamma, q, n)
    return rho


# ---------------------------------------------------------------------
# Descriptive amplitude-damping controls
# ---------------------------------------------------------------------

def excitation_content(rho: np.ndarray, n: int) -> float:
    diag = np.real(np.diag(rho))
    idx = np.arange(2 ** n, dtype=np.uint64)

    counts = np.zeros(2 ** n, dtype=float)
    for q in range(n):
        bit = (idx >> np.uint64(n - 1 - q)) & np.uint64(1)
        counts += bit.astype(float)

    return float(np.dot(diag, counts))


def regularized_relative_entropy_from_entropy(
    rho: np.ndarray,
    entropy_rho: float,
    n: int,
    eps: float = REL_ENT_EPS,
) -> float:
    """
    D_eps(rho || sigma_eps), where
      sigma_eps = (1-eps)|0...0><0...0| + eps I / 2^N.

    sigma_eps is diagonal with one eigenvalue lambda0 and a common
    orthogonal eigenvalue lambda_perp, so no matrix logarithm is needed.
    """
    d = 2 ** n
    lambda0 = (1.0 - eps) + eps / d
    lambda_perp = eps / d

    p0 = float(np.real(rho[0, 0]))
    cross_entropy_term = (
        p0 * np.log(lambda0)
        + (1.0 - p0) * np.log(lambda_perp)
    )

    return float(-entropy_rho - cross_entropy_term)


# ---------------------------------------------------------------------
# Mutual-information sector — same construction as Paper 29
# ---------------------------------------------------------------------

def pair_mutual_information(
    rho_a: np.ndarray,
    a_size: int,
    pairs=PAIRS,
) -> np.ndarray:
    values = []
    for i, j in pairs:
        rho_i = reduced_density(rho_a, (i,), a_size)
        rho_j = reduced_density(rho_a, (j,), a_size)
        rho_ij = reduced_density(rho_a, (i, j), a_size)
        values.append(entropy(rho_i) + entropy(rho_j) - entropy(rho_ij))
    return np.asarray(values, dtype=float)


# ---------------------------------------------------------------------
# Modular sector — copied without scientific change from Paper 29
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


def normalize_vector(vec: np.ndarray, name: str):
    norm = float(np.linalg.norm(vec))
    if norm <= STEP_EPS:
        return None, norm, f"{name}_UNDEFINED_ZERO_NORM"
    return vec / norm, norm, "VALID"


def analyze_state(rho_global: np.ndarray, n: int) -> dict:
    global_diag = global_density_diagnostics(rho_global)

    rho_a = reduced_density(rho_global, A_KEEP, n)
    validate_small_density(rho_a)

    w = pair_mutual_information(rho_a, A_SIZE, PAIRS)
    w_hat, m_w, w_status = normalize_vector(w, "W")

    k_raw, k_tilde, modular_diag = modular_hamiltonian_normalized(rho_a)
    a_coeff = modular_coefficients(rho_a, k_tilde, A_SIZE, PAIRS)
    v = np.sqrt(a_coeff)
    v_hat, m_mod, v_status = normalize_vector(v, "MOD")

    entropy_global = global_diag["entropy_global"]
    e_exc = excitation_content(rho_global, n)
    d_eps = regularized_relative_entropy_from_entropy(
        rho_global,
        entropy_global,
        n,
        REL_ENT_EPS,
    )

    return {
        "rho_A": rho_a,
        "w": w,
        "w_hat": w_hat,
        "w_status": w_status,
        "M_W": m_w,
        "K_raw": k_raw,
        "K_tilde": k_tilde,
        "A_coeff": a_coeff,
        "v": v,
        "v_hat": v_hat,
        "v_status": v_status,
        "M_mod": m_mod,
        "entropy_global": entropy_global,
        "E_exc": e_exc,
        "D_eps": d_eps,
        "global_diag": global_diag,
        "modular_diag": modular_diag,
    }


# ---------------------------------------------------------------------
# Paper 28 / Paper 29 s=0 recovery for N8 and N10
# ---------------------------------------------------------------------

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

PAPER28_RECOVERY_ATOL = 1e-10


def paper28_metric(v: np.ndarray, g: np.ndarray):
    nv = np.linalg.norm(v)
    dot = float(v @ g)
    scale = max(0.0, dot / float(g @ g))
    cosine = dot / (nv * np.linalg.norm(g))
    residual = float(np.linalg.norm(v - scale * g) / nv)
    return float(cosine), float(residual)


def paper28_recovery(result: dict, n: int) -> dict:
    if n not in PAPER28_EXPECTED:
        return {
            "applicable": False,
            "passed": None,
            "reason": "NO_FROZEN_PAPER28_REFERENCE_FOR_THIS_N",
        }

    g = np.array(
        [1.0 if j == i + 1 else 0.0 for i, j in PAIRS],
        dtype=float,
    )

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
        "applicable": True,
        "passed": bool(passed),
        "atol": PAPER28_RECOVERY_ATOL,
        "observed": observed,
        "expected": expected,
        "delta": deltas,
    }


# ---------------------------------------------------------------------
# Transition construction
# ---------------------------------------------------------------------

def transition_step(cur: dict, nxt: dict) -> dict:
    if cur["w_hat"] is None or nxt["w_hat"] is None:
        return {
            "d_W": np.nan,
            "d_mod": np.nan,
            "status": "DIRECTION_UNDEFINED_W_ZERO_NORM",
        }

    if cur["v_hat"] is None or nxt["v_hat"] is None:
        return {
            "d_W": np.nan,
            "d_mod": np.nan,
            "status": "DIRECTION_UNDEFINED_MOD_ZERO_NORM",
        }

    dw = nxt["w_hat"] - cur["w_hat"]
    dv = nxt["v_hat"] - cur["v_hat"]

    d_w = float(np.linalg.norm(dw))
    d_mod = float(np.linalg.norm(dv))

    if d_w <= STEP_EPS or d_mod <= STEP_EPS:
        return {
            "d_W": d_w,
            "d_mod": d_mod,
            "status": "DIRECTION_UNDEFINED_ZERO_STEP",
        }

    return {
        "d_W": d_w,
        "d_mod": d_mod,
        "status": "VALID",
    }


# ---------------------------------------------------------------------
# Exact Paper 29 nonparametric detrending implementation,
# with Q replaced by preregistered gamma
# ---------------------------------------------------------------------

def rankdata(a):
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)

    i = 0
    while i < len(a):
        j = i + 1
        while j < len(a) and a[order[j]] == a[order[i]]:
            j += 1
        rank = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = rank
        i = j

    return ranks


def pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def local_linear_loo_predict(x, y, bandwidth):
    """
    Exact Paper 29 LOO local-linear Gaussian-kernel procedure.

    At x_i, point i is excluded. A weighted local line in
    dx = x_j - x_i is fitted; the intercept predicts y at x_i.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if bandwidth <= 0.0:
        raise ValueError("bandwidth must be positive")

    n = len(x)
    pred = np.empty(n, dtype=float)
    effective_weight_sums = np.empty(n, dtype=float)

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False

        dx = x[mask] - x[i]
        yy = y[mask]
        w = np.exp(-0.5 * (dx / bandwidth) ** 2)

        effective_weight_sums[i] = float(np.sum(w))

        if np.sum(w) <= WEIGHT_FLOOR:
            j_local = int(np.argmin(np.abs(dx)))
            pred[i] = float(yy[j_local])
            continue

        x_design = np.column_stack([np.ones_like(dx), dx])
        sqrt_w = np.sqrt(w)
        xw = x_design * sqrt_w[:, None]
        yw = yy * sqrt_w

        beta, *_ = np.linalg.lstsq(xw, yw, rcond=None)
        pred[i] = float(beta[0])

    return pred, effective_weight_sums


# ---------------------------------------------------------------------
# Compensation / quasi-invariant metrics
# ---------------------------------------------------------------------

def zscore(x):
    x = np.asarray(x, dtype=float)
    mu = float(np.mean(x))
    sd = float(np.std(x))

    if not np.isfinite(sd) or sd <= 0.0:
        raise RuntimeError("ZERO_OR_INVALID_RESIDUAL_STD")

    return (x - mu) / sd, mu, sd


def optimal_alpha(z_w, z_m):
    den = float(np.dot(z_w, z_w))
    if den <= 0.0:
        raise RuntimeError("ZERO_ALPHA_DENOMINATOR")
    return float(-np.dot(z_w, z_m) / den)


def quasi_invariant_metrics(z_w, z_m, alpha):
    c = z_m + float(alpha) * z_w
    var_mod = float(np.var(z_m))
    var_c = float(np.var(c))

    if var_mod <= 0.0:
        raise RuntimeError("ZERO_STANDARDIZED_MOD_VARIANCE")

    return {
        "alpha": float(alpha),
        "variance_C": var_c,
        "variance_ratio_R_C": float(var_c / var_mod),
        "variance_suppression_fraction": float(1.0 - var_c / var_mod),
        "rms_C": float(np.sqrt(np.mean(c ** 2))),
        "max_abs_C": float(np.max(np.abs(c))),
        "mean_C": float(np.mean(c)),
    }


# ---------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------

def write_json(path: Path, data: dict) -> None:
    path.write_text(
        json.dumps(data, indent=2, sort_keys=True, allow_nan=True),
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------
# Raw trajectory
# ---------------------------------------------------------------------

def run_trajectory(n: int) -> dict:
    rho0 = initial_density(n)

    state_rows = []
    payloads = []

    identity_max_abs = None
    previous_exc = None
    excitation_increase_count = 0
    excitation_max_increase = 0.0

    for idx, (s, gamma) in enumerate(zip(S_GRID, GAMMA_GRID)):
        rho_g = amplitude_damp_all(rho0, float(gamma), n)

        if idx == 0:
            identity_max_abs = float(np.max(np.abs(rho_g - rho0)))
            if identity_max_abs > DENSITY_TRACE_TOL:
                raise RuntimeError(
                    f"AMPLITUDE_DAMPING_IDENTITY_FAIL_N{n}:{identity_max_abs:.17e}"
                )

        result = analyze_state(rho_g, n)
        payloads.append(result)

        e_exc = float(result["E_exc"])
        if previous_exc is not None:
            delta_exc = e_exc - previous_exc
            excitation_max_increase = max(excitation_max_increase, delta_exc)
            if delta_exc > DENSITY_TRACE_TOL:
                excitation_increase_count += 1
        previous_exc = e_exc

        row = {
            "N": n,
            "state_index": idx,
            "s": float(s),
            "gamma": float(gamma),
            "entropy_global": float(result["entropy_global"]),
            "E_exc": e_exc,
            "D_eps": float(result["D_eps"]),
            "M_W": float(result["M_W"]),
            "M_mod": float(result["M_mod"]),
            "W_status": result["w_status"],
            "MOD_status": result["v_status"],
            **result["global_diag"],
            **result["modular_diag"],
        }

        for pair_idx, pair in enumerate(PAIRS):
            label = f"{pair[0]}{pair[1]}"
            row[f"W_{label}"] = float(result["w"][pair_idx])
            row[f"W_hat_{label}"] = (
                float(result["w_hat"][pair_idx])
                if result["w_hat"] is not None else np.nan
            )
            row[f"A_{label}"] = float(result["A_coeff"][pair_idx])
            row[f"v_{label}"] = float(result["v"][pair_idx])
            row[f"v_hat_{label}"] = (
                float(result["v_hat"][pair_idx])
                if result["v_hat"] is not None else np.nan
            )

        state_rows.append(row)

        # Release the full global damped matrix before the next trajectory point.
        del rho_g

    if excitation_increase_count > 0:
        raise RuntimeError(
            f"EXCITATION_MONOTONICITY_FAIL_N{n}:"
            f"count={excitation_increase_count},"
            f"max_increase={excitation_max_increase:.17e}"
        )

    recovery = paper28_recovery(payloads[0], n)

    if recovery.get("applicable") and not recovery["passed"]:
        raise RuntimeError(f"PAPER28_S0_RECOVERY_FAIL_N{n}")

    transition_rows = []
    undefined_count = 0

    for k in range(len(GAMMA_GRID) - 1):
        step = transition_step(payloads[k], payloads[k + 1])

        if step["status"] != "VALID":
            undefined_count += 1

        transition_rows.append({
            "N": n,
            "n": k,
            "s_from": float(S_GRID[k]),
            "s_to": float(S_GRID[k + 1]),
            "gamma_from": float(GAMMA_GRID[k]),
            "gamma_to": float(GAMMA_GRID[k + 1]),
            "delta_gamma": float(GAMMA_GRID[k + 1] - GAMMA_GRID[k]),
            "d_W": float(step["d_W"]),
            "d_mod": float(step["d_mod"]),
            "status": step["status"],
        })

    trajectory_summary = {
        "N": n,
        "A_size": A_SIZE,
        "n_states": len(state_rows),
        "n_transitions": len(transition_rows),
        "n_undefined_transitions": int(undefined_count),
        "gamma_min": float(GAMMA_GRID.min()),
        "gamma_max": float(GAMMA_GRID.max()),
        "identity_channel_max_abs_error": identity_max_abs,
        "excitation_increase_count": int(excitation_increase_count),
        "excitation_max_increase": float(excitation_max_increase),
        "paper28_s0_recovery": recovery,
    }

    return {
        "states": state_rows,
        "transitions": transition_rows,
        "summary": trajectory_summary,
    }


# ---------------------------------------------------------------------
# Primary preregistered compensation analysis
# ---------------------------------------------------------------------

def analyze_compensation(n: int, transition_rows: list[dict]) -> dict:
    if len(transition_rows) != len(GAMMA_GRID) - 1:
        raise RuntimeError(f"UNEXPECTED_TRANSITION_COUNT_N{n}")

    statuses = [r["status"] for r in transition_rows]
    if any(status != "VALID" for status in statuses):
        return {
            "N": n,
            "analysis_status": "PRIMARY_INVALID_UNDEFINED_TRANSITIONS",
            "n_undefined": int(sum(status != "VALID" for status in statuses)),
            "primary_supported": False,
            "bandwidth_results": [],
            "residual_rows": [],
        }

    gamma_x = np.asarray(
        [float(r["gamma_from"]) for r in transition_rows],
        dtype=float,
    )
    d_w = np.asarray([float(r["d_W"]) for r in transition_rows], dtype=float)
    d_mod = np.asarray([float(r["d_mod"]) for r in transition_rows], dtype=float)

    # Paper 29 used the range of the predictor values actually entering
    # the transition regression. We retain that exact convention here.
    gamma_range = float(np.max(gamma_x) - np.min(gamma_x))
    if gamma_range <= 0.0:
        raise RuntimeError(f"NONPOSITIVE_GAMMA_RANGE_N{n}")

    bandwidth_results = []
    residual_rows = []

    for frac in BANDWIDTH_FRACTIONS:
        h = float(frac * gamma_range)

        fit_w, weight_w = local_linear_loo_predict(gamma_x, d_w, h)
        fit_m, weight_m = local_linear_loo_predict(gamma_x, d_mod, h)

        resid_w = d_w - fit_w
        resid_m = d_mod - fit_m

        rp = pearson(resid_w, resid_m)
        rs = spearman(resid_w, resid_m)

        z_w, mean_rw, std_rw = zscore(resid_w)
        z_m, mean_rm, std_rm = zscore(resid_m)

        alpha = optimal_alpha(z_w, z_m)
        qinv = quasi_invariant_metrics(z_w, z_m, alpha)
        qinv_alpha1 = quasi_invariant_metrics(z_w, z_m, 1.0)

        endpoint_p1 = bool(rp < 0.0)
        endpoint_p2 = bool(ALPHA_MIN <= alpha <= ALPHA_MAX)
        endpoint_p3 = bool(qinv["variance_ratio_R_C"] <= RC_MAX)
        bw_supported = bool(endpoint_p1 and endpoint_p2 and endpoint_p3)

        bandwidth_results.append({
            "bandwidth_fraction_gamma_range": float(frac),
            "bandwidth_absolute": h,
            "pearson_residuals": rp,
            "spearman_residuals": rs,
            "mean_resid_dW": mean_rw,
            "std_resid_dW": std_rw,
            "mean_resid_dmod": mean_rm,
            "std_resid_dmod": std_rm,
            "min_effective_weight_sum_dW": float(np.min(weight_w)),
            "min_effective_weight_sum_dmod": float(np.min(weight_m)),
            "alpha_optimal": alpha,
            "quasi_invariant_optimal": qinv,
            "alpha_equal_one_control": qinv_alpha1,
            "P1_negative_pearson": endpoint_p1,
            "P2_alpha_in_window": endpoint_p2,
            "P3_R_C_le_0p20": endpoint_p3,
            "bandwidth_supported": bw_supported,
        })

        for i in range(len(gamma_x)):
            residual_rows.append({
                "N": n,
                "bandwidth_fraction_gamma_range": float(frac),
                "bandwidth_absolute": h,
                "n": int(i),
                "gamma": float(gamma_x[i]),
                "d_W": float(d_w[i]),
                "d_mod": float(d_mod[i]),
                "fit_dW": float(fit_w[i]),
                "fit_dmod": float(fit_m[i]),
                "resid_dW": float(resid_w[i]),
                "resid_dmod": float(resid_m[i]),
                "z_resid_dW": float(z_w[i]),
                "z_resid_dmod": float(z_m[i]),
            })

    all_supported = all(r["bandwidth_supported"] for r in bandwidth_results)

    return {
        "N": n,
        "analysis_status": "PREREGISTERED_PRIMARY_ANALYSIS_COMPLETED",
        "gamma_predictor_min": float(np.min(gamma_x)),
        "gamma_predictor_max": float(np.max(gamma_x)),
        "gamma_predictor_range": gamma_range,
        "raw_pearson_dW_dmod": pearson(d_w, d_mod),
        "raw_spearman_dW_dmod": spearman(d_w, d_mod),
        "bandwidth_results": bandwidth_results,
        "primary_supported": bool(all_supported),
        "residual_rows": residual_rows,
    }


# ---------------------------------------------------------------------
# Cross-size transfer — preregistered secondary endpoint
# ---------------------------------------------------------------------

def residual_block_for_bandwidth(analysis: dict, frac: float):
    rows = [
        r for r in analysis["residual_rows"]
        if abs(float(r["bandwidth_fraction_gamma_range"]) - frac) < 1e-12
    ]
    if len(rows) != len(GAMMA_GRID) - 1:
        raise RuntimeError(
            f"UNEXPECTED_RESIDUAL_ROW_COUNT_N{analysis['N']}_BW{frac}"
        )

    z_w = np.asarray([float(r["z_resid_dW"]) for r in rows], dtype=float)
    z_m = np.asarray([float(r["z_resid_dmod"]) for r in rows], dtype=float)
    return z_w, z_m


def alpha_for_bandwidth(analysis: dict, frac: float) -> float:
    matches = [
        r for r in analysis["bandwidth_results"]
        if abs(float(r["bandwidth_fraction_gamma_range"]) - frac) < 1e-12
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"UNEXPECTED_ALPHA_MATCH_COUNT_N{analysis['N']}_BW{frac}"
        )
    return float(matches[0]["alpha_optimal"])


def cross_size_transfer(analyses: dict[int, dict]) -> list[dict]:
    if any(not a["bandwidth_results"] for a in analyses.values()):
        return []

    out = []

    for frac in BANDWIDTH_FRACTIONS:
        alpha12 = alpha_for_bandwidth(analyses[12], frac)
        alpha8 = alpha_for_bandwidth(analyses[8], frac)
        alpha10 = alpha_for_bandwidth(analyses[10], frac)

        zw12, zm12 = residual_block_for_bandwidth(analyses[12], frac)
        zw8, zm8 = residual_block_for_bandwidth(analyses[8], frac)
        zw10, zm10 = residual_block_for_bandwidth(analyses[10], frac)

        out.append({
            "bandwidth_fraction_gamma_range": float(frac),
            "alpha_N12_applied_to_N8": quasi_invariant_metrics(zw8, zm8, alpha12),
            "alpha_N12_applied_to_N10": quasi_invariant_metrics(zw10, zm10, alpha12),
            "alpha_N8_applied_to_N12": quasi_invariant_metrics(zw12, zm12, alpha8),
            "alpha_N10_applied_to_N12": quasi_invariant_metrics(zw12, zm12, alpha10),
        })

    return out


# ---------------------------------------------------------------------
# Frozen protocol manifest
# ---------------------------------------------------------------------

def frozen_manifest() -> dict:
    return {
        "paper": 30,
        "analysis_status": "PROSPECTIVE_PREREGISTERED",
        "primary_system": {
            "N": N_PRIMARY,
            "A": list(A_KEEP),
            "role": "NEW_SIZE_PROSPECTIVE_VALIDATION",
        },
        "replication_systems": [
            {"N": n, "A": list(A_KEEP), "role": "SAME_SIZE_CHANNEL_REPLICATION"}
            for n in N_REPLICATIONS
        ],
        "tfim": {
            "J": J,
            "h": H_FIELD,
            "p": P,
            "dt": DT,
            "nominal_t": P * DT,
            "state": "time-reversal mixture of Strang-evolved |+>^N",
        },
        "channel": {
            "name": "independent local amplitude damping",
            "s_grid": [float(x) for x in S_GRID],
            "gamma_definition": "gamma(s)=1-exp(-s)",
            "gamma_grid": [float(x) for x in GAMMA_GRID],
            "evaluation": "each gamma evaluated directly from same prepared rho0",
        },
        "progress_variable": {
            "primary": "gamma",
            "secondary_regularized_relative_entropy_eps": REL_ENT_EPS,
        },
        "pairs": [list(p) for p in PAIRS],
        "detrending": {
            "method": "leave-one-out local-linear Gaussian-kernel regression",
            "predictor": "gamma_from for the 30 transitions",
            "bandwidth_fractions_of_predictor_range": list(BANDWIDTH_FRACTIONS),
            "weight_floor": WEIGHT_FLOOR,
        },
        "primary_endpoints_N12": {
            "P1": "Pearson(resid_dW,resid_dmod)<0 at all 5 bandwidths",
            "P2": f"{ALPHA_MIN}<=alpha_N12<={ALPHA_MAX} at all 5 bandwidths",
            "P3": f"R_C<={RC_MAX} at all 5 bandwidths",
            "success_rule": "all P1/P2/P3 at all fixed bandwidths",
        },
        "numerics": {
            "state_norm_tol": STATE_NORM_TOL,
            "density_trace_tol": DENSITY_TRACE_TOL,
            "density_positivity_tol": DENSITY_POSITIVITY_TOL,
            "hermiticity_tol": HERMITICITY_TOL,
            "entropy_eig_cutoff": ENTROPY_EIG_CUTOFF,
            "zero_vector_step_eps": STEP_EPS,
            "strict_modular_positivity": True,
            "paper28_recovery_atol": PAPER28_RECOVERY_ATOL,
        },
        "interpretation_boundary": (
            "Support establishes finite-size generalization of the residual "
            "compensation signature to amplitude damping at a new system size. "
            "It does not establish a fundamental conservation law, literal "
            "information transfer, black-hole horizon, thermodynamic-limit "
            "universality, or universal alpha constant."
        ),
    }


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not PREREG_PATH.exists():
        raise RuntimeError(f"MISSING_PREREGISTRATION:{PREREG_PATH}")

    script_path = Path(__file__).resolve()

    provenance = {
        "script_sha256": sha256_file(script_path),
        "preregistration_sha256": sha256_file(PREREG_PATH),
        "script_path": str(script_path),
        "preregistration_path": str(PREREG_PATH),
    }

    write_json(OUT_DIR / "frozen_manifest.json", frozen_manifest())
    write_json(OUT_DIR / "provenance.json", provenance)

    trajectories = {}
    analyses = {}

    # Run all systems without printing scientific results mid-run.
    # Primary N12 first, then direct replication systems.
    for n in (N_PRIMARY, *N_REPLICATIONS):
        trajectory = run_trajectory(n)
        analysis = analyze_compensation(n, trajectory["transitions"])

        trajectories[n] = trajectory
        analyses[n] = analysis

        role = "primary" if n == N_PRIMARY else "replication"

        write_csv(
            OUT_DIR / f"{role}_states_N{n}_A4.csv",
            trajectory["states"],
        )
        write_csv(
            OUT_DIR / f"{role}_transitions_N{n}_A4.csv",
            trajectory["transitions"],
        )
        write_csv(
            OUT_DIR / f"{role}_residuals_N{n}_A4.csv",
            analysis["residual_rows"],
        )

        analysis_to_write = dict(analysis)
        analysis_to_write.pop("residual_rows", None)

        write_json(
            OUT_DIR / f"{role}_summary_N{n}_A4.json",
            {
                "trajectory": trajectory["summary"],
                "compensation": analysis_to_write,
            },
        )

    transfer = cross_size_transfer(analyses)

    primary_analysis = analyses[N_PRIMARY]
    primary_supported = bool(primary_analysis["primary_supported"])

    final_summary = {
        "analysis_status": "PREREGISTERED_PAPER30_COMPLETED",
        "primary_N": N_PRIMARY,
        "primary_classification": (
            "SUPPORTED" if primary_supported else "NOT_SUPPORTED"
        ),
        "primary_success_rule": (
            "N12 must satisfy P1 negative residual Pearson, "
            "P2 alpha in [0.8,1.2], and P3 R_C<=0.20 "
            "at all five frozen bandwidths."
        ),
        "primary_summary_file": "primary_summary_N12_A4.json",
        "replication_summary_files": [
            "replication_summary_N8_A4.json",
            "replication_summary_N10_A4.json",
        ],
        "cross_size_transfer": transfer,
        "regularized_relative_entropy": {
            "epsilon": REL_ENT_EPS,
            "role": "SECONDARY_DIAGNOSTIC_ONLY",
        },
        "interpretation_boundary": frozen_manifest()["interpretation_boundary"],
    }

    write_json(OUT_DIR / "paper30_summary.json", final_summary)

    print("PAPER30_EXECUTION_COMPLETED")
    print(OUT_DIR)


if __name__ == "__main__":
    main()

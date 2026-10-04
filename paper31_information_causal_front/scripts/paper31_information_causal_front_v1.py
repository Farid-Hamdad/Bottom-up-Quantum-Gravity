#!/usr/bin/env python3
"""
Paper 31 — Information Causal Front and Emergent Propagation Speed

Prospectively frozen implementation candidate.

IMPORTANT
---------
- Implements the frozen Paper 31 preregistration.
- Do not modify primary observable, perturbation, time grid, front levels,
  fit rules, tolerances, or classification rules after scientific unblinding.
- Preserve null, negative, failed, and undefined outcomes.
- The Paper 22 comparison remains gated unless an independent unit mapping
  has been established before scientific interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

J = 1.0
H_FIELD = 1.0
DT = 0.35
N = 16
A_KEEP = (0, 1, 2, 3)
A_SIZE = len(A_KEEP)
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
P_GRID = tuple(range(1, 9))
T_GRID = tuple(float(p * DT) for p in P_GRID)
PERTURBATION_SITES = tuple(range(4, N))
DISTANCES = tuple(r - 3 for r in PERTURBATION_SITES)
FRONT_LEVELS = (0.10, 0.25, 0.50)

STATE_NORM_TOL = 1e-12
DENSITY_TRACE_TOL = 1e-12
HERMITICITY_TOL = 1e-12
DENSITY_POSITIVITY_TOL = 1e-12
SIGNAL_ZERO_TOL = 1e-14
CIRCUIT_SUPPORT_TOL = 1e-12
ENTROPY_EIG_CUTOFF = 1e-15
MIN_VALID_FRONT_POINTS = 4
MIN_FRONT_R2 = 0.95
MAX_RELATIVE_VELOCITY_SPREAD = 0.15
C_EDGE_PAPER22 = 1.0
CROSS_SECTOR_REL_TOL = 0.10
UNIT_COMPATIBILITY_GATE = False
UNIT_COMPATIBILITY_NOTE = (
    "No independent common spatial/temporal unit mapping between Paper 31 "
    "lattice-time units and Paper 22 calibrated graph units is frozen in "
    "this implementation. Cross-sector equality claims are therefore gated."
)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"X": X, "Y": Y, "Z": Z}

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "paper31_information_causal_front_v1"
PREREG_PATH = ROOT / "notes" / "PREREGISTRATION_FROZEN.md"


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
    idx = np.arange(2 ** n, dtype=np.uint64)
    energy = np.zeros(2 ** n, dtype=float)
    for q in range(n - 1):
        bq = (idx >> np.uint64(n - 1 - q)) & np.uint64(1)
        br = (idx >> np.uint64(n - 2 - q)) & np.uint64(1)
        zq = 1.0 - 2.0 * bq.astype(float)
        zr = 1.0 - 2.0 * br.astype(float)
        energy += zq * zr
    return state * np.exp(1j * J * dt * energy)


def strang_step(state: np.ndarray, dt: float, n: int) -> np.ndarray:
    state = x_half(state, dt, n)
    state = zz_layer(state, dt, n)
    state = x_half(state, dt, n)
    return state


def validate_state(state: np.ndarray, label: str) -> None:
    norm = float(np.real(np.vdot(state, state)))
    if abs(norm - 1.0) > STATE_NORM_TOL:
        raise RuntimeError(f"STATE_NORM_DRIFT:{label}:{norm:.17e}")


def evolve_forward(initial: np.ndarray, p: int, n: int) -> np.ndarray:
    state = initial.copy()
    for _ in range(p):
        state = strang_step(state, DT, n)
    validate_state(state, f"forward_P{p}")
    return state


def evolution_snapshots(initial: np.ndarray, n: int) -> dict[int, np.ndarray]:
    state = initial.copy()
    out = {}
    for p in P_GRID:
        state = strang_step(state, DT, n)
        validate_state(state, f"snapshot_P{p}")
        out[p] = state.copy()
    return out


def apply_initial_z_perturbation(state: np.ndarray, q: int, n: int) -> np.ndarray:
    return one_gate(state, Z, q, n)


def hermitize(rho: np.ndarray) -> np.ndarray:
    return 0.5 * (rho + rho.conj().T)


def reduced_density_pure(state: np.ndarray, keep: tuple[int, ...], n: int) -> np.ndarray:
    keep = tuple(keep)
    rest = tuple(q for q in range(n) if q not in keep)
    tensor = np.transpose(state.reshape([2] * n), keep + rest)
    dk = 2 ** len(keep)
    dr = 2 ** len(rest)
    psi_matrix = tensor.reshape(dk, dr)
    return hermitize(psi_matrix @ psi_matrix.conj().T)


def reduced_density_mixed(rho: np.ndarray, keep: tuple[int, ...], n: int) -> np.ndarray:
    keep = tuple(keep)
    rest = tuple(q for q in range(n) if q not in keep)
    perm = keep + rest + tuple(n + q for q in keep) + tuple(n + q for q in rest)
    tensor = np.transpose(rho.reshape([2] * (2 * n)), perm)
    dk = 2 ** len(keep)
    dr = 2 ** len(rest)
    out = np.einsum("abcb->ac", tensor.reshape(dk, dr, dk, dr), optimize=True)
    return hermitize(out)


def density_diagnostics(rho: np.ndarray, label: str) -> dict:
    tr = float(np.real(np.trace(rho)))
    herm_err = float(np.max(np.abs(rho - rho.conj().T)))
    eigvals = np.real(np.linalg.eigvalsh(hermitize(rho)))
    lam_min = float(eigvals.min())
    if abs(tr - 1.0) > DENSITY_TRACE_TOL:
        raise RuntimeError(f"TRACE_FAIL:{label}:{tr:.17e}")
    if herm_err > HERMITICITY_TOL:
        raise RuntimeError(f"HERMITICITY_FAIL:{label}:{herm_err:.17e}")
    if lam_min < -DENSITY_POSITIVITY_TOL:
        raise RuntimeError(f"POSITIVITY_FAIL:{label}:{lam_min:.17e}")
    return {
        "trace": tr,
        "hermiticity_error": herm_err,
        "lambda_min": lam_min,
        "rank_above_entropy_cutoff": int(np.sum(eigvals > ENTROPY_EIG_CUTOFF)),
    }


def trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    vals = np.real(np.linalg.eigvalsh(hermitize(rho - sigma)))
    value = 0.5 * float(np.sum(np.abs(vals)))
    if value < -DENSITY_TRACE_TOL or value > 1.0 + DENSITY_TRACE_TOL:
        raise RuntimeError(f"TRACE_DISTANCE_RANGE_FAIL:{value:.17e}")
    return float(np.clip(value, 0.0, 1.0))


def frobenius_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    return float(np.linalg.norm(rho - sigma, ord="fro"))


def entropy(rho: np.ndarray) -> float:
    vals = np.real(np.linalg.eigvalsh(hermitize(rho)))
    vals = np.clip(vals, 0.0, None)
    vals = vals[vals > ENTROPY_EIG_CUTOFF]
    return float(-np.sum(vals * np.log(vals)))


def pair_mutual_information(rho_a: np.ndarray, a_size: int, pairs=PAIRS) -> np.ndarray:
    values = []
    for i, j in pairs:
        rho_i = reduced_density_mixed(rho_a, (i,), a_size)
        rho_j = reduced_density_mixed(rho_a, (j,), a_size)
        rho_ij = reduced_density_mixed(rho_a, (i, j), a_size)
        values.append(entropy(rho_i) + entropy(rho_j) - entropy(rho_ij))
    return np.asarray(values, dtype=float)


def local_ops(a_size: int) -> dict[tuple[int, str], np.ndarray]:
    ops = {}
    for q in range(a_size):
        for label, sigma in PAULI.items():
            op = np.array([[1.0 + 0.0j]])
            for k in range(a_size):
                op = np.kron(op, sigma if k == q else I2)
            ops[(q, label)] = op
    return ops


def comm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def modular_sector_or_undefined(rho_a: np.ndarray) -> dict:
    vals, vecs = np.linalg.eigh(hermitize(rho_a))
    vals = np.real(vals)
    if vals.min() <= ENTROPY_EIG_CUTOFF:
        return {
            "status": "UNDEFINED_RANK_DEFICIENT",
            "lambda_min": float(vals.min()),
            "rank": int(np.sum(vals > ENTROPY_EIG_CUTOFF)),
            "v": None,
        }
    kappa = -np.log(vals)
    kappa_mean = float(kappa.mean())
    kappa_std = float(kappa.std(ddof=0))
    if kappa_std <= 0.0:
        return {
            "status": "UNDEFINED_MODULAR_STD_ZERO",
            "lambda_min": float(vals.min()),
            "rank": int(len(vals)),
            "v": None,
        }
    kappa_norm = (kappa - kappa_mean) / kappa_std
    k_tilde = hermitize(vecs @ np.diag(kappa_norm) @ vecs.conj().T)
    ops = local_ops(A_SIZE)
    coeffs = []
    for i, j in PAIRS:
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
    v = np.sqrt(np.clip(coeffs, 0.0, None))
    return {
        "status": "VALID",
        "lambda_min": float(vals.min()),
        "rank": int(len(vals)),
        "kappa_mean": kappa_mean,
        "kappa_std": kappa_std,
        "v": v,
    }


def subsystem_metrics(rho_a: np.ndarray, label: str) -> dict:
    return {
        "density_diag": density_diagnostics(rho_a, label),
        "w": pair_mutual_information(rho_a, A_SIZE, PAIRS),
        "modular": modular_sector_or_undefined(rho_a),
    }


def outward_crossing(distances: np.ndarray, profile: np.ndarray, q: float) -> float | None:
    crossings = []
    for i in range(len(distances) - 1):
        d0, d1 = float(distances[i]), float(distances[i + 1])
        f0, f1 = float(profile[i]), float(profile[i + 1])
        if not np.isfinite(f0) or not np.isfinite(f1):
            continue
        if f0 == q:
            crossings.append(d0)
        product = (f0 - q) * (f1 - q)
        if product < 0.0:
            frac = (q - f0) / (f1 - f0)
            crossings.append(d0 + frac * (d1 - d0))
        elif f1 == q:
            crossings.append(d1)
    return None if not crossings else float(max(crossings))


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2:
        return {"status": "INSUFFICIENT_POINTS", "n": int(len(x)), "slope": np.nan, "intercept": np.nan, "R2": np.nan}
    design = np.column_stack([x, np.ones_like(x)])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 0.0 and ss_res <= 1e-30 else (np.nan if ss_tot <= 0.0 else float(1.0 - ss_res / ss_tot))
    return {"status": "VALID", "n": int(len(x)), "slope": float(beta[0]), "intercept": float(beta[1]), "R2": r2}


def recovery_forward_branch(reference_snapshots: dict[int, np.ndarray]) -> dict:
    direct = evolve_forward(plus_state(N), 6, N)
    delta = float(np.max(np.abs(direct - reference_snapshots[6])))
    return {"passed": bool(delta <= STATE_NORM_TOL), "max_abs_state_difference": delta, "tolerance": STATE_NORM_TOL}


def theoretical_support_control(raw_rows: list[dict]) -> dict:
    outside = [r for r in raw_rows if int(r["distance"]) > int(r["P"])]
    max_signal = max((float(r["trace_distance"]) for r in outside), default=0.0)
    return {
        "n_outside_support_points": int(len(outside)),
        "max_trace_distance_outside_support": float(max_signal),
        "tolerance": CIRCUIT_SUPPORT_TOL,
        "passed": bool(max_signal <= CIRCUIT_SUPPORT_TOL),
    }


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_raw_signal_grid() -> tuple[list[dict], dict[int, dict]]:
    ref_initial = plus_state(N)
    validate_state(ref_initial, "reference_initial")
    reference_snapshots = evolution_snapshots(ref_initial, N)
    reference_metrics = {}
    for p in P_GRID:
        rho_a = reduced_density_pure(reference_snapshots[p], A_KEEP, N)
        reference_metrics[p] = {"rho_A": rho_a, **subsystem_metrics(rho_a, f"reference_P{p}")}

    rows = []
    for site in PERTURBATION_SITES:
        distance = site - 3
        pert_initial = apply_initial_z_perturbation(ref_initial, site, N)
        validate_state(pert_initial, f"pert_initial_site{site}")
        pert_snapshots = evolution_snapshots(pert_initial, N)
        for p in P_GRID:
            rho_ref = reference_metrics[p]["rho_A"]
            rho_pert = reduced_density_pure(pert_snapshots[p], A_KEEP, N)
            pert_metrics = subsystem_metrics(rho_pert, f"pert_site{site}_P{p}")
            td = trace_distance(rho_pert, rho_ref)
            fd = frobenius_distance(rho_pert, rho_ref)
            w_ref = reference_metrics[p]["w"]
            w_pert = pert_metrics["w"]
            delta_w_vec = w_pert - w_ref
            mod_ref = reference_metrics[p]["modular"]
            mod_pert = pert_metrics["modular"]
            if mod_ref["status"] == "VALID" and mod_pert["status"] == "VALID":
                delta_v_vec = mod_pert["v"] - mod_ref["v"]
                delta_v_norm = float(np.linalg.norm(delta_v_vec))
                mod_status = "VALID"
            else:
                delta_v_vec = np.full(len(PAIRS), np.nan)
                delta_v_norm = np.nan
                mod_status = f"REF_{mod_ref['status']}__PERT_{mod_pert['status']}"
            row = {
                "P": int(p), "t": float(p * DT), "site": int(site), "distance": int(distance),
                "trace_distance": td, "frobenius_distance": fd,
                "delta_W_norm": float(np.linalg.norm(delta_w_vec)), "delta_v_norm": delta_v_norm,
                "modular_status": mod_status,
                "ref_modular_status": mod_ref["status"], "pert_modular_status": mod_pert["status"],
                "ref_rank": int(reference_metrics[p]["density_diag"]["rank_above_entropy_cutoff"]),
                "pert_rank": int(pert_metrics["density_diag"]["rank_above_entropy_cutoff"]),
            }
            for idx, pair in enumerate(PAIRS):
                tag = f"{pair[0]}{pair[1]}"
                row[f"delta_W_{tag}"] = float(delta_w_vec[idx])
                row[f"delta_v_{tag}"] = float(delta_v_vec[idx])
            rows.append(row)
    return rows, reference_metrics


def analyze_front(raw_rows: list[dict]) -> dict:
    global_max = max(float(r["trace_distance"]) for r in raw_rows)
    if global_max <= SIGNAL_ZERO_TOL:
        return {"status": "SIGNAL_UNDEFINED_GLOBAL_MAX_TOO_SMALL", "global_signal_max": float(global_max), "front_rows": [], "front_fits": [], "v_info": np.nan, "sigma_v": np.nan, "internal_coherence": False}
    distances = np.asarray(DISTANCES, dtype=float)
    front_rows = []
    for p in P_GRID:
        time_rows = sorted([r for r in raw_rows if int(r["P"]) == p], key=lambda r: int(r["distance"]))
        profile = np.asarray([float(r["trace_distance"]) / global_max for r in time_rows], dtype=float)
        for row, value in zip(time_rows, profile):
            row["normalized_global_signal"] = float(value)
        for q in FRONT_LEVELS:
            crossing = outward_crossing(distances, profile, q)
            front_rows.append({
                "P": int(p), "t": float(p * DT), "q": float(q),
                "front_position": float(crossing) if crossing is not None else np.nan,
                "status": "VALID" if crossing is not None else "UNDEFINED_NO_CROSSING",
            })
    front_fits = []
    for q in FRONT_LEVELS:
        valid = [r for r in front_rows if float(r["q"]) == q and r["status"] == "VALID"]
        fit = linear_fit(np.asarray([r["t"] for r in valid]), np.asarray([r["front_position"] for r in valid]))
        fit.update({
            "q": float(q),
            "enough_points": bool(fit["n"] >= MIN_VALID_FRONT_POINTS),
            "positive_velocity": bool(np.isfinite(fit["slope"]) and fit["slope"] > 0.0),
            "R2_pass": bool(np.isfinite(fit["R2"]) and fit["R2"] >= MIN_FRONT_R2),
        })
        front_fits.append(fit)
    slopes = np.asarray([float(f["slope"]) for f in front_fits], dtype=float)
    if np.all(np.isfinite(slopes)):
        v_info = float(np.mean(slopes))
        sigma_v = float(np.std(slopes, ddof=0))
        rel_spread = float(sigma_v / v_info) if v_info > 0.0 else np.inf
    else:
        v_info = sigma_v = rel_spread = np.nan
    coherent = bool(
        all(f["enough_points"] and f["positive_velocity"] and f["R2_pass"] for f in front_fits)
        and np.isfinite(rel_spread)
        and rel_spread <= MAX_RELATIVE_VELOCITY_SPREAD
    )
    return {"status": "COMPLETED", "global_signal_max": float(global_max), "front_rows": front_rows, "front_fits": front_fits, "v_info": v_info, "sigma_v": sigma_v, "relative_velocity_spread": rel_spread, "internal_coherence": coherent}


def analyze_circuit_support(raw_rows: list[dict]) -> dict:
    support_rows = []
    for p in P_GRID:
        rows = [r for r in raw_rows if int(r["P"]) == p and float(r["trace_distance"]) > CIRCUIT_SUPPORT_TOL]
        furthest = max((int(r["distance"]) for r in rows), default=None)
        support_rows.append({"P": int(p), "t": float(p * DT), "furthest_distance_above_1e-12": furthest})
    valid = [r for r in support_rows if r["furthest_distance_above_1e-12"] is not None]
    fit = linear_fit(np.asarray([r["t"] for r in valid]), np.asarray([r["furthest_distance_above_1e-12"] for r in valid]))
    return {"support_rows": support_rows, "fit": fit, "v_circuit_estimate": float(fit["slope"]) if fit["status"] == "VALID" else np.nan}


def classify(front: dict, unit_gate: bool) -> dict:
    if not front["internal_coherence"]:
        return {"internal_front_classification": "NOT_INTERNALLY_COHERENT", "cross_sector_classification": "NOT_EVALUABLE", "reason": "Internal front-consistency endpoint failed."}
    if not unit_gate:
        return {"internal_front_classification": "INTERNALLY_COHERENT", "cross_sector_classification": "UNIT_COMPATIBILITY_NOT_ESTABLISHED", "reason": UNIT_COMPATIBILITY_NOTE}
    ratio = float(front["v_info"] / C_EDGE_PAPER22)
    compatible = abs(ratio - 1.0) <= CROSS_SECTOR_REL_TOL
    return {
        "internal_front_classification": "INTERNALLY_COHERENT",
        "cross_sector_classification": "CROSS_SECTOR_COMPATIBLE" if compatible else "CROSS_SECTOR_NOT_COMPATIBLE",
        "v_info_over_c_edge": ratio,
        "relative_difference_from_one": abs(ratio - 1.0),
        "tolerance": CROSS_SECTOR_REL_TOL,
    }


def frozen_manifest() -> dict:
    return {
        "paper": 31,
        "analysis_status": "PROSPECTIVE_PREREGISTERED",
        "system": {"N": N, "A": list(A_KEEP), "boundary": "open nearest-neighbor chain"},
        "dynamics": {"J": J, "h": H_FIELD, "dt": DT, "P_grid": list(P_GRID), "t_grid": list(T_GRID), "initial_state": "|+>^N", "evolution": "forward pure-state Strang branch"},
        "perturbation": {"operator": "Pauli Z", "sites": list(PERTURBATION_SITES), "distances_from_A_boundary": list(DISTANCES)},
        "primary_observable": "trace distance 0.5*||rho_A_pert-rho_A_ref||_1",
        "normalization": "single global maximum over complete preregistered d,t grid",
        "front_levels": list(FRONT_LEVELS),
        "front_fit": "r_q(t)=v_q*t+b_q by OLS",
        "v_info": "mean of v_0.10, v_0.25, v_0.50",
        "internal_endpoint": {"minimum_valid_points_per_q": MIN_VALID_FRONT_POINTS, "minimum_R2": MIN_FRONT_R2, "positive_velocity_required": True, "maximum_relative_velocity_spread": MAX_RELATIVE_VELOCITY_SPREAD},
        "circuit_support_threshold": CIRCUIT_SUPPORT_TOL,
        "unit_compatibility_gate": UNIT_COMPATIBILITY_GATE,
        "unit_compatibility_note": UNIT_COMPATIBILITY_NOTE,
        "paper22": {"c_edge_calibrated": C_EDGE_PAPER22, "comparison_relative_tolerance": CROSS_SECTOR_REL_TOL},
        "numerics": {"state_norm_tol": STATE_NORM_TOL, "density_trace_tol": DENSITY_TRACE_TOL, "hermiticity_tol": HERMITICITY_TOL, "density_positivity_tol": DENSITY_POSITIVITY_TOL, "signal_zero_tol": SIGNAL_ZERO_TOL, "entropy_eig_cutoff": ENTROPY_EIG_CUTOFF},
        "secondary_modular_policy": "No regularization. Rank-deficient rho_A is explicitly marked UNDEFINED_RANK_DEFICIENT.",
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not PREREG_PATH.exists():
        raise RuntimeError(f"MISSING_PREREGISTRATION:{PREREG_PATH}")
    script_path = Path(__file__).resolve()
    write_json(OUT_DIR / "frozen_manifest.json", frozen_manifest())
    write_json(OUT_DIR / "provenance.json", {"script_path": str(script_path), "script_sha256": sha256_file(script_path), "preregistration_path": str(PREREG_PATH), "preregistration_sha256": sha256_file(PREREG_PATH)})
    raw_rows, reference_metrics = build_raw_signal_grid()
    recovery_forward = recovery_forward_branch(evolution_snapshots(plus_state(N), N))
    support_control = theoretical_support_control(raw_rows)
    if not recovery_forward["passed"]:
        raise RuntimeError("FORWARD_BRANCH_RECOVERY_FAIL")
    if not support_control["passed"]:
        raise RuntimeError("OUTSIDE_SUPPORT_NONZERO_SIGNAL_FAIL")
    front = analyze_front(raw_rows)
    circuit = analyze_circuit_support(raw_rows)
    classification = classify(front, UNIT_COMPATIBILITY_GATE)
    write_csv(OUT_DIR / "paper31_signal_grid.csv", raw_rows)
    write_csv(OUT_DIR / "paper31_front_crossings.csv", front["front_rows"])
    write_json(OUT_DIR / "paper31_front_analysis.json", {k: v for k, v in front.items() if k != "front_rows"})
    write_json(OUT_DIR / "paper31_circuit_support.json", circuit)
    reference_summary = []
    for p in P_GRID:
        mod = reference_metrics[p]["modular"]
        reference_summary.append({"P": int(p), "t": float(p * DT), "rank": int(reference_metrics[p]["density_diag"]["rank_above_entropy_cutoff"]), "modular_status": mod["status"], "lambda_min": float(mod["lambda_min"])})
    write_csv(OUT_DIR / "paper31_reference_diagnostics.csv", reference_summary)
    final_summary = {
        "analysis_status": "PREREGISTERED_PAPER31_COMPLETED",
        "primary_system_N": N,
        "A": list(A_KEEP),
        "v_info": front["v_info"],
        "sigma_v": front["sigma_v"],
        "relative_velocity_spread": front.get("relative_velocity_spread", np.nan),
        "front_fits": front["front_fits"],
        "internal_coherence": front["internal_coherence"],
        "v_circuit_estimate": circuit["v_circuit_estimate"],
        "unit_compatibility_gate": UNIT_COMPATIBILITY_GATE,
        "c_edge_paper22": C_EDGE_PAPER22,
        "classification": classification,
        "recovery": {"forward_branch": recovery_forward, "outside_exact_support": support_control},
        "interpretation_boundary": "Paper31 may establish a reproducible finite-speed information front distinct from the formal circuit-support boundary. A direct equality claim with Paper22 c_edge is prohibited unless the independently frozen unit-compatibility gate passes.",
    }
    write_json(OUT_DIR / "paper31_summary.json", final_summary)
    print("PAPER31_EXECUTION_COMPLETED")
    print(OUT_DIR)


if __name__ == "__main__":
    main()

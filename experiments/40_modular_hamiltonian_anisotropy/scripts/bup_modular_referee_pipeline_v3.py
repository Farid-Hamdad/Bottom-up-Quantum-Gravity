#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bup_modular_referee_pipeline_v3.py

Referee-oriented modular pipeline with:
- scaling in N and subsystem size
- anisotropy scan
- comparison of state classes:
    * exact ground state of XXZ/Heisenberg chain
    * thermal Gibbs state exp(-beta H)/Z
    * Haar-random pure state
- automatic plots
- explicit fit quality, locality, and anisotropy amplification diagnostics

Main question
-------------
Does the effective modular Hamiltonian K_A:
1. remain approximately few-body?
2. remain local under local perturbations?
3. inherit or amplify microscopic anisotropy?
4. behave differently for structured states (ground / thermal) vs random states?

Typical usage
-------------
python bup_modular_referee_pipeline_v3.py \
  --n-list 8,10 \
  --subsystem-sizes 3,4,5 \
  --jxy-list 1.0 \
  --jz-list 1.0,1.1,1.25,1.5,2.0 \
  --state-kinds ground,thermal,random \
  --beta 2.0 \
  --periodic \
  --perturb-site 2 \
  --perturb-strength 0.2 \
  --output-dir results_bup_modular_referee_v3

Outputs
-------
results_dir/
  aggregate_rows.csv
  aggregate_summary.json
  plots/
    amplification_vs_eta_micro.png
    locality_vs_eta_micro.png
    fit_r2_vs_eta_micro.png
    amplification_by_state_kind.png
    locality_by_state_kind.png
  N_<N>_A_<m>_state_<kind>_eta_<tag>/
    summary.json
    coefficients_baseline.csv
    coefficients_perturbed.csv
    delta_by_distance.csv
    xxz_diagnostics.json
"""

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

I2 = np.eye(2, dtype=np.complex128)
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)

PAULI_DICT = {"I": I2, "X": X, "Y": Y, "Z": Z}


@dataclass
class FitResult:
    K_fit: np.ndarray
    err_fro: float
    rel_err: float
    r2_matrix: float
    coeff_table: List[dict]
    rank: int
    singular_values: List[float]


def kron_many(ops: List[np.ndarray]) -> np.ndarray:
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def apply_to_site(op: np.ndarray, site: int, n: int) -> np.ndarray:
    ops = [I2] * n
    ops[site] = op
    return kron_many(ops)


def apply_two_sites(op1: np.ndarray, i: int, op2: np.ndarray, j: int, n: int) -> np.ndarray:
    ops = [I2] * n
    ops[i] = op1
    ops[j] = op2
    return kron_many(ops)


def build_xxz_hamiltonian(n: int, periodic: bool, Jxy: float, Jz: float) -> np.ndarray:
    H = np.zeros((2**n, 2**n), dtype=np.complex128)
    edges = [(i, i + 1) for i in range(n - 1)]
    if periodic and n > 2:
        edges.append((n - 1, 0))
    for i, j in edges:
        H += Jxy * apply_two_sites(X, i, X, j, n)
        H += Jxy * apply_two_sites(Y, i, Y, j, n)
        H += Jz  * apply_two_sites(Z, i, Z, j, n)
    H = 0.5 * (H + H.conj().T)
    return H


def build_local_field(n: int, site: int, axis: str = "Z", strength: float = 0.0) -> np.ndarray:
    if strength == 0.0:
        return np.zeros((2**n, 2**n), dtype=np.complex128)
    op = {"X": X, "Y": Y, "Z": Z}[axis.upper()]
    return strength * apply_to_site(op, site, n)


def ground_state(H: np.ndarray) -> Tuple[float, np.ndarray]:
    vals, vecs = np.linalg.eigh(H)
    idx = int(np.argmin(np.real(vals)))
    psi = vecs[:, idx]
    psi = psi / np.linalg.norm(psi)
    return float(np.real(vals[idx])), psi


def thermal_density_matrix(H: np.ndarray, beta: float) -> np.ndarray:
    vals, vecs = np.linalg.eigh(H)
    weights = np.exp(-beta * (np.real(vals) - np.min(np.real(vals))))
    Zpart = np.sum(weights)
    rho = vecs @ np.diag(weights / Zpart) @ vecs.conj().T
    rho = 0.5 * (rho + rho.conj().T)
    return rho


def haar_random_state(dim: int, rng: np.random.Generator) -> np.ndarray:
    x = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    x = x / np.linalg.norm(x)
    return x


def reduced_density_matrix_pure(state: np.ndarray, keep: List[int], n: int) -> np.ndarray:
    keep = sorted(keep)
    rest = [i for i in range(n) if i not in keep]
    psi = state.reshape((2,) * n)
    psi = np.transpose(psi, axes=keep + rest)
    dim_keep = 2 ** len(keep)
    dim_rest = 2 ** len(rest)
    psi = psi.reshape(dim_keep, dim_rest)
    rho = psi @ psi.conj().T
    rho = 0.5 * (rho + rho.conj().T)
    return rho


def reduced_density_matrix_mixed(rho: np.ndarray, keep: List[int], n: int) -> np.ndarray:
    keep = sorted(keep)
    rest = [i for i in range(n) if i not in keep]
    dim_keep = 2 ** len(keep)
    dim_rest = 2 ** len(rest)

    rho_t = rho.reshape((2,) * (2 * n))
    axes = keep + rest + [n + i for i in keep] + [n + i for i in rest]
    rho_t = np.transpose(rho_t, axes=axes)
    rho_t = rho_t.reshape(dim_keep, dim_rest, dim_keep, dim_rest)
    rho_red = np.einsum("a b c b -> a c", rho_t)
    rho_red = 0.5 * (rho_red + rho_red.conj().T)
    return rho_red


def modular_hamiltonian(rho: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    vals, vecs = np.linalg.eigh(rho)
    vals = np.maximum(np.real(vals), eps)
    K = vecs @ np.diag(-np.log(vals)) @ vecs.conj().T
    K = 0.5 * (K + K.conj().T)
    return K


def von_neumann_entropy(rho: np.ndarray, eps: float = 1e-12) -> float:
    vals = np.linalg.eigvalsh(rho)
    vals = np.maximum(np.real(vals), eps)
    vals = vals / np.sum(vals)
    return float(-np.sum(vals * np.log(vals)))


def build_operator_basis(m: int, include_cross: bool = True) -> Tuple[List[np.ndarray], List[str]]:
    basis = []
    labels = []

    basis.append(kron_many([I2] * m))
    labels.append("I")

    for i in range(m):
        for p in ["X", "Y", "Z"]:
            ops = [I2] * m
            ops[i] = PAULI_DICT[p]
            basis.append(kron_many(ops))
            labels.append(f"{p}_{i}")

    pair_ops = [("X", "X"), ("Y", "Y"), ("Z", "Z")]
    if include_cross:
        pair_ops += [
            ("X", "Y"), ("X", "Z"),
            ("Y", "X"), ("Y", "Z"),
            ("Z", "X"), ("Z", "Y"),
        ]

    for i, j in combinations(range(m), 2):
        for p, q in pair_ops:
            ops = [I2] * m
            ops[i] = PAULI_DICT[p]
            ops[j] = PAULI_DICT[q]
            basis.append(kron_many(ops))
            labels.append(f"{p}_{i}{q}_{j}")

    return basis, labels


def fit_effective_hamiltonian(K: np.ndarray, basis: List[np.ndarray], labels: List[str]) -> FitResult:
    d = K.shape[0]
    M = len(basis)
    A = np.zeros((d * d, M), dtype=np.complex128)
    b = K.reshape(-1)

    for mu, O in enumerate(basis):
        A[:, mu] = O.reshape(-1)

    coeffs, residuals, rank, svals = np.linalg.lstsq(A, b, rcond=None)

    K_fit = np.zeros_like(K, dtype=np.complex128)
    for c, O in zip(coeffs, basis):
        K_fit += c * O
    K_fit = 0.5 * (K_fit + K_fit.conj().T)

    err_fro = float(np.linalg.norm(K - K_fit, ord="fro"))
    norm_K = float(np.linalg.norm(K, ord="fro"))
    rel_err = err_fro / norm_K if norm_K > 0 else float("nan")
    r2_matrix = 1.0 - (err_fro**2 / norm_K**2) if norm_K > 0 else float("nan")

    coeff_table = []
    for label, c in zip(labels, coeffs):
        c_real = float(np.real_if_close(c))
        coeff_table.append({"term": label, "coeff": c_real, "abs_coeff": abs(c_real)})
    coeff_table = sorted(coeff_table, key=lambda x: x["abs_coeff"], reverse=True)

    return FitResult(
        K_fit=K_fit,
        err_fro=err_fro,
        rel_err=rel_err,
        r2_matrix=float(r2_matrix),
        coeff_table=coeff_table,
        rank=int(rank),
        singular_values=[float(np.real_if_close(x)) for x in svals],
    )


def coeff_map(coeff_table: List[dict]) -> Dict[str, float]:
    return {row["term"]: float(row["coeff"]) for row in coeff_table}


def write_coeff_table(path: Path, coeff_table: List[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["term", "coeff", "abs_coeff"])
        w.writeheader()
        for row in coeff_table:
            w.writerow(row)


def extract_pair_type(term: str):
    m = re.fullmatch(r'([XYZ])_(\d+)([XYZ])_(\d+)', str(term))
    if not m:
        return None
    p, i, q, j = m.groups()
    return p, int(i), q, int(j)


def pair_distance_from_site(term: str, site: int) -> float:
    if term == "I":
        return math.inf
    matches = re.findall(r'([XYZ])_(\d+)', str(term))
    if not matches:
        return math.inf
    qubits = [int(idx) for _, idx in matches]
    return float(min(abs(q - site) for q in qubits))


def locality_from_coeffs(coeffs_base: Dict[str, float], coeffs_pert: Dict[str, float], site: int, noise_floor: float = 1e-12):
    keys = sorted(set(coeffs_base) | set(coeffs_pert))
    rows = []
    for k in keys:
        d = abs(coeffs_pert.get(k, 0.0) - coeffs_base.get(k, 0.0))
        r = pair_distance_from_site(k, site)
        if math.isfinite(r):
            rows.append({"term": k, "distance": r, "delta": d})
    df = pd.DataFrame(rows)
    grouped = (
        df.groupby("distance", as_index=False)
          .agg(mean_delta=("delta", "mean"), max_delta=("delta", "max"), n_terms=("delta", "size"))
          .sort_values("distance")
    )

    delta_near = float(grouped.loc[grouped["distance"] == 0.0, "mean_delta"].iloc[0]) if np.any(grouped["distance"] == 0.0) else 0.0
    far = grouped[grouped["distance"] > 0.0]
    delta_far = float(far["mean_delta"].mean()) if not far.empty else 0.0
    compact_locality = bool(np.all(far["mean_delta"].to_numpy(dtype=float) < noise_floor)) if not far.empty else True
    xi_upper_bound = float(far["distance"].min()) if compact_locality and not far.empty else float("nan")
    locality_score = delta_near - delta_far
    return grouped, {
        "delta_near": delta_near,
        "delta_far": delta_far,
        "locality_score": locality_score,
        "compact_locality": compact_locality,
        "xi_upper_bound": xi_upper_bound,
    }


def xxz_diagnostics(coeff_table: List[dict], eta_micro: float):
    coeffs = coeff_map(coeff_table)
    xx, yy, zz, cross = [], [], [], []

    for term, c in coeffs.items():
        parsed = extract_pair_type(term)
        if parsed is None:
            continue
        p, i, q, j = parsed
        ac = abs(c)
        if p == "X" and q == "X":
            xx.append(ac)
        elif p == "Y" and q == "Y":
            yy.append(ac)
        elif p == "Z" and q == "Z":
            zz.append(ac)
        else:
            cross.append(ac)

    mean_xx = float(np.mean(xx)) if xx else 0.0
    mean_yy = float(np.mean(yy)) if yy else 0.0
    mean_zz = float(np.mean(zz)) if zz else 0.0
    mean_cross = float(np.mean(cross)) if cross else 0.0

    mean_xy = 0.5 * (mean_xx + mean_yy)
    eta_eff = float(mean_zz / mean_xy) if mean_xy > 0 else float("nan")
    amplification = float(eta_eff / eta_micro) if eta_micro > 0 and np.isfinite(eta_eff) else float("nan")

    xxyy_gap = abs(mean_xx - mean_yy)
    xxz_anisotropy = abs(mean_xy - mean_zz)

    regime = "undetermined"
    scale = max(mean_xx, mean_yy, mean_zz, mean_cross, 1e-12)
    if xxyy_gap < 0.1 * scale and mean_cross < 0.5 * scale:
        if xxz_anisotropy < 0.1 * scale:
            regime = "heisenberg_like"
        else:
            regime = "xxz_like"
    else:
        regime = "generic_anisotropic"

    return {
        "mean_abs_xx": mean_xx,
        "mean_abs_yy": mean_yy,
        "mean_abs_zz": mean_zz,
        "mean_abs_cross": mean_cross,
        "xxyy_gap": float(xxyy_gap),
        "xxz_anisotropy": float(xxz_anisotropy),
        "eta_micro": float(eta_micro),
        "eta_eff": eta_eff,
        "amplification": amplification,
        "regime": regime,
    }


def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_str_list(s: str) -> List[str]:
    return [str(x.strip()) for x in s.split(",") if x.strip()]


def eta_tag(x: float) -> str:
    return str(x).replace("-", "m").replace(".", "p")


def build_state_object(state_kind: str, H: np.ndarray, n: int, beta: float, rng: np.random.Generator):
    if state_kind == "ground":
        E0, psi = ground_state(H)
        return {"kind": "pure", "energy": E0, "state": psi}
    elif state_kind == "thermal":
        rho = thermal_density_matrix(H, beta=beta)
        # energy expectation
        E = float(np.real(np.trace(rho @ H)))
        return {"kind": "mixed", "energy": E, "state": rho}
    elif state_kind == "random":
        psi = haar_random_state(2**n, rng)
        E = float(np.real(np.vdot(psi, H @ psi)))
        return {"kind": "pure", "energy": E, "state": psi}
    else:
        raise ValueError(f"Unknown state_kind: {state_kind}")


def reduced_density_from_stateobj(obj, keep: List[int], n: int):
    if obj["kind"] == "pure":
        return reduced_density_matrix_pure(obj["state"], keep, n)
    elif obj["kind"] == "mixed":
        return reduced_density_matrix_mixed(obj["state"], keep, n)
    else:
        raise ValueError("Unknown state object kind")


def make_line_plot(df: pd.DataFrame, xcol: str, ycol: str, hue: str, outpath: Path, title: str, xlabel: str, ylabel: str):
    if df.empty:
        return
    plt.figure(figsize=(7, 5))
    for key, sub in df.groupby(hue):
        sub = sub.sort_values(xcol)
        plt.plot(sub[xcol], sub[ycol], marker='o', label=str(key))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def make_bar_plot(df: pd.DataFrame, category: str, value: str, outpath: Path, title: str, ylabel: str):
    if df.empty:
        return
    agg = df.groupby(category, as_index=False)[value].mean()
    plt.figure(figsize=(7, 5))
    plt.bar(agg[category], agg[value])
    plt.title(title)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Referee-grade modular pipeline with anisotropy scans and state-class comparisons")
    parser.add_argument("--n-list", type=str, default="8,10")
    parser.add_argument("--subsystem-sizes", type=str, default="3,4,5")
    parser.add_argument("--jxy-list", type=str, default="1.0")
    parser.add_argument("--jz-list", type=str, default="1.0,1.1,1.25,1.5,2.0")
    parser.add_argument("--state-kinds", type=str, default="ground,thermal,random")
    parser.add_argument("--beta", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--periodic", action="store_true")
    parser.add_argument("--include-cross", action="store_true", default=True)
    parser.add_argument("--perturb-site", type=int, default=2)
    parser.add_argument("--perturb-strength", type=float, default=0.2)
    parser.add_argument("--perturb-axis", type=str, default="Z", choices=["X", "Y", "Z"])
    parser.add_argument("--output-dir", type=str, default="results_bup_modular_referee_v3")
    args = parser.parse_args()

    n_list = parse_int_list(args.n_list)
    subsystem_sizes = parse_int_list(args.subsystem_sizes)
    jxy_list = parse_float_list(args.jxy_list)
    jz_list = parse_float_list(args.jz_list)
    state_kinds = parse_str_list(args.state_kinds)

    outroot = Path(args.output_dir)
    outroot.mkdir(parents=True, exist_ok=True)
    plots_dir = outroot / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    aggregate_rows = []

    print("=" * 92)
    print("BUP MODULAR REFEREE PIPELINE — v3")
    print("=" * 92)
    print(f"N list            = {n_list}")
    print(f"subsystem sizes   = {subsystem_sizes}")
    print(f"Jxy list          = {jxy_list}")
    print(f"Jz list           = {jz_list}")
    print(f"state kinds       = {state_kinds}")
    print(f"beta              = {args.beta}")
    print(f"seed              = {args.seed}")
    print(f"periodic          = {args.periodic}")
    print(f"perturb site      = {args.perturb_site}")
    print(f"perturb strength  = {args.perturb_strength}")
    print(f"perturb axis      = {args.perturb_axis}")
    print(f"include_cross     = {args.include_cross}")
    print(f"output dir        = {outroot.resolve()}")
    print()

    for N in n_list:
        for Jxy in jxy_list:
            for Jz in jz_list:
                eta_micro = float(Jz / Jxy) if Jxy != 0 else float("nan")
                H_base = build_xxz_hamiltonian(N, periodic=args.periodic, Jxy=Jxy, Jz=Jz)
                H_pert = H_base + build_local_field(N, args.perturb_site, axis=args.perturb_axis, strength=args.perturb_strength)

                for state_kind in state_kinds:
                    base_obj = build_state_object(state_kind, H_base, N, args.beta, rng)
                    pert_obj = build_state_object(state_kind, H_pert, N, args.beta, rng)

                    for m in subsystem_sizes:
                        if m >= N:
                            continue

                        A = list(range(m))
                        case_dir = outroot / f"N_{N}_A_{m}_state_{state_kind}_eta_{eta_tag(eta_micro)}"
                        case_dir.mkdir(parents=True, exist_ok=True)

                        rho_A_base = reduced_density_from_stateobj(base_obj, A, N)
                        rho_A_pert = reduced_density_from_stateobj(pert_obj, A, N)
                        K_A_base = modular_hamiltonian(rho_A_base)
                        K_A_pert = modular_hamiltonian(rho_A_pert)

                        S_A_base = von_neumann_entropy(rho_A_base)
                        S_A_pert = von_neumann_entropy(rho_A_pert)

                        basis, labels = build_operator_basis(m, include_cross=args.include_cross)
                        fit_base = fit_effective_hamiltonian(K_A_base, basis, labels)
                        fit_pert = fit_effective_hamiltonian(K_A_pert, basis, labels)

                        write_coeff_table(case_dir / "coefficients_baseline.csv", fit_base.coeff_table)
                        write_coeff_table(case_dir / "coefficients_perturbed.csv", fit_pert.coeff_table)

                        coeffs_base = coeff_map(fit_base.coeff_table)
                        coeffs_pert = coeff_map(fit_pert.coeff_table)

                        delta_by_distance, locality = locality_from_coeffs(coeffs_base, coeffs_pert, site=args.perturb_site)
                        delta_by_distance.to_csv(case_dir / "delta_by_distance.csv", index=False)

                        xxz_base = xxz_diagnostics(fit_base.coeff_table, eta_micro=eta_micro)
                        xxz_pert = xxz_diagnostics(fit_pert.coeff_table, eta_micro=eta_micro)

                        with open(case_dir / "xxz_diagnostics.json", "w", encoding="utf-8") as f:
                            json.dump({"baseline": xxz_base, "perturbed": xxz_pert}, f, indent=2)

                        summary = {
                            "N": N,
                            "subsystem": A,
                            "subsystem_size": m,
                            "periodic": args.periodic,
                            "state_kind": state_kind,
                            "Jxy": Jxy,
                            "Jz": Jz,
                            "eta_micro": eta_micro,
                            "perturb_site": args.perturb_site,
                            "perturb_strength": args.perturb_strength,
                            "perturb_axis": args.perturb_axis,
                            "energy_baseline": base_obj["energy"],
                            "energy_perturbed": pert_obj["energy"],
                            "entropy_S_A_baseline": S_A_base,
                            "entropy_S_A_perturbed": S_A_pert,
                            "fit_baseline": {
                                "rank": fit_base.rank,
                                "fro_error": fit_base.err_fro,
                                "relative_error": fit_base.rel_err,
                                "r2_matrix": fit_base.r2_matrix,
                                "top_terms": fit_base.coeff_table[:20],
                            },
                            "fit_perturbed": {
                                "rank": fit_pert.rank,
                                "fro_error": fit_pert.err_fro,
                                "relative_error": fit_pert.rel_err,
                                "r2_matrix": fit_pert.r2_matrix,
                                "top_terms": fit_pert.coeff_table[:20],
                            },
                            "xxz_diagnostics_baseline": xxz_base,
                            "xxz_diagnostics_perturbed": xxz_pert,
                            "locality": locality,
                        }
                        with open(case_dir / "summary.json", "w", encoding="utf-8") as f:
                            json.dump(summary, f, indent=2)

                        row = {
                            "N": N,
                            "A_size": m,
                            "state_kind": state_kind,
                            "Jxy": Jxy,
                            "Jz": Jz,
                            "eta_micro": eta_micro,
                            "energy_baseline": base_obj["energy"],
                            "energy_perturbed": pert_obj["energy"],
                            "S_A_baseline": S_A_base,
                            "S_A_perturbed": S_A_pert,
                            "fit_rel_err_baseline": fit_base.rel_err,
                            "fit_r2_baseline": fit_base.r2_matrix,
                            "fit_rel_err_perturbed": fit_pert.rel_err,
                            "fit_r2_perturbed": fit_pert.r2_matrix,
                            "xxz_regime_baseline": xxz_base["regime"],
                            "xxz_regime_perturbed": xxz_pert["regime"],
                            "mean_abs_xx_baseline": xxz_base["mean_abs_xx"],
                            "mean_abs_yy_baseline": xxz_base["mean_abs_yy"],
                            "mean_abs_zz_baseline": xxz_base["mean_abs_zz"],
                            "mean_abs_cross_baseline": xxz_base["mean_abs_cross"],
                            "eta_eff_baseline": xxz_base["eta_eff"],
                            "amplification_baseline": xxz_base["amplification"],
                            "delta_near": locality["delta_near"],
                            "delta_far": locality["delta_far"],
                            "locality_score": locality["locality_score"],
                            "compact_locality": locality["compact_locality"],
                            "xi_upper_bound": locality["xi_upper_bound"],
                        }
                        aggregate_rows.append(row)

                        print(
                            f"[OK] state={state_kind:7s} | N={N} | |A|={m} | eta_micro={eta_micro:.3f} | "
                            f"R²={fit_base.r2_matrix:.5f} | regime={xxz_base['regime']} | "
                            f"eta_eff={xxz_base['eta_eff']:.3f} | amp={xxz_base['amplification']:.3f}"
                        )

    agg_csv = outroot / "aggregate_rows.csv"
    if aggregate_rows:
        with open(agg_csv, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(aggregate_rows[0].keys()))
            w.writeheader()
            for row in aggregate_rows:
                w.writerow(row)

    df = pd.DataFrame(aggregate_rows)

    aggregate_summary = {
        "n_cases": int(len(df)),
        "n_list": n_list,
        "subsystem_sizes": subsystem_sizes,
        "jxy_list": jxy_list,
        "jz_list": jz_list,
        "state_kinds": state_kinds,
        "beta": args.beta,
        "periodic": args.periodic,
        "perturb_site": args.perturb_site,
        "perturb_strength": args.perturb_strength,
        "mean_fit_r2_baseline": float(df["fit_r2_baseline"].mean()) if not df.empty else float("nan"),
        "mean_fit_r2_perturbed": float(df["fit_r2_perturbed"].mean()) if not df.empty else float("nan"),
        "mean_locality_score": float(df["locality_score"].mean()) if not df.empty else float("nan"),
        "fraction_compact_locality": float(df["compact_locality"].mean()) if not df.empty else float("nan"),
        "regimes_baseline": sorted(list(df["xxz_regime_baseline"].dropna().unique())) if not df.empty else [],
        "mean_amplification_baseline": float(df["amplification_baseline"].replace([np.inf, -np.inf], np.nan).dropna().mean()) if not df.empty else float("nan"),
    }
    with open(outroot / "aggregate_summary.json", "w", encoding="utf-8") as f:
        json.dump(aggregate_summary, f, indent=2)

    # Plots
    if not df.empty:
        make_line_plot(
            df.groupby(["state_kind", "eta_micro"], as_index=False)["amplification_baseline"].mean(),
            xcol="eta_micro", ycol="amplification_baseline", hue="state_kind",
            outpath=plots_dir / "amplification_vs_eta_micro.png",
            title="Amplification vs microscopic anisotropy",
            xlabel="eta_micro = Jz / Jxy", ylabel="eta_eff / eta_micro"
        )
        make_line_plot(
            df.groupby(["state_kind", "eta_micro"], as_index=False)["locality_score"].mean(),
            xcol="eta_micro", ycol="locality_score", hue="state_kind",
            outpath=plots_dir / "locality_vs_eta_micro.png",
            title="Locality score vs microscopic anisotropy",
            xlabel="eta_micro = Jz / Jxy", ylabel="mean locality score"
        )
        make_line_plot(
            df.groupby(["state_kind", "eta_micro"], as_index=False)["fit_r2_baseline"].mean(),
            xcol="eta_micro", ycol="fit_r2_baseline", hue="state_kind",
            outpath=plots_dir / "fit_r2_vs_eta_micro.png",
            title="Fit R² vs microscopic anisotropy",
            xlabel="eta_micro = Jz / Jxy", ylabel="mean matrix R²"
        )
        make_bar_plot(
            df, category="state_kind", value="amplification_baseline",
            outpath=plots_dir / "amplification_by_state_kind.png",
            title="Mean amplification by state class",
            ylabel="mean amplification"
        )
        make_bar_plot(
            df, category="state_kind", value="locality_score",
            outpath=plots_dir / "locality_by_state_kind.png",
            title="Mean locality score by state class",
            ylabel="mean locality score"
        )

    print()
    print("Files written")
    print("-" * 92)
    print(outroot / "aggregate_summary.json")
    print(outroot / "aggregate_rows.csv")
    print(plots_dir)
    print("=" * 92)


if __name__ == "__main__":
    main()

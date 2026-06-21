#!/usr/bin/env python3
"""Classical-shadow style pairwise MI tests for the Zoller/Joshi raw data.

This script compares the modular profile 1/beta_j with pairwise mutual
information densities reconstructed from raw Pauli measurements.

It reports three densities for each analysed subsystem A:
  - all_chain: sum_k I(j:k)
  - inside_A:  sum_{k in A} I(j:k)
  - outside_A: sum_{k not in A} I(j:k)

The outside_A channel is important because beta_j is a modular/entanglement
Hamiltonian profile for a subsystem, so its boundary structure is not expected
to be identical to the purely internal pairwise MI.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from analyze_zoller_bup import load_mat_v5
from raw_pairwise_mi_test import (
    AXES,
    PAULI,
    decode_measurements,
    load_beta_profiles,
    pearson,
    single_rho,
    spearman,
)


def project_density(rho: np.ndarray) -> np.ndarray:
    """Project a Hermitian trace-one matrix to the PSD trace-one cone."""
    herm = (rho + rho.conj().T) / 2.0
    vals, vecs = np.linalg.eigh(herm)
    vals = np.clip(np.real(vals), 0.0, None)
    if vals.sum() <= 0:
        vals[:] = 1.0 / len(vals)
    else:
        vals /= vals.sum()
    return (vecs @ np.diag(vals) @ vecs.conj().T)


def entropy(rho: np.ndarray, physical: bool = True) -> float:
    if physical:
        rho = project_density(rho)
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    vals = np.clip(np.real(vals), 1e-12, None)
    vals /= vals.sum()
    return float(-np.sum(vals * np.log(vals)))


def reconstruct_raw_observables(raw_path: Path, estimator: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return single Bloch vectors, pair correlations and coverage counts.

    estimator:
      - conditional: coefficient is the conditional mean for the requested
        Pauli bases. This is exact when that basis combination exists.
      - shadow_uniform: coefficient is the global classical-shadow estimator
        with factors 3 and 9, assuming uniform Pauli bases. Missing basis
        combinations contribute zero.
    """
    data = load_mat_v5(raw_path)
    axes, outcomes = decode_measurements(data["combdatasorted"], data["combtomohistory"].astype(np.uint8))
    n_rows, n_shots, n_sites = outcomes.shape
    total = n_rows * n_shots

    r = np.zeros((n_sites, 3), dtype=float)
    corr = np.full((n_sites, n_sites, 3, 3), np.nan, dtype=float)
    coverage = np.zeros((n_sites, n_sites, 3, 3), dtype=int)

    for i in range(n_sites):
        for a in range(3):
            mask = axes[:, i] == a
            if estimator == "conditional":
                r[i, a] = outcomes[mask, :, i].mean()
            elif estimator == "shadow_uniform":
                r[i, a] = 3.0 * outcomes[mask, :, i].sum() / total
            else:
                raise ValueError(estimator)

    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            for a in range(3):
                for b in range(3):
                    mask = (axes[:, i] == a) & (axes[:, j] == b)
                    count = int(mask.sum())
                    coverage[i, j, a, b] = coverage[j, i, b, a] = count
                    if count == 0:
                        if estimator == "shadow_uniform":
                            corr[i, j, a, b] = 0.0
                            corr[j, i, b, a] = 0.0
                        continue
                    vals = outcomes[mask, :, i] * outcomes[mask, :, j]
                    if estimator == "conditional":
                        value = vals.mean()
                    else:
                        value = 9.0 * vals.sum() / total
                    corr[i, j, a, b] = value
                    corr[j, i, b, a] = value
    return r, corr, coverage


def pair_rho_from_coeffs(ri: np.ndarray, rj: np.ndarray, corr: np.ndarray) -> np.ndarray | None:
    if np.isnan(corr).any():
        return None
    rho = np.kron(PAULI["I"], PAULI["I"]).astype(complex)
    for a, label in enumerate(AXES):
        rho += ri[a] * np.kron(PAULI[label], PAULI["I"])
        rho += rj[a] * np.kron(PAULI["I"], PAULI[label])
    for a, label_a in enumerate(AXES):
        for b, label_b in enumerate(AXES):
            rho += corr[a, b] * np.kron(PAULI[label_a], PAULI[label_b])
    return rho / 4.0


def mi_matrix(raw_path: Path, estimator: str, physical: bool = True) -> tuple[np.ndarray, np.ndarray]:
    r, corr, coverage = reconstruct_raw_observables(raw_path, estimator)
    n_sites = r.shape[0]
    singles = np.array([entropy(single_rho(r[i]), physical=physical) for i in range(n_sites)])
    mi = np.full((n_sites, n_sites), np.nan, dtype=float)
    np.fill_diagonal(mi, 0.0)

    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            rho_ij = pair_rho_from_coeffs(r[i], r[j], corr[i, j])
            if rho_ij is None:
                continue
            sij = entropy(rho_ij, physical=physical)
            val = max(0.0, singles[i] + singles[j] - sij)
            mi[i, j] = mi[j, i] = val
    return mi, coverage


def normalized(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    total = np.nansum(v)
    return v / total if total > 0 else v


def edge_to_center(v: np.ndarray) -> float:
    v = np.asarray(v, dtype=float)
    return float(((v[0] + v[-1]) / 2.0) / v[len(v) // 2])


def compare_density(state: str, estimator: str, density_name: str, density: np.ndarray, profile: dict) -> dict:
    sites0 = profile["sites"] - 1
    y = normalized(density[sites0])
    x = profile["inv_beta_norm"]
    return {
        "state": state,
        "estimator": estimator,
        "density": density_name,
        "cell_index": profile["cell_index"],
        "n_sites": int(len(sites0)),
        "kind": profile["kind"],
        "pearson": pearson(y, x),
        "spearman": spearman(y, x),
        "mae_norm": float(np.nanmean(np.abs(y - x))),
        "density_edge_to_center": edge_to_center(y),
        "invbeta_edge_to_center": edge_to_center(x),
    }


def run_one(state: str, raw_path: Path, beta_path: Path, estimator: str) -> tuple[list[dict], dict]:
    mi, coverage = mi_matrix(raw_path, estimator=estimator, physical=True)
    profiles = load_beta_profiles(beta_path, group="databulk")

    rows = []
    all_density = np.nansum(mi, axis=1)
    for prof in profiles:
        sites0 = prof["sites"] - 1
        in_mask = np.zeros(mi.shape[0], dtype=bool)
        in_mask[sites0] = True

        inside = np.zeros(mi.shape[0], dtype=float)
        outside = np.zeros(mi.shape[0], dtype=float)
        for site in sites0:
            inside[site] = np.nansum(mi[site, in_mask])
            outside[site] = np.nansum(mi[site, ~in_mask])

        rows.append(compare_density(state, estimator, "all_chain_sumI", all_density, prof))
        rows.append(compare_density(state, estimator, "inside_A_sumI", inside, prof))
        rows.append(compare_density(state, estimator, "outside_A_sumI", outside, prof))

    exp = [r for r in rows if r["kind"] == "experimental"]
    summary = {}
    for density_name in sorted(set(r["density"] for r in exp)):
        sub = [r for r in exp if r["density"] == density_name]
        summary[density_name] = {
            "n": len(sub),
            "mean_pearson": float(np.nanmean([r["pearson"] for r in sub])),
            "mean_spearman": float(np.nanmean([r["spearman"] for r in sub])),
            "mean_mae_norm": float(np.nanmean([r["mae_norm"] for r in sub])),
            "mean_density_edge_to_center": float(np.nanmean([r["density_edge_to_center"] for r in sub])),
            "mean_invbeta_edge_to_center": float(np.nanmean([r["invbeta_edge_to_center"] for r in sub])),
        }
    return rows, summary


def main() -> None:
    outdir = Path("/workspace/bup_zoller_tomography/results_shadow_mi")
    outdir.mkdir(parents=True, exist_ok=True)

    jobs = [
        ("GS_Delta_1", Path("/workspace/.cache/03-dataRawParam_GS_Delta_1.mat"), Path("/workspace/.cache/03-GroundStateBetas.mat")),
        ("ES_Delta_1", Path("/workspace/.cache/01-dataRawParam_ES_Delta_1.mat"), Path("/workspace/.cache/01-ExcitedStateBetas.mat")),
    ]
    estimators = ["conditional", "shadow_uniform"]

    all_rows = []
    summary = {}
    for estimator in estimators:
        summary[estimator] = {}
        for state, raw_path, beta_path in jobs:
            rows, state_summary = run_one(state, raw_path, beta_path, estimator)
            all_rows.extend(rows)
            summary[estimator][state] = state_summary

    csv_path = outdir / "bup_shadow_MI_vs_inv_beta_tests.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    summary_path = outdir / "bup_shadow_MI_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Pairwise mutual-information test from Zoller/Joshi raw bitstring data.

This reconstructs every two-site reduced density matrix by Pauli linear
inversion from the 243 measurement bases, then tests the BuP relation

    sum_k I(j:k)  ~  1 / beta_j

on the sites covered by the analysed beta profiles.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from analyze_zoller_bup import load_mat_v5


PAULI = {
    "I": np.array([[1, 0], [0, 1]], dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
AXIS_INDEX = {"X": 0, "Y": 1, "Z": 2}
AXES = ["X", "Y", "Z"]


def basis_axis_sign(code: int) -> tuple[int, int]:
    # 1,2,3,4,5,6 denote X,-X,Y,-Y,Z,-Z.
    if code == 1:
        return 0, +1
    if code == 2:
        return 0, -1
    if code == 3:
        return 1, +1
    if code == 4:
        return 1, -1
    if code == 5:
        return 2, +1
    if code == 6:
        return 2, -1
    raise ValueError(f"unknown basis code: {code}")


def entropy_from_eigs(eigs: np.ndarray) -> float:
    eigs = np.real(eigs)
    eigs = np.clip(eigs, 1e-12, None)
    eigs = eigs / eigs.sum()
    return float(-np.sum(eigs * np.log(eigs)))


def entropy_density(rho: np.ndarray) -> float:
    eigs = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    return entropy_from_eigs(eigs)


def decode_measurements(raw_data: np.ndarray, history: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return axes/signs and canonical outcomes.

    outcomes[row, shot, site] is the eigenvalue of the canonical Pauli axis
    X/Y/Z selected by history[row, site].
    """
    bitstrings = raw_data[:, 1:].astype(np.uint64)
    n_rows, n_shots = bitstrings.shape
    n_sites = history.shape[1]

    axes = np.zeros((n_rows, n_sites), dtype=np.int8)
    signs = np.ones((n_rows, n_sites), dtype=np.int8)
    for code in range(1, 7):
        mask = history == code
        if np.any(mask):
            axis, sign = basis_axis_sign(code)
            axes[mask] = axis
            signs[mask] = sign

    outcomes = np.empty((n_rows, n_shots, n_sites), dtype=np.int8)
    for site in range(n_sites):
        bits = ((bitstrings >> np.uint64(site)) & np.uint64(1)).astype(np.int8)
        # bit=0 -> + along measured basis; bit=1 -> - along measured basis.
        outcomes[:, :, site] = signs[:, site][:, None] * (1 - 2 * bits)

    return axes, outcomes


def estimate_single(axes: np.ndarray, outcomes: np.ndarray) -> np.ndarray:
    n_sites = axes.shape[1]
    r = np.zeros((n_sites, 3), dtype=float)
    for site in range(n_sites):
        for a in range(3):
            mask = axes[:, site] == a
            r[site, a] = outcomes[mask, :, site].mean()
    return r


def estimate_pair_corr(axes: np.ndarray, outcomes: np.ndarray, i: int, j: int) -> np.ndarray | None:
    corr = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            mask = (axes[:, i] == a) & (axes[:, j] == b)
            if not np.any(mask):
                return None
            vals = outcomes[mask, :, i] * outcomes[mask, :, j]
            corr[a, b] = vals.mean()
    return corr


def single_rho(r: np.ndarray) -> np.ndarray:
    return 0.5 * (PAULI["I"] + r[0] * PAULI["X"] + r[1] * PAULI["Y"] + r[2] * PAULI["Z"])


def pair_rho(ri: np.ndarray, rj: np.ndarray, corr: np.ndarray) -> np.ndarray:
    rho = np.kron(PAULI["I"], PAULI["I"]).astype(complex)
    for a, label_a in enumerate(AXES):
        rho += ri[a] * np.kron(PAULI[label_a], PAULI["I"])
        rho += rj[a] * np.kron(PAULI["I"], PAULI[label_a])
    for a, label_a in enumerate(AXES):
        for b, label_b in enumerate(AXES):
            rho += corr[a, b] * np.kron(PAULI[label_a], PAULI[label_b])
    return rho / 4.0


def mutual_information_matrix(raw_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = load_mat_v5(raw_path)
    raw = data["combdatasorted"]
    history = data["combtomohistory"].astype(np.uint8)
    axes, outcomes = decode_measurements(raw, history)
    n_sites = history.shape[1]

    r = estimate_single(axes, outcomes)
    s1 = np.array([entropy_density(single_rho(r[j])) for j in range(n_sites)])

    mi = np.full((n_sites, n_sites), np.nan, dtype=float)
    coverage = np.zeros((n_sites, n_sites), dtype=int)
    np.fill_diagonal(mi, 0.0)
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            corr = estimate_pair_corr(axes, outcomes, i, j)
            if corr is None:
                continue
            sij = entropy_density(pair_rho(r[i], r[j], corr))
            val = max(0.0, s1[i] + s1[j] - sij)
            mi[i, j] = mi[j, i] = val
            coverage[i, j] = coverage[j, i] = 1
    return mi, s1, coverage


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    def ranks(v):
        order = np.argsort(v)
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(v), dtype=float)
        return r

    return pearson(ranks(np.asarray(x)), ranks(np.asarray(y)))


def load_beta_profiles(beta_path: Path, group: str = "databulk") -> list[dict]:
    mat = load_mat_v5(beta_path)
    profiles = []
    for idx, (_, arr) in enumerate(mat[group]["cells"]):
        if not isinstance(arr, np.ndarray) or arr.shape[1] < 2:
            continue
        sites = arr[:, 0].astype(int)
        beta = arr[:, 1].astype(float)
        mask = beta > 1e-9
        if mask.sum() < 3:
            continue
        profiles.append(
            {
                "cell_index": idx,
                "n_rows": int(arr.shape[0]),
                "kind": "experimental" if arr.shape[1] >= 4 else "mps_or_fit",
                "sites": sites[mask],
                "beta": beta[mask],
                "inv_beta_norm": (1.0 / beta[mask]) / np.sum(1.0 / beta[mask]),
            }
        )
    return profiles


def compare_profiles(state: str, mi: np.ndarray, profiles: list[dict]) -> list[dict]:
    rho_sum = np.nansum(mi, axis=1)
    rows = []
    for prof in profiles:
        # Site labels in the paper are 1-based ion numbers.
        idx = prof["sites"] - 1
        mi_local = rho_sum[idx]
        mi_norm = mi_local / mi_local.sum() if mi_local.sum() > 0 else mi_local
        inv_beta = prof["inv_beta_norm"]
        rows.append(
            {
                "state": state,
                "cell_index": prof["cell_index"],
                "n_sites": int(len(idx)),
                "kind": prof["kind"],
                "test": "site_sumI_vs_invbeta",
                "pearson_sumI_invbeta": pearson(mi_norm, inv_beta),
                "spearman_sumI_invbeta": spearman(mi_norm, inv_beta),
                "mae_norm": float(np.mean(np.abs(mi_norm - inv_beta))),
                "sumI_edge_to_center": float(((mi_local[0] + mi_local[-1]) / 2.0) / mi_local[len(mi_local) // 2]),
                "invbeta_edge_to_center": float(((inv_beta[0] + inv_beta[-1]) / 2.0) / inv_beta[len(inv_beta) // 2]),
            }
        )
    return rows


def compare_bond_profiles(state: str, mi: np.ndarray, profiles: list[dict]) -> list[dict]:
    """Compare beta_j, which multiplies h_j, with nearest-neighbour MI I(j:j+1)."""
    rows = []
    for prof in profiles:
        sites = prof["sites"]
        inv_beta = prof["inv_beta_norm"]
        bond_mi = []
        keep_inv_beta = []
        for site, invb in zip(sites, inv_beta):
            if site < mi.shape[0] and np.isfinite(mi[site - 1, site]):
                bond_mi.append(mi[site - 1, site])
                keep_inv_beta.append(invb)
        bond_mi = np.asarray(bond_mi, dtype=float)
        keep_inv_beta = np.asarray(keep_inv_beta, dtype=float)
        if bond_mi.size == 0:
            continue
        bond_norm = bond_mi / bond_mi.sum() if bond_mi.sum() > 0 else bond_mi
        rows.append(
            {
                "state": state,
                "cell_index": prof["cell_index"],
                "n_sites": int(bond_mi.size),
                "kind": prof["kind"],
                "test": "nearest_neighbor_MI_vs_invbeta",
                "pearson_sumI_invbeta": pearson(bond_norm, keep_inv_beta),
                "spearman_sumI_invbeta": spearman(bond_norm, keep_inv_beta),
                "mae_norm": float(np.mean(np.abs(bond_norm - keep_inv_beta))),
                "sumI_edge_to_center": float(((bond_mi[0] + bond_mi[-1]) / 2.0) / bond_mi[len(bond_mi) // 2]),
                "invbeta_edge_to_center": float(((keep_inv_beta[0] + keep_inv_beta[-1]) / 2.0) / keep_inv_beta[len(keep_inv_beta) // 2]),
            }
        )
    return rows


def write_matrix_csv(path: Path, matrix: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["site"] + list(range(1, matrix.shape[1] + 1)))
        for i, row in enumerate(matrix, start=1):
            writer.writerow([i] + [f"{x:.10g}" for x in row])


def main() -> None:
    outdir = Path("/workspace/bup_zoller_tomography/results_raw_mi")
    outdir.mkdir(parents=True, exist_ok=True)

    jobs = [
        {
            "state": "GS_Delta_1",
            "raw": Path("/workspace/.cache/03-dataRawParam_GS_Delta_1.mat"),
            "beta": Path("/workspace/.cache/03-GroundStateBetas.mat"),
        },
        {
            "state": "ES_Delta_1",
            "raw": Path("/workspace/.cache/01-dataRawParam_ES_Delta_1.mat"),
            "beta": Path("/workspace/.cache/01-ExcitedStateBetas.mat"),
        },
    ]

    all_rows = []
    summary = {}
    for job in jobs:
        mi, s1, coverage = mutual_information_matrix(job["raw"])
        rho_sum = np.nansum(mi, axis=1)
        write_matrix_csv(outdir / f"{job['state']}_pairwise_MI.csv", np.nan_to_num(mi, nan=-1.0))
        write_matrix_csv(outdir / f"{job['state']}_pairwise_MI_coverage.csv", coverage.astype(float))
        np.savetxt(outdir / f"{job['state']}_sumI_per_site.csv", np.column_stack([np.arange(1, 52), rho_sum]), delimiter=",", header="site,sum_pairwise_MI", comments="")
        np.savetxt(outdir / f"{job['state']}_single_site_entropy.csv", np.column_stack([np.arange(1, 52), s1]), delimiter=",", header="site,S_i", comments="")

        profiles = load_beta_profiles(job["beta"], group="databulk")
        rows = compare_profiles(job["state"], mi, profiles)
        rows += compare_bond_profiles(job["state"], mi, profiles)
        all_rows.extend(rows)
        exp_rows = [r for r in rows if r["kind"] == "experimental" and r["test"] == "site_sumI_vs_invbeta"]
        bond_exp_rows = [r for r in rows if r["kind"] == "experimental" and r["test"] == "nearest_neighbor_MI_vs_invbeta"]
        summary[job["state"]] = {
            "site_sumI_mean_pearson_experimental": float(np.nanmean([r["pearson_sumI_invbeta"] for r in exp_rows])),
            "site_sumI_mean_spearman_experimental": float(np.nanmean([r["spearman_sumI_invbeta"] for r in exp_rows])),
            "site_sumI_mean_mae_norm_experimental": float(np.nanmean([r["mae_norm"] for r in exp_rows])),
            "site_sumI_mean_edge_to_center_experimental": float(np.nanmean([r["sumI_edge_to_center"] for r in exp_rows])),
            "nearest_neighbor_mean_pearson_experimental": float(np.nanmean([r["pearson_sumI_invbeta"] for r in bond_exp_rows])),
            "nearest_neighbor_mean_spearman_experimental": float(np.nanmean([r["spearman_sumI_invbeta"] for r in bond_exp_rows])),
            "nearest_neighbor_mean_mae_norm_experimental": float(np.nanmean([r["mae_norm"] for r in bond_exp_rows])),
            "nearest_neighbor_mean_edge_to_center_experimental": float(np.nanmean([r["sumI_edge_to_center"] for r in bond_exp_rows])),
            "mean_invbeta_edge_to_center_experimental": float(np.nanmean([r["invbeta_edge_to_center"] for r in exp_rows])),
            "n_experimental_profiles": len(exp_rows),
        }

    csv_path = outdir / "bup_raw_MI_vs_inv_beta_profile_tests.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    (outdir / "bup_raw_MI_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()

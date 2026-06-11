#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP Paper 13 — Microscopic closure scan v1

Goal
----
Robustness scan for the microscopic closure hypothesis:

    |Psi_gal> -> I_ij -> rho_ent^micro(R) ? proportional to Sigma(R)

This script scans the two most important microscopic control parameters:

    xi : correlation length in J_ij
    h0 : transverse field amplitude

For each pair (xi, h0), it:
    1. builds a small quantum disk,
    2. constructs J_ij ~ sqrt(Sigma_i Sigma_j) exp(-r_ij / xi),
    3. diagonalizes the transverse-Ising-like Hamiltonian,
    4. computes pairwise mutual information I_ij,
    5. computes rho_ent^micro(i) = sum_j I_ij,
    6. radially coarse-grains rho_ent^micro(R),
    7. compares it to Sigma(R),
    8. stores correlation, RMSE, KL and fitted R_ent.

Recommended run
---------------
python3 papers/paper13_micro_closure_bup/scripts/bup_paper13_micro_closure_scan_v1.py \
    --n-rings 3 \
    --n-theta 4 \
    --r-max 10.0 \
    --rd 3.0 \
    --xi-list 0.8 1.2 1.6 2.0 2.5 3.0 4.0 \
    --h0-list 0.1 0.2 0.3 0.5 0.8 1.0 \
    --j0 1.0 \
    --output-dir papers/paper13_micro_closure_bup/results/scan_v1

Notes
-----
Exact diagonalization scales as 2^N. Keep N <= 14 for comfortable runs.
For the first robustness scan, N = 3 rings x 4 theta = 12 qubits is recommended.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import eigsh
from scipy.optimize import curve_fit
from scipy.stats import pearsonr


EPS = 1e-12


@dataclass
class DiskData:
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    theta: np.ndarray
    ring_id: np.ndarray
    sigma: np.ndarray
    dist: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="BuP Paper 13 microscopic closure robustness scan v1"
    )

    parser.add_argument("--n-rings", type=int, default=3,
                        help="Number of radial rings. Default: 3")
    parser.add_argument("--n-theta", type=int, default=4,
                        help="Number of angular cells per ring. Default: 4")
    parser.add_argument("--r-max", type=float, default=10.0,
                        help="Disk maximum radius in arbitrary/kpc units. Default: 10")
    parser.add_argument("--rd", type=float, default=3.0,
                        help="Target exponential disk scale length R_d. Default: 3")
    parser.add_argument("--sigma0", type=float, default=1.0,
                        help="Central target surface density. Default: 1")
    parser.add_argument("--xi-list", type=float, nargs="+",
                        default=[0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0],
                        help="List of correlation lengths xi to scan.")
    parser.add_argument("--h0-list", type=float, nargs="+",
                        default=[0.1, 0.2, 0.3, 0.5, 0.8, 1.0],
                        help="List of transverse fields h0 to scan.")
    parser.add_argument("--j0", type=float, default=1.0,
                        help="Overall coupling amplitude. Default: 1")
    parser.add_argument("--max-qubits", type=int, default=14,
                        help="Safety limit for exact diagonalization. Default: 14")
    parser.add_argument("--output-dir", type=str, default="results_paper13_micro_closure_scan_v1",
                        help="Output directory.")
    parser.add_argument("--save-all-radial", action="store_true",
                        help="Save one radial profile CSV per parameter point.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Disable figure generation.")

    return parser.parse_args()


def make_quantum_disk(
    n_rings: int,
    n_theta: int,
    r_max: float,
    rd: float,
    sigma0: float,
) -> DiskData:
    xs, ys, rs, ths, ring_ids, sigmas = [], [], [], [], [], []
    dr = r_max / n_rings

    for k in range(n_rings):
        radius = (k + 0.5) * dr
        for m in range(n_theta):
            theta = 2.0 * np.pi * m / n_theta
            x = radius * np.cos(theta)
            y = radius * np.sin(theta)
            sigma = sigma0 * np.exp(-radius / rd)

            xs.append(x)
            ys.append(y)
            rs.append(radius)
            ths.append(theta)
            ring_ids.append(k)
            sigmas.append(sigma)

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    r = np.asarray(rs, dtype=float)
    theta = np.asarray(ths, dtype=float)
    ring_id = np.asarray(ring_ids, dtype=int)
    sigma = np.asarray(sigmas, dtype=float)

    coords = np.column_stack([x, y])
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(diff ** 2, axis=-1))

    return DiskData(
        x=x,
        y=y,
        r=r,
        theta=theta,
        ring_id=ring_id,
        sigma=sigma,
        dist=dist,
    )


def build_coupling_matrix(disk: DiskData, j0: float, xi: float) -> np.ndarray:
    J = j0 * np.sqrt(np.outer(disk.sigma, disk.sigma)) * np.exp(-disk.dist / xi)
    np.fill_diagonal(J, 0.0)
    return J


def pauli_sparse() -> Tuple[csr_matrix, csr_matrix, csr_matrix]:
    X = csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]], dtype=float))
    Z = csr_matrix(np.array([[1.0, 0.0], [0.0, -1.0]], dtype=float))
    I = identity(2, format="csr", dtype=float)
    return X, Z, I


def one_site_operator(op: csr_matrix, site: int, n: int) -> csr_matrix:
    _, _, I = pauli_sparse()

    result = None
    for k in range(n):
        factor = op if k == site else I
        result = factor if result is None else kron(result, factor, format="csr")
    return result


def two_site_xx_operator(site_i: int, site_j: int, n: int) -> csr_matrix:
    X, _, I = pauli_sparse()

    result = None
    for k in range(n):
        factor = X if (k == site_i or k == site_j) else I
        result = factor if result is None else kron(result, factor, format="csr")
    return result


def precompute_operators(n: int) -> Tuple[List[csr_matrix], Dict[Tuple[int, int], csr_matrix]]:
    """
    Precompute Z_i and X_i X_j operators once to speed up the scan.
    """

    _, Z, _ = pauli_sparse()

    z_ops = [one_site_operator(Z, i, n) for i in range(n)]

    xx_ops = {}
    for i in range(n):
        for j in range(i + 1, n):
            xx_ops[(i, j)] = two_site_xx_operator(i, j, n)

    return z_ops, xx_ops


def build_hamiltonian_from_precomputed(
    J: np.ndarray,
    h0: float,
    z_ops: List[csr_matrix],
    xx_ops: Dict[Tuple[int, int], csr_matrix],
) -> csr_matrix:
    n = J.shape[0]
    dim = 2 ** n
    H = csr_matrix((dim, dim), dtype=float)

    for i in range(n):
        H += -h0 * z_ops[i]

    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > EPS:
                H += -J[i, j] * xx_ops[(i, j)]

    return H


def ground_state(H: csr_matrix) -> Tuple[float, np.ndarray]:
    vals, vecs = eigsh(H, k=1, which="SA")
    e0 = float(vals[0])
    psi0 = np.asarray(vecs[:, 0], dtype=complex)
    psi0 = psi0 / np.linalg.norm(psi0)
    return e0, psi0


def reduced_density_matrix_pure_state(
    psi: np.ndarray,
    keep: List[int],
    n: int,
) -> np.ndarray:
    keep = list(keep)
    trace = [i for i in range(n) if i not in keep]

    psi_tensor = psi.reshape([2] * n)
    perm = keep + trace
    psi_perm = np.transpose(psi_tensor, axes=perm)

    dim_keep = 2 ** len(keep)
    dim_trace = 2 ** (n - len(keep))

    psi_matrix = psi_perm.reshape(dim_keep, dim_trace)
    rho = psi_matrix @ psi_matrix.conj().T
    rho = 0.5 * (rho + rho.conj().T)

    return rho


def von_neumann_entropy(rho: np.ndarray, base: float = 2.0) -> float:
    eigvals = np.linalg.eigvalsh(rho)
    eigvals = np.real(eigvals)
    eigvals = eigvals[eigvals > EPS]

    if len(eigvals) == 0:
        return 0.0

    logs = np.log(eigvals) / np.log(base)
    return float(-np.sum(eigvals * logs))


def mutual_information_matrix(psi: np.ndarray, n: int) -> Tuple[np.ndarray, np.ndarray]:
    S1 = np.zeros(n, dtype=float)

    for i in range(n):
        rho_i = reduced_density_matrix_pure_state(psi, [i], n)
        S1[i] = von_neumann_entropy(rho_i)

    MI = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            rho_ij = reduced_density_matrix_pure_state(psi, [i, j], n)
            Sij = von_neumann_entropy(rho_ij)
            mij = S1[i] + S1[j] - Sij

            if mij < 0 and abs(mij) < 1e-10:
                mij = 0.0

            MI[i, j] = mij
            MI[j, i] = mij

    return MI, S1


def radial_coarse_grain(values: np.ndarray, disk: DiskData, n_rings: int) -> pd.DataFrame:
    rows = []

    for k in range(n_rings):
        mask = disk.ring_id == k

        rows.append({
            "ring_id": k,
            "R_mean": float(np.mean(disk.r[mask])),
            "n_cells": int(np.sum(mask)),
            "sigma_mean": float(np.mean(disk.sigma[mask])),
            "rho_ent_mean": float(np.mean(values[mask])),
            "rho_ent_sum": float(np.sum(values[mask])),
        })

    return pd.DataFrame(rows)


def normalize_profile(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    total = np.sum(y)
    if total <= EPS:
        return np.zeros_like(y)
    return y / total


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = normalize_profile(p) + EPS
    q = normalize_profile(q) + EPS

    p = p / np.sum(p)
    q = q / np.sum(q)

    return float(np.sum(p * np.log(p / q)))


def exp_profile(R: np.ndarray, A: float, R_scale: float) -> np.ndarray:
    return A * np.exp(-R / R_scale)


def fit_exponential_scale(R: np.ndarray, y: np.ndarray, rd_guess: float) -> Tuple[float, float]:
    R = np.asarray(R, dtype=float)
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, EPS)

    try:
        popt, _ = curve_fit(
            exp_profile,
            R,
            y,
            p0=(float(np.max(y)), float(rd_guess)),
            bounds=([0.0, EPS], [np.inf, np.inf]),
            maxfev=10000,
        )
        return float(popt[0]), float(popt[1])
    except Exception:
        return float("nan"), float("nan")


def compute_diagnostics(radial_df: pd.DataFrame, rd_target: float) -> Dict[str, float]:
    sigma = radial_df["sigma_mean"].to_numpy(dtype=float)
    rho = radial_df["rho_ent_mean"].to_numpy(dtype=float)
    R = radial_df["R_mean"].to_numpy(dtype=float)

    sigma_n = normalize_profile(sigma)
    rho_n = normalize_profile(rho)

    if len(R) >= 2 and np.std(sigma_n) > EPS and np.std(rho_n) > EPS:
        corr = float(pearsonr(sigma_n, rho_n).statistic)
    else:
        corr = float("nan")

    rmse = float(np.sqrt(np.mean((rho_n - sigma_n) ** 2)))
    kl_rho_sigma = kl_divergence(rho, sigma)
    kl_sigma_rho = kl_divergence(sigma, rho)

    _, R_ent = fit_exponential_scale(R, rho, rd_target)
    _, R_sigma_fit = fit_exponential_scale(R, sigma, rd_target)

    if np.isfinite(R_ent):
        rd_rel_error = float(abs(R_ent - rd_target) / rd_target)
    else:
        rd_rel_error = float("nan")

    return {
        "corr_rho_sigma": corr,
        "rmse_normalized": rmse,
        "kl_rho_ent_to_sigma": kl_rho_sigma,
        "kl_sigma_to_rho_ent": kl_sigma_rho,
        "R_ent_fit": R_ent,
        "R_sigma_fit": R_sigma_fit,
        "R_d_target": float(rd_target),
        "R_d_relative_error": rd_rel_error,
    }


def closure_verdict(diag: Dict[str, float]) -> str:
    corr = diag["corr_rho_sigma"]
    rmse = diag["rmse_normalized"]
    rd_err = diag["R_d_relative_error"]

    if np.isfinite(corr) and np.isfinite(rmse) and np.isfinite(rd_err):
        if corr >= 0.90 and rmse <= 0.10 and rd_err <= 0.20:
            return "strong_closure"
        if corr >= 0.80 and rmse <= 0.15 and rd_err <= 0.35:
            return "moderate_closure"
        return "weak_or_failed_closure"

    return "insufficient_diagnostics"


def pivot_metric(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    return df.pivot(index="h0", columns="xi", values=metric).sort_index().sort_index(axis=1)


def save_heatmap(df: pd.DataFrame, metric: str, output_path: Path, title: str) -> None:
    pivot = pivot_metric(df, metric)

    x_vals = pivot.columns.to_numpy(dtype=float)
    y_vals = pivot.index.to_numpy(dtype=float)
    Z = pivot.to_numpy(dtype=float)

    plt.figure(figsize=(8, 5.5))
    im = plt.imshow(
        Z,
        origin="lower",
        aspect="auto",
        extent=[x_vals.min(), x_vals.max(), y_vals.min(), y_vals.max()],
    )
    plt.colorbar(im, label=metric)
    plt.xlabel("xi")
    plt.ylabel("h0")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_verdict_map(df: pd.DataFrame, output_path: Path) -> None:
    verdict_order = {
        "weak_or_failed_closure": 0,
        "insufficient_diagnostics": 0,
        "moderate_closure": 1,
        "strong_closure": 2,
    }

    tmp = df.copy()
    tmp["verdict_code"] = tmp["verdict"].map(verdict_order).fillna(0)

    pivot = pivot_metric(tmp, "verdict_code")

    x_vals = pivot.columns.to_numpy(dtype=float)
    y_vals = pivot.index.to_numpy(dtype=float)
    Z = pivot.to_numpy(dtype=float)

    plt.figure(figsize=(8, 5.5))
    im = plt.imshow(
        Z,
        origin="lower",
        aspect="auto",
        extent=[x_vals.min(), x_vals.max(), y_vals.min(), y_vals.max()],
        vmin=0,
        vmax=2,
    )
    cbar = plt.colorbar(im, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["weak/failed", "moderate", "strong"])
    plt.xlabel("xi")
    plt.ylabel("h0")
    plt.title("Microscopic closure verdict map")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def run_one_point(
    disk: DiskData,
    n: int,
    n_rings: int,
    rd: float,
    j0: float,
    xi: float,
    h0: float,
    z_ops: List[csr_matrix],
    xx_ops: Dict[Tuple[int, int], csr_matrix],
) -> Tuple[Dict[str, float], pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    J = build_coupling_matrix(disk, j0=j0, xi=xi)
    H = build_hamiltonian_from_precomputed(J, h0=h0, z_ops=z_ops, xx_ops=xx_ops)

    e0, psi0 = ground_state(H)
    MI, S1 = mutual_information_matrix(psi0, n)

    rho_ent_site = np.sum(MI, axis=1)
    radial_df = radial_coarse_grain(rho_ent_site, disk, n_rings)

    diagnostics = compute_diagnostics(radial_df, rd_target=rd)
    verdict = closure_verdict(diagnostics)

    upper = MI[np.triu_indices(n, k=1)]

    row = {
        "xi": float(xi),
        "h0": float(h0),
        "j0": float(j0),
        "ground_state_energy": float(e0),
        "mean_single_site_entropy": float(np.mean(S1)),
        "max_single_site_entropy": float(np.max(S1)),
        "mean_mutual_information": float(np.mean(upper)),
        "max_mutual_information": float(np.max(MI)),
        "corr_rho_sigma": diagnostics["corr_rho_sigma"],
        "rmse_normalized": diagnostics["rmse_normalized"],
        "kl_rho_ent_to_sigma": diagnostics["kl_rho_ent_to_sigma"],
        "kl_sigma_to_rho_ent": diagnostics["kl_sigma_to_rho_ent"],
        "R_ent_fit": diagnostics["R_ent_fit"],
        "R_sigma_fit": diagnostics["R_sigma_fit"],
        "R_d_target": diagnostics["R_d_target"],
        "R_d_relative_error": diagnostics["R_d_relative_error"],
        "verdict": verdict,
    }

    return row, radial_df, J, MI, rho_ent_site


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n = args.n_rings * args.n_theta
    if n > args.max_qubits:
        raise ValueError(
            f"Requested N={n} qubits, but max-qubits={args.max_qubits}. "
            f"Exact diagonalization scales as 2^N. Reduce --n-rings or --n-theta, "
            f"or increase --max-qubits if you know what you are doing."
        )

    disk = make_quantum_disk(
        n_rings=args.n_rings,
        n_theta=args.n_theta,
        r_max=args.r_max,
        rd=args.rd,
        sigma0=args.sigma0,
    )

    n_total = len(args.xi_list) * len(args.h0_list)

    print("\n=== BuP Paper 13 — microscopic closure scan v1 ===")
    print(f"N qubits/cells        : {n}")
    print(f"Hilbert dimension     : {2 ** n}")
    print(f"n_rings, n_theta      : {args.n_rings}, {args.n_theta}")
    print(f"R_d target            : {args.rd:.6g}")
    print(f"j0                    : {args.j0:.6g}")
    print(f"xi values             : {args.xi_list}")
    print(f"h0 values             : {args.h0_list}")
    print(f"Total runs            : {n_total}")
    print("\nPrecomputing operators...")

    z_ops, xx_ops = precompute_operators(n)

    radial_dir = output_dir / "radial_profiles_by_point"
    if args.save_all_radial:
        radial_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    best_by_rmse = None
    best_payload = None

    t0 = time.time()
    run_idx = 0

    for xi in args.xi_list:
        for h0 in args.h0_list:
            run_idx += 1
            print(f"\n[{run_idx}/{n_total}] Running xi={xi:.6g}, h0={h0:.6g}...")

            try:
                row, radial_df, J, MI, rho_ent_site = run_one_point(
                    disk=disk,
                    n=n,
                    n_rings=args.n_rings,
                    rd=args.rd,
                    j0=args.j0,
                    xi=float(xi),
                    h0=float(h0),
                    z_ops=z_ops,
                    xx_ops=xx_ops,
                )

                rows.append(row)

                print(
                    f"  corr={row['corr_rho_sigma']:.6f} | "
                    f"rmse={row['rmse_normalized']:.6f} | "
                    f"R_ent={row['R_ent_fit']:.6f} | "
                    f"rd_err={row['R_d_relative_error']:.6f} | "
                    f"{row['verdict']}"
                )

                if args.save_all_radial:
                    radial_name = f"radial_xi_{xi:.6g}_h0_{h0:.6g}.csv".replace(".", "p")
                    radial_df.to_csv(radial_dir / radial_name, index=False)

                if best_by_rmse is None or row["rmse_normalized"] < best_by_rmse["rmse_normalized"]:
                    best_by_rmse = row
                    best_payload = {
                        "radial_df": radial_df.copy(),
                        "J": J.copy(),
                        "MI": MI.copy(),
                        "rho_ent_site": rho_ent_site.copy(),
                    }

            except Exception as exc:
                print(f"  FAILED: {exc}")

                rows.append({
                    "xi": float(xi),
                    "h0": float(h0),
                    "j0": float(args.j0),
                    "ground_state_energy": float("nan"),
                    "mean_single_site_entropy": float("nan"),
                    "max_single_site_entropy": float("nan"),
                    "mean_mutual_information": float("nan"),
                    "max_mutual_information": float("nan"),
                    "corr_rho_sigma": float("nan"),
                    "rmse_normalized": float("nan"),
                    "kl_rho_ent_to_sigma": float("nan"),
                    "kl_sigma_to_rho_ent": float("nan"),
                    "R_ent_fit": float("nan"),
                    "R_sigma_fit": float("nan"),
                    "R_d_target": float(args.rd),
                    "R_d_relative_error": float("nan"),
                    "verdict": "failed",
                    "error": str(exc),
                })

    scan_df = pd.DataFrame(rows)
    scan_df.to_csv(output_dir / "scan_results.csv", index=False)

    n_success = int(np.sum(scan_df["verdict"] != "failed"))
    n_strong = int(np.sum(scan_df["verdict"] == "strong_closure"))
    n_moderate = int(np.sum(scan_df["verdict"] == "moderate_closure"))
    n_failed_or_weak = int(np.sum(scan_df["verdict"].isin(["weak_or_failed_closure", "failed", "insufficient_diagnostics"])))

    summary = {
        "experiment": "BuP Paper 13 microscopic closure scan v1",
        "n_qubits": int(n),
        "hilbert_dimension": int(2 ** n),
        "n_rings": int(args.n_rings),
        "n_theta": int(args.n_theta),
        "r_max": float(args.r_max),
        "rd_target": float(args.rd),
        "sigma0": float(args.sigma0),
        "j0": float(args.j0),
        "xi_list": [float(x) for x in args.xi_list],
        "h0_list": [float(x) for x in args.h0_list],
        "n_total_runs": int(n_total),
        "n_successful_runs": n_success,
        "n_strong_closure": n_strong,
        "n_moderate_closure": n_moderate,
        "n_weak_failed_or_insufficient": n_failed_or_weak,
        "fraction_strong_closure": float(n_strong / max(n_total, 1)),
        "fraction_moderate_or_strong": float((n_strong + n_moderate) / max(n_total, 1)),
        "best_by_rmse": None if best_by_rmse is None else {
            k: (float(v) if isinstance(v, (int, float, np.floating)) and np.isfinite(v) else str(v))
            for k, v in best_by_rmse.items()
        },
        "runtime_seconds": float(time.time() - t0),
        "outputs": {
            "scan_results": "scan_results.csv",
            "scan_summary": "scan_summary.json",
            "fig_corr_xi_h0": "fig_corr_xi_h0.png",
            "fig_rmse_xi_h0": "fig_rmse_xi_h0.png",
            "fig_rd_error_xi_h0": "fig_rd_error_xi_h0.png",
            "fig_verdict_map": "fig_verdict_map.png",
            "best_radial_profile": "best_radial_profile.csv",
            "best_mutual_information_matrix": "best_mutual_information_matrix.csv",
            "best_coupling_matrix": "best_coupling_matrix.csv",
        },
    }

    with open(output_dir / "scan_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    if best_payload is not None and best_by_rmse is not None:
        best_payload["radial_df"].to_csv(output_dir / "best_radial_profile.csv", index=False)
        pd.DataFrame(best_payload["MI"]).to_csv(output_dir / "best_mutual_information_matrix.csv", index=False)
        pd.DataFrame(best_payload["J"]).to_csv(output_dir / "best_coupling_matrix.csv", index=False)

    if not args.no_plots and len(scan_df) > 0:
        valid_df = scan_df.replace([np.inf, -np.inf], np.nan)

        try:
            save_heatmap(
                valid_df,
                "corr_rho_sigma",
                output_dir / "fig_corr_xi_h0.png",
                r"Closure correlation: corr($\rho_{\rm ent}^{micro}$, $\Sigma$)",
            )
        except Exception as exc:
            print(f"Could not save correlation heatmap: {exc}")

        try:
            save_heatmap(
                valid_df,
                "rmse_normalized",
                output_dir / "fig_rmse_xi_h0.png",
                "Normalized RMSE map",
            )
        except Exception as exc:
            print(f"Could not save RMSE heatmap: {exc}")

        try:
            save_heatmap(
                valid_df,
                "R_d_relative_error",
                output_dir / "fig_rd_error_xi_h0.png",
                "Relative scale-length error map",
            )
        except Exception as exc:
            print(f"Could not save R_d error heatmap: {exc}")

        try:
            save_verdict_map(
                valid_df,
                output_dir / "fig_verdict_map.png",
            )
        except Exception as exc:
            print(f"Could not save verdict map: {exc}")

    print("\n=== Scan summary ===")
    print(f"Total runs                 : {n_total}")
    print(f"Successful runs            : {n_success}")
    print(f"Strong closure             : {n_strong}")
    print(f"Moderate closure           : {n_moderate}")
    print(f"Weak/failed/insufficient   : {n_failed_or_weak}")
    print(f"Fraction strong            : {summary['fraction_strong_closure']:.3f}")
    print(f"Fraction moderate+strong   : {summary['fraction_moderate_or_strong']:.3f}")

    if best_by_rmse is not None:
        print("\nBest point by RMSE:")
        print(f"  xi                    : {best_by_rmse['xi']}")
        print(f"  h0                    : {best_by_rmse['h0']}")
        print(f"  corr                  : {best_by_rmse['corr_rho_sigma']:.6f}")
        print(f"  rmse                  : {best_by_rmse['rmse_normalized']:.6f}")
        print(f"  R_ent                 : {best_by_rmse['R_ent_fit']:.6f}")
        print(f"  R_d relative error    : {best_by_rmse['R_d_relative_error']:.6f}")
        print(f"  verdict               : {best_by_rmse['verdict']}")

    print("\nFiles written:")
    expected = [
        "scan_results.csv",
        "scan_summary.json",
        "fig_corr_xi_h0.png",
        "fig_rmse_xi_h0.png",
        "fig_rd_error_xi_h0.png",
        "fig_verdict_map.png",
        "best_radial_profile.csv",
        "best_mutual_information_matrix.csv",
        "best_coupling_matrix.csv",
    ]
    for name in expected:
        if (output_dir / name).exists():
            print(f"  - {output_dir / name}")

    print("\nDone.")


if __name__ == "__main__":
    main()

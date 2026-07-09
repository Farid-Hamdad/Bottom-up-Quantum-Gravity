#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Circle spectrum convergence v1

Goal
----
Test the continuum-limit normalization of the unnormalized graph Laplacian
on the unit circle S^1.

We build:

    W_ij = exp[-d_S1(theta_i, theta_j)^2 / (4 epsilon)]

with diagonal removed, and:

    L = D - W.

For S^1, the continuum spectrum of -d^2/dtheta^2 is:

    lambda_k = k^2,

with multiplicity 2 for k >= 1.

We estimate a normalization c_{N,epsilon} using the first non-zero spectral
band and test whether:

    c_{N,epsilon} lambda_k(L) ~ k^2.

Outputs
-------
paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1/
  circle_spectrum_convergence.csv
  normalization_scan.csv
  summary.json
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd


OUTDIR = Path("paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1")
OUTDIR.mkdir(parents=True, exist_ok=True)


def periodic_distance(theta: np.ndarray) -> np.ndarray:
    """
    Pairwise geodesic distance on the unit circle.

    theta is in [0, 2pi).
    distance is min(|theta_i-theta_j|, 2pi-|theta_i-theta_j|).
    """
    diff = np.abs(theta[:, None] - theta[None, :])
    return np.minimum(diff, 2.0 * np.pi - diff)


def build_circle_laplacian(N: int, epsilon: float) -> np.ndarray:
    """
    Build unnormalized Gaussian graph Laplacian on S^1.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)
    dist = periodic_distance(theta)

    W = np.exp(-(dist ** 2) / (4.0 * epsilon))
    np.fill_diagonal(W, 0.0)

    D = np.diag(W.sum(axis=1))
    L = D - W
    return L


def continuum_circle_eigenvalues(num_modes: int) -> np.ndarray:
    """
    Continuum eigenvalues of -d^2/dtheta^2 on S^1.

    Ordering:
        0,
        1,1,
        4,4,
        9,9,
        ...
    """
    vals = [0.0]
    k = 1
    while len(vals) < num_modes:
        vals.extend([float(k * k), float(k * k)])
        k += 1
    return np.array(vals[:num_modes], dtype=float)


def spectral_error(normalized_vals: np.ndarray, target_vals: np.ndarray, start: int = 1) -> dict:
    """
    Compute low-mode errors, excluding the zero mode by default.
    """
    a = normalized_vals[start:]
    b = target_vals[start:]

    abs_err = np.abs(a - b)
    rel_err = abs_err / np.maximum(np.abs(b), 1e-12)

    return {
        "mean_abs_error": float(np.mean(abs_err)),
        "max_abs_error": float(np.max(abs_err)),
        "mean_rel_error": float(np.mean(rel_err)),
        "max_rel_error": float(np.max(rel_err)),
    }


def run_single(N: int, epsilon: float, num_modes: int = 11) -> dict:
    """
    Build L, diagonalize, estimate normalization using the k=1 band,
    then compare low modes to continuum.
    """
    L = build_circle_laplacian(N=N, epsilon=epsilon)

    eigvals = np.linalg.eigvalsh(L)
    eigvals = np.sort(np.real(eigvals))

    target = continuum_circle_eigenvalues(num_modes)

    # First non-zero band on S^1 has multiplicity 2:
    # indices 1 and 2 should correspond to k=1, lambda=1.
    first_band_mean = float(np.mean(eigvals[1:3]))
    c_empirical = 1.0 / first_band_mean if first_band_mean > 0 else np.nan

    normalized = c_empirical * eigvals[:num_modes]
    errors = spectral_error(normalized, target, start=1)

    # Analytic scaling expected for D=1:
    # c ~ const / (N * epsilon^(3/2))
    scaling_core = N * (epsilon ** 1.5)
    prefactor_estimate = c_empirical * scaling_core

    row = {
        "N": N,
        "epsilon": epsilon,
        "num_modes": num_modes,
        "lambda0_graph": float(eigvals[0]),
        "lambda1_graph": float(eigvals[1]),
        "lambda2_graph": float(eigvals[2]),
        "lambda1_band_mean": first_band_mean,
        "c_empirical": float(c_empirical),
        "scaling_core_N_eps_3_over_2": float(scaling_core),
        "prefactor_estimate_c_times_N_eps_3_over_2": float(prefactor_estimate),
        **errors,
    }

    # Add mode-by-mode comparison.
    for i in range(num_modes):
        row[f"graph_lambda_{i}"] = float(eigvals[i])
        row[f"norm_lambda_{i}"] = float(normalized[i])
        row[f"target_lambda_{i}"] = float(target[i])
        row[f"abs_error_{i}"] = float(abs(normalized[i] - target[i]))

    return row


def main() -> None:
    Ns = [64, 128, 256, 512]
    epsilons = [0.0025, 0.005, 0.01, 0.02, 0.04]
    num_modes = 11

    rows = []

    print("=" * 80)
    print("Paper 26 — Circle spectrum convergence v1")
    print("=" * 80)

    for N in Ns:
        for epsilon in epsilons:
            # Avoid too narrow kernels for small N, which can become numerically poor.
            dx = 2.0 * np.pi / N
            if np.sqrt(epsilon) < 0.35 * dx:
                continue

            row = run_single(N=N, epsilon=epsilon, num_modes=num_modes)
            rows.append(row)

            print(
                f"N={N:4d} eps={epsilon:7.4f} "
                f"c={row['c_empirical']:.6e} "
                f"pref={row['prefactor_estimate_c_times_N_eps_3_over_2']:.6e} "
                f"mean_rel={row['mean_rel_error']:.6e} "
                f"max_rel={row['max_rel_error']:.6e}"
            )

    df = pd.DataFrame(rows)

    spectrum_csv = OUTDIR / "circle_spectrum_convergence.csv"
    normalization_csv = OUTDIR / "normalization_scan.csv"
    summary_json = OUTDIR / "summary.json"

    df.to_csv(spectrum_csv, index=False)

    norm_cols = [
        "N",
        "epsilon",
        "c_empirical",
        "scaling_core_N_eps_3_over_2",
        "prefactor_estimate_c_times_N_eps_3_over_2",
        "mean_rel_error",
        "max_rel_error",
    ]
    df[norm_cols].to_csv(normalization_csv, index=False)

    best_idx = df["mean_rel_error"].idxmin()
    best = df.loc[best_idx].to_dict()

    summary = {
        "experiment": "paper26_circle_spectrum_convergence_v1",
        "geometry": "unit circle S^1",
        "kernel": "periodic Gaussian exp(-d^2/(4 epsilon))",
        "laplacian": "unnormalized L = D - W, diagonal removed",
        "normalization": "empirical c estimated from first nonzero spectral band",
        "expected_scaling": "c_{N,epsilon} proportional to 1/(N epsilon^(3/2)) for D=1",
        "num_rows": int(len(df)),
        "Ns": Ns,
        "epsilons": epsilons,
        "num_modes": num_modes,
        "best_by_mean_rel_error": {
            "N": int(best["N"]),
            "epsilon": float(best["epsilon"]),
            "c_empirical": float(best["c_empirical"]),
            "prefactor_estimate": float(best["prefactor_estimate_c_times_N_eps_3_over_2"]),
            "mean_rel_error": float(best["mean_rel_error"]),
            "max_rel_error": float(best["max_rel_error"]),
        },
        "outputs": {
            "circle_spectrum_convergence_csv": str(spectrum_csv),
            "normalization_scan_csv": str(normalization_csv),
            "summary_json": str(summary_json),
        },
    }

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("-" * 80)
    print(f"[OK] wrote {spectrum_csv}")
    print(f"[OK] wrote {normalization_csv}")
    print(f"[OK] wrote {summary_json}")
    print("-" * 80)
    print("Best row:")
    print(
        f"N={summary['best_by_mean_rel_error']['N']}, "
        f"epsilon={summary['best_by_mean_rel_error']['epsilon']}, "
        f"mean_rel_error={summary['best_by_mean_rel_error']['mean_rel_error']:.6e}, "
        f"max_rel_error={summary['best_by_mean_rel_error']['max_rel_error']:.6e}"
    )


if __name__ == "__main__":
    main()

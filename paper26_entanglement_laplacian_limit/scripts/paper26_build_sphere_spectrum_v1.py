#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Sphere spectrum convergence v1

Goal
----
Test the continuum-limit normalization of the unnormalized graph Laplacian
on the unit sphere S^2.

Kernel:
    W_ij = exp[-d_S2(x_i,x_j)^2 / (4 epsilon)]

Graph Laplacian:
    L = D - W

Continuum spectrum of -Delta on unit S^2:
    lambda_l = l(l+1), multiplicity 2l+1.

General normalization:
    c_{N,epsilon} = 1 / [rho (4pi)^(D/2) epsilon^(D/2+1)]

For unit S^2:
    D = 2
    Vol(S^2) = 4pi
    rho = N/(4pi)

Therefore:
    c_{N,epsilon}^{S2} = 1/(N epsilon^2)

Outputs
-------
paper26_entanglement_laplian_limit/results/paper26_laplacian_limit_v1/
  sphere_spectrum_convergence.csv
  sphere_prefactor_theory.csv
  sphere_summary.json
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd


OUTDIR = Path("paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1")
OUTDIR.mkdir(parents=True, exist_ok=True)


def fibonacci_sphere(N: int) -> np.ndarray:
    """
    Quasi-uniform points on the unit sphere using the Fibonacci lattice.
    Returns array shape (N, 3).
    """
    i = np.arange(N, dtype=float)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))

    z = 1.0 - 2.0 * (i + 0.5) / N
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = golden_angle * i

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    pts = np.column_stack([x, y, z])
    return pts


def sphere_geodesic_distances(pts: np.ndarray) -> np.ndarray:
    """
    Pairwise geodesic distances on unit S^2.
    """
    cosang = pts @ pts.T
    cosang = np.clip(cosang, -1.0, 1.0)
    return np.arccos(cosang)


def build_sphere_laplacian(N: int, epsilon: float) -> np.ndarray:
    pts = fibonacci_sphere(N)
    dist = sphere_geodesic_distances(pts)

    W = np.exp(-(dist ** 2) / (4.0 * epsilon))
    np.fill_diagonal(W, 0.0)

    D = np.diag(W.sum(axis=1))
    L = D - W
    return L


def continuum_sphere_eigenvalues(num_modes: int) -> np.ndarray:
    vals = []
    ell = 0
    while len(vals) < num_modes:
        lam = float(ell * (ell + 1))
        multiplicity = 2 * ell + 1
        vals.extend([lam] * multiplicity)
        ell += 1
    return np.array(vals[:num_modes], dtype=float)


def spectral_error(normalized: np.ndarray, target: np.ndarray, start: int = 1) -> dict:
    a = normalized[start:]
    b = target[start:]

    abs_err = np.abs(a - b)
    rel_err = abs_err / np.maximum(np.abs(b), 1e-12)

    return {
        "mean_abs_error": float(np.mean(abs_err)),
        "max_abs_error": float(np.max(abs_err)),
        "mean_rel_error": float(np.mean(rel_err)),
        "max_rel_error": float(np.max(rel_err)),
    }


def run_single(N: int, epsilon: float, num_modes: int = 25) -> dict:
    L = build_sphere_laplacian(N=N, epsilon=epsilon)

    eigvals = np.linalg.eigvalsh(L)
    eigvals = np.sort(np.real(eigvals))

    target = continuum_sphere_eigenvalues(num_modes=num_modes)

    # First non-zero band on S^2:
    # ell = 1, lambda = 2, multiplicity 3.
    first_band_mean = float(np.mean(eigvals[1:4]))

    c_empirical = 2.0 / first_band_mean
    c_theory = 1.0 / (N * epsilon * epsilon)

    norm_emp = c_empirical * eigvals[:num_modes]
    norm_theory = c_theory * eigvals[:num_modes]

    err_emp = spectral_error(norm_emp, target, start=1)
    err_theory = spectral_error(norm_theory, target, start=1)

    scaling_core = N * epsilon * epsilon
    pref_emp = c_empirical * scaling_core
    pref_theory = c_theory * scaling_core

    row = {
        "N": N,
        "epsilon": epsilon,
        "num_modes": num_modes,
        "sqrt_epsilon": float(np.sqrt(epsilon)),
        "lambda0_graph": float(eigvals[0]),
        "lambda1_graph": float(eigvals[1]),
        "lambda3_graph": float(eigvals[3]),
        "lambda1_band_mean": first_band_mean,
        "c_empirical": float(c_empirical),
        "c_theory_1_over_N_eps2": float(c_theory),
        "scaling_core_N_eps2": float(scaling_core),
        "prefactor_empirical_c_times_N_eps2": float(pref_emp),
        "prefactor_theory_1": float(pref_theory),
        "prefactor_rel_error_to_1": float(abs(pref_emp - 1.0)),
        "mean_rel_error_empirical_c": err_emp["mean_rel_error"],
        "max_rel_error_empirical_c": err_emp["max_rel_error"],
        "mean_rel_error_theory_c": err_theory["mean_rel_error"],
        "max_rel_error_theory_c": err_theory["max_rel_error"],
    }

    for i in range(num_modes):
        row[f"graph_lambda_{i}"] = float(eigvals[i])
        row[f"norm_emp_lambda_{i}"] = float(norm_emp[i])
        row[f"norm_theory_lambda_{i}"] = float(norm_theory[i])
        row[f"target_lambda_{i}"] = float(target[i])
        row[f"abs_error_emp_{i}"] = float(abs(norm_emp[i] - target[i]))
        row[f"abs_error_theory_{i}"] = float(abs(norm_theory[i] - target[i]))

    return row


def main() -> None:
    Ns = [256, 384, 512, 768, 1024]
    epsilons = [0.01, 0.02, 0.04, 0.08]
    num_modes = 25

    rows = []

    print("=" * 80)
    print("Paper 26 — Sphere spectrum convergence v1")
    print("=" * 80)

    for N in Ns:
        # Typical spacing on unit sphere is approximately sqrt(4pi/N).
        spacing = np.sqrt(4.0 * np.pi / N)

        for epsilon in epsilons:
            # Avoid kernels much narrower than point spacing.
            if np.sqrt(epsilon) / spacing < 0.8:
                continue

            row = run_single(N=N, epsilon=epsilon, num_modes=num_modes)
            rows.append(row)

            print(
                f"N={N:4d} eps={epsilon:6.3f} "
                f"sqrt_eps/spacing={np.sqrt(epsilon)/spacing:.2f} "
                f"pref_emp={row['prefactor_empirical_c_times_N_eps2']:.6f} "
                f"target=1.000000 "
                f"pref_abs_err={row['prefactor_rel_error_to_1']:.3e} "
                f"mean_rel_theory={row['mean_rel_error_theory_c']:.3e} "
                f"max_rel_theory={row['max_rel_error_theory_c']:.3e}"
            )

    df = pd.DataFrame(rows)

    spectrum_csv = OUTDIR / "sphere_spectrum_convergence.csv"
    prefactor_csv = OUTDIR / "sphere_prefactor_theory.csv"
    summary_json = OUTDIR / "sphere_summary.json"

    df.to_csv(spectrum_csv, index=False)

    pref_cols = [
        "N",
        "epsilon",
        "sqrt_epsilon",
        "c_empirical",
        "c_theory_1_over_N_eps2",
        "prefactor_empirical_c_times_N_eps2",
        "prefactor_theory_1",
        "prefactor_rel_error_to_1",
        "mean_rel_error_empirical_c",
        "max_rel_error_empirical_c",
        "mean_rel_error_theory_c",
        "max_rel_error_theory_c",
    ]
    df[pref_cols].to_csv(prefactor_csv, index=False)

    best_pref = df.loc[df["prefactor_rel_error_to_1"].idxmin()].to_dict()
    best_theory = df.loc[df["mean_rel_error_theory_c"].idxmin()].to_dict()

    summary = {
        "experiment": "paper26_sphere_spectrum_convergence_v1",
        "geometry": "unit sphere S^2",
        "kernel": "geodesic Gaussian exp(-d^2/(4 epsilon))",
        "laplacian": "unnormalized L = D - W, diagonal removed",
        "sampling": "Fibonacci quasi-uniform sphere",
        "expected_scaling": "c_{N,epsilon} = 1/(N epsilon^2) for unit S^2",
        "theory_prefactor": 1.0,
        "num_rows": int(len(df)),
        "Ns": Ns,
        "epsilons": epsilons,
        "num_modes": num_modes,
        "best_by_prefactor_error": {
            "N": int(best_pref["N"]),
            "epsilon": float(best_pref["epsilon"]),
            "empirical_prefactor": float(best_pref["prefactor_empirical_c_times_N_eps2"]),
            "absolute_error_to_1": float(best_pref["prefactor_rel_error_to_1"]),
            "mean_rel_error_theory_c": float(best_pref["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_pref["max_rel_error_theory_c"]),
        },
        "best_by_theory_spectral_error": {
            "N": int(best_theory["N"]),
            "epsilon": float(best_theory["epsilon"]),
            "empirical_prefactor": float(best_theory["prefactor_empirical_c_times_N_eps2"]),
            "absolute_error_to_1": float(best_theory["prefactor_rel_error_to_1"]),
            "mean_rel_error_theory_c": float(best_theory["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_theory["max_rel_error_theory_c"]),
        },
        "outputs": {
            "sphere_spectrum_convergence_csv": str(spectrum_csv),
            "sphere_prefactor_theory_csv": str(prefactor_csv),
            "sphere_summary_json": str(summary_json),
        },
    }

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("-" * 80)
    print(f"[OK] wrote {spectrum_csv}")
    print(f"[OK] wrote {prefactor_csv}")
    print(f"[OK] wrote {summary_json}")
    print("-" * 80)
    print("Best prefactor row:")
    print(summary["best_by_prefactor_error"])
    print("Best theory-spectrum row:")
    print(summary["best_by_theory_spectral_error"])


if __name__ == "__main__":
    main()

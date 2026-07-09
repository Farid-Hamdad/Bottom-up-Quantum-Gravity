#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Torus spectrum convergence v1

Goal
----
Test the D=2 continuum-limit normalization of the unnormalized graph Laplacian
on the flat torus T^2 = [0,2pi)^2.

Kernel:

    W_ij = exp[-d_T2(x_i,x_j)^2 / (4 epsilon)]

Graph Laplacian:

    L = D - W

Continuum spectrum of -Delta on the flat square torus:

    lambda_{m,n} = m^2 + n^2.

Expected normalization:

    c_{N,epsilon}^{T2} = 2*pi / (N epsilon^2)

for the current convention and uniform sampling on area (2pi)^2.

Outputs
-------
paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1/
  torus_spectrum_convergence.csv
  torus_prefactor_theory.csv
  torus_summary.json
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd


OUTDIR = Path("paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1")
OUTDIR.mkdir(parents=True, exist_ok=True)


def periodic_delta(a: np.ndarray) -> np.ndarray:
    return np.minimum(np.abs(a), 2.0 * np.pi - np.abs(a))


def make_torus_grid(n_side: int) -> np.ndarray:
    theta = np.linspace(0.0, 2.0 * np.pi, n_side, endpoint=False)
    X, Y = np.meshgrid(theta, theta, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel()])
    return pts


def torus_pairwise_distance_squared(pts: np.ndarray) -> np.ndarray:
    dx = periodic_delta(pts[:, 0][:, None] - pts[:, 0][None, :])
    dy = periodic_delta(pts[:, 1][:, None] - pts[:, 1][None, :])
    return dx * dx + dy * dy


def build_torus_laplacian(n_side: int, epsilon: float) -> np.ndarray:
    pts = make_torus_grid(n_side)
    dist2 = torus_pairwise_distance_squared(pts)

    W = np.exp(-dist2 / (4.0 * epsilon))
    np.fill_diagonal(W, 0.0)

    D = np.diag(W.sum(axis=1))
    L = D - W
    return L


def continuum_torus_eigenvalues(num_modes: int, max_k: int = 8) -> np.ndarray:
    vals = []
    for m in range(-max_k, max_k + 1):
        for n in range(-max_k, max_k + 1):
            vals.append(float(m * m + n * n))
    vals = sorted(vals)
    return np.array(vals[:num_modes], dtype=float)


def spectral_error(normalized_vals: np.ndarray, target_vals: np.ndarray, start: int = 1) -> dict:
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


def run_single(n_side: int, epsilon: float, num_modes: int = 20) -> dict:
    N = n_side * n_side

    L = build_torus_laplacian(n_side=n_side, epsilon=epsilon)

    eigvals = np.linalg.eigvalsh(L)
    eigvals = np.sort(np.real(eigvals))

    target = continuum_torus_eigenvalues(num_modes=num_modes)

    # The first non-zero band on T^2 has multiplicity 4:
    # (±1,0), (0,±1), eigenvalue 1.
    first_band_mean = float(np.mean(eigvals[1:5]))
    c_empirical = 1.0 / first_band_mean if first_band_mean > 0 else np.nan

    c_theory = 2.0 * np.pi / (N * epsilon * epsilon)

    normalized_emp = c_empirical * eigvals[:num_modes]
    normalized_theory = c_theory * eigvals[:num_modes]

    err_emp = spectral_error(normalized_emp, target, start=1)
    err_theory = spectral_error(normalized_theory, target, start=1)

    scaling_core = N * epsilon * epsilon
    prefactor_emp = c_empirical * scaling_core
    prefactor_theory = c_theory * scaling_core

    row = {
        "n_side": n_side,
        "N": N,
        "epsilon": epsilon,
        "num_modes": num_modes,
        "lambda0_graph": float(eigvals[0]),
        "lambda1_graph": float(eigvals[1]),
        "lambda4_graph": float(eigvals[4]),
        "lambda1_band_mean": first_band_mean,
        "c_empirical": float(c_empirical),
        "c_theory_2pi_over_N_eps2": float(c_theory),
        "scaling_core_N_eps2": float(scaling_core),
        "prefactor_empirical_c_times_N_eps2": float(prefactor_emp),
        "prefactor_theory_2pi": float(prefactor_theory),
        "prefactor_rel_error_to_2pi": float(abs(prefactor_emp - 2.0 * np.pi) / (2.0 * np.pi)),
        "mean_abs_error_empirical_c": err_emp["mean_abs_error"],
        "max_abs_error_empirical_c": err_emp["max_abs_error"],
        "mean_rel_error_empirical_c": err_emp["mean_rel_error"],
        "max_rel_error_empirical_c": err_emp["max_rel_error"],
        "mean_abs_error_theory_c": err_theory["mean_abs_error"],
        "max_abs_error_theory_c": err_theory["max_abs_error"],
        "mean_rel_error_theory_c": err_theory["mean_rel_error"],
        "max_rel_error_theory_c": err_theory["max_rel_error"],
    }

    for i in range(num_modes):
        row[f"graph_lambda_{i}"] = float(eigvals[i])
        row[f"norm_emp_lambda_{i}"] = float(normalized_emp[i])
        row[f"norm_theory_lambda_{i}"] = float(normalized_theory[i])
        row[f"target_lambda_{i}"] = float(target[i])
        row[f"abs_error_emp_{i}"] = float(abs(normalized_emp[i] - target[i]))
        row[f"abs_error_theory_{i}"] = float(abs(normalized_theory[i] - target[i]))

    return row


def main() -> None:
    n_sides = [16, 20, 24, 28]
    epsilons = [0.01, 0.02, 0.04, 0.08]
    num_modes = 20

    rows = []

    print("=" * 80)
    print("Paper 26 — Torus spectrum convergence v1")
    print("=" * 80)

    for n_side in n_sides:
        for epsilon in epsilons:
            dx = 2.0 * np.pi / n_side
            if np.sqrt(epsilon) < 0.30 * dx:
                continue

            row = run_single(n_side=n_side, epsilon=epsilon, num_modes=num_modes)
            rows.append(row)

            print(
                f"n_side={n_side:3d} N={row['N']:4d} eps={epsilon:7.4f} "
                f"pref_emp={row['prefactor_empirical_c_times_N_eps2']:.6f} "
                f"2pi={2*np.pi:.6f} "
                f"pref_rel={row['prefactor_rel_error_to_2pi']:.3e} "
                f"mean_rel_theory={row['mean_rel_error_theory_c']:.3e} "
                f"max_rel_theory={row['max_rel_error_theory_c']:.3e}"
            )

    df = pd.DataFrame(rows)

    spectrum_csv = OUTDIR / "torus_spectrum_convergence.csv"
    prefactor_csv = OUTDIR / "torus_prefactor_theory.csv"
    summary_json = OUTDIR / "torus_summary.json"

    df.to_csv(spectrum_csv, index=False)

    pref_cols = [
        "n_side",
        "N",
        "epsilon",
        "c_empirical",
        "c_theory_2pi_over_N_eps2",
        "prefactor_empirical_c_times_N_eps2",
        "prefactor_theory_2pi",
        "prefactor_rel_error_to_2pi",
        "mean_rel_error_empirical_c",
        "max_rel_error_empirical_c",
        "mean_rel_error_theory_c",
        "max_rel_error_theory_c",
    ]
    df[pref_cols].to_csv(prefactor_csv, index=False)

    best_pref = df.loc[df["prefactor_rel_error_to_2pi"].idxmin()].to_dict()
    best_theory = df.loc[df["mean_rel_error_theory_c"].idxmin()].to_dict()

    summary = {
        "experiment": "paper26_torus_spectrum_convergence_v1",
        "geometry": "flat square torus T^2 = [0,2pi)^2",
        "kernel": "periodic Gaussian exp(-d^2/(4 epsilon))",
        "laplacian": "unnormalized L = D - W, diagonal removed",
        "expected_scaling": "c_{N,epsilon} = 2*pi/(N epsilon^2) for D=2 and this convention",
        "num_rows": int(len(df)),
        "n_sides": n_sides,
        "epsilons": epsilons,
        "num_modes": num_modes,
        "theory_prefactor": float(2.0 * np.pi),
        "best_by_prefactor_error": {
            "n_side": int(best_pref["n_side"]),
            "N": int(best_pref["N"]),
            "epsilon": float(best_pref["epsilon"]),
            "empirical_prefactor": float(best_pref["prefactor_empirical_c_times_N_eps2"]),
            "relative_error_to_2pi": float(best_pref["prefactor_rel_error_to_2pi"]),
            "mean_rel_error_theory_c": float(best_pref["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_pref["max_rel_error_theory_c"]),
        },
        "best_by_theory_spectral_error": {
            "n_side": int(best_theory["n_side"]),
            "N": int(best_theory["N"]),
            "epsilon": float(best_theory["epsilon"]),
            "empirical_prefactor": float(best_theory["prefactor_empirical_c_times_N_eps2"]),
            "relative_error_to_2pi": float(best_theory["prefactor_rel_error_to_2pi"]),
            "mean_rel_error_theory_c": float(best_theory["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_theory["max_rel_error_theory_c"]),
        },
        "outputs": {
            "torus_spectrum_convergence_csv": str(spectrum_csv),
            "torus_prefactor_theory_csv": str(prefactor_csv),
            "torus_summary_json": str(summary_json),
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

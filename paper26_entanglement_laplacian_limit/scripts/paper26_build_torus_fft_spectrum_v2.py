#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Torus FFT spectrum convergence v2

This script corrects and improves the first T^2 test.

Geometry:
    T^2 = [0, 2pi)^2

Kernel:
    W_ij = exp[-d_T2(x_i,x_j)^2 / (4 epsilon)]

Laplacian:
    L = D - W

For a regular periodic grid, the graph Laplacian is block-circulant.
Its eigenvectors are Fourier modes, and its eigenvalues are obtained by FFT.

Analytic continuum expansion:
    L_{N,epsilon} f ≈ -(N/pi) epsilon^2 Delta f

Therefore:
    c_{N,epsilon}^{T2} = pi / (N epsilon^2)

Target:
    c_{N,epsilon} lambda_{m,n}(L) ≈ m^2+n^2
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd


OUTDIR = Path("paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1")
OUTDIR.mkdir(parents=True, exist_ok=True)


def periodic_grid_offsets(n_side: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Periodic coordinate offsets from the origin on [0, 2pi)^2.
    """
    dx = 2.0 * np.pi / n_side
    idx = np.arange(n_side)

    # Minimal periodic integer offsets: 0,1,2,...,n/2,-n/2+1,...
    offsets = np.where(idx <= n_side // 2, idx, idx - n_side)
    x = offsets * dx

    X, Y = np.meshgrid(x, x, indexing="ij")
    return X, Y


def torus_kernel_first_row(n_side: int, epsilon: float) -> np.ndarray:
    """
    First row of the periodic Gaussian kernel as an n_side x n_side array.
    """
    X, Y = periodic_grid_offsets(n_side)
    dist2 = X * X + Y * Y

    W0 = np.exp(-dist2 / (4.0 * epsilon))
    W0[0, 0] = 0.0  # diagonal removed
    return W0


def graph_laplacian_eigenvalues_fft(n_side: int, epsilon: float) -> np.ndarray:
    """
    Eigenvalues of L = D - W using FFT.

    For a circulant convolution matrix W:
        eigenvalues(W) = fft2(first row kernel)
        eigenvalues(L) = degree - eigenvalues(W)
    """
    W0 = torus_kernel_first_row(n_side, epsilon)
    degree = float(W0.sum())

    W_hat = np.fft.fft2(W0).real
    L_eigs = degree - W_hat

    eigvals = np.sort(L_eigs.ravel())
    return eigvals


def continuum_torus_eigenvalues(num_modes: int, max_k: int = 12) -> np.ndarray:
    vals = []
    for m in range(-max_k, max_k + 1):
        for n in range(-max_k, max_k + 1):
            vals.append(float(m * m + n * n))
    vals = np.array(sorted(vals), dtype=float)
    return vals[:num_modes]


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


def run_single(n_side: int, epsilon: float, num_modes: int = 30) -> dict:
    N = n_side * n_side

    eigvals = graph_laplacian_eigenvalues_fft(n_side, epsilon)
    target = continuum_torus_eigenvalues(num_modes)

    # First non-zero band has multiplicity 4:
    # (±1,0), (0,±1), eigenvalue 1.
    first_band_mean = float(np.mean(eigvals[1:5]))

    c_empirical = 1.0 / first_band_mean
    c_theory = np.pi / (N * epsilon * epsilon)

    norm_emp = c_empirical * eigvals[:num_modes]
    norm_theory = c_theory * eigvals[:num_modes]

    err_emp = spectral_error(norm_emp, target, start=1)
    err_theory = spectral_error(norm_theory, target, start=1)

    scaling_core = N * epsilon * epsilon
    pref_emp = c_empirical * scaling_core
    pref_theory = c_theory * scaling_core

    row = {
        "n_side": n_side,
        "N": N,
        "epsilon": epsilon,
        "num_modes": num_modes,
        "dx": float(2.0 * np.pi / n_side),
        "sqrt_epsilon_over_dx": float(np.sqrt(epsilon) / (2.0 * np.pi / n_side)),
        "lambda0_graph": float(eigvals[0]),
        "lambda1_graph": float(eigvals[1]),
        "lambda4_graph": float(eigvals[4]),
        "lambda1_band_mean": first_band_mean,
        "c_empirical": float(c_empirical),
        "c_theory_pi_over_N_eps2": float(c_theory),
        "scaling_core_N_eps2": float(scaling_core),
        "prefactor_empirical_c_times_N_eps2": float(pref_emp),
        "prefactor_theory_pi": float(pref_theory),
        "prefactor_rel_error_to_pi": float(abs(pref_emp - np.pi) / np.pi),
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
    n_sides = [32, 48, 64, 96, 128]
    epsilons = [0.01, 0.02, 0.04, 0.08, 0.12]
    num_modes = 30

    rows = []

    print("=" * 80)
    print("Paper 26 — Torus FFT spectrum convergence v2")
    print("=" * 80)

    for n_side in n_sides:
        dx = 2.0 * np.pi / n_side

        for epsilon in epsilons:
            # Need the kernel to be resolved by the grid.
            # sqrt(epsilon)/dx should not be too small.
            if np.sqrt(epsilon) / dx < 1.2:
                continue

            row = run_single(n_side=n_side, epsilon=epsilon, num_modes=num_modes)
            rows.append(row)

            print(
                f"n_side={n_side:3d} N={row['N']:5d} eps={epsilon:6.3f} "
                f"sqrt_eps/dx={row['sqrt_epsilon_over_dx']:.2f} "
                f"pref_emp={row['prefactor_empirical_c_times_N_eps2']:.6f} "
                f"pi={np.pi:.6f} "
                f"pref_rel={row['prefactor_rel_error_to_pi']:.3e} "
                f"mean_rel_theory={row['mean_rel_error_theory_c']:.3e} "
                f"max_rel_theory={row['max_rel_error_theory_c']:.3e}"
            )

    df = pd.DataFrame(rows)

    spectrum_csv = OUTDIR / "torus_fft_spectrum_convergence_v2.csv"
    prefactor_csv = OUTDIR / "torus_fft_prefactor_theory_v2.csv"
    summary_json = OUTDIR / "torus_fft_summary_v2.json"

    df.to_csv(spectrum_csv, index=False)

    pref_cols = [
        "n_side",
        "N",
        "epsilon",
        "dx",
        "sqrt_epsilon_over_dx",
        "c_empirical",
        "c_theory_pi_over_N_eps2",
        "prefactor_empirical_c_times_N_eps2",
        "prefactor_theory_pi",
        "prefactor_rel_error_to_pi",
        "mean_rel_error_empirical_c",
        "max_rel_error_empirical_c",
        "mean_rel_error_theory_c",
        "max_rel_error_theory_c",
    ]
    df[pref_cols].to_csv(prefactor_csv, index=False)

    best_pref = df.loc[df["prefactor_rel_error_to_pi"].idxmin()].to_dict()
    best_theory = df.loc[df["mean_rel_error_theory_c"].idxmin()].to_dict()

    summary = {
        "experiment": "paper26_torus_fft_spectrum_convergence_v2",
        "geometry": "flat square torus T^2 = [0,2pi)^2",
        "kernel": "periodic Gaussian exp(-d^2/(4 epsilon))",
        "laplacian": "unnormalized L = D - W, diagonal removed",
        "method": "FFT eigenvalues on regular periodic grid",
        "corrected_expected_scaling": "c_{N,epsilon} = pi/(N epsilon^2)",
        "theory_prefactor": float(np.pi),
        "num_rows": int(len(df)),
        "n_sides": n_sides,
        "epsilons": epsilons,
        "num_modes": num_modes,
        "best_by_prefactor_error": {
            "n_side": int(best_pref["n_side"]),
            "N": int(best_pref["N"]),
            "epsilon": float(best_pref["epsilon"]),
            "sqrt_epsilon_over_dx": float(best_pref["sqrt_epsilon_over_dx"]),
            "empirical_prefactor": float(best_pref["prefactor_empirical_c_times_N_eps2"]),
            "relative_error_to_pi": float(best_pref["prefactor_rel_error_to_pi"]),
            "mean_rel_error_theory_c": float(best_pref["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_pref["max_rel_error_theory_c"]),
        },
        "best_by_theory_spectral_error": {
            "n_side": int(best_theory["n_side"]),
            "N": int(best_theory["N"]),
            "epsilon": float(best_theory["epsilon"]),
            "sqrt_epsilon_over_dx": float(best_theory["sqrt_epsilon_over_dx"]),
            "empirical_prefactor": float(best_theory["prefactor_empirical_c_times_N_eps2"]),
            "relative_error_to_pi": float(best_theory["prefactor_rel_error_to_pi"]),
            "mean_rel_error_theory_c": float(best_theory["mean_rel_error_theory_c"]),
            "max_rel_error_theory_c": float(best_theory["max_rel_error_theory_c"]),
        },
        "outputs": {
            "torus_fft_spectrum_convergence_csv": str(spectrum_csv),
            "torus_fft_prefactor_theory_csv": str(prefactor_csv),
            "torus_fft_summary_json": str(summary_json),
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

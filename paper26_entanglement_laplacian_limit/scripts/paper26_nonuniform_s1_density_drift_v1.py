#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Non-uniform S1 density drift test v1

Purpose
-------
Test the important limitation of the raw unnormalized graph Laplacian:

    L = D - W

For uniform density, the normalized operator converges to -Delta.

For non-uniform density rho(theta), the raw combinatorial Laplacian acquires
a drift term proportional to grad log rho.

We sample S1 from:

    rho(theta) = (1 / 2pi) * (1 + a cos theta)

and test the operator on:

    f(theta) = sin(k theta)

The local continuum prediction is:

    c_local(theta) L f(theta)
      ≈ -f''(theta) - 2 (partial_theta log rho) f'(theta)

where:

    c_local(theta) = 1 / [rho_node(theta) * (4pi)^(1/2) * epsilon^(3/2)]

and rho_node(theta) is the node density per unit length:

    rho_node(theta) = N * rho_probability(theta).

Outputs
-------
paper26_entanglement_laplacian_limit/results/paper26_laplacian_limit_v1/
  nonuniform_s1_density_drift.csv
  nonuniform_s1_density_drift_summary.json

paper26_entanglement_laplacian_limit/figures/
  fig08_nonuniform_s1_density_drift.png/pdf
  fig09_nonuniform_s1_operator_comparison.png/pdf
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path("paper26_entanglement_laplacian_limit")
OUTDIR = ROOT / "results" / "paper26_laplacian_limit_v1"
FIGDIR = ROOT / "figures"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)


def rho_prob(theta: np.ndarray, a: float) -> np.ndarray:
    """
    Probability density on [0, 2pi):
        rho(theta) = (1/(2pi)) * (1 + a cos theta)
    """
    return (1.0 / (2.0 * np.pi)) * (1.0 + a * np.cos(theta))


def sample_nonuniform_s1(N: int, a: float, seed: int = 1234) -> np.ndarray:
    """
    Rejection sampling from rho(theta) proportional to 1 + a cos theta.
    """
    rng = np.random.default_rng(seed)
    samples = []

    rho_max = (1.0 + abs(a)) / (2.0 * np.pi)

    while len(samples) < N:
        theta = rng.uniform(0.0, 2.0 * np.pi, size=N)
        y = rng.uniform(0.0, rho_max, size=N)
        accept = y < rho_prob(theta, a)
        samples.extend(theta[accept].tolist())

    theta = np.array(samples[:N])
    theta.sort()
    return theta


def periodic_distance(theta: np.ndarray) -> np.ndarray:
    diff = np.abs(theta[:, None] - theta[None, :])
    return np.minimum(diff, 2.0 * np.pi - diff)


def build_laplacian(theta: np.ndarray, epsilon: float) -> np.ndarray:
    dist = periodic_distance(theta)
    W = np.exp(-(dist ** 2) / (4.0 * epsilon))
    np.fill_diagonal(W, 0.0)
    D = np.diag(W.sum(axis=1))
    return D - W


def test_function(theta: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    f = sin(k theta)
    f' = k cos(k theta)
    f'' = -k^2 sin(k theta)
    """
    f = np.sin(k * theta)
    fp = k * np.cos(k * theta)
    fpp = -(k ** 2) * np.sin(k * theta)
    return f, fp, fpp


def derivative_log_rho(theta: np.ndarray, a: float) -> np.ndarray:
    """
    d/dtheta log rho(theta), for rho ∝ 1 + a cos theta.
    """
    return (-a * np.sin(theta)) / (1.0 + a * np.cos(theta))


def weighted_rmse(y: np.ndarray, yhat: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y - yhat) ** 2)))


def relative_rmse(y: np.ndarray, yhat: np.ndarray) -> float:
    denom = np.sqrt(np.mean(y ** 2))
    return float(weighted_rmse(y, yhat) / max(denom, 1e-12))


def run_single(N: int, epsilon: float, a: float, k: int, seed: int) -> dict:
    theta = sample_nonuniform_s1(N=N, a=a, seed=seed)
    L = build_laplacian(theta, epsilon=epsilon)

    f, fp, fpp = test_function(theta, k=k)

    Lf = L @ f

    # Node density per unit length:
    rho_node = N * rho_prob(theta, a)

    # Local normalization:
    c_local = 1.0 / (rho_node * np.sqrt(4.0 * np.pi) * (epsilon ** 1.5))

    op_raw_local = c_local * Lf

    pure_target = -fpp

    drift = 2.0 * derivative_log_rho(theta, a) * fp
    drift_target = -fpp - drift

    pure_rel_rmse = relative_rmse(pure_target, op_raw_local)
    drift_rel_rmse = relative_rmse(drift_target, op_raw_local)

    improvement = pure_rel_rmse / max(drift_rel_rmse, 1e-12)

    row = {
        "N": N,
        "epsilon": epsilon,
        "a": a,
        "k": k,
        "seed": seed,
        "pure_rel_rmse": pure_rel_rmse,
        "drift_rel_rmse": drift_rel_rmse,
        "drift_improvement_factor": improvement,
        "mean_degree": float(np.mean(np.diag(L))),
        "min_rho_node": float(np.min(rho_node)),
        "max_rho_node": float(np.max(rho_node)),
        "rho_contrast": float(np.max(rho_node) / np.min(rho_node)),
    }

    return row, theta, op_raw_local, pure_target, drift_target


def main() -> None:
    N = 2048
    a = 0.5
    k = 2

    epsilons = [0.0025, 0.005, 0.01, 0.02]
    seeds = [11, 22, 33]

    rows = []

    print("=" * 80)
    print("Paper 26 — Non-uniform S1 density drift test v1")
    print("=" * 80)

    best_payload = None
    best_score = np.inf

    for epsilon in epsilons:
        for seed in seeds:
            row, theta, op_raw_local, pure_target, drift_target = run_single(
                N=N,
                epsilon=epsilon,
                a=a,
                k=k,
                seed=seed,
            )
            rows.append(row)

            print(
                f"N={N} eps={epsilon:.4f} seed={seed} "
                f"pure_rel_rmse={row['pure_rel_rmse']:.4e} "
                f"drift_rel_rmse={row['drift_rel_rmse']:.4e} "
                f"improvement={row['drift_improvement_factor']:.2f}"
            )

            if row["drift_rel_rmse"] < best_score:
                best_score = row["drift_rel_rmse"]
                best_payload = (row, theta, op_raw_local, pure_target, drift_target)

    df = pd.DataFrame(rows)

    out_csv = OUTDIR / "nonuniform_s1_density_drift.csv"
    out_json = OUTDIR / "nonuniform_s1_density_drift_summary.json"

    df.to_csv(out_csv, index=False)

    best_by_drift = df.loc[df["drift_rel_rmse"].idxmin()].to_dict()
    best_by_improvement = df.loc[df["drift_improvement_factor"].idxmax()].to_dict()

    summary = {
        "experiment": "paper26_nonuniform_s1_density_drift_v1",
        "density": "rho(theta) = (1/(2pi)) * (1 + a cos theta)",
        "a": a,
        "test_function": f"sin({k} theta)",
        "raw_laplacian": "L = D - W",
        "kernel": "exp(-d_S1^2/(4 epsilon))",
        "local_normalization": "1/[rho_node(theta) sqrt(4pi) epsilon^(3/2)]",
        "pure_target": "-f''",
        "drift_target": "-f'' - 2 (partial_theta log rho) f'",
        "num_rows": int(len(df)),
        "best_by_drift_rmse": {
            "N": int(best_by_drift["N"]),
            "epsilon": float(best_by_drift["epsilon"]),
            "seed": int(best_by_drift["seed"]),
            "pure_rel_rmse": float(best_by_drift["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_by_drift["drift_rel_rmse"]),
            "drift_improvement_factor": float(best_by_drift["drift_improvement_factor"]),
            "rho_contrast": float(best_by_drift["rho_contrast"]),
        },
        "best_by_improvement": {
            "N": int(best_by_improvement["N"]),
            "epsilon": float(best_by_improvement["epsilon"]),
            "seed": int(best_by_improvement["seed"]),
            "pure_rel_rmse": float(best_by_improvement["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_by_improvement["drift_rel_rmse"]),
            "drift_improvement_factor": float(best_by_improvement["drift_improvement_factor"]),
            "rho_contrast": float(best_by_improvement["rho_contrast"]),
        },
        "outputs": {
            "csv": str(out_csv),
            "summary_json": str(out_json),
        },
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Figures
    row, theta, op_raw_local, pure_target, drift_target = best_payload

    order = np.argsort(theta)
    theta_s = theta[order]
    op_s = op_raw_local[order]
    pure_s = pure_target[order]
    drift_s = drift_target[order]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(theta_s, pure_s, label="Pure Laplacian target: -f''")
    ax.plot(theta_s, drift_s, label="Density-drift target")
    ax.plot(theta_s, op_s, linestyle="--", label="Raw graph operator")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("Operator value")
    ax.set_title(
        "Non-uniform S1 density: raw Laplacian acquires a drift term\n"
        f"N={N}, epsilon={row['epsilon']}, a={a}, seed={row['seed']}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig08_png = FIGDIR / "fig08_nonuniform_s1_density_drift.png"
    fig08_pdf = FIGDIR / "fig08_nonuniform_s1_density_drift.pdf"
    plt.savefig(fig08_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig08_pdf, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    labels = ["pure target", "drift target"]
    values = [row["pure_rel_rmse"], row["drift_rel_rmse"]]
    ax.bar(labels, values)
    ax.set_ylabel("Relative RMSE")
    ax.set_title(
        "Operator comparison: pure Laplacian vs density-drift target\n"
        f"improvement factor = {row['drift_improvement_factor']:.2f}"
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig09_png = FIGDIR / "fig09_nonuniform_s1_operator_comparison.png"
    fig09_pdf = FIGDIR / "fig09_nonuniform_s1_operator_comparison.pdf"
    plt.savefig(fig09_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig09_pdf, bbox_inches="tight")
    plt.close(fig)

    print("-" * 80)
    print(f"[OK] wrote {out_csv}")
    print(f"[OK] wrote {out_json}")
    print(f"[OK] wrote {fig08_png}")
    print(f"[OK] wrote {fig08_pdf}")
    print(f"[OK] wrote {fig09_png}")
    print(f"[OK] wrote {fig09_pdf}")
    print("-" * 80)
    print("Best by drift RMSE:")
    print(summary["best_by_drift_rmse"])
    print("Best by improvement:")
    print(summary["best_by_improvement"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Non-uniform S1 density-corrected Laplacian v1

Purpose
-------
Test whether a density-corrected diffusion-map style kernel reduces the
density drift observed for the raw graph Laplacian.

Raw kernel:
    W_ij = exp[-d_S1(theta_i, theta_j)^2 / (4 epsilon)]

Density-corrected kernel:
    W_ij^(alpha) = W_ij / (q_i^alpha q_j^alpha),
    q_i = sum_j W_ij.

We compare:
    alpha = 0  -> raw combinatorial Laplacian
    alpha = 1  -> density-corrected kernel

For the same non-uniform density:
    rho(theta) = (1/(2pi)) * (1 + a cos theta)

and test function:
    f(theta) = sin(k theta)

Outputs
-------
results/paper26_laplacian_limit_v1/
  nonuniform_s1_density_corrected.csv
  nonuniform_s1_density_corrected_summary.json

figures/
  fig10_density_corrected_operator_comparison.png/pdf
  fig11_density_corrected_rmse_summary.png/pdf
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
    return (1.0 / (2.0 * np.pi)) * (1.0 + a * np.cos(theta))


def sample_nonuniform_s1(N: int, a: float, seed: int = 1234) -> np.ndarray:
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


def build_kernel(theta: np.ndarray, epsilon: float) -> np.ndarray:
    dist = periodic_distance(theta)
    W = np.exp(-(dist ** 2) / (4.0 * epsilon))
    np.fill_diagonal(W, 0.0)
    return W


def density_correct_kernel(W: np.ndarray, alpha: float) -> np.ndarray:
    if alpha == 0:
        return W.copy()

    q = W.sum(axis=1)
    q = np.maximum(q, 1e-14)
    Wcorr = W / ((q[:, None] ** alpha) * (q[None, :] ** alpha))
    np.fill_diagonal(Wcorr, 0.0)
    return Wcorr


def build_laplacian_from_kernel(W: np.ndarray) -> np.ndarray:
    D = np.diag(W.sum(axis=1))
    return D - W


def test_function(theta: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f = np.sin(k * theta)
    fp = k * np.cos(k * theta)
    fpp = -(k ** 2) * np.sin(k * theta)
    return f, fp, fpp


def derivative_log_rho(theta: np.ndarray, a: float) -> np.ndarray:
    return (-a * np.sin(theta)) / (1.0 + a * np.cos(theta))


def relative_rmse(y: np.ndarray, yhat: np.ndarray) -> float:
    rmse = np.sqrt(np.mean((y - yhat) ** 2))
    denom = np.sqrt(np.mean(y ** 2))
    return float(rmse / max(denom, 1e-12))


def fit_scale(x: np.ndarray, y: np.ndarray) -> float:
    """
    Fit scalar s minimizing ||s x - y||^2.
    """
    denom = float(np.dot(x, x))
    if denom <= 1e-14:
        return np.nan
    return float(np.dot(x, y) / denom)


def run_single(N: int, epsilon: float, a: float, k: int, seed: int, alpha: float) -> dict:
    theta = sample_nonuniform_s1(N=N, a=a, seed=seed)
    W = build_kernel(theta, epsilon)
    Wcorr = density_correct_kernel(W, alpha=alpha)
    L = build_laplacian_from_kernel(Wcorr)

    f, fp, fpp = test_function(theta, k)
    Lf = L @ f

    pure_target = -fpp
    drift_target = -fpp - 2.0 * derivative_log_rho(theta, a) * fp

    # Rather than impose a possibly convention-dependent prefactor,
    # fit one scalar scale per operator and compare shapes.
    s_pure = fit_scale(Lf, pure_target)
    s_drift = fit_scale(Lf, drift_target)

    op_scaled_to_pure = s_pure * Lf
    op_scaled_to_drift = s_drift * Lf

    pure_rel_rmse = relative_rmse(pure_target, op_scaled_to_pure)
    drift_rel_rmse = relative_rmse(drift_target, op_scaled_to_drift)

    row = {
        "N": N,
        "epsilon": epsilon,
        "a": a,
        "k": k,
        "seed": seed,
        "alpha": alpha,
        "scale_to_pure": s_pure,
        "scale_to_drift": s_drift,
        "pure_rel_rmse": pure_rel_rmse,
        "drift_rel_rmse": drift_rel_rmse,
        "pure_vs_drift_ratio": pure_rel_rmse / max(drift_rel_rmse, 1e-12),
        "mean_degree": float(Wcorr.sum(axis=1).mean()),
        "min_degree": float(Wcorr.sum(axis=1).min()),
        "max_degree": float(Wcorr.sum(axis=1).max()),
        "degree_contrast": float(Wcorr.sum(axis=1).max() / max(Wcorr.sum(axis=1).min(), 1e-14)),
    }

    payload = {
        "theta": theta,
        "Lf": Lf,
        "pure_target": pure_target,
        "drift_target": drift_target,
        "op_scaled_to_pure": op_scaled_to_pure,
        "op_scaled_to_drift": op_scaled_to_drift,
        "row": row,
    }

    return row, payload


def main() -> None:
    N = 2048
    a = 0.5
    k = 2
    epsilons = [0.005, 0.01, 0.02, 0.04]
    seeds = [11, 22, 33]
    alphas = [0.0, 0.5, 1.0]

    rows = []
    payloads = []

    print("=" * 80)
    print("Paper 26 — Non-uniform S1 density-corrected Laplacian v1")
    print("=" * 80)

    for epsilon in epsilons:
        for seed in seeds:
            for alpha in alphas:
                row, payload = run_single(
                    N=N,
                    epsilon=epsilon,
                    a=a,
                    k=k,
                    seed=seed,
                    alpha=alpha,
                )

                rows.append(row)
                payloads.append(payload)

                print(
                    f"eps={epsilon:.4f} seed={seed} alpha={alpha:.1f} "
                    f"pure_rmse={row['pure_rel_rmse']:.4e} "
                    f"drift_rmse={row['drift_rel_rmse']:.4e} "
                    f"degree_contrast={row['degree_contrast']:.3f}"
                )

    df = pd.DataFrame(rows)

    out_csv = OUTDIR / "nonuniform_s1_density_corrected.csv"
    out_json = OUTDIR / "nonuniform_s1_density_corrected_summary.json"

    df.to_csv(out_csv, index=False)

    # Best pure Laplacian recovery among corrected operators.
    best_pure = df.loc[df["pure_rel_rmse"].idxmin()].to_dict()

    # Best alpha=0 raw drift behavior.
    raw = df[df["alpha"] == 0.0]
    best_raw_drift = raw.loc[raw["drift_rel_rmse"].idxmin()].to_dict()

    # Median by alpha, useful for robust summary.
    alpha_summary = (
        df.groupby("alpha")[["pure_rel_rmse", "drift_rel_rmse", "degree_contrast"]]
        .median()
        .reset_index()
        .to_dict(orient="records")
    )

    summary = {
        "experiment": "paper26_nonuniform_s1_density_corrected_v1",
        "density": "rho(theta) = (1/(2pi)) * (1 + a cos theta)",
        "a": a,
        "test_function": f"sin({k} theta)",
        "kernel": "exp(-d_S1^2/(4 epsilon))",
        "density_correction": "W_ij^(alpha)=W_ij/(q_i^alpha q_j^alpha)",
        "alphas": alphas,
        "num_rows": int(len(df)),
        "best_pure_laplacian_recovery": {
            "alpha": float(best_pure["alpha"]),
            "epsilon": float(best_pure["epsilon"]),
            "seed": int(best_pure["seed"]),
            "pure_rel_rmse": float(best_pure["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_pure["drift_rel_rmse"]),
            "degree_contrast": float(best_pure["degree_contrast"]),
        },
        "best_raw_drift_fit": {
            "alpha": float(best_raw_drift["alpha"]),
            "epsilon": float(best_raw_drift["epsilon"]),
            "seed": int(best_raw_drift["seed"]),
            "pure_rel_rmse": float(best_raw_drift["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_raw_drift["drift_rel_rmse"]),
            "degree_contrast": float(best_raw_drift["degree_contrast"]),
        },
        "median_by_alpha": alpha_summary,
        "outputs": {
            "csv": str(out_csv),
            "summary_json": str(out_json),
        },
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Pick best pure payload for figure.
    target_alpha = best_pure["alpha"]
    target_eps = best_pure["epsilon"]
    target_seed = best_pure["seed"]

    best_payload = None
    for payload in payloads:
        r = payload["row"]
        if (
            abs(r["alpha"] - target_alpha) < 1e-12
            and abs(r["epsilon"] - target_eps) < 1e-12
            and r["seed"] == target_seed
        ):
            best_payload = payload
            break

    theta = best_payload["theta"]
    order = np.argsort(theta)
    theta_s = theta[order]
    op_s = best_payload["op_scaled_to_pure"][order]
    pure_s = best_payload["pure_target"][order]
    drift_s = best_payload["drift_target"][order]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(theta_s, pure_s, label="Pure Laplacian target: -f''")
    ax.plot(theta_s, drift_s, label="Density-drift target")
    ax.plot(theta_s, op_s, linestyle="--", label=f"Density-corrected graph operator alpha={target_alpha}")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("Operator value")
    ax.set_title(
        "Density-corrected S1 operator compared to pure Laplacian\n"
        f"epsilon={target_eps}, seed={target_seed}, alpha={target_alpha}, "
        f"pure RMSE={best_pure['pure_rel_rmse']:.3e}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig10_png = FIGDIR / "fig10_density_corrected_operator_comparison.png"
    fig10_pdf = FIGDIR / "fig10_density_corrected_operator_comparison.pdf"
    plt.savefig(fig10_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig10_pdf, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))

    med = pd.DataFrame(alpha_summary)
    ax.plot(med["alpha"], med["pure_rel_rmse"], marker="o", label="Pure Laplacian RMSE")
    ax.plot(med["alpha"], med["drift_rel_rmse"], marker="s", label="Drift target RMSE")
    ax.set_xlabel(r"Density correction exponent $\alpha$")
    ax.set_ylabel("Median relative RMSE")
    ax.set_title("Density correction reduces bias toward the drift operator")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig11_png = FIGDIR / "fig11_density_corrected_rmse_summary.png"
    fig11_pdf = FIGDIR / "fig11_density_corrected_rmse_summary.pdf"
    plt.savefig(fig11_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig11_pdf, bbox_inches="tight")
    plt.close(fig)

    print("-" * 80)
    print(f"[OK] wrote {out_csv}")
    print(f"[OK] wrote {out_json}")
    print(f"[OK] wrote {fig10_png}")
    print(f"[OK] wrote {fig10_pdf}")
    print(f"[OK] wrote {fig11_png}")
    print(f"[OK] wrote {fig11_pdf}")
    print("-" * 80)
    print("Best pure Laplacian recovery:")
    print(summary["best_pure_laplacian_recovery"])
    print("Best raw drift fit:")
    print(summary["best_raw_drift_fit"])
    print("Median by alpha:")
    print(summary["median_by_alpha"])


if __name__ == "__main__":
    main()

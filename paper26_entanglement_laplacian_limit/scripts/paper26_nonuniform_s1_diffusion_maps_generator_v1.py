#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Non-uniform S1 diffusion-maps generator test v1

Purpose
-------
Separate two conventions:

1. Symmetric combinatorial corrected Laplacian:
      L = D - W_alpha

2. Diffusion-maps Markov generator:
      K_alpha = W / (q_i^alpha q_j^alpha)
      P_alpha = row_normalize(K_alpha)
      G_alpha = (I - P_alpha) / epsilon

The previous test showed that alpha=1/2 worked best for the symmetric
combinatorial construction. This script tests whether alpha=1 behaves better
for the Markov diffusion-maps generator.

Density:
    rho(theta) = (1/(2pi)) * (1 + a cos theta), a=0.5

Test function:
    f(theta) = sin(2 theta)

Targets:
    pure  = -f''
    drift = -f'' - 2 (partial_theta log rho) f'

Outputs
-------
results/paper26_laplacian_limit_v1/
  nonuniform_s1_diffusion_maps_generator.csv
  nonuniform_s1_diffusion_maps_generator_summary.json

figures/
  fig12_diffusion_maps_alpha_summary.png/pdf
  fig13_diffusion_maps_best_operator.png/pdf
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


def sample_nonuniform_s1(N: int, a: float, seed: int) -> np.ndarray:
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


def diffusion_maps_operator(W: np.ndarray, epsilon: float, alpha: float) -> np.ndarray:
    """
    Build G_alpha = (I - P_alpha) / epsilon.

    q_i = sum_j W_ij
    K_alpha_ij = W_ij / (q_i^alpha q_j^alpha)
    P_alpha = row-normalized K_alpha
    """
    q = W.sum(axis=1)
    q = np.maximum(q, 1e-14)

    K = W / ((q[:, None] ** alpha) * (q[None, :] ** alpha))
    np.fill_diagonal(K, 0.0)

    d = K.sum(axis=1)
    d = np.maximum(d, 1e-14)

    P = K / d[:, None]

    N = W.shape[0]
    G = (np.eye(N) - P) / epsilon

    return G


def test_function(theta: np.ndarray, k: int):
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
    denom = float(np.dot(x, x))
    if denom <= 1e-14:
        return np.nan
    return float(np.dot(x, y) / denom)


def run_single(N: int, epsilon: float, a: float, k: int, seed: int, alpha: float):
    theta = sample_nonuniform_s1(N=N, a=a, seed=seed)
    W = build_kernel(theta, epsilon)

    G = diffusion_maps_operator(W, epsilon=epsilon, alpha=alpha)

    f, fp, fpp = test_function(theta, k=k)

    Gf = G @ f

    pure_target = -fpp
    drift_target = -fpp - 2.0 * derivative_log_rho(theta, a) * fp

    # Fit scalar scale because generator prefactors may differ by convention.
    s_pure = fit_scale(Gf, pure_target)
    s_drift = fit_scale(Gf, drift_target)

    op_scaled_pure = s_pure * Gf
    op_scaled_drift = s_drift * Gf

    pure_rmse = relative_rmse(pure_target, op_scaled_pure)
    drift_rmse = relative_rmse(drift_target, op_scaled_drift)

    q = W.sum(axis=1)

    row = {
        "N": N,
        "epsilon": epsilon,
        "a": a,
        "k": k,
        "seed": seed,
        "alpha": alpha,
        "pure_rel_rmse": pure_rmse,
        "drift_rel_rmse": drift_rmse,
        "pure_vs_drift_ratio": pure_rmse / max(drift_rmse, 1e-12),
        "scale_to_pure": s_pure,
        "scale_to_drift": s_drift,
        "q_min": float(q.min()),
        "q_max": float(q.max()),
        "q_contrast": float(q.max() / max(q.min(), 1e-14)),
    }

    payload = {
        "theta": theta,
        "op_scaled_pure": op_scaled_pure,
        "pure_target": pure_target,
        "drift_target": drift_target,
        "row": row,
    }

    return row, payload


def main():
    N = 2048
    a = 0.5
    k = 2

    epsilons = [0.005, 0.01, 0.02, 0.04, 0.08]
    seeds = [11, 22, 33]
    alphas = [0.0, 0.5, 1.0]

    rows = []
    payloads = []

    print("=" * 80)
    print("Paper 26 — Non-uniform S1 diffusion-maps generator test v1")
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
                    f"ratio={row['pure_vs_drift_ratio']:.2f}"
                )

    df = pd.DataFrame(rows)

    out_csv = OUTDIR / "nonuniform_s1_diffusion_maps_generator.csv"
    out_json = OUTDIR / "nonuniform_s1_diffusion_maps_generator_summary.json"

    df.to_csv(out_csv, index=False)

    best_pure = df.loc[df["pure_rel_rmse"].idxmin()].to_dict()
    best_drift = df.loc[df["drift_rel_rmse"].idxmin()].to_dict()

    alpha_summary = (
        df.groupby("alpha")[["pure_rel_rmse", "drift_rel_rmse", "pure_vs_drift_ratio"]]
        .median()
        .reset_index()
        .to_dict(orient="records")
    )

    summary = {
        "experiment": "paper26_nonuniform_s1_diffusion_maps_generator_v1",
        "operator": "G_alpha=(I-P_alpha)/epsilon, P_alpha=row_normalized(W/(q_i^alpha q_j^alpha))",
        "density": "rho(theta)=(1/(2pi))*(1+a cos theta)",
        "a": a,
        "test_function": f"sin({k} theta)",
        "alphas": alphas,
        "num_rows": int(len(df)),
        "best_pure_laplacian_recovery": {
            "alpha": float(best_pure["alpha"]),
            "epsilon": float(best_pure["epsilon"]),
            "seed": int(best_pure["seed"]),
            "pure_rel_rmse": float(best_pure["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_pure["drift_rel_rmse"]),
            "pure_vs_drift_ratio": float(best_pure["pure_vs_drift_ratio"]),
        },
        "best_drift_recovery": {
            "alpha": float(best_drift["alpha"]),
            "epsilon": float(best_drift["epsilon"]),
            "seed": int(best_drift["seed"]),
            "pure_rel_rmse": float(best_drift["pure_rel_rmse"]),
            "drift_rel_rmse": float(best_drift["drift_rel_rmse"]),
            "pure_vs_drift_ratio": float(best_drift["pure_vs_drift_ratio"]),
        },
        "median_by_alpha": alpha_summary,
        "outputs": {
            "csv": str(out_csv),
            "summary_json": str(out_json),
        },
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Figure 12: alpha summary.
    med = pd.DataFrame(alpha_summary)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(med["alpha"], med["pure_rel_rmse"], marker="o", label="Pure Laplacian RMSE")
    ax.plot(med["alpha"], med["drift_rel_rmse"], marker="s", label="Drift target RMSE")
    ax.set_xlabel(r"Diffusion maps exponent $\alpha$")
    ax.set_ylabel("Median relative RMSE")
    ax.set_title("Diffusion-maps Markov generator on non-uniform S1")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig12_png = FIGDIR / "fig12_diffusion_maps_alpha_summary.png"
    fig12_pdf = FIGDIR / "fig12_diffusion_maps_alpha_summary.pdf"
    plt.savefig(fig12_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig12_pdf, bbox_inches="tight")
    plt.close(fig)

    # Figure 13: best pure operator.
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
    op_s = best_payload["op_scaled_pure"][order]
    pure_s = best_payload["pure_target"][order]
    drift_s = best_payload["drift_target"][order]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(theta_s, pure_s, label="Pure Laplacian target: -f''")
    ax.plot(theta_s, drift_s, label="Density-drift target")
    ax.plot(theta_s, op_s, linestyle="--", label=f"Markov generator alpha={target_alpha}")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("Operator value")
    ax.set_title(
        "Best diffusion-maps Markov generator recovery\n"
        f"epsilon={target_eps}, seed={target_seed}, alpha={target_alpha}, "
        f"pure RMSE={best_pure['pure_rel_rmse']:.3e}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig13_png = FIGDIR / "fig13_diffusion_maps_best_operator.png"
    fig13_pdf = FIGDIR / "fig13_diffusion_maps_best_operator.pdf"
    plt.savefig(fig13_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig13_pdf, bbox_inches="tight")
    plt.close(fig)

    print("-" * 80)
    print(f"[OK] wrote {out_csv}")
    print(f"[OK] wrote {out_json}")
    print(f"[OK] wrote {fig12_png}")
    print(f"[OK] wrote {fig12_pdf}")
    print(f"[OK] wrote {fig13_png}")
    print(f"[OK] wrote {fig13_pdf}")
    print("-" * 80)
    print("Best pure Laplacian recovery:")
    print(summary["best_pure_laplacian_recovery"])
    print("Best drift recovery:")
    print(summary["best_drift_recovery"])
    print("Median by alpha:")
    print(summary["median_by_alpha"])


if __name__ == "__main__":
    main()

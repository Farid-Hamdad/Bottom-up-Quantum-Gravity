#!/usr/bin/env python3
"""
Paper 29 — Post-hoc quasi-invariant compensation analysis.

Goal
----
Test whether the standardized nonparametric residuals admit a quasi-invariant

    C_alpha = z(epsilon_mod) + alpha * z(epsilon_W)

with a common alpha across N8/A4 and N10/A4.

For each bandwidth:
- compute the variance-minimizing alpha separately for N8 and N10;
- compute the pooled alpha using both sizes together;
- evaluate variance suppression relative to Var[z(epsilon_mod)] = 1;
- test cross-size transfer using alpha learned on one size and applied to the other;
- compare alpha to the natural complementarity value alpha = 1.

Reads only frozen nonparametric residual CSVs.

Status: POST_HOC_EXPLORATORY
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "posthoc_analysis" / "nonparametric_coupling_output"
OUT = ROOT / "posthoc_analysis" / "quasi_invariant_output"

BANDWIDTHS = [0.05, 0.10, 0.20, 0.30, 0.40]


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def subset(rows, frac):
    return [
        r for r in rows
        if abs(float(r["bandwidth_fraction_Q_range"]) - frac) < 1e-12
    ]


def arr(rows, key):
    return np.asarray([float(r[key]) for r in rows], dtype=float)


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = float(np.std(x))
    if sd <= 0.0:
        raise RuntimeError("ZERO_STD")
    return (x - np.mean(x)) / sd


def optimal_alpha(z_w, z_m):
    """
    Minimize mean[(z_m + alpha z_w)^2].
    Since both are centered, alpha* = -Cov(z_w,z_m)/Var(z_w).
    """
    den = float(np.dot(z_w, z_w))
    if den == 0.0:
        raise RuntimeError("ZERO_DENOMINATOR")
    return float(-np.dot(z_w, z_m) / den)


def evaluate_alpha(z_w, z_m, alpha):
    c = z_m + alpha * z_w
    var_c = float(np.var(c))
    rms_c = float(np.sqrt(np.mean(c ** 2)))
    baseline_var = float(np.var(z_m))
    suppression = float(var_c / baseline_var) if baseline_var > 0.0 else float("nan")
    explained = float(1.0 - suppression) if np.isfinite(suppression) else float("nan")
    return {
        "alpha": float(alpha),
        "variance_C": var_c,
        "rms_C": rms_c,
        "variance_ratio_to_zmod": suppression,
        "variance_suppression_fraction": explained,
        "max_abs_C": float(np.max(np.abs(c))),
        "mean_C": float(np.mean(c)),
        "C": c,
    }


def prepare(rows_all, frac):
    rows = subset(rows_all, frac)
    if len(rows) != 30:
        raise RuntimeError(f"UNEXPECTED_ROW_COUNT_{frac}_{len(rows)}")
    ew = arr(rows, "resid_dW")
    em = arr(rows, "resid_dmod")
    return zscore(ew), zscore(em)


def stripped(d):
    return {k: v for k, v in d.items() if k != "C"}


def plot_bandwidth(frac, e8, e10, pooled_alpha):
    c8 = e8["zM"] + pooled_alpha * e8["zW"]
    c10 = e10["zM"] + pooled_alpha * e10["zW"]

    plt.figure(figsize=(8, 5))
    plt.plot(np.arange(len(c8)), c8, marker="o", label="N8/A4")
    plt.plot(np.arange(len(c10)), c10, marker="o", label="N10/A4")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("transition index")
    plt.ylabel("C_alpha")
    plt.title(f"Paper 29 post-hoc quasi-invariant — h={frac:.2f} ΔQ, alpha={pooled_alpha:.4f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"quasi_invariant_bw_{frac:.2f}.png", dpi=180)
    plt.close()


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    rows8_all = read_csv(SOURCE / "nonparametric_residuals_N8_A4.csv")
    rows10_all = read_csv(SOURCE / "nonparametric_residuals_N10_A4.csv")

    results = []

    for frac in BANDWIDTHS:
        zw8, zm8 = prepare(rows8_all, frac)
        zw10, zm10 = prepare(rows10_all, frac)

        alpha8 = optimal_alpha(zw8, zm8)
        alpha10 = optimal_alpha(zw10, zm10)

        zw_pool = np.concatenate([zw8, zw10])
        zm_pool = np.concatenate([zm8, zm10])
        alpha_pool = optimal_alpha(zw_pool, zm_pool)

        eval8_opt = evaluate_alpha(zw8, zm8, alpha8)
        eval10_opt = evaluate_alpha(zw10, zm10, alpha10)
        eval8_pool = evaluate_alpha(zw8, zm8, alpha_pool)
        eval10_pool = evaluate_alpha(zw10, zm10, alpha_pool)

        eval8_alpha1 = evaluate_alpha(zw8, zm8, 1.0)
        eval10_alpha1 = evaluate_alpha(zw10, zm10, 1.0)

        transfer_8_to_10 = evaluate_alpha(zw10, zm10, alpha8)
        transfer_10_to_8 = evaluate_alpha(zw8, zm8, alpha10)

        results.append({
            "bandwidth_fraction_Q_range": float(frac),
            "alpha_optimal": {
                "N8_A4": alpha8,
                "N10_A4": alpha10,
                "pooled": alpha_pool,
                "absolute_difference_N8_N10": float(abs(alpha8 - alpha10)),
                "relative_difference_to_mean": float(
                    abs(alpha8 - alpha10) / (0.5 * (abs(alpha8) + abs(alpha10)))
                ),
                "distance_pooled_from_one": float(abs(alpha_pool - 1.0)),
            },
            "within_size_optimal": {
                "N8_A4": stripped(eval8_opt),
                "N10_A4": stripped(eval10_opt),
            },
            "pooled_alpha_applied": {
                "N8_A4": stripped(eval8_pool),
                "N10_A4": stripped(eval10_pool),
            },
            "alpha_equal_one_control": {
                "N8_A4": stripped(eval8_alpha1),
                "N10_A4": stripped(eval10_alpha1),
            },
            "cross_size_transfer": {
                "alpha_N8_applied_to_N10": stripped(transfer_8_to_10),
                "alpha_N10_applied_to_N8": stripped(transfer_10_to_8),
            },
        })

        plot_bandwidth(
            frac,
            {"zW": zw8, "zM": zm8},
            {"zW": zw10, "zM": zm10},
            alpha_pool,
        )

    pooled_alphas = np.asarray(
        [r["alpha_optimal"]["pooled"] for r in results],
        dtype=float,
    )
    var8_pool = np.asarray(
        [r["pooled_alpha_applied"]["N8_A4"]["variance_ratio_to_zmod"] for r in results],
        dtype=float,
    )
    var10_pool = np.asarray(
        [r["pooled_alpha_applied"]["N10_A4"]["variance_ratio_to_zmod"] for r in results],
        dtype=float,
    )
    var8_one = np.asarray(
        [r["alpha_equal_one_control"]["N8_A4"]["variance_ratio_to_zmod"] for r in results],
        dtype=float,
    )
    var10_one = np.asarray(
        [r["alpha_equal_one_control"]["N10_A4"]["variance_ratio_to_zmod"] for r in results],
        dtype=float,
    )

    payload = {
        "analysis_status": "POST_HOC_EXPLORATORY",
        "question": (
            "Do standardized residuals admit a size-transferable quasi-invariant "
            "C_alpha = z(epsilon_mod) + alpha z(epsilon_W) with alpha near 1?"
        ),
        "source": (
            "Frozen nonparametric residual CSV files from "
            "posthoc_analysis/nonparametric_coupling_output."
        ),
        "bandwidth_fractions_Q_range": BANDWIDTHS,
        "bandwidth_results": results,
        "robustness_summary": {
            "mean_pooled_alpha": float(np.mean(pooled_alphas)),
            "std_pooled_alpha": float(np.std(pooled_alphas)),
            "max_distance_pooled_alpha_from_one": float(np.max(np.abs(pooled_alphas - 1.0))),
            "mean_variance_ratio_pooled_alpha_N8": float(np.mean(var8_pool)),
            "mean_variance_ratio_pooled_alpha_N10": float(np.mean(var10_pool)),
            "max_variance_ratio_pooled_alpha_N8": float(np.max(var8_pool)),
            "max_variance_ratio_pooled_alpha_N10": float(np.max(var10_pool)),
            "mean_variance_ratio_alpha_one_N8": float(np.mean(var8_one)),
            "mean_variance_ratio_alpha_one_N10": float(np.mean(var10_one)),
            "max_variance_ratio_alpha_one_N8": float(np.max(var8_one)),
            "max_variance_ratio_alpha_one_N10": float(np.max(var10_one)),
        },
        "interpretation_boundary": (
            "A low-variance C_alpha with alpha near 1 would support an exploratory "
            "quasi-invariant compensation signature in standardized residual coordinates. "
            "It would not establish a fundamental conservation law or literal information transfer."
        ),
    }

    (OUT / "quasi_invariant_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("PAPER29_POSTHOC_QUASI_INVARIANT_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

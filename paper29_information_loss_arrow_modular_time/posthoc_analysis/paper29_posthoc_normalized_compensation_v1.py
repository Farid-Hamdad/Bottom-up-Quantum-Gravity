#!/usr/bin/env python3
"""
Paper 29 — Post-hoc normalized compensation analysis.

Question
--------
After standardizing the nonparametric residuals within each system and bandwidth,
does the compensation relation approach a universal form

    z(epsilon_mod) = a + b * z(epsilon_W)

with b close to -1?

This script reads ONLY the already frozen nonparametric residual CSV files.
It does not recompute the smoothing stage and does not modify any frozen result.

Status: POST_HOC_EXPLORATORY
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "posthoc_analysis" / "nonparametric_coupling_output"
OUT = ROOT / "posthoc_analysis" / "normalized_compensation_output"

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
    mu = float(np.mean(x))
    sd = float(np.std(x))
    if sd <= 0.0:
        raise RuntimeError("ZERO_STD")
    return (x - mu) / sd, mu, sd


def fit_line(x, y):
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    resid = y - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else float("nan")
    return float(beta[0]), float(beta[1]), pred, r2


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])


def analyze_dataset(label, rows_all, frac):
    rows = subset(rows_all, frac)
    if len(rows) != 30:
        raise RuntimeError(f"UNEXPECTED_ROW_COUNT_{label}_{frac}_{len(rows)}")

    ew = arr(rows, "resid_dW")
    em = arr(rows, "resid_dmod")

    zw, mu_w, sd_w = zscore(ew)
    zm, mu_m, sd_m = zscore(em)

    intercept, slope, pred, r2 = fit_line(zw, zm)

    # With separate z-scoring, slope equals Pearson r for OLS with intercept.
    corr = pearson(zw, zm)

    return {
        "label": label,
        "n": int(len(rows)),
        "bandwidth_fraction_Q_range": float(frac),
        "mean_epsilon_W": mu_w,
        "std_epsilon_W": sd_w,
        "mean_epsilon_mod": mu_m,
        "std_epsilon_mod": sd_m,
        "scale_ratio_std_mod_over_std_W": float(sd_m / sd_w),
        "normalized_intercept": intercept,
        "normalized_slope": slope,
        "pearson_standardized": corr,
        "normalized_r2": r2,
        "distance_of_slope_from_minus_one": float(abs(slope + 1.0)),
        "zW": zw,
        "zM": zm,
        "pred": pred,
    }


def external_standardized_prediction(source, target):
    # Apply source's normalized law zM = a + b zW to target standardized coordinates.
    pred = source["normalized_intercept"] + source["normalized_slope"] * target["zW"]
    y = target["zM"]

    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else float("nan")
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))

    return {
        "r2_external": r2,
        "rmse_standardized": rmse,
        "pearson_observed_predicted": pearson(y, pred),
    }


def plot_bandwidth(frac, n8, n10):
    plt.figure(figsize=(7, 5))

    plt.scatter(n8["zW"], n8["zM"], label="N8/A4")
    plt.scatter(n10["zW"], n10["zM"], label="N10/A4")

    xmin = min(float(np.min(n8["zW"])), float(np.min(n10["zW"])))
    xmax = max(float(np.max(n8["zW"])), float(np.max(n10["zW"])))
    xx = np.linspace(xmin, xmax, 200)

    plt.plot(xx, n8["normalized_intercept"] + n8["normalized_slope"] * xx, label="fit N8")
    plt.plot(xx, n10["normalized_intercept"] + n10["normalized_slope"] * xx, label="fit N10")
    plt.plot(xx, -xx, linestyle="--", label="ideal y=-x")

    plt.axhline(0.0, linewidth=1)
    plt.axvline(0.0, linewidth=1)
    plt.xlabel("z(epsilon_W)")
    plt.ylabel("z(epsilon_mod)")
    plt.title(f"Paper 29 post-hoc normalized compensation — h={frac:.2f} ΔQ")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"normalized_compensation_bw_{frac:.2f}.png", dpi=180)
    plt.close()


def stripped(d):
    return {
        k: v for k, v in d.items()
        if k not in ("zW", "zM", "pred")
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    rows8_all = read_csv(SOURCE / "nonparametric_residuals_N8_A4.csv")
    rows10_all = read_csv(SOURCE / "nonparametric_residuals_N10_A4.csv")

    results = []

    for frac in BANDWIDTHS:
        n8 = analyze_dataset("N8_A4", rows8_all, frac)
        n10 = analyze_dataset("N10_A4", rows10_all, frac)

        transfer_8_to_10 = external_standardized_prediction(n8, n10)
        transfer_10_to_8 = external_standardized_prediction(n10, n8)

        slope_abs_diff = abs(n8["normalized_slope"] - n10["normalized_slope"])
        slope_mean_abs = 0.5 * (
            abs(n8["normalized_slope"]) + abs(n10["normalized_slope"])
        )

        results.append({
            "bandwidth_fraction_Q_range": float(frac),
            "N8_A4": stripped(n8),
            "N10_A4": stripped(n10),
            "slope_comparison": {
                "absolute_difference": float(slope_abs_diff),
                "relative_difference_to_mean_abs_slope": float(
                    slope_abs_diff / slope_mean_abs
                ) if slope_mean_abs > 0.0 else float("nan"),
                "both_negative": bool(
                    n8["normalized_slope"] < 0.0 and n10["normalized_slope"] < 0.0
                ),
                "both_within_0p2_of_minus_one": bool(
                    abs(n8["normalized_slope"] + 1.0) <= 0.2
                    and abs(n10["normalized_slope"] + 1.0) <= 0.2
                ),
            },
            "cross_size_transfer_standardized": {
                "fit_N8_predict_N10": transfer_8_to_10,
                "fit_N10_predict_N8": transfer_10_to_8,
            },
        })

        plot_bandwidth(frac, n8, n10)

    slopes8 = np.asarray([r["N8_A4"]["normalized_slope"] for r in results], dtype=float)
    slopes10 = np.asarray([r["N10_A4"]["normalized_slope"] for r in results], dtype=float)
    scale_ratios8 = np.asarray(
        [r["N8_A4"]["scale_ratio_std_mod_over_std_W"] for r in results], dtype=float
    )
    scale_ratios10 = np.asarray(
        [r["N10_A4"]["scale_ratio_std_mod_over_std_W"] for r in results], dtype=float
    )
    r2_8_to_10 = np.asarray(
        [r["cross_size_transfer_standardized"]["fit_N8_predict_N10"]["r2_external"]
         for r in results],
        dtype=float,
    )
    r2_10_to_8 = np.asarray(
        [r["cross_size_transfer_standardized"]["fit_N10_predict_N8"]["r2_external"]
         for r in results],
        dtype=float,
    )

    payload = {
        "analysis_status": "POST_HOC_EXPLORATORY",
        "question": (
            "After within-system standardization, does the residual compensation relation "
            "approach a universal z(epsilon_mod) ≈ -z(epsilon_W) law?"
        ),
        "source": (
            "Frozen nonparametric residual CSV files from "
            "posthoc_analysis/nonparametric_coupling_output."
        ),
        "bandwidth_fractions_Q_range": BANDWIDTHS,
        "bandwidth_results": results,
        "robustness_summary": {
            "mean_normalized_slope_N8": float(np.mean(slopes8)),
            "mean_normalized_slope_N10": float(np.mean(slopes10)),
            "std_normalized_slope_N8": float(np.std(slopes8)),
            "std_normalized_slope_N10": float(np.std(slopes10)),
            "max_distance_from_minus_one_N8": float(np.max(np.abs(slopes8 + 1.0))),
            "max_distance_from_minus_one_N10": float(np.max(np.abs(slopes10 + 1.0))),
            "mean_scale_ratio_std_mod_over_std_W_N8": float(np.mean(scale_ratios8)),
            "mean_scale_ratio_std_mod_over_std_W_N10": float(np.mean(scale_ratios10)),
            "mean_cross_R2_standardized_N8_to_N10": float(np.mean(r2_8_to_10)),
            "min_cross_R2_standardized_N8_to_N10": float(np.min(r2_8_to_10)),
            "mean_cross_R2_standardized_N10_to_N8": float(np.mean(r2_10_to_8)),
            "min_cross_R2_standardized_N10_to_N8": float(np.min(r2_10_to_8)),
            "negative_slope_count_N8": int(np.sum(slopes8 < 0.0)),
            "negative_slope_count_N10": int(np.sum(slopes10 < 0.0)),
        },
        "interpretation_boundary": (
            "A slope near -1 after standardization would support an exploratory universal "
            "complementarity pattern in standardized residual coordinates. It would not "
            "establish conservation of information or a fundamental transfer law."
        ),
    }

    (OUT / "normalized_compensation_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("PAPER29_POSTHOC_NORMALIZED_COMPENSATION_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Paper 29 — Post-hoc compensation-law analysis.

Purpose
-------
Test whether the nonparametric residuals support a quantitative compensation law

    epsilon_mod = a + b * epsilon_W

and whether that law transfers between N8/A4 and N10/A4.

This script reads ONLY the already frozen post-hoc nonparametric residual CSV files.
It does not recompute the nonparametric regressions and does not modify any frozen result.

Status: POST_HOC_EXPLORATORY
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "posthoc_analysis" / "nonparametric_coupling_output"
OUT = ROOT / "posthoc_analysis" / "compensation_law_output"

BANDWIDTHS = [0.05, 0.10, 0.20, 0.30, 0.40]


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def rows_for_bandwidth(rows, frac):
    out = []
    for r in rows:
        if abs(float(r["bandwidth_fraction_Q_range"]) - frac) < 1e-12:
            out.append(r)
    return out


def arr(rows, key):
    return np.asarray([float(r[key]) for r in rows], dtype=float)


def ols_with_intercept(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    resid = y - pred

    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    return {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
        "pred": pred,
        "resid": resid,
        "r2": r2,
        "rmse": float(np.sqrt(np.mean(resid ** 2))),
    }


def slope_through_origin(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    den = float(np.dot(x, x))
    if den == 0.0:
        return float("nan")
    return float(np.dot(x, y) / den)


def r2_external(y, pred):
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def normalized_rmse(y, pred):
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    sd = float(np.std(y))
    return float(rmse / sd) if sd > 0 else float("nan")


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])


def fit_dataset(rows):
    x = arr(rows, "resid_dW")
    y = arr(rows, "resid_dmod")
    fit = ols_with_intercept(x, y)

    return {
        "n": int(len(x)),
        "intercept": fit["intercept"],
        "slope": fit["slope"],
        "compensation_coefficient_c_minus_slope": float(-fit["slope"]),
        "slope_through_origin": slope_through_origin(x, y),
        "pearson": pearson(x, y),
        "r2": fit["r2"],
        "rmse": fit["rmse"],
        "std_resid_dW": float(np.std(x)),
        "std_resid_dmod": float(np.std(y)),
        "x": x,
        "y": y,
        "pred": fit["pred"],
    }


def external_prediction(source_fit, target_rows):
    x = arr(target_rows, "resid_dW")
    y = arr(target_rows, "resid_dmod")
    pred = source_fit["intercept"] + source_fit["slope"] * x

    return {
        "r2_external": r2_external(y, pred),
        "normalized_rmse": normalized_rmse(y, pred),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
        "pearson_observed_predicted": pearson(y, pred),
    }


def pooled_fit(rows8, rows10):
    x = np.concatenate([arr(rows8, "resid_dW"), arr(rows10, "resid_dW")])
    y = np.concatenate([arr(rows8, "resid_dmod"), arr(rows10, "resid_dmod")])
    fit = ols_with_intercept(x, y)
    return {
        "n": int(len(x)),
        "intercept": fit["intercept"],
        "slope": fit["slope"],
        "compensation_coefficient_c_minus_slope": float(-fit["slope"]),
        "r2": fit["r2"],
        "rmse": fit["rmse"],
        "pearson": pearson(x, y),
    }


def plot_bandwidth(frac, f8, f10):
    plt.figure(figsize=(7, 5))

    plt.scatter(f8["x"], f8["y"], label="N8/A4")
    plt.scatter(f10["x"], f10["y"], label="N10/A4")

    xmin = min(float(np.min(f8["x"])), float(np.min(f10["x"])))
    xmax = max(float(np.max(f8["x"])), float(np.max(f10["x"])))
    xx = np.linspace(xmin, xmax, 200)

    plt.plot(xx, f8["intercept"] + f8["slope"] * xx, label="fit N8")
    plt.plot(xx, f10["intercept"] + f10["slope"] * xx, label="fit N10")

    plt.axhline(0.0, linewidth=1)
    plt.axvline(0.0, linewidth=1)
    plt.xlabel("epsilon_W")
    plt.ylabel("epsilon_mod")
    plt.title(f"Paper 29 post-hoc compensation law — h={frac:.2f} ΔQ")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"compensation_law_bw_{frac:.2f}.png", dpi=180)
    plt.close()


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    rows8_all = read_csv(SOURCE / "nonparametric_residuals_N8_A4.csv")
    rows10_all = read_csv(SOURCE / "nonparametric_residuals_N10_A4.csv")

    bandwidth_results = []

    for frac in BANDWIDTHS:
        rows8 = rows_for_bandwidth(rows8_all, frac)
        rows10 = rows_for_bandwidth(rows10_all, frac)

        if len(rows8) != 30 or len(rows10) != 30:
            raise RuntimeError(
                f"UNEXPECTED_ROW_COUNT_BW_{frac}: N8={len(rows8)} N10={len(rows10)}"
            )

        f8 = fit_dataset(rows8)
        f10 = fit_dataset(rows10)

        n8_to_n10 = external_prediction(f8, rows10)
        n10_to_n8 = external_prediction(f10, rows8)
        pooled = pooled_fit(rows8, rows10)

        slope8 = f8["slope"]
        slope10 = f10["slope"]
        mean_abs_slope = 0.5 * (abs(slope8) + abs(slope10))

        slope_relative_difference = (
            abs(slope8 - slope10) / mean_abs_slope
            if mean_abs_slope > 0.0
            else float("nan")
        )

        bandwidth_results.append({
            "bandwidth_fraction_Q_range": float(frac),
            "N8_A4": {
                k: v for k, v in f8.items()
                if k not in ("x", "y", "pred")
            },
            "N10_A4": {
                k: v for k, v in f10.items()
                if k not in ("x", "y", "pred")
            },
            "slope_comparison": {
                "absolute_difference": float(abs(slope8 - slope10)),
                "relative_difference_to_mean_abs_slope": float(slope_relative_difference),
                "same_negative_sign": bool(slope8 < 0.0 and slope10 < 0.0),
                "ratio_N10_over_N8": float(slope10 / slope8) if slope8 != 0.0 else float("nan"),
            },
            "cross_size_transfer": {
                "fit_N8_predict_N10": n8_to_n10,
                "fit_N10_predict_N8": n10_to_n8,
            },
            "pooled_fit": pooled,
        })

        plot_bandwidth(frac, f8, f10)

    slope8_all = np.asarray(
        [r["N8_A4"]["slope"] for r in bandwidth_results], dtype=float
    )
    slope10_all = np.asarray(
        [r["N10_A4"]["slope"] for r in bandwidth_results], dtype=float
    )
    rel_diff_all = np.asarray(
        [r["slope_comparison"]["relative_difference_to_mean_abs_slope"]
         for r in bandwidth_results],
        dtype=float,
    )
    cross_r2_8to10 = np.asarray(
        [r["cross_size_transfer"]["fit_N8_predict_N10"]["r2_external"]
         for r in bandwidth_results],
        dtype=float,
    )
    cross_r2_10to8 = np.asarray(
        [r["cross_size_transfer"]["fit_N10_predict_N8"]["r2_external"]
         for r in bandwidth_results],
        dtype=float,
    )

    robustness = {
        "negative_slope_count_N8": int(np.sum(slope8_all < 0.0)),
        "negative_slope_count_N10": int(np.sum(slope10_all < 0.0)),
        "mean_slope_N8": float(np.mean(slope8_all)),
        "mean_slope_N10": float(np.mean(slope10_all)),
        "std_slope_N8": float(np.std(slope8_all)),
        "std_slope_N10": float(np.std(slope10_all)),
        "mean_relative_slope_difference": float(np.mean(rel_diff_all)),
        "max_relative_slope_difference": float(np.max(rel_diff_all)),
        "mean_cross_R2_N8_to_N10": float(np.mean(cross_r2_8to10)),
        "min_cross_R2_N8_to_N10": float(np.min(cross_r2_8to10)),
        "mean_cross_R2_N10_to_N8": float(np.mean(cross_r2_10to8)),
        "min_cross_R2_N10_to_N8": float(np.min(cross_r2_10to8)),
    }

    payload = {
        "analysis_status": "POST_HOC_EXPLORATORY",
        "question": (
            "Do the nonparametric residuals obey a stable compensation law "
            "epsilon_mod = a + b epsilon_W, and does that law transfer between N8 and N10?"
        ),
        "source": (
            "Frozen nonparametric residual CSV files from "
            "posthoc_analysis/nonparametric_coupling_output."
        ),
        "bandwidth_fractions_Q_range": BANDWIDTHS,
        "bandwidth_results": bandwidth_results,
        "robustness_summary": robustness,
        "interpretation_boundary": (
            "A stable negative slope and cross-size transfer would support an exploratory "
            "quantitative compensation signature, not establish conservation of information "
            "or a fundamental transfer law."
        ),
    }

    (OUT / "compensation_law_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("PAPER29_POSTHOC_COMPENSATION_LAW_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

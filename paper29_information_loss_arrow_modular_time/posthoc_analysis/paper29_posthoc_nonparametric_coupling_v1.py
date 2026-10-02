#!/usr/bin/env python3
"""
Paper 29 — Post-hoc nonparametric residual coupling analysis.

Question:
Does the residual association between d_W and d_mod survive after removing
a smooth, nonparametric dependence on Q(s)?

Method:
- Frozen Paper 29 result files are read only.
- Separate leave-one-out local-linear Gaussian-kernel regressions are fitted:
      d_W   ~ f(Q)
      d_mod ~ g(Q)
- Bandwidths are fixed as fractions of the observed Q range:
      0.05, 0.10, 0.20, 0.30, 0.40
- Residual Pearson and Spearman correlations are reported for each bandwidth.
- The same analysis is run independently for N8/A4 and N10/A4.

Status:
POST_HOC_EXPLORATORY.
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "paper29_arrow_of_modular_time_v1"
OUT = ROOT / "posthoc_analysis" / "nonparametric_coupling_output"

BANDWIDTH_FRACTIONS = [0.05, 0.10, 0.20, 0.30, 0.40]
WEIGHT_FLOOR = 1e-14


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fcol(rows, key):
    return np.asarray([float(r[key]) for r in rows], dtype=float)


def save_json(path: Path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def rankdata(a):
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    i = 0
    while i < len(a):
        j = i + 1
        while j < len(a) and a[order[j]] == a[order[i]]:
            j += 1
        rank = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = rank
        i = j
    return ranks


def pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def local_linear_loo_predict(x, y, bandwidth):
    """
    Leave-one-out local-linear Gaussian-kernel regression.

    At x_i, point i is excluded from the fit. A weighted local line in
    dx = x_j - x_i is fitted, and the intercept is the prediction at x_i.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if bandwidth <= 0.0:
        raise ValueError("bandwidth must be positive")

    n = len(x)
    pred = np.empty(n, dtype=float)
    effective_weight_sums = np.empty(n, dtype=float)

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False

        dx = x[mask] - x[i]
        yy = y[mask]
        w = np.exp(-0.5 * (dx / bandwidth) ** 2)

        effective_weight_sums[i] = float(np.sum(w))

        if np.sum(w) <= WEIGHT_FLOOR:
            # Deterministic fallback: nearest neighbor excluding self.
            j_local = int(np.argmin(np.abs(dx)))
            pred[i] = float(yy[j_local])
            continue

        X = np.column_stack([np.ones_like(dx), dx])
        sqrt_w = np.sqrt(w)
        Xw = X * sqrt_w[:, None]
        yw = yy * sqrt_w

        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        pred[i] = float(beta[0])

    return pred, effective_weight_sums


def write_rows(path: Path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def analyze(label, states_name, trans_name):
    states = read_csv(RESULTS / states_name)
    trans = read_csv(RESULTS / trans_name)

    q = fcol(states, "Q")[:-1]
    s = fcol(trans, "s_from")
    d_w = fcol(trans, "d_W")
    d_mod = fcol(trans, "d_mod")

    q_range = float(np.max(q) - np.min(q))
    if q_range <= 0.0:
        raise RuntimeError(f"NONPOSITIVE_Q_RANGE_{label}")

    raw = {
        "pearson_dW_dmod": pearson(d_w, d_mod),
        "spearman_dW_dmod": spearman(d_w, d_mod),
        "Q_min": float(np.min(q)),
        "Q_max": float(np.max(q)),
        "Q_range": q_range,
    }

    bandwidth_results = []
    residual_csv_rows = []

    for frac in BANDWIDTH_FRACTIONS:
        h = frac * q_range

        fit_w, weight_w = local_linear_loo_predict(q, d_w, h)
        fit_m, weight_m = local_linear_loo_predict(q, d_mod, h)

        resid_w = d_w - fit_w
        resid_m = d_mod - fit_m

        r_p = pearson(resid_w, resid_m)
        r_s = spearman(resid_w, resid_m)

        bandwidth_results.append({
            "bandwidth_fraction_Q_range": float(frac),
            "bandwidth_absolute": float(h),
            "pearson_residuals": r_p,
            "spearman_residuals": r_s,
            "residual_std_dW": float(np.std(resid_w)),
            "residual_std_dmod": float(np.std(resid_m)),
            "min_effective_weight_sum_dW": float(np.min(weight_w)),
            "min_effective_weight_sum_dmod": float(np.min(weight_m)),
        })

        for i in range(len(q)):
            residual_csv_rows.append({
                "label": label,
                "bandwidth_fraction_Q_range": float(frac),
                "bandwidth_absolute": float(h),
                "n": int(i),
                "s": float(s[i]),
                "Q": float(q[i]),
                "d_W": float(d_w[i]),
                "d_mod": float(d_mod[i]),
                "fit_dW": float(fit_w[i]),
                "fit_dmod": float(fit_m[i]),
                "resid_dW": float(resid_w[i]),
                "resid_dmod": float(resid_m[i]),
            })

        plt.figure(figsize=(7, 5))
        plt.scatter(resid_w, resid_m)
        plt.axhline(0.0, linewidth=1)
        plt.axvline(0.0, linewidth=1)
        plt.xlabel("Residual d_W")
        plt.ylabel("Residual d_mod")
        plt.title(f"Paper 29 post-hoc — LOO local-linear residuals, {label}, h={frac:.2f} ΔQ")
        plt.tight_layout()
        plt.savefig(
            OUT / f"residual_scatter_{label}_bw_{frac:.2f}.png",
            dpi=180,
        )
        plt.close()

    write_rows(OUT / f"nonparametric_residuals_{label}.csv", residual_csv_rows)

    pearsons = np.asarray(
        [r["pearson_residuals"] for r in bandwidth_results],
        dtype=float,
    )
    spearmans = np.asarray(
        [r["spearman_residuals"] for r in bandwidth_results],
        dtype=float,
    )

    summary = {
        "label": label,
        "n_transitions": int(len(trans)),
        "raw": raw,
        "bandwidth_results": bandwidth_results,
        "sign_robustness": {
            "pearson_negative_count": int(np.sum(pearsons < 0.0)),
            "pearson_positive_count": int(np.sum(pearsons > 0.0)),
            "spearman_negative_count": int(np.sum(spearmans < 0.0)),
            "spearman_positive_count": int(np.sum(spearmans > 0.0)),
            "pearson_min": float(np.min(pearsons)),
            "pearson_max": float(np.max(pearsons)),
            "spearman_min": float(np.min(spearmans)),
            "spearman_max": float(np.max(spearmans)),
        },
    }

    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    n8 = analyze(
        "N8_A4",
        "primary_states_N8_A4.csv",
        "primary_transitions_N8_A4.csv",
    )
    n10 = analyze(
        "N10_A4",
        "replication_states_N10_A4.csv",
        "replication_transitions_N10_A4.csv",
    )

    payload = {
        "analysis_status": "POST_HOC_EXPLORATORY",
        "method": (
            "Leave-one-out local-linear Gaussian-kernel regression of d_W and d_mod "
            "against Q, evaluated at five fixed bandwidth fractions of the observed Q range."
        ),
        "bandwidth_fractions_Q_range": BANDWIDTH_FRACTIONS,
        "datasets": [n8, n10],
        "interpretation_boundary": (
            "This is an exploratory robustness analysis. It does not modify or replace "
            "the preregistered Paper 29 endpoints."
        ),
    }

    save_json(OUT / "nonparametric_coupling_summary.json", payload)

    print("PAPER29_POSTHOC_NONPARAMETRIC_COUPLING_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

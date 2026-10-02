#!/usr/bin/env python3
"""
Paper 29 post-hoc residual coupling analysis.

Purpose:
Test whether the coupling between d_W and d_mod survives after removing
their shared dependence on Q(s) and, as a control, on s itself.

This is POST_HOC_EXPLORATORY and reads only frozen Paper 29 outputs.
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "paper29_arrow_of_modular_time_v1"
OUT = ROOT / "posthoc_analysis" / "residual_coupling_output"


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


def linear_residuals(y, x):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    resid = y - fitted
    return resid, fitted, beta


def quadratic_residuals(y, x):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones_like(x), x, x * x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    resid = y - fitted
    return resid, fitted, beta


def rank_residuals(y, x):
    ry = rankdata(y)
    rx = rankdata(x)
    resid, fitted, beta = linear_residuals(ry, rx)
    return resid, fitted, beta


def analyze(label, states_name, trans_name):
    states = read_csv(RESULTS / states_name)
    trans = read_csv(RESULTS / trans_name)

    q = fcol(states, "Q")[:-1]
    s = fcol(trans, "s_from")
    d_w = fcol(trans, "d_W")
    d_mod = fcol(trans, "d_mod")

    # Raw association.
    raw = {
        "pearson_dW_dmod": pearson(d_w, d_mod),
        "spearman_dW_dmod": spearman(d_w, d_mod),
    }

    # Linear residualization on Q.
    rw_q_lin, fw_q_lin, bw_q_lin = linear_residuals(d_w, q)
    rm_q_lin, fm_q_lin, bm_q_lin = linear_residuals(d_mod, q)

    # Quadratic residualization on Q: exploratory robustness check.
    rw_q_quad, fw_q_quad, bw_q_quad = quadratic_residuals(d_w, q)
    rm_q_quad, fm_q_quad, bm_q_quad = quadratic_residuals(d_mod, q)

    # Rank-based residualization on Q: partial-rank style check.
    rw_q_rank, _, _ = rank_residuals(d_w, q)
    rm_q_rank, _, _ = rank_residuals(d_mod, q)

    # Linear residualization on s as an explicit common-parameter control.
    rw_s_lin, fw_s_lin, bw_s_lin = linear_residuals(d_w, s)
    rm_s_lin, fm_s_lin, bm_s_lin = linear_residuals(d_mod, s)

    summary = {
        "label": label,
        "n_transitions": int(len(trans)),
        "raw": raw,
        "control_Q_linear": {
            "pearson_residuals": pearson(rw_q_lin, rm_q_lin),
            "spearman_residuals": spearman(rw_q_lin, rm_q_lin),
            "beta_dW_on_Q": [float(x) for x in bw_q_lin],
            "beta_dmod_on_Q": [float(x) for x in bm_q_lin],
            "residual_std_dW": float(np.std(rw_q_lin)),
            "residual_std_dmod": float(np.std(rm_q_lin)),
        },
        "control_Q_quadratic": {
            "pearson_residuals": pearson(rw_q_quad, rm_q_quad),
            "spearman_residuals": spearman(rw_q_quad, rm_q_quad),
            "beta_dW_on_Q_Q2": [float(x) for x in bw_q_quad],
            "beta_dmod_on_Q_Q2": [float(x) for x in bm_q_quad],
            "residual_std_dW": float(np.std(rw_q_quad)),
            "residual_std_dmod": float(np.std(rm_q_quad)),
        },
        "control_Q_rank_linear": {
            "pearson_rank_residuals": pearson(rw_q_rank, rm_q_rank),
            "spearman_rank_residuals": spearman(rw_q_rank, rm_q_rank),
        },
        "control_s_linear": {
            "pearson_residuals": pearson(rw_s_lin, rm_s_lin),
            "spearman_residuals": spearman(rw_s_lin, rm_s_lin),
            "beta_dW_on_s": [float(x) for x in bw_s_lin],
            "beta_dmod_on_s": [float(x) for x in bm_s_lin],
            "residual_std_dW": float(np.std(rw_s_lin)),
            "residual_std_dmod": float(np.std(rm_s_lin)),
        },
    }

    # Residual-residual plots.
    plt.figure(figsize=(7, 5))
    plt.scatter(rw_q_lin, rm_q_lin)
    plt.axhline(0.0, linewidth=1)
    plt.axvline(0.0, linewidth=1)
    plt.xlabel("Residual d_W after linear Q control")
    plt.ylabel("Residual d_mod after linear Q control")
    plt.title(f"Paper 29 post-hoc — residual coupling vs Q, {label}")
    plt.tight_layout()
    plt.savefig(OUT / f"residual_coupling_Q_linear_{label}.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.scatter(rw_q_quad, rm_q_quad)
    plt.axhline(0.0, linewidth=1)
    plt.axvline(0.0, linewidth=1)
    plt.xlabel("Residual d_W after quadratic Q control")
    plt.ylabel("Residual d_mod after quadratic Q control")
    plt.title(f"Paper 29 post-hoc — residual coupling vs Q², {label}")
    plt.tight_layout()
    plt.savefig(OUT / f"residual_coupling_Q_quadratic_{label}.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.scatter(rw_s_lin, rm_s_lin)
    plt.axhline(0.0, linewidth=1)
    plt.axvline(0.0, linewidth=1)
    plt.xlabel("Residual d_W after linear s control")
    plt.ylabel("Residual d_mod after linear s control")
    plt.title(f"Paper 29 post-hoc — residual coupling vs s, {label}")
    plt.tight_layout()
    plt.savefig(OUT / f"residual_coupling_s_linear_{label}.png", dpi=180)
    plt.close()

    # Save row-level residuals for audit.
    rows = []
    for i in range(len(s)):
        rows.append({
            "n": i,
            "s": float(s[i]),
            "Q": float(q[i]),
            "d_W": float(d_w[i]),
            "d_mod": float(d_mod[i]),
            "resid_dW_Q_linear": float(rw_q_lin[i]),
            "resid_dmod_Q_linear": float(rm_q_lin[i]),
            "resid_dW_Q_quadratic": float(rw_q_quad[i]),
            "resid_dmod_Q_quadratic": float(rm_q_quad[i]),
            "resid_rank_dW_Q": float(rw_q_rank[i]),
            "resid_rank_dmod_Q": float(rm_q_rank[i]),
            "resid_dW_s_linear": float(rw_s_lin[i]),
            "resid_dmod_s_linear": float(rm_s_lin[i]),
        })

    csv_path = OUT / f"residuals_{label}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

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
        "question": (
            "Does d_W remain associated with d_mod after removing shared dependence on Q(s), "
            "with s used as an explicit common-parameter control?"
        ),
        "datasets": [n8, n10],
        "interpretation_boundary": (
            "Residual associations are exploratory and do not retroactively modify preregistered endpoints."
        ),
    }

    save_json(OUT / "residual_coupling_summary.json", payload)
    print("PAPER29_POSTHOC_RESIDUAL_COUPLING_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

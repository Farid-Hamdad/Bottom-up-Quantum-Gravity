#!/usr/bin/env python3
"""
Paper 29 post-hoc exploratory analysis.

IMPORTANT:
- This script reads only the frozen Paper 29 outputs.
- It does not modify the preregistered script, preregistration, or frozen results.
- Every output produced here is explicitly post-hoc / exploratory.
"""

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "paper29_arrow_of_modular_time_v1"
OUT = ROOT / "posthoc_analysis" / "output"

PAIRS = ["01", "02", "03", "12", "13", "23"]


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return rows


def fcol(rows, key):
    return np.asarray([float(r[key]) for r in rows], dtype=float)


def normalize01(x):
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if hi == lo:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def save_json(path: Path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def dataset(label, state_file, transition_file):
    states = read_csv(RESULTS / state_file)
    trans = read_csv(RESULTS / transition_file)

    s = fcol(states, "s")
    q = fcol(states, "Q")

    s_from = fcol(trans, "s_from")
    sigma = fcol(trans, "Sigma")
    d_w = fcol(trans, "d_W")
    d_mod = fcol(trans, "d_mod")
    chi = fcol(trans, "chi")

    # 1) chi(s)
    plt.figure(figsize=(8, 5))
    plt.plot(s_from, chi, marker="o")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("s")
    plt.ylabel("chi")
    plt.title(f"Paper 29 post-hoc — directional cosine, {label}")
    plt.tight_layout()
    plt.savefig(OUT / f"chi_vs_s_{label}.png", dpi=180)
    plt.close()

    # 2) normalized magnitudes of the three preregistered transition quantities
    plt.figure(figsize=(8, 5))
    plt.plot(s_from, normalize01(sigma), marker="o", label="Sigma")
    plt.plot(s_from, normalize01(d_w), marker="o", label="d_W")
    plt.plot(s_from, normalize01(d_mod), marker="o", label="d_mod")
    plt.xlabel("s")
    plt.ylabel("min-max normalized value")
    plt.title(f"Paper 29 post-hoc — transition amplitudes, {label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"normalized_transition_amplitudes_{label}.png", dpi=180)
    plt.close()

    # 3) pairwise contributions to the dot product defining chi
    contribution_rows = []
    contribution_matrix = np.zeros((len(trans), len(PAIRS)), dtype=float)

    for k, row in enumerate(trans):
        for p_i, pair in enumerate(PAIRS):
            dw = float(row[f"delta_W_hat_{pair}"])
            dv = float(row[f"delta_v_hat_{pair}"])
            contribution_matrix[k, p_i] = dw * dv

        total_dot = float(np.sum(contribution_matrix[k]))
        dwn = float(row["d_W"])
        dvn = float(row["d_mod"])
        denom = dwn * dvn
        reconstructed_chi = total_dot / denom if denom > 0 else float("nan")

        out_row = {
            "n": int(float(row["n"])),
            "s_from": float(row["s_from"]),
            "s_to": float(row["s_to"]),
            "chi_frozen": float(row["chi"]),
            "chi_reconstructed": reconstructed_chi,
            "dot_total": total_dot,
        }
        for p_i, pair in enumerate(PAIRS):
            out_row[f"dot_contribution_{pair}"] = float(contribution_matrix[k, p_i])
        contribution_rows.append(out_row)

    write_csv(OUT / f"pair_dot_contributions_{label}.csv", contribution_rows)

    # mean signed pair contributions across the irreversible trajectory
    mean_pair_contrib = np.mean(contribution_matrix, axis=0)
    plt.figure(figsize=(8, 5))
    plt.bar(PAIRS, mean_pair_contrib)
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("pair")
    plt.ylabel("mean delta_W_hat * delta_v_hat")
    plt.title(f"Paper 29 post-hoc — mean directional contribution by pair, {label}")
    plt.tight_layout()
    plt.savefig(OUT / f"mean_pair_directional_contributions_{label}.png", dpi=180)
    plt.close()

    # 4) common-scalar diagnostics against Q(s)
    # Transition quantities are associated with the start of each step.
    q_from = q[:-1]

    def pearson(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        return float(np.corrcoef(a, b)[0, 1])

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

    def spearman(a, b):
        return pearson(rankdata(a), rankdata(b))

    scalar_summary = {
        "label": label,
        "n_transitions": int(len(trans)),
        "chi_min": float(np.min(chi)),
        "chi_max": float(np.max(chi)),
        "chi_mean": float(np.mean(chi)),
        "fraction_chi_positive": float(np.mean(chi > 0.0)),
        "pearson_Qfrom_dW": pearson(q_from, d_w),
        "pearson_Qfrom_dmod": pearson(q_from, d_mod),
        "pearson_Qfrom_Sigma": pearson(q_from, sigma),
        "spearman_Qfrom_dW": spearman(q_from, d_w),
        "spearman_Qfrom_dmod": spearman(q_from, d_mod),
        "spearman_Qfrom_Sigma": spearman(q_from, sigma),
        "pearson_dW_dmod": pearson(d_w, d_mod),
        "spearman_dW_dmod": spearman(d_w, d_mod),
        "mean_pair_dot_contributions": {
            pair: float(mean_pair_contrib[i]) for i, pair in enumerate(PAIRS)
        },
        "sum_mean_pair_dot_contributions": float(np.sum(mean_pair_contrib)),
    }

    # Scatter views to test whether d_W and d_mod look like monotonic functions of Q.
    plt.figure(figsize=(8, 5))
    plt.plot(q_from, d_w, marker="o", label="d_W")
    plt.plot(q_from, d_mod, marker="o", label="d_mod")
    plt.xlabel("Q(s) at transition start")
    plt.ylabel("step magnitude")
    plt.title(f"Paper 29 post-hoc — step magnitude vs Q, {label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"step_magnitudes_vs_Q_{label}.png", dpi=180)
    plt.close()

    return scalar_summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    summaries = [
        dataset(
            "N8_A4",
            "primary_states_N8_A4.csv",
            "primary_transitions_N8_A4.csv",
        ),
        dataset(
            "N10_A4",
            "replication_states_N10_A4.csv",
            "replication_transitions_N10_A4.csv",
        ),
    ]

    combined = {
        "analysis_status": "POST_HOC_EXPLORATORY",
        "source_results_commit_note": "Reads frozen Paper 29 result files only.",
        "datasets": summaries,
        "interpretation_boundary": (
            "These diagnostics are exploratory and must not be reclassified as preregistered tests."
        ),
    }

    save_json(OUT / "posthoc_summary.json", combined)
    print("PAPER29_POSTHOC_ANALYSIS_COMPLETED")
    print(OUT)


if __name__ == "__main__":
    main()

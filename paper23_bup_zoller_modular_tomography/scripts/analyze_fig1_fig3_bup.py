#!/usr/bin/env python3
"""BuP-oriented diagnostics from Zoller/Joshi figure 1 and 3 analysed data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from analyze_zoller_bup import load_mat_v5


ROOT = Path("/workspace/.cache")
OUT = Path("/workspace/bup_zoller_tomography/results_fig1_fig3")


def fit_linear(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[mask], dtype=float)
    y = np.asarray(y[mask], dtype=float)
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {"slope": float(a), "intercept": float(b), "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")}


def fit_const(x, y):
    mask = np.isfinite(y)
    y = np.asarray(y[mask], dtype=float)
    pred = np.full_like(y, np.mean(y))
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {"mean": float(np.mean(y)), "r2": 0.0 if ss_tot > 0 else float("nan"), "rmse": float(np.sqrt(np.mean((y - pred) ** 2)))}


def fit_entropy_file(path: Path, label: str):
    d = load_mat_v5(path)
    rows = []
    summary = {}
    for state in ["GS", "ES"]:
        arr = d[f"{state}ExpAvg"]
        L = arr[:, 0].astype(float)
        S = arr[:, 1].astype(float)
        mask = np.isfinite(S) & (L >= 2)
        L, S = L[mask], S[mask]
        lin_L = fit_linear(L, S)
        lin_log = fit_linear(np.log(L), S)
        const = fit_const(L, S)
        summary[f"{label}_{state}"] = {
            "n": int(len(L)),
            "S_first": float(S[0]),
            "S_last": float(S[-1]),
            "delta_S": float(S[-1] - S[0]),
            "linear_L": lin_L,
            "linear_logL": lin_log,
            "const": const,
            "best_scaling": "volume_like_linear_L" if lin_L["r2"] > lin_log["r2"] and lin_L["slope"] > 0.15 else "area_or_log_like",
        }
        for l, s in zip(L, S):
            rows.append({"delta": label, "state": state, "L_A": float(l), "S_vN": float(s)})
    return rows, summary


def analyze_beta_vs_j():
    files = {
        "Delta1_ES": ROOT / "01-BetaVsJDelta1p0ES.mat",
        "Delta1_GS": ROOT / "02-BetaVsJDelta1p0GS.mat",
        "Delta1p7_ES": ROOT / "03-BetaVsJDelta1p7ES.mat",
        "Delta1p7_GS": ROOT / "04-BetaVsJDelta1p7GS.mat",
    }
    summary = {}
    rows = []
    for label, path in files.items():
        d = load_mat_v5(path)
        j = d["indexj"].ravel().astype(float)
        beta = d["meanall"].ravel().astype(float)
        mask = beta > 1e-9
        jj = j[mask]
        bb = beta[mask]
        span = jj.max() - jj.min()
        u = (jj - jj.min()) / span if span else jj
        cft = u * (1 - u)
        tri = np.minimum(u, 1 - u)
        flat_cv = float(np.std(bb) / np.mean(bb))
        summary[label] = {
            "n": int(len(bb)),
            "beta_max": float(np.max(bb)),
            "beta_cv": flat_cv,
            "pearson_cft_parabola": float(np.corrcoef(bb, cft)[0, 1]),
            "pearson_triangle": float(np.corrcoef(bb, tri)[0, 1]),
            "T_edge_to_center": float(((1 / bb[0] + 1 / bb[-1]) / 2) / (1 / bb[len(bb) // 2])),
        }
        for x, y in zip(j, beta):
            rows.append({"dataset": label, "j": float(x), "beta_mean": float(y)})
    return rows, summary


def analyze_mi_and_fidelity():
    dmi = load_mat_v5(ROOT / "10-MutualInformation.mat")
    avg = dmi["dataMIAvg"]
    theory = dmi["dataMITheory"]
    dist, mi = avg[:, 0].astype(float), avg[:, 1].astype(float)
    mask = (dist >= 1) & np.isfinite(mi) & (mi > 0)
    exp_decay = fit_linear(dist[mask], np.log(mi[mask]))
    power_decay = fit_linear(np.log(dist[mask]), np.log(mi[mask]))

    dfid = load_mat_v5(ROOT / "09-FidelitiesLinksWoutLinks.mat")
    with_links = dfid["dataWithLinksAvg"]
    without_links = dfid["dataWithoutLinksAvg"]
    common = min(len(with_links), len(without_links))
    gain = with_links[:common, 1] - without_links[:common, 1]

    return {
        "mutual_information_decay": {
            "MI_at_distance_0": float(mi[dist == 0][0]) if np.any(dist == 0) else None,
            "MI_at_max_distance": float(mi[-1]),
            "exp_decay_log_slope": exp_decay,
            "power_decay_log_slope": power_decay,
        },
        "fidelity_links": {
            "mean_with_links": float(np.mean(with_links[:, 1])),
            "mean_without_links": float(np.mean(without_links[:, 1])),
            "mean_gain": float(np.mean(gain)),
            "min_gain": float(np.min(gain)),
            "max_gain": float(np.max(gain)),
        },
    }


def analyze_beta_mats():
    d = load_mat_v5(ROOT / "07-betaMats.mat")
    summary = {}
    for name, mat in d.items():
        if not isinstance(mat, np.ndarray):
            continue
        # Treat the two 5-site intervals as blocks 0:5 and 5:10.
        intra = np.concatenate([mat[:5, :5].ravel(), mat[5:, 5:].ravel()])
        inter = mat[:5, 5:].ravel()
        intra_abs = np.abs(intra[np.abs(intra) > 1e-12])
        inter_abs = np.abs(inter[np.abs(inter) > 1e-12])
        summary[name] = {
            "mean_abs_intra": float(np.mean(intra_abs)),
            "mean_abs_inter": float(np.mean(inter_abs)),
            "inter_over_intra": float(np.mean(inter_abs) / np.mean(intra_abs)),
            "max_abs_inter": float(np.max(inter_abs)),
        }
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    entropy_rows = []
    entropy_summary = {}
    for label, path in [("Delta1", ROOT / "05-EntropyDelta1p0.mat"), ("Delta1p7", ROOT / "06-EntropyDelta1p7.mat")]:
        rows, summ = fit_entropy_file(path, label)
        entropy_rows.extend(rows)
        entropy_summary.update(summ)

    with (OUT / "entropy_scaling_data.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["delta", "state", "L_A", "S_vN"])
        writer.writeheader()
        writer.writerows(entropy_rows)

    beta_rows, beta_summary = analyze_beta_vs_j()
    with (OUT / "beta_vs_j_profiles.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["dataset", "j", "beta_mean"])
        writer.writeheader()
        writer.writerows(beta_rows)

    summary = {
        "entropy_scaling": entropy_summary,
        "beta_vs_j": beta_summary,
        **analyze_mi_and_fidelity(),
        "beta_mats": analyze_beta_mats(),
    }
    (OUT / "fig1_fig3_bup_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

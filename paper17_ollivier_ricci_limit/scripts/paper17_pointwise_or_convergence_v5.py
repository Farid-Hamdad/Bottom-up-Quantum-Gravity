#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 17 — Pointwise Ollivier-Ricci convergence test v5

This script tests the remaining pointwise-convergence question:

    y_ij = kappa_ij^OR / r_rms(i,j)^2
    y_ij ≈ alpha_N + beta_N K(x_ij)

on the compact-support variable-curvature conformal torus protocol that gave
strong high-resolution signals in v4b.

Recommended run:

cd ~/bottomup

python3 papers/paper17_ollivier_ricci_limit/scripts/paper17_pointwise_or_convergence_v5.py \
  --N-values 784 1024 1225 1600 \
  --amplitude 0.20 \
  --epsilon-fixed 0.003 \
  --measure-radius-fixed 0.08 \
  --edge-max-length 0.06 \
  --idleness 0.5 \
  --n-edges-per-stratum 250 \
  --graph-neighbor-radius 2 \
  --k-threshold-quantile 0.35 \
  --max-support 12 \
  --seeds 11 22 33 \
  --n-curvature-bins 8 \
  --output-dir papers/paper17_ollivier_ricci_limit/results/pointwise_or_convergence_v5
"""

from __future__ import annotations

import argparse
import json
import math
import heapq
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.stats import spearmanr, pearsonr


def make_periodic_grid(N: int):
    m = int(round(math.sqrt(N)))
    xs = np.arange(m) / m
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    return coords, m, m * m


def conformal_factor(coords: np.ndarray, amplitude: float):
    return amplitude * np.cos(2 * np.pi * coords[:, 0]) * np.cos(2 * np.pi * coords[:, 1])


def gaussian_curvature(coords: np.ndarray, amplitude: float):
    f = conformal_factor(coords, amplitude)
    base = np.cos(2 * np.pi * coords[:, 0]) * np.cos(2 * np.pi * coords[:, 1])
    return 8.0 * np.pi * np.pi * amplitude * base * np.exp(-2.0 * f)


def grid_neighbor_edges(coords: np.ndarray, m: int, amplitude: float, neighbor_radius: int):
    f = conformal_factor(coords, amplitude)
    idx = np.arange(m * m).reshape(m, m)
    adj = [[] for _ in range(m * m)]
    offsets = [(dy, dx) for dy in range(-neighbor_radius, neighbor_radius + 1)
               for dx in range(-neighbor_radius, neighbor_radius + 1)
               if not (dx == 0 and dy == 0)]
    for iy in range(m):
        for ix in range(m):
            i = idx[iy, ix]
            for dy, dx in offsets:
                jy = (iy + dy) % m
                jx = (ix + dx) % m
                j = idx[jy, jx]
                if j <= i:
                    continue
                ddx = min(abs(dx) / m, 1.0 - abs(dx) / m)
                ddy = min(abs(dy) / m, 1.0 - abs(dy) / m)
                d0 = math.sqrt(ddx * ddx + ddy * ddy)
                length = math.exp(0.5 * (f[i] + f[j])) * d0
                adj[i].append((j, length))
                adj[j].append((i, length))
    return adj


def all_pairs_dijkstra(adj):
    N = len(adj)
    dist = np.full((N, N), np.inf, dtype=float)
    for s in range(N):
        d = np.full(N, np.inf, dtype=float)
        d[s] = 0.0
        pq = [(0.0, s)]
        while pq:
            du, u = heapq.heappop(pq)
            if du != d[u]:
                continue
            for v, w in adj[u]:
                nd = du + w
                if nd < d[v]:
                    d[v] = nd
                    heapq.heappush(pq, (nd, v))
        dist[s] = d
    return dist


def midpoint_periodic(coords, i, j):
    a = coords[i].copy()
    b = coords[j].copy()
    delta = b - a
    delta = (delta + 0.5) % 1.0 - 0.5
    return (a + 0.5 * delta) % 1.0


def heat_kernel(dist, eps):
    K = np.exp(-(dist * dist) / (4.0 * eps))
    np.fill_diagonal(K, 0.0)
    return K


def local_measure_radius(i, K, dist, radius, idleness, max_support=None):
    row = K[i].copy()
    row[i] = 0.0
    neigh = np.where((dist[i] > 0) & (dist[i] <= radius) & np.isfinite(dist[i]))[0]
    if max_support is not None and len(neigh) > max_support:
        weights_tmp = row[neigh]
        idx = np.argpartition(weights_tmp, -max_support)[-max_support:]
        neigh = neigh[idx]
    weights = row[neigh]
    if len(neigh) == 0 or weights.sum() <= 0:
        return np.array([i], dtype=int), np.array([1.0], dtype=float), 0.0, 0.0
    p_neigh = (1.0 - idleness) * weights / weights.sum()
    support = np.concatenate([[i], neigh])
    prob = np.concatenate([[idleness], p_neigh])
    prob = prob / prob.sum()
    d = dist[i, neigh]
    mean_r = float(np.average(d, weights=weights))
    rms_r = float(np.sqrt(np.average(d * d, weights=weights)))
    return support.astype(int), prob.astype(float), mean_r, rms_r


def wasserstein_1(supp_a, prob_a, supp_b, prob_b, dist):
    na, nb = len(supp_a), len(supp_b)
    C = dist[np.ix_(supp_a, supp_b)]
    c = C.ravel()
    A_eq, b_eq = [], []
    for a in range(na):
        row = np.zeros(na * nb)
        row[a * nb:(a + 1) * nb] = 1.0
        A_eq.append(row)
        b_eq.append(prob_a[a])
    for b in range(nb):
        row = np.zeros(na * nb)
        row[b::nb] = 1.0
        A_eq.append(row)
        b_eq.append(prob_b[b])
    res = linprog(c, A_eq=np.asarray(A_eq), b_eq=np.asarray(b_eq), bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.fun)


def build_candidate_edges(coords, dist, amplitude, edge_max_length):
    N = dist.shape[0]
    iu = np.triu_indices(N, k=1)
    dvals = dist[iu]
    mask = np.isfinite(dvals) & (dvals > 0) & (dvals <= edge_max_length)
    rows = []
    for i, j, d in zip(iu[0][mask], iu[1][mask], dvals[mask]):
        mid = midpoint_periodic(coords, int(i), int(j))
        K_mid = float(gaussian_curvature(mid.reshape(1, 2), amplitude)[0])
        rows.append((int(i), int(j), float(d), K_mid))
    if not rows:
        raise RuntimeError("No candidate edges. Increase --edge-max-length or --graph-neighbor-radius.")
    return pd.DataFrame(rows, columns=["i", "j", "edge_length", "K_mid"])


def stratified_sample_edges(candidates, n_per_stratum, k_threshold_quantile, rng):
    absK = candidates["K_mid"].abs()
    threshold = float(absK.quantile(k_threshold_quantile))
    neg = candidates[candidates["K_mid"] < -threshold].copy()
    neu = candidates[candidates["K_mid"].abs() <= threshold].copy()
    pos = candidates[candidates["K_mid"] > threshold].copy()
    samples = []
    for label, sub in [("negative", neg), ("neutral", neu), ("positive", pos)]:
        if len(sub) == 0:
            continue
        n = min(n_per_stratum, len(sub))
        w = 1.0 / np.maximum(sub["edge_length"].values, 1e-12)
        w = w / w.sum()
        idx = rng.choice(sub.index.values, size=n, replace=False, p=w)
        s = sub.loc[idx].copy()
        s["stratum"] = label
        s["K_threshold"] = threshold
        samples.append(s)
    if not samples:
        raise RuntimeError("No stratified samples produced.")
    return pd.concat(samples, ignore_index=True)


def compute_case(N_input, seed, args):
    rng = np.random.default_rng(seed)
    coords, m, N = make_periodic_grid(N_input)
    adj = grid_neighbor_edges(coords, m, args.amplitude, args.graph_neighbor_radius)
    dist = all_pairs_dijkstra(adj)
    eps = float(args.epsilon_fixed)
    Kheat = heat_kernel(dist, eps)
    candidates = build_candidate_edges(coords, dist, args.amplitude, args.edge_max_length)
    sampled = stratified_sample_edges(candidates, args.n_edges_per_stratum, args.k_threshold_quantile, rng)
    cache, rows = {}, []
    for edge_idx, row in sampled.iterrows():
        i, j = int(row["i"]), int(row["j"])
        if i not in cache:
            cache[i] = local_measure_radius(i, Kheat, dist, args.measure_radius_fixed, args.idleness, args.max_support)
        if j not in cache:
            cache[j] = local_measure_radius(j, Kheat, dist, args.measure_radius_fixed, args.idleness, args.max_support)
        si, pi, mean_ri, rms_ri = cache[i]
        sj, pj, mean_rj, rms_rj = cache[j]
        dij = float(dist[i, j])
        if dij <= 0:
            continue
        W1 = wasserstein_1(si, pi, sj, pj, dist)
        kappa = 1.0 - W1 / dij
        mean_r = 0.5 * (mean_ri + mean_rj)
        rms_r = 0.5 * (rms_ri + rms_rj)
        rows.append({
            "N_input": int(N_input), "N_actual": int(N), "m_grid": int(m), "seed": int(seed),
            "amplitude": float(args.amplitude), "epsilon_fixed": eps,
            "measure_radius_fixed": float(args.measure_radius_fixed), "edge_max_length": float(args.edge_max_length),
            "edge_index": int(edge_idx), "stratum": str(row["stratum"]), "K_threshold": float(row["K_threshold"]),
            "i": i, "j": j, "edge_length": dij, "measure_mean_radius": mean_r, "measure_rms_radius": rms_r,
            "support_i": int(len(si)), "support_j": int(len(sj)), "K_mid": float(row["K_mid"]),
            "W1": float(W1), "kappa_or": float(kappa),
            "kappa_over_epsilon": float(kappa / eps),
            "kappa_over_l2": float(kappa / (dij * dij)),
            "kappa_over_measure_mean_r2": float(kappa / (mean_r * mean_r)) if mean_r > 0 else np.nan,
            "kappa_over_measure_rms_r2": float(kappa / (rms_r * rms_r)) if rms_r > 0 else np.nan,
        })
    return rows


def fit_affine(K_mid, metric_values):
    x = pd.to_numeric(pd.Series(K_mid), errors="coerce").replace([np.inf, -np.inf], np.nan).values
    y = pd.to_numeric(pd.Series(metric_values), errors="coerce").replace([np.inf, -np.inf], np.nan).values
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return None
    X = np.column_stack([np.ones(ok.sum()), x[ok]])
    beta, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    pred = X @ beta
    residual = y[ok] - pred
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y[ok] - np.mean(y[ok])) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    sp = spearmanr(x[ok], y[ok])
    pe = pearsonr(x[ok], y[ok])
    centered_y = y[ok] - beta[0]
    sign_mask = np.abs(x[ok]) > np.quantile(np.abs(x[ok]), 0.20)
    sign_accuracy = float((np.sign(centered_y[sign_mask]) == np.sign(x[ok][sign_mask])).mean()) if sign_mask.sum() else np.nan
    metric_std = float(np.std(y[ok], ddof=0))
    return {
        "n": int(ok.sum()), "intercept": float(beta[0]), "slope": float(beta[1]), "r2": r2,
        "spearman": float(sp.statistic), "spearman_p": float(sp.pvalue),
        "pearson": float(pe.statistic), "pearson_p": float(pe.pvalue),
        "residual_mean": float(np.mean(residual)), "residual_std": float(np.std(residual, ddof=0)),
        "residual_mae": float(np.mean(np.abs(residual))), "metric_std": metric_std,
        "normalized_residual_std": float(np.std(residual, ddof=0) / metric_std) if metric_std > 0 else np.nan,
        "normalized_residual_mae": float(np.mean(np.abs(residual)) / metric_std) if metric_std > 0 else np.nan,
        "sign_accuracy": sign_accuracy,
    }


def summarize_seed(df, metric):
    fit = fit_affine(df["K_mid"], df[metric])
    if fit is None:
        return None
    neg = df[df["stratum"] == "negative"]
    pos = df[df["stratum"] == "positive"]
    def mean_metric(s):
        return float(pd.to_numeric(s[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).mean()) if len(s) else np.nan
    fit.update({
        "metric": metric,
        "negative_n": int(len(neg)), "positive_n": int(len(pos)),
        "negative_mean": mean_metric(neg), "positive_mean": mean_metric(pos),
        "delta_positive_minus_negative": mean_metric(pos) - mean_metric(neg),
        "edge_length_mean": float(df["edge_length"].mean()),
        "measure_rms_radius_mean": float(df["measure_rms_radius"].mean()),
        "support_mean": float(0.5 * (df["support_i"].mean() + df["support_j"].mean())),
    })
    return fit


def curvature_bin_summary(df, metric, n_bins):
    d = df.copy()
    d[metric] = pd.to_numeric(d[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)
    d = d[np.isfinite(d[metric]) & np.isfinite(d["K_mid"])].copy()
    if len(d) < n_bins:
        return pd.DataFrame()
    d["K_bin"] = pd.qcut(d["K_mid"], q=n_bins, duplicates="drop")
    rows = []
    for b, sub in d.groupby("K_bin", observed=False):
        if len(sub) == 0:
            continue
        rows.append({
            "metric": metric, "K_bin": str(b), "K_bin_mean": float(sub["K_mid"].mean()),
            "metric_mean": float(sub[metric].mean()), "metric_std": float(sub[metric].std(ddof=0)),
            "n_edges": int(len(sub)),
        })
    return pd.DataFrame(rows)


def trend_vs_logN(N_values, y_values):
    xN = np.asarray(N_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    ok = np.isfinite(xN) & np.isfinite(y) & (xN > 0)
    if ok.sum() < 3:
        return np.nan, np.nan
    X = np.column_stack([np.ones(ok.sum()), np.log(xN[ok])])
    beta, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((y[ok] - pred) ** 2))
    ss_tot = float(np.sum((y[ok] - np.mean(y[ok])) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    return float(beta[1]), r2


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-values", nargs="+", type=int, default=[784, 1024, 1225, 1600])
    p.add_argument("--amplitude", type=float, default=0.20)
    p.add_argument("--epsilon-fixed", type=float, default=0.003)
    p.add_argument("--measure-radius-fixed", type=float, default=0.08)
    p.add_argument("--edge-max-length", type=float, default=0.06)
    p.add_argument("--idleness", type=float, default=0.5)
    p.add_argument("--n-edges-per-stratum", type=int, default=250)
    p.add_argument("--graph-neighbor-radius", type=int, default=2)
    p.add_argument("--k-threshold-quantile", type=float, default=0.35)
    p.add_argument("--max-support", type=int, default=12)
    p.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    p.add_argument("--n-curvature-bins", type=int, default=8)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    if args.max_support == 0:
        args.max_support = None
    outdir = Path(args.output_dir)
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for N in args.N_values:
        for seed in args.seeds:
            print(f"[run] N={N} seed={seed}")
            rows.extend(compute_case(N, seed, args))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "paper17_pointwise_or_v5_rows.csv", index=False)
    metrics = ["kappa_over_epsilon", "kappa_over_l2", "kappa_over_measure_mean_r2", "kappa_over_measure_rms_r2"]
    summary_rows, bin_rows = [], []
    for (N, seed), sub in df.groupby(["N_actual", "seed"]):
        for metric in metrics:
            s = summarize_seed(sub, metric)
            if s is None:
                continue
            s.update({
                "N_actual": int(N), "seed": int(seed), "amplitude": float(args.amplitude),
                "epsilon_fixed": float(args.epsilon_fixed), "measure_radius_fixed": float(args.measure_radius_fixed),
                "edge_max_length": float(args.edge_max_length), "n_edges": int(len(sub)),
            })
            summary_rows.append(s)
            bdf = curvature_bin_summary(sub, metric, args.n_curvature_bins)
            if len(bdf):
                bdf["N_actual"] = int(N)
                bdf["seed"] = int(seed)
                bin_rows.extend(bdf.to_dict(orient="records"))
    summary = pd.DataFrame(summary_rows).sort_values(["metric", "N_actual", "seed"])
    summary.to_csv(outdir / "paper17_pointwise_or_v5_summary_by_seed.csv", index=False)
    bins = pd.DataFrame(bin_rows)
    bins.to_csv(outdir / "paper17_pointwise_or_v5_bins.csv", index=False)
    agg_rows = []
    for (N, metric), sub in summary.groupby(["N_actual", "metric"]):
        agg_rows.append({
            "N_actual": int(N), "metric": metric, "n_seeds": int(len(sub)),
            "slope_mean": float(sub["slope"].mean()), "slope_std": float(sub["slope"].std(ddof=0)),
            "r2_mean": float(sub["r2"].mean()), "r2_std": float(sub["r2"].std(ddof=0)),
            "spearman_mean": float(sub["spearman"].mean()), "spearman_std": float(sub["spearman"].std(ddof=0)),
            "pearson_mean": float(sub["pearson"].mean()), "pearson_std": float(sub["pearson"].std(ddof=0)),
            "residual_std_mean": float(sub["residual_std"].mean()),
            "residual_mae_mean": float(sub["residual_mae"].mean()),
            "normalized_residual_std_mean": float(sub["normalized_residual_std"].mean()),
            "normalized_residual_mae_mean": float(sub["normalized_residual_mae"].mean()),
            "sign_accuracy_mean": float(sub["sign_accuracy"].mean()), "sign_accuracy_std": float(sub["sign_accuracy"].std(ddof=0)),
            "delta_pos_minus_neg_mean": float(sub["delta_positive_minus_negative"].mean()),
            "delta_pos_minus_neg_std": float(sub["delta_positive_minus_negative"].std(ddof=0)),
            "positive_delta_fraction": float((sub["delta_positive_minus_negative"] > 0).mean()),
            "support_mean": float(sub["support_mean"].mean()),
            "edge_length_mean": float(sub["edge_length_mean"].mean()),
            "measure_rms_radius_mean": float(sub["measure_rms_radius_mean"].mean()),
        })
    agg = pd.DataFrame(agg_rows).sort_values(["metric", "N_actual"])
    agg.to_csv(outdir / "paper17_pointwise_or_v5_aggregate.csv", index=False)
    scaling_rows = []
    for metric, sub in agg.groupby("metric"):
        for observable in ["r2_mean", "spearman_mean", "normalized_residual_std_mean", "sign_accuracy_mean", "slope_mean"]:
            trend, trend_r2 = trend_vs_logN(sub["N_actual"].values, sub[observable].values)
            scaling_rows.append({
                "metric": metric, "observable": observable, "logN_trend": trend, "trend_r2": trend_r2,
                "first_value": float(sub[observable].iloc[0]), "last_value": float(sub[observable].iloc[-1]),
                "delta_last_minus_first": float(sub[observable].iloc[-1] - sub[observable].iloc[0]),
            })
    scaling = pd.DataFrame(scaling_rows)
    scaling.to_csv(outdir / "paper17_pointwise_or_v5_scaling.csv", index=False)
    primary = "kappa_over_measure_rms_r2"
    primary_agg = agg[agg["metric"] == primary]
    for obs, ylabel in [("r2_mean", "OLS R²"), ("spearman_mean", "Spearman"), ("normalized_residual_std_mean", "normalized residual std"), ("sign_accuracy_mean", "sign accuracy"), ("delta_pos_minus_neg_mean", "delta positive-negative"), ("slope_mean", "OLS slope")]:
        plt.figure()
        yerr_col = obs.replace("_mean", "_std")
        yerr = primary_agg[yerr_col] if yerr_col in primary_agg.columns else None
        plt.errorbar(primary_agg["N_actual"], primary_agg[obs], yerr=yerr, marker="o")
        if obs in ["spearman_mean", "delta_pos_minus_neg_mean", "slope_mean"]:
            plt.axhline(0, linestyle="--")
        plt.xlabel("N")
        plt.ylabel(ylabel)
        plt.title(f"Pointwise convergence diagnostic: {primary}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_primary_{obs}.png", dpi=160)
        plt.close()
    if len(bins):
        for N, subN in bins[bins["metric"] == primary].groupby("N_actual"):
            plt.figure()
            for seed, sub in subN.groupby("seed"):
                s = sub.sort_values("K_bin_mean")
                plt.plot(s["K_bin_mean"], s["metric_mean"], marker="o", alpha=0.5, label=f"seed={seed}")
            plt.xlabel("mean analytic curvature K in bin")
            plt.ylabel(f"mean {primary}")
            plt.title(f"Binned pointwise relation, N={N}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(figdir / f"fig_bins_{primary}_N{N}.png", dpi=160)
            plt.close()
    payload = {
        "experiment": "Paper 17 pointwise OR convergence v5",
        "args": vars(args),
        "aggregate": agg.to_dict(orient="records"),
        "scaling": scaling.to_dict(orient="records"),
        "primary_metric": primary,
    }
    with open(outdir / "paper17_pointwise_or_v5_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    lines = []
    lines.append("# Paper 17 — Pointwise OR convergence v5")
    lines.append("")
    lines.append("The pointwise test fits:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"y_{ij}=\alpha_N+\beta_N K(x_{ij})+\eta_{ij},\qquad y_{ij}=\kappa_{ij}^{OR}/r_{\rm rms}^2.")
    lines.append(r"\]")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append("| N | metric | slope | R2 | Spearman | norm resid std | sign acc | delta pos-neg | positive delta frac |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in agg.iterrows():
        lines.append(f"| {int(r['N_actual'])} | {r['metric']} | {r['slope_mean']:.6g} | {r['r2_mean']:.6g} | {r['spearman_mean']:.6g} | {r['normalized_residual_std_mean']:.6g} | {r['sign_accuracy_mean']:.6g} | {r['delta_pos_minus_neg_mean']:.6g} | {r['positive_delta_fraction']:.3g} |")
    lines.append("")
    lines.append("## Scaling")
    lines.append("")
    lines.append("| metric | observable | first | last | delta | logN trend | trend R2 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, r in scaling.iterrows():
        lines.append(f"| {r['metric']} | {r['observable']} | {r['first_value']:.6g} | {r['last_value']:.6g} | {r['delta_last_minus_first']:.6g} | {r['logN_trend']:.6g} | {r['trend_r2']:.6g} |")
    (outdir / "paper17_pointwise_or_v5_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("=" * 110)
    print("Paper 17 — Pointwise OR convergence v5")
    print("=" * 110)
    print(agg.to_string(index=False))
    print("\n" + "=" * 110)
    print("SCALING DIAGNOSTICS")
    print("=" * 110)
    print(scaling.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "paper17_pointwise_or_v5_rows.csv")
    print(outdir / "paper17_pointwise_or_v5_summary_by_seed.csv")
    print(outdir / "paper17_pointwise_or_v5_aggregate.csv")
    print(outdir / "paper17_pointwise_or_v5_bins.csv")
    print(outdir / "paper17_pointwise_or_v5_scaling.csv")
    print(outdir / "paper17_pointwise_or_v5_summary.json")
    print(outdir / "paper17_pointwise_or_v5_summary.md")
    print(figdir)


if __name__ == "__main__":
    main()

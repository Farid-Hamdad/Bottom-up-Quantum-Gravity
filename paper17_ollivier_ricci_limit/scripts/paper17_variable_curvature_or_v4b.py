#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 17 — Variable-curvature local Ollivier-Ricci test v4b

Purpose
-------
v4 showed a stable positive local Ricci signal for N=400,576,784 but a sign
reversal at N=1024. The suspected issue is that the numerical protocol used
kNN-dependent scales:

    epsilon_N ~ d_kNN^2,
    measure support from fixed k nearest neighbors,

so the physical scale of the kernel/measure changes with N.

v4b fixes the physical scales directly:

    --epsilon-fixed
    --measure-radius-fixed
    --edge-max-length

This tests local OR curvature at a fixed physical resolution while increasing N.

Geometry
--------
Conformal torus:

    g = exp(2 f) (dx^2 + dy^2),
    f(x,y) = A cos(2pi x) cos(2pi y).

In 2D:

    Ric(u,u) = K(x),

with

    K = - exp(-2f) Delta_flat f
      = 8 pi^2 A cos(2pi x) cos(2pi y) exp(-2f).

Main test
---------
For the metric

    kappa_over_measure_rms_r2,

we want, across N and seeds:

    Spearman(metric, K_mid) > 0
    OLS slope > 0
    delta_positive_minus_negative > 0
    positive_delta_fraction = 1

Recommended run
---------------
cd ~/bottomup

python3 papers/paper17_ollivier_ricci_limit/scripts/paper17_variable_curvature_or_v4b.py \
  --N-values 400 576 784 1024 1600 \
  --amplitude 0.10 \
  --epsilon-fixed 0.006 \
  --measure-radius-fixed 0.14 \
  --edge-max-length 0.10 \
  --idleness 0.5 \
  --n-edges-per-stratum 300 \
  --graph-neighbor-radius 3 \
  --k-threshold-quantile 0.35 \
  --seeds 11 22 33 \
  --output-dir papers/paper17_ollivier_ricci_limit/results/variable_curvature_or_v4b_fixed_epsilon

Notes
-----
- N=1600 may be slower because all-pairs Dijkstra is used.
- If N=1600 is too slow, first run N=400 576 784 1024.
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


# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------

def make_periodic_grid(N: int):
    m = int(round(math.sqrt(N)))
    xs = np.arange(m) / m
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    return coords, m, m*m


def conformal_factor(coords: np.ndarray, amplitude: float):
    return amplitude * np.cos(2*np.pi*coords[:, 0]) * np.cos(2*np.pi*coords[:, 1])


def gaussian_curvature(coords: np.ndarray, amplitude: float):
    f = conformal_factor(coords, amplitude)
    base = np.cos(2*np.pi*coords[:, 0]) * np.cos(2*np.pi*coords[:, 1])
    return 8.0*np.pi*np.pi*amplitude*base*np.exp(-2.0*f)


def grid_neighbor_edges(coords: np.ndarray, m: int, amplitude: float, neighbor_radius: int):
    f = conformal_factor(coords, amplitude)
    idx = np.arange(m*m).reshape(m, m)
    adj = [[] for _ in range(m*m)]

    offsets = []
    for dy in range(-neighbor_radius, neighbor_radius+1):
        for dx in range(-neighbor_radius, neighbor_radius+1):
            if dx == 0 and dy == 0:
                continue
            offsets.append((dy, dx))

    for iy in range(m):
        for ix in range(m):
            i = idx[iy, ix]
            for dy, dx in offsets:
                jy = (iy + dy) % m
                jx = (ix + dx) % m
                j = idx[jy, jx]
                if j <= i:
                    continue

                ddx = min(abs(dx)/m, 1.0 - abs(dx)/m)
                ddy = min(abs(dy)/m, 1.0 - abs(dy)/m)
                d0 = math.sqrt(ddx*ddx + ddy*ddy)
                fmid = 0.5*(f[i] + f[j])
                length = math.exp(fmid) * d0

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
    return (a + 0.5*delta) % 1.0


# -----------------------------------------------------------------------------
# OR machinery
# -----------------------------------------------------------------------------

def heat_kernel(dist, eps):
    K = np.exp(-(dist*dist)/(4.0*eps))
    np.fill_diagonal(K, 0.0)
    return K


def local_measure_radius(i, K, dist, radius, idleness, max_support=None):
    row = K[i].copy()
    row[i] = 0.0

    neigh = np.where((dist[i] > 0) & (dist[i] <= radius) & np.isfinite(dist[i]))[0]

    if max_support is not None and len(neigh) > max_support:
        # Keep strongest kernel weights inside radius.
        weights_tmp = row[neigh]
        idx = np.argpartition(weights_tmp, -max_support)[-max_support:]
        neigh = neigh[idx]

    weights = row[neigh]
    if len(neigh) == 0 or weights.sum() <= 0:
        return np.array([i], dtype=int), np.array([1.0], dtype=float), 0.0, 0.0

    p_neigh = (1.0-idleness)*weights/weights.sum()
    support = np.concatenate([[i], neigh])
    prob = np.concatenate([[idleness], p_neigh])
    prob = prob/prob.sum()

    d = dist[i, neigh]
    mean_r = float(np.average(d, weights=weights))
    rms_r = float(np.sqrt(np.average(d*d, weights=weights)))
    return support.astype(int), prob.astype(float), mean_r, rms_r


def wasserstein_1(supp_a, prob_a, supp_b, prob_b, dist):
    na, nb = len(supp_a), len(supp_b)
    C = dist[np.ix_(supp_a, supp_b)]
    c = C.ravel()

    A_eq = []
    b_eq = []

    for a in range(na):
        row = np.zeros(na*nb)
        row[a*nb:(a+1)*nb] = 1.0
        A_eq.append(row)
        b_eq.append(prob_a[a])

    for b in range(nb):
        row = np.zeros(na*nb)
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
        raise RuntimeError(
            f"No candidate edges with edge_max_length={edge_max_length}. "
            "Increase --edge-max-length or --graph-neighbor-radius."
        )

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
    sampled = stratified_sample_edges(
        candidates,
        args.n_edges_per_stratum,
        args.k_threshold_quantile,
        rng,
    )

    cache = {}
    rows = []

    for edge_idx, row in sampled.iterrows():
        i = int(row["i"])
        j = int(row["j"])

        if i not in cache:
            cache[i] = local_measure_radius(
                i, Kheat, dist,
                radius=args.measure_radius_fixed,
                idleness=args.idleness,
                max_support=args.max_support,
            )
        if j not in cache:
            cache[j] = local_measure_radius(
                j, Kheat, dist,
                radius=args.measure_radius_fixed,
                idleness=args.idleness,
                max_support=args.max_support,
            )

        si, pi, mean_ri, rms_ri = cache[i]
        sj, pj, mean_rj, rms_rj = cache[j]

        dij = float(dist[i, j])
        if dij <= 0:
            continue

        W1 = wasserstein_1(si, pi, sj, pj, dist)
        kappa = 1.0 - W1/dij

        mean_r = 0.5*(mean_ri + mean_rj)
        rms_r = 0.5*(rms_ri + rms_rj)

        rows.append({
            "N_input": int(N_input),
            "N_actual": int(N),
            "m_grid": int(m),
            "seed": int(seed),
            "amplitude": float(args.amplitude),
            "epsilon_fixed": eps,
            "measure_radius_fixed": float(args.measure_radius_fixed),
            "edge_max_length": float(args.edge_max_length),
            "graph_neighbor_radius": int(args.graph_neighbor_radius),
            "idleness": float(args.idleness),
            "edge_index": int(edge_idx),
            "stratum": str(row["stratum"]),
            "K_threshold": float(row["K_threshold"]),
            "i": i,
            "j": j,
            "edge_length": dij,
            "measure_mean_radius": mean_r,
            "measure_rms_radius": rms_r,
            "support_i": int(len(si)),
            "support_j": int(len(sj)),
            "K_mid": float(row["K_mid"]),
            "W1": float(W1),
            "kappa_or": float(kappa),
            "kappa_over_epsilon": float(kappa/eps),
            "kappa_over_l2": float(kappa/(dij*dij)),
            "kappa_over_measure_mean_r2": float(kappa/(mean_r*mean_r)) if mean_r > 0 else np.nan,
            "kappa_over_measure_rms_r2": float(kappa/(rms_r*rms_r)) if rms_r > 0 else np.nan,
        })

    return rows


def corr_summary(sub, metric):
    x = pd.to_numeric(sub[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)
    y = pd.to_numeric(sub["K_mid"], errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)

    if ok.sum() < 5:
        return {
            "spearman": np.nan, "spearman_p": np.nan,
            "pearson": np.nan, "pearson_p": np.nan,
            "ols_intercept": np.nan, "ols_slope": np.nan, "ols_r2": np.nan,
        }

    sp = spearmanr(x[ok], y[ok])
    pe = pearsonr(x[ok], y[ok])

    X = np.column_stack([np.ones(ok.sum()), y[ok].values])
    beta, *_ = np.linalg.lstsq(X, x[ok].values, rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((x[ok].values - pred)**2))
    ss_tot = float(np.sum((x[ok].values - np.mean(x[ok].values))**2))
    r2 = float(1.0 - ss_res/ss_tot) if ss_tot > 0 else np.nan

    return {
        "spearman": float(sp.statistic),
        "spearman_p": float(sp.pvalue),
        "pearson": float(pe.statistic),
        "pearson_p": float(pe.pvalue),
        "ols_intercept": float(beta[0]),
        "ols_slope": float(beta[1]),
        "ols_r2": r2,
    }


def summarize_case(df, metric):
    c = corr_summary(df, metric)

    neg = df[df["stratum"] == "negative"]
    neu = df[df["stratum"] == "neutral"]
    pos = df[df["stratum"] == "positive"]

    def mean_metric(s):
        return float(pd.to_numeric(s[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).mean()) if len(s) else np.nan

    return {
        "metric": metric,
        **c,
        "negative_n": int(len(neg)),
        "neutral_n": int(len(neu)),
        "positive_n": int(len(pos)),
        "negative_mean": mean_metric(neg),
        "neutral_mean": mean_metric(neu),
        "positive_mean": mean_metric(pos),
        "delta_positive_minus_negative": mean_metric(pos)-mean_metric(neg),
        "K_mid_mean": float(df["K_mid"].mean()),
        "K_mid_std": float(df["K_mid"].std(ddof=0)),
        "edge_length_mean": float(df["edge_length"].mean()),
        "edge_length_std": float(df["edge_length"].std(ddof=0)),
        "measure_rms_radius_mean": float(df["measure_rms_radius"].mean()),
        "support_mean": float(0.5*(df["support_i"].mean()+df["support_j"].mean())),
    }


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-values", nargs="+", type=int, default=[400, 576, 784, 1024, 1600])
    p.add_argument("--amplitude", type=float, default=0.10)
    p.add_argument("--epsilon-fixed", type=float, default=0.006)
    p.add_argument("--measure-radius-fixed", type=float, default=0.14)
    p.add_argument("--edge-max-length", type=float, default=0.10)
    p.add_argument("--idleness", type=float, default=0.5)
    p.add_argument("--n-edges-per-stratum", type=int, default=300)
    p.add_argument("--graph-neighbor-radius", type=int, default=3)
    p.add_argument("--k-threshold-quantile", type=float, default=0.35)
    p.add_argument("--max-support", type=int, default=80,
                   help="Cap support size for LP speed; use 0 for no cap.")
    p.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
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
    df.to_csv(outdir / "paper17_variable_curvature_or_v4b_rows.csv", index=False)

    metrics = [
        "kappa_over_epsilon",
        "kappa_over_l2",
        "kappa_over_measure_mean_r2",
        "kappa_over_measure_rms_r2",
    ]

    summary_rows = []
    for (N, seed), sub in df.groupby(["N_actual", "seed"]):
        for metric in metrics:
            s = summarize_case(sub, metric)
            s.update({
                "N_actual": int(N),
                "seed": int(seed),
                "amplitude": float(args.amplitude),
                "epsilon_fixed": float(args.epsilon_fixed),
                "measure_radius_fixed": float(args.measure_radius_fixed),
                "edge_max_length": float(args.edge_max_length),
                "n_edges": int(len(sub)),
            })
            summary_rows.append(s)

    summary = pd.DataFrame(summary_rows).sort_values(["metric", "N_actual", "seed"])
    summary.to_csv(outdir / "paper17_variable_curvature_or_v4b_summary_by_seed.csv", index=False)

    aggregate_rows = []
    for (N, metric), sub in summary.groupby(["N_actual", "metric"]):
        aggregate_rows.append({
            "N_actual": int(N),
            "metric": metric,
            "n_seeds": int(len(sub)),
            "spearman_mean": float(sub["spearman"].mean()),
            "spearman_std": float(sub["spearman"].std(ddof=0)),
            "pearson_mean": float(sub["pearson"].mean()),
            "pearson_std": float(sub["pearson"].std(ddof=0)),
            "ols_slope_mean": float(sub["ols_slope"].mean()),
            "ols_slope_std": float(sub["ols_slope"].std(ddof=0)),
            "ols_r2_mean": float(sub["ols_r2"].mean()),
            "ols_r2_std": float(sub["ols_r2"].std(ddof=0)),
            "delta_pos_minus_neg_mean": float(sub["delta_positive_minus_negative"].mean()),
            "delta_pos_minus_neg_std": float(sub["delta_positive_minus_negative"].std(ddof=0)),
            "positive_delta_fraction": float((sub["delta_positive_minus_negative"] > 0).mean()),
            "edge_length_mean": float(sub["edge_length_mean"].mean()),
            "measure_rms_radius_mean": float(sub["measure_rms_radius_mean"].mean()),
            "support_mean": float(sub["support_mean"].mean()),
        })

    agg = pd.DataFrame(aggregate_rows).sort_values(["metric", "N_actual"])
    agg.to_csv(outdir / "paper17_variable_curvature_or_v4b_aggregate.csv", index=False)

    # Figures
    for metric in ["kappa_over_measure_rms_r2", "kappa_over_l2", "kappa_over_epsilon"]:
        sub = agg[agg["metric"] == metric]

        plt.figure()
        plt.errorbar(sub["N_actual"], sub["spearman_mean"], yerr=sub["spearman_std"], marker="o")
        plt.axhline(0, linestyle="--")
        plt.xlabel("N")
        plt.ylabel("Spearman mean ± std")
        plt.title(f"v4b fixed-scale Spearman: {metric}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_spearman_{metric}.png", dpi=160)
        plt.close()

        plt.figure()
        plt.errorbar(sub["N_actual"], sub["delta_pos_minus_neg_mean"], yerr=sub["delta_pos_minus_neg_std"], marker="o")
        plt.axhline(0, linestyle="--")
        plt.xlabel("N")
        plt.ylabel("delta positive-negative mean ± std")
        plt.title(f"v4b fixed-scale sign-bin delta: {metric}")
        plt.tight_layout()
        plt.savefig(figdir / f"fig_delta_{metric}.png", dpi=160)
        plt.close()

    # Example scatter for first N and seed
    first_N = sorted(df["N_actual"].unique())[0]
    first_seed = sorted(df["seed"].unique())[0]
    ex = df[(df["N_actual"] == first_N) & (df["seed"] == first_seed)]
    colors = {"negative": "tab:blue", "neutral": "tab:gray", "positive": "tab:red"}
    for metric in ["kappa_over_measure_rms_r2", "kappa_over_l2"]:
        plt.figure()
        for label, sub in ex.groupby("stratum"):
            plt.scatter(sub["K_mid"], sub[metric], s=12, alpha=0.6, label=label, c=colors.get(label, None))
        plt.xlabel("analytic K_mid")
        plt.ylabel(metric)
        plt.title(f"v4b scatter example: N={first_N}, seed={first_seed}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(figdir / f"fig_scatter_example_{metric}.png", dpi=160)
        plt.close()

    payload = {
        "experiment": "Paper 17 variable-curvature OR v4b fixed physical scale",
        "args": vars(args),
        "summary_by_seed": summary.to_dict(orient="records"),
        "aggregate": agg.to_dict(orient="records"),
        "interpretation": {
            "main_change": "epsilon, measure radius, and edge length are fixed in physical units",
            "main_test": "stability of Spearman and sign-bin delta across N",
            "positive_condition": "positive_delta_fraction=1 and positive Spearman for kappa_over_measure_rms_r2",
        },
    }
    with open(outdir / "paper17_variable_curvature_or_v4b_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    lines = []
    lines.append("# Paper 17 — Variable-curvature local OR v4b")
    lines.append("")
    lines.append("Fixed physical scale test:")
    lines.append("")
    lines.append(r"\[")
    lines.append(r"\epsilon=\mathrm{const},\qquad r_{\rm meas}=\mathrm{const},\qquad \ell_{\rm edge}<\ell_{\max}.")
    lines.append(r"\]")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append("| N | metric | Spearman mean | Spearman std | slope mean | R2 mean | delta pos-neg mean | delta std | positive delta fraction | support mean |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in agg.iterrows():
        lines.append(
            f"| {int(r['N_actual'])} | {r['metric']} | {r['spearman_mean']:.6g} | {r['spearman_std']:.6g} | "
            f"{r['ols_slope_mean']:.6g} | {r['ols_r2_mean']:.6g} | "
            f"{r['delta_pos_minus_neg_mean']:.6g} | {r['delta_pos_minus_neg_std']:.6g} | "
            f"{r['positive_delta_fraction']:.3g} | {r['support_mean']:.3g} |"
        )

    (outdir / "paper17_variable_curvature_or_v4b_summary.md").write_text("\n".join(lines), encoding="utf-8")

    print("="*110)
    print("Paper 17 — Variable-curvature local OR v4b")
    print("="*110)
    print(agg.to_string(index=False))

    print("\nFiles written:")
    print(outdir / "paper17_variable_curvature_or_v4b_rows.csv")
    print(outdir / "paper17_variable_curvature_or_v4b_summary_by_seed.csv")
    print(outdir / "paper17_variable_curvature_or_v4b_aggregate.csv")
    print(outdir / "paper17_variable_curvature_or_v4b_summary.json")
    print(outdir / "paper17_variable_curvature_or_v4b_summary.md")
    print(figdir)


if __name__ == "__main__":
    main()

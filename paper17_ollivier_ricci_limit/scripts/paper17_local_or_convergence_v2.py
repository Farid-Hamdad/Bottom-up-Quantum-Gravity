#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 17 — Local Ollivier-Ricci convergence test v2

v1 showed that kappa/epsilon alone is not the right local diagnostic.
v2 tests local normalizations:
  kappa / edge_length^2
  kappa / measure_radius^2
  kappa / epsilon_N

Clean intrinsic Ricci targets:
  flat_torus2d: Ric(u,u)=0
  sphere S^2:  Ric(u,u)=1

Recommended run:
python3 papers/paper17_ollivier_ricci_limit/scripts/paper17_local_or_convergence_v2.py \
  --geometries flat_torus2d sphere \
  --N-values 128 256 512 \
  --kernel-factor 0.25 \
  --k-mode sqrt \
  --measure-k-mode sqrt \
  --idleness 0.5 \
  --n-edges 400 \
  --edge-percentile 5 \
  --n-bins 5 \
  --seed 123 \
  --output-dir papers/paper17_ollivier_ricci_limit/results/local_or_convergence_v2
"""

from __future__ import annotations

import argparse, json, math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.stats import spearmanr


def make_flat_torus2d(N: int):
    m = int(round(math.sqrt(N)))
    grid = np.arange(m) / m
    xx, yy = np.meshgrid(grid, grid, indexing="ij")
    x = np.column_stack([xx.ravel(), yy.ravel()])
    return x, {"D": 2, "N_actual": m*m, "ricci_target": 0.0, "volume": 1.0}


def make_sphere(N: int):
    i = np.arange(N)
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    z = 1.0 - 2.0 * (i + 0.5) / N
    r = np.sqrt(np.maximum(0.0, 1.0 - z*z))
    theta = 2.0 * np.pi * i / phi
    x = np.column_stack([r*np.cos(theta), r*np.sin(theta), z])
    return x, {"D": 2, "N_actual": N, "ricci_target": 1.0, "volume": 4*np.pi}


def torus_distance(x):
    dx = np.abs(x[:, None, :] - x[None, :, :])
    dx = np.minimum(dx, 1.0 - dx)
    return np.sqrt(np.sum(dx*dx, axis=-1))


def sphere_distance(x):
    dots = np.clip(x @ x.T, -1.0, 1.0)
    return np.arccos(dots)


def build_geometry(name, N):
    if name == "flat_torus2d":
        x, meta = make_flat_torus2d(N)
        return x, meta, torus_distance(x)
    if name == "sphere":
        x, meta = make_sphere(N)
        return x, meta, sphere_distance(x)
    raise ValueError(name)


def choose_k(N, mode, fixed):
    if fixed is not None:
        return int(fixed)
    if mode == "sqrt":
        return max(4, int(round(math.sqrt(N))))
    if mode == "log":
        return max(4, int(round(math.log(N)**2)))
    raise ValueError(mode)


def epsilon_from_knn(dist, k, kernel_factor):
    sd = np.sort(dist, axis=1)
    kth = sd[:, min(k, dist.shape[1]-1)]
    eps_knn = float(np.median(kth*kth))
    return eps_knn, float(kernel_factor * eps_knn)


def heat_kernel(dist, eps):
    K = np.exp(-(dist*dist)/(4.0*eps))
    np.fill_diagonal(K, 0.0)
    return K


def local_measure(i, K, dist, measure_k, idleness):
    row = K[i].copy()
    row[i] = 0.0
    if measure_k >= len(row)-1:
        neigh = np.where(row > 0)[0]
    else:
        neigh = np.argpartition(row, -measure_k)[-measure_k:]
        neigh = neigh[row[neigh] > 0]

    weights = row[neigh]
    if weights.sum() <= 0:
        return np.array([i], int), np.array([1.0], float), 0.0, 0.0

    p_neigh = (1.0-idleness) * weights / weights.sum()
    support = np.concatenate([[i], neigh])
    prob = np.concatenate([[idleness], p_neigh])
    prob = prob / prob.sum()

    d = dist[i, neigh]
    mean_r = float(np.average(d, weights=weights)) if len(d) else 0.0
    rms_r = float(np.sqrt(np.average(d*d, weights=weights))) if len(d) else 0.0
    return support.astype(int), prob.astype(float), mean_r, rms_r


def wasserstein_1(supp_a, prob_a, supp_b, prob_b, dist):
    na, nb = len(supp_a), len(supp_b)
    C = dist[np.ix_(supp_a, supp_b)]
    c = C.ravel()
    A_eq, b_eq = [], []
    for a in range(na):
        row = np.zeros(na*nb)
        row[a*nb:(a+1)*nb] = 1.0
        A_eq.append(row); b_eq.append(prob_a[a])
    for b in range(nb):
        row = np.zeros(na*nb)
        row[b::nb] = 1.0
        A_eq.append(row); b_eq.append(prob_b[b])
    res = linprog(c, A_eq=np.asarray(A_eq), b_eq=np.asarray(b_eq), bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.fun)


def sample_edges(dist, n_edges, rng, edge_percentile):
    N = dist.shape[0]
    iu = np.triu_indices(N, k=1)
    dvals = dist[iu]
    thr = np.percentile(dvals, edge_percentile)
    mask = dvals <= thr
    cand = np.column_stack([iu[0][mask], iu[1][mask]])
    cd = dvals[mask]
    if len(cand) == 0:
        raise RuntimeError("No edge candidates.")
    if len(cand) <= n_edges:
        idx = np.arange(len(cand))
    else:
        w = 1.0 / np.maximum(cd, 1e-12)
        w /= w.sum()
        idx = rng.choice(len(cand), size=n_edges, replace=False, p=w)
    return [(int(cand[t,0]), int(cand[t,1])) for t in idx]


def compute_case(geom, N_input, args, rng):
    x, meta, dist = build_geometry(geom, N_input)
    N = len(x)
    k = choose_k(N, args.k_mode, args.k_fixed)
    mk = choose_k(N, args.measure_k_mode, args.measure_k_fixed)
    eps_knn, eps = epsilon_from_knn(dist, k, args.kernel_factor)
    K = heat_kernel(dist, eps)
    edges = sample_edges(dist, args.n_edges, rng, args.edge_percentile)

    cache = {}
    rows = []
    for idx, (i, j) in enumerate(edges):
        if i not in cache:
            cache[i] = local_measure(i, K, dist, mk, args.idleness)
        if j not in cache:
            cache[j] = local_measure(j, K, dist, mk, args.idleness)
        si, pi, mean_ri, rms_ri = cache[i]
        sj, pj, mean_rj, rms_rj = cache[j]
        dij = float(dist[i, j])
        if dij <= 0:
            continue
        W1 = wasserstein_1(si, pi, sj, pj, dist)
        kappa = 1.0 - W1/dij
        mean_r = 0.5*(mean_ri+mean_rj)
        rms_r = 0.5*(rms_ri+rms_rj)
        rows.append({
            "geometry": geom, "N_input": int(N_input), "N_actual": int(N),
            "D": int(meta["D"]), "ricci_target": float(meta["ricci_target"]),
            "volume": float(meta["volume"]), "k": int(k), "measure_k": int(mk),
            "kernel_factor": float(args.kernel_factor), "epsilon_knn": float(eps_knn),
            "epsilon_N": float(eps), "idleness": float(args.idleness),
            "edge_index": int(idx), "i": int(i), "j": int(j),
            "edge_length": dij, "measure_mean_radius": mean_r, "measure_rms_radius": rms_r,
            "W1": W1, "kappa_or": kappa,
            "kappa_over_epsilon": kappa/eps,
            "kappa_over_l2": kappa/(dij*dij),
            "kappa_over_measure_mean_r2": kappa/(mean_r*mean_r) if mean_r > 0 else np.nan,
            "kappa_over_measure_rms_r2": kappa/(rms_r*rms_r) if rms_r > 0 else np.nan,
            "support_i": len(si), "support_j": len(sj),
        })
    return rows


def summarize(sub):
    out = {"n_edges": int(len(sub))}
    for col in ["kappa_or","kappa_over_epsilon","kappa_over_l2","kappa_over_measure_mean_r2","kappa_over_measure_rms_r2","edge_length","measure_rms_radius"]:
        x = pd.to_numeric(sub[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        out[f"mean_{col}"] = float(x.mean())
        out[f"median_{col}"] = float(x.median())
        out[f"std_{col}"] = float(x.std(ddof=0))
        out[f"q10_{col}"] = float(x.quantile(0.10))
        out[f"q90_{col}"] = float(x.quantile(0.90))
    out["positive_fraction"] = float((sub["kappa_or"] > 0).mean())
    return out


def safe_spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or len(np.unique(y[ok])) < 2:
        return np.nan, np.nan
    r = spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--geometries", nargs="+", default=["flat_torus2d","sphere"], choices=["flat_torus2d","sphere"])
    p.add_argument("--N-values", nargs="+", type=int, default=[128,256,512])
    p.add_argument("--kernel-factor", type=float, default=0.25)
    p.add_argument("--k-mode", default="sqrt", choices=["sqrt","log"])
    p.add_argument("--k-fixed", type=int, default=None)
    p.add_argument("--measure-k-mode", default="sqrt", choices=["sqrt","log"])
    p.add_argument("--measure-k-fixed", type=int, default=None)
    p.add_argument("--idleness", type=float, default=0.5)
    p.add_argument("--n-edges", type=int, default=400)
    p.add_argument("--edge-percentile", type=float, default=5.0)
    p.add_argument("--n-bins", type=int, default=5)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows = []
    for geom in args.geometries:
        for N in args.N_values:
            print(f"[run] geometry={geom} N={N}")
            rows.extend(compute_case(geom, N, args, rng))

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "paper17_local_or_v2_rows.csv", index=False)

    summary_rows = []
    for (geom, N), sub in df.groupby(["geometry","N_actual"]):
        s = summarize(sub)
        s.update({
            "geometry": geom, "N_actual": int(N),
            "epsilon_N": float(sub["epsilon_N"].iloc[0]),
            "k": int(sub["k"].iloc[0]), "measure_k": int(sub["measure_k"].iloc[0]),
            "idleness": float(sub["idleness"].iloc[0]),
            "ricci_target": float(sub["ricci_target"].iloc[0]),
        })
        summary_rows.append(s)
    summary = pd.DataFrame(summary_rows).sort_values(["geometry","N_actual"])
    summary.to_csv(outdir / "paper17_local_or_v2_summary.csv", index=False)

    bin_rows = []
    if {"flat_torus2d","sphere"}.issubset(set(df["geometry"])):
        for N in sorted(set(df["N_actual"])):
            dN = df[df["N_actual"] == N].copy()
            if set(dN["geometry"]) != {"flat_torus2d","sphere"}:
                continue
            lo, hi = dN["edge_length"].quantile(0.02), dN["edge_length"].quantile(0.98)
            bins = np.linspace(lo, hi, args.n_bins+1)
            dN["length_bin"] = pd.cut(dN["edge_length"], bins=bins, include_lowest=True)
            for b, sb in dN.groupby("length_bin", observed=False):
                flat = sb[sb["geometry"]=="flat_torus2d"]
                sph = sb[sb["geometry"]=="sphere"]
                if len(flat) < 3 or len(sph) < 3:
                    continue
                row = {"N_actual": int(N), "length_bin": str(b), "bin_edge_length_mean": float(sb["edge_length"].mean()), "flat_n": int(len(flat)), "sphere_n": int(len(sph))}
                for metric in ["kappa_over_epsilon","kappa_over_l2","kappa_over_measure_mean_r2","kappa_over_measure_rms_r2"]:
                    row[f"flat_mean_{metric}"] = float(flat[metric].mean())
                    row[f"sphere_mean_{metric}"] = float(sph[metric].mean())
                    row[f"delta_sphere_minus_flat_{metric}"] = float(sph[metric].mean()-flat[metric].mean())
                bin_rows.append(row)
    bins_df = pd.DataFrame(bin_rows)
    bins_df.to_csv(outdir / "paper17_local_or_v2_binned_delta.csv", index=False)

    discr_rows = []
    for metric in ["kappa_over_epsilon","kappa_over_l2","kappa_over_measure_mean_r2","kappa_over_measure_rms_r2"]:
        r, p = safe_spearman(df[metric].replace([np.inf,-np.inf], np.nan), df["ricci_target"])
        discr_rows.append({"metric": metric, "spearman_vs_ricci_target": r, "p_value": p})
    discr = pd.DataFrame(discr_rows)
    discr.to_csv(outdir / "paper17_local_or_v2_discrimination.csv", index=False)

    for metric in ["kappa_over_epsilon","kappa_over_l2","kappa_over_measure_rms_r2"]:
        plt.figure()
        for geom, sub in df.groupby("geometry"):
            vals = sub[metric].replace([np.inf,-np.inf], np.nan).dropna()
            plt.hist(vals, bins=40, alpha=0.45, label=geom)
        plt.xlabel(metric); plt.ylabel("count"); plt.title(f"Local OR metric: {metric}")
        plt.legend(); plt.tight_layout()
        plt.savefig(figdir / f"fig_hist_{metric}.png", dpi=160); plt.close()

    if len(bins_df):
        for metric in ["kappa_over_l2","kappa_over_measure_rms_r2"]:
            plt.figure()
            for N, sub in bins_df.groupby("N_actual"):
                plt.plot(sub["bin_edge_length_mean"], sub[f"delta_sphere_minus_flat_{metric}"], marker="o", label=f"N={N}")
            plt.axhline(0, linestyle="--")
            plt.xlabel("mean edge length bin"); plt.ylabel(f"sphere-flat delta: {metric}")
            plt.title(f"Binned sphere-minus-flat: {metric}")
            plt.legend(); plt.tight_layout()
            plt.savefig(figdir / f"fig_binned_delta_{metric}.png", dpi=160); plt.close()

    payload = {
        "experiment": "Paper 17 local OR convergence v2",
        "args": vars(args),
        "summary": summary.to_dict(orient="records"),
        "binned_delta": bins_df.to_dict(orient="records") if len(bins_df) else [],
        "discrimination": discr.to_dict(orient="records"),
    }
    with open(outdir / "paper17_local_or_v2_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    lines = []
    lines.append("# Paper 17 — Local OR convergence v2")
    lines.append("")
    lines.append("v2 tests local normalizations: kappa/d^2 and kappa/measure_radius^2.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("| geometry | N | eps | edges | mean kappa/eps | mean kappa/l2 | mean kappa/rms_r2 | positive fraction |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in summary.iterrows():
        lines.append(f"| {r['geometry']} | {int(r['N_actual'])} | {r['epsilon_N']:.6g} | {int(r['n_edges'])} | {r['mean_kappa_over_epsilon']:.6g} | {r['mean_kappa_over_l2']:.6g} | {r['mean_kappa_over_measure_rms_r2']:.6g} | {r['positive_fraction']:.6g} |")
    lines.append("")
    lines.append("## Discrimination")
    lines.append("")
    lines.append("| metric | Spearman vs Ricci target | p-value |")
    lines.append("|---|---:|---:|")
    for _, r in discr.iterrows():
        lines.append(f"| {r['metric']} | {r['spearman_vs_ricci_target']:.6g} | {r['p_value']:.3e} |")
    (outdir / "paper17_local_or_v2_summary.md").write_text("\n".join(lines), encoding="utf-8")

    print("="*110)
    print("Paper 17 — Local Ollivier-Ricci convergence v2")
    print("="*110)
    print(summary.to_string(index=False))
    print("\n" + "="*110)
    print("EDGE-LEVEL DISCRIMINATION")
    print("="*110)
    print(discr.to_string(index=False))
    if len(bins_df):
        print("\n" + "="*110)
        print("BINNED SPHERE MINUS FLAT")
        print("="*110)
        print(bins_df.to_string(index=False))
    print("\nFiles written:")
    print(outdir / "paper17_local_or_v2_rows.csv")
    print(outdir / "paper17_local_or_v2_summary.csv")
    print(outdir / "paper17_local_or_v2_binned_delta.csv")
    print(outdir / "paper17_local_or_v2_discrimination.csv")
    print(outdir / "paper17_local_or_v2_summary.json")
    print(outdir / "paper17_local_or_v2_summary.md")
    print(figdir)


if __name__ == "__main__":
    main()

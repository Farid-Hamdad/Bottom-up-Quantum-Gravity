#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP gravitational-wave mass-gap spectral-origin test v15.

Question:
    Does the fitted edge-Hessian mass gap m_eff^2 correlate with the smallest
    non-zero eigenvalue lambda_1(L_ent) of the entanglement Laplacian?

This is a post-processing script. It reads one or more mass-gap summary CSVs
from v13/v14, reconstructs the entanglement graph W_ij for each row, computes
several Laplacian gaps, and fits/correlates them against m_eff^2.

It computes:
    - lambda1_comb       : second eigenvalue of combinatorial L = D - W
    - lambda1_norm       : second eigenvalue of normalized L = I - D^{-1/2} W D^{-1/2}
    - lambda1_rw         : second eigenvalue of random-walk equivalent spectrum
    - lambda1_comb_phys  : lambda1_comb / a^2, optional physical scaling

Main fits:
    m_eff^2 = alpha + beta * lambda1_*

Outputs:
    results/v15_laplacian_gap_table.csv
    results/v15_laplacian_gap_fits.json
    results/v15_laplacian_gap_summary.md
    figures/fig_v15_meff2_vs_lambda1_*.png

Example:
python3 correlate_mass_gap_laplacian_v15.py \
  --summary-csv papers/paper22_bup_gravitational_waves/results/mass_gap_v14_fixed_density/results/mass_gap_v14_summary.csv \
  --output-dir papers/paper22_bup_gravitational_waves/results/mass_gap_v15_laplacian_corr
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh


def build_grid(n_side: int, extent: float) -> np.ndarray:
    xs = np.linspace(-extent, extent, n_side)
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def pairwise_dist(coords: np.ndarray) -> np.ndarray:
    d = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(d * d, axis=-1))


def build_W(coords: np.ndarray, ell: float, cutoff: float) -> np.ndarray:
    R = pairwise_dist(coords)
    W = np.exp(-R / ell)
    W[R > cutoff] = 0.0
    np.fill_diagonal(W, 0.0)
    return 0.5 * (W + W.T)


def smallest_nonzero(vals: np.ndarray, tol: float = 1e-9) -> float:
    vals = np.sort(np.real(vals[np.isfinite(vals)]))
    nz = vals[vals > tol]
    return float(nz[0]) if len(nz) else float("nan")


def laplacian_gaps(W: np.ndarray, dense_threshold: int = 900) -> Dict[str, float]:
    n = W.shape[0]
    deg = W.sum(axis=1)
    L = np.diag(deg) - W

    # Combinatorial lambda_1
    if n <= dense_threshold:
        vals_comb = eigh(L, eigvals_only=True)
    else:
        vals_comb = eigsh(sparse.csr_matrix(L), k=8, which="SM", return_eigenvectors=False, tol=1e-10)
    lam_comb = smallest_nonzero(vals_comb)

    # Normalized Laplacian. Isolated nodes are handled by setting invsqrt=0.
    invsqrt = np.zeros_like(deg)
    ok = deg > 1e-15
    invsqrt[ok] = 1.0 / np.sqrt(deg[ok])
    Lnorm = np.eye(n) - (invsqrt[:, None] * W * invsqrt[None, :])
    # For isolated nodes, normalized diagonal convention is 0 rather than 1.
    if np.any(~ok):
        Lnorm[~ok, ~ok] = 0.0
    if n <= dense_threshold:
        vals_norm = eigh(Lnorm, eigvals_only=True)
    else:
        vals_norm = eigsh(sparse.csr_matrix(Lnorm), k=8, which="SM", return_eigenvectors=False, tol=1e-10)
    lam_norm = smallest_nonzero(vals_norm)

    # Random-walk Laplacian has same nonzero spectrum as normalized L for undirected graphs.
    # We still report it explicitly for naming clarity.
    return {
        "lambda1_comb": lam_comb,
        "lambda1_norm": lam_norm,
        "lambda1_rw": lam_norm,
        "degree_min": float(np.min(deg)),
        "degree_mean": float(np.mean(deg)),
        "degree_max": float(np.max(deg)),
    }


def linfit(x: np.ndarray, y: np.ndarray) -> Dict[str, object]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return {"ok": False, "n_used": int(ok.sum())}
    X = np.column_stack([np.ones(ok.sum()), x[ok]])
    coeff, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    pred = X @ coeff
    resid = y[ok] - pred
    ss_res = float(np.sum(resid * resid))
    ss_tot = float(np.sum((y[ok] - y[ok].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rmse = float(np.sqrt(np.mean(resid * resid)))
    corr = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() >= 3 else float("nan")
    return {
        "ok": True,
        "n_used": int(ok.sum()),
        "alpha": float(coeff[0]),
        "beta": float(coeff[1]),
        "r2": float(r2),
        "rmse": rmse,
        "pearson_corr": corr,
    }


def plot_relation(df: pd.DataFrame, xkey: str, outpath: Path, title: str):
    x = df[xkey].to_numpy(float)
    y = df["meff2"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    plt.figure(figsize=(7.0, 4.6))
    if ok.sum():
        plt.scatter(x[ok], y[ok])
    if ok.sum() >= 3:
        fit = linfit(x, y)
        xs = np.linspace(np.min(x[ok]), np.max(x[ok]), 200)
        ys = fit["alpha"] + fit["beta"] * xs
        plt.plot(xs, ys, linestyle="--", label=f"R²={fit['r2']:.3f}, r={fit['pearson_corr']:.3f}")
        plt.legend()
    plt.xlabel(xkey)
    plt.ylabel(r"$m_{\rm eff}^2$")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def infer_n_side(row: pd.Series) -> int:
    if "N_side" in row and pd.notna(row["N_side"]):
        return int(row["N_side"])
    n_nodes = int(row["N_nodes"])
    n_side = int(round(math.sqrt(n_nodes)))
    if n_side * n_side != n_nodes:
        raise ValueError(f"Cannot infer N_side from N_nodes={n_nodes}")
    return n_side


def read_summaries(paths: List[str]) -> pd.DataFrame:
    frames = []
    for p in paths:
        path = Path(p)
        df = pd.read_csv(path)
        df["source_csv"] = str(path)
        frames.append(df)
    if not frames:
        raise ValueError("No summary CSV provided")
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-csv", action="append", required=True,
                    help="Mass-gap summary CSV from v13/v14. Can be passed multiple times.")
    ap.add_argument("--ell", type=float, default=0.22)
    ap.add_argument("--cutoff", type=float, default=0.45)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    resdir = outdir / "results"
    figdir = outdir / "figures"
    resdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    df0 = read_summaries(args.summary_csv)
    rows = []
    print("=" * 100)
    print("Mass-gap v15 Laplacian-correlation test")
    print("=" * 100)
    for idx, row in df0.iterrows():
        n_side = infer_n_side(row)
        extent = float(row["extent"])
        a = float(row["a"]) if "a" in row and pd.notna(row["a"]) else 2.0 * extent / (n_side - 1)
        Lphys = float(row["L"]) if "L" in row and pd.notna(row["L"]) else 2.0 * extent
        coords = build_grid(n_side, extent)
        W = build_W(coords, args.ell, args.cutoff)
        gaps = laplacian_gaps(W)
        out = row.to_dict()
        out.update({
            "N_side": n_side,
            "extent": extent,
            "a": a,
            "L": Lphys,
            "ell": args.ell,
            "cutoff": args.cutoff,
            **gaps,
            "lambda1_comb_phys": gaps["lambda1_comb"] / (a * a) if a > 0 else np.nan,
            "lambda1_norm_phys": gaps["lambda1_norm"] / (a * a) if a > 0 else np.nan,
        })
        rows.append(out)
        print(f"N={int(row['N_nodes']):4d} extent={extent:.6g} a={a:.6g} L={Lphys:.6g} "
              f"meff2={float(row['meff2']):.6g} "
              f"lambda1_comb={gaps['lambda1_comb']:.6g} lambda1_norm={gaps['lambda1_norm']:.6g}")

    df = pd.DataFrame(rows)
    out_csv = resdir / "v15_laplacian_gap_table.csv"
    df.to_csv(out_csv, index=False)

    keys = ["lambda1_comb", "lambda1_norm", "lambda1_rw", "lambda1_comb_phys", "lambda1_norm_phys"]
    fits = {key: linfit(df[key].to_numpy(float), df["meff2"].to_numpy(float)) for key in keys}

    # Also fit per family if available.
    if "family" in df.columns:
        fits_by_family = {}
        for fam, g in df.groupby("family"):
            fits_by_family[str(fam)] = {key: linfit(g[key].to_numpy(float), g["meff2"].to_numpy(float)) for key in keys}
        fits["by_family"] = fits_by_family

    with open(resdir / "v15_laplacian_gap_fits.json", "w", encoding="utf-8") as f:
        json.dump(fits, f, indent=2)

    for key in keys:
        plot_relation(df, key, figdir / f"fig_v15_meff2_vs_{key}.png", f"Mass gap vs {key}")

    # Markdown summary.
    best_key = None
    best_r2 = -np.inf
    for key in keys:
        fit = fits.get(key, {})
        if isinstance(fit, dict) and fit.get("ok") and fit.get("r2", -np.inf) > best_r2:
            best_key = key
            best_r2 = fit["r2"]
    md = []
    md.append("# v15 — Mass gap vs entanglement-Laplacian spectral gap")
    md.append("")
    md.append(f"Input CSVs: {', '.join(args.summary_csv)}")
    md.append("")
    md.append("## Best linear relation")
    if best_key is not None:
        fit = fits[best_key]
        md.append(f"Best predictor: `{best_key}`")
        md.append(f"- alpha = {fit['alpha']:.12g}")
        md.append(f"- beta = {fit['beta']:.12g}")
        md.append(f"- R2 = {fit['r2']:.6f}")
        md.append(f"- Pearson r = {fit['pearson_corr']:.6f}")
        md.append(f"- RMSE = {fit['rmse']:.6g}")
    md.append("")
    md.append("## Interpretation")
    md.append("If `meff2` is tightly correlated with `lambda1_*`, the fitted mass gap has a direct spectral origin in the BuP entanglement Laplacian.")
    md.append("The strongest claim should use the predictor with the highest R2 and stable sign across families.")
    (resdir / "v15_laplacian_gap_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("\nFits:")
    print(json.dumps(fits, indent=2))
    print(f"\nWrote outputs to: {outdir}")


if __name__ == "__main__":
    main()

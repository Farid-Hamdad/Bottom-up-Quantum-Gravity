#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP gravitational-wave mass-gap relative spectral-gap test v16.

Purpose
-------
Post-process the v15 table and test whether the residual edge-Hessian mass gap
is better described as a two-parameter spectral law

    m_eff^2 = m0^2 - beta * lambda1

or as a constrained relative-gap law

    m_eff^2 = beta * (1 - lambda1)

and its generalized form

    m_eff^2 = beta * (1 - lambda1 / lambda_star).

For the generalized model with free lambda_star, the fit is algebraically
identical to the free linear model, with

    lambda_star = m0^2 / beta.

The main theoretical question is therefore whether lambda_star is compatible
with 1, equivalently whether m0^2 ~= beta.

Inputs
------
- A v15 table CSV, typically:
  results/mass_gap_v15_laplacian_corr_combined/results/v15_laplacian_gap_table.csv

Outputs
-------
results/v16_relative_gap_model_comparison.csv
results/v16_relative_gap_fits.json
results/v16_relative_gap_summary.md
figures/fig_v16_*.png

Example
-------
python3 test_relative_spectral_gap_v16.py \
  --v15-table papers/paper22_bup_gravitational_waves/results/mass_gap_v15_laplacian_corr_combined/results/v15_laplacian_gap_table.csv \
  --lambda-key lambda1_norm \
  --output-dir papers/paper22_bup_gravitational_waves/results/mass_gap_v16_relative_spectral
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def finite_xy(df: pd.DataFrame, xkey: str, ykey: str = "meff2") -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    sub = df[[xkey, ykey] + (["family"] if "family" in df.columns else [])].copy()
    sub[xkey] = pd.to_numeric(sub[xkey], errors="coerce")
    sub[ykey] = pd.to_numeric(sub[ykey], errors="coerce")
    sub = sub[np.isfinite(sub[xkey]) & np.isfinite(sub[ykey])]
    return sub[xkey].to_numpy(float), sub[ykey].to_numpy(float), sub


def metrics(y: np.ndarray, yhat: np.ndarray, k: int) -> Dict[str, float]:
    n = len(y)
    resid = y - yhat
    sse = float(np.sum(resid * resid))
    rmse = float(math.sqrt(sse / max(n, 1)))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - sse / sst) if sst > 0 else float("nan")
    # Gaussian likelihood criteria up to additive constants.
    eps = 1e-300
    aic = float(n * math.log(max(sse / max(n, 1), eps)) + 2 * k) if n > 0 else float("nan")
    bic = float(n * math.log(max(sse / max(n, 1), eps)) + k * math.log(max(n, 1))) if n > 0 else float("nan")
    return {"n": n, "k_params": k, "sse": sse, "rmse": rmse, "r2": r2, "aic": aic, "bic": bic}


def fit_free_linear(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    # y = alpha + slope*x = m0_squared - beta_positive*x
    X = np.column_stack([np.ones_like(x), x])
    coeff, *_ = np.linalg.lstsq(X, y, rcond=None)
    alpha = float(coeff[0])
    slope = float(coeff[1])
    yhat = X @ coeff
    beta_positive = float(-slope)
    lam_star = float(alpha / beta_positive) if beta_positive != 0 else float("nan")
    out = {
        "ok": True,
        "model": "free_linear",
        "formula": "m2 = m0_squared - beta * lambda1",
        "m0_squared": alpha,
        "slope": slope,
        "beta": beta_positive,
        "lambda_star_free": lam_star,
        "m0_over_beta": lam_star,
    }
    out.update(metrics(y, yhat, k=2))
    return out


def fit_constrained_one_minus(x: np.ndarray, y: np.ndarray, lambda_star: float = 1.0) -> Dict[str, float]:
    # y = beta * (1 - x/lambda_star), one parameter beta.
    z = 1.0 - x / float(lambda_star)
    denom = float(np.dot(z, z))
    beta = float(np.dot(z, y) / denom) if denom > 0 else float("nan")
    yhat = beta * z
    out = {
        "ok": True,
        "model": f"relative_gap_lambda_star_{lambda_star:g}",
        "formula": "m2 = beta * (1 - lambda1/lambda_star)",
        "lambda_star": float(lambda_star),
        "beta": beta,
        "m0_squared_constrained": beta,
        "m0_over_beta": 1.0,
    }
    out.update(metrics(y, yhat, k=1))
    return out


def fit_zero_intercept_lambda(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    # Naive direct proportionality y = gamma * x, included as a negative control.
    denom = float(np.dot(x, x))
    gamma = float(np.dot(x, y) / denom) if denom > 0 else float("nan")
    yhat = gamma * x
    out = {"ok": True, "model": "direct_proportional", "formula": "m2 = gamma * lambda1", "gamma": gamma}
    out.update(metrics(y, yhat, k=1))
    return out


def bootstrap_free_ratio(x: np.ndarray, y: np.ndarray, n_boot: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(y)
    vals = []
    m0s = []
    betas = []
    if n < 3 or n_boot <= 0:
        return {"ok": False, "reason": "not enough points or disabled"}
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        xb = x[idx]
        yb = y[idx]
        try:
            fit = fit_free_linear(xb, yb)
            ratio = fit.get("m0_over_beta", float("nan"))
            beta = fit.get("beta", float("nan"))
            m0 = fit.get("m0_squared", float("nan"))
            if np.isfinite(ratio) and np.isfinite(beta) and beta > 0:
                vals.append(ratio)
                m0s.append(m0)
                betas.append(beta)
        except Exception:
            pass
    if not vals:
        return {"ok": False, "reason": "all bootstrap fits failed"}
    arr = np.array(vals, dtype=float)
    return {
        "ok": True,
        "n_boot_requested": int(n_boot),
        "n_boot_used": int(len(arr)),
        "ratio_mean": float(np.mean(arr)),
        "ratio_median": float(np.median(arr)),
        "ratio_p05": float(np.percentile(arr, 5)),
        "ratio_p16": float(np.percentile(arr, 16)),
        "ratio_p84": float(np.percentile(arr, 84)),
        "ratio_p95": float(np.percentile(arr, 95)),
        "m0_mean": float(np.mean(m0s)),
        "beta_mean": float(np.mean(betas)),
        "frac_ratio_contains_1_rough": float(np.mean(np.abs(arr - 1.0) < 0.1)),
    }


def fit_suite(df: pd.DataFrame, xkey: str, lambda_stars: List[float], bootstrap: int, seed: int) -> Dict[str, object]:
    x, y, sub = finite_xy(df, xkey)
    if len(y) < 3:
        return {"ok": False, "reason": f"not enough finite points for {xkey}", "n": int(len(y))}

    models: List[Dict[str, float]] = []
    models.append(fit_free_linear(x, y))
    models.append(fit_constrained_one_minus(x, y, lambda_star=1.0))
    for ls in lambda_stars:
        if abs(ls - 1.0) > 1e-12:
            models.append(fit_constrained_one_minus(x, y, lambda_star=ls))
    models.append(fit_zero_intercept_lambda(x, y))

    # Sort by AIC.
    models_sorted = sorted(models, key=lambda d: d.get("aic", float("inf")))
    best_aic = models_sorted[0]["aic"]
    best_bic = min(m["bic"] for m in models_sorted)
    for m in models_sorted:
        m["delta_aic"] = float(m["aic"] - best_aic)
        m["delta_bic"] = float(m["bic"] - best_bic)

    boot = bootstrap_free_ratio(x, y, bootstrap, seed)
    return {
        "ok": True,
        "lambda_key": xkey,
        "n_used": int(len(y)),
        "models_by_aic": models_sorted,
        "free_ratio_bootstrap": boot,
    }


def plot_models(df: pd.DataFrame, xkey: str, result: Dict[str, object], outpath: Path, title: str) -> None:
    x, y, sub = finite_xy(df, xkey)
    if len(y) < 3 or not result.get("ok"):
        return
    xs = np.linspace(float(np.min(x)), float(np.max(x)), 250)
    plt.figure(figsize=(7, 5))
    if "family" in sub.columns:
        for fam, g in sub.groupby("family"):
            plt.scatter(g[xkey], g["meff2"], label=str(fam))
    else:
        plt.scatter(x, y, label="data")

    models = result["models_by_aic"]
    for m in models[:3]:
        name = m["model"]
        if name == "free_linear":
            ys = m["m0_squared"] - m["beta"] * xs
            label = f"free: m0/beta={m['m0_over_beta']:.3g}, R2={m['r2']:.3f}"
            plt.plot(xs, ys, label=label)
        elif name.startswith("relative_gap"):
            lam_star = m["lambda_star"]
            ys = m["beta"] * (1.0 - xs / lam_star)
            label = f"rel λ*={lam_star:.3g}: R2={m['r2']:.3f}"
            plt.plot(xs, ys, linestyle="--", label=label)
        elif name == "direct_proportional":
            ys = m["gamma"] * xs
            plt.plot(xs, ys, linestyle=":", label=f"direct: R2={m['r2']:.3f}")
    plt.xlabel(xkey)
    plt.ylabel("m_eff²")
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description="v16 relative spectral-gap test for BuP edge-Hessian mass gap")
    p.add_argument("--v15-table", required=True, help="Path to v15_laplacian_gap_table.csv")
    p.add_argument("--lambda-key", action="append", default=None,
                   help="Lambda column to test. Can be repeated. Default: lambda1_norm and lambda1_comb")
    p.add_argument("--lambda-star-list", default="1,2",
                   help="Fixed lambda_star values for constrained model m2=beta*(1-lambda/lambda_star). Default: 1,2")
    p.add_argument("--bootstrap", type=int, default=2000, help="Bootstrap count for m0/beta ratio. Default: 2000")
    p.add_argument("--seed", type=int, default=123, help="Random seed")
    p.add_argument("--output-dir", required=True, help="Output directory")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    resdir = outdir / "results"
    figdir = outdir / "figures"
    resdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.v15_table)
    lambda_keys = args.lambda_key or ["lambda1_norm", "lambda1_comb", "lambda1_norm_phys", "lambda1_comb_phys"]
    lambda_keys = [k for k in lambda_keys if k in df.columns]
    lambda_stars = [float(x.strip()) for x in args.lambda_star_list.split(",") if x.strip()]

    all_results: Dict[str, object] = {}
    comparison_rows = []

    for key in lambda_keys:
        result = fit_suite(df, key, lambda_stars=lambda_stars, bootstrap=args.bootstrap, seed=args.seed)
        all_results[key] = result
        if result.get("ok"):
            for m in result["models_by_aic"]:
                row = {"lambda_key": key}
                row.update(m)
                comparison_rows.append(row)
            plot_models(df, key, result, figdir / f"fig_v16_models_{key}.png", f"v16 relative spectral-gap test: {key}")

    # Per-family tests, useful because v15 showed strong family-specific trends.
    by_family = {}
    if "family" in df.columns:
        for fam, g in df.groupby("family"):
            by_family[str(fam)] = {}
            for key in lambda_keys:
                by_family[str(fam)][key] = fit_suite(g, key, lambda_stars=lambda_stars, bootstrap=args.bootstrap, seed=args.seed)
    all_results["by_family"] = by_family

    comp = pd.DataFrame(comparison_rows)
    comp.to_csv(resdir / "v16_relative_gap_model_comparison.csv", index=False)
    with open(resdir / "v16_relative_gap_fits.json", "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2)

    # Human summary.
    md = []
    md.append("# v16 — Relative spectral-gap test")
    md.append("")
    md.append(f"Input table: `{args.v15_table}`")
    md.append("")
    md.append("The tested theoretical constraint is:")
    md.append("")
    md.append("```math")
    md.append("m_{\\rm eff}^2 = \\beta(1-\\lambda_1)")
    md.append("```")
    md.append("")
    md.append("compared with the free two-parameter law:")
    md.append("")
    md.append("```math")
    md.append("m_{\\rm eff}^2 = m_0^2 - \\beta\\lambda_1.")
    md.append("```")
    md.append("")
    for key in lambda_keys:
        r = all_results.get(key, {})
        if not isinstance(r, dict) or not r.get("ok"):
            continue
        free = next((m for m in r["models_by_aic"] if m["model"] == "free_linear"), None)
        rel1 = next((m for m in r["models_by_aic"] if m["model"] == "relative_gap_lambda_star_1"), None)
        best = r["models_by_aic"][0]
        md.append(f"## Predictor `{key}`")
        md.append("")
        if free:
            md.append(f"Free model: m0² = `{free['m0_squared']:.12g}`, beta = `{free['beta']:.12g}`, m0²/beta = `{free['m0_over_beta']:.6g}`, R² = `{free['r2']:.6f}`, RMSE = `{free['rmse']:.6g}`.")
        if rel1:
            md.append(f"Constrained λ*=1 model: beta = `{rel1['beta']:.12g}`, R² = `{rel1['r2']:.6f}`, RMSE = `{rel1['rmse']:.6g}`, ΔAIC = `{rel1['delta_aic']:.6g}`.")
        md.append(f"Best by AIC: `{best['model']}` with AIC `{best['aic']:.6g}` and BIC `{best['bic']:.6g}`.")
        boot = r.get("free_ratio_bootstrap", {})
        if isinstance(boot, dict) and boot.get("ok"):
            md.append(f"Bootstrap m0²/beta median = `{boot['ratio_median']:.6g}`; 16–84% = `[{boot['ratio_p16']:.6g}, {boot['ratio_p84']:.6g}]`; 5–95% = `[{boot['ratio_p05']:.6g}, {boot['ratio_p95']:.6g}]`.")
        md.append("")
    md.append("## Reading rule")
    md.append("")
    md.append("If the constrained λ*=1 model has nearly the same RMSE/AIC as the free model and the bootstrap interval for m0²/beta contains 1, then the relation m0²≈beta is supported. If not, the data support the generalized relative law with λ*=m0²/beta rather than exact λ*=1.")
    (resdir / "v16_relative_gap_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Mass-gap v16 relative spectral-gap test")
    print("=" * 100)
    if len(comp):
        cols = [c for c in ["lambda_key", "model", "r2", "rmse", "aic", "bic", "delta_aic", "m0_squared", "beta", "m0_over_beta", "lambda_star"] if c in comp.columns]
        print(comp[cols].to_string(index=False))
    print("\nWrote outputs to:", outdir)


if __name__ == "__main__":
    main()

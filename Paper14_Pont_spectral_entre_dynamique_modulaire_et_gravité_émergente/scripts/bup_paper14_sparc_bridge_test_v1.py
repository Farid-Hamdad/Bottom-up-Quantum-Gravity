#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP Paper 14 — SPARC bridge test v1

Goal
----
Test the Paper 14 bridge on SPARC-like galaxy rotation files.

For each galaxy:
1. Read SPARC rotmod file.
2. Reconstruct a baryonic surface-density proxy Sigma(R).
3. Build a 2D disk entanglement graph W_ij.
4. Compute L_ent = D - W and graph invariants.
5. Compute beta_grav, beta_mod_bridge, beta_mod_ML.
6. Compare bridge prediction and ML prediction.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.optimize import curve_fit
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

EPS = 1e-12


def parse_args():
    p = argparse.ArgumentParser(description="BuP Paper 14 SPARC bridge test v1")
    p.add_argument("--rotmod", nargs="*", default=None, help="One or more SPARC *_rotmod.dat files.")
    p.add_argument("--rotmod-dir", default=None, help="Directory containing SPARC *_rotmod.dat files.")
    p.add_argument("--training-csv", default=None, help="Paper 14 training CSV for ML beta_mod prediction.")
    p.add_argument("--candidate", default="schur")
    p.add_argument("--family", default="pos")
    p.add_argument("--min-fit-r2", type=float, default=0.90)
    p.add_argument("--rmt-min", type=float, default=None)
    p.add_argument("--rmt-max", type=float, default=None)
    p.add_argument("--n-rings", type=int, default=12)
    p.add_argument("--n-theta", type=int, default=12)
    p.add_argument("--lambda-corr", type=float, default=None)
    p.add_argument("--lambda-corr-factor", type=float, default=1.0)
    p.add_argument("--sigma-source", default="auto", choices=["auto", "sb", "vbar"])
    p.add_argument("--heat-t-min", type=float, default=0.05)
    p.add_argument("--heat-t-max", type=float, default=20.0)
    p.add_argument("--heat-n", type=int, default=60)
    p.add_argument("--walk-t-min", type=float, default=0.05)
    p.add_argument("--walk-t-max", type=float, default=20.0)
    p.add_argument("--walk-n", type=int, default=60)
    p.add_argument("--max-files", type=int, default=None)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def galaxy_name_from_path(path: Path) -> str:
    name = path.name
    name = re.sub(r"_rotmod\.dat$", "", name)
    name = re.sub(r"\.dat$", "", name)
    return name


def read_rotmod(path: Path) -> pd.DataFrame:
    arr = np.genfromtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    names = ["R", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    cols = {names[i]: arr[:, i] for i in range(min(arr.shape[1], len(names)))}
    df = pd.DataFrame(cols).replace([np.inf, -np.inf], np.nan).dropna(subset=["R"])
    df = df[df["R"] > 0].sort_values("R").reset_index(drop=True)
    for c in names:
        if c not in df.columns:
            df[c] = np.nan
    return df


def baryonic_velocity_proxy_sigma(df: pd.DataFrame):
    R = df["R"].to_numpy(float)
    vgas = np.nan_to_num(df["Vgas"].to_numpy(float), nan=0.0)
    vdisk = np.nan_to_num(df["Vdisk"].to_numpy(float), nan=0.0)
    vbul = np.nan_to_num(df["Vbul"].to_numpy(float), nan=0.0)
    vbar2 = np.abs(vgas) * vgas + vdisk**2 + vbul**2
    vbar2 = np.maximum(vbar2, 0.0)
    return vbar2 / np.maximum(R, EPS)


def sigma_profile_from_rotmod(df: pd.DataFrame, source="auto"):
    R = df["R"].to_numpy(float)
    sb = df["SBdisk"].to_numpy(float)
    has_sb = np.isfinite(sb).sum() >= max(3, len(df)//2) and np.nanmax(sb) > 0
    if source == "sb" or (source == "auto" and has_sb):
        sigma = np.nan_to_num(df["SBdisk"].to_numpy(float), nan=0.0) + np.nan_to_num(df["SBbul"].to_numpy(float), nan=0.0)
        used = "SBdisk+SBbul"
    else:
        sigma = baryonic_velocity_proxy_sigma(df)
        used = "vbar2_over_R"
    sigma = np.maximum(sigma, 0.0)
    if np.nanmax(sigma) <= EPS:
        sigma = np.ones_like(R)
        used = "flat_fallback"
    sigma = np.maximum(sigma, 1e-6 * np.nanmax(sigma))
    return R, sigma, used


def exp_profile(R, A, Rd, C):
    return A * np.exp(-R / np.maximum(Rd, EPS)) + C


def fit_Rd(R, sigma):
    if len(R) < 4:
        return float(np.nanmedian(R)) if len(R) else 3.0
    A0 = max(float(np.nanmax(sigma) - np.nanmin(sigma)), EPS)
    Rd0 = max(float((np.nanmax(R) - np.nanmin(R)) / 3.0), EPS)
    C0 = max(float(np.nanmin(sigma)), 0.0)
    try:
        popt, _ = curve_fit(exp_profile, R, sigma, p0=[A0, Rd0, C0], bounds=([0.0, 0.05, 0.0], [np.inf, 100.0, np.inf]), maxfev=20000)
        return float(popt[1])
    except Exception:
        mask = sigma > 0
        if mask.sum() >= 3:
            coef = np.polyfit(R[mask], np.log(sigma[mask]), 1)
            if coef[0] < 0:
                return float(-1.0 / coef[0])
        return Rd0


def build_disk_graph(R_data, sigma_data, n_rings, n_theta, lambda_corr):
    r_min = max(float(np.nanmin(R_data)), 1e-3)
    r_max = float(np.nanmax(R_data))
    rings = np.linspace(r_min, r_max, n_rings)
    coords, R_nodes, sigma_nodes = [], [], []
    for r in rings:
        sig_r = float(np.interp(r, R_data, sigma_data))
        for k in range(n_theta):
            th = 2.0 * np.pi * k / n_theta
            coords.append([r * np.cos(th), r * np.sin(th)])
            R_nodes.append(r)
            sigma_nodes.append(sig_r)
    coords = np.asarray(coords, float)
    R_nodes = np.asarray(R_nodes, float)
    sigma_nodes = np.asarray(sigma_nodes, float)
    sigma_nodes = sigma_nodes / max(float(np.mean(sigma_nodes)), EPS)
    sigma_nodes = np.maximum(sigma_nodes, EPS)
    dist = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=-1))
    W = np.sqrt(np.outer(sigma_nodes, sigma_nodes)) * np.exp(-dist / max(lambda_corr, EPS))
    np.fill_diagonal(W, 0.0)
    return W, coords, R_nodes, sigma_nodes


def laplacian(W):
    W = np.array(W, float, copy=True)
    np.fill_diagonal(W, 0.0)
    return np.diag(W.sum(axis=1)) - W


def normalized_laplacian(W):
    W = np.array(W, float, copy=True)
    np.fill_diagonal(W, 0.0)
    deg = W.sum(axis=1)
    inv = np.zeros_like(deg)
    m = deg > EPS
    inv[m] = 1.0 / np.sqrt(deg[m])
    return np.eye(W.shape[0]) - inv[:, None] * W * inv[None, :]


def graph_shortest_paths_from_weights(W):
    n = W.shape[0]
    dist = np.full((n, n), np.inf, dtype=float)
    np.fill_diagonal(dist, 0.0)
    m = W > EPS
    dist[m] = 1.0 / W[m]
    np.fill_diagonal(dist, 0.0)
    for k in range(n):
        dist = np.minimum(dist, dist[:, [k]] + dist[[k], :])
    return dist


def clustering_binary(W):
    A = (W > EPS).astype(float)
    np.fill_diagonal(A, 0.0)
    vals = []
    for i in range(A.shape[0]):
        neigh = np.where(A[i] > 0)[0]
        k = len(neigh)
        if k < 2:
            vals.append(0.0)
        else:
            sub = A[np.ix_(neigh, neigh)]
            vals.append(2.0 * (sub.sum() / 2.0) / (k * (k - 1)))
    return float(np.mean(vals))


def fit_slope_loglog(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > EPS) & (y > EPS)
    x, y = x[m], y[m]
    if len(x) < 6:
        return np.nan, np.nan
    lx, ly = np.log(x), np.log(y)
    order = np.argsort(lx)
    lx, ly = lx[order], ly[order]
    a, b = int(0.2 * len(lx)), int(0.8 * len(lx))
    if b <= a + 3:
        a, b = 0, len(lx)
    coef = np.polyfit(lx[a:b], ly[a:b], 1)
    pred = np.polyval(coef, lx[a:b])
    ss_res = np.sum((ly[a:b] - pred) ** 2)
    ss_tot = np.sum((ly[a:b] - ly[a:b].mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > EPS else np.nan
    return float(coef[0]), float(r2)


def spectral_dimension(Ln, t_min, t_max, n_t):
    eig = np.maximum(np.linalg.eigvalsh(Ln), 0.0)
    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)
    Z = np.array([np.sum(np.exp(-tt * eig)) for tt in t]) / len(eig)
    slope, r2 = fit_slope_loglog(t, Z)
    return float(-2.0 * slope) if np.isfinite(slope) else np.nan, float(r2)


def walk_dimension(Ln, dist, t_min, t_max, n_t):
    eig, U = eigh(Ln)
    eig = np.maximum(eig, 0.0)
    t = np.logspace(np.log10(t_min), np.log10(t_max), n_t)
    D2 = dist**2
    finite = np.isfinite(D2)
    if not np.all(finite):
        max_f = np.nanmax(D2[finite]) if finite.any() else 1.0
        D2 = np.where(finite, D2, max_f)
    msd = []
    for tt in t:
        K = (U * np.exp(-tt * eig)) @ U.T
        K = np.maximum(K, 0.0)
        P = K / np.maximum(K.sum(axis=1, keepdims=True), EPS)
        msd.append(float(np.mean(np.sum(P * D2, axis=1))))
    slope, r2 = fit_slope_loglog(t, np.asarray(msd))
    dw = 2.0 / slope if np.isfinite(slope) and abs(slope) > EPS else np.nan
    return float(dw), float(r2)


def delta3_proxy_from_spectrum(vals):
    x = np.sort(np.asarray(vals, float))
    if len(x) < 5 or np.std(x) < EPS:
        return np.nan
    u = (x - x.min()) / (x.max() - x.min() + EPS)
    N = np.arange(1, len(u)+1, dtype=float) / len(u)
    X = np.vstack([u, np.ones_like(u)]).T
    coef, *_ = np.linalg.lstsq(X, N, rcond=None)
    fit = X @ coef
    return float(np.mean((N - fit)**2))


def graph_invariants(W, args):
    W = np.array(W, float, copy=True)
    np.fill_diagonal(W, 0.0)
    n = W.shape[0]
    E = float(np.count_nonzero(np.triu(W > EPS, 1)))
    rho = 2.0 * E / max(n * (n - 1), 1)
    deg_w = W.sum(axis=1)
    deg_b = (W > EPS).sum(axis=1)
    L = laplacian(W)
    Ln = normalized_laplacian(W)
    eigLn = np.sort(np.maximum(np.linalg.eigvalsh(Ln), 0.0))
    eigL = np.sort(np.maximum(np.linalg.eigvalsh(L), 0.0))
    dist = graph_shortest_paths_from_weights(W)
    finite = dist[np.isfinite(dist) & (dist > 0)]
    ds, ds_r2 = spectral_dimension(Ln, args.heat_t_min, args.heat_t_max, args.heat_n)
    dw, dw_r2 = walk_dimension(Ln, dist, args.walk_t_min, args.walk_t_max, args.walk_n)
    alpha_eff = float(2.0 * ds / dw + dw - 4.0) if np.isfinite(ds) and np.isfinite(dw) and abs(dw) > EPS else np.nan
    beta_grav = float((alpha_eff + 1.0) / 2.0) if np.isfinite(alpha_eff) else np.nan
    return {
        "N_graph": int(n), "n_edges_input": E, "mi_rho_graph": rho,
        "mi_mean_degree_weighted": float(np.mean(deg_w)), "mi_std_degree_weighted": float(np.std(deg_w)),
        "mi_mean_degree_binary": float(np.mean(deg_b)), "mi_std_degree_binary": float(np.std(deg_b)),
        "mi_lambda2": float(eigL[1]) if len(eigL) > 1 else np.nan,
        "mi_norm_lambda2": float(eigLn[1]) if len(eigLn) > 1 else np.nan,
        "mi_clustering_binary": clustering_binary(W),
        "mi_effective_diameter_90": float(np.percentile(finite, 90)) if finite.size else np.nan,
        "mi_mean_shortest_path": float(np.mean(finite)) if finite.size else np.nan,
        "mi_laplacian_trace": float(np.trace(L)),
        "mi_ds_graph": ds, "mi_ds_fit_r2": ds_r2,
        "mi_dw_graph": dw, "mi_dw_fit_r2": dw_r2,
        "mi_ds_over_dw": float(ds / dw) if np.isfinite(ds) and np.isfinite(dw) and abs(dw) > EPS else np.nan,
        "mi_alpha_eff_graph": alpha_eff, "beta_grav": beta_grav,
        "delta3_proxy": delta3_proxy_from_spectrum(eigLn),
        "mean_MI": float(np.mean(W[np.triu_indices(n, 1)])), "max_MI": float(np.max(W)),
    }

DEFAULT_FEATURES = [
    "mi_norm_lambda2", "mi_ds_over_dw", "delta3_proxy", "mi_ds_graph", "mi_dw_graph",
    "mi_lambda2", "mi_alpha_eff_graph", "n_edges_input", "mean_MI", "max_MI",
    "mi_rho_graph", "mi_mean_degree_weighted", "mi_mean_degree_binary",
    "mi_clustering_binary", "mi_effective_diameter_90", "mi_mean_shortest_path",
]


def train_beta_model(training_csv, args):
    if training_csv is None:
        return None, [], {}
    path = Path(training_csv)
    if not path.exists():
        raise FileNotFoundError(f"Training CSV not found: {path}")
    df = pd.read_csv(path)
    if "fit_ok" in df.columns:
        df = df[df["fit_ok"] == 1].copy()
    if "family" in df.columns:
        df = df[df["family"] == args.family].copy()
    if "candidate_L" in df.columns and args.candidate is not None:
        df = df[df["candidate_L"] == args.candidate].copy()
    if "r2" in df.columns:
        df = df[pd.to_numeric(df["r2"], errors="coerce") >= args.min_fit_r2].copy()
    if args.rmt_min is not None and "K_gap_ratio" in df.columns:
        df = df[pd.to_numeric(df["K_gap_ratio"], errors="coerce") >= args.rmt_min].copy()
    if args.rmt_max is not None and "K_gap_ratio" in df.columns:
        df = df[pd.to_numeric(df["K_gap_ratio"], errors="coerce") <= args.rmt_max].copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df[np.isfinite(pd.to_numeric(df["beta"], errors="coerce"))].copy()
    features = [f for f in DEFAULT_FEATURES if f in df.columns]
    if len(features) < 3:
        raise ValueError(f"Not enough training features found: {features}")
    X, y = df[features].copy(), pd.to_numeric(df["beta"], errors="coerce").to_numpy(float)
    model = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("gb", GradientBoostingRegressor(n_estimators=300, learning_rate=0.035, max_depth=3, random_state=42)),
    ])
    model.fit(X, y)
    pred = model.predict(X)
    return model, features, {
        "training_csv": str(path), "n_train": int(len(df)), "features": features,
        "train_r2_in_sample": float(r2_score(y, pred)) if len(df) > 2 else np.nan,
        "target_beta_mean": float(np.mean(y)), "target_beta_std": float(np.std(y)),
    }


def fallback_beta_mod_ml(inv):
    lam2n = inv.get("mi_norm_lambda2", np.nan)
    if not np.isfinite(lam2n):
        return np.nan
    return float(0.31186451112106184 + 1.716023029675049 * lam2n)


def verdict(consistency):
    if not np.isfinite(consistency):
        return "insufficient_diagnostics"
    if consistency < 0.20:
        return "strong_bridge"
    if consistency < 0.35:
        return "moderate_bridge"
    return "failed_bridge"


def collect_rotmods(args):
    files = []
    if args.rotmod:
        files.extend([Path(p) for p in args.rotmod])
    if args.rotmod_dir:
        files.extend(sorted(Path(args.rotmod_dir).glob("*rotmod*.dat")))
        if not files:
            files.extend(sorted(Path(args.rotmod_dir).glob("*.dat")))
    files = [p for p in files if p.exists()]
    if args.max_files is not None:
        files = files[:args.max_files]
    return files


def main():
    args = parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = collect_rotmods(args)
    if not files:
        raise SystemExit("No rotmod files found. Use --rotmod or --rotmod-dir.")
    model, ml_features, train_info = train_beta_model(args.training_csv, args)
    print("=" * 90)
    print("BuP Paper 14 — SPARC bridge test v1")
    print("=" * 90)
    print(f"rotmod files      : {len(files)}")
    print(f"training csv      : {args.training_csv}")
    print(f"ML features       : {ml_features if ml_features else 'fallback lambda2_norm only'}")
    print(f"output            : {out}")
    rows = []
    for idx, path in enumerate(files, start=1):
        gal = galaxy_name_from_path(path)
        print(f"\n[{idx}/{len(files)}] {gal}")
        try:
            df = read_rotmod(path)
            R, sigma, sigma_used = sigma_profile_from_rotmod(df, args.sigma_source)
            Rd = fit_Rd(R, sigma)
            lambda_corr = args.lambda_corr if args.lambda_corr is not None else args.lambda_corr_factor * Rd
            lambda_corr = max(float(lambda_corr), 1e-3)
            W, coords, R_nodes, sigma_nodes = build_disk_graph(R, sigma, args.n_rings, args.n_theta, lambda_corr)
            inv = graph_invariants(W, args)
            beta_grav = inv["beta_grav"]
            beta_mod_bridge = float(0.531 + 1.726 * beta_grav) if np.isfinite(beta_grav) else np.nan
            if model is not None:
                Xrow = pd.DataFrame([{f: inv.get(f, np.nan) for f in ml_features}])
                beta_mod_ML = float(model.predict(Xrow)[0])
                ml_mode = "trained_gradient_boosting"
            else:
                beta_mod_ML = fallback_beta_mod_ml(inv)
                ml_mode = "fallback_lambda2_linear"
            consistency = abs(beta_mod_bridge - beta_mod_ML) / max(abs(beta_mod_ML), EPS) if np.isfinite(beta_mod_bridge) and np.isfinite(beta_mod_ML) else np.nan
            row = {
                "galaxy": gal, "rotmod": str(path), "n_points_rotmod": int(len(df)),
                "sigma_source": sigma_used, "Rd_fit": float(Rd), "lambda_corr": float(lambda_corr),
                "n_rings": args.n_rings, "n_theta": args.n_theta, "ml_mode": ml_mode,
                "beta_grav": beta_grav, "beta_mod_bridge": beta_mod_bridge, "beta_mod_ML": beta_mod_ML,
                "bridge_consistency": consistency, "verdict": verdict(consistency),
            }
            row.update(inv)
            rows.append(row)
            print(f"  ds={inv['mi_ds_graph']:.4f} dw={inv['mi_dw_graph']:.4f} alpha={inv['mi_alpha_eff_graph']:.4f}")
            print(f"  beta_grav={beta_grav:.4f} beta_bridge={beta_mod_bridge:.4f} beta_ML={beta_mod_ML:.4f}")
            print(f"  consistency={consistency:.3f} verdict={row['verdict']}")
            gal_dir = out / "galaxies" / gal
            gal_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1], "R": R_nodes, "sigma_node": sigma_nodes}).to_csv(gal_dir / "graph_nodes.csv", index=False)
            np.savetxt(gal_dir / "W_matrix.csv", W, delimiter=",")
        except Exception as exc:
            print("  FAILED:", exc)
            rows.append({"galaxy": gal, "rotmod": str(path), "error": str(exc), "verdict": "failed_runtime"})
    res = pd.DataFrame(rows)
    res.to_csv(out / "sparc_bridge_results.csv", index=False)
    valid = res[np.isfinite(pd.to_numeric(res.get("bridge_consistency", np.nan), errors="coerce"))].copy()
    summary = {
        "experiment": "BuP Paper 14 SPARC bridge test v1", "n_files": int(len(files)),
        "n_valid": int(len(valid)), "training": train_info,
        "parameters": {"n_rings": args.n_rings, "n_theta": args.n_theta, "lambda_corr": args.lambda_corr,
                       "lambda_corr_factor": args.lambda_corr_factor, "sigma_source": args.sigma_source},
        "counts_by_verdict": res["verdict"].value_counts(dropna=False).to_dict() if "verdict" in res else {},
    }
    if len(valid):
        summary.update({
            "mean_consistency": float(valid["bridge_consistency"].mean()),
            "median_consistency": float(valid["bridge_consistency"].median()),
            "fraction_strong": float((valid["verdict"] == "strong_bridge").mean()),
            "fraction_moderate_or_strong": float(valid["verdict"].isin(["strong_bridge", "moderate_bridge"]).mean()),
            "mean_beta_grav": float(valid["beta_grav"].mean()),
            "mean_beta_mod_bridge": float(valid["beta_mod_bridge"].mean()),
            "mean_beta_mod_ML": float(valid["beta_mod_ML"].mean()),
        })
    with open(out / "sparc_bridge_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if len(valid):
        plt.figure(figsize=(6, 6))
        plt.scatter(valid["beta_mod_ML"], valid["beta_mod_bridge"], s=45, alpha=0.8)
        mn = float(min(valid["beta_mod_ML"].min(), valid["beta_mod_bridge"].min()))
        mx = float(max(valid["beta_mod_ML"].max(), valid["beta_mod_bridge"].max()))
        plt.plot([mn, mx], [mn, mx], "--")
        plt.xlabel(r"$\beta_{\rm mod}^{ML}$")
        plt.ylabel(r"$\beta_{\rm mod}^{bridge}$")
        plt.title("SPARC bridge test")
        plt.tight_layout()
        plt.savefig(out / "fig_beta_bridge_vs_ml.png", dpi=220)
        plt.close()
        plt.figure(figsize=(7, 5))
        plt.hist(valid["bridge_consistency"], bins=20)
        plt.axvline(0.20, linestyle="--", label="strong threshold")
        plt.axvline(0.35, linestyle="--", label="moderate threshold")
        plt.xlabel("bridge consistency error")
        plt.ylabel("count")
        plt.title("Bridge consistency over galaxies")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out / "fig_bridge_consistency_hist.png", dpi=220)
        plt.close()
    print("\n" + "=" * 90)
    print("SUMMARY")
    print("=" * 90)
    print(json.dumps(summary, indent=2))
    print("\nFiles written:")
    for p in sorted(out.iterdir()):
        if p.is_file():
            print("  -", p)
    print("\nDone.")


if __name__ == "__main__":
    main()

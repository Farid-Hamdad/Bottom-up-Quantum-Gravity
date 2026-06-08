#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS residual test — Paper 9 potential v1

Test:
    S_flux -> Phi_BuP = L_ent^+ S_flux -> C_obs

Target:
    C_obs = log(theta_E_obs / theta_E_baryon)

Purpose:
    Check whether a continuous BuP potential feature explains the residual
    better than the discrete LOW/HIGH regime, and whether it adds information
    beyond logM.

Input CSV must contain at least:
    name, theta_E_obs_arcsec, theta_E_baryon_arcsec, Re_arcsec, logMs

Optional:
    sersic_n_final, sersic_n_measured, sersic_n, q_SIE
"""

import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error


def read_csv_auto(path):
    return pd.read_csv(path, sep=None, engine="python")


def sersic_b_n(n):
    n = float(max(n, 0.3))
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n * n)


def sersic_sigma(R, Re, n):
    R = np.asarray(R, dtype=float)
    Re = max(float(Re), 1e-12)
    n = max(float(n), 0.3)
    b = sersic_b_n(n)
    x = np.maximum(R / Re, 1e-12)
    return np.exp(-b * (x ** (1.0 / n) - 1.0))


def make_polar_nodes(Re, q, n_rings=12, n_theta=16, rmin_factor=0.03, rmax_factor=5.0):
    r = np.logspace(np.log10(rmin_factor * Re), np.log10(rmax_factor * Re), n_rings)
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)

    nodes = []
    for rr in r:
        for tt in th:
            x = rr * np.cos(tt)
            y = rr * q * np.sin(tt)
            Rell = np.sqrt(x * x + (y / max(q, 1e-6)) ** 2)
            nodes.append((x, y, Rell, rr, tt))

    return np.array(nodes, dtype=float)


def build_W(nodes, sigma, xi):
    x = nodes[:, 0]
    y = nodes[:, 1]
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    d2 = dx * dx + dy * dy

    sigma = np.maximum(np.asarray(sigma, dtype=float), 1e-300)

    W = np.sqrt(sigma[:, None] * sigma[None, :]) * np.exp(-d2 / max(2.0 * xi * xi, 1e-12))
    np.fill_diagonal(W, 0.0)

    return W


def laplacian(W, normalized=True):
    W = np.asarray(W, dtype=float)
    d = np.sum(W, axis=1)

    if not normalized:
        return np.diag(d) - W

    invsqrt = 1.0 / np.sqrt(np.maximum(d, 1e-12))
    L = np.eye(W.shape[0]) - (invsqrt[:, None] * W * invsqrt[None, :])
    return L


def pseudo_inverse_potential(L, S, ridge=1e-8):
    S = np.asarray(S, dtype=float)
    S0 = S - np.mean(S)

    # Eigen pseudo-inverse, stable for Laplacian null mode
    vals, vecs = np.linalg.eigh(L)
    inv = np.zeros_like(vals)
    mask = vals > ridge
    inv[mask] = 1.0 / vals[mask]

    Lplus = (vecs * inv[None, :]) @ vecs.T
    phi = Lplus @ S0
    phi = phi - np.mean(phi)

    return phi



def radial_gradient_feature(nodes, phi):
    """
    Safe radial gradient on a polar grid.

    The previous version sorted all nodes by R and called np.gradient(phi, R).
    But many nodes share the same radius for different theta values, producing
    zero radial spacings and NaN/Inf gradients.

    This version:
      1. groups nodes by unique radial shells,
      2. averages phi on each shell,
      3. computes dphi/dR across shells,
      4. maps the shell gradient back to all nodes.
    """
    R = np.asarray(nodes[:, 2], dtype=float)
    phi = np.asarray(phi, dtype=float)

    # Round only for grouping identical polar shells robustly
    R_round = np.round(R, 12)
    unique_R = np.unique(R_round)

    if len(unique_R) < 3:
        return np.zeros_like(phi, dtype=float)

    shell_R = []
    shell_phi = []

    for rr in unique_R:
        mask = R_round == rr
        vals = phi[mask]
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            shell_phi.append(0.0)
        else:
            shell_phi.append(float(np.mean(vals)))
        shell_R.append(float(np.mean(R[mask])))

    shell_R = np.asarray(shell_R, dtype=float)
    shell_phi = np.asarray(shell_phi, dtype=float)

    order = np.argsort(shell_R)
    shell_R = shell_R[order]
    shell_phi = shell_phi[order]

    # Ensure strictly increasing radii
    ok = np.isfinite(shell_R) & np.isfinite(shell_phi)
    shell_R = shell_R[ok]
    shell_phi = shell_phi[ok]

    if len(shell_R) < 3:
        return np.zeros_like(phi, dtype=float)

    grad_shell = np.gradient(shell_phi, shell_R)
    grad_shell = np.nan_to_num(grad_shell, nan=0.0, posinf=0.0, neginf=0.0)

    # Map back to each node
    grad_map = {}
    for rr, gg in zip(np.round(shell_R, 12), grad_shell):
        grad_map[rr] = gg

    gout = np.zeros_like(phi, dtype=float)
    for i, rr in enumerate(R_round):
        gout[i] = grad_map.get(rr, 0.0)

    return gout


def safe_mean(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if len(x) else np.nan


def safe_std(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.std(x)) if len(x) else np.nan


def get_sersic_n(row, default=4.0):
    for c in ["sersic_n_final", "sersic_n_measured", "sersic_n", "n"]:
        if c in row and pd.notna(row[c]):
            return float(row[c])
    return float(default)


def get_q(row):
    for c in ["q_SIE", "q", "axis_ratio"]:
        if c in row and pd.notna(row[c]):
            q = float(row[c])
            if 0.1 <= q <= 1.5:
                return q
    return 1.0


def compute_phi_features(row, args):
    Re = float(row["Re_arcsec"])
    nser = get_sersic_n(row, args.default_sersic_n)
    q = get_q(row)

    nodes = make_polar_nodes(
        Re,
        q,
        n_rings=args.n_rings,
        n_theta=args.n_theta,
        rmin_factor=args.rmin_factor,
        rmax_factor=args.rmax_factor,
    )

    R = nodes[:, 2]
    sigma = sersic_sigma(R, Re, nser)
    sigma_norm = sigma / np.maximum(np.sum(sigma), 1e-12)

    xi = args.xi_factor * Re
    W = build_W(nodes, sigma_norm, xi=xi)

    deg = np.sum(W, axis=1)
    rho_ent = deg / np.maximum(np.sum(deg), 1e-12)

    # Paper 9 proxy sources
    # S1: mismatch energy between entanglement density and baryonic normalized density
    S1 = (rho_ent - sigma_norm) ** 2

    # S2: gradient-supported source
    grad_rho = radial_gradient_feature(nodes, rho_ent)
    S2 = S1 + args.grad_weight * np.abs(grad_rho)

    # S3: flux-like source, minimal proxy
    # flux proxy = weighted radial transport imbalance
    x = nodes[:, 0]
    y = nodes[:, 1]
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    dist = np.sqrt(dx * dx + dy * dy) + 1e-12
    radial_dir = R[:, None] - R[None, :]
    flux_proxy = np.sum(W * np.abs(radial_dir) / dist, axis=1)
    flux_proxy = flux_proxy / np.maximum(np.sum(flux_proxy), 1e-12)

    S3 = S2 + args.flux_weight * flux_proxy

    L = laplacian(W, normalized=True)

    features = {
        "Re_arcsec_used": Re,
        "sersic_n_used": nser,
        "q_used": q,
        "xi_used": xi,
    }

    for label, S in [("S1", S1), ("S2", S2), ("S3", S3)]:
        phi = pseudo_inverse_potential(L, S, ridge=args.ridge)
        grad_phi = radial_gradient_feature(nodes, phi)

        inner = R <= Re
        mid = (R > Re) & (R <= 2.0 * Re)
        outer = R > 2.0 * Re

        features[f"{label}_sum"] = float(np.sum(S))
        features[f"{label}_std"] = safe_std(S)
        features[f"{label}_phi_rms"] = float(np.sqrt(np.mean(phi * phi)))
        features[f"{label}_phi_abs_mean"] = safe_mean(np.abs(phi))
        features[f"{label}_phi_std"] = safe_std(phi)
        features[f"{label}_phi_inner_mean"] = safe_mean(phi[inner])
        features[f"{label}_phi_mid_mean"] = safe_mean(phi[mid])
        features[f"{label}_phi_outer_mean"] = safe_mean(phi[outer])
        features[f"{label}_phi_inner_minus_outer"] = features[f"{label}_phi_inner_mean"] - features[f"{label}_phi_outer_mean"]
        features[f"{label}_grad_rms"] = float(np.sqrt(np.mean(grad_phi * grad_phi)))
        features[f"{label}_grad_abs_mean"] = safe_mean(np.abs(grad_phi))
        features[f"{label}_force_inner_mean"] = safe_mean(np.abs(grad_phi[inner]))
        features[f"{label}_force_outer_mean"] = safe_mean(np.abs(grad_phi[outer]))
        features[f"{label}_force_inner_minus_outer"] = features[f"{label}_force_inner_mean"] - features[f"{label}_force_outer_mean"]

    return features


def metrics_y(y, pred):
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)

    if len(y) > 2 and np.std(pred) > 1e-12:
        pr = pearsonr(y, pred)[0]
        sr = spearmanr(y, pred).correlation
    else:
        pr = np.nan
        sr = np.nan

    return {
        "pearson": float(pr) if np.isfinite(pr) else np.nan,
        "spearman": float(sr) if np.isfinite(sr) else np.nan,
        "rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "mae": float(mean_absolute_error(y, pred)),
        "bias": float(np.mean(pred - y)),
    }


def theta_improvement(df, C_pred):
    theta_obs = df["theta_E_obs_arcsec"].values.astype(float)
    theta_bar = df["theta_E_baryon_arcsec"].values.astype(float)
    theta_pred = theta_bar * np.exp(C_pred)

    log_obs = np.log(theta_obs + 1e-12)
    log_bar = np.log(theta_bar + 1e-12)
    log_pred = np.log(theta_pred + 1e-12)

    rmse_bar = float(np.sqrt(mean_squared_error(log_obs, log_bar)))
    rmse_pred = float(np.sqrt(mean_squared_error(log_obs, log_pred)))

    return {
        "rmse_log_model": rmse_pred,
        "rmse_log_baryon": rmse_bar,
        "improvement_pct": 100.0 * (rmse_bar - rmse_pred) / rmse_bar,
    }



def standardize(train, test, features):
    Xtr = [np.ones(len(train))]
    Xte = [np.ones(len(test))]
    stats = {}

    for f in features:
        xtr = train[f].values.astype(float)
        xte = test[f].values.astype(float)

        # Keep only finite values for training statistics
        finite = np.isfinite(xtr)

        if np.sum(finite) < 3:
            mu = 0.0
            sd = 1.0
        else:
            mu = float(np.mean(xtr[finite]))
            sd = float(np.std(xtr[finite]))

        if not np.isfinite(sd) or sd < 1e-12:
            sd = 1.0

        # Impute non-finite values by train mean
        xtr = np.where(np.isfinite(xtr), xtr, mu)
        xte = np.where(np.isfinite(xte), xte, mu)

        stats[f] = {"mean": mu, "std": sd}

        Xtr.append((xtr - mu) / sd)
        Xte.append((xte - mu) / sd)

    Xtr = np.column_stack(Xtr)
    Xte = np.column_stack(Xte)

    Xtr = np.nan_to_num(Xtr, nan=0.0, posinf=0.0, neginf=0.0)
    Xte = np.nan_to_num(Xte, nan=0.0, posinf=0.0, neginf=0.0)

    return Xtr, Xte, stats


def fit_predict(train, test, features):
    ytr = train["C_obs"].values.astype(float)
    Xtr, Xte, stats = standardize(train, test, features)

    beta = np.linalg.lstsq(Xtr, ytr, rcond=None)[0]
    pred = Xte @ beta

    coeffs = {"intercept": float(beta[0])}
    for f, b in zip(features, beta[1:]):
        coeffs[f] = float(b)

    return pred, coeffs, stats


def full_model(df, features):
    pred, coeffs, stats = fit_predict(df, df, features)
    return pred, coeffs, stats


def loo_model(df, features):
    preds = np.zeros(len(df), dtype=float)

    for i in range(len(df)):
        train = df.drop(df.index[i])
        test = df.iloc[[i]]
        pred, _, _ = fit_predict(train, test, features)
        preds[i] = pred[0]

    return preds


def evaluate_model(df, model_name, features, outdir):
    model_dir = os.path.join(outdir, model_name)
    os.makedirs(model_dir, exist_ok=True)

    C_obs = df["C_obs"].values.astype(float)

    C_full, coeffs, stats = full_model(df, features)
    C_loo = loo_model(df, features)

    m_full = metrics_y(C_obs, C_full)
    m_loo = metrics_y(C_obs, C_loo)
    th_full = theta_improvement(df, C_full)
    th_loo = theta_improvement(df, C_loo)

    pred = df.copy()
    pred["C_pred_full"] = C_full
    pred["C_pred_LOO"] = C_loo
    pred["theta_pred_full"] = df["theta_E_baryon_arcsec"].values * np.exp(C_full)
    pred["theta_pred_LOO"] = df["theta_E_baryon_arcsec"].values * np.exp(C_loo)
    pred.to_csv(os.path.join(model_dir, "predictions.csv"), index=False)

    summary = {
        "model": model_name,
        "features": features,
        "n_features": len(features),
        "coefficients": coeffs,
        "feature_stats": stats,
        "metrics_C_full": m_full,
        "metrics_C_LOO": m_loo,
        "theta_full": th_full,
        "theta_LOO": th_loo,
        "files": {"predictions": "predictions.csv"},
    }

    with open(os.path.join(model_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def bootstrap(df, features, n_boot, seed):
    rng = np.random.default_rng(seed)
    n = len(df)
    vals = []

    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = df.iloc[idx].copy()
        C_full, _, _ = full_model(sample, features)
        th = theta_improvement(sample, C_full)
        vals.append(th["improvement_pct"])

    vals = np.asarray(vals, dtype=float)
    return {
        "n_boot": int(n_boot),
        "prob_improvement": float(np.mean(vals > 0)),
        "improvement_mean": float(np.mean(vals)),
        "improvement_q05": float(np.quantile(vals, 0.05)),
        "improvement_q50": float(np.quantile(vals, 0.50)),
        "improvement_q95": float(np.quantile(vals, 0.95)),
    }


def added_value_permutation(df, base_features, add_features, n_perm, seed):
    rng = np.random.default_rng(seed)

    base = evaluate_model(df, "__tmp_base", base_features, "/tmp")
    full = evaluate_model(df, "__tmp_full", base_features + add_features, "/tmp")

    obs_delta = full["theta_LOO"]["improvement_pct"] - base["theta_LOO"]["improvement_pct"]

    count = 0
    deltas = []

    for _ in range(n_perm):
        sh = df.copy()
        for f in add_features:
            sh[f] = rng.permutation(sh[f].values)

        base_sh = evaluate_model(sh, "__tmp_base_sh", base_features, "/tmp")
        full_sh = evaluate_model(sh, "__tmp_full_sh", base_features + add_features, "/tmp")

        delta = full_sh["theta_LOO"]["improvement_pct"] - base_sh["theta_LOO"]["improvement_pct"]
        deltas.append(delta)

        if delta >= obs_delta:
            count += 1

    return {
        "obs_added_LOO_improvement_pct": float(obs_delta),
        "perm_p_one_sided": float((1.0 + count) / (1.0 + n_perm)),
        "perm_delta_mean": float(np.mean(deltas)),
        "perm_delta_q05": float(np.quantile(deltas, 0.05)),
        "perm_delta_q50": float(np.quantile(deltas, 0.50)),
        "perm_delta_q95": float(np.quantile(deltas, 0.95)),
        "n_perm": int(n_perm),
    }


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--output-dir", required=True)

    ap.add_argument("--n-rings", type=int, default=12)
    ap.add_argument("--n-theta", type=int, default=16)
    ap.add_argument("--rmin-factor", type=float, default=0.03)
    ap.add_argument("--rmax-factor", type=float, default=5.0)
    ap.add_argument("--xi-factor", type=float, default=0.35)
    ap.add_argument("--default-sersic-n", type=float, default=4.0)

    ap.add_argument("--grad-weight", type=float, default=0.5)
    ap.add_argument("--flux-weight", type=float, default=0.5)
    ap.add_argument("--ridge", type=float, default=1e-8)

    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = read_csv_auto(args.input_csv)

    required = ["name", "theta_E_obs_arcsec", "theta_E_baryon_arcsec", "Re_arcsec", "logMs"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    df = df.dropna(subset=required).copy()
    df = df[(df["theta_E_obs_arcsec"] > 0) & (df["theta_E_baryon_arcsec"] > 0)].copy()

    df["C_obs"] = np.log(df["theta_E_obs_arcsec"].astype(float) / df["theta_E_baryon_arcsec"].astype(float))

    print("Computing Paper 9 Phi_BuP features...")
    rows = []
    for i, row in df.iterrows():
        name = row["name"]
        print(f"[{len(rows)+1}/{len(df)}] {name}")
        feats = compute_phi_features(row, args)
        feats["name"] = name
        rows.append(feats)

    feats_df = pd.DataFrame(rows)

    full = df.merge(feats_df, on="name", how="left")

    features_path = os.path.join(args.output_dir, "paper9_phi_features.csv")
    full.to_csv(features_path, index=False)

    phi_features = []
    for c in full.columns:
        if not (c.startswith("S1_") or c.startswith("S2_") or c.startswith("S3_")):
            continue
        if not pd.api.types.is_numeric_dtype(full[c]):
            continue

        arr = full[c].values.astype(float)
        finite = np.isfinite(arr)

        # Keep only informative finite features
        if np.sum(finite) < max(10, int(0.5 * len(arr))):
            continue

        std = np.std(arr[finite])
        if np.isfinite(std) and std > 1e-12:
            phi_features.append(c)

    # Candidate subsets
    s1_features = [c for c in phi_features if c.startswith("S1_")]
    s2_features = [c for c in phi_features if c.startswith("S2_")]
    s3_features = [c for c in phi_features if c.startswith("S3_")]

    # Single-feature screening
    screening = []
    for f in phi_features:
        s = evaluate_model(full, f"single_{f}", [f], args.output_dir)
        screening.append({
            "feature": f,
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
        })

    screen_df = pd.DataFrame(screening).sort_values("loo_improvement_pct", ascending=False)
    screen_df.to_csv(os.path.join(args.output_dir, "single_phi_feature_screen.csv"), index=False)

    best_single = screen_df.iloc[0]["feature"] if len(screen_df) else None

    models = {
        "baseline_constant": [],
        "logM_only": ["logMs"],
    }

    if best_single:
        models["Phi_best_single"] = [best_single]
        models["logM_plus_Phi_best"] = ["logMs", best_single]

    # Keep compact multi-feature models to avoid overfit
    compact = []
    for name in [
        "S1_phi_rms", "S1_phi_inner_minus_outer", "S1_grad_rms",
        "S2_phi_rms", "S2_phi_inner_minus_outer", "S2_grad_rms",
        "S3_phi_rms", "S3_phi_inner_minus_outer", "S3_grad_rms",
    ]:
        if name in full.columns and np.nanstd(full[name].values.astype(float)) > 1e-12:
            compact.append(name)

    if compact:
        models["Phi_compact"] = compact
        models["logM_plus_Phi_compact"] = ["logMs"] + compact

    summaries = []
    for name, feats in models.items():
        print("Running model:", name, feats)
        summaries.append(evaluate_model(full, name, feats, args.output_dir))

    rows = []
    for s in summaries:
        rows.append({
            "model": s["model"],
            "features": ",".join(s["features"]),
            "n_features": s["n_features"],
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
            "full_rmse_log": s["theta_full"]["rmse_log_model"],
            "full_pearson_C": s["metrics_C_full"]["pearson"],
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "loo_spearman_C": s["metrics_C_LOO"]["spearman"],
        })

    comp = pd.DataFrame(rows).sort_values("loo_improvement_pct", ascending=False)
    comp.to_csv(os.path.join(args.output_dir, "paper9_phi_model_summary.csv"), index=False)

    boot = {}
    for name, feats in models.items():
        if name in ["logM_only", "Phi_best_single", "logM_plus_Phi_best"]:
            boot[name] = bootstrap(full, feats, args.n_boot, args.seed)

    perm_added = {}
    if best_single:
        perm_added = added_value_permutation(
            full,
            base_features=["logMs"],
            add_features=[best_single],
            n_perm=args.n_perm,
            seed=args.seed,
        )

    summary = {
        "title": "BuP SLACS Paper 9 potential residual test v1",
        "input_csv": args.input_csv,
        "n": int(len(full)),
        "target": "C_obs = log(theta_E_obs/theta_E_baryon)",
        "settings": vars(args),
        "n_phi_features": int(len(phi_features)),
        "best_single_phi_feature": best_single,
        "best_by_LOO": comp.iloc[0].to_dict() if len(comp) else {},
        "bootstrap": boot,
        "added_value_permutation_logM_plus_bestPhi": perm_added,
        "files": {
            "features": "paper9_phi_features.csv",
            "single_feature_screen": "single_phi_feature_screen.csv",
            "model_summary": "paper9_phi_model_summary.csv",
        },
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 100)
    print("Paper 9 Phi_BuP model summary")
    print("=" * 100)
    print(comp.to_string(index=False))

    print("\nBest single Phi feature:")
    print(best_single)

    print("\nWrote:", os.path.join(args.output_dir, "summary.json"))


if __name__ == "__main__":
    main()

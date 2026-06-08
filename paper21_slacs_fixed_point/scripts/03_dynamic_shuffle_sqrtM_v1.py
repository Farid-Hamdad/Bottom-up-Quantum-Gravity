#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS Paper 9 — dynamic sqrtM shuffle control v1

True dynamic control:
    sqrt(Mstar)_true or sqrt(Mstar)_shuffle is injected BEFORE W_ij construction.

Pipeline:
    amplitude -> W_ij -> L_raw = D-W -> S_zscore -> Phi = L_raw^+ S -> feature

Canonical dynamic feature:
    rawL_zS_amp_sqrtMstar_phi_mid_mean

Target:
    C_obs = log(theta_E_obs / theta_E_baryon)
"""

import argparse
import json
import os
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error


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


def get_mstar_rel(row):
    if "Mstar_rel" in row and pd.notna(row["Mstar_rel"]):
        return float(row["Mstar_rel"])
    if "logMs" in row and pd.notna(row["logMs"]):
        return float(10 ** (float(row["logMs"]) - 11.0))
    return 1.0


def make_polar_nodes(Re, q, n_rings, n_theta, rmin_factor, rmax_factor):
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


def build_W(nodes, sigma_amp, xi):
    x = nodes[:, 0]
    y = nodes[:, 1]

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    d2 = dx * dx + dy * dy

    sigma_amp = np.maximum(np.asarray(sigma_amp, dtype=float), 1e-300)

    W = np.sqrt(sigma_amp[:, None] * sigma_amp[None, :]) * np.exp(
        -d2 / max(2.0 * xi * xi, 1e-12)
    )

    np.fill_diagonal(W, 0.0)
    return W


def laplacian_raw(W):
    d = np.sum(W, axis=1)
    return np.diag(d) - W


def pseudo_inverse_potential(L, S, ridge=1e-8):
    S = np.asarray(S, dtype=float)
    S0 = S - np.mean(S)

    vals, vecs = np.linalg.eigh(L)

    inv = np.zeros_like(vals)
    mask = vals > ridge
    inv[mask] = 1.0 / vals[mask]

    Lplus = (vecs * inv[None, :]) @ vecs.T

    phi = Lplus @ S0
    phi -= np.mean(phi)

    return np.nan_to_num(phi, nan=0.0, posinf=0.0, neginf=0.0)


def radial_gradient(nodes, field):
    R = np.asarray(nodes[:, 2], dtype=float)
    field = np.asarray(field, dtype=float)

    R_round = np.round(R, 12)
    shells = np.unique(R_round)

    if len(shells) < 3:
        return np.zeros_like(field)

    shell_R = []
    shell_field = []

    for rr in shells:
        mask = R_round == rr
        shell_R.append(float(np.mean(R[mask])))
        shell_field.append(float(np.mean(field[mask])))

    shell_R = np.asarray(shell_R)
    shell_field = np.asarray(shell_field)

    order = np.argsort(shell_R)
    shell_R = shell_R[order]
    shell_field = shell_field[order]

    grad_shell = np.gradient(shell_field, shell_R)
    grad_shell = np.nan_to_num(grad_shell, nan=0.0, posinf=0.0, neginf=0.0)

    grad_map = {round(float(rr), 12): float(gg) for rr, gg in zip(shell_R, grad_shell)}

    out = np.zeros_like(field)
    for i, rr in enumerate(R_round):
        out[i] = grad_map.get(round(float(rr), 12), 0.0)

    return out


def safe_mean(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if len(x) else np.nan


def zscale(v):
    v = np.asarray(v, dtype=float)
    v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    s = np.std(v)
    if not np.isfinite(s) or s < 1e-12:
        return np.zeros_like(v)
    return v / s


def compute_phi_feature(row, sqrtM_amp, args, prefix):
    Re = float(row["Re_arcsec"])
    nser = get_sersic_n(row, args.default_sersic_n)
    q = get_q(row)

    nodes = make_polar_nodes(
        Re,
        q,
        args.n_rings,
        args.n_theta,
        args.rmin_factor,
        args.rmax_factor,
    )

    R = nodes[:, 2]

    sigma_shape = sersic_sigma(R, Re, nser)
    sigma_shape = sigma_shape / np.maximum(np.sum(sigma_shape), 1e-12)

    sigma_amp = sigma_shape * max(float(sqrtM_amp), 1e-12)

    xi = args.xi_factor * Re

    W = build_W(nodes, sigma_amp, xi)
    L = laplacian_raw(W)

    deg = np.sum(W, axis=1)

    # rawL_zS mode from v2:
    # delta normalized, source z-scored, amplitude enters through raw L.
    rho_field = deg / np.maximum(np.sum(deg), 1e-12)
    sigma_field = sigma_shape

    delta = rho_field - sigma_field

    T00 = delta ** 2
    grad_delta = radial_gradient(nodes, delta)
    Taa = delta * grad_delta

    S = zscale(T00) - 0.5 * zscale(Taa)

    phi = pseudo_inverse_potential(L, S, ridge=args.ridge)
    force = np.abs(radial_gradient(nodes, phi))

    inner = R <= Re
    mid = (R > Re) & (R <= 2.0 * Re)
    outer = R > 2.0 * Re

    return {
        f"{prefix}_phi_mid_mean": safe_mean(phi[mid]),
        f"{prefix}_phi_inner_minus_outer": safe_mean(phi[inner]) - safe_mean(phi[outer]),
        f"{prefix}_force_mid_mean": safe_mean(force[mid]),
        f"{prefix}_force_inner_minus_outer": safe_mean(force[inner]) - safe_mean(force[outer]),
        f"{prefix}_grad_rms": float(np.sqrt(np.mean(force * force))),
        f"{prefix}_amp_used": float(sqrtM_amp),
    }


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

    for f in features:
        xtr = train[f].values.astype(float)
        xte = test[f].values.astype(float)

        finite = np.isfinite(xtr)
        if np.sum(finite) < 3:
            mu = 0.0
            sd = 1.0
        else:
            mu = float(np.mean(xtr[finite]))
            sd = float(np.std(xtr[finite]))

        if not np.isfinite(sd) or sd < 1e-12:
            sd = 1.0

        xtr = np.where(np.isfinite(xtr), xtr, mu)
        xte = np.where(np.isfinite(xte), xte, mu)

        Xtr.append((xtr - mu) / sd)
        Xte.append((xte - mu) / sd)

    return np.column_stack(Xtr), np.column_stack(Xte)


def fit_predict(train, test, features):
    ytr = train["C_obs"].values.astype(float)
    Xtr, Xte = standardize(train, test, features)

    beta = np.linalg.lstsq(Xtr, ytr, rcond=None)[0]
    return Xte @ beta


def evaluate_model(df, features):
    y = df["C_obs"].values.astype(float)

    C_full = fit_predict(df, df, features)

    C_loo = np.zeros(len(df))
    for i in range(len(df)):
        train = df.drop(df.index[i])
        test = df.iloc[[i]]
        C_loo[i] = fit_predict(train, test, features)[0]

    return {
        "features": features,
        "n_features": len(features),
        "metrics_C_full": metrics_y(y, C_full),
        "metrics_C_LOO": metrics_y(y, C_loo),
        "theta_full": theta_improvement(df, C_full),
        "theta_LOO": theta_improvement(df, C_loo),
    }


def bootstrap_model(df, features, n_boot, seed):
    rng = np.random.default_rng(seed)
    n = len(df)
    vals = []

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = df.iloc[idx].copy()
        C_full = fit_predict(sample, sample, features)
        vals.append(theta_improvement(sample, C_full)["improvement_pct"])

    vals = np.asarray(vals, dtype=float)

    return {
        "n_boot": int(n_boot),
        "prob_improvement": float(np.mean(vals > 0)),
        "improvement_mean": float(np.mean(vals)),
        "improvement_q05": float(np.quantile(vals, 0.05)),
        "improvement_q50": float(np.quantile(vals, 0.50)),
        "improvement_q95": float(np.quantile(vals, 0.95)),
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
    ap.add_argument("--ridge", type=float, default=1e-8)

    ap.add_argument("--n-boot", type=int, default=1000)
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
    df["sqrtM_proxy"] = np.sqrt(np.maximum(df.apply(get_mstar_rel, axis=1).values.astype(float), 1e-12))

    rng = np.random.default_rng(args.seed)
    shuffled_amp = rng.permutation(df["sqrtM_proxy"].values)

    rows = []

    print("Computing true/shuffled dynamic sqrtM control...")
    for k, ((idx, row), amp_shuf) in enumerate(zip(df.iterrows(), shuffled_amp), start=1):
        print(f"[{k}/{len(df)}] {row['name']}")

        amp_true = float(row["sqrtM_proxy"])

        out = {
            "name": row["name"],
            "logMs": float(row["logMs"]),
            "sqrtM_proxy": amp_true,
            "sqrtM_shuffle": float(amp_shuf),
        }

        out.update(compute_phi_feature(row, amp_true, args, "Phi_true"))
        out.update(compute_phi_feature(row, amp_shuf, args, "Phi_shuffle"))

        rows.append(out)

    feats = pd.DataFrame(rows)
    full = df.merge(feats, on=["name", "logMs", "sqrtM_proxy"], how="left")

    full.to_csv(os.path.join(args.output_dir, "dynamic_shuffle_sqrtM_features.csv"), index=False)

    candidate_features = [
        "Phi_true_phi_mid_mean",
        "Phi_true_phi_inner_minus_outer",
        "Phi_true_force_mid_mean",
        "Phi_true_force_inner_minus_outer",
        "Phi_true_grad_rms",
        "Phi_shuffle_phi_mid_mean",
        "Phi_shuffle_phi_inner_minus_outer",
        "Phi_shuffle_force_mid_mean",
        "Phi_shuffle_force_inner_minus_outer",
        "Phi_shuffle_grad_rms",
        "sqrtM_proxy",
        "sqrtM_shuffle",
        "logMs",
    ]

    results = []

    for f in candidate_features:
        if f not in full.columns:
            continue
        arr = full[f].values.astype(float)
        finite = np.isfinite(arr)
        if np.sum(finite) <= 10 or np.std(arr[finite]) <= 1e-12:
            continue

        s = evaluate_model(full, [f])

        if f.startswith("Phi_true"):
            mode = "Phi_true"
        elif f.startswith("Phi_shuffle"):
            mode = "Phi_shuffle"
        elif f == "sqrtM_proxy":
            mode = "sqrtM_true"
        elif f == "sqrtM_shuffle":
            mode = "sqrtM_shuffle"
        elif f == "logMs":
            mode = "logM"
        else:
            mode = "other"

        results.append({
            "model": f,
            "mode": mode,
            "features": f,
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "loo_spearman_C": s["metrics_C_LOO"]["spearman"],
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
            "full_pearson_C": s["metrics_C_full"]["pearson"],
        })

    # Combined models
    combo_models = {
        "Phi_true_plus_sqrtM": ["sqrtM_proxy", "Phi_true_phi_mid_mean"],
        "Phi_shuffle_plus_sqrtM": ["sqrtM_proxy", "Phi_shuffle_phi_mid_mean"],
        "Phi_true_plus_shuffleM": ["sqrtM_shuffle", "Phi_true_phi_mid_mean"],
        "Phi_shuffle_plus_shuffleM": ["sqrtM_shuffle", "Phi_shuffle_phi_mid_mean"],
        "logM_plus_Phi_true": ["logMs", "Phi_true_phi_mid_mean"],
    }

    for name, feats_model in combo_models.items():
        s = evaluate_model(full, feats_model)
        results.append({
            "model": name,
            "mode": "combined",
            "features": ",".join(feats_model),
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "loo_spearman_C": s["metrics_C_LOO"]["spearman"],
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
            "full_pearson_C": s["metrics_C_full"]["pearson"],
        })

    res = pd.DataFrame(results).sort_values("loo_improvement_pct", ascending=False)
    res.to_csv(os.path.join(args.output_dir, "dynamic_shuffle_sqrtM_summary.csv"), index=False)

    by_mode = (
        res.groupby("mode")["loo_improvement_pct"]
        .agg(["max", "median", "mean"])
        .reset_index()
        .sort_values("max", ascending=False)
    )
    by_mode.to_csv(os.path.join(args.output_dir, "dynamic_shuffle_sqrtM_by_mode.csv"), index=False)

    boot = {}
    for name, feats_model in {
        "logM_only": ["logMs"],
        "sqrtM_only": ["sqrtM_proxy"],
        "Phi_true_phi_mid_mean": ["Phi_true_phi_mid_mean"],
        "Phi_shuffle_phi_mid_mean": ["Phi_shuffle_phi_mid_mean"],
        "Phi_true_plus_sqrtM": ["sqrtM_proxy", "Phi_true_phi_mid_mean"],
        "Phi_shuffle_plus_sqrtM": ["sqrtM_proxy", "Phi_shuffle_phi_mid_mean"],
    }.items():
        boot[name] = bootstrap_model(full, feats_model, args.n_boot, args.seed)

    corr = {}
    for f in ["Phi_true_phi_mid_mean", "Phi_shuffle_phi_mid_mean"]:
        if f in full.columns:
            phi = full[f].values.astype(float)
            sm = full["sqrtM_proxy"].values.astype(float)
            corr[f + "_vs_sqrtM"] = {
                "pearson": float(pearsonr(phi, sm)[0]),
                "spearman": float(spearmanr(phi, sm).correlation),
            }

    summary = {
        "title": "BuP SLACS Paper 9 dynamic sqrtM shuffle control v1",
        "input_csv": args.input_csv,
        "n": int(len(full)),
        "target": "C_obs = log(theta_E_obs/theta_E_baryon)",
        "settings": vars(args),
        "best_by_LOO": res.iloc[0].to_dict() if len(res) else {},
        "by_mode": by_mode.to_dict(orient="records"),
        "correlations": corr,
        "bootstrap": boot,
        "files": {
            "features": "dynamic_shuffle_sqrtM_features.csv",
            "summary": "dynamic_shuffle_sqrtM_summary.csv",
            "by_mode": "dynamic_shuffle_sqrtM_by_mode.csv",
        },
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 100)
    print("DYNAMIC SQRTM SHUFFLE CONTROL SUMMARY")
    print("=" * 100)
    print(res.to_string(index=False))

    print("\nBY MODE")
    print(by_mode.to_string(index=False))

    print("\nCORRELATIONS")
    print(json.dumps(corr, indent=2))

    print("\nWrote:", os.path.join(args.output_dir, "summary.json"))


if __name__ == "__main__":
    main()

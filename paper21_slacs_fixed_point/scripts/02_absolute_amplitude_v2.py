#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP SLACS Paper 9 — absolute amplitude test v2

Purpose:
    Test whether absolute stellar amplitude truly enters the BuP dynamics
    once we remove the two amplitude-killing operations:

    1) normalized Laplacian L_norm
    2) z-scored source S

We compare:

    norm_zS     : normalized Laplacian + z-scored source
    rawL_zS     : raw Laplacian D-W + z-scored source
    norm_rawS   : normalized Laplacian + raw source
    rawL_rawS   : raw Laplacian D-W + raw source

Amplitude modes:

    shape_norm
    amp_sqrtMstar
    amp_Mstar
    amp_logMstar

Target:
    C_obs = log(theta_E_obs/theta_E_baryon)
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


def laplacian(W, mode):
    W = np.asarray(W, dtype=float)
    d = np.sum(W, axis=1)

    if mode == "norm":
        invsqrt = 1.0 / np.sqrt(np.maximum(d, 1e-12))
        return np.eye(W.shape[0]) - (invsqrt[:, None] * W * invsqrt[None, :])

    if mode == "raw":
        return np.diag(d) - W

    raise ValueError(f"Unknown Laplacian mode: {mode}")


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


def compute_source(delta, use_zscore):
    T00 = delta ** 2
    grad_delta = radial_gradient(current_nodes, delta)
    Taa = delta * grad_delta

    if use_zscore:
        return zscale(T00) - 0.5 * zscale(Taa)

    return T00 - 0.5 * Taa


# Global holder used only to avoid passing nodes into compute_source repeatedly
current_nodes = None


def compute_features_for_mode(row, args, amp_mode, variant):
    global current_nodes

    Re = float(row["Re_arcsec"])
    nser = get_sersic_n(row, args.default_sersic_n)
    q = get_q(row)
    logMs = float(row["logMs"])
    Mstar_rel = get_mstar_rel(row)

    nodes = make_polar_nodes(
        Re,
        q,
        args.n_rings,
        args.n_theta,
        args.rmin_factor,
        args.rmax_factor,
    )

    current_nodes = nodes

    R = nodes[:, 2]

    sigma_shape = sersic_sigma(R, Re, nser)
    sigma_shape = sigma_shape / np.maximum(np.sum(sigma_shape), 1e-12)

    if amp_mode == "shape_norm":
        amp_factor = 1.0
    elif amp_mode == "amp_sqrtMstar":
        amp_factor = np.sqrt(max(Mstar_rel, 1e-12))
    elif amp_mode == "amp_Mstar":
        amp_factor = max(Mstar_rel, 1e-12)
    elif amp_mode == "amp_logMstar":
        amp_factor = max(logMs - 10.0, 1e-3)
    else:
        raise ValueError(f"Unknown amp mode: {amp_mode}")

    sigma_amp = sigma_shape * amp_factor

    xi = args.xi_factor * Re

    W = build_W(nodes, sigma_amp, xi)

    if variant.startswith("norm"):
        L = laplacian(W, "norm")
    elif variant.startswith("rawL"):
        L = laplacian(W, "raw")
    else:
        raise ValueError(f"Unknown variant: {variant}")

    use_zscore = variant.endswith("zS")

    deg = np.sum(W, axis=1)

    # Two options:
    # - normalized delta for shape-only behavior
    # - amplitude-carrying delta for raw-amplitude behavior
    if "rawS" in variant:
        rho_field = deg
        sigma_field = sigma_amp
    else:
        rho_field = deg / np.maximum(np.sum(deg), 1e-12)
        sigma_field = sigma_shape

    delta = rho_field - sigma_field

    S = compute_source(delta, use_zscore=use_zscore)

    phi = pseudo_inverse_potential(L, S, ridge=args.ridge)
    force = np.abs(radial_gradient(nodes, phi))

    inner = R <= Re
    mid = (R > Re) & (R <= 2.0 * Re)
    outer = R > 2.0 * Re

    prefix = f"{variant}_{amp_mode}"

    return {
        f"{prefix}_force_inner_mean": safe_mean(force[inner]),
        f"{prefix}_force_mid_mean": safe_mean(force[mid]),
        f"{prefix}_force_outer_mean": safe_mean(force[outer]),
        f"{prefix}_force_inner_minus_outer": safe_mean(force[inner]) - safe_mean(force[outer]),
        f"{prefix}_force_mid_minus_outer": safe_mean(force[mid]) - safe_mean(force[outer]),
        f"{prefix}_grad_rms": float(np.sqrt(np.mean(force * force))),
        f"{prefix}_phi_mid_mean": safe_mean(phi[mid]),
        f"{prefix}_phi_inner_minus_outer": safe_mean(phi[inner]) - safe_mean(phi[outer]),
        f"{prefix}_S_std": float(np.std(S)),
        f"{prefix}_amp_factor": float(amp_factor),
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


def parse_feature_meta(feature):
    parts = feature.split("_")

    # variant examples:
    # norm_zS_amp_Mstar_force_mid_mean
    # rawL_rawS_amp_sqrtMstar_force_mid_mean
    variant = "_".join(parts[0:2])

    rest = "_".join(parts[2:])

    amp_mode = None
    for m in ["shape_norm", "amp_sqrtMstar", "amp_Mstar", "amp_logMstar"]:
        if rest.startswith(m):
            amp_mode = m
            break

    return variant, amp_mode


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

    variants = ["norm_zS", "rawL_zS", "norm_rawS", "rawL_rawS"]
    amp_modes = ["shape_norm", "amp_sqrtMstar", "amp_Mstar", "amp_logMstar"]

    rows = []

    print("Computing absolute-amplitude v2 features...")
    for k, (_, row) in enumerate(df.iterrows(), start=1):
        print(f"[{k}/{len(df)}] {row['name']}")

        out = {
            "name": row["name"],
            "logMs": float(row["logMs"]),
            "Mstar_rel_used": get_mstar_rel(row),
        }

        for variant in variants:
            for amp_mode in amp_modes:
                out.update(compute_features_for_mode(row, args, amp_mode, variant))

        rows.append(out)

    feats = pd.DataFrame(rows)
    full = df.merge(feats, on=["name", "logMs"], how="left")

    full.to_csv(os.path.join(args.output_dir, "absolute_amplitude_v2_features.csv"), index=False)

    candidate_features = []

    for variant in variants:
        for amp_mode in amp_modes:
            prefix = f"{variant}_{amp_mode}"
            for suffix in [
                "force_inner_mean",
                "force_mid_mean",
                "force_outer_mean",
                "force_inner_minus_outer",
                "force_mid_minus_outer",
                "grad_rms",
                "phi_mid_mean",
                "phi_inner_minus_outer",
                "S_std",
                "amp_factor",
            ]:
                f = f"{prefix}_{suffix}"
                if f in full.columns:
                    arr = full[f].values.astype(float)
                    finite = np.isfinite(arr)
                    if np.sum(finite) > 10 and np.std(arr[finite]) > 1e-12:
                        candidate_features.append(f)

    rows_screen = []

    for f in candidate_features:
        print("Evaluating", f)
        s = evaluate_model(full, [f])

        variant, amp_mode = parse_feature_meta(f)

        rows_screen.append({
            "feature": f,
            "variant": variant,
            "amp_mode": amp_mode,
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "loo_spearman_C": s["metrics_C_LOO"]["spearman"],
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
            "full_pearson_C": s["metrics_C_full"]["pearson"],
        })

    screen = pd.DataFrame(rows_screen).sort_values("loo_improvement_pct", ascending=False)
    screen.to_csv(os.path.join(args.output_dir, "absolute_amplitude_v2_screen.csv"), index=False)

    by_variant_amp = (
        screen.groupby(["variant", "amp_mode"])["loo_improvement_pct"]
        .agg(["max", "median", "mean"])
        .reset_index()
        .sort_values("max", ascending=False)
    )

    by_variant_amp.to_csv(os.path.join(args.output_dir, "absolute_amplitude_v2_by_variant_amp.csv"), index=False)

    by_variant = (
        screen.groupby("variant")["loo_improvement_pct"]
        .agg(["max", "median", "mean"])
        .reset_index()
        .sort_values("max", ascending=False)
    )

    by_variant.to_csv(os.path.join(args.output_dir, "absolute_amplitude_v2_by_variant.csv"), index=False)

    models = {
        "logM_only": ["logMs"],
    }

    best_feature = screen.iloc[0]["feature"] if len(screen) else None

    if best_feature:
        models["best_amp_v2_single"] = [best_feature]
        models["logM_plus_best_amp_v2"] = ["logMs", best_feature]

    canonical_candidates = [
        "norm_zS_shape_norm_force_inner_minus_outer",
        "rawL_zS_amp_Mstar_force_inner_minus_outer",
        "norm_rawS_amp_Mstar_force_inner_minus_outer",
        "rawL_rawS_amp_Mstar_force_inner_minus_outer",
        "rawL_rawS_amp_sqrtMstar_force_inner_minus_outer",
        "rawL_rawS_amp_Mstar_force_mid_mean",
        "rawL_rawS_amp_sqrtMstar_force_mid_mean",
        "rawL_rawS_amp_Mstar_amp_factor",
    ]

    for f in canonical_candidates:
        if f in full.columns:
            models[f"canonical_{f}"] = [f]

    comp_rows = []

    for name, feats_model in models.items():
        s = evaluate_model(full, feats_model)
        comp_rows.append({
            "model": name,
            "features": ",".join(feats_model),
            "n_features": len(feats_model),
            "full_improvement_pct": s["theta_full"]["improvement_pct"],
            "full_rmse_log": s["theta_full"]["rmse_log_model"],
            "full_pearson_C": s["metrics_C_full"]["pearson"],
            "loo_improvement_pct": s["theta_LOO"]["improvement_pct"],
            "loo_rmse_log": s["theta_LOO"]["rmse_log_model"],
            "loo_pearson_C": s["metrics_C_LOO"]["pearson"],
            "loo_spearman_C": s["metrics_C_LOO"]["spearman"],
        })

    comp = pd.DataFrame(comp_rows).sort_values("loo_improvement_pct", ascending=False)
    comp.to_csv(os.path.join(args.output_dir, "absolute_amplitude_v2_model_summary.csv"), index=False)

    boot = {}
    for name, feats_model in models.items():
        if name in ["logM_only", "best_amp_v2_single", "logM_plus_best_amp_v2"]:
            boot[name] = bootstrap_model(full, feats_model, args.n_boot, args.seed)

    summary = {
        "title": "BuP SLACS Paper 9 absolute amplitude test v2",
        "input_csv": args.input_csv,
        "n": int(len(full)),
        "target": "C_obs = log(theta_E_obs/theta_E_baryon)",
        "settings": vars(args),
        "best_feature": best_feature,
        "best_by_LOO": comp.iloc[0].to_dict() if len(comp) else {},
        "by_variant": by_variant.to_dict(orient="records"),
        "by_variant_amp": by_variant_amp.to_dict(orient="records"),
        "bootstrap": boot,
        "files": {
            "features": "absolute_amplitude_v2_features.csv",
            "screen": "absolute_amplitude_v2_screen.csv",
            "by_variant_amp": "absolute_amplitude_v2_by_variant_amp.csv",
            "by_variant": "absolute_amplitude_v2_by_variant.csv",
            "model_summary": "absolute_amplitude_v2_model_summary.csv",
        },
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 100)
    print("ABSOLUTE AMPLITUDE V2 SUMMARY")
    print("=" * 100)
    print(comp.to_string(index=False))

    print("\nBY VARIANT")
    print(by_variant.to_string(index=False))

    print("\nBY VARIANT + AMP")
    print(by_variant_amp.head(20).to_string(index=False))

    print("\nTOP FEATURES")
    print(screen.head(40).to_string(index=False))

    print("\nWrote:", os.path.join(args.output_dir, "summary.json"))


if __name__ == "__main__":
    main()

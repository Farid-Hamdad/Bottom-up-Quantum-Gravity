#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 10 — Sparse Newtonian model comparison v4.

Goal:
    Test whether the point-source Green function on a 3D lattice is better
    described by:

        Newton fixed model:
            Phi(r) = C + A/r

    or by:

        Free power-law model:
            Phi(r) = C + A r^{-alpha}

    The fit is restricted to an intermediate physical window:

        r in [f_min L, f_max L]

    default:
        f_min = 0.10
        f_max = 0.40

This version uses sparse matrices and conjugate gradient, so it can test
larger 3D lattices than the dense v2/v3 scripts.

Outputs:
    - newton_compare_summary.csv
    - newton_compare_summary.json
    - profiles/profile_L*.csv
    - figures/*.png
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse.linalg import cg
from scipy.optimize import curve_fit


# ============================================================
# Utilities
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def idx3(x, y, z, L):
    return x * L * L + y * L + z


def coords_from_index(i, L):
    x = i // (L * L)
    rem = i % (L * L)
    y = rem // L
    z = rem % L
    return x, y, z


def center_index(L):
    c = L // 2
    return idx3(c, c, c, L)


# ============================================================
# Sparse 3D lattice Laplacian
# ============================================================

def build_sparse_laplacian_3d(L, periodic=True):
    """
    Build sparse combinatorial Laplacian for a 3D cubic lattice.

    Nodes: N = L^3.
    Nearest-neighbor edges.
    """
    N = L ** 3

    rows = []
    cols = []
    data = []

    degree = np.zeros(N, dtype=float)

    directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = idx3(x, y, z, L)

                for dx, dy, dz in directions:
                    x2 = x + dx
                    y2 = y + dy
                    z2 = z + dz

                    if periodic:
                        x2 %= L
                        y2 %= L
                        z2 %= L
                    else:
                        if x2 >= L or y2 >= L or z2 >= L:
                            continue

                    j = idx3(x2, y2, z2, L)

                    # undirected edge i-j
                    rows.extend([i, j])
                    cols.extend([j, i])
                    data.extend([-1.0, -1.0])
                    degree[i] += 1.0
                    degree[j] += 1.0

    # diagonal
    rows.extend(range(N))
    cols.extend(range(N))
    data.extend(degree)

    Lmat = sparse.csr_matrix((data, (rows, cols)), shape=(N, N))
    return Lmat


# ============================================================
# Solve point-source Poisson
# ============================================================

def solve_point_source_sparse(Lmat, source_index, rtol=1e-10, atol=1e-12, maxiter=20000):
    """
    Solve:

        L Phi = delta_i0 - 1/N

    Since the RHS has zero mean, the singular Laplacian system is compatible.
    CG returns a solution in the orthogonal complement of the constant mode
    up to numerical precision. We remove the mean afterward.
    """
    N = Lmat.shape[0]

    b = np.full(N, -1.0 / N, dtype=float)
    b[source_index] += 1.0

    try:
        phi, info = cg(Lmat, b, rtol=rtol, atol=atol, maxiter=maxiter)
    except TypeError:
        # Compatibility with older scipy versions
        phi, info = cg(Lmat, b, tol=rtol, maxiter=maxiter)

    phi = np.asarray(phi, dtype=float)
    phi -= np.mean(phi)

    return phi, info


# ============================================================
# Torus distance and radial profile
# ============================================================

def torus_distance_to_center(L, periodic=True):
    """
    Compute r_i from each node to central node.

    If periodic=True, use torus distance:
        dx = min(|x-x0|, L-|x-x0|)
    """
    N = L ** 3
    c = L // 2
    r = np.zeros(N, dtype=float)

    for i in range(N):
        x, y, z = coords_from_index(i, L)

        dx = abs(x - c)
        dy = abs(y - c)
        dz = abs(z - c)

        if periodic:
            dx = min(dx, L - dx)
            dy = min(dy, L - dy)
            dz = min(dz, L - dz)

        r[i] = np.sqrt(dx * dx + dy * dy + dz * dz)

    return r


def radial_profile(phi, r, bins=40):
    mask = r > 0
    rr = r[mask]
    yy = phi[mask]

    r_min = float(np.min(rr))
    r_max = float(np.max(rr))

    edges = np.linspace(r_min, r_max, bins + 1)

    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (rr >= a) & (rr < b)
        if b == edges[-1]:
            m = (rr >= a) & (rr <= b)

        vals = yy[m]
        dists = rr[m]

        rows.append({
            "r_min": float(a),
            "r_max": float(b),
            "r_mid": float(0.5 * (a + b)),
            "r_mean": float(np.mean(dists)) if vals.size else np.nan,
            "phi_mean": float(np.mean(vals)) if vals.size else np.nan,
            "phi_std": float(np.std(vals)) if vals.size else np.nan,
            "count": int(vals.size)
        })

    return pd.DataFrame(rows)


# ============================================================
# Models and fit metrics
# ============================================================

def model_newton(r, C, A):
    return C + A / r


def model_free(r, C, A, alpha):
    return C + A * np.power(r, -alpha)


def r2_score(y, yhat):
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 1e-15:
        return np.nan
    return 1.0 - ss_res / ss_tot


def information_criteria(y, yhat, p):
    """
    AIC/BIC for Gaussian residuals up to constants:
        AIC = n log(RSS/n) + 2p
        BIC = n log(RSS/n) + p log(n)
    """
    y = np.asarray(y, dtype=float)
    yhat = np.asarray(yhat, dtype=float)

    n = len(y)
    rss = float(np.sum((y - yhat) ** 2))
    rss = max(rss, 1e-300)

    aic = n * np.log(rss / n) + 2 * p
    bic = n * np.log(rss / n) + p * np.log(n)

    return rss, aic, bic


def fit_models(profile, L, rmin_frac=0.10, rmax_frac=0.40, min_count=5):
    df = profile.copy()
    df = df[df["count"] >= min_count]
    df = df[np.isfinite(df["r_mean"]) & np.isfinite(df["phi_mean"])]

    r_low = rmin_frac * L
    r_high = rmax_frac * L

    df = df[(df["r_mean"] >= r_low) & (df["r_mean"] <= r_high)]

    if len(df) < 4:
        return None, df

    r = df["r_mean"].values.astype(float)
    y = df["phi_mean"].values.astype(float)

    # Newton fixed model: C + A/r
    X = np.vstack([np.ones_like(r), 1.0 / r]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    C_newt, A_newt = beta
    yhat_newt = model_newton(r, C_newt, A_newt)

    rss_newt, aic_newt, bic_newt = information_criteria(y, yhat_newt, p=2)
    r2_newt = r2_score(y, yhat_newt)

    # Free model: C + A r^-alpha
    try:
        p0 = [C_newt, A_newt, 1.0]
        popt, _ = curve_fit(
            model_free,
            r,
            y,
            p0=p0,
            bounds=([-np.inf, -np.inf, -5.0], [np.inf, np.inf, 5.0]),
            maxfev=50000
        )
        C_free, A_free, alpha_free = popt
        yhat_free = model_free(r, C_free, A_free, alpha_free)

        rss_free, aic_free, bic_free = information_criteria(y, yhat_free, p=3)
        r2_free = r2_score(y, yhat_free)

    except Exception as e:
        C_free, A_free, alpha_free = np.nan, np.nan, np.nan
        yhat_free = np.full_like(y, np.nan)
        rss_free, aic_free, bic_free, r2_free = np.nan, np.nan, np.nan, np.nan

    delta_aic_free_minus_newt = aic_free - aic_newt
    delta_bic_free_minus_newt = bic_free - bic_newt

    if np.isfinite(delta_bic_free_minus_newt):
        if delta_bic_free_minus_newt > 2:
            preferred = "newton"
        elif delta_bic_free_minus_newt < -2:
            preferred = "free"
        else:
            preferred = "comparable"
    else:
        preferred = "undefined"

    result = {
        "n_fit": int(len(df)),
        "r_window_min": float(r_low),
        "r_window_max": float(r_high),

        "C_newton": float(C_newt),
        "A_newton": float(A_newt),
        "R2_newton": float(r2_newt),
        "RSS_newton": float(rss_newt),
        "AIC_newton": float(aic_newt),
        "BIC_newton": float(bic_newt),

        "C_free": float(C_free),
        "A_free": float(A_free),
        "alpha_free": float(alpha_free),
        "R2_free": float(r2_free),
        "RSS_free": float(rss_free),
        "AIC_free": float(aic_free),
        "BIC_free": float(bic_free),

        "delta_AIC_free_minus_newton": float(delta_aic_free_minus_newt),
        "delta_BIC_free_minus_newton": float(delta_bic_free_minus_newt),
        "preferred_by_BIC": preferred
    }

    fit_df = pd.DataFrame({
        "r": r,
        "phi": y,
        "phi_newton": yhat_newt,
        "phi_free": yhat_free
    })

    return result, fit_df


# ============================================================
# Figures
# ============================================================

def make_figures(profile, fit_df, result, L, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    # Full radial profile
    df = profile[profile["count"] > 0].copy()

    plt.figure(figsize=(8, 5))
    plt.errorbar(df["r_mean"], df["phi_mean"], yerr=df["phi_std"], fmt="o-", capsize=3)
    plt.axvspan(result["r_window_min"], result["r_window_max"], alpha=0.15, label="fit window")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\Phi(r)$")
    plt.title(f"Point-source Green profile — L={L}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_profile_L{L}.png"), dpi=250)
    plt.close()

    # Model comparison in fit window
    plt.figure(figsize=(8, 5))
    plt.plot(fit_df["r"], fit_df["phi"], "o", label="data")
    plt.plot(fit_df["r"], fit_df["phi_newton"], "-", label=fr"Newton $\alpha=1$, $R^2={result['R2_newton']:.4f}$")
    plt.plot(fit_df["r"], fit_df["phi_free"], "--", label=fr"free $\alpha={result['alpha_free']:.3f}$, $R^2={result['R2_free']:.4f}$")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\Phi(r)$")
    plt.title(f"Newton fixed vs free power-law — L={L}")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_newton_vs_free_L{L}.png"), dpi=250)
    plt.close()

    # Residuals
    plt.figure(figsize=(8, 5))
    plt.plot(fit_df["r"], fit_df["phi"] - fit_df["phi_newton"], "o-", label="Newton residual")
    plt.plot(fit_df["r"], fit_df["phi"] - fit_df["phi_free"], "s-", label="free residual")
    plt.axhline(0, linewidth=1)
    plt.xlabel(r"$r$")
    plt.ylabel("residual")
    plt.title(f"Residuals — L={L}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_residuals_L{L}.png"), dpi=250)
    plt.close()


def make_summary_figure(summary_df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    plt.figure(figsize=(8, 5))
    plt.plot(summary_df["L"], summary_df["alpha_free"], "o-", label=r"free $\alpha$")
    plt.axhline(1.0, linestyle="--", linewidth=1, label="Newton alpha=1")
    plt.xlabel("L")
    plt.ylabel(r"$\alpha$")
    plt.title("Free exponent vs lattice size")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_alpha_free_vs_L.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(summary_df["L"], summary_df["delta_BIC_free_minus_newton"], "o-")
    plt.axhline(0, linewidth=1)
    plt.axhline(2, linestyle="--", linewidth=1, label="Newton preferred")
    plt.axhline(-2, linestyle="--", linewidth=1, label="Free preferred")
    plt.xlabel("L")
    plt.ylabel(r"$\Delta BIC = BIC_{\rm free}-BIC_{\rm Newt}$")
    plt.title("Model comparison: Newton fixed vs free")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_delta_bic_vs_L.png"), dpi=250)
    plt.close()


# ============================================================
# One L run
# ============================================================

def run_one_L(L, args):
    print("=" * 100)
    print(f"Running L={L}, N={L**3}")
    print("=" * 100)

    outL = os.path.join(args.output_dir, f"L{L}")
    ensure_dir(outL)

    Lmat = build_sparse_laplacian_3d(L, periodic=args.periodic)
    i0 = center_index(L)

    phi, info = solve_point_source_sparse(
        Lmat,
        i0,
        rtol=args.rtol,
        atol=args.atol,
        maxiter=args.maxiter
    )

    if info != 0:
        print(f"WARNING: CG info={info} for L={L}")

    r = torus_distance_to_center(L, periodic=args.periodic)

    profile = radial_profile(phi, r, bins=args.bins)
    profile_path = os.path.join(outL, f"profile_L{L}.csv")
    profile.to_csv(profile_path, index=False)

    result, fit_df = fit_models(
        profile,
        L=L,
        rmin_frac=args.rmin_frac,
        rmax_frac=args.rmax_frac,
        min_count=args.min_count
    )

    if result is None:
        raise RuntimeError(f"Not enough fit points for L={L}.")

    fit_path = os.path.join(outL, f"fit_points_L{L}.csv")
    fit_df.to_csv(fit_path, index=False)

    make_figures(profile, fit_df, result, L, outL)

    row = {
        "L": L,
        "N": L ** 3,
        "periodic": bool(args.periodic),
        "cg_info": int(info),
        "profile_csv": profile_path,
        "fit_points_csv": fit_path,
    }
    row.update(result)

    print(
        f"L={L} "
        f"alpha_free={row['alpha_free']:.6g} "
        f"R2_newton={row['R2_newton']:.6g} "
        f"R2_free={row['R2_free']:.6g} "
        f"dBIC={row['delta_BIC_free_minus_newton']:.6g} "
        f"preferred={row['preferred_by_BIC']}"
    )

    return row


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--L-list", type=int, nargs="+", default=[13, 17, 21])
    parser.add_argument("--output-dir", default="papers/paper10_limite_newtonienne_propagateur/results/sparse_newton_compare_v4")

    parser.add_argument("--periodic", action="store_true")

    parser.add_argument("--bins", type=int, default=50)
    parser.add_argument("--min-count", type=int, default=10)

    parser.add_argument("--rmin-frac", type=float, default=0.10)
    parser.add_argument("--rmax-frac", type=float, default=0.40)

    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--maxiter", type=int, default=50000)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    rows = []
    for L in args.L_list:
        rows.append(run_one_L(L, args))

    df = pd.DataFrame(rows)
    summary_csv = os.path.join(args.output_dir, "newton_compare_summary.csv")
    df.to_csv(summary_csv, index=False)

    summary = {
        "title": "BuP Paper 10 — Sparse Newtonian model comparison v4",
        "description": "Compares fixed Newton model C+A/r against free power law C+A r^-alpha.",
        "parameters": vars(args),
        "results": rows,
        "files": {
            "summary_csv": "newton_compare_summary.csv",
            "figures_dir": "figures/",
            "per_L_dirs": "L*/"
        }
    }

    summary_json = os.path.join(args.output_dir, "newton_compare_summary.json")
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    make_summary_figure(df, args.output_dir)

    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(df[[
        "L",
        "N",
        "alpha_free",
        "R2_newton",
        "R2_free",
        "AIC_newton",
        "AIC_free",
        "BIC_newton",
        "BIC_free",
        "delta_AIC_free_minus_newton",
        "delta_BIC_free_minus_newton",
        "preferred_by_BIC"
    ]].to_string(index=False))

    print("=" * 100)
    print("Files:")
    print(summary_csv)
    print(summary_json)
    print("DONE")


if __name__ == "__main__":
    main()

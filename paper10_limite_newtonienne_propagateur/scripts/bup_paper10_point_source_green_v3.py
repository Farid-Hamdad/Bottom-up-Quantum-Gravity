#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BuP Paper 10 — Point-source Green function test v3 with torus distance.

Goal:
    Correctly test the Newtonian limit of the graph Green function.

Instead of averaging all entries of L^+, this script solves a point-source
Poisson problem:

    L Phi = delta_i0 - 1/N

where i0 is the central node.

Then it measures the radial potential Phi(r) around the source and fits:

    Phi(r) = C + A r^{-alpha}

For a 3D lattice, the Newtonian target is:

    alpha ≈ 1

For 2D, the Green function is logarithmic:
    Phi(r) ≈ C + A log(r)

For 1D, the Green function grows approximately linearly with r.

Outputs:
    synthetic_point_source_summary.csv
    synthetic_point_source_summary.json
    profiles/profile_<graph>.csv
    figures/*.png
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.linalg import expm
from scipy.stats import linregress


# ============================================================
# Utilities
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def lattice_index(coords, shape):
    idx = 0
    mul = 1
    for c, s in zip(reversed(coords), reversed(shape)):
        idx += c * mul
        mul *= s
    return idx


def build_lattice_adjacency(dim, Lsize, periodic=False):
    shape = tuple([Lsize] * dim)
    N = Lsize ** dim
    A = np.zeros((N, N), dtype=float)

    for coords in np.ndindex(shape):
        i = lattice_index(coords, shape)

        for axis in range(dim):
            for step in [-1, 1]:
                nb = list(coords)
                nb[axis] += step

                if periodic:
                    nb[axis] %= Lsize
                else:
                    if nb[axis] < 0 or nb[axis] >= Lsize:
                        continue

                j = lattice_index(tuple(nb), shape)
                A[i, j] = 1.0

    A = np.maximum(A, A.T)
    np.fill_diagonal(A, 0.0)
    return A


def lattice_coordinates(dim, Lsize):
    shape = tuple([Lsize] * dim)
    coords = np.array(list(np.ndindex(shape)), dtype=float)
    return coords


def center_node(dim, Lsize):
    shape = tuple([Lsize] * dim)
    c = tuple([Lsize // 2] * dim)
    return lattice_index(c, shape)


def combinatorial_laplacian(A):
    deg = np.sum(A, axis=1)
    return np.diag(deg) - A


def solve_point_source(L, source_index, eig_tol=1e-10):
    """
    Solve L Phi = delta_i0 - 1/N via spectral pseudo-inverse.
    """
    N = L.shape[0]

    b = np.zeros(N)
    b[source_index] = 1.0
    b -= 1.0 / N

    evals, evecs = np.linalg.eigh(L)
    inv = np.zeros_like(evals)
    good = evals > eig_tol
    inv[good] = 1.0 / evals[good]

    Phi = evecs @ (inv * (evecs.T @ b))
    Phi -= np.mean(Phi)

    return Phi, evals


def radial_profile_point_source(
    Phi,
    coords,
    source_index,
    bins=16,
    r_min=None,
    r_max=None,
    periodic=False,
    Lsize=None
):
    """
    Radial profile around the point source.

    If periodic=True, use torus distance:

        r = sqrt(sum_a min(|dx_a|, L - |dx_a|)^2)

    This is essential when the adjacency has periodic boundary conditions.
    """
    x0 = coords[source_index]

    diff = np.abs(coords - x0)

    if periodic:
        if Lsize is None:
            raise ValueError("Lsize must be provided when periodic=True.")
        diff = np.minimum(diff, Lsize - diff)

    r = np.linalg.norm(diff, axis=1)

    mask = r > 0
    r = r[mask]
    phi = Phi[mask]

    if r_min is None:
        r_min = float(np.min(r))
    if r_max is None:
        r_max = float(np.max(r))

    edges = np.linspace(r_min, r_max, bins + 1)

    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        if b == edges[-1]:
            m = (r >= a) & (r <= b)

        vals = phi[m]
        rr = r[m]

        rows.append({
            "r_min": float(a),
            "r_max": float(b),
            "r_mid": float(0.5 * (a + b)),
            "r_mean": float(np.mean(rr)) if vals.size else np.nan,
            "phi_mean": float(np.mean(vals)) if vals.size else np.nan,
            "phi_std": float(np.std(vals)) if vals.size else np.nan,
            "phi_abs_mean": float(np.mean(np.abs(vals))) if vals.size else np.nan,
            "count": int(vals.size)
        })

    return pd.DataFrame(rows)


# ============================================================
# Fits
# ============================================================

def model_power(r, C, A, alpha):
    return C + A * np.power(r, -alpha)


def model_log(r, C, A):
    return C + A * np.log(r)


def model_linear(r, C, A):
    return C + A * r


def r2_score(y, yhat):
    y = np.asarray(y, dtype=float)
    yhat = np.asarray(yhat, dtype=float)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 1e-15:
        return np.nan
    return 1.0 - ss_res / ss_tot


def fit_point_source_profile(profile, dim, min_count=5, drop_first=1, drop_last=1):
    df = profile.copy()
    df = df[df["count"] >= min_count]
    df = df[np.isfinite(df["r_mean"]) & np.isfinite(df["phi_mean"])]
    df = df[df["r_mean"] > 0]

    if drop_first > 0 and len(df) > drop_first:
        df = df.iloc[drop_first:]

    if drop_last > 0 and len(df) > drop_last:
        df = df.iloc[:-drop_last]

    if len(df) < 4:
        return {
            "fit_type": "insufficient_data",
            "alpha": np.nan,
            "C": np.nan,
            "A": np.nan,
            "r2": np.nan,
            "n_fit": int(len(df))
        }

    r = df["r_mean"].values.astype(float)
    y = df["phi_mean"].values.astype(float)

    result = {
        "n_fit": int(len(df))
    }

    # Power-law fit C + A r^-alpha
    try:
        C0 = float(np.median(y[-max(2, len(y)//3):]))
        A0 = float(y[0] - C0)
        alpha0 = 1.0 if dim == 3 else max(dim - 2, 0.1)

        popt, pcov = curve_fit(
            model_power,
            r,
            y,
            p0=[C0, A0, alpha0],
            bounds=([-np.inf, -np.inf, -5.0], [np.inf, np.inf, 5.0]),
            maxfev=20000
        )

        yhat = model_power(r, *popt)
        result.update({
            "fit_type": "power",
            "C": float(popt[0]),
            "A": float(popt[1]),
            "alpha": float(popt[2]),
            "r2": float(r2_score(y, yhat)),
        })
    except Exception as e:
        result.update({
            "fit_type": "power_failed",
            "C": np.nan,
            "A": np.nan,
            "alpha": np.nan,
            "r2": np.nan,
            "error": str(e)
        })

    # Log fit for 2D comparison
    try:
        popt_log, _ = curve_fit(model_log, r, y, p0=[np.mean(y), 1.0], maxfev=20000)
        yhat_log = model_log(r, *popt_log)
        result["log_C"] = float(popt_log[0])
        result["log_A"] = float(popt_log[1])
        result["log_r2"] = float(r2_score(y, yhat_log))
    except Exception:
        result["log_C"] = np.nan
        result["log_A"] = np.nan
        result["log_r2"] = np.nan

    # Linear fit for 1D comparison
    try:
        lr = linregress(r, y)
        yhat_lin = lr.intercept + lr.slope * r
        result["linear_C"] = float(lr.intercept)
        result["linear_A"] = float(lr.slope)
        result["linear_r2"] = float(r2_score(y, yhat_lin))
    except Exception:
        result["linear_C"] = np.nan
        result["linear_A"] = np.nan
        result["linear_r2"] = np.nan

    return result


# ============================================================
# Spectral dimension
# ============================================================

def heat_trace_spectral_dimension(L, tau_min, tau_max, tau_points):
    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_points)
    Z = []

    for t in taus:
        K = expm(-t * L)
        Z.append(float(np.trace(K)))

    Z = np.asarray(Z)
    logt = np.log(taus)
    logZ = np.log(np.maximum(Z, 1e-300))

    slope = np.gradient(logZ, logt)
    ds = -2.0 * slope

    n = len(ds)
    lo = max(1, n // 4)
    hi = min(n - 1, 3 * n // 4)
    ds_eff = float(np.nanmedian(ds[lo:hi]))

    return taus, Z, ds, ds_eff


# ============================================================
# Figures
# ============================================================

def make_graph_figures(profile, fit, dim, graph_name, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)
    safe = graph_name.replace(" ", "_")

    df = profile.copy()
    df = df[df["count"] > 0]
    df = df[np.isfinite(df["r_mean"]) & np.isfinite(df["phi_mean"])]

    r = df["r_mean"].values
    y = df["phi_mean"].values

    # Radial potential
    plt.figure(figsize=(8, 5))
    plt.errorbar(df["r_mean"], df["phi_mean"], yerr=df["phi_std"], fmt="o-", capsize=3)
    plt.xlabel(r"$r$")
    plt.ylabel(r"$\Phi(r)$")
    plt.title(f"Point-source Green profile — {graph_name}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_point_source_profile_{safe}.png"), dpi=250)
    plt.close()

    # Fit
    plt.figure(figsize=(8, 5))
    plt.plot(r, y, "o", label="data")

    rr = np.linspace(np.min(r), np.max(r), 300)

    if np.isfinite(fit.get("alpha", np.nan)):
        yy = model_power(rr, fit["C"], fit["A"], fit["alpha"])
        plt.plot(rr, yy, "-", label=fr"power $\alpha={fit['alpha']:.3f}$, $R^2={fit['r2']:.3f}$")

    if np.isfinite(fit.get("log_r2", np.nan)):
        yy_log = model_log(rr, fit["log_C"], fit["log_A"])
        plt.plot(rr, yy_log, "--", label=fr"log fit $R^2={fit['log_r2']:.3f}$")

    if np.isfinite(fit.get("linear_r2", np.nan)):
        yy_lin = model_linear(rr, fit["linear_C"], fit["linear_A"])
        plt.plot(rr, yy_lin, ":", label=fr"linear fit $R^2={fit['linear_r2']:.3f}$")

    plt.xlabel(r"$r$")
    plt.ylabel(r"$\Phi(r)$")
    plt.title(f"Point-source Green fits — {graph_name}")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_point_source_fit_{safe}.png"), dpi=250)
    plt.close()

    # Log-log absolute shifted profile for visual check
    plt.figure(figsize=(8, 5))

    if np.isfinite(fit.get("C", np.nan)):
        y_shift = np.abs(y - fit["C"])
        plt.loglog(r, y_shift, "o-", label=r"$|\Phi(r)-C|$")
        if np.isfinite(fit.get("alpha", np.nan)):
            yy = np.abs(fit["A"]) * rr ** (-fit["alpha"])
            plt.loglog(rr, yy, "-", label=fr"$r^{{-{fit['alpha']:.3f}}}$")
    else:
        plt.loglog(r, np.abs(y), "o-", label=r"$|\Phi(r)|$")

    plt.xlabel(r"$r$")
    plt.ylabel(r"$|\Phi-C|$")
    plt.title(f"Log-log Newtonian check — {graph_name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"fig_point_source_loglog_{safe}.png"), dpi=250)
    plt.close()


def make_summary_figures(summary_df, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    ensure_dir(fig_dir)

    plt.figure(figsize=(8, 5))
    plt.plot(summary_df["dimension"], summary_df["alpha"], "o-", label=r"measured $\alpha$")
    plt.plot(summary_df["dimension"], summary_df["expected_alpha"], "s--", label=r"expected $d-2$")
    plt.axhline(1.0, linestyle=":", linewidth=1, label=r"Newtonian $\alpha=1$")
    plt.xlabel("dimension géométrique")
    plt.ylabel(r"exposant $\alpha$")
    plt.title("Point-source Green exponent")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_point_source_alpha_vs_dimension.png"), dpi=250)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(summary_df["dimension"], summary_df["ds_eff"], "o-", label=r"$d_s$ measured")
    plt.plot(summary_df["dimension"], summary_df["dimension"], "s--", label="geometric dimension")
    plt.xlabel("dimension géométrique")
    plt.ylabel(r"$d_s$")
    plt.title("Spectral dimension calibration")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fig_point_source_ds_vs_dimension.png"), dpi=250)
    plt.close()


# ============================================================
# One graph
# ============================================================

def run_graph(dim, Lsize, args):
    graph_name = f"{dim}D_L{Lsize}"

    print(f"===== {graph_name} =====")

    A = build_lattice_adjacency(dim, Lsize, periodic=args.periodic)
    L = combinatorial_laplacian(A)
    coords = lattice_coordinates(dim, Lsize)
    i0 = center_node(dim, Lsize)

    Phi, evals = solve_point_source(L, i0, eig_tol=args.eig_tol)

    profile = radial_profile_point_source(
        Phi,
        coords,
        i0,
        bins=args.bins,
        periodic=args.periodic,
        Lsize=Lsize
    )

    fit = fit_point_source_profile(
        profile,
        dim=dim,
        min_count=args.min_count,
        drop_first=args.drop_first,
        drop_last=args.drop_last
    )

    taus, Z, ds_curve, ds_eff = heat_trace_spectral_dimension(
        L,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
        tau_points=args.tau_points
    )

    profiles_dir = os.path.join(args.output_dir, "profiles")
    ensure_dir(profiles_dir)

    profile_path = os.path.join(profiles_dir, f"profile_point_source_{graph_name}.csv")
    profile.to_csv(profile_path, index=False)

    heat_path = os.path.join(profiles_dir, f"spectral_dimension_{graph_name}.csv")
    pd.DataFrame({
        "tau": taus,
        "Z": Z,
        "d_s": ds_curve
    }).to_csv(heat_path, index=False)

    make_graph_figures(profile, fit, dim, graph_name, args.output_dir)

    expected_alpha = dim - 2

    row = {
        "graph": graph_name,
        "dimension": dim,
        "Lsize": Lsize,
        "N": int(A.shape[0]),
        "periodic": bool(args.periodic),
        "source_index": int(i0),
        "ds_eff": float(ds_eff),
        "expected_alpha": float(expected_alpha),
        "alpha": float(fit.get("alpha", np.nan)),
        "alpha_minus_expected": float(fit.get("alpha", np.nan) - expected_alpha) if np.isfinite(fit.get("alpha", np.nan)) else np.nan,
        "distance_to_newtonian_alpha_1": float(abs(fit.get("alpha", np.nan) - 1.0)) if np.isfinite(fit.get("alpha", np.nan)) else np.nan,
        "power_r2": float(fit.get("r2", np.nan)),
        "log_r2": float(fit.get("log_r2", np.nan)),
        "linear_r2": float(fit.get("linear_r2", np.nan)),
        "fit_type": fit.get("fit_type", ""),
        "fit_C": float(fit.get("C", np.nan)),
        "fit_A": float(fit.get("A", np.nan)),
        "n_fit": int(fit.get("n_fit", 0)),
        "lambda_2": float(np.sort(evals)[1]) if len(evals) > 1 else np.nan,
        "lambda_max": float(np.max(evals)),
        "zero_modes_estimate": int(np.sum(evals < args.eig_tol)),
        "profile_csv": profile_path,
        "spectral_dimension_csv": heat_path
    }

    print(
        f"N={row['N']} ds={row['ds_eff']:.4g} "
        f"alpha={row['alpha']:.4g} expected={expected_alpha:.4g} "
        f"R2power={row['power_r2']:.4g} R2log={row['log_r2']:.4g} R2lin={row['linear_r2']:.4g}"
    )

    return row


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--output-dir", default="papers/paper10_limite_newtonienne_propagateur/results/point_source_green_v3")

    parser.add_argument("--L1", type=int, default=120)
    parser.add_argument("--L2", type=int, default=31)
    parser.add_argument("--L3", type=int, default=11)

    parser.add_argument("--periodic", action="store_true")

    parser.add_argument("--bins", type=int, default=18)
    parser.add_argument("--min-count", type=int, default=4)
    parser.add_argument("--drop-first", type=int, default=1)
    parser.add_argument("--drop-last", type=int, default=2)

    parser.add_argument("--tau-min", type=float, default=0.01)
    parser.add_argument("--tau-max", type=float, default=100.0)
    parser.add_argument("--tau-points", type=int, default=60)

    parser.add_argument("--eig-tol", type=float, default=1e-10)

    args = parser.parse_args()
    ensure_dir(args.output_dir)

    print("=" * 100)
    print("BuP Paper 10 — Point-source Green function test v3 — torus distance")
    print("=" * 100)
    print(f"output: {args.output_dir}")
    print("-" * 100)

    rows = []
    rows.append(run_graph(1, args.L1, args))
    rows.append(run_graph(2, args.L2, args))
    rows.append(run_graph(3, args.L3, args))

    df = pd.DataFrame(rows)
    summary_csv = os.path.join(args.output_dir, "point_source_green_summary.csv")
    df.to_csv(summary_csv, index=False)

    summary = {
        "title": "BuP Paper 10 — Point-source Green function test v3 with torus distance",
        "description": "Tests Newtonian limit using L Phi = delta_i0 - 1/N on 1D, 2D, 3D lattices with torus distance for periodic BC.",
        "parameters": {
            "L1": args.L1,
            "L2": args.L2,
            "L3": args.L3,
            "periodic": args.periodic,
            "bins": args.bins,
            "min_count": args.min_count,
            "drop_first": args.drop_first,
            "drop_last": args.drop_last,
            "tau_min": args.tau_min,
            "tau_max": args.tau_max,
            "tau_points": args.tau_points
        },
        "results": rows,
        "files": {
            "summary_csv": "point_source_green_summary.csv",
            "profiles_dir": "profiles/",
            "figures_dir": "figures/"
        }
    }

    summary_json = os.path.join(args.output_dir, "point_source_green_summary.json")
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    make_summary_figures(df, args.output_dir)

    print("-" * 100)
    print("SUMMARY")
    print("-" * 100)
    print(df[[
        "graph",
        "N",
        "dimension",
        "ds_eff",
        "alpha",
        "expected_alpha",
        "alpha_minus_expected",
        "power_r2",
        "log_r2",
        "linear_r2"
    ]].to_string(index=False))
    print("-" * 100)
    print("Files:")
    print(summary_csv)
    print(summary_json)
    print("DONE")


if __name__ == "__main__":
    main()

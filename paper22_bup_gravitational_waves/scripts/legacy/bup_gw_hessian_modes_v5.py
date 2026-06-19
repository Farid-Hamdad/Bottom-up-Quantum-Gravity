#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 21 — BuP gravitational-wave prototype v5

Hessian modes of an effective BuP action.

v5 moves from kinematics to dynamics. We linearize an effective BuP action
around a background entanglement graph W0 and compute Hessian eigenmodes.

Controlled reduction:
    W_ij(phi) = W_ij^0 * exp((phi_i + phi_j)/2)

Then:
    H_ab = d^2 S_BuP[phi] / dphi_a dphi_b |_{phi=0}

Action:
    S_spec   = Tr[(L_norm + mu I)^(-beta)]
    S_loc    = lambda * sum_ij W_ij d_ij^2
    S_smooth = eta * sum_ij W0_ij (phi_i - phi_j)^2
    S_mass   = m2 * sum_i phi_i^2

Run:
python3 bup_gw_hessian_modes_v5.py \
  --N-side 11 \
  --beta 1.0 \
  --mu 0.05 \
  --lambda-locality 0.25 \
  --eta-smooth 1.0 \
  --m2 1e-3 \
  --output-dir results/hessian_modes_v5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh


def build_grid(n_side: int, extent: float):
    xs = np.linspace(-extent, extent, n_side)
    yy, xx = np.meshgrid(xs, xs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel()])


def pairwise_dist(coords):
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def base_weights(coords, ell, cutoff):
    R = pairwise_dist(coords)
    W = np.exp(-R / ell)
    W[R > cutoff] = 0.0
    np.fill_diagonal(W, 0.0)
    return W, R


def normalized_laplacian(W):
    deg = W.sum(axis=1)
    invsqrt = np.zeros_like(deg)
    mask = deg > 0
    invsqrt[mask] = 1.0 / np.sqrt(deg[mask])
    S = (invsqrt[:, None] * W) * invsqrt[None, :]
    return np.eye(W.shape[0]) - S


def combinatorial_laplacian(W):
    return np.diag(W.sum(axis=1)) - W


def deform_W(W0, phi):
    W = W0 * np.exp(0.5 * (phi[:, None] + phi[None, :]))
    np.fill_diagonal(W, 0.0)
    return W


def spectral_action(W, beta, mu):
    L = normalized_laplacian(W)
    vals = np.linalg.eigvalsh(L + mu * np.eye(L.shape[0]))
    vals = np.clip(vals, 1e-12, None)
    return float(np.sum(vals ** (-beta)))


def action_phi(phi, W0, R, beta, mu, lambda_locality, eta_smooth, m2):
    W = deform_W(W0, phi)
    S_spec = spectral_action(W, beta=beta, mu=mu)
    S_loc = 0.5 * lambda_locality * float(np.sum(W * (R ** 2)))
    dphi = phi[:, None] - phi[None, :]
    S_smooth = 0.5 * eta_smooth * float(np.sum(W0 * dphi * dphi))
    S_mass = 0.5 * m2 * float(np.sum(phi * phi))
    return S_spec + S_loc + S_smooth + S_mass


def finite_difference_hessian(W0, R, beta, mu, lambda_locality, eta_smooth, m2, eps):
    n = W0.shape[0]
    H = np.zeros((n, n), dtype=float)
    phi0 = np.zeros(n, dtype=float)
    S0 = action_phi(phi0, W0, R, beta, mu, lambda_locality, eta_smooth, m2)

    for i in range(n):
        ep = np.zeros(n); ep[i] = eps
        Sp = action_phi(ep, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
        Sm = action_phi(-ep, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
        H[i, i] = (Sp - 2*S0 + Sm) / (eps * eps)

    for i in range(n):
        ei = np.zeros(n); ei[i] = eps
        for j in range(i + 1, n):
            ej = np.zeros(n); ej[j] = eps
            Spp = action_phi(ei + ej, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
            Spm = action_phi(ei - ej, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
            Smp = action_phi(-ei + ej, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
            Smm = action_phi(-ei - ej, W0, R, beta, mu, lambda_locality, eta_smooth, m2)
            val = (Spp - Spm - Smp + Smm) / (4 * eps * eps)
            H[i, j] = val
            H[j, i] = val

    return 0.5 * (H + H.T)


def mode_laplacian_energy(mode, L):
    denom = float(mode @ mode)
    if denom <= 0:
        return np.nan
    return float(mode @ L @ mode / denom)


def sign_fix(mode):
    idx = int(np.argmax(np.abs(mode)))
    return mode if mode[idx] >= 0 else -mode


def fit_dispersion(k2, omega2):
    x = np.asarray(k2, dtype=float)
    y = np.asarray(omega2, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan
    X = np.column_stack([x[ok], np.ones(ok.sum())])
    slope, intercept = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    pred = X @ np.array([slope, intercept])
    ss_res = np.sum((y[ok] - pred) ** 2)
    ss_tot = np.sum((y[ok] - y[ok].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(slope), float(intercept), float(r2)


def plot_mode(coords, mode, outpath, title):
    plt.figure(figsize=(5, 4.5))
    plt.scatter(coords[:, 0], coords[:, 1], c=mode, s=45)
    plt.colorbar(label="phi")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=170)
    plt.close()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--N-side", type=int, default=11)
    p.add_argument("--extent", type=float, default=1.0)
    p.add_argument("--ell", type=float, default=0.22)
    p.add_argument("--cutoff", type=float, default=0.45)
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--mu", type=float, default=0.05)
    p.add_argument("--lambda-locality", type=float, default=0.25)
    p.add_argument("--eta-smooth", type=float, default=1.0)
    p.add_argument("--m2", type=float, default=1e-3)
    p.add_argument("--fd-eps", type=float, default=1e-4)
    p.add_argument("--n-modes-save", type=int, default=12)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    coords = build_grid(args.N_side, args.extent)
    W0, R = base_weights(coords, args.ell, args.cutoff)
    n = W0.shape[0]

    L_comb = combinatorial_laplacian(W0)
    L_norm = normalized_laplacian(W0)
    lap_vals = np.linalg.eigvalsh(L_comb)

    H = finite_difference_hessian(
        W0, R,
        beta=args.beta,
        mu=args.mu,
        lambda_locality=args.lambda_locality,
        eta_smooth=args.eta_smooth,
        m2=args.m2,
        eps=args.fd_eps,
    )

    h_vals, h_vecs = eigh(H)
    idx = np.argsort(h_vals)
    h_vals = h_vals[idx]
    h_vecs = h_vecs[:, idx]

    rows = []
    for a in range(n):
        mode = sign_fix(h_vecs[:, a])
        k2_eff = mode_laplacian_energy(mode, L_comb)
        k2_norm = mode_laplacian_energy(mode, L_norm)
        omega2 = float(h_vals[a])
        rows.append({
            "mode": a,
            "omega2_hessian": omega2,
            "omega_hessian": float(np.sqrt(max(omega2, 0.0))) if omega2 >= 0 else np.nan,
            "k2_eff_comb_laplacian": float(k2_eff),
            "k_eff": float(np.sqrt(max(k2_eff, 0.0))) if np.isfinite(k2_eff) else np.nan,
            "k2_norm_laplacian": float(k2_norm),
            "min_component": float(mode.min()),
            "max_component": float(mode.max()),
            "mean_component": float(mode.mean()),
            "std_component": float(mode.std()),
        })
        if a < args.n_modes_save:
            plot_mode(
                coords,
                mode,
                outdir / f"fig_mode_{a:03d}.png",
                f"Hessian mode {a}: omega2={omega2:.4e}, k2={k2_eff:.4e}",
            )

    modes_df = pd.DataFrame(rows)
    modes_df.to_csv(outdir / "eigenmodes.csv", index=False)

    fit_df = modes_df[(modes_df["mode"] >= 1) & (modes_df["omega2_hessian"] > 0)].copy()
    low = fit_df.head(min(25, len(fit_df))).copy()
    slope, intercept, r2_disp = fit_dispersion(low["k2_eff_comb_laplacian"], low["omega2_hessian"])
    c_graph = np.sqrt(slope) if slope > 0 else np.nan
    low.to_csv(outdir / "dispersion.csv", index=False)

    summary = {
        "N_side": args.N_side,
        "N_nodes": int(n),
        "N_edges_nonzero_directed": int(np.sum(W0 > 0)),
        "N_edges_undirected": int(np.sum(np.triu(W0 > 0, k=1))),
        "beta": args.beta,
        "mu": args.mu,
        "lambda_locality": args.lambda_locality,
        "eta_smooth": args.eta_smooth,
        "m2": args.m2,
        "fd_eps": args.fd_eps,
        "hessian_min_eigenvalue": float(np.min(h_vals)),
        "hessian_max_eigenvalue": float(np.max(h_vals)),
        "hessian_negative_modes": int(np.sum(h_vals < -1e-8)),
        "hessian_near_zero_modes": int(np.sum(np.abs(h_vals) <= 1e-8)),
        "hessian_positive_modes": int(np.sum(h_vals > 1e-8)),
        "dispersion_fit_modes_used": int(len(low)),
        "dispersion_slope_c2": float(slope),
        "dispersion_intercept_mass2": float(intercept),
        "dispersion_r2": float(r2_disp),
        "c_graph_from_hessian": float(c_graph),
        "interpretation": "Hessian modes support wave dynamics if Hessian is mostly positive and omega^2 correlates with graph k^2 for low modes.",
    }

    pd.DataFrame([summary]).to_csv(outdir / "summary.csv", index=False)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    np.savetxt(outdir / "hessian_matrix.csv", H, delimiter=",")
    np.savetxt(outdir / "hessian_eigenvalues.csv", h_vals, delimiter=",")
    np.savetxt(outdir / "laplacian_eigenvalues.csv", lap_vals, delimiter=",")

    plt.figure(figsize=(7, 4))
    plt.plot(np.arange(len(h_vals)), h_vals, marker=".", linewidth=1)
    plt.axhline(0, linestyle="--")
    plt.xlabel("mode index")
    plt.ylabel(r"Hessian eigenvalue $\omega^2$")
    plt.title("BuP Hessian spectrum")
    plt.tight_layout()
    plt.savefig(outdir / "fig_eigenvalues.png", dpi=170)
    plt.close()

    plt.figure(figsize=(6, 4.5))
    plt.scatter(low["k2_eff_comb_laplacian"], low["omega2_hessian"], label="low modes")
    if np.isfinite(slope):
        x = np.linspace(low["k2_eff_comb_laplacian"].min(), low["k2_eff_comb_laplacian"].max(), 100)
        plt.plot(x, slope*x + intercept, label=f"fit: c^2={slope:.3g}, R2={r2_disp:.3f}")
    plt.xlabel(r"graph $k^2$ from Laplacian energy")
    plt.ylabel(r"Hessian $\omega^2$")
    plt.title("Effective dispersion relation")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_dispersion.png", dpi=170)
    plt.close()

    md = [
        "# Paper 21 — BuP Hessian modes v5",
        "",
        "The effective BuP action is linearized around a background entanglement graph.",
        "",
        r"\[H_{ab}=\left.\frac{\partial^2 S_{\rm BuP}}{\partial\phi_a\partial\phi_b}\right|_{\phi=0}.\]",
        "",
        r"\[H\phi_n=\omega_n^2\phi_n.\]",
        "",
        "## Summary",
        "",
    ]
    for k, v in summary.items():
        md.append(f"- `{k}`: `{v}`")
    (outdir / "summary.md").write_text("\n".join(md), encoding="utf-8")

    print("=" * 100)
    print("Paper 21 — BuP gravitational-wave prototype v5 Hessian modes")
    print("=" * 100)
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nFiles written:")
    for name in [
        "summary.csv",
        "summary.json",
        "summary.md",
        "eigenmodes.csv",
        "dispersion.csv",
        "hessian_matrix.csv",
        "hessian_eigenvalues.csv",
        "laplacian_eigenvalues.csv",
        "fig_eigenvalues.png",
        "fig_dispersion.png",
        "fig_mode_*.png",
    ]:
        print(outdir / name)


if __name__ == "__main__":
    main()

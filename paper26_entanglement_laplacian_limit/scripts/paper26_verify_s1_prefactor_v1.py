#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — S1 prefactor verification v1

Checks whether the empirical prefactor

    c_empirical * N * epsilon^(3/2)

converges to sqrt(pi), as predicted by the continuum expansion of the
unnormalized Gaussian graph Laplacian on the unit circle S^1.

For the convention

    W_ij = exp[-d(theta_i, theta_j)^2 / (4 epsilon)]
    L = D - W

and uniform sampling on S^1 of length 2pi, the continuum expansion gives

    L f(theta) ~ -(N / sqrt(pi)) epsilon^(3/2) f''(theta).

Therefore

    c_{N,epsilon} = sqrt(pi) / (N epsilon^(3/2)).
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path("paper26_entanglement_laplacian_limit")
RESULTS = ROOT / "results" / "paper26_laplian_limit_v1"
# fallback for the actual folder name
if not RESULTS.exists():
    RESULTS = ROOT / "results" / "paper26_laplacian_limit_v1"

FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

norm_csv = RESULTS / "normalization_scan.csv"
df = pd.read_csv(norm_csv)

sqrt_pi = float(np.sqrt(np.pi))

df["theory_prefactor_sqrt_pi"] = sqrt_pi
df["prefactor_minus_sqrt_pi"] = (
    df["prefactor_estimate_c_times_N_eps_3_over_2"] - sqrt_pi
)
df["prefactor_rel_error_to_sqrt_pi"] = (
    df["prefactor_minus_sqrt_pi"].abs() / sqrt_pi
)

out_csv = RESULTS / "prefactor_theory_s1.csv"
df.to_csv(out_csv, index=False)

best_idx = df["prefactor_rel_error_to_sqrt_pi"].idxmin()
best = df.loc[best_idx].to_dict()

summary = {
    "experiment": "paper26_s1_prefactor_verification_v1",
    "theory": "c_{N,epsilon} = sqrt(pi)/(N epsilon^(3/2))",
    "theory_prefactor": sqrt_pi,
    "best_by_prefactor_error": {
        "N": int(best["N"]),
        "epsilon": float(best["epsilon"]),
        "empirical_prefactor": float(best["prefactor_estimate_c_times_N_eps_3_over_2"]),
        "theory_prefactor": sqrt_pi,
        "relative_error_to_sqrt_pi": float(best["prefactor_rel_error_to_sqrt_pi"]),
        "mean_rel_spectral_error": float(best["mean_rel_error"]),
        "max_rel_spectral_error": float(best["max_rel_error"]),
    },
    "mean_prefactor_rel_error": float(df["prefactor_rel_error_to_sqrt_pi"].mean()),
    "max_prefactor_rel_error": float(df["prefactor_rel_error_to_sqrt_pi"].max()),
    "output_csv": str(out_csv),
}

out_json = RESULTS / "prefactor_theory_s1_summary.json"
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

# Figure
fig, ax = plt.subplots(figsize=(10, 6))

for N, group in df.groupby("N"):
    group = group.sort_values("epsilon")
    ax.plot(
        group["epsilon"],
        group["prefactor_estimate_c_times_N_eps_3_over_2"],
        marker="o",
        label=f"N={int(N)}",
    )

ax.axhline(sqrt_pi, linestyle="--", linewidth=1.5, label=r"Theory $\sqrt{\pi}$")

ax.set_xscale("log")
ax.set_xlabel(r"$\epsilon$")
ax.set_ylabel(r"$c_{N,\epsilon}N\epsilon^{3/2}$")
ax.set_title(r"S1 prefactor verification: $c_{N,\epsilon}N\epsilon^{3/2}\to\sqrt{\pi}$")
ax.grid(True, alpha=0.3)
ax.legend()

png = FIGURES / "fig04_s1_prefactor_to_sqrt_pi.png"
pdf = FIGURES / "fig04_s1_prefactor_to_sqrt_pi.pdf"

plt.savefig(png, dpi=300, bbox_inches="tight")
plt.savefig(pdf, bbox_inches="tight")

print(f"[OK] wrote {out_csv}")
print(f"[OK] wrote {out_json}")
print(f"[OK] wrote {png}")
print(f"[OK] wrote {pdf}")
print()
print("Theory prefactor sqrt(pi):", sqrt_pi)
print("Best empirical prefactor:", summary["best_by_prefactor_error"])

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 21 — Make figures

Generates:
  fig1_pipeline_paper9_to_slacs.png
  fig2_transition_window_summary.png
  fig3_dynamic_shuffle_control.png
  fig4_sersic_measured_vs_n4.png
  fig5_fixed_point_Cobs_zero.png
  fig6_alpha_eff_diffusion_window.png
  fig7_triple_convergence.png
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


BASE = "papers/paper21_slacs_fixed_point"
FIG = f"{BASE}/figures"
os.makedirs(FIG, exist_ok=True)


def savefig(name):
    path = f"{FIG}/{name}"
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()
    print("Wrote:", path)


# ------------------------------------------------------------
# Fig 1 — Pipeline schematic
# ------------------------------------------------------------

plt.figure(figsize=(10, 3))

steps = [
    r"$S_{\rm flux}$",
    r"$L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}$",
    r"$\Phi_{\rm BuP}$",
    r"$C_{\rm obs}$",
    r"SLACS fixed point"
]

x = np.arange(len(steps))
y = np.zeros(len(steps))

plt.scatter(x, y, s=600)

for i, label in enumerate(steps):
    plt.text(i, 0.08, label, ha="center", va="bottom", fontsize=11)
    if i < len(steps) - 1:
        plt.arrow(i + 0.18, 0, 0.55, 0, length_includes_head=True,
                  head_width=0.035, head_length=0.08)

plt.ylim(-0.2, 0.35)
plt.axis("off")
plt.title("Paper 9 prediction tested on SLACS")
savefig("fig1_pipeline_paper9_to_slacs.png")


# ------------------------------------------------------------
# Fig 2 — Transition window summary
# ------------------------------------------------------------

summary_path = f"{BASE}/results/compare_sersic_measured_transition_v1/compare_sersic_transition_top_windows.csv"

if os.path.exists(summary_path):
    df = pd.read_csv(summary_path)

    # Prefer measured_n_only
    g = df[df["tag"] == "measured_n_only"].copy()
    if len(g):
        g = g.sort_values("mean_delta_logM_minus_Phi", ascending=False).head(20)

        plt.figure(figsize=(8, 5))
        plt.scatter(g["center"], g["mean_delta_logM_minus_Phi"], s=70)
        plt.axvspan(11.545, 11.645, alpha=0.2)
        plt.axhline(0, linestyle="--")
        plt.xlabel(r"$\log M_\star$")
        plt.ylabel(r"$\langle |\epsilon_{\log M}|-|\epsilon_{\Phi}| \rangle$")
        plt.title(r"SLACS transition window: $\Phi_{\rm BuP}$ vs $\log M_\star$")
        savefig("fig2_transition_window_summary.png")


# ------------------------------------------------------------
# Fig 3 — Dynamic shuffle control
# ------------------------------------------------------------

shuffle_path = f"{BASE}/results/dynamic_shuffle_sqrtM_v1/dynamic_shuffle_sqrtM_by_mode.csv"

if os.path.exists(shuffle_path):
    df = pd.read_csv(shuffle_path)

    plt.figure(figsize=(8, 5))
    plt.bar(df["mode"].astype(str), df["max"])
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Max LOO improvement (%)")
    plt.title("Dynamic shuffle control")
    savefig("fig3_dynamic_shuffle_control.png")


# ------------------------------------------------------------
# Fig 4 — Measured Sérsic vs forced n=4
# ------------------------------------------------------------

sersic_summary = f"{BASE}/results/compare_sersic_measured_transition_v1/compare_sersic_transition_summary.csv"

if os.path.exists(sersic_summary):
    df = pd.read_csv(sersic_summary)

    keep = df[df["tag"].isin(["measured_n_only", "matched_forced_n4", "measured_plus_fallback"])].copy()

    if len(keep):
        plt.figure(figsize=(8, 5))
        plt.bar(keep["tag"], keep["Phi_LOO_improvement"])
        plt.xticks(rotation=25, ha="right")
        plt.ylabel(r"$\Phi_{\rm BuP}$ LOO improvement (%)")
        plt.title("Effect of measured Sérsic indices")
        savefig("fig4_sersic_measured_vs_n4.png")


# ------------------------------------------------------------
# Fig 5 — Fixed point Cobs = 0
# ------------------------------------------------------------

data_path = f"{BASE}/data/slacs_v5_measured_n_only.csv"

if os.path.exists(data_path):
    df = pd.read_csv(data_path)

    if {"logMs", "theta_E_obs_arcsec", "theta_E_baryon_arcsec"}.issubset(df.columns):
        df = df.dropna(subset=["logMs", "theta_E_obs_arcsec", "theta_E_baryon_arcsec"]).copy()
        df = df[(df["theta_E_obs_arcsec"] > 0) & (df["theta_E_baryon_arcsec"] > 0)]
        df["C_obs"] = np.log(df["theta_E_obs_arcsec"] / df["theta_E_baryon_arcsec"])

        x = df["logMs"].values
        y = df["C_obs"].values

        coeff = np.polyfit(x, y, deg=1)
        xx = np.linspace(np.min(x), np.max(x), 200)
        yy = np.polyval(coeff, xx)

        root = -coeff[1] / coeff[0] if abs(coeff[0]) > 1e-12 else np.nan

        plt.figure(figsize=(8, 5))
        plt.scatter(x, y, s=45)
        plt.plot(xx, yy)
        plt.axhline(0, linestyle="--")
        if np.isfinite(root):
            plt.axvline(root, linestyle="--")
            plt.text(root + 0.01, np.nanmin(y), rf"$C_{{obs}}=0$: {root:.3f}", rotation=90, va="bottom")
        plt.axvspan(11.545, 11.645, alpha=0.2)
        plt.xlabel(r"$\log M_\star$")
        plt.ylabel(r"$C_{\rm obs}$")
        plt.title(r"Observed fixed point: $C_{\rm obs}=0$")
        savefig("fig5_fixed_point_Cobs_zero.png")


# ------------------------------------------------------------
# Fig 6 — Alpha_eff diffusion-window scan
# ------------------------------------------------------------

scan_path = f"{BASE}/results/dw_window_alpha_scan_v1/measured_n_only_scan_window_results.csv"

if os.path.exists(scan_path):
    df = pd.read_csv(scan_path)

    plt.figure(figsize=(8, 5))
    sc = plt.scatter(df["t_high"], df["mean_alpha"], c=df["t_low"], s=35)
    plt.axhline(1.0, linestyle="--")
    plt.xlabel(r"$t_{\max}$")
    plt.ylabel(r"$\langle \alpha_{\rm eff} \rangle$")
    plt.title(r"Diffusion-window scan: $\alpha_{\rm eff}=1$")
    plt.colorbar(sc, label=r"$t_{\min}$")
    savefig("fig6_alpha_eff_diffusion_window.png")


# ------------------------------------------------------------
# Fig 7 — Triple convergence
# ------------------------------------------------------------

plt.figure(figsize=(8, 5))

labels = [
    r"$C_{\rm obs}=0$",
    r"$\Phi_{\rm BuP}>\log M_\star$",
    r"$\alpha_{\rm eff}\simeq1$"
]

x = [11.60, 11.595, 11.595]
y = [3, 2, 1]

plt.scatter(x, y, s=140)

for xi, yi, lab in zip(x, y, labels):
    plt.text(xi + 0.006, yi, lab, va="center", fontsize=11)

plt.axvspan(11.545, 11.645, alpha=0.2)
plt.xlabel(r"$\log M_\star$")
plt.yticks([])
plt.title("SLACS fixed point: triple convergence")
plt.xlim(11.45, 11.75)
plt.ylim(0.5, 3.5)
savefig("fig7_triple_convergence.png")


print("Done. Figures directory:", FIG)

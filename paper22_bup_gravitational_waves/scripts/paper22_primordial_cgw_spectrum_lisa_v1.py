#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 22 — BuP primordial gravitational-wave spectrum prototype v6.

Exploratory phenomenological simulation of a reduced primordial graph speed
c_GW(z)=c alpha(z), with alpha(z)<1 in a strongly entangled early phase.

This is not a full Boltzmann/tensor-mode solver. It is a controlled prototype
to test whether alpha(z)<1 can leave a PTA/LISA-band spectral imprint.

Outputs:
  summary.json
  spectrum.csv
  band_summary.csv
  fig_spectrum_lisa_pta.png
  fig_alpha_of_z.png
  fig_transfer.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def alpha_of_z(z, alpha_early, z_transition, sharpness):
    z = np.asarray(z, dtype=float)
    return alpha_early + (1.0 - alpha_early) / (1.0 + (z / z_transition) ** sharpness)


def z_eff_of_f(f, f_min, z_min, z_max, z_power):
    f = np.asarray(f, dtype=float)
    z = z_min * (f / f_min) ** z_power
    return np.clip(z, z_min, z_max)


def omega_baseline(f, omega_ref, f_ref, n_t, running):
    x = f / f_ref
    return omega_ref * np.exp(n_t * np.log(x) + 0.5 * running * np.log(x) ** 2)


def lisa_toy_sensitivity(f):
    # Rough visual guide only, not an official LISA sensitivity curve.
    f0 = 3e-3
    floor = 1.0e-13
    low_wall = (f0 / f) ** 4
    high_wall = (f / f0) ** 2
    return floor * (1.0 + 0.08 * low_wall + 0.15 * high_wall)


def transfer_alpha(alpha, p_accum, tau_damp):
    return alpha ** (-p_accum) * np.exp(-tau_damp * (1.0 / alpha - 1.0))


def transfer_mass(f, f_mass, mass_power):
    if f_mass <= 0:
        return np.ones_like(f)
    return 1.0 / (1.0 + (f_mass / f) ** mass_power)


def summarize_band(df, name, f_lo, f_hi):
    sub = df[(df["f_Hz"] >= f_lo) & (df["f_Hz"] <= f_hi)]
    if len(sub) == 0:
        return {
            "band": name,
            "f_lo": f_lo,
            "f_hi": f_hi,
            "n": 0,
            "Omega_GR_median": np.nan,
            "Omega_BuP_median": np.nan,
            "ratio_median": np.nan,
            "ratio_min": np.nan,
            "ratio_max": np.nan,
            "alpha_median": np.nan,
            "z_eff_median": np.nan,
        }
    ratio = sub["Omega_BuP"] / sub["Omega_GR"]
    return {
        "band": name,
        "f_lo": f_lo,
        "f_hi": f_hi,
        "n": int(len(sub)),
        "Omega_GR_median": float(np.median(sub["Omega_GR"])),
        "Omega_BuP_median": float(np.median(sub["Omega_BuP"])),
        "ratio_median": float(np.median(ratio)),
        "ratio_min": float(np.min(ratio)),
        "ratio_max": float(np.max(ratio)),
        "alpha_median": float(np.median(sub["alpha_eff"])),
        "z_eff_median": float(np.median(sub["z_eff"])),
    }


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True)

    p.add_argument("--f-min", type=float, default=1e-10)
    p.add_argument("--f-max", type=float, default=1e1)
    p.add_argument("--n-f", type=int, default=1400)

    p.add_argument("--omega-ref", type=float, default=1e-12)
    p.add_argument("--f-ref", type=float, default=1e-3)
    p.add_argument("--n-t", type=float, default=0.0)
    p.add_argument("--running", type=float, default=0.0)

    p.add_argument("--alpha-early", type=float, default=0.35)
    p.add_argument("--z-transition", type=float, default=1e13)
    p.add_argument("--sharpness", type=float, default=3.0)

    p.add_argument("--z-min", type=float, default=1.0)
    p.add_argument("--z-max", type=float, default=1e15)
    p.add_argument("--z-power", type=float, default=2.0)

    p.add_argument("--p-accum", type=float, default=1.0)
    p.add_argument("--tau-damp", type=float, default=0.35)
    p.add_argument("--f-mass", type=float, default=1e-6)
    p.add_argument("--mass-power", type=float, default=2.0)
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    f = np.logspace(np.log10(args.f_min), np.log10(args.f_max), args.n_f)
    z_eff = z_eff_of_f(f, args.f_min, args.z_min, args.z_max, args.z_power)
    alpha = alpha_of_z(z_eff, args.alpha_early, args.z_transition, args.sharpness)

    omega_gr = omega_baseline(f, args.omega_ref, args.f_ref, args.n_t, args.running)
    t_alpha = transfer_alpha(alpha, args.p_accum, args.tau_damp)
    t_mass = transfer_mass(f, args.f_mass, args.mass_power)
    omega_bup = omega_gr * t_alpha * t_mass
    lisa_guide = lisa_toy_sensitivity(f)

    df = pd.DataFrame({
        "f_Hz": f,
        "z_eff": z_eff,
        "alpha_eff": alpha,
        "T_alpha": t_alpha,
        "T_mass": t_mass,
        "T_total": t_alpha * t_mass,
        "Omega_GR": omega_gr,
        "Omega_BuP": omega_bup,
        "Omega_BuP_over_GR": omega_bup / omega_gr,
        "LISA_toy_Omega_sensitivity": lisa_guide,
    })
    df.to_csv(outdir / "spectrum.csv", index=False)

    bands = [
        ("PTA_nHz", 1e-9, 1e-7),
        ("LISA_mHz", 1e-4, 1e-1),
        ("ground_10_100Hz", 10.0, 100.0),
    ]
    band_df = pd.DataFrame([summarize_band(df, *b) for b in bands])
    band_df.to_csv(outdir / "band_summary.csv", index=False)

    lisa_row = band_df[band_df["band"] == "LISA_mHz"].iloc[0].to_dict()
    pta_row = band_df[band_df["band"] == "PTA_nHz"].iloc[0].to_dict()

    summary = {
        "status": "exploratory_phenomenological_model",
        "alpha_early": args.alpha_early,
        "z_transition": args.z_transition,
        "sharpness": args.sharpness,
        "p_accum": args.p_accum,
        "tau_damp": args.tau_damp,
        "f_mass": args.f_mass,
        "mass_power": args.mass_power,
        "omega_ref": args.omega_ref,
        "f_ref": args.f_ref,
        "n_t": args.n_t,
        "running": args.running,
        "PTA_ratio_median": float(pta_row["ratio_median"]),
        "LISA_ratio_median": float(lisa_row["ratio_median"]),
        "PTA_alpha_median": float(pta_row["alpha_median"]),
        "LISA_alpha_median": float(lisa_row["alpha_median"]),
        "interpretation": (
            "A reduced primordial graph speed changes the observed SGWB through "
            "accumulation, damping and a possible geometric mass/gap filter. "
            "A LISA-band curvature or enhancement/suppression relative to GR is "
            "the target observable."
        ),
    }
    with open(outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    plt.figure(figsize=(8, 5))
    plt.loglog(df["f_Hz"], df["Omega_GR"], label="GR baseline")
    plt.loglog(df["f_Hz"], df["Omega_BuP"], label="BuP variable-speed model")
    plt.loglog(df["f_Hz"], df["LISA_toy_Omega_sensitivity"], linestyle="--", label="LISA toy guide")
    plt.axvspan(1e-9, 1e-7, alpha=0.12, label="PTA band")
    plt.axvspan(1e-4, 1e-1, alpha=0.10, label="LISA band")
    plt.xlabel("frequency today f [Hz]")
    plt.ylabel("Omega_GW(f)")
    plt.title("BuP primordial SGWB: reduced early graph speed")
    y_min = max(1e-22, min(df["Omega_BuP"].min(), df["Omega_GR"].min()) * 0.5)
    y_max = max(df["Omega_BuP"].max(), df["Omega_GR"].max(), 1e-13) * 5
    plt.ylim(y_min, y_max)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "fig_spectrum_lisa_pta.png", dpi=180)
    plt.close()

    z_plot = np.logspace(0, 15, 1000)
    a_plot = alpha_of_z(z_plot, args.alpha_early, args.z_transition, args.sharpness)
    plt.figure(figsize=(7, 4.5))
    plt.semilogx(z_plot, a_plot)
    plt.axvline(args.z_transition, linestyle="--", label="z_transition")
    plt.xlabel("redshift z")
    plt.ylabel("alpha(z)=c_GW(z)/c")
    plt.title("BuP graph-speed transition")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "fig_alpha_of_z.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.semilogx(df["f_Hz"], df["T_alpha"], label="T_alpha")
    plt.semilogx(df["f_Hz"], df["T_mass"], label="T_mass")
    plt.semilogx(df["f_Hz"], df["T_total"], label="T_total")
    plt.axvspan(1e-9, 1e-7, alpha=0.12, label="PTA")
    plt.axvspan(1e-4, 1e-1, alpha=0.10, label="LISA")
    plt.xlabel("frequency today f [Hz]")
    plt.ylabel("transfer factor")
    plt.title("BuP transfer functions")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / "fig_transfer.png", dpi=180)
    plt.close()

    print("=" * 100)
    print("Paper 22 — v6 primordial SGWB with reduced early graph speed")
    print("=" * 100)
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nBand summary:")
    print(band_df.to_string(index=False))

    print("\nFiles written:")
    for name in [
        "summary.json",
        "spectrum.csv",
        "band_summary.csv",
        "fig_spectrum_lisa_pta.png",
        "fig_alpha_of_z.png",
        "fig_transfer.png",
    ]:
        print(outdir / name)


if __name__ == "__main__":
    main()

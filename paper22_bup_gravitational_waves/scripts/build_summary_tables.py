#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build summary tables for Paper 22 — BuP gravitational waves.

This script collects key JSON/CSV outputs from v4, v5 and v6 and writes
publication-ready summary CSV tables into:

papers/paper22_bup_gravitational_waves/tables/
"""

from pathlib import Path
import json
import glob
import pandas as pd


ROOT = Path("papers/paper22_bup_gravitational_waves")
RESULTS = ROOT / "results"
TABLES = ROOT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def build_v4_speed_table():
    rows = []

    candidates = [
        ("plus", RESULTS / "v4_dynamic_source_plus" / "summary.json"),
        ("cross", RESULTS / "v4_dynamic_source_cross" / "summary.json"),
    ]

    for pol, path in candidates:
        if not path.exists():
            print(f"[skip] missing {path}")
            continue

        d = load_json(path)

        rows.append({
            "polarization": pol,
            "speed_from_arrival": d.get("speed_from_arrival"),
            "speed_from_arrival_over_input": d.get("speed_from_arrival_over_input"),
            "arrival_fit_r2": d.get("arrival_fit_r2"),
            "speed_from_peak": d.get("speed_from_peak"),
            "speed_from_peak_over_input": d.get("speed_from_peak_over_input"),
            "peak_fit_r2": d.get("peak_fit_r2"),
            "mean_quadrupole_fit_r2": d.get("mean_quadrupole_fit_r2"),
            "median_quadrupole_fit_r2": d.get("median_quadrupole_fit_r2"),
        })

    if rows:
        df = pd.DataFrame(rows)
        out = TABLES / "table_v4_speed.csv"
        df.to_csv(out, index=False)
        print(f"written: {out}")
    else:
        print("[skip] no v4 rows")


def build_v5_eta_table():
    rows = []

    pattern = str(RESULTS / "v5_eta_scan" / "*" / "summary.json")
    for path in glob.glob(pattern):
        d = load_json(path)
        rows.append({
            "eta_smooth": d.get("eta_smooth"),
            "negative_modes": d.get("hessian_negative_modes"),
            "dispersion_r2": d.get("dispersion_r2"),
            "c_graph": d.get("c_graph_from_hessian"),
            "c2": d.get("dispersion_slope_c2"),
            "mass2": d.get("dispersion_intercept_mass2"),
        })

    if rows:
        df = pd.DataFrame(rows).sort_values("eta_smooth")
        out = TABLES / "table_v5_eta_scan.csv"
        df.to_csv(out, index=False)
        print(f"written: {out}")
    else:
        print("[skip] no v5 eta rows")


def build_v6_alpha_table():
    rows = []

    pattern = str(RESULTS / "v6_alpha_scan" / "*" / "summary.json")
    for path in glob.glob(pattern):
        d = load_json(path)
        rows.append({
            "alpha_early": d.get("alpha_early"),
            "PTA_ratio": d.get("PTA_ratio_median"),
            "LISA_ratio": d.get("LISA_ratio_median"),
            "PTA_alpha": d.get("PTA_alpha_median"),
            "LISA_alpha": d.get("LISA_alpha_median"),
            "f_mass": d.get("f_mass"),
            "tau_damp": d.get("tau_damp"),
        })

    if rows:
        df = pd.DataFrame(rows).sort_values("alpha_early")
        out = TABLES / "table_v6_alpha_scan.csv"
        df.to_csv(out, index=False)
        print(f"written: {out}")
    else:
        print("[skip] no v6 alpha rows")


def build_v6_mass_gap_table():
    rows = []

    candidates = [
        ("0", RESULTS / "v6_primordial_cgw_lisa_no_mass" / "summary.json"),
        ("1e-11", RESULTS / "v6_primordial_cgw_lisa_mass_1e-11" / "summary.json"),
        ("1e-6", RESULTS / "v6_primordial_cgw_lisa" / "summary.json"),
    ]

    for label, path in candidates:
        if not path.exists():
            print(f"[skip] missing {path}")
            continue

        d = load_json(path)
        rows.append({
            "f_mass_Hz": d.get("f_mass"),
            "PTA_ratio": d.get("PTA_ratio_median"),
            "LISA_ratio": d.get("LISA_ratio_median"),
            "PTA_alpha": d.get("PTA_alpha_median"),
            "LISA_alpha": d.get("LISA_alpha_median"),
            "verdict": (
                "viable_PTA_unchanged_LISA_enhanced"
                if d.get("PTA_ratio_median", 0) > 0.9
                else "not_viable_PTA_suppressed"
            ),
        })

    if rows:
        df = pd.DataFrame(rows)
        out = TABLES / "table_v6_mass_gap_comparison.csv"
        df.to_csv(out, index=False)
        print(f"written: {out}")
    else:
        print("[skip] no v6 mass-gap rows")


def main():
    build_v4_speed_table()
    build_v5_eta_table()
    build_v6_alpha_table()
    build_v6_mass_gap_table()


if __name__ == "__main__":
    main()

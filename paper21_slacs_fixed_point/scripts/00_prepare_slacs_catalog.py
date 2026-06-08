#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 21 — Prepare SLACS catalogues

Creates:
  data/slacs_v5_bup_final.csv
  data/slacs_v5_measured_n_only.csv
  data/slacs_v5_matched_forced_n4.csv
  data/slacs_v5_measured_plus_fallback.csv

Input expected:
  original SLACS BuP CSV from Paper 13 / SLACS Hbup direct experiment.

This script standardizes the Sérsic-index catalogues used in Paper 21.
"""

import argparse
import os
import numpy as np
import pandas as pd


def read_csv_auto(path):
    return pd.read_csv(path, sep=None, engine="python")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--input-csv",
        default="papers/paper13_matiere_intrication_critique/experiments/slacs_Hbup_direct/data/slacs_v5_bup_final.csv",
        help="Original SLACS BuP CSV"
    )

    ap.add_argument(
        "--output-dir",
        default="papers/paper21_slacs_fixed_point/data",
        help="Paper 21 data directory"
    )

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = read_csv_auto(args.input_csv)

    # Save original copy
    out_base = os.path.join(args.output_dir, "slacs_v5_bup_final.csv")
    df.to_csv(out_base, index=False)

    print("Input:", args.input_csv)
    print("N total:", len(df))
    print("Columns:")
    print(df.columns.tolist())

    # Find measured Sérsic column
    candidates = [
        "sersic_n_measured",
        "n_sersic",
        "sersic_n",
        "n",
        "sersic_n_final",
    ]

    ncol = None

    for c in candidates:
        if c in df.columns:
            x = pd.to_numeric(df[c], errors="coerce")
            valid = x.notna() & np.isfinite(x) & (x > 0.3) & (x < 12.0)
            print(f"{c}: valid {valid.sum()}/{len(df)}")
            if valid.sum() > 0 and ncol is None:
                ncol = c

    if ncol is None:
        raise RuntimeError("No usable measured Sérsic-index column found.")

    print("\nUsing Sérsic column:", ncol)

    df["_n_measured_tmp"] = pd.to_numeric(df[ncol], errors="coerce")
    valid = (
        df["_n_measured_tmp"].notna()
        & np.isfinite(df["_n_measured_tmp"])
        & (df["_n_measured_tmp"] > 0.3)
        & (df["_n_measured_tmp"] < 12.0)
    )

    # 1. Measured n only
    measured = df[valid].copy()
    measured["sersic_n_measured"] = measured["_n_measured_tmp"]
    measured["sersic_n_final"] = measured["_n_measured_tmp"]
    measured["n_source"] = "measured"

    # 2. Same galaxies, forced n=4
    forced = df[valid].copy()
    forced["sersic_n_measured"] = forced["_n_measured_tmp"]
    forced["sersic_n_final"] = 4.0
    forced["n_source"] = "forced_4_matched"

    # 3. Full catalogue: measured n if available, fallback n=4
    fallback = df.copy()
    fallback["sersic_n_measured"] = fallback["_n_measured_tmp"]
    fallback["sersic_n_final"] = np.where(valid, fallback["_n_measured_tmp"], 4.0)
    fallback["n_source"] = np.where(valid, "measured", "default_4")

    for d in (measured, forced, fallback):
        d.drop(columns=["_n_measured_tmp"], inplace=True, errors="ignore")

    out_measured = os.path.join(args.output_dir, "slacs_v5_measured_n_only.csv")
    out_forced = os.path.join(args.output_dir, "slacs_v5_matched_forced_n4.csv")
    out_fallback = os.path.join(args.output_dir, "slacs_v5_measured_plus_fallback.csv")

    measured.to_csv(out_measured, index=False)
    forced.to_csv(out_forced, index=False)
    fallback.to_csv(out_fallback, index=False)

    print("\nWrote:")
    print(out_base, len(df))
    print(out_measured, len(measured))
    print(out_forced, len(forced))
    print(out_fallback, len(fallback))

    print("\nMeasured n summary:")
    print(measured["sersic_n_final"].describe())

    preview_cols = [c for c in ["name", "name_key", "logMs", "Re_arcsec", "sersic_n_final", "n_source"] if c in measured.columns]
    print("\nPreview:")
    print(measured[preview_cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()


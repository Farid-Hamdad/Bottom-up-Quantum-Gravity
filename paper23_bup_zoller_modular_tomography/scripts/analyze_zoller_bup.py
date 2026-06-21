#!/usr/bin/env python3
"""Minimal BuP tomography diagnostics for the Zoller/Joshi 2023 data files.

The script intentionally avoids scipy so it can run in this workspace with only
NumPy/Pandas available. It implements the small subset of MATLAB v5 parsing
needed for the uploaded analysed-data files.
"""

from __future__ import annotations

import csv
import json
import math
import struct
import zlib
from pathlib import Path

import numpy as np


MI_COMPRESSED = 15
MI_MATRIX = 14
MX_CELL_CLASS = 1
MX_CHAR_CLASS = 4
MX_DOUBLE_CLASS = 6

DTYPES = {
    1: "i1",
    2: "u1",
    3: "<i2",
    4: "<u2",
    5: "<i4",
    6: "<u4",
    7: "<f4",
    9: "<f8",
    12: "<i8",
    13: "<u8",
}


def read_tag(buf: bytes, pos: int, align: bool = True):
    raw = struct.unpack_from("<I", buf, pos)[0]
    small_n = raw >> 16
    small_type = raw & 0xFFFF
    if small_n:
        return small_type, small_n, buf[pos + 4 : pos + 4 + small_n], pos + 8

    dtype, nbytes = struct.unpack_from("<II", buf, pos)
    data_start = pos + 8
    next_pos = data_start + nbytes
    if align:
        next_pos += (-nbytes) % 8
    return dtype, nbytes, buf[data_start : data_start + nbytes], next_pos


def parse_matrix(data: bytes):
    pos = 0
    _, _, flags, pos = read_tag(data, pos)
    matlab_class = struct.unpack_from("<I", flags + b"\0" * 4, 0)[0] & 0xFF

    _, _, dims_raw, pos = read_tag(data, pos)
    dims = np.frombuffer(dims_raw, dtype="<i4").tolist()

    _, _, name_raw, pos = read_tag(data, pos)
    name = name_raw.decode("latin1") if name_raw else ""

    if matlab_class == MX_DOUBLE_CLASS:
        dtype, _, payload, _ = read_tag(data, pos)
        arr = np.frombuffer(payload, dtype=np.dtype(DTYPES[dtype])).copy()
        if dims:
            arr = arr.reshape(tuple(dims), order="F")
        return name, arr

    if matlab_class == MX_CHAR_CLASS:
        dtype, _, payload, _ = read_tag(data, pos)
        if dtype in DTYPES:
            arr = np.frombuffer(payload, dtype=np.dtype(DTYPES[dtype])).copy()
            text = "".join(chr(int(x)) for x in arr.ravel(order="F") if int(x) != 0)
        else:
            text = payload.decode("latin1", errors="ignore")
        return name, text

    if matlab_class == MX_CELL_CLASS:
        cells = []
        while pos + 8 <= len(data):
            dtype, _, payload, pos = read_tag(data, pos)
            if dtype == MI_MATRIX:
                if payload:
                    cells.append(parse_matrix(payload))
                else:
                    cells.append(("", np.array([])))
        return name, {"type": "cell", "dims": dims, "cells": cells}

    return name, {"type": f"class_{matlab_class}", "dims": dims}


def load_mat_v5(path: Path) -> dict:
    buf = path.read_bytes()
    out = {}
    pos = 128
    while pos + 8 <= len(buf):
        dtype, nbytes, payload, next_pos = read_tag(buf, pos, align=False)
        if dtype == MI_COMPRESSED:
            decompressed = zlib.decompress(payload)
            sub_pos = 0
            while sub_pos + 8 <= len(decompressed):
                sub_dtype, _, sub_payload, sub_pos = read_tag(decompressed, sub_pos)
                if sub_dtype == MI_MATRIX:
                    name, value = parse_matrix(sub_payload)
                    out[name] = value
        elif dtype == MI_MATRIX:
            name, value = parse_matrix(payload)
            out[name] = value
        pos = next_pos
    return out


def fit_r2(y: np.ndarray, xcols: list[np.ndarray]) -> float:
    x = np.vstack(xcols).T
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def profile_metrics(beta: np.ndarray) -> dict:
    positive = beta[beta > 1e-9]
    n = len(beta)
    interior_idx = np.where(beta > 1e-9)[0]
    if positive.size < 2:
        return {}

    first = positive[0]
    last = positive[-1]
    edge_mean = float((first + last) / 2.0)
    center = float(positive[len(positive) // 2])
    beta_edge_to_center = edge_mean / center if center else float("nan")
    rho = 1.0 / positive
    rho /= rho.sum()
    rho_edge_fraction = float(rho[0] + rho[-1]) if rho.size >= 2 else float("nan")
    rho_center_fraction = float(rho[len(rho) // 2])
    rho_edge_to_center = float(((rho[0] + rho[-1]) / 2.0) / rho_center_fraction)

    local_pos = interior_idx.astype(float)
    local_pos -= local_pos.min()
    span = max(float(local_pos.max()), 1.0)
    u = local_pos / span
    parabola = u * (1.0 - u)
    triangle = np.minimum(u, 1.0 - u)

    return {
        "n_coeff": int(n),
        "n_positive": int(positive.size),
        "beta_min_positive": float(np.min(positive)),
        "beta_max": float(np.max(positive)),
        "beta_mean": float(np.mean(positive)),
        "beta_cv": float(np.std(positive) / np.mean(positive)),
        "beta_edge_to_center": float(beta_edge_to_center),
        "rho_edge_fraction": rho_edge_fraction,
        "rho_center_fraction": rho_center_fraction,
        "rho_edge_to_center": rho_edge_to_center,
        "r2_flat": fit_r2(positive, [np.ones_like(positive)]),
        "r2_parabola": fit_r2(positive, [np.ones_like(positive), parabola]),
        "r2_triangle": fit_r2(positive, [np.ones_like(positive), triangle]),
    }


def extract_profiles(label: str, mat: dict) -> list[dict]:
    rows = []
    for group_name, group in mat.items():
        if not isinstance(group, dict) or group.get("type") != "cell":
            continue
        for idx, (_, arr) in enumerate(group["cells"]):
            if not isinstance(arr, np.ndarray) or arr.ndim != 2 or arr.shape[1] < 2:
                continue
            beta = arr[:, 1].astype(float)
            metrics = profile_metrics(beta)
            if not metrics:
                continue
            rows.append(
                {
                    "state": label,
                    "group": group_name,
                    "cell_index": idx,
                    "n_rows": int(arr.shape[0]),
                    "n_cols": int(arr.shape[1]),
                    "kind": "experimental" if arr.shape[1] >= 4 else "mps_or_fit",
                    **metrics,
                }
            )
    return rows


def write_svg_profile(path: Path, title: str, profiles: list[tuple[str, np.ndarray]]):
    width, height = 760, 420
    margin = 55
    all_y = np.concatenate([p[1][p[1] > 1e-9] for p in profiles])
    ymax = float(all_y.max()) * 1.05

    def sx(i, n):
        return margin + (width - 2 * margin) * i / max(n - 1, 1)

    def sy(y):
        return height - margin - (height - 2 * margin) * y / ymax

    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="Arial" font-size="18">{title}</text>',
        f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#333"/>',
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="#333"/>',
        f'<text x="{width/2}" y="{height-15}" text-anchor="middle" font-family="Arial" font-size="13">index local</text>',
        f'<text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" font-family="Arial" font-size="13">beta_j</text>',
    ]
    for pi, (name, beta) in enumerate(profiles):
        pts = " ".join(f"{sx(i, len(beta)):.1f},{sy(float(y)):.1f}" for i, y in enumerate(beta))
        color = colors[pi % len(colors)]
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        parts.append(f'<text x="{width-margin-120}" y="{margin+18*pi}" font-family="Arial" font-size="13" fill="{color}">{name}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    root = Path("/workspace/.cache")
    outdir = Path("/workspace/bup_zoller_tomography/results")
    outdir.mkdir(parents=True, exist_ok=True)

    gs = load_mat_v5(root / "03-GroundStateBetas.mat")
    es = load_mat_v5(root / "01-ExcitedStateBetas.mat")
    rows = extract_profiles("ground", gs) + extract_profiles("excited", es)

    metrics_csv = outdir / "bup_zoller_profile_metrics.csv"
    with metrics_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Export normalized rho=1/beta profiles for the bulk experimental cells.
    rho_csv = outdir / "bup_zoller_rho_inv_beta_profiles.csv"
    with rho_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["state", "group", "cell_index", "local_index", "site", "beta", "rho_inv_beta_norm"])
        writer.writeheader()
        for state, mat in [("ground", gs), ("excited", es)]:
            group = mat["databulk"]
            for idx, (_, arr) in enumerate(group["cells"]):
                if arr.shape[1] < 4:
                    continue
                beta = arr[:, 1].astype(float)
                mask = beta > 1e-9
                rho = np.zeros_like(beta, dtype=float)
                rho[mask] = 1.0 / beta[mask]
                if rho.sum() > 0:
                    rho /= rho.sum()
                for local_i, (site, b, r) in enumerate(zip(arr[:, 0], beta, rho)):
                    writer.writerow(
                        {
                            "state": state,
                            "group": "databulk",
                            "cell_index": idx,
                            "local_index": local_i,
                            "site": float(site),
                            "beta": float(b),
                            "rho_inv_beta_norm": float(r),
                        }
                    )

    def avg(where):
        subset = [r for r in rows if where(r)]
        return {
            "count": len(subset),
            "mean_beta_cv": float(np.mean([r["beta_cv"] for r in subset])),
            "mean_r2_parabola": float(np.mean([r["r2_parabola"] for r in subset])),
            "mean_r2_triangle": float(np.mean([r["r2_triangle"] for r in subset])),
            "mean_rho_edge_to_center": float(np.mean([r["rho_edge_to_center"] for r in subset])),
            "mean_rho_edge_fraction": float(np.mean([r["rho_edge_fraction"] for r in subset])),
        }

    summary = {
        "interpretation": "rho_ent is tested here as the normalized operational proxy rho_j = (1/beta_j)/sum_k(1/beta_k), using analysed beta profiles.",
        "ground_bulk_experimental": avg(lambda r: r["state"] == "ground" and r["group"] == "databulk" and r["kind"] == "experimental"),
        "excited_bulk_experimental": avg(lambda r: r["state"] == "excited" and r["group"] == "databulk" and r["kind"] == "experimental"),
        "ground_bulk_all": avg(lambda r: r["state"] == "ground" and r["group"] == "databulk"),
        "excited_bulk_all": avg(lambda r: r["state"] == "excited" and r["group"] == "databulk"),
    }
    (outdir / "bup_zoller_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Visual sanity plots without matplotlib.
    gs_profile = gs["databulk"]["cells"][4][1][:, 1].astype(float)
    es_profile = es["databulk"]["cells"][5][1][:, 1].astype(float)
    write_svg_profile(outdir / "beta_profiles_ground_vs_excited.svg", "Zoller data: beta_j profiles", [("ground bulk L=13", gs_profile), ("excited bulk L=13", es_profile)])

    for state, profile in [("ground", gs_profile), ("excited", es_profile)]:
        mask = profile > 1e-9
        rho = np.zeros_like(profile)
        rho[mask] = 1.0 / profile[mask]
        rho /= rho.sum()
        write_svg_profile(outdir / f"rho_inv_beta_{state}.svg", f"BuP proxy rho_ent ∝ 1/beta_j ({state})", [(f"{state} rho", rho)])

    print(json.dumps(summary, indent=2))
    print(f"Wrote {metrics_csv}")
    print(f"Wrote {rho_csv}")


if __name__ == "__main__":
    main()

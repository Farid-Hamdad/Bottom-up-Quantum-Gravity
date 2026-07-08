#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 25 — Figure 03
Correction tensor hierarchy.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = Path("paper25_variational_foundation/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(15, 9))
ax.set_xlim(0, 15)
ax.set_ylim(0, 9)
ax.axis("off")


def box(x, y, w, h, text, fontsize=11):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.4,
        edgecolor="black",
        facecolor="white"
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize)
    return patch


def arrow(x1, y1, x2, y2):
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.3,
        color="black"
    )
    ax.add_patch(arr)


ax.text(
    7.5,
    8.45,
    "Correction tensor hierarchy in Bottom-Up Quantum Gravity",
    ha="center",
    va="center",
    fontsize=16,
    fontweight="bold"
)

# Top level
box(
    5.0, 7.2, 5.0, 0.75,
    r"$\mathcal{H}_{\mu\nu}^{\rm full}$"
    "\n"
    r"$=\mathcal{H}_{\mu\nu}^{\rm static}+\mathcal{H}_{\mu\nu}^{\rm dyn}$",
    12
)

# Static and dynamic
box(
    1.0, 5.9, 5.2, 0.75,
    r"$\mathcal{H}_{\mu\nu}^{\rm static}$"
    "\n"
    "Paper 20 — smooth-limit deviations",
    11
)

box(
    8.8, 5.9, 5.2, 0.75,
    r"$\mathcal{H}_{\mu\nu}^{\rm dyn}$"
    "\n"
    "Paper 22 — dynamical wave sector",
    11
)

arrow(6.2, 7.2, 3.6, 6.65)
arrow(8.8, 7.2, 11.4, 6.65)

# Static sectors
static_sectors = [
    (0.5, 4.55, r"$\mathcal{H}^{\rm spec}$", "spectral / heat-kernel residuals", "controlled residual"),
    (4.0, 4.55, r"$\mathcal{H}^{\rm curv}$", "OR-to-Ricci calibration", "controlled calibration"),
    (7.5, 4.55, r"$\mathcal{H}^{\rm source}$", "modular source residuals", "controlled residual"),
    (11.0, 4.55, r"$\mathcal{H}^{\rm dim}$", "dimension flow / SPARC / cosmology", "proxy sector"),
    (0.5, 2.95, r"$\mathcal{H}^{\rm nonlocal}$", "long-range entanglement / lensing", "proxy sector"),
    (4.0, 2.95, r"$\mathcal{H}^{\rm topo}$", "global topology / holonomies", "open sector"),
    (7.5, 2.95, r"$\mathcal{H}^{\rm finite}$", "finite-N / resolution effects", "scaling sector"),
]

for x, y, title, desc, status in static_sectors:
    box(x, y, 3.0, 0.95, title + "\n" + desc + "\n" + status, 9.5)
    arrow(3.6, 5.9, x + 1.5, y + 0.95)

# Dynamic detail
box(
    11.0, 2.95, 3.0, 0.95,
    r"$H_{\rm edge}$"
    "\nedge Hessian modes"
    "\nwave propagation",
    9.5
)
arrow(11.4, 5.9, 12.5, 3.9)

# Notes
box(
    1.0, 1.15, 6.2, 0.9,
    "Controlled sectors: spec / curv / source\n"
    "tested on synthetic smooth geometries",
    10
)

box(
    7.8, 1.15, 6.2, 0.9,
    "Phenomenological proxy sectors: dim / nonlocal\n"
    "tested on SPARC, lensing and cosmological regimes",
    10
)

ax.text(
    7.5,
    0.35,
    "Open objective: replace heterogeneous proxy diagnostics by a common tensor norm for all correction sectors.",
    ha="center",
    va="center",
    fontsize=12
)

png_path = OUTDIR / "fig03_correction_tensor_hierarchy.png"
pdf_path = OUTDIR / "fig03_correction_tensor_hierarchy.pdf"

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")

print(f"[OK] wrote {png_path}")
print(f"[OK] wrote {pdf_path}")

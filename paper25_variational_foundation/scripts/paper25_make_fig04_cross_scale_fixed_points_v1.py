#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 25 — Figure 04
Cross-scale fixed-point structure of BuP.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = Path("paper25_variational_foundation/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(15, 8.5))
ax.set_xlim(0, 15)
ax.set_ylim(0, 8.5)
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
        linewidth=1.25,
        color="black"
    )
    ax.add_patch(arr)


ax.text(
    7.5,
    8.05,
    "Cross-scale fixed-point structure in Bottom-Up Quantum Gravity",
    ha="center",
    va="center",
    fontsize=16,
    fontweight="bold"
)

# Central fixed point
box(
    5.15, 5.9, 4.7, 1.15,
    r"BuP fixed-point sector"
    "\n"
    r"$\alpha_{\rm eff}\simeq1,\quad m_{\rm eff}^2\simeq0,\quad v_g\simeq1$",
    12
)

# Surrounding scales
box(
    0.6, 6.0, 3.3, 0.95,
    "Solar system / Mercury"
    "\n"
    r"$|d_s^\odot-3|\ll1$"
    "\nfuture constraint",
    10
)

box(
    0.6, 4.35, 3.3, 0.95,
    "SPARC galaxies"
    "\n"
    r"$\alpha_{\rm eff}(r)$ organizes $V(r)$"
    "\nPapers 12--14",
    10
)

box(
    0.6, 2.7, 3.3, 0.95,
    "Cosmology"
    "\n"
    r"$d(z)\to G_{\rm eff}(z)$"
    "\nPapers 2--4",
    10
)

box(
    11.1, 6.0, 3.3, 0.95,
    "SLACS strong lensing"
    "\n"
    r"$\langle\alpha_{\rm eff}\rangle\simeq1.014$"
    "\nPaper 21",
    10
)

box(
    11.1, 4.35, 3.3, 0.95,
    "Gravitational waves"
    "\n"
    r"$m_{\rm eff}^2\simeq0,\quad v_g\simeq1$"
    "\nPaper 22",
    10
)

box(
    11.1, 2.7, 3.3, 0.95,
    "Modular tomography"
    "\n"
    r"$W_{ij}=I(i:j)$ observable"
    "\nPaper 23",
    10
)

# Arrows to center
arrow(3.9, 6.48, 5.15, 6.48)
arrow(3.9, 4.83, 5.15, 6.15)
arrow(3.9, 3.18, 5.15, 5.95)

arrow(11.1, 6.48, 9.85, 6.48)
arrow(11.1, 4.83, 9.85, 6.15)
arrow(11.1, 3.18, 9.85, 5.95)

# Underlying law
box(
    4.8, 3.65, 5.4, 0.9,
    r"Effective propagator law"
    "\n"
    r"$\alpha_{\rm eff}=\frac{2d_s}{d_w}+d_w-4$",
    12
)

arrow(7.5, 4.55, 7.5, 5.9)

# Correction sectors
box(
    4.8, 1.8, 5.4, 1.0,
    r"Correction-tensor interpretation"
    "\n"
    r"$\mathcal{H}^{\rm dim}+\mathcal{H}^{\rm nonlocal}+\mathcal{H}^{\rm dyn}$",
    12
)

arrow(7.5, 3.65, 7.5, 2.8)

# Footer
ax.text(
    7.5,
    0.55,
    "Paper 25 objective: turn these cross-scale consistencies into a controlled tensorial correction framework.",
    ha="center",
    va="center",
    fontsize=12
)

png_path = OUTDIR / "fig04_cross_scale_fixed_points.png"
pdf_path = OUTDIR / "fig04_cross_scale_fixed_points.pdf"

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")

print(f"[OK] wrote {png_path}")
print(f"[OK] wrote {pdf_path}")

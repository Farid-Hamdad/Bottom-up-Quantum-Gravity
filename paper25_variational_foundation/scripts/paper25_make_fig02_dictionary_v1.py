#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 25 — Figure 02
Discrete-to-continuum dictionary of Bottom-Up Quantum Gravity.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = Path("paper25_variational_foundation/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(15, 8))
ax.set_xlim(0, 15)
ax.set_ylim(0, 8)
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
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fontsize)
    return patch


def arrow(x1, y1, x2, y2):
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.4,
        color="black"
    )
    ax.add_patch(arr)


ax.text(
    7.5, 7.55,
    "BuP discrete-to-continuum dictionary",
    ha="center",
    va="center",
    fontsize=16,
    fontweight="bold"
)

# Headers
box(0.8, 6.55, 3.2, 0.55, "Discrete BuP object", 12)
box(5.6, 6.55, 3.2, 0.55, "Continuum target", 12)
box(10.4, 6.55, 3.7, 0.55, "Status", 12)

rows = [
    (
        r"$W_{ij}=I(i:j)$",
        r"$g_{\mu\nu}^{\rm ent}$",
        "fundamental postulate / reconstruction map"
    ),
    (
        r"$L_{\rm ent}=D-W$",
        r"$-\Delta_g$",
        "Papers 15--16: spectral convergence"
    ),
    (
        r"$L_{\rm ent}^{+}$",
        r"$(-\Delta_g)^{-1}$",
        "Papers 9--10: weak-field Green function"
    ),
    (
        r"$\kappa^{OR}_{ij}$",
        r"$R_{\mu\nu}u^\mu u^\nu$",
        "Paper 17: mean affine calibration"
    ),
    (
        r"$\delta\langle K_A\rangle$",
        r"$T_{\mu\nu}^{\rm ent}$",
        "Paper 18: modular source precursor"
    ),
    (
        r"$\alpha_{\rm eff}$",
        "weak-field / galaxies / lensing",
        "Papers 11--12--21: effective scaling law"
    ),
    (
        r"$\mathcal{H}_{\mu\nu}^{\rm BuP}$",
        "deviation from Einstein limit",
        "Papers 20--22: correction hierarchy"
    ),
]

y0 = 5.6
dy = 0.8

for i, (left, mid, right) in enumerate(rows):
    y = y0 - i * dy
    box(0.8, y, 3.2, 0.55, left, 11)
    box(5.6, y, 3.2, 0.55, mid, 11)
    box(10.4, y, 3.7, 0.55, right, 10)
    arrow(4.0, y + 0.275, 5.6, y + 0.275)
    arrow(8.8, y + 0.275, 10.4, y + 0.275)

ax.text(
    7.5,
    0.35,
    r"$\alpha_{\rm eff}=2d_s/d_w+d_w-4,\qquad "
    r"\mathcal{H}_{\mu\nu}^{\rm full}=\mathcal{H}_{\mu\nu}^{\rm static}+\mathcal{H}_{\mu\nu}^{\rm dyn}$",
    ha="center",
    va="center",
    fontsize=13
)

png_path = OUTDIR / "fig02_discrete_to_continuum_dictionary.png"
pdf_path = OUTDIR / "fig02_discrete_to_continuum_dictionary.pdf"

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")

print(f"[OK] wrote {png_path}")
print(f"[OK] wrote {pdf_path}")

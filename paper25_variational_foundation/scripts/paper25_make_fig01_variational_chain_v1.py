#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 25 — Figure 01
Variational chain of Bottom-Up Quantum Gravity.
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
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.5,
        edgecolor="black",
        facecolor="white"
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        wrap=True
    )
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


# Core variable
box(0.6, 4.0, 2.5, 1.0, r"$W_{ij}=I(i:j)$" "\n" "mutual-information graph", 12)

# First layer
box(4.0, 6.8, 2.6, 0.9, r"$L_{\rm ent}=D-W$" "\n" "spectral geometry", 11)
box(4.0, 5.2, 2.6, 0.9, r"$\kappa^{OR}_{ij}[W]$" "\n" "discrete Ricci sector", 11)
box(4.0, 3.6, 2.6, 0.9, r"$K_A=-\log\rho_A$" "\n" "modular source", 11)
box(4.0, 2.0, 2.6, 0.9, r"$L_{\rm ent}^{+}$" "\n" "Green function", 11)

# Second layer
box(7.6, 6.8, 2.8, 0.9, r"$-\Delta_g$" "\n" "Laplace--Beltrami limit", 11)
box(7.6, 5.2, 2.8, 0.9, r"$R_{\mu\nu}u^\mu u^\nu$" "\n" "Ricci signal", 11)
box(7.6, 3.6, 2.8, 0.9, r"$T_{\mu\nu}^{\rm ent}$" "\n" "stress-energy precursor", 11)
box(7.6, 2.0, 2.8, 0.9, r"$\Phi_{\rm BuP}$" "\n" "weak-field potential", 11)

# Geometry and phenomenology
box(11.0, 6.0, 3.0, 1.0, r"$g_{\mu\nu}^{\rm ent}$" "\n" "emergent geometry", 11)
box(11.0, 4.4, 3.0, 1.0, r"$G_{\mu\nu}+\Lambda g_{\mu\nu}$" "\n" "Einstein side", 11)
box(11.0, 2.8, 3.0, 1.0, r"$8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}$" "\n" "source side", 11)
box(11.0, 1.2, 3.0, 1.0, r"$\mathcal{H}_{\mu\nu}^{\rm BuP}$" "\n" "correction tensor", 11)

# Arrows from core
arrow(3.1, 4.5, 4.0, 7.25)
arrow(3.1, 4.5, 4.0, 5.65)
arrow(3.1, 4.5, 4.0, 4.05)
arrow(3.1, 4.5, 4.0, 2.45)

# Arrows first to second layer
arrow(6.6, 7.25, 7.6, 7.25)
arrow(6.6, 5.65, 7.6, 5.65)
arrow(6.6, 4.05, 7.6, 4.05)
arrow(6.6, 2.45, 7.6, 2.45)

# Arrows second to final
arrow(10.4, 7.25, 11.0, 6.55)
arrow(10.4, 5.65, 11.0, 4.95)
arrow(10.4, 4.05, 11.0, 3.3)
arrow(10.4, 2.45, 11.0, 1.75)

# Main equation title
ax.text(
    7.5,
    8.55,
    "Bottom-Up Quantum Gravity as a variational theory of entanglement equilibrium",
    ha="center",
    va="center",
    fontsize=15,
    fontweight="bold"
)

ax.text(
    7.5,
    0.35,
    r"$\frac{\delta S_{\rm BuP}[W,\rho]}{\delta W_{ij}}=0"
    r"\quad\Longrightarrow\quad"
    r"$G_{\mu\nu}[g^{\rm ent}]+\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}"
    r"=8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}+\mathcal{H}_{\mu\nu}^{\rm BuP}$",
    ha="center",
    va="center",
    fontsize=13
)

png_path = OUTDIR / "fig01_variational_chain.png"
pdf_path = OUTDIR / "fig01_variational_chain.pdf"

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")

print(f"[OK] wrote {png_path}")
print(f"[OK] wrote {pdf_path}")

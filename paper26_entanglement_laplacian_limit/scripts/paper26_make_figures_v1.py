#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Figures v1

Generates:
  fig01_limit_chain.png/pdf
  fig02_circle_spectrum_convergence.png/pdf
  fig03_normalization_scaling.png/pdf
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path("paper26_entanglement_laplacian_limit")
RESULTS = ROOT / "results" / "paper26_laplacian_limit_v1"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def savefig(name: str):
    png = FIGURES / f"{name}.png"
    pdf = FIGURES / f"{name}.pdf"
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"[OK] wrote {png}")
    print(f"[OK] wrote {pdf}")


def make_fig01_limit_chain():
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 5)
    ax.axis("off")

    def box(x, y, w, h, text, fontsize=11):
        patch = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            linewidth=1.4,
            edgecolor="black",
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize)
        return patch

    def arrow(x1, y1, x2, y2):
        arr = FancyArrowPatch(
            (x1, y1), (x2, y2),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.4,
            color="black",
        )
        ax.add_patch(arr)

    ax.text(
        7,
        4.55,
        "Paper 26 — Continuum-limit chain of the entanglement Laplacian",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )

    boxes = [
        (0.4, 2.45, 2.2, 0.9, r"$W_{ij}=I(i:j)$" + "\nmutual-information graph"),
        (3.0, 2.45, 2.2, 0.9, r"$L_{\rm ent}=D-W$" + "\nentanglement Laplacian"),
        (5.6, 2.45, 2.2, 0.9, r"$c_{N,\epsilon}L_N$" + "\nnormalized graph operator"),
        (8.2, 2.45, 2.2, 0.9, r"$-\Delta_g$" + "\nLaplace--Beltrami limit"),
        (10.8, 2.45, 2.6, 0.9, r"$(-\Delta_g)^{-1}$" + "\nGreen-function sector"),
    ]

    for x, y, w, h, text in boxes:
        box(x, y, w, h, text, 10.5)

    arrow(2.6, 2.9, 3.0, 2.9)
    arrow(5.2, 2.9, 5.6, 2.9)
    arrow(7.8, 2.9, 8.2, 2.9)
    arrow(10.4, 2.9, 10.8, 2.9)

    box(
        3.7,
        0.85,
        6.6,
        0.85,
        r"Circle result: $c_{N,\epsilon}\lambda_k(L_N)\simeq k^2$"
        "\n"
        r"with $c_{N,\epsilon}\propto 1/(N\epsilon^{3/2})$ for $D=1$",
        11,
    )

    savefig("fig01_limit_chain")
    plt.close(fig)


def make_fig02_circle_spectrum():
    df = pd.read_csv(RESULTS / "circle_spectrum_convergence.csv")

    # Use best row by mean_rel_error.
    best = df.loc[df["mean_rel_error"].idxmin()]
    num_modes = int(best["num_modes"])

    modes = list(range(num_modes))
    norm_vals = [best[f"norm_lambda_{i}"] for i in modes]
    target_vals = [best[f"target_lambda_{i}"] for i in modes]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(modes, target_vals, marker="o", label="Continuum target")
    ax.plot(modes, norm_vals, marker="s", linestyle="--", label="Normalized graph spectrum")

    ax.set_xlabel("Spectral index")
    ax.set_ylabel("Eigenvalue")
    ax.set_title(
        "Circle spectrum convergence\n"
        f"N={int(best['N'])}, epsilon={best['epsilon']}, "
        f"mean relative error={best['mean_rel_error']:.3e}"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    savefig("fig02_circle_spectrum_convergence")
    plt.close(fig)


def make_fig03_normalization_scaling():
    df = pd.read_csv(RESULTS / "normalization_scan.csv")

    fig, ax = plt.subplots(figsize=(10, 6))

    for N, group in df.groupby("N"):
        group = group.sort_values("epsilon")
        ax.plot(
            group["epsilon"],
            group["prefactor_estimate_c_times_N_eps_3_over_2"],
            marker="o",
            label=f"N={int(N)}",
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$c_{N,\epsilon}N\epsilon^{3/2}$")
    ax.set_title(
        "Normalization scaling on the circle\n"
        r"Expected stability of $c_{N,\epsilon}N\epsilon^{3/2}$"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    savefig("fig03_normalization_scaling")
    plt.close(fig)


def main():
    make_fig01_limit_chain()
    make_fig02_circle_spectrum()
    make_fig03_normalization_scaling()


if __name__ == "__main__":
    main()

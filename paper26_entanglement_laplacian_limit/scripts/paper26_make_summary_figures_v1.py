#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 26 — Summary figures v1

Generates:
  fig05_general_normalization_law.png/pdf
  fig06_prefactor_validation_summary.png/pdf
  fig07_sphere_spectrum_validation.png/pdf
"""

from pathlib import Path
import json
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


def make_fig05_general_law():
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
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
            linewidth=1.3,
            color="black",
        )
        ax.add_patch(arr)

    ax.text(
        7,
        5.55,
        "Paper 26 — General normalization law for the graph Laplacian",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
    )

    box(
        0.7, 3.55, 3.2, 0.9,
        r"$W_{ij}=\exp[-d_g(x_i,x_j)^2/(4\epsilon)]$"
        "\nlocal Gaussian graph",
        10.5
    )

    box(
        5.2, 3.55, 3.4, 0.9,
        r"$L_{N,\epsilon}=D-W$"
        "\nunnormalized Laplacian",
        11
    )

    box(
        10.0, 3.55, 3.2, 0.9,
        r"$c_{N,\epsilon}L_{N,\epsilon}\to-\Delta_g$"
        "\ncontinuum limit",
        11
    )

    arrow(3.9, 4.0, 5.2, 4.0)
    arrow(8.6, 4.0, 10.0, 4.0)

    box(
        2.0, 1.65, 10.0, 1.05,
        r"$c_{N,\epsilon}="
        r"\frac{1}{\rho(4\pi)^{D/2}\epsilon^{D/2+1}},"
        r"\quad \rho=\frac{N}{\mathrm{Vol}(M)}$",        
        15
    )

    box(
        1.0, 0.35, 3.6, 0.8,
        r"$S^1:\ c=\frac{\sqrt{\pi}}{N\epsilon^{3/2}}$",
        12
    )
    box(
        5.2, 0.35, 3.6, 0.8,
        r"$T^2:\ c=\frac{\pi}{N\epsilon^2}$",
        12
    )
    box(
        9.4, 0.35, 3.6, 0.8,
        r"$S^2:\ c=\frac{1}{N\epsilon^2}$",
        12
    )

    savefig("fig05_general_normalization_law")
    plt.close(fig)


def make_fig06_prefactor_summary():
    s1 = json.loads((RESULTS / "prefactor_theory_s1_summary.json").read_text())
    t2 = json.loads((RESULTS / "torus_fft_summary_v2.json").read_text())
    s2 = json.loads((RESULTS / "sphere_summary.json").read_text())

    labels = ["S1", "T2", "S2"]
    empirical = [
        s1["best_by_prefactor_error"]["empirical_prefactor"],
        t2["best_by_prefactor_error"]["empirical_prefactor"],
        s2["best_by_prefactor_error"]["empirical_prefactor"],
    ]
    theory = [
        s1["theory_prefactor"],
        t2["theory_prefactor"],
        s2["theory_prefactor"],
    ]

    rel_errors = [
        s1["best_by_prefactor_error"]["relative_error_to_sqrt_pi"],
        t2["best_by_prefactor_error"]["relative_error_to_pi"],
        s2["best_by_prefactor_error"]["absolute_error_to_1"],
    ]

    x = range(len(labels))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x, theory, marker="o", label="Theory prefactor")
    ax.plot(x, empirical, marker="s", linestyle="--", label="Empirical prefactor")

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Prefactor")
    ax.set_title("Prefactor validation across geometries")

    for i, err in enumerate(rel_errors):
        ax.text(i, empirical[i], f"err={100*err:.2f}%", ha="center", va="bottom", fontsize=9)

    ax.grid(True, alpha=0.3)
    ax.legend()

    savefig("fig06_prefactor_validation_summary")
    plt.close(fig)


def make_fig07_sphere_spectrum():
    df = pd.read_csv(RESULTS / "sphere_spectrum_convergence.csv")
    best = df.loc[df["mean_rel_error_theory_c"].idxmin()]
    num_modes = int(best["num_modes"])

    modes = list(range(num_modes))
    target = [best[f"target_lambda_{i}"] for i in modes]
    theory_norm = [best[f"norm_theory_lambda_{i}"] for i in modes]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(modes, target, marker="o", label="Continuum S2 spectrum")
    ax.plot(modes, theory_norm, marker="s", linestyle="--", label="Graph spectrum, theory normalization")

    ax.set_xlabel("Spectral index")
    ax.set_ylabel("Eigenvalue")
    ax.set_title(
        "Sphere spectrum validation\n"
        f"N={int(best['N'])}, epsilon={best['epsilon']}, "
        f"mean relative error={best['mean_rel_error_theory_c']:.3e}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend()

    savefig("fig07_sphere_spectrum_validation")
    plt.close(fig)


def main():
    make_fig05_general_law()
    make_fig06_prefactor_summary()
    make_fig07_sphere_spectrum()


if __name__ == "__main__":
    main()

# Paper 11 — The Dimensional Flow of BuP Gravity

## Spectral Dimension, Walk Dimension and the Effective Gravitational Exponent

**Author:** Farid Hamdad
**Project:** Bottom-Up Quantum Gravity
**Year:** 2026

---

## Overview

This folder contains the numerical and theoretical material associated with **Paper 11** of the Bottom-Up Quantum Gravity program.

Paper 11 introduces the dimensional law connecting the spectral properties of an entanglement graph to an effective gravitational exponent:

[
\alpha_{\rm eff}
================

\frac{2d_s}{d_w}
+
d_w
---

4.

]

Here:

* (d_s) is the spectral dimension of the entanglement graph;
* (d_w) is the walk dimension;
* (\alpha_{\rm eff}) is the effective gravitational exponent controlling the emergent large-scale response.

This paper is the bridge between the microscopic graph program and the phenomenological tests performed later in the BuP sequence.

---

## Role in the BuP program

Paper 11 provides the dimensional backbone for later papers:

```text
Paper 11:
    W_ij → L_ent → d_s, d_w → alpha_eff

Paper 12:
    Sigma(R) → W_ij → L_ent → alpha_eff → V(r)

Paper 14:
    SPARC galaxy rotation curves and LOW/HIGH regimes

Paper 21:
    SLACS strong-lensing fixed point and alpha_eff ≈ 1
```

Thus, Paper 11 does not primarily fit a galaxy catalogue. It establishes the theoretical and numerical law that later papers test observationally.

---

## Central equation

The central result is:

[
\boxed{
\alpha_{\rm eff}
================

\frac{2d_s}{d_w}
+
d_w
---

4
}
]

This formula combines two graph-diffusion quantities:

[
P(t)\sim t^{-d_s/2},
]

and

[
\langle r^2(t)\rangle\sim t^{2/d_w}.
]

The BuP gravitational exponent is therefore not inserted by hand. It is inferred from diffusion on the entanglement graph.

---

## Newtonian fixed point

The Newtonian or baryonic fixed point corresponds to:

[
\alpha_{\rm eff}=1.
]

Therefore,

[
\frac{2d_s}{d_w}+d_w-4=1.
]

Equivalently,

[
2d_s=d_w(5-d_w),
]

or

[
d_s=\frac{d_w(5-d_w)}{2}.
]

For Brownian diffusion,

[
d_w=2,
]

this gives:

[
d_s=3.
]

Thus, ordinary three-dimensional Newtonian behavior appears as a special fixed point of the BuP dimensional flow.

---

## Main numerical outputs

The expected outputs are:

```text
results/
  alpha_eff_table.csv
  finite_size_summary.csv
  alpha_predictions_vs_N.csv
  paper11_summary.json
```

and the figures:

```text
figures/
  fig1_pipeline_dimensional_flow.png
  fig2_ds_vs_N.png
  fig3_dw_vs_N.png
  fig4_alpha_predictions_vs_N.png
  fig5_alpha_fixed_point_curve.png
  fig6_interpretation_regimes.png
```

---

## Interpretation

Paper 11 shows that the effective gravitational behavior is controlled by the pair:

[
(d_s,d_w).
]

Different regimes correspond to different gravitational responses:

| Regime                | Condition                  | Interpretation                      |
| --------------------- | -------------------------- | ----------------------------------- |
| Newtonian fixed point | (\alpha_{\rm eff}=1)       | ordinary baryonic/Newtonian scaling |
| sub-Newtonian         | (\alpha_{\rm eff}<1)       | softened or under-coupled regime    |
| super-Newtonian       | (\alpha_{\rm eff}>1)       | enhanced effective response         |
| transition regime     | (\alpha_{\rm eff}\approx1) | crossover between graph phases      |

---

## Reproducibility

Run:

```bash
cd papers/paper11_dimensional_flow
bash scripts/run_all.sh
```

This generates all numerical tables and figures.

---

## Status

Paper 11 should be read as a theoretical and numerical derivation of the BuP dimensional law. Its observational consequences are tested later in Paper 12, Paper 14 and Paper 21.

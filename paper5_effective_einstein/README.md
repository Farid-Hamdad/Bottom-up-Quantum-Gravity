# BuP Paper 5 — Effective Einstein Equation from Entanglement Geometry

**Title FR:** Équation d'Einstein effective émergente dans la Gravité Quantique Bottom-Up (BuP)

**Title EN:** Effective Emergent Einstein Equation in Bottom-Up Quantum Gravity (BuP)

**Author:** Farid Hamdad  
**Year:** 2026

---

## Core equation

$$
G_{\mu\nu}[g^{\rm ent}]
=
8\pi G_{\rm eff}(d_s)T_{\mu\nu}^{\rm matter}
+
8\pi G\,T_{\mu\nu}^{\rm ent}[d_s]
$$

No cosmological constant is postulated, and no dark-energy fluid is added.  
The late-time acceleration is interpreted as a geometric consequence of the effective dimension \(d_{\rm bg}(z)<3\).

---

## Numerical internal consistency test

Discrete curvature proxies are compared with spectral-dimension stress proxies on BuP mutual-information graphs (\(N=9,12,16,20\)).

Most stable relation found for \(N=20\), across \(k=3,5,7\):

$$
G_{ij}^{\rm proxy}
=
-\kappa_{ij}
+
\tfrac{1}{2}R_{\rm edge}
$$

This curvature proxy is monotonically correlated with spectral proxies constructed from:

$$
(\nabla d_s)^2,
\qquad
\Delta d_s.
$$

| \(k\) | Dominant spectral proxy | \(\langle r \rangle\) | \(\langle\rho\rangle\) | Spearman sign |
|---|---|---:|---:|---:|
| 3 | \((\nabla d_s)^2+|\Delta d_s|\) | 0.707 | 0.598 | 100% |
| 5 | \((\nabla d_s)^2+|\Delta d_s|\) | 0.558 | 0.566 | 100% |
| 7 | \((\nabla d_s)^2-\Delta d_s\) | 0.995 | 0.417 | 100% |

The Spearman correlation sign is positive in **100% of tested MI matrices** for all tested \(k\)-values.

This provides a first internal numerical consistency test for the proposed BuP effective Einstein equation. It does not constitute a full derivation of the continuum tensor equation.

---

## Structure

```text
paper/
└── main_fr.tex                         ← Full French manuscript with numerical section

figures/
├── fig_best_G_vs_T_N20.png              ← Best curvature-stress proxy correlation
└── fig_top10_stability_N20.png          ← Proxy stability ranking

scripts/
├── bup_einstein_tensor_correlation_v1.py
└── bup_einstein_tensor_correlation_v2.py

results/
├── results_einstein_corr_v2/
├── results_einstein_corr_N20_v2/
├── results_einstein_corr_N20_k3_v2/
├── results_einstein_corr_N20_k5_v2/
└── results_einstein_corr_N20_k7_v2/

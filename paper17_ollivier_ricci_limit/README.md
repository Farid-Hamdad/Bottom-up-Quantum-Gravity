# Paper 17 — Ollivier--Ricci Curvature and the Ricci Limit

**From discrete Ollivier--Ricci curvature to the continuum Ricci tensor**

---

## Status

Paper 17 studies the second analytic arrow required for the Einstein limit of BuP:

\[
\kappa_{ij}^{OR}
\longrightarrow
R_{\mu\nu}u^\mu u^\nu.
\]

Paper 15 showed a first positive relative signal:

\[
\bar\kappa_{\rm sphere}
>
\bar\kappa_{\rm flat}.
\]

Paper 17 strengthens this by testing how the Ollivier--Ricci curvature scales with the local graph scale and whether it can be renormalized into a Ricci proxy.

---

## 1. Starting point

The discrete curvature used in BuP is the Ollivier--Ricci curvature:

\[
\kappa_{ij}^{OR}
=
1-
\frac{W_1(m_i,m_j)}{d(i,j)},
\]

where:

- \(W_1(m_i,m_j)\) is the Wasserstein-1 distance between neighborhood measures,
- \(d(i,j)\) is the graph/geodesic distance between nodes,
- \(m_i\) and \(m_j\) are local probability measures around nodes \(i\) and \(j\).

The expected continuum target is:

\[
\kappa_{ij}^{OR}
\sim
\epsilon_N
R_{\mu\nu}u^\mu u^\nu
+
\text{finite-scale bias}.
\]

Thus, the natural normalized quantity is:

\[
\frac{\kappa^{OR}}{\epsilon_N}.
\]

---

## 2. Link with Paper 15

Paper 15 tested Ollivier--Ricci curvature on controlled geometries.

The key result was:

\[
\bar\kappa_{\rm flat}=-0.001671\simeq0,
\]

\[
\bar\kappa_{\rm sphere}=0.008854,
\]

so:

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.010525.
\]

This provided the first relative numerical signal for:

\[
\kappa_{ij}^{OR}
\to
R_{\mu\nu}u^\mu u^\nu.
\]

Paper 17 now asks a sharper question:

\[
\frac{\kappa^{OR}}{\epsilon_N}
=
B_N
+
C_N R_{\mu\nu}u^\mu u^\nu
+
o(1)?
\]

---

## 3. Main numerical results

Paper 17 currently has two positive tests.

| Step | Quantity | \(N\) | Value | Reference | Status |
|---|---|---:|---:|---|---|
| v1 relative OR signal | sphere - flat raw mean kappa | 256 | 0.013383 | \(>0\) | positive |
| v1 relative OR signal | sphere - flat mean kappa/epsilon | 256 | 0.197719 | \(>0\) | positive |
| v1 relative OR signal | sphere - flat mean kappa/\(\ell^2\) | 256 | 0.557041 | \(>0\) | positive |
| v2 affine calibration | \(B_N\) flat baseline for kappa/epsilon | 512 | -0.301281 | finite-\(N\) bias | identified |
| v2 affine calibration | sphere-flat delta for kappa/epsilon | 512 | 0.391933 | \(>0\) | positive |
| v2 affine calibration | \(A_N\) epsilon calibration | 512 | 2.551455 | \(1/\Delta\) | identified |

---

## 4. Test v1 — Relative sphere-minus-flat signal

The v1 test compares three quantities:

\[
\kappa^{OR},
\]

\[
\frac{\kappa^{OR}}{\epsilon_N},
\]

\[
\frac{\kappa^{OR}}{\ell_{ij}^2}.
\]

At \(N=256\), the sphere-minus-flat signals are:

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.013383,
\]

\[
\Delta
\left\langle
\frac{\kappa}{\epsilon}
\right\rangle_{\rm sphere-flat}
=
0.197719,
\]

\[
\Delta
\left\langle
\frac{\kappa}{\ell^2}
\right\rangle_{\rm sphere-flat}
=
0.557041.
\]

All three are positive.

This confirms that the Ollivier--Ricci signal is not an artifact of one specific normalization. The sphere consistently carries a larger curvature signal than the flat torus.

---

## 5. Test v2 — Affine Ricci calibration

The v2 test introduces the calibrated Ricci proxy:

\[
\widehat{R}_{OR}
=
A_N
\left(
\frac{\kappa^{OR}}{\epsilon}
-
B_N
\right).
\]

The flat torus is used as the zero-Ricci reference:

\[
B_N=
\left\langle
\frac{\kappa^{OR}}{\epsilon}
\right\rangle_{\rm flat}.
\]

The unit sphere is used as the reference:

\[
R_{\mu\nu}u^\mu u^\nu=1.
\]

Thus:

\[
A_N
=
\frac{1}{
\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm sphere}
-
\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm flat}
}.
\]

The calibration constants are:

| \(N_{\rm input}\) | \(B_N\) flat baseline | sphere-flat \(\Delta\) | \(A_N\) |
|---:|---:|---:|---:|
| 128 | -0.2103 | 0.3053 | 3.2760 |
| 256 | -0.3079 | 0.4005 | 2.4968 |
| 512 | -0.3013 | 0.3919 | 2.5515 |

At \(N=512\):

\[
B_N=-0.301281,
\]

\[
C_N=\Delta_{\rm sphere-flat}=0.391933,
\]

\[
A_N=\frac{1}{C_N}=2.551455.
\]

Thus:

\[
\widehat{R}_{OR}
=
2.551455
\left(
\frac{\kappa^{OR}}{\epsilon}
+
0.301281
\right).
\]

By construction:

\[
\left\langle
\widehat{R}_{OR}
\right\rangle_{\rm flat}
=
0,
\]

\[
\left\langle
\widehat{R}_{OR}
\right\rangle_{\rm sphere}
=
1.
\]

The important finite-scale relation is:

\[
\boxed{
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N
+
C_N R_{\mu\nu}u^\mu u^\nu
}
\]

with, at \(N=512\):

\[
B_N\simeq-0.301,
\qquad
C_N\simeq0.392.
\]

---

## 6. Interpretation

Paper 17 v2 identifies a finite-scale affine renormalization of the Ollivier--Ricci signal.

The result supports a mean-level Ricci proxy:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N
+
C_N R_{\mu\nu}u^\mu u^\nu.
\]

This is not yet a pointwise convergence theorem.

The calibration works on geometric averages:

\[
\langle \widehat{R}_{OR}\rangle_{\rm flat}=0,
\]

\[
\langle \widehat{R}_{OR}\rangle_{\rm sphere}=1.
\]

However, the distributional dispersion remains significant, especially on the flat torus. Therefore, Paper 17 currently establishes a **mean Ricci proxy**, not yet a local pointwise Ricci estimator.

---

## 7. Limits

The current result has several limits:

1. The calibration is affine and empirical.
2. The result is mean-level, not pointwise.
3. The flat torus distribution has large dispersion.
4. The \(k\)-NN construction may introduce anisotropy.
5. The exact continuum coefficient \(C_N\) is not analytically derived.
6. Larger \(N\) tests require more efficient optimal transport.
7. More geometries are needed: sphere radii, hyperbolic surfaces, product manifolds.

---

## 8. Folder structure

```text
papers/paper17_ollivier_ricci_limit/
  README.md
  paper17_ollivier_ricci_limit.tex

  scripts/
    paper17_or_ricci_scaling_v1.py
    paper17_or_ricci_scaling_v2.py
    paper17_build_ricci_summary_v1.py

  results/
    or_ricci_scaling_v1/
    or_ricci_scaling_v2/
    paper17_ricci_summary_v1/

  figures/
    # final figures copied from selected result folders

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

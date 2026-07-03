# Paper 17 — Ricci Summary v1

This summary gathers the first Paper 17 tests of the continuum Ricci limit of Ollivier--Ricci curvature.

The target relation is:

\[
\kappa_{ij}^{OR}
\longrightarrow
R_{\mu\nu}u^\mu u^\nu.
\]

The finite-scale diagnostic suggested by the tests is:

\[
\frac{\kappa^{OR}}{\epsilon}
=
B_N+C_N R_{\mu\nu}u^\mu u^\nu+o(1).
\]

---

## 1. Key results

| Step | Quantity | N | Value | Reference | Status |
|---|---|---:|---:|---|---|
| v1 relative OR signal | sphere - flat raw mean kappa | 256 | 0.013383 | > 0 | positive |
| v1 relative OR signal | sphere - flat mean kappa/epsilon | 256 | 0.197719 | > 0 | positive |
| v1 relative OR signal | sphere - flat mean kappa/l^2 | 256 | 0.557041 | > 0 | positive |
| v2 affine calibration | B_N flat baseline for kappa/epsilon | 512 | -0.301281 | finite-N bias | identified |
| v2 affine calibration | sphere-flat delta for kappa/epsilon | 512 | 0.391933 | > 0 | positive |
| v2 affine calibration | A_N epsilon calibration | 512 | 2.551455 | 1/delta | identified |

---

## 2. v1 — Relative sphere-minus-flat signal

The first v1 test confirmed that the sphere has a positive Ollivier--Ricci signal relative to the flat torus.

At \(N=256\):

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.013383,
\]

\[
\Delta\left\langle\frac{\kappa}{\epsilon}\right\rangle_{\rm sphere-flat}=0.197719,
\]

\[
\Delta\left\langle\frac{\kappa}{\ell^2}\right\rangle_{\rm sphere-flat}=0.557041.
\]

Interpretation: the sign of the Ricci signal is correct, but the raw normalizations are not yet unbiased.

---

## 3. v2 — Affine Ricci calibration

The v2 test introduces an affine calibration:

\[
\widehat{R}_{OR}
=
A_N\left(
\frac{\kappa^{OR}}{\epsilon}-B_N
\right).
\]

The flat torus is used as the Ricci-zero reference:

\[
B_N=
\left\langle\frac{\kappa^{OR}}{\epsilon}\right\rangle_{\rm flat}.
\]

The unit sphere is used as the \(R_{\mu
u}u^\mu u^
u=1\) reference:

\[
A_N=
\frac{1}{
\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm sphere}
-
\left\langle\kappa^{OR}/\epsilon\right\rangle_{\rm flat}
}.
\]

The calibration constants are:

| \(N_{
m input}\) | \(B_N\) flat baseline | sphere-flat \(\Delta\) | \(A_N\) |
|---:|---:|---:|---:|
| 128 | -0.2103 | 0.3053 | 3.2760 |
| 256 | -0.3079 | 0.4005 | 2.4968 |
| 512 | -0.3013 | 0.3919 | 2.5515 |

At \(N=512\), this gives:

\[
\widehat{R}_{OR}=2.551
\left(
\frac{\kappa^{OR}}{\epsilon}+0.301
\right).
\]

Between \(N=256\) and \(N=512\), the calibration is relatively stable:

| Quantity | N=256 | N=512 | absolute change |
|---|---:|---:|---:|
| \(B_N\) | -0.3079 | -0.3013 | 0.0066 |
| sphere-flat \(\Delta\) | 0.4005 | 0.3919 | 0.0086 |
| \(A_N\) | 2.4968 | 2.5515 | 0.0546 |

---

## 4. Interpretation

Paper 17 v2 identifies a finite-scale affine renormalization of Ollivier--Ricci curvature.

The result supports the mean-level relation:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_N R_{\mu\nu}u^\mu u^\nu.
\]

This is not yet a pointwise convergence theorem. The calibrated means are correct by construction, while the distributional dispersion remains large, especially on the flat torus.

The next step is to reduce this dispersion by testing radius graphs, edge-length bins, trimmed means, and larger \(N\).

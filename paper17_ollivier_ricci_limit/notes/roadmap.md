# Paper 17 — Roadmap

## Objective

Paper 17 studies the continuum Ricci limit of Ollivier--Ricci curvature:

\[
\kappa_{ij}^{OR}
\longrightarrow
R_{\mu\nu}u^\mu u^\nu.
\]

The numerical target is to identify a finite-scale relation of the form:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N
+
C_N R_{\mu\nu}u^\mu u^\nu.
\]

## Current numerical status

### v1 — Relative sphere-minus-flat signal

At \(N=256\):

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.013383.
\]

\[
\Delta\left\langle\frac{\kappa}{\epsilon}\right\rangle_{\rm sphere-flat}=0.197719.
\]

\[
\Delta\left\langle\frac{\kappa}{\ell^2}\right\rangle_{\rm sphere-flat}=0.557041.
\]

All three relative signals are positive.

## v2 — Affine Ricci calibration

The calibrated proxy is:

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

This maps the flat torus to mean Ricci \(0\) and the unit sphere to mean Ricci \(1\).

## Interpretation

Paper 17 currently establishes a mean Ricci proxy, not a pointwise convergence theorem.

## Next steps

1. Reduce dispersion on the flat torus.
2. Test radius graphs instead of \(k\)-NN graphs.
3. Bin edges by length \(\ell_{ij}\).
4. Use trimmed means or robust estimators.
5. Increase \(N\) using optimized optimal transport.
6. Test spheres of different radii:
   \[
   R_{\mu\nu}u^\mu u^\nu=\frac{1}{r^2}.
   \]
7. Test negatively curved geometries.
8. Derive \(B_N\) and \(C_N\) analytically.

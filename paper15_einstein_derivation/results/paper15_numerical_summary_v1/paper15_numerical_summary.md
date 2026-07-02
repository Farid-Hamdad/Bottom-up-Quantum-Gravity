# Paper 15 — Numerical Summary v1

This file summarizes the first three numerical tests of Paper 15.

The target chain is:

\[
W_{ij}\to L_\epsilon\to\Delta_g,
\qquad
\kappa_{ij}^{\rm OR}\to R_{\mu\nu}u^\mu u^\nu,
\qquad
\delta W_{\rm loc}\to\delta\kappa(r).
\]

---

## 1. Key results table

| step                          | quantity                                 |     value | target_or_reference         | error_or_delta               | status            |
|:------------------------------|:-----------------------------------------|----------:|:----------------------------|:-----------------------------|:------------------|
| Step B — spectral convergence | circle: d_s                              |  1.0021   | 1                           | 0.0021                       | positive          |
| Step B — spectral convergence | grid2d: d_s                              |  1.9938   | 2                           | 0.0062                       | positive          |
| Step B — spectral convergence | interval: d_s                            |  0.9445   | 1                           | 0.0555                       | positive          |
| Step B — spectral convergence | sphere: d_s                              |  1.8964   | 2                           | 0.1036                       | positive          |
| Step C — Ricci signal         | sphere mean kappa                        |  0.008854 | flat torus mean = -0.001671 | delta = 0.010525             | positive_relative |
| Step E — source response      | flat_torus2d: near/far response, s=-0.30 |  7.75     | near/far > 1                | Spearman(phi,|dkappa|)=0.583 | positive          |
| Step E — source response      | sphere: near/far response, s=-0.30       | 11.51     | near/far > 1                | Spearman(phi,|dkappa|)=0.752 | positive          |

---

## 2. Step B — Spectral convergence

The continuum-scaled Laplacian

\[
L_\epsilon=\frac{D-W}{\epsilon}
\]

recovers the expected heat-trace spectral dimension on controlled geometries.

| Geometry | Target dimension | Measured \(d_s\) | Error |
|---|---:|---:|---:|
| circle | 1 | 1.0021 | 0.0021 |
| grid2d | 2 | 1.9938 | 0.0062 |
| interval | 1 | 0.9445 | 0.0555 |
| sphere | 2 | 1.8964 | 0.1036 |

Interpretation: this is the first positive numerical evidence for

\[
L_N\longrightarrow \Delta_g.
\]

---

## 3. Step C — Ollivier--Ricci curvature signal

The flat periodic torus gives

\[
\bar\kappa_{\rm flat}=-0.001671,
\]

while the sphere gives

\[
\bar\kappa_{\rm sphere}=0.008854.
\]

The relative curvature excess is

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.010525.
\]

Interpretation: this provides a first relative numerical signal for

\[
\kappa_{ij}^{\rm OR}\longrightarrow R_{\mu\nu}u^\mu u^\nu.
\]

---

## 4. Step E — Source-response test

A smooth radial perturbation of the entanglement weights is applied:

\[
\phi_i=\exp\left[-\frac{d(i,\mathrm{source})^2}{2\sigma^2}\right],
\]

\[
W'_{ij}=W_{ij}\left[1+s\frac{\phi_i+\phi_j}{2}\right].
\]

For \(s=-0.30\), the key source-response results are:

| Geometry | near/far response | Spearman \((\phi,|\Delta\kappa|)\) | p-value |
|---|---:|---:|---:|
| flat_torus2d | 7.75 | 0.583 | 4.12e-83 |
| sphere | 11.51 | 0.752 | 5.49e-165 |

Interpretation: a local entanglement defect generates a localized curvature response:

\[
\delta W_{\rm loc}\longrightarrow \delta\kappa(r).
\]

This connects Paper 15 directly to the Paper 8 construction of
\(T_{\mu\nu}^{\rm eff}\) from local entanglement perturbations.

---

## 5. Current status

Paper 15 now has three positive numerical pillars:

1. \(L_\epsilon=(D-W)/\epsilon\) recovers the expected spectral dimension.
2. Ollivier--Ricci curvature gives a positive sphere-minus-flat signal.
3. A radial entanglement defect produces a localized curvature response.

These results do not yet prove the continuum Einstein equation, but they support the three main arrows required for the effective derivation.

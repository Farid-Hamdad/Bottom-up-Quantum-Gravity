# Paper 19 — Einstein Limit Summary v1

This file assembles the numerical pillars needed for the effective Einstein limit of Bottom-Up Quantum Gravity.

The target equation is

\[
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}.
\]

In the smooth low-energy limit, the correction tensor is expected to vanish:

\[
\mathcal{H}_{\mu\nu}\to0.
\]

---

## 1. Three pillars

### Paper 16 — Spectral geometry

Paper 16 supports

\[
c_NL_N\to-\Delta_g.
\]

Key results:

| Geometry | Mean relative spectral error |
|---|---:|
| \(S^1\) | 0.005072 |
| \(T^2\) | 0.026627 |

Interpretation: the entanglement Laplacian reconstructs the low spectrum of the Laplace--Beltrami operator.

### Paper 17 — Ricci curvature

Paper 17 supports

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_NR_{\mu\nu}u^\mu u^\nu.
\]

At \(N=512\):

\[
B_N=-0.301281,
\qquad
C_N=0.391933,
\qquad
A_N=2.551455.
\]

Interpretation: Ollivier--Ricci curvature carries a calibrated mean Ricci signal.

### Paper 18 — Modular source

Paper 18 supports

\[
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r).
\]

Graph modular first-law slopes:

| Geometry | slope \(\delta S_A\) vs \(\delta\langle K_A
angle\) |
|---|---:|
| flat torus | 0.989481 |
| sphere | 0.989718 |

Modular source to curvature response:

| Geometry | \(R^2(|\delta K|,\langle|\Delta\kappa|
angle_{
m near})\) |
|---|---:|
| flat torus | 0.967229 |
| sphere | 0.982252 |

Interpretation: the modular response behaves as an effective source for curvature.

---

## 2. Assembly

The BuP variational equation is

\[
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0.
\]

The action is

\[
S_{\rm BuP}[W]
=
\mathrm{Tr}\,L(W)^{-\beta[W]}
+
\sum_{ij}W_{ij}\kappa_{ij}[W]
+
\lambda\sum_{ij}W_{ij}d_{ij}^2
+
S_{\rm topo}[W].
\]

Using Papers 16--18, the continuum identifications are:

\[
L(W)\to-\Delta_g,
\qquad
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
\qquad
\delta\langle K_A\rangle\to T_{\mu\nu}^{\rm ent}.
\]

Therefore, the expected continuum equation is

\[
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}.
\]

---

## 3. Key results table

| Paper | Pillar | Arrow | Test | Quantity | Geometry | Value | Target | Status |
|---|---|---|---|---|---|---:|---|---|
| Paper 16 | spectral geometry | `L_N -> Delta_g` | Circle spectrum | mean relative spectral error | Circle | 0.005072 | small | positive |
| Paper 16 | spectral geometry | `L_N -> Delta_g` | Circle spectrum | lambda1 scaled | Circle | 1.006212 | 1.0 | positive |
| Paper 16 | spectral geometry | `L_N -> Delta_g` | Flat torus spectrum | mean relative spectral error | Flat torus | 0.026627 | small | positive |
| Paper 16 | spectral geometry | `L_N -> Delta_g` | Flat torus spectrum | lambda1 scaled | Flat torus | 40.961270 | 39.47841760435743 | positive |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v1 relative OR signal | sphere - flat raw mean kappa | flat_torus2d / sphere | 0.013383 | > 0 | positive |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v1 relative OR signal | sphere - flat mean kappa/epsilon | flat_torus2d / sphere | 0.197719 | > 0 | positive |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v1 relative OR signal | sphere - flat mean kappa/l^2 | flat_torus2d / sphere | 0.557041 | > 0 | positive |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v2 affine calibration | B_N flat baseline for kappa/epsilon | flat_torus2d / sphere | -0.301281 | finite-N bias | identified |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v2 affine calibration | sphere-flat delta for kappa/epsilon | flat_torus2d / sphere | 0.391933 | > 0 | positive |
| Paper 17 | curvature | `kappa_OR -> Ricci` | v2 affine calibration | A_N epsilon calibration | flat_torus2d / sphere | 2.551455 | 1/delta | identified |
| Paper 18 | source modular first law | `delta W_loc -> delta S_A ~= delta<K_A>` | v1 graph modular first law | slope deltaS vs delta<K> | flat_torus2d | 0.989481 | 1 | positive |
| Paper 18 | source modular first law | `delta W_loc -> delta S_A ~= delta<K_A>` | v1 graph modular first law | R2 deltaS vs delta<K> | flat_torus2d | 0.996590 | near 1 | positive |
| Paper 18 | source modular first law | `delta W_loc -> delta S_A ~= delta<K_A>` | v1 graph modular first law | slope deltaS vs delta<K> | sphere | 0.989718 | 1 | positive |
| Paper 18 | source modular first law | `delta W_loc -> delta S_A ~= delta<K_A>` | v1 graph modular first law | R2 deltaS vs delta<K> | sphere | 0.996592 | near 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | near/far localization ratio | flat_torus2d | 6.293007 | > 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | R2 |deltaK| -> near |delta kappa| | flat_torus2d | 0.967229 | near 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | signed Pearson deltaK -> near delta kappa | flat_torus2d | -0.999168 | |Pearson| near 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | near/far localization ratio | sphere | 9.533529 | > 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | R2 |deltaK| -> near |delta kappa| | sphere | 0.982252 | near 1 | positive |
| Paper 18 | source to curvature | `delta<K_A> -> delta kappa(r)` | v2 modular source curvature | signed Pearson deltaK -> near delta kappa | sphere | -0.999428 | |Pearson| near 1 | positive |

---

## 4. Interpretation

Paper 19 does not introduce an isolated new numerical test. It assembles the three validated continuum arrows:

\[
L_N\to\Delta_g,
\qquad
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
\qquad
\delta W_{\rm loc}\to T_{\mu\nu}^{\rm ent}.
\]

Together, these support the claim that the BuP entanglement equilibrium admits an effective Einstein limit.

The remaining task is analytic: derive the coefficients, control the correction tensor \(\mathcal{H}_{\mu
u}\), and specify the low-energy/smooth regime in which it vanishes.

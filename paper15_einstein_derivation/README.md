# Paper 15 — From Entanglement Graphs to Einstein Equations

**Dérivation effective de la gravité continue dans Bottom-Up Quantum Gravity**

---

## Status

Paper 15 studies whether the discrete BuP variational equation on the entanglement network,

\[
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0,
\]

admits a smooth continuum limit of Einstein type:

\[
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

The goal is not to postulate general relativity, but to recover it as a stable continuum limit of an equilibrium condition on the entanglement graph.

---

## 1. Starting point

The fundamental object is the mutual-information matrix:

\[
W_{ij}=I(i:j).
\]

From it, one constructs the entanglement Laplacian:

\[
L_{\rm ent}=D-W,
\qquad
D_{ii}=\sum_j W_{ij}.
\]

Paper 15 starts from the minimal BuP action:

\[
S_{\rm BuP}[W]
=
\mathrm{Tr}\,L(W)^{-\beta[W]}
+
\sum_{ij}W_{ij}\kappa_{ij}[W]
+
\lambda\sum_{ij}W_{ij}d_{ij}^{2}
+
S_{\rm topo}[W].
\]

The four terms are:

| Term | Meaning | Role |
|---|---|---|
| \(\mathrm{Tr}\,L^{-\beta}\) | spectral action | global geometry |
| \(\sum W_{ij}\kappa_{ij}\) | discrete curvature | local Ricci response |
| \(\lambda\sum W_{ij}d_{ij}^2\) | locality cost | suppresses nonlocal complete graphs |
| \(S_{\rm topo}[W]\) | topology | global connectivity constraints |

Paper 14 fixed the spectral exponent:

\[
\beta
=
F(\lambda_2,d_s,d_w,\Delta_3),
\]

so the spectral term is no longer free.

---

## 2. Central objective

Paper 15 tests the chain:

\[
W_{ij}
\to
L_\epsilon
\to
\Delta_g,
\]

\[
\kappa_{ij}^{\rm OR}
\to
R_{\mu\nu}u^\mu u^\nu,
\]

\[
\delta W_{\rm loc}
\to
\delta\kappa(r)
\to
T_{\mu\nu}^{\rm eff}.
\]

If these limits hold, the continuum action should take the schematic form:

\[
S_{\rm cont}[g]
=
\int d^Dx\sqrt{|g|}
\left[
\frac{1}{16\pi G_{\rm eff}}R
+
\Lambda_{\rm ent}
+
\mathcal{L}_{\rm ent-source}
+
\mathcal{H}
\right],
\]

where \(\mathcal{H}\) contains higher-curvature or nonlocal spectral corrections.

In the smooth low-energy limit,

\[
\mathcal{H}_{\mu\nu}\to0,
\]

and one expects:

\[
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

---

## 3. Numerical status — first three positive tests

Paper 15 currently has three positive numerical pillars.

---

### Step B — Spectral convergence

The first task is to test:

\[
L_N\longrightarrow\Delta_g.
\]

A first attempt using the unscaled normalized graph Laplacian,

\[
L_{\rm norm}=I-D^{-1/2}WD^{-1/2},
\]

failed on 2D geometries, with both grid and sphere collapsing toward:

\[
d_s\simeq1.2.
\]

The corrected version uses the continuum-rescaled Laplacian:

\[
L_\epsilon=\frac{D-W}{\epsilon}.
\]

Using \(k=\sqrt{N}\), \(\epsilon=0.5\epsilon_{\rm knn}\), and \(N=1024\), the measured spectral dimensions are:

| Geometry | Target dimension | Measured \(d_s\) | Error |
|---|---:|---:|---:|
| circle | 1 | 1.0021 | 0.0021 |
| interval | 1 | 0.9445 | 0.0555 |
| grid2d | 2 | 1.9938 | 0.0062 |
| sphere | 2 | 1.8964 | 0.1036 |

This provides the first positive numerical support for:

\[
L_N\longrightarrow\Delta_g.
\]

---

### Step C — Ollivier--Ricci curvature signal

The second task is to test whether discrete Ollivier--Ricci curvature detects the continuum Ricci signal:

\[
\kappa_{ij}^{\rm OR}
\longrightarrow
R_{\mu\nu}u^\mu u^\nu.
\]

A first version using a flat grid with boundary was contaminated by edge effects. The corrected version uses a flat periodic torus as the zero-curvature reference and geodesic distances for the sphere.

At the final resolution:

\[
\bar\kappa_{\rm flat}=-0.001671\simeq0,
\]

while:

\[
\bar\kappa_{\rm sphere}=0.008854.
\]

The relative curvature excess is:

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.010525.
\]

This gives a first relative numerical signal for:

\[
\kappa_{ij}^{\rm OR}
\to
R_{\mu\nu}u^\mu u^\nu.
\]

The result is still qualitative: it distinguishes flat periodic geometry from positive curvature, but does not yet prove pointwise convergence to the Ricci tensor.

---

### Step E — Source-response test

The third task is to test:

\[
\delta W_{\rm loc}
\longrightarrow
\delta\kappa(r).
\]

The source-response v2 test uses a smooth radial entanglement perturbation:

\[
\phi_i
=
\exp\left[
-\frac{d(i,\mathrm{source})^2}{2\sigma^2}
\right],
\]

and modifies the weights as:

\[
W'_{ij}
=
W_{ij}
\left[
1+s\frac{\phi_i+\phi_j}{2}
\right].
\]

The curvature response is:

\[
\Delta\kappa_{ij}
=
\kappa_{ij}^{\rm after}
-
\kappa_{ij}^{\rm before}.
\]

For \(s=-0.30\), the response is strongly localized:

| Geometry | near/far response | Spearman \((\phi_{\rm edge},|\Delta\kappa|)\) | p-value |
|---|---:|---:|---:|
| flat_torus2d | 7.75 | 0.583 | \(4.12\times10^{-83}\) |
| sphere | 11.51 | 0.752 | \(5.49\times10^{-165}\) |

Thus:

\[
\delta W_{\rm loc}
\longrightarrow
\delta\kappa(r)
\]

is numerically supported on controlled geometries.

---

## 4. Connection with Paper 7 and Paper 8

The Step E result is not isolated. It confirms the same coupling already observed in Paper 7 and Paper 8:

\[
\delta W_{\rm loc}
\longrightarrow
\delta\kappa.
\]

| Paper | Test | Framework | Signal |
|---|---|---|---|
| Paper 7 | direct curvature response | quantum MI graphs, \(N=16\) | \(\Delta\kappa_{\rm edge}=0.076\), positive fraction \(=100\%\) |
| Paper 8 | reconstructed \(T_{\mu\nu}^{\rm eff}\) | source from \(\delta W_{\rm loc}\) | Spearman \(\rho=0.741\) |
| Paper 15 Step E | radial source response | controlled geometries | Spearman \(\rho=0.752\), near/far \(=11.51\) on sphere |

The same signal appears from three independent angles:

\[
\text{local curvature response}
\quad
\leftrightarrow
\quad
\text{effective stress-energy reconstruction}
\quad
\leftrightarrow
\quad
\text{controlled source-response geometry}.
\]

Paper 15 is stronger than the previous tests in one sense: the geometry is controlled, and the flat/curved reference cases are known.

---

## 5. Relation to Paper 8 source construction

Paper 8 constructed an effective stress-energy tensor from local entanglement perturbations:

\[
T_{\mu\nu}^{\rm eff}
=
T_{\mu\nu}^{\rm matter}[\delta W^{\rm loc}]
+
T_{\mu\nu}^{\rm ent}[d_s].
\]

The best source candidate was:

\[
S_{\rm flux}
=
T_{00}
-
\frac{1}{2}T_{aa}
+
\frac{1}{2}T_{\rm grad}
+
F.
\]

Paper 8 found:

\[
\rho_{\rm Spearman}=0.741,
\qquad
p=1.84\times10^{-4},
\]

between the reconstructed source and curvature fluctuations.

It also showed that the matter-like sector is not conserved in isolation:

\[
\nabla^\mu T_{\mu\nu}^{\rm matter}\neq0.
\]

In BuP, this is interpreted as exchange with the entanglement background:

\[
J_\nu^{\rm exchange}
=
\nabla^\mu
\left[
G_{\rm eff}(d_s)T_{\mu\nu}^{\rm matter}
\right],
\]

with compensation:

\[
\nabla^\mu
\left[
G\,T_{\mu\nu}^{\rm ent}
\right]
=
-
J_\nu^{\rm exchange}.
\]

Thus Step E of Paper 15 is not an open guess: it is the controlled-geometries version of the Paper 8 source mechanism.

---

## 6. Current interpretation

The current result can be summarized as:

\[
\boxed{
\text{Paper 15 does not yet prove Einstein, but validates the three arrows required for the effective derivation.}
}
\]

The three validated arrows are:

\[
W_{ij}\to L_\epsilon\to\Delta_g,
\]

\[
\kappa_{ij}^{\rm OR}\to \text{Ricci signal},
\]

\[
\delta W_{\rm loc}\to\delta\kappa(r).
\]

Together, they support the continuum target:

\[
\frac{\delta S_{\rm BuP}}{\delta W_{ij}}=0
\quad
\xrightarrow[N\to\infty]{}
\quad
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

---

## 7. Folder structure

```text
papers/paper15_einstein_derivation/
  README.md
  paper15_einstein_derivation.tex

  scripts/
    paper15_spectral_convergence_v1.py
    paper15_spectral_convergence_v2.py
    paper15_ricci_convergence_v1.py
    paper15_ricci_convergence_v2.py
    paper15_source_response_v1.py
    paper15_source_response_v2.py
    paper15_build_numerical_summary_v1.py

  results/
    spectral_convergence_v1/
    spectral_convergence_v2/
    spectral_convergence_v2_best_unnormalized/
    ricci_convergence_v1/
    ricci_convergence_v2/
    source_response_v1/
    source_response_v2/
    paper15_numerical_summary_v1/

  figures/
    # final figures copied from selected result folders

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

# Paper 14 — Modular/Spectral Action and Galactic Phenomenology

## Status

This folder contains the numerical and phenomenological tests for Paper 14 of the Bottom-Up Quantum Gravity program.

Paper 14 connects three levels of the BuP framework:

1. the modular Hamiltonian \(K_A=-\log\rho_A\),
2. the entanglement Laplacian \(L_{\rm ent}\),
3. the effective gravitational response tested on SPARC rotation curves.

The central result is that the same spectral structure controls both the modular exponent \(\beta_{\rm mod}\) and the gravitational exponent \(\beta_{\rm grav}\), and that the resulting correlation length \(\lambda_{\rm corr}\) can be predicted on galaxies without fitting it individually.

---

## 1. Scientific motivation

The starting question of Paper 14 is:

\[
\text{Can modular quantum dynamics and emergent gravity be derived from the same entanglement spectrum?}
\]

In previous papers, BuP had already introduced:

\[
L_{\rm ent}
\]

as the Laplacian of the entanglement graph, and the effective gravitational exponent:

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
\]

Paper 14 asks whether the modular Hamiltonian,

\[
K_A=-\log\rho_A,
\]

can also be represented as a spectral function of the same entanglement Laplacian.

---

## 2. Core hypothesis

The working hypothesis is:

\[
K_A \sim f(L_{\rm ent}).
\]

More specifically, the tests compare the spectrum of \(K_A\) with spectral functions of candidate Laplacians:

\[
L_A^{\rm induced},
\qquad
L_A^{\rm Schur},
\qquad
L_A^{\rm MI-normalized}.
\]

The effective fit takes the form:

\[
K_A \sim (L_A)^{\beta_{\rm mod}}.
\]

The exponent \(\beta_{\rm mod}\) is then compared to the gravitational exponent:

\[
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
\]

---

## 3. Modular spectral test

The first step was to test whether \(K_A\) can be represented by a spectral function of the entanglement Laplacian.

The best \(N=16\) tests gave:

\[
R^2 \simeq 0.99.
\]

This shows that the modular Hamiltonian admits an effective spectral representation in terms of the entanglement graph.

The positive spectral family works well, while the negative family collapses to near-zero explanatory power.

---

## 4. Failure of the naive \(C\to\beta\) relation

A natural first hypothesis was that the modular topological constant

\[
C_{\rm modular}
=
d_A g_2^{\rm plateau}
\]

could directly predict \(\beta_{\rm mod}\).

This hypothesis failed.

Across the clean \(N=9\), \(N=16\), and optimal \(N=16\) datasets, the direct linear correlation between \(C_{\rm modular}\) and \(\beta_{\rm mod}\) remains weak:

\[
R^2 \approx 0.
\]

This failure is important. It shows that \(\beta_{\rm mod}\) is not controlled by a single global plateau constant. The relevant structure is multivariate and spectral.

---

## 5. Spectral features and SFF diagnostics

Paper 14 then extracted additional spectral-form-factor features:

\[
t_{\rm dip},
\qquad
t_{\rm ramp},
\qquad
{\rm slope}_{\rm ramp},
\qquad
\Delta_3.
\]

The feature \(\Delta_3\) is more informative than \(C_{\rm modular}\), but still insufficient by itself.

This led to a multivariate model using:

\[
C_{\rm modular},
\quad
K_{\rm gap},
\quad
\Delta_3,
\quad
\lambda_2,
\quad
d_s,
\quad
d_w,
\quad
d_s/d_w,
\quad
\langle I\rangle,
\quad
I_{\max}.
\]

---

## 6. Predicting \(\beta_{\rm mod}\)

The decisive improvement comes from adding the graph and diffusion invariants:

\[
d_s,
\qquad
d_w,
\qquad
\lambda_2,
\qquad
d_s/d_w.
\]

On the Schur/RMT subset, the best models reach:

\[
R^2 \simeq 0.96
\]

in k-fold validation, and remain strong in leave-one-regime-out tests.

The most important features are:

\[
\lambda_2^{\rm norm},
\qquad
d_s/d_w,
\qquad
\Delta_3,
\qquad
d_s.
\]

This establishes that:

\[
\beta_{\rm mod}
=
F(\lambda_2,d_s,d_w,\Delta_3,\ldots).
\]

---

## 7. Modular--gravitational bridge

Using the BuP gravitational relation:

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4,
\]

we define:

\[
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
\]

The measured relation is:

\[
\beta_{\rm mod}
\simeq
0.531
+
1.726\,\beta_{\rm grav}.
\]

This is the central bridge of Paper 14.

It means that modular dynamics and effective gravity are two spectral projections of the same entanglement Laplacian.

---

## 8. SPARC bridge test

The bridge was then tested on galaxy-scale entanglement graphs reconstructed from SPARC baryonic profiles.

The pipeline is:

\[
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s,d_w,\lambda_2,\Delta_3)
\rightarrow
\beta_{\rm grav},\beta_{\rm mod}.
\]

For each galaxy, the bridge compares:

\[
\beta_{\rm mod}^{\rm bridge}
\]

against:

\[
\beta_{\rm mod}^{\rm ML}.
\]

The bridge error is:

\[
\epsilon_{\rm bridge}
=
\frac{
|\beta_{\rm mod}^{\rm ML}-\beta_{\rm mod}^{\rm bridge}|
}{
|\beta_{\rm mod}^{\rm ML}|
}.
\]

---

## 9. Predictive \(\lambda_{\rm corr}\) test

The first SPARC tests scanned:

\[
f_\lambda
=
\lambda_{\rm corr}/R_d
\in
\{0.5,1.0,1.5,2.0,3.0\}.
\]

The decisive test was then performed without fitting \(\lambda_{\rm corr}\) galaxy by galaxy.

A leave-one-galaxy-out model was trained on 174 galaxies and used to predict:

\[
\widehat f_\lambda
\]

for the held-out galaxy.

The best model was a Random Forest classifier.

Result:

\[
173/175
\]

galaxies pass the bridge with predicted \(\lambda_{\rm corr}\), without individual adjustment.

The success rates are:

\[
98.86\%
\]

strong or moderate bridge, and:

\[
72.0\%
\]

strong bridge.

Thus:

\[
\lambda_{\rm corr}
\]

is not merely a fitted scale. It is predictable from the graph invariants.

---

## 10. Unified galaxy taxonomy

Paper 14 also builds a unified taxonomy combining:

1. the LOW/HIGH dimension phase,
2. the fit-quality category,
3. the optimal coherence class \(f_\lambda^{\rm opt}\).

The LOW/HIGH split is:

\[
N_{\rm HIGH}=88,
\qquad
N_{\rm LOW}=87.
\]

with:

\[
d_{\min}^{\rm HIGH}=2.487107,
\qquad
d_{\min}^{\rm LOW}=2.274448.
\]

The Paper 14 coherence classes are:

```text
short_0p5Rd          : 89 galaxies
standard_1Rd         : 38 galaxies
extended_1p5_2Rd     : 39 galaxies
very_extended_3Rd    : 9 galaxies

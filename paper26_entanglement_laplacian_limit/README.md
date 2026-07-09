# Paper 26 — Normalization and Continuum Limit of the Entanglement Laplacian

## Subtitle

From Mutual-Information Graphs to the Laplace--Beltrami Operator

---

## Purpose

Paper 26 attacks the central mathematical debt identified in Paper 25:

\[
c_N L_N \longrightarrow -\Delta_g.
\]

The goal is to determine the normalization under which the discrete BuP entanglement Laplacian converges to the continuum Laplace--Beltrami operator.

---

## Central object

BuP starts from the mutual-information graph

\[
W_{ij}=I(i:j),
\]

and defines the entanglement Laplacian

\[
L_{\rm ent}=D-W,
\qquad
D_{ii}=\sum_j W_{ij}.
\]

Paper 26 studies when and how a normalized version satisfies

\[
c_{N,\epsilon} L_{\rm ent}^{(N,\epsilon)} f
\to
-\Delta_g f.
\]

---

## Why this matters

This limit is the mathematical foundation behind the chain

\[
W_{ij}
\to
L_{\rm ent}
\to
-\Delta_g
\to
(-\Delta_g)^{-1}
\to
\Phi_{\rm BuP}
\to
\text{weak-field gravity}.
\]

If this convergence is controlled, then the BuP Green-function sector becomes much stronger:

\[
L_{\rm ent}^{+}
\to
(-\Delta_g)^{-1}.
\]

---

## Theorem target

Let \(M\) be a compact smooth Riemannian manifold and let \(x_i\in M\) be sampled from a density \(\rho\). Define a local kernel graph

\[
W_{ij}^{(\epsilon)}
=
\exp\left[
-\frac{d_g(x_i,x_j)^2}{4\epsilon}
\right].
\]

The target statement is:

\[
\lim_{N\to\infty,\epsilon\to0}
\left\|
c_{N,\epsilon} L_{N,\epsilon} f
+
\Delta_g f
\right\|_{L^2(M)}
=
0,
\]

for smooth test functions \(f\).

---

## BuP interpretation

For the physical BuP graph,

\[
W_{ij}=I(i:j),
\]

Paper 26 introduces the local-kernel hypothesis:

\[
I(i:j)
\simeq
F(d_g(x_i,x_j)),
\]

with a local regime in which \(F\) behaves like a diffusion kernel.

Thus, Paper 26 separates two layers:

1. mathematical layer: kernel graph Laplacian convergence;
2. physical layer: mutual-information graph as an effective local kernel.

---

## Numerical programme

The first numerical checks are performed on analytically controlled geometries:

| Geometry | Continuum spectrum |
|---|---|
| Circle \(S^1\) | \(\lambda_k=k^2\) |
| Flat torus \(T^2\) | \(\lambda_{m,n}=m^2+n^2\) up to geometric factors |
| Sphere \(S^2\) | \(\lambda_\ell=\ell(\ell+1)\) |

The first target is the circle because its spectrum is simple and exact.

---

## Expected outputs

```text
paper26_entanglement_laplacian_limit/
  README.md
  paper26_entanglement_laplacian_limit.tex

  scripts/
    paper26_build_circle_spectrum_v1.py
    paper26_scan_normalization_v1.py
    paper26_make_figures_v1.py

  results/
    paper26_laplacian_limit_v1/
      circle_spectrum_convergence.csv
      normalization_scan.csv
      summary.json

  figures/
    fig01_limit_chain.png
    fig02_circle_spectrum_convergence.png
    fig03_normalization_scaling.png

  notes/
    roadmap.md
    analytic_derivation.md
    referee_notes.md
    reproducibility.md
cat > paper26_entanglement_laplacian_limit/README.md <<'MD'
# Paper 26 — Normalization and Continuum Limit of the Entanglement Laplacian

## Subtitle

From Mutual-Information Graphs to the Laplace--Beltrami Operator

---

## Purpose

Paper 26 attacks the central mathematical debt identified in Paper 25:

\[
c_N L_N \longrightarrow -\Delta_g.
\]

The goal is to determine the normalization under which the discrete BuP entanglement Laplacian converges to the continuum Laplace--Beltrami operator.

---

## Central object

BuP starts from the mutual-information graph

\[
W_{ij}=I(i:j),
\]

and defines the entanglement Laplacian

\[
L_{\rm ent}=D-W,
\qquad
D_{ii}=\sum_j W_{ij}.
\]

Paper 26 studies when and how a normalized version satisfies

\[
c_{N,\epsilon} L_{\rm ent}^{(N,\epsilon)} f
\to
-\Delta_g f.
\]

---

## Why this matters

This limit is the mathematical foundation behind the chain

\[
W_{ij}
\to
L_{\rm ent}
\to
-\Delta_g
\to
(-\Delta_g)^{-1}
\to
\Phi_{\rm BuP}
\to
\text{weak-field gravity}.
\]

If this convergence is controlled, then the BuP Green-function sector becomes much stronger:

\[
L_{\rm ent}^{+}
\to
(-\Delta_g)^{-1}.
\]

---

## Theorem target

Let \(M\) be a compact smooth Riemannian manifold and let \(x_i\in M\) be sampled from a density \(\rho\). Define a local kernel graph

\[
W_{ij}^{(\epsilon)}
=
\exp\left[
-\frac{d_g(x_i,x_j)^2}{4\epsilon}
\right].
\]

The target statement is:

\[
\lim_{N\to\infty,\epsilon\to0}
\left\|
c_{N,\epsilon} L_{N,\epsilon} f
+
\Delta_g f
\right\|_{L^2(M)}
=
0,
\]

for smooth test functions \(f\).

---

## BuP interpretation

For the physical BuP graph,

\[
W_{ij}=I(i:j),
\]

Paper 26 introduces the local-kernel hypothesis:

\[
I(i:j)
\simeq
F(d_g(x_i,x_j)),
\]

with a local regime in which \(F\) behaves like a diffusion kernel.

Thus, Paper 26 separates two layers:

1. mathematical layer: kernel graph Laplacian convergence;
2. physical layer: mutual-information graph as an effective local kernel.

---

## Numerical programme

The first numerical checks are performed on analytically controlled geometries:

| Geometry | Continuum spectrum |
|---|---|
| Circle \(S^1\) | \(\lambda_k=k^2\) |
| Flat torus \(T^2\) | \(\lambda_{m,n}=m^2+n^2\) up to geometric factors |
| Sphere \(S^2\) | \(\lambda_\ell=\ell(\ell+1)\) |

The first target is the circle because its spectrum is simple and exact.

---

## Expected outputs

```text
paper26_entanglement_laplacian_limit/
  README.md
  paper26_entanglement_laplacian_limit.tex

  scripts/
    paper26_build_circle_spectrum_v1.py
    paper26_scan_normalization_v1.py
    paper26_make_figures_v1.py

  results/
    paper26_laplacian_limit_v1/
      circle_spectrum_convergence.csv
      normalization_scan.csv
      summary.json

  figures/
    fig01_limit_chain.png
    fig02_circle_spectrum_convergence.png
    fig03_normalization_scaling.png

  notes/
    roadmap.md
    analytic_derivation.md
    referee_notes.md
    reproducibility.md
Status

Paper 26 is a mathematical foundation paper.

It does not claim that the full mutual-information graph automatically converges to a geometric Laplacian. Instead, it aims to identify the exact assumptions under which the convergence holds and to isolate the additional BuP physical hypothesis required for W
ij
	​

=I(i:j).

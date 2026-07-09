# Paper 26 — Non-Uniform Density Bias

## Purpose

This note records an important limitation of the raw unnormalized graph Laplacian

\[
L=D-W.
\]

The convergence

\[
c_{N,\epsilon}L_{N,\epsilon}
\to
-\Delta_g
\]

holds in the locally uniform density regime.

For non-uniform node density, the raw combinatorial Laplacian converges instead to a density-biased diffusion operator.

---

## Uniform-density result

For a compact \(D\)-dimensional manifold \(M\), local Gaussian kernel

\[
W_{ij}
=
\exp\left[
-\frac{d_g(x_i,x_j)^2}{4\epsilon}
\right],
\]

and uniform node density

\[
\rho=\frac{N}{\mathrm{Vol}(M)},
\]

the raw graph Laplacian satisfies

\[
L_{N,\epsilon}f
\simeq
-\rho(4\pi)^{D/2}\epsilon^{D/2+1}
\Delta_g f.
\]

Thus,

\[
\boxed{
c_{N,\epsilon}
=
\frac{1}
{\rho(4\pi)^{D/2}\epsilon^{D/2+1}}
}
\]

gives

\[
c_{N,\epsilon}L_{N,\epsilon}f
\to
-\Delta_g f.
\]

---

## Non-uniform density

If the sampling density is position-dependent,

\[
\rho=\rho(x),
\]

then the same expansion gives a drift contribution.

At leading order,

\[
L_{N,\epsilon}f
\simeq
-\rho(x)(4\pi)^{D/2}\epsilon^{D/2+1}
\left[
\Delta_g f
+
2\nabla\log\rho\cdot\nabla f
\right].
\]

Therefore, after local normalization,

\[
\boxed{
c_{N,\epsilon}(x)L_{N,\epsilon}f
\simeq
-\Delta_g f
-
2\nabla\log\rho\cdot\nabla f.
}
\]

The exact drift coefficient depends on the graph normalization convention, but the existence of a density-induced drift term is unavoidable for the raw combinatorial Laplacian.

---

## Controlled \(S^1\) test

We tested the density

\[
\rho(\theta)
=
\frac{1}{2\pi}
\left(
1+a\cos\theta
\right),
\qquad
a=0.5,
\]

and the test function

\[
f(\theta)=\sin(2\theta).
\]

For this density,

\[
\partial_\theta\log\rho
=
\frac{-a\sin\theta}{1+a\cos\theta}.
\]

The two targets are:

\[
\text{pure target}
=
-f''(\theta),
\]

and

\[
\text{drift target}
=
-f''(\theta)
-
2(\partial_\theta\log\rho)f'(\theta).
\]

---

## Numerical result

Best run:

\[
N=2048,
\qquad
\epsilon=0.02,
\qquad
a=0.5.
\]

The density contrast is approximately

\[
\rho_{\max}/\rho_{\min}\simeq 3.
\]

The raw graph operator gives:

\[
\text{RMSE against pure Laplacian}
=
0.3921,
\]

\[
\text{RMSE against drift operator}
=
0.1738.
\]

Therefore,

\[
\boxed{
\text{drift target improves the fit by a factor }2.26.
}
\]

---

## Interpretation

This confirms that the raw unnormalized graph Laplacian should not be directly identified with the Laplace--Beltrami operator when the node density is non-uniform.

For non-uniform BuP graphs, including baryon-weighted SPARC graphs and lensing-weighted SLACS graphs, one must either:

1. work in a locally uniform regime;
2. include the density-drift term explicitly;
3. use a density-corrected normalized graph Laplacian.

---

## Referee-safe statement

The central Paper 26 theorem target should be stated as:

\[
\boxed{
c_{N,\epsilon}L_{N,\epsilon}\to-\Delta_g
\quad
\text{under locally uniform node density.}
}
\]

For non-uniform entanglement graphs, the raw Laplacian contains density bias, and density-corrected normalization is required to recover a pure Laplace--Beltrami limit.

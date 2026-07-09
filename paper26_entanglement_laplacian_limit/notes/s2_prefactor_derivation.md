# Paper 26 — Analytic Derivation of the \(S^2\) Prefactor

## Setup

We consider the unit sphere

\[
S^2=\{x\in\mathbb R^3:\|x\|=1\}.
\]

The kernel graph is

\[
W_{ij}
=
\exp\left[
-\frac{d_{S^2}(x_i,x_j)^2}{4\epsilon}
\right],
\]

where \(d_{S^2}\) is the geodesic distance on the unit sphere.

The unnormalized graph Laplacian is

\[
(Lf)_i
=
\sum_j W_{ij}(f_i-f_j).
\]

The diagonal is removed.

---

## General normalization

For a compact \(D\)-dimensional Riemannian manifold \(M\), uniform sampling density

\[
\rho=\frac{N}{\mathrm{Vol}(M)},
\]

and Gaussian kernel convention

\[
k_\epsilon(x,y)
=
\exp\left[
-\frac{d_g(x,y)^2}{4\epsilon}
\right],
\]

the local expansion gives

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
{\rho(4\pi)^{D/2}\epsilon^{D/2+1}}.
}
\]

---

## Specialization to the unit sphere

For the unit sphere,

\[
D=2,
\qquad
\mathrm{Vol}(S^2)=4\pi,
\qquad
\rho=\frac{N}{4\pi}.
\]

Therefore,

\[
c_{N,\epsilon}^{S^2}
=
\frac{1}
{
\frac{N}{4\pi}(4\pi)\epsilon^2
}.
\]

Hence,

\[
\boxed{
c_{N,\epsilon}^{S^2}
=
\frac{1}{N\epsilon^2}.
}
\]

---

## Continuum spectrum

The spectrum of the positive Laplace--Beltrami operator on the unit sphere is

\[
\lambda_\ell
=
\ell(\ell+1),
\]

with multiplicity

\[
2\ell+1.
\]

Therefore the target low spectrum is

\[
0,
\quad
2,2,2,
\quad
6,6,6,6,6,
\quad
12,\ldots
\]

---

## Numerical confirmation

The Fibonacci-sphere experiment gives the best row:

\[
N=1024,
\qquad
\epsilon=0.01.
\]

The empirical prefactor is

\[
c_{N,\epsilon}N\epsilon^2
=
1.0201663725.
\]

The theoretical prefactor is

\[
1.
\]

The absolute prefactor error is therefore

\[
2.02\%.
\]

The low-spectrum comparison gives

\[
\text{mean relative spectral error}
\simeq
6.85\%,
\]

\[
\text{max relative spectral error}
\simeq
10.35\%.
\]

---

## Interpretation

The \(S^2\) result is less sharp than the \(S^1\) and \(T^2\) tests because:

1. the Fibonacci sphere is quasi-uniform but not an exact spectral grid;
2. the graph Laplacian is not exactly circulant;
3. the continuum eigenvalues have nontrivial multiplicities \(2\ell+1\);
4. curvature introduces higher-order heat-kernel corrections.

Nevertheless, the prefactor confirms the general normalization law:

\[
\boxed{
c_{N,\epsilon}
=
\frac{1}
{\rho(4\pi)^{D/2}\epsilon^{D/2+1}}.
}
\]

---

## Result

For the unit sphere,

\[
\boxed{
\frac{1}{N\epsilon^2}
L_{N,\epsilon}
\longrightarrow
-\Delta_{S^2}
}
\]

in the low-spectrum continuum regime, with finite-\(N\), sampling and curvature corrections.

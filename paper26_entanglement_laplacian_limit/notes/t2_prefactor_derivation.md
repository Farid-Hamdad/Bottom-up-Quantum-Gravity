# Paper 26 — Analytic Derivation of the \(T^2\) Prefactor

## Setup

We consider the flat square torus

\[
T^2=[0,2\pi)^2
\]

with periodic coordinates

\[
x=(\theta_1,\theta_2).
\]

The kernel graph is

\[
W_{ij}
=
\exp\left[
-\frac{d_{T^2}(x_i,x_j)^2}{4\epsilon}
\right],
\]

and the unnormalized graph Laplacian is

\[
(Lf)_i
=
\sum_j W_{ij}(f_i-f_j).
\]

The diagonal is removed.

---

## Continuum approximation

For a uniform \(n_{\rm side}\times n_{\rm side}\) grid,

\[
N=n_{\rm side}^2,
\]

and

\[
\sum_j
\approx
\frac{N}{(2\pi)^2}
\int_{T^2}d^2y.
\]

For small \(\epsilon\), the kernel is local, so we use local coordinates

\[
u=y-x.
\]

Then

\[
(Lf)(x)
\approx
\frac{N}{(2\pi)^2}
\int_{\mathbb R^2}
e^{-|u|^2/(4\epsilon)}
\left[
f(x)-f(x+u)
\right]d^2u.
\]

Using Taylor expansion,

\[
f(x+u)
=
f(x)
+
u^a\partial_a f(x)
+
\frac12 u^a u^b\partial_a\partial_b f(x)
+
O(|u|^3).
\]

Thus,

\[
f(x)-f(x+u)
=
-u^a\partial_a f(x)
-
\frac12 u^a u^b\partial_a\partial_b f(x)
+
O(|u|^3).
\]

The odd term vanishes by symmetry:

\[
\int_{\mathbb R^2}
u^a e^{-|u|^2/(4\epsilon)}d^2u
=
0.
\]

By rotational symmetry,

\[
\int_{\mathbb R^2}
u^a u^b
e^{-|u|^2/(4\epsilon)}
d^2u
=
\delta^{ab}
\frac{1}{2}
\int_{\mathbb R^2}
|u|^2 e^{-|u|^2/(4\epsilon)}d^2u.
\]

Equivalently, using one-dimensional Gaussian moments,

\[
\int_{\mathbb R^2}
u_1^2
e^{-|u|^2/(4\epsilon)}
d^2u
=
8\pi\epsilon^2.
\]

Therefore,

\[
(Lf)(x)
\approx
-
\frac{N}{(2\pi)^2}
\frac12
\left(
8\pi\epsilon^2
\right)
\Delta f(x).
\]

Hence

\[
(Lf)(x)
\approx
-
\frac{N}{\pi}
\epsilon^2
\Delta f(x).
\]

---

## Normalization

To obtain

\[
c_{N,\epsilon}L_{N,\epsilon}
\to
-\Delta_{T^2},
\]

we need

\[
\boxed{
c_{N,\epsilon}^{T^2}
=
\frac{\pi}{N\epsilon^2}.
}
\]

---

## Numerical confirmation

The FFT torus experiment gives:

\[
c_{N,\epsilon}N\epsilon^2
\simeq
3.1573267968,
\]

while

\[
\pi
=
3.1415926536.
\]

The relative prefactor error is

\[
5.01\times10^{-3},
\]

or approximately \(0.50\%\).

The low-spectrum comparison gives:

\[
\text{mean relative spectral error}
\simeq
2.46\%,
\]

\[
\text{max relative spectral error}
\simeq
4.84\%.
\]

---

## Result

For the kernel convention

\[
W_{ij}
=
\exp\left[
-\frac{d_{T^2}(x_i,x_j)^2}{4\epsilon}
\right],
\]

and the unnormalized graph Laplacian

\[
L=D-W,
\]

one obtains

\[
\boxed{
\frac{\pi}{N\epsilon^2}
L_{N,\epsilon}
\longrightarrow
-\Delta_{T^2}.
}
\]

This confirms the \(D=2\) case of the general scaling law

\[
c_{N,\epsilon}^{(D)}
\propto
\frac{1}{N\epsilon^{D/2+1}}.
\]

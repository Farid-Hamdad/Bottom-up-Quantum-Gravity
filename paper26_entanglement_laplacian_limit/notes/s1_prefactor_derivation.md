# Paper 26 — Analytic Derivation of the \(S^1\) Prefactor

## Setup

We consider the unit circle \(S^1\) parameterized by

\[
\theta\in[0,2\pi).
\]

The kernel graph is defined by

\[
W_{ij}
=
\exp\left[
-\frac{d_{S^1}(\theta_i,\theta_j)^2}{4\epsilon}
\right],
\]

where

\[
d_{S^1}(\theta,\phi)
=
\min\left(
|\theta-\phi|,
2\pi-|\theta-\phi|
\right).
\]

The unnormalized graph Laplacian is

\[
(Lf)_i
=
\sum_j W_{ij}(f_i-f_j).
\]

The diagonal is removed.

---

## Continuum approximation

For uniformly sampled points on \(S^1\),

\[
\theta_j=\frac{2\pi j}{N},
\]

the sum is approximated by

\[
\sum_j
\approx
\frac{N}{2\pi}
\int_0^{2\pi}d\phi.
\]

For small \(\epsilon\), the kernel is local, so locally we set

\[
u=\phi-\theta.
\]

Then

\[
(Lf)(\theta)
\approx
\frac{N}{2\pi}
\int_{-\infty}^{+\infty}
e^{-u^2/(4\epsilon)}
\left[
f(\theta)-f(\theta+u)
\right]du.
\]

Using the Taylor expansion

\[
f(\theta+u)
=
f(\theta)
+
u f'(\theta)
+
\frac{u^2}{2}f''(\theta)
+
O(u^3),
\]

we obtain

\[
f(\theta)-f(\theta+u)
=
-u f'(\theta)
-
\frac{u^2}{2}f''(\theta)
+
O(u^3).
\]

The odd term vanishes by symmetry:

\[
\int_{-\infty}^{+\infty}
u e^{-u^2/(4\epsilon)}du
=
0.
\]

Thus,

\[
(Lf)(\theta)
\approx
-
\frac{N}{2\pi}
\frac{f''(\theta)}{2}
\int_{-\infty}^{+\infty}
u^2 e^{-u^2/(4\epsilon)}du.
\]

The Gaussian moment is

\[
\int_{-\infty}^{+\infty}
u^2 e^{-u^2/(4\epsilon)}du
=
4\sqrt{\pi}\epsilon^{3/2}.
\]

Therefore,

\[
(Lf)(\theta)
\approx
-
\frac{N}{2\pi}
\frac{1}{2}
4\sqrt{\pi}\epsilon^{3/2}
f''(\theta).
\]

Hence

\[
(Lf)(\theta)
\approx
-
\frac{N}{\sqrt{\pi}}
\epsilon^{3/2}
f''(\theta).
\]

---

## Normalization

Since on \(S^1\),

\[
-\Delta_{S^1}
=
-\frac{d^2}{d\theta^2},
\]

we need

\[
c_{N,\epsilon}L_{N,\epsilon}
\to
-\Delta_{S^1}.
\]

Thus,

\[
\boxed{
c_{N,\epsilon}^{S^1}
=
\frac{\sqrt{\pi}}{N\epsilon^{3/2}}.
}
\]

---

## Numerical confirmation

The numerical scan gives:

\[
c_{N,\epsilon}N\epsilon^{3/2}
\simeq
1.77467,
\]

while

\[
\sqrt{\pi}
=
1.7724538509.
\]

The best relative error is

\[
1.25\times10^{-3},
\]

or approximately \(0.125\%\).

The low-spectrum comparison gives:

\[
\text{mean relative spectral error}
\simeq
1.23\%,
\]

\[
\text{max relative spectral error}
\simeq
2.94\%.
\]

---

## Result

For the kernel convention

\[
W_{ij}
=
\exp\left[
-\frac{d_{S^1}(\theta_i,\theta_j)^2}{4\epsilon}
\right],
\]

and the unnormalized graph Laplacian

\[
L=D-W,
\]

one obtains

\[
\boxed{
\frac{\sqrt{\pi}}{N\epsilon^{3/2}}
L_{N,\epsilon}
\longrightarrow
-\Delta_{S^1}.
}
\]

This is the first controlled analytic normalization result of Paper 26.

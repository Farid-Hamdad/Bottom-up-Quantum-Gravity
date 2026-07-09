# Paper 26 — Analytic Derivation Notes

## Kernel operator

Consider

\[
K_\epsilon f(x)
=
\int_M
\exp\left[
-\frac{d_g(x,y)^2}{4\epsilon}
\right]
f(y)\rho(y)\,d{\rm vol}_g(y).
\]

For small \(\epsilon\), heat-kernel asymptotics suggest an expansion of the form

\[
K_\epsilon f(x)
=
A_\epsilon \rho(x) f(x)
+
B_\epsilon \epsilon
\Delta_g(\rho f)(x)
+
O(\epsilon^2).
\]

---

## Unnormalized graph Laplacian

The unnormalized graph Laplacian is

\[
L f_i
=
\sum_j W_{ij}(f_i-f_j).
\]

In the large-\(N\) limit,

\[
L f(x)
\approx
N
\int_M
k_\epsilon(x,y)
(f(x)-f(y))
\rho(y)d{\rm vol}_g(y).
\]

For uniform sampling \(\rho=\text{constant}\), one expects

\[
L f(x)
\approx
-C_D N \epsilon^{D/2+1}\Delta_g f(x).
\]

Therefore,

\[
c_{N,\epsilon}
\sim
\frac{1}{C_D N\epsilon^{D/2+1}}.
\]

---

## Expected scaling

For dimension \(D\),

\[
c_{N,\epsilon}
\propto
\frac{1}{N\epsilon^{D/2+1}}.
\]

For the circle \(D=1\),

\[
c_{N,\epsilon}
\propto
\frac{1}{N\epsilon^{3/2}}.
\]

For the torus or sphere \(D=2\),

\[
c_{N,\epsilon}
\propto
\frac{1}{N\epsilon^{2}}.
\]

The exact prefactor depends on:

- kernel convention,
- volume convention,
- graph normalization,
- sampling density,
- whether the diagonal is removed,
- whether the Laplacian is normalized or unnormalized.

---

## Working target

Paper 26 should determine the precise prefactor for the convention used in the numerical scripts.

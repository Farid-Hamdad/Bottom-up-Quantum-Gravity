# Paper 16 — Spectral Summary v1

This summary gathers the first eigenvalue-level convergence tests for Paper 16.

The target statement is:

\[
c_N L_N
=
c_N\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
-\Delta_g.
\]

The scalar \(c_N\) is currently fitted numerically. A main theoretical task of Paper 16 is to derive this normalization analytically.

---

## 1. Key results

| Test | Limit | \(N\) | Modes | Mean rel. error | Median rel. error | Max rel. error | \(\lambda_1^{scaled}\) | \(\lambda_1^{target}\) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Circle spectrum | `c_N L_N -> -Delta_{S^1}` | 1024 | 12 | 0.0051 | 0.0051 | 0.0107 | 1.0062 | 1.0000 |
| Flat torus spectrum | `c_N L_N -> -Delta_{T^2}` | 1024 | 20 | 0.0266 | 0.0263 | 0.0522 | 40.9613 | 39.4784 |

---

## 2. Circle spectrum

For the unit circle \(S^1\), the analytic spectrum is

\[
\lambda_m=m^2,
\qquad
m=0,1,1,2,2,3,3,\ldots
\]

At \(N=1024\), the first 12 nonzero modes are reconstructed with:

\[
\mathrm{mean\ relative\ error}=0.0051,
\qquad
\mathrm{median\ relative\ error}=0.0051,
\qquad
\mathrm{max\ relative\ error}=0.0107.
\]

The first nonzero eigenvalue is

\[
\lambda_1^{\rm scaled}=1.0062,
\qquad
\lambda_1^{S^1}=1.
\]

This gives eigenvalue-level support for

\[
c_N L_N\longrightarrow -\Delta_{S^1}.
\]

---

## 3. Flat torus spectrum

For the unit flat torus \(T^2=[0,1)^2\), the analytic spectrum is

\[
\lambda_{m,n}=4\pi^2(m^2+n^2),
\qquad
(m,n)\in\mathbb{Z}^2\setminus\{(0,0)\}.
\]

At \(N=1024\), the first 20 nonzero modes are reconstructed with:

\[
\mathrm{mean\ relative\ error}=0.0266,
\qquad
\mathrm{median\ relative\ error}=0.0263,
\qquad
\mathrm{max\ relative\ error}=0.0522.
\]

The first nonzero eigenvalue is

\[
\lambda_1^{\rm scaled}=40.9613,
\qquad
\lambda_1^{T^2}=39.4784.
\]

This gives eigenvalue-level support for

\[
c_N L_N\longrightarrow -\Delta_{T^2}.
\]

---

## 4. Interpretation

Paper 15 showed heat-trace dimension convergence. Paper 16 now strengthens this by testing individual eigenvalues.

The current results show:

- \(S^1\): low-spectrum convergence with mean relative error \(0.0051\).
- \(T^2\): low-spectrum convergence with mean relative error \(0.0266\).

This supports the spectral continuum limit:

\[
L_N\to \Delta_g
\]

at the level of low eigenmodes, up to a scalar normalization \(c_N\).

---

## 5. Remaining theoretical issue

The current tests fit \(c_N\) numerically. The next task is to derive \(c_N\) analytically from the kernel normalization, the sampling density, and the choice of graph Laplacian.

A target theorem should specify conditions under which

\[
c_N\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
-\Delta_g
\]

in a controlled convergence mode: pointwise on smooth test functions, quadratic-form convergence, heat-kernel convergence, or spectral convergence.

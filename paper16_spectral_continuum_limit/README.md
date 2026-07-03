# Paper 16 — Spectral Continuum Limit

**From the BuP entanglement Laplacian to the Laplace--Beltrami operator**

---

## Status

Paper 16 studies the first analytic arrow required for the continuum Einstein limit of BuP:

\[
L_N
=
\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
\Delta_g.
\]

Paper 15 provided numerical evidence at the heat-trace level:

\[
d_s(t)
=
-2\frac{d\log{\rm Tr}(e^{-tL})}{d\log t}.
\]

Paper 16 strengthens this by testing the convergence of individual eigenvalues:

\[
\lambda_k(L_N)
\longrightarrow
\lambda_k(\Delta_g).
\]

The current results support the low-spectrum convergence on two controlled geometries:

\[
S^1,
\qquad
T^2.
\]

---

## 1. Starting point

The fundamental object is the mutual-information matrix:

\[
W_{ij}=I(i:j).
\]

From it, one constructs the graph degree matrix:

\[
D_{ii}=\sum_j W_{ij},
\]

and the unscaled entanglement Laplacian:

\[
L_N^0=D_N-W_N.
\]

The continuum-scaled operator used in Paper 15 and Paper 16 is:

\[
L_N
=
\frac{D_N-W_N}{\epsilon_N}.
\]

The target continuum operator is the Laplace--Beltrami operator:

\[
\Delta_g
=
\frac{1}{\sqrt{|g|}}
\partial_\mu
\left(
\sqrt{|g|}g^{\mu\nu}\partial_\nu
\right).
\]

Because the raw graph Laplacian still has an unknown multiplicative normalization, the spectral tests currently fit one scalar factor \(c_N\):

\[
c_N L_N
\longrightarrow
-\Delta_g.
\]

One of the main theoretical tasks of Paper 16 is to derive \(c_N\) analytically rather than fitting it.

---

## 2. Link with Paper 15

Paper 15 validated the heat-trace spectral dimension using the rescaled Laplacian:

\[
L_\epsilon=\frac{D-W}{\epsilon}.
\]

The main Paper 15 results were:

| Geometry | Target dimension | Measured \(d_s\) | Error |
|---|---:|---:|---:|
| circle | 1 | 1.0021 | 0.0021 |
| interval | 1 | 0.9445 | 0.0555 |
| grid2d | 2 | 1.9938 | 0.0062 |
| sphere | 2 | 1.8964 | 0.1036 |

This supported:

\[
L_N\to\Delta_g
\]

at the heat-trace / dimension level.

Paper 16 now tests the stronger statement:

\[
\lambda_k(c_NL_N)
\simeq
\lambda_k(-\Delta_g).
\]

---

## 3. Main numerical results

Paper 16 currently has two positive eigenvalue-level tests.

| Test | Limit | \(N\) | Modes compared | Mean rel. error | Median rel. error | Max rel. error | \(\lambda_1^{\rm scaled}\) | \(\lambda_1^{\rm target}\) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Circle spectrum | \(c_NL_N\to-\Delta_{S^1}\) | 1024 | 12 | 0.0051 | 0.0051 | 0.0107 | 1.0062 | 1.0000 |
| Flat torus spectrum | \(c_NL_N\to-\Delta_{T^2}\) | 1024 | 20 | 0.0266 | 0.0263 | 0.0522 | 40.9613 | 39.4784 |

Thus, Paper 16 strengthens Paper 15:

\[
d_s\to D
\quad
\Longrightarrow
\quad
\lambda_k(L_N)\to\lambda_k(\Delta_g).
\]

---

## 4. Test 1 — Circle spectrum

For the unit circle \(S^1\), the analytic positive spectrum of the Laplace--Beltrami operator is:

\[
\lambda_m=m^2,
\]

with sine/cosine degeneracy:

\[
1,1,4,4,9,9,16,16,\ldots
\]

The graph is built from \(N\) equally spaced points on the circle, using the intrinsic geodesic distance:

\[
d(\theta_i,\theta_j)
=
\min(|\theta_i-\theta_j|,2\pi-|\theta_i-\theta_j|).
\]

The weights are:

\[
W_{ij}
=
\exp\left[
-\frac{d(\theta_i,\theta_j)^2}{4\epsilon_N}
\right]
\]

on a local \(k\)-NN graph, and the operator is:

\[
L_N=\frac{D_N-W_N}{\epsilon_N}.
\]

A scalar normalization \(c_N\) is fitted so that:

\[
c_N\lambda_k^{\rm graph}
\simeq
\lambda_k^{S^1}.
\]

The errors decrease systematically with \(N\):

| \(N\) | Mean rel. error | Median rel. error | Max rel. error |
|---:|---:|---:|---:|
| 64 | 0.0937 | 0.0952 | 0.1852 |
| 128 | 0.0492 | 0.0480 | 0.1003 |
| 256 | 0.0214 | 0.0215 | 0.0446 |
| 512 | 0.0109 | 0.0094 | 0.0249 |
| 1024 | 0.0051 | 0.0051 | 0.0107 |

At \(N=1024\):

\[
\lambda_1^{\rm scaled}=1.0062,
\]

while the analytic value is:

\[
\lambda_1^{S^1}=1.
\]

This provides eigenvalue-level evidence for:

\[
c_NL_N\longrightarrow-\Delta_{S^1}.
\]

---

## 5. Test 2 — Flat torus spectrum

For the unit flat torus:

\[
T^2=[0,1)^2,
\]

with periodic boundary conditions, the analytic spectrum is:

\[
\lambda_{m,n}
=
4\pi^2(m^2+n^2),
\qquad
(m,n)\in\mathbb{Z}^2\setminus\{(0,0)\}.
\]

The graph is built from a periodic \(m\times m\) grid using the intrinsic periodic distance.

Again, the operator is:

\[
L_N=\frac{D_N-W_N}{\epsilon_N},
\]

and a scalar \(c_N\) is fitted:

\[
c_N\lambda_k^{\rm graph}
\simeq
\lambda_k^{T^2}.
\]

The errors decrease with resolution:

| \(N\) | Mean rel. error | Median rel. error | Max rel. error |
|---:|---:|---:|---:|
| 64 | 0.0951 | 0.1340 | 0.1467 |
| 144 | 0.0757 | 0.0437 | 0.1367 |
| 256 | 0.0604 | 0.0574 | 0.1850 |
| 576 | 0.0334 | 0.0455 | 0.0504 |
| 1024 | 0.0266 | 0.0263 | 0.0522 |

At \(N=1024\):

\[
\lambda_1^{\rm scaled}=40.9613,
\]

while the analytic value is:

\[
\lambda_1^{T^2}=4\pi^2=39.4784.
\]

This provides eigenvalue-level evidence for:

\[
c_NL_N\longrightarrow-\Delta_{T^2}.
\]

---

## 6. Interpretation

The Paper 16 results show that the BuP rescaled graph Laplacian does not only reproduce the global heat-trace dimension. It also reconstructs the low spectrum of the corresponding Laplace--Beltrami operator.

Current evidence:

\[
S^1:
\quad
{\rm mean\ relative\ error}=0.0051,
\]

\[
T^2:
\quad
{\rm mean\ relative\ error}=0.0266.
\]

Thus, the spectral continuum limit is supported in dimension 1 and dimension 2:

\[
c_N\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
-\Delta_g.
\]

The remaining theoretical problem is to identify the correct normalization \(c_N\) and prove convergence under controlled assumptions.

---

## 7. Target theorem

A possible target theorem for Paper 16 is:

> Let \(M\) be a compact \(D\)-dimensional Riemannian manifold sampled by points \(x_i\). Suppose the BuP weights satisfy locally:
>
> \[
> W_{ij}
> =
> \exp\left[
> -\frac{d_g(x_i,x_j)^2}{4\epsilon_N}
> \right]
> +
> o(1).
> \]
>
> If:
>
> \[
> \epsilon_N\to0,
> \qquad
> N\epsilon_N^{D/2}\to\infty,
> \]
>
> then there exists a normalization \(c_N\) such that:
>
> \[
> c_N\frac{D_N-W_N}{\epsilon_N}
> \longrightarrow
> -\Delta_g
> \]
>
> in a controlled spectral or heat-kernel sense.

The exact convergence mode remains to be chosen.

---

## 8. Possible convergence modes

Paper 16 may use one or several of the following:

1. pointwise convergence on smooth test functions;
2. convergence of quadratic forms;
3. Mosco convergence;
4. eigenvalue convergence;
5. heat-kernel convergence;
6. heat-trace convergence;
7. convergence of spectral dimension.

---

## 9. Open theoretical questions

The main open questions are:

1. What is the analytic expression for \(c_N\)?
2. How does \(c_N\) depend on \(D\), \(\epsilon_N\), \(k_N\), and sampling density?
3. Does \(W_{ij}=I(i:j)\) naturally approximate a heat kernel?
4. How is the entanglement distance \(d_{\rm ent}\) related to \(d_g\)?
5. How should non-uniform sampling be handled?
6. How should boundaries be treated?
7. Can the proof be extended from synthetic manifolds to genuine quantum mutual-information graphs?

---

## 10. Folder structure

```text
papers/paper16_spectral_continuum_limit/
  README.md
  paper16_spectral_continuum_limit.tex

  scripts/
    paper16_circle_spectrum_convergence_v1.py
    paper16_flat_torus_spectrum_convergence_v1.py
    paper16_build_spectral_summary_v1.py

  results/
    circle_spectrum_convergence_v1/
    flat_torus_spectrum_convergence_v1/
    paper16_spectral_summary_v1/

  figures/
    # final figures copied from selected result folders

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

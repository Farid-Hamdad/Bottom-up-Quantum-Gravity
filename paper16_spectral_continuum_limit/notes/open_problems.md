# Paper 16 — Open Problems

## 1. Normalization \(c_N\)

The current spectral tests fit a scalar normalization:

\[
c_N\lambda_k^{graph}\simeq\lambda_k^{target}.
\]

The main theoretical open problem is to derive \(c_N\) analytically.

Expected dependencies:

\[
c_N=C(D,\rho_N,\epsilon_N,k_N,\gamma).
\]

## 2. Convergence mode

Possible convergence modes:

1. pointwise convergence on smooth functions;
2. convergence of quadratic forms;
3. Mosco convergence;
4. eigenvalue convergence;
5. heat-kernel convergence;
6. heat-trace convergence.

## 3. Mutual-information origin of \(W_{ij}\)

The synthetic tests use Gaussian weights:

\[
W_{ij}=\exp\left[-\frac{d_g(i,j)^2}{4\epsilon_N}\right].
\]

For true BuP, one needs:

\[
W_{ij}=I(i:j).
\]

Open question:

under what conditions does mutual information approximate a local heat kernel?

## 4. Non-uniform sampling

Real entanglement graphs may not correspond to uniform sampling.

Need to study:

- density correction;
- alpha-normalized diffusion maps;
- renormalized graph Laplacians.

## 5. Boundaries

Current tests use compact geometries without boundary:

\[
S^1,\quad T^2.
\]

Boundary conditions remain open.

## 6. High modes

Current convergence is low-spectrum.

High modes are more sensitive to:

- graph discreteness;
- kernel cutoff;
- kNN anisotropy;
- numerical noise.

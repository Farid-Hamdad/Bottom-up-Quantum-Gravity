# Paper 17 — Open Problems

## 1. Pointwise convergence

The current result is mean-level:

\[
\left\langle
\frac{\kappa^{OR}}{\epsilon}
\right\rangle
\simeq
B_N+C_N
\left\langle
R_{\mu\nu}u^\mu u^\nu
\right\rangle.
\]

The open problem is the stronger statement:

\[
\frac{\kappa_{ij}^{OR}}{\epsilon}
\to
C R_{\mu\nu}u^\mu u^\nu
\]

edge by edge or locally after averaging.

## 2. Flat torus dispersion

The flat torus mean is calibrated to zero, but the distribution remains broad.

Need to understand whether this comes from:

- kNN anisotropy;
- finite-size effects;
- edge-length mixing;
- Wasserstein numerical noise;
- graph topology irregularity.

## 3. Affine constants

The current constants are empirical:

\[
B_N\simeq -0.301,
\qquad
C_N\simeq 0.392.
\]

Need analytic derivation of:

\[
B_N,
\qquad
C_N.
\]

## 4. Graph construction

Current tests use \(k\)-NN graphs.

Need to test:

- fixed-radius graphs;
- heat-kernel full graphs with cutoff;
- diffusion maps normalization;
- density-corrected measures.

## 5. Other curvatures

Need to compare Ollivier--Ricci with:

- Forman--Ricci curvature;
- Bakry--Émery curvature;
- curvature from Laplacian heat kernel;
- Regge-like curvature proxies.

## 6. Other geometries

Need to test:

- spheres of radius \(r\);
- hyperbolic surfaces;
- product manifolds;
- non-constant curvature manifolds.

## 7. True MI graphs

Current tests use synthetic geometric kernels.

The BuP target is:

\[
W_{ij}=I(i:j).
\]

Need to test whether true mutual-information graphs exhibit the same Ricci calibration.

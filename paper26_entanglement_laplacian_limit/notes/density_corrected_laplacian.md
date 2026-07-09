# Paper 26 — Density-Corrected Laplacian Test

## Purpose

This note records the controlled non-uniform \(S^1\) test comparing the raw graph Laplacian with density-corrected kernels.

The goal is to verify that the density bias of the raw combinatorial Laplacian is structural and can be reduced by an appropriate graph normalization.

---

## Raw graph

The raw kernel is

\[
W_{ij}
=
\exp\left[
-\frac{d_{S^1}(\theta_i,\theta_j)^2}{4\epsilon}
\right].
\]

The raw graph Laplacian is

\[
L=D-W.
\]

For non-uniform density

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

the raw operator is expected to contain a density-drift term.

---

## Density correction

Define

\[
q_i=\sum_j W_{ij}.
\]

The corrected kernel is

\[
W_{ij}^{(\alpha)}
=
\frac{W_{ij}}
{q_i^\alpha q_j^\alpha}.
\]

The tested values were

\[
\alpha=0,
\qquad
\alpha=\frac12,
\qquad
\alpha=1.
\]

The case \(\alpha=0\) is the raw graph.

The case \(\alpha=1/2\) is a symmetric degree-flattening correction:

\[
W_{ij}^{(1/2)}
=
\frac{W_{ij}}{\sqrt{q_iq_j}}.
\]

---

## Test function

The test function is

\[
f(\theta)=\sin(2\theta).
\]

The pure Laplacian target is

\[
-f''(\theta).
\]

The density-drift target is

\[
-f''(\theta)
-
2(\partial_\theta\log\rho)f'(\theta).
\]

---

## Results

Median relative RMSE by \(\alpha\):

| \(\alpha\) | pure RMSE | drift RMSE | degree contrast |
|---:|---:|---:|---:|
| 0 | 0.4588 | 0.3803 | 3.3045 |
| 0.5 | 0.2690 | 0.2916 | 1.0622 |
| 1 | 0.3931 | 0.5163 | 3.0197 |

Best pure Laplacian recovery:

\[
\alpha=0.5,
\qquad
\epsilon=0.04,
\qquad
{\rm RMSE}_{pure}=0.16098.
\]

The degree contrast is reduced to

\[
1.0537,
\]

showing that the symmetric correction nearly flattens the graph degree distribution.

---

## Interpretation

The raw graph Laplacian is closer to the density-drift operator than to the pure Laplacian.

The symmetric correction \(\alpha=1/2\) reduces this density bias and gives the best recovery of the pure Laplace--Beltrami operator in this test.

The \(\alpha=1\) correction is not optimal for this particular symmetric combinatorial Laplacian construction. This is not a contradiction with diffusion-map theory: the optimal correction depends on whether one uses a raw combinatorial, random-walk, or symmetric normalized graph operator.

---

## Referee-safe conclusion

The Paper 26 theorem target must distinguish three regimes:

1. uniform density:
   \[
   c_{N,\epsilon}L_{N,\epsilon}\to-\Delta_g;
   \]

2. non-uniform density with raw combinatorial Laplacian:
   \[
   c_{N,\epsilon}(x)L_{N,\epsilon}
   \to
   -\Delta_g
   -
   2\nabla\log\rho\cdot\nabla;
   \]

3. density-corrected graph:
   a suitable normalization can reduce the drift and recover the pure Laplace--Beltrami operator more accurately.

The non-uniform \(S^1\) tests show that the density bias is measurable and that the symmetric correction \(W_{ij}/\sqrt{q_iq_j}\) is effective in this controlled setting.

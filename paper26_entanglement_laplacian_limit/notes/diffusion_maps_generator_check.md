# Paper 26 — Diffusion Maps Generator Check

## Purpose

This note resolves the apparent discrepancy observed in the density-corrected
combinatorial Laplacian test.

The previous finite-\(N\) symmetric combinatorial test found that

\[
W_{ij}^{(1/2)}
=
\frac{W_{ij}}{\sqrt{q_iq_j}}
\]

gave the best recovery of the pure Laplacian among the tested combinatorial
operators

\[
L^{(\alpha)}=D^{(\alpha)}-W^{(\alpha)}.
\]

This does not contradict Coifman--Lafon diffusion maps, because the standard
density-corrected diffusion-maps construction uses a row-normalized Markov
operator, not the raw symmetric combinatorial Laplacian.

---

## Diffusion maps convention

Define

\[
q_i=\sum_j W_{ij}.
\]

Then

\[
K_{ij}^{(\alpha)}
=
\frac{W_{ij}}
{q_i^\alpha q_j^\alpha}.
\]

The Markov matrix is

\[
P_{ij}^{(\alpha)}
=
\frac{K_{ij}^{(\alpha)}}
{\sum_j K_{ij}^{(\alpha)}}.
\]

The associated generator is

\[
G^{(\alpha)}
=
\frac{I-P^{(\alpha)}}{\epsilon}.
\]

This is the convention in which \(\alpha=1\) is expected to remove sampling
density bias and recover the pure Laplace--Beltrami operator.

---

## Controlled non-uniform \(S^1\) test

The density is

\[
\rho(\theta)
=
\frac{1}{2\pi}
\left(
1+a\cos\theta
\right),
\qquad
a=0.5.
\]

The test function is

\[
f(\theta)=\sin(2\theta).
\]

Two targets are compared:

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

## Results

Best pure Laplacian recovery:

\[
\alpha=1,
\qquad
\epsilon=0.08,
\qquad
{\rm RMSE}_{pure}=0.0441.
\]

The corresponding drift error is

\[
{\rm RMSE}_{drift}=0.3161.
\]

Median errors by \(\alpha\):

| \(\alpha\) | pure RMSE | drift RMSE | pure/drift ratio |
|---:|---:|---:|---:|
| 0 | 0.3757 | 0.2521 | 1.4077 |
| 0.5 | 0.2343 | 0.2858 | 0.9182 |
| 1 | 0.1216 | 0.3680 | 0.3492 |

---

## Interpretation

The apparent preference for \(\alpha=1/2\) in the previous test was specific to
the symmetric combinatorial operator

\[
D^{(\alpha)}-W^{(\alpha)}.
\]

When the proper diffusion-maps Markov generator

\[
G^{(\alpha)}
=
\frac{I-P^{(\alpha)}}{\epsilon}
\]

is tested, the expected density-corrected choice

\[
\alpha=1
\]

indeed gives the best recovery of the pure Laplace--Beltrami operator.

---

## Referee-safe conclusion

There is no contradiction with the diffusion-maps theory.

Paper 26 should explicitly distinguish:

1. the raw combinatorial Laplacian \(D-W\), which exhibits density drift;
2. the symmetric corrected combinatorial Laplacian, where \(\alpha=1/2\) works
   well as a finite-\(N\) degree-flattening correction;
3. the diffusion-maps Markov generator, where \(\alpha=1\) recovers the pure
   Laplace--Beltrami operator most accurately.

The final theorem target for density-corrected convergence should be stated for
the Markov generator convention, not for the raw combinatorial operator.

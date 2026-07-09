# Paper 26 — Referee Notes

## Main claim

Paper 26 should not claim that all mutual-information graphs automatically converge to Laplace--Beltrami operators.

The safe claim is:

\[
\boxed{
\text{Under a local diffusion-kernel hypothesis, the BuP entanglement Laplacian admits a controlled continuum limit.}
}
\]

---

## Expected objection 1

“The graph is built from mutual information, not from a Gaussian kernel.”

### Response

Correct. Paper 26 separates:

1. the mathematical kernel-laplacian theorem;
2. the physical BuP hypothesis that local mutual information behaves as a decaying geometric kernel.

---

## Expected objection 2

“The normalization \(c_{N,\epsilon}\) is convention-dependent.”

### Response

Correct. Paper 26 explicitly states the convention used and derives the corresponding scaling.

---

## Expected objection 3

“Unnormalized graph Laplacians suffer density bias.”

### Response

Correct. Paper 26 first treats uniform sampling, then records density correction as an open or secondary problem.

---

## Expected objection 4

“Spectral convergence is weaker than operator convergence.”

### Response

Correct. The first numerical evidence is spectral. The theorem target concerns operator convergence on smooth test functions.

---

## Referee-safe wording

Use:

- continuum-limit target;
- controlled kernel regime;
- local mutual-information hypothesis;
- convention-dependent normalization;
- spectral evidence;
- operator-level theorem target.

Avoid:

- complete proof for arbitrary \(I(i:j)\);
- universal convergence;
- exact gravity derivation from this paper alone.

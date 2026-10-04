# Paper 30 — Causal Compensation Landscape

## IMPLEMENTATION CLARIFICATION — REFLECTION CONTROL

### 1. Status

This clarification is frozen before execution of any previously unknown
landscape cell.

It does not change the preregistered grid, observables, thresholds,
support criterion, or interpretation boundary.

### 2. Local ordering under spatial reflection

For a contiguous four-qubit subsystem

A_a = (a, a+1, a+2, a+3),

the spatially mirrored subsystem is

A_m = (m, m+1, m+2, m+3),

with

m = N - 4 - a.

Spatial reflection reverses the internal local ordering.

Therefore the state of A_m must be transformed by the four-qubit
reversal permutation before matrix-element comparison with A_a.

In local indices:

0 <-> 3
1 <-> 2.

The right-subsystem reduced density matrix is compared after the
corresponding ket and bra permutation.

### 3. Pair-observable mapping

Under the same local reflection, the six unordered local pairs map as:

01 <-> 23
02 <-> 13
03 <-> 03
12 <-> 12.

Accordingly, mirrored W_ij, modular A_ij, and v_ij values are compared
after applying this pair permutation.

Scalar permutation-invariant quantities such as M_W and M_mod are
compared directly.

The five Pearson residual correlations, alpha_optimal values, R_C
values, and support booleans are also compared directly because the
transition norms are invariant under the same fixed pair permutation.

### 4. Frozen tolerance

The preregistered reflection tolerance remains unchanged:

1e-10.

Historical support booleans must agree exactly.

### 5. Interpretation

This clarification specifies the mathematically correct coordinate map
for the already-preregistered statement that mirrored subsystems should
agree.

It introduces no new physical hypothesis and uses no unknown landscape
outcome.

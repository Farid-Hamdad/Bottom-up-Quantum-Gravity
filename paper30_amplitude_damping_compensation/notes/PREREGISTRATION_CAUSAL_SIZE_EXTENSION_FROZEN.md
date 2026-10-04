# Paper 30 — Causal-size extension preregistration

## Status

Prospective extension motivated after completion of the frozen Paper 30 analysis.

The original Paper 30 preregistration, implementation, results, and classification remain unchanged.

This extension tests whether the previously observed near-identity of N=10 and N=12 local compensation observables is explained by finite causal reach and boundary placement.

## 1. Scientific question

Does the Paper 30 residual compensation structure remain present once the observed subsystem A becomes genuinely sensitive to global system size?

The extension separates two possible causes of the original N10/N12 numerical equivalence:

1. A was located at the boundary of the open chain.
2. The finite Trotter depth P=6 limited the causal region influencing A.

## 2. Frozen quantities inherited unchanged from Paper 30

- TFIM parameters: J=1, h=1.
- Strang time step: dt=0.35.
- Initial state: |+>^N.
- Time-reversal mixture construction.
- Independent local amplitude-damping channel.
- s grid: 0.0 to 3.0 in steps of 0.1.
- gamma(s)=1-exp(-s).
- A size: 4 qubits.
- Pair ordering: (01,02,03,12,13,23).
- Mutual-information sector W.
- Normalized modular-Hamiltonian sector.
- Transition definitions d_W and d_mod.
- LOO local-linear Gaussian-kernel detrending against gamma.
- Bandwidth fractions: 0.05, 0.10, 0.20, 0.30, 0.40.
- Alpha window: 0.8 <= alpha <= 1.2.
- Compensation threshold: R_C <= 0.20.
- All original numerical safeguards.

No compensation metric or threshold may be modified after unblinding.

## 3. Systems

Only N=10 and N=12 are used in this extension.

They are compared under four frozen configurations.

### C0 — original control
- P=6
- N10: A=(0,1,2,3)
- N12: A=(0,1,2,3)

Purpose: reproduce the original Paper 30 local equivalence.

### C1 — centered subsystem
- P=6
- N10: A=(3,4,5,6)
- N12: A=(4,5,6,7)

Purpose: expose A symmetrically to both directions and test boundary-placement effects.

### C2 — enlarged causal depth
- P=8
- N10: A=(0,1,2,3)
- N12: A=(0,1,2,3)

Purpose: enlarge the causal reach while retaining the original edge geometry.

### C3 — centered and enlarged
- P=8
- N10: A=(3,4,5,6)
- N12: A=(4,5,6,7)

Purpose: combined test with both original restrictions removed.

## 4. Primary size-sensitivity observable

For each configuration and each gamma, compute

D_A^(10,12)(gamma)
    = 1/2 || rho_A^(N10)(gamma) - rho_A^(N12)(gamma) ||_1.

The two reduced states are compared in their common four-qubit local ordering.

Define

D_A,max = max_gamma D_A^(10,12)(gamma).

Numerical equivalence:
D_A,max <= 1e-10.

Size sensitivity established:
D_A,max > 1e-10.

The threshold is frozen before execution.

## 5. Secondary size-sensitivity diagnostics

For each gamma also record absolute N10/N12 differences in:

- all six W_ij,
- all six modular coefficients A_ij,
- all six v_ij,
- M_W,
- M_mod,
- normalized W direction,
- normalized modular direction.

These are descriptive diagnostics and do not replace D_A,max as the primary size-sensitivity test.

## 6. Compensation analysis

For each N and each configuration independently, run the unchanged Paper 30 compensation pipeline.

At all five frozen bandwidths evaluate:

P1: Pearson(resid_dW,resid_dmod) < 0.
P2: 0.8 <= alpha <= 1.2.
P3: R_C <= 0.20.

A configuration is compensation-supported for a given N only if P1, P2, and P3 pass at all five bandwidths.

## 7. Main prospective interpretation

The key question is conditional:

If size sensitivity is established in C1, C2, or C3, does the compensation signature remain supported for both N=10 and N=12?

Classification:

- CAUSAL_SIZE_EXPOSURE_NOT_ESTABLISHED:
  none of C1/C2/C3 has D_A,max > 1e-10.

- COMPENSATION_ROBUST_AFTER_CAUSAL_SIZE_EXPOSURE:
  at least one of C1/C2/C3 has D_A,max > 1e-10 and, for every size-sensitive configuration, both N10 and N12 satisfy all frozen Paper 30 compensation endpoints.

- COMPENSATION_NOT_ROBUST_AFTER_CAUSAL_SIZE_EXPOSURE:
  size sensitivity is established, but at least one size-sensitive configuration fails the frozen compensation criteria for N10 or N12.

No additional category may be introduced after unblinding.

## 8. Control requirement

C0 must reproduce the original Paper 30 N10/N12 local observables and compensation results within numerical tolerance.

Failure of this recovery control invalidates the extension run until the implementation discrepancy is resolved.

## 9. Interpretation boundary

This extension tests finite-size sensitivity, boundary placement, causal reach, and robustness of the previously observed information-sector compensation.

It does not establish:

- a fundamental conservation law,
- literal transfer of information between sectors,
- thermodynamic-limit universality,
- a universal alpha constant,
- relativistic signal propagation,
- equality with the gravitational-wave speed of Paper 22,
- or a fundamental spacetime causal cone.

Any additional structure discovered after execution is exploratory unless separately preregistered.

# Paper 30 — Causal Compensation Landscape

## IMPLEMENTATION AMENDMENT — FROZEN BEFORE LANDSCAPE UNBLINDING

### 1. Status

This amendment is frozen after the landscape preregistration and before any new landscape cell has been executed.

The previously known C0-C3 configurations were used only as recovery anchors, as explicitly allowed by the frozen preregistration.

No landscape result for any new (N,P,a) cell has been inspected.

### 2. Reason for amendment

The frozen preregistration specified an algebraically reduced implementation in which amplitude damping is applied directly to the four-qubit reduced state rho_A.

That implementation was validated at the density-matrix level to approximately 10^-16, but the frozen C0-C3 anchor recovery requirement is stricter: all archived local compensation quantities must reproduce within absolute tolerance 1e-10.

Direct local damping produced transition discrepancies of order 10^-12, which were amplified by the nonparametric detrending into compensation-metric discrepancies up to approximately 10^-8 for some anchors.

Therefore the direct-local implementation does not satisfy the preregistered frozen-anchor recovery requirement and will not be used for the landscape run.

### 3. Replacement implementation

The landscape implementation will use an exact compressed representation of the global-order amplitude-damping calculation.

For subsystem A and complement B, retain only matrix elements contributing to the final partial trace:

    B[a,r,a_prime] = rho[(a,r),(a_prime,r)]

where a and a_prime index the four-qubit subsystem and r is the common complement basis index.

The two TFIM branches are first transformed into the (A,B) basis ordering and combined into the same time-reversal mixture as Paper 30.

Amplitude damping is then applied in the exact original global qubit order q = 0,1,...,N-1.

For q in A, the full ket/bra four-block update is applied.

For q outside A, only the diagonal complement blocks required by the eventual partial trace are propagated:

    B[...,r_q=0,...] <- B[...,r_q=0,...] + gamma B[...,r_q=1,...]
    B[...,r_q=1,...] <- (1-gamma) B[...,r_q=1,...]

The reduced state rho_A is then obtained by summing over the complement index r.

This representation is mathematically equivalent to the full global density-matrix channel followed by partial trace, but preserves the floating-point update ordering relevant to frozen-anchor recovery.

### 4. Validation performed before freezing this amendment

For the known recovery anchor C1 / N=10 / P=6 / A=(3,4,5,6), the compressed global-order implementation reproduced the original full-global calculation exactly for representative gamma values:

- gamma index 0: max |Delta rho_A| = 0
- gamma index 10: max |Delta rho_A| = 0
- gamma index 30: max |Delta rho_A| = 0

Across the complete trajectory it also reproduced exactly:

- archived local-state fields,
- d_W and d_mod transitions,
- transition status,
- all five frozen Pearson residual correlations.

The implementation is not accepted for the landscape until all frozen C0-C3 anchors for N10 and N12 pass the original 1e-10 recovery tolerance.

### 5. Scientific invariants unchanged

This amendment changes no scientific parameter, observable, threshold, bandwidth, system size, circuit depth, subsystem position, support criterion, output definition, or interpretation boundary.

In particular, the following remain exactly as frozen:

- N in {8,10,12,14,16};
- P in {1,...,12};
- every contiguous four-qubit subsystem;
- the s and gamma grids;
- W and modular-sector definitions;
- transition_step;
- nonparametric detrending;
- Pearson/Spearman calculations;
- alpha optimization;
- R_C;
- the historical Paper 30 support criterion;
- reflection-symmetry tolerance 1e-10;
- frozen-anchor recovery tolerance 1e-10.

### 6. Interpretation

The amendment is an implementation-level numerical-equivalence correction required by the preregistered recovery gate.

It does not constitute a new physical hypothesis and does not use any previously unknown landscape outcome.

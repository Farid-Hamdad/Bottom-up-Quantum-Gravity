# Paper 30 — Causal Compensation Landscape

## POST-UNBLINDING AMENDMENT — RANK-DEFICIENT MODULAR SECTOR

### 1. Status and timing

This amendment is frozen after the first official landscape execution attempt and before any continuation of the landscape scan.

The frozen implementation commit was:

54c8495 Freeze Paper 30 causal compensation landscape implementation

The frozen script SHA256 was:

7701d483e1a0fc5f5783e33a74955e82da7af374ea588a49c8fbdea24509090f

The first official run passed all eight C0-C3 recovery anchors exactly and then stopped at the first landscape cell encountered in the preregistered iteration order.

No later unknown landscape cell was inspected before this amendment.

### 2. Observed first-run condition

The first stopping cell was:

- N = 8
- P = 1
- A_start = 0
- A = (0, 1, 2, 3)
- state_index = 0
- s = 0
- gamma = 0

The reduced state produced by the compressed implementation and by the original full-global implementation was exactly identical:

max_abs_rho_delta = 0

Its numerical spectrum had twelve eigenvalues at approximately machine-zero scale and four positive eigenvalues above 1e-2.

The numerical rank was 4 of 16 for thresholds from 1e-15 through 1e-10.

The minimum eigenvalue was approximately:

-3.32503514717640009e-17

This is far inside the inherited density positivity tolerance 1e-12 and therefore does not indicate a physically meaningful negative density eigenvalue.

The two time-reversal branches each had Schmidt rank 2 at tolerance 1e-12 across the A versus complement bipartition.

The rank deficiency is therefore treated as a structural property of this finite-depth reduced state, not as a compressed-trajectory implementation error.

### 3. Consequence for the modular observable

The inherited Paper 30 modular definition is retained unchanged:

K_A = -log(rho_A)

with the inherited strict gate that rho_A must be positive definite.

No eigenvalue clipping, eigenvalue floor, pseudologarithm, support-restricted logarithm, regulator, or modified positivity threshold will be introduced after unblinding.

If rho_A is not positive definite at any gamma-grid state required for a cell, the normalized modular sector is undefined for that cell under the frozen observable definition.

### 4. Frozen handling of undefined cells

A cell encountering the inherited POSITIVITY_GATE_FAIL is retained in the 540-cell landscape and assigned:

analysis_status = MODULAR_SECTOR_UNDEFINED_RANK_DEFICIENT

The cell is not removed from the grid.

For that cell:

- geometric descriptors remain recorded;
- the first failing state index, s, gamma, lambda_min, lambda_max, and numerical rank are recorded;
- mutual-information quantities already well-defined may be retained as diagnostics;
- modular A_ij, v_ij, modular transitions, residual compensation correlations, alpha, R_C, and bandwidth support values are not fabricated;
- PAPER30_SUPPORTED is not interpreted as a valid compensation failure.

The same principle applies to MODULAR_STD_ZERO if it is encountered later: the cell is retained and marked as modular-sector undefined with its exact failure reason.

No other RuntimeError is silently converted into an invalid cell. Unexpected errors remain fatal and require separate investigation.

### 5. Aggregate reporting

All preregistered cells remain in the denominator of the landscape inventory.

Aggregate outputs will explicitly report, for every relevant grouping:

- total cell count;
- valid compensation-analysis cell count;
- modular-sector undefined cell count;
- supported valid cell count;
- fraction valid;
- fraction supported among all preregistered cells;
- fraction supported among valid cells when at least one valid cell exists.

The original question of where Paper 30 compensation is supported will therefore be reported together with the domain on which the frozen modular observable is actually defined.

Undefined cells are not relabeled as ordinary compensation failures.

### 6. Reflection control

For a mirrored pair, reflection comparison of compensation quantities is performed only when both cells have valid modular analyses.

If both mirrored cells are undefined for the same inherited modular reason at corresponding reflected states, this is recorded as an exact qualitative reflection-status match.

If only one member is undefined, or failure reasons or corresponding failure locations disagree beyond the frozen reflection mapping, the reflection control fails and requires investigation.

The frozen numerical reflection tolerance remains 1e-10 for quantities that are defined.

### 7. Scientific interpretation

This amendment does not change the landscape grid, TFIM preparation, amplitude-damping channel, W observable, modular-Hamiltonian definition, detrending procedure, bandwidths, historical support criterion, reflection tolerance, anchor tolerance, or any numerical threshold.

It only specifies how to preserve and report preregistered cells for which the already-frozen modular observable is mathematically undefined.

The existence and distribution of rank-deficient cells is itself a result of the preregistered finite-size and finite-depth landscape and must be reported rather than regularized away.

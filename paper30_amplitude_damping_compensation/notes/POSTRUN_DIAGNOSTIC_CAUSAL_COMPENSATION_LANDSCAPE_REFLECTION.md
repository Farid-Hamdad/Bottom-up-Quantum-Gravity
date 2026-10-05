# Paper 30 — Causal Compensation Landscape

## Post-run diagnostic of the frozen reflection control

### 1. Status

The Attempt 002 landscape results are frozen at commit:

23acc83 Freeze Paper 30 landscape attempt 002 results

The preregistered reflection tolerance remains unchanged at:

1e-10

The frozen Attempt 002 result remains:

reflection_control_passed = False

No threshold, support criterion, observable definition, landscape cell, or frozen result is changed by this diagnostic.

### 2. Attempt 002 landscape result

The complete preregistered landscape contains 540 cells.

Observed:

- valid compensation-analysis cells: 485 / 540
- modular-sector undefined cells: 55 / 540
- Paper 30 supported cells among valid cells: 68 / 485
- reflection checks: 300
- reflection checks passed: 242
- reflection checks failed: 58
- maximum frozen reflection discrepancy: 1.3766901396650155e-08

All observed reflection failures retained exact Paper 30 support-status agreement between mirrored partners.

### 3. Rank-deficient modular domain

The 55 modular-sector undefined cells form a structured depth pattern.

For every N in {8,10,12,14,16}:

- P = 1: every four-qubit position is modular-sector undefined;
- P = 2: only the two edge positions are modular-sector undefined;
- P >= 3: all preregistered four-qubit positions are modular-sector defined.

The failures occur through the inherited POSITIVITY_GATE_FAIL at gamma = 0 and correspond to rank-deficient reduced density matrices, not to a compressed-trajectory mismatch.

This is reported as a domain-of-definition property of the frozen modular observable K_A = -log(rho_A). No eigenvalue floor, pseudologarithm, or post-unblinding regularization is introduced.

### 4. Reflection-failure decomposition

The 58 failed reflection checks separate into two numerically distinct classes.

#### Class I — spectral-logarithm amplification

Five failed pairs have extreme gamma=0 conditioning:

condition number approximately 1.442e9.

They are the edge-reflection pairs at P = 3 for:

N = 8, 10, 12, 14, 16.

For the worst case N=16, P=3, A=0 <-> 12:

- max reflected rho_A discrepancy: approximately 5.27e-16 at gamma=0;
- lambda_min: approximately 4.036e-10;
- condition number: approximately 1.442e9;
- K_raw discrepancy: approximately 5.93e-08;
- K_tilde discrepancy: approximately 6.58e-09;
- A_ij discrepancy: approximately 1.38e-08;
- M_mod discrepancy: approximately 7.43e-09.

Thus machine-scale reflected-state differences are strongly amplified by the spectral map -log(rho_A) when rho_A is extremely ill-conditioned.

#### Class II — normalization / detrending / standardization amplification

The remaining 53 failed pairs have gamma=0 condition numbers below 1e6.

Across these 53 pairs:

- max reflected raw-W discrepancy:
  min 6.661e-16,
  median 1.776e-15,
  max 3.775e-15;

- max reflected normalized-W discrepancy:
  min 5.131e-13,
  median 1.480e-12,
  max 3.409e-12;

- max reflected d_W discrepancy:
  min 6.094e-13,
  median 1.128e-12,
  max 4.067e-12;

- max reflected d_mod discrepancy:
  min 4.366e-15,
  median 1.370e-14,
  max 2.361e-14;

- max reflected z-residual W discrepancy:
  min 1.855e-09,
  median 6.575e-09,
  max 2.744e-08;

- max reflected z-residual modular discrepancy:
  min 1.962e-12,
  median 6.754e-12,
  max 5.252e-11.

All 53 non-extreme failures satisfy:

max reflected raw-W discrepancy < 1e-12

and:

max reflected z-residual W discrepancy > 1e-10.

None satisfies:

max reflected z-residual modular discrepancy > 1e-10.

The median amplification ratio

max_delta_zW / max_delta_residual_W

is approximately 4.197e3, while the median inverse minimum residual standard deviation is approximately 4.327e3.

A representative case, N=14, P=7, A=4 <-> 6, gives:

- raw W discrepancy: 1.332e-15;
- normalized W discrepancy: 3.409e-12;
- d_W discrepancy: 3.242e-12;
- residual d_W discrepancy: 3.121e-12;
- z-residual W discrepancy: 1.984e-08;
- residual standard deviation: approximately 1.457e-04;
- observed z amplification: approximately 6.359e3;
- inverse residual standard deviation: approximately 6.864e3;
- Pearson discrepancy: approximately 1.111e-09;
- R_C discrepancy: approximately 2.052e-09;
- Paper 30 bandwidth support status: exact mirror match.

This establishes a second numerical amplification mechanism distinct from the ill-conditioned modular logarithm.

### 5. Interpretation of the frozen reflection control

The preregistered reflection control remains formally failed because 58 of 300 comparisons exceed the frozen tolerance 1e-10.

The post-run diagnostic does not convert this result into a PASS.

However, the numerical evidence does not support a physical left-right asymmetry of the underlying reduced states.

The failed checks are attributable to:

1. spectral amplification through -log(rho_A) in five extremely ill-conditioned P=3 edge cases; and
2. amplification of machine-scale W differences through normalization, transition construction, local detrending, and residual standardization in the other 53 cases.

The exact Paper 30 support classification is reflection-consistent for all failed pairs.

### 6. Scientific consequence

The causal-compensation landscape should therefore be reported with two simultaneous statements:

- the strict preregistered reflection control failed numerically and interpretation under that gate remains formally blocked;
- the post-run audit finds no evidence that this failure represents a genuine physical breaking of reflection symmetry.

No tolerance is relaxed and no result is reclassified.

The substantive landscape observations remain descriptive:

- the frozen modular observable has a depth-dependent domain of definition;
- valid compensation is present only in a minority of valid cells;
- the original Paper 30 compensation regime is therefore conditional on finite-size, depth, subsystem placement, and analysis scale rather than universal.

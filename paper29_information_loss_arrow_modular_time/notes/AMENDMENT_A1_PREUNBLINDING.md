# Paper 29 — Amendment A1 (Pre-unblinding technical correction)

## Status

Technical amendment after the first execution attempt stopped at the Paper 28 recovery gate.

No Paper 29 primary, secondary, control, or replication result was produced or inspected.

## Failure observed

The first official Paper 29 execution stopped with:

`RuntimeError: PAPER28_S0_RECOVERY_FAIL_N8`

The failure occurred at the mandatory s=0 recovery check before any Paper 29 result files were generated.

## Diagnosis

The Paper 29 implementation used:

`PAPER28_RECOVERY_ATOL = 1e-12`

However, the frozen Paper 28 public reproduction verifier uses a numerical tolerance of:

`1e-10`

for non-integer reproduced metrics.

The observed N8/A4 recovery differences were numerical roundoff-level differences and were within the official Paper 28 tolerance.

The largest relevant discrepancy was approximately:

`1.14e-11`

for the condition number.

## Correction

The Paper 29 recovery tolerance was changed from:

`1e-12`

to:

`1e-10`

No other line of the implementation was modified.

## Scientific status

This amendment does not change:

- the irreversible channel;
- the s grid;
- the reference state;
- the subsystem;
- the pair ordering;
- the geometry definition;
- the modular observable;
- the primary endpoint;
- the secondary metrics;
- the zero-step rule;
- the controls;
- the replication protocol;
- the interpretation policy.

This is strictly a technical alignment of the Paper 29 recovery gate with the published Paper 28 verification tolerance.

## Script hashes

Original frozen script SHA256:

`f05efa730ceded76277e05414ceb1d5418b568da7e83e491d97d5f35e721d1ba`

Amended script SHA256:

`28030768ccb38b67af83538a37a07e933f0a95658b4ef91a5d1a24df9f45de48`

## Unblinding status

No Paper 29 scientific result had been observed when this amendment was made.

The next authorized action is to commit this amendment and the amended script before re-execution.

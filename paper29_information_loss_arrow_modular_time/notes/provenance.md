# Paper 29 — Provenance and Freeze Record

## Repository context

Public theoretical repository:

`Farid-Hamdad/Bottom-up-Quantum-Gravity`

Branch:

`bup_cosmology`

Paper directory:

`paper29_information_loss_arrow_modular_time/`

## Parent scientific benchmark

Paper 29 directly extends:

`paper28_modular_time_geometry/`

Paper 28 public commit:

`47cfa760899713fb9b5ba8c63e059f1315159124`

Commit message:

`Add Paper 28 modular time geometry benchmark`

Paper 28 provides the frozen TFIM state preparation, modular normalization, pair ordering conventions, and exact short-time double-commutator construction used by Paper 29.

## Paper 29 frozen implementation

Script:

`scripts/paper29_arrow_of_modular_time_v1.py`

Frozen SHA256:

`f05efa730ceded76277e05414ceb1d5418b568da7e83e491d97d5f35e721d1ba`

The SHA256 was verified after manual transfer to the local Mac repository and matched the audited implementation byte-for-byte.

## Pre-unblinding state

At freeze time:

- the Paper 29 preregistration was finalized;
- the implementation was created;
- the script SHA256 was frozen;
- no Paper 29 simulation result had been inspected;
- no primary endpoint had been computed;
- no Paper 29 result file had been generated.

## Execution policy

Before first execution:

1. stage only the Paper 29 preregistration, provenance file, and frozen script;
2. commit them in Git;
3. verify the committed script SHA256;
4. only then execute the benchmark;
5. retain all positive, null, negative, failed, and undefined outcomes.

No primary metric, grid, normalization, pair selection, control, or zero-step rule may be altered after first execution.

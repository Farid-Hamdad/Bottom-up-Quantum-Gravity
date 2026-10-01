# Paper 28 — Provenance

Paper 28 is assembled from frozen deterministic TFIM modular-flow campaigns performed before manuscript construction.

## Frozen N10 audit chain

- preregistration commit: `c5b537a5b263e9eef50aa260a847650bdc330479`
- primary script freeze: `aed4186d2011db0caae2ad5b610cd7eacc7d64ea`
- primary result freeze: `db370627e37658eb7be709f82789803d770bc3ae`
- secondary curve generator freeze: `42c9594d44e8ac98effd7f0e0b366d251cc0c186`
- secondary raw-data freeze: `aaada1293ace117ba4390d5ead43ef9d2cd52206`
- descriptive-analysis script freeze: `786fd35d6c089c89b7cbd8c490b0d132bb24ae49`
- descriptive-result freeze: `413a9ce70a20df5d8121857029b33effbd44dea0`

The public reproduction script in this paper is a clean standalone implementation of the frozen mathematical protocol. It does not require the private research repository.

## Frozen public reference values

The file `results/paper28_modular_time_v1/frozen_reference_metrics.csv` records the reference values transcribed from the frozen N6/N8/N10 audit chain. The verification script recomputes the benchmark independently and checks those values to numerical tolerance.

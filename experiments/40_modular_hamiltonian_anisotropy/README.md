# Modular Hamiltonian Anisotropy — BuP v3

## Objective

Test whether modular Hamiltonians encode a state-dependent anisotropy
and extract an effective scaling law:

    η_eff ≈ f(η_micro)

## Setup

- System size: N = 12
- Subsystem: |A| = N/2
- States tested:
  - Ground state
  - Thermal state
  - Random state

## Key Result

Ground states exhibit strong anisotropy amplification:

    η_eff > η_micro

while random/thermal states remain near isotropic.

## Outputs

- `aggregate_rows.csv` → raw results
- `aggregate_summary.json` → aggregated stats
- `figures/` → scaling plots
- `tables/` → summary tables

## Reproduce

```bash
bash run_experiment.sh


---
 SCRIPT DE RUN

`run_experiment.sh`

```bash
#!/bin/bash

OUTDIR="results"

mkdir -p $OUTDIR

python analyze_bup_v3_results.py \
  --output-dir $OUTDIR

python make_bup_v3_figures_tables_fixed.py \
  --csv $OUTDIR/aggregate_rows.csv \
  --output-dir $OUTDIR

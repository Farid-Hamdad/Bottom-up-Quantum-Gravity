# Experiment 40 — Micro dimension coherence across sizes and couplings

This experiment collects the numerical evidence that the microscopic emergent
dimension in BuP remains close to

\[
d_s \approx 2.61 \pm 0.03
\]

across:

- system sizes $N = 6, 8, 9, 12, 16$,
- non-local couplings $\lambda \in [0,1]$.

## Purpose

The goal is to test whether the emergent dimension is:

1. stable in the non-local coupling $\lambda$,
2. stable in system size $N$,
3. numerically compatible with the galactic SPARC scale
   $d_{\mathrm{eff}} = 2.644 \pm 0.295$,
4. compatible with the late-time cosmological range inferred from
   DESI 2024 + Pantheon+.

## Files

### Scripts

- `scripts/bup_ds_fss_scan_v1.py`
- `scripts/bup_plot_micro_coherence_v1.py`

### Results

- `results/raw_scan.csv`
- `results/fss_extrapolated.csv`
- `results/d_lambda_summary.csv`
- `results/fig_micro_coherence.svg`

## Main result

The raw data show a stable microscopic plateau:

\[
d_s(\lambda, N) \approx 2.61 \pm 0.03
\]

without a significant drift over the explored values of $N$ and $\lambda$.

## Important note

The finite-size-scaling extrapolation is exploratory and not used as the main conclusion.

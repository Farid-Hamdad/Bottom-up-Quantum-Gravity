# Paper 28 — Reproducibility

## Environment

Required:

- Python 3.10+
- NumPy
- Matplotlib (figures only)
- a LaTeX distribution with `pdflatex` for the manuscript PDF

## Full numerical reproduction

```bash
cd paper28_modular_time_geometry
./scripts/run_all.sh
```

This performs three steps:

1. exact finite-system TFIM reconstruction and modular-flow calculation;
2. comparison against frozen reference metrics;
3. regeneration of the publication figures.

The calculation uses no randomness.

## Frozen protocol

- open TFIM chain
- J = h = 1
- initial state |+>^N
- Strang splitting with p = 6 and dt = 0.35
- time-reversal mixture of +t and -t evolved states
- K_A = -log rho_A
- K_A normalized by the unweighted mean and population standard deviation of its eigenvalues
- primary profile v_ij = sqrt(A_ij)
- secondary profile C_ij(tau)
- tau grid: 1024 equally spaced points on [0,120]
- all 1023 tau > 0 points retained in the descriptive summary
- positive-scale least-squares residual comparison against W and bare chain adjacency

No adaptive time window, smoothing, pair selection, regularization or post-hoc threshold is used.

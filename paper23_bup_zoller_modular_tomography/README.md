# Paper 23 — BuP/Zoller Modular Tomography

This folder contains the scientific paper connecting the BuP tomographic principle to the Zoller/Joshi 2023 trapped-ion entanglement data.

Core result:

```math
\rho_{\rm ent}^{\rm cut}(j)
=
\sum_{k\notin A} I(j:k)
\propto
T_{\rm mod}(j)
=
\frac{1}{\beta_j}.
```

## Contents

```text
main.tex
main.pdf
scripts/
results/
  modular_profiles/
  raw_mi/
  fig1_fig3/
figures/
```

## Scripts

- `analyze_zoller_bup.py`: extracts modular profiles and tests `1/beta_j`.
- `raw_pairwise_mi_test.py`: reconstructs pairwise mutual information from raw Pauli data.
- `shadow_mi_test.py`: compares global, internal, and cut mutual-information densities.
- `analyze_fig1_fig3_bup.py`: analyzes entropy scaling, MI decay, and disjoint-subsystem diagnostics.

## Main Numerical Claims

- Ground states: parabolic modular profiles compatible with CFT form `j(L-j)`.
- Excited states: volume-law entropy with linear fits `R^2 ≈ 0.995`.
- Naive global density `sum_k I(j:k)` does not match `1/beta_j`.
- Cut density `sum_{k notin A} I(j:k)` correlates better with `1/beta_j`.
- Disjoint-subsystem links improve reconstruction fidelity.

## Data Note

The raw `.mat` files from Joshi et al. are not included here by default. Keep derived CSV/JSON summaries in `results/`, and document the original data source in the paper.

# Modular Hamiltonian Anisotropy — BuP v3

## Objective

This experiment tests whether modular Hamiltonians encode a **state-dependent geometric anisotropy** in the Bottom-Up (BuP) framework.

The goal is to extract an effective scaling law:

    η_eff ≈ f(η_micro)

and determine how this behavior depends on the underlying quantum state.

---

## Physical Motivation

In the BuP framework:

- Geometry is not postulated but **emerges from entanglement**
- The modular Hamiltonian \( K_A \) acts as a **local probe of this structure**
- Anisotropy in \( K_A \) reflects **effective interaction geometry**

This experiment tests whether:

- Ground states exhibit **anisotropy amplification**
- Thermal and random states remain **near isotropic**

---

## Setup

- System sizes: `N = 8, 10, 12`
- Subsystem sizes: `|A| = 3, 4, 5, 6`
- State types:
  - `ground`
  - `thermal`
  - `random`
- Parameters:
  - \( J_{xy} = 1.0 \)
  - \( J_z \in [1.0, 2.0] \)
  - inverse temperature: `β = 2.0`
- Boundary conditions: periodic
- Local perturbation:
  - site: `2`
  - strength: `0.2`

---

## Scripts Used

All scripts are located in the `jauge/` directory:

- `bup_modular_referee_pipeline_v3.py`  
  → runs the full simulation

- `analyze_bup_v3_results.py`  
  → aggregates and fits scaling laws

- `make_bup_v3_figures_tables_fixed.py`  
  → generates figures and tables

---

## Reproducibility

Run the full pipeline from the project root:

```bash
python jauge/bup_modular_referee_pipeline_v3.py \
  --n-list 8,10,12 \
  --subsystem-sizes 3,4,5,6 \
  --jxy-list 1.0 \
  --jz-list 1.0,1.1,1.25,1.5,2.0 \
  --state-kinds ground,thermal,random \
  --beta 2.0 \
  --periodic \
  --perturb-site 2 \
  --perturb-strength 0.2 \
  --output-dir experiments/40_modular_hamiltonian_anisotropy/results

python jauge/analyze_bup_v3_results.py \
  --csv experiments/40_modular_hamiltonian_anisotropy/results/aggregate_rows.csv \
  --output-dir experiments/40_modular_hamiltonian_anisotropy/results

python jauge/make_bup_v3_figures_tables_fixed.py \
  --csv experiments/40_modular_hamiltonian_anisotropy/results/aggregate_rows.csv \
  --output-dir experiments/40_modular_hamiltonian_anisotropy/results

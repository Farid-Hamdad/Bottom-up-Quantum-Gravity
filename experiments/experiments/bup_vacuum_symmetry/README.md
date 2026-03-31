# BuP Vacuum — Effective Symmetry Selection

This experiment tests whether the ground state of a BuP-inspired Hamiltonian
exhibits an emergent symmetry structure.

We consider an XXZ-type Hamiltonian with weak random local fields:

H = -Jzz Σ ZiZj - Jxy Σ (XiXj + YiYj) + Σ hi · σi

The goal is to determine whether the vacuum remains isotropic (SO(3)-like)
or selects a preferred axis (U(1)-like).

---

## Method

1. Compute the ground state |ψ₀⟩ of H
2. Extract 2-qubit reduced density matrices
3. Build correlation tensors T_ij = Tr(ρ_ij σ_a ⊗ σ_b)
4. Project T_ij → SO(3) rotations via polar decomposition
5. Map rotations → Lie algebra vectors
6. Analyze alignment of these vectors:

- strong alignment → U(1)-like
- isotropic distribution → SO(3)-like

---

## Observables

- Order parameter (OP): mean alignment with dominant axis
- Singular ratio: λ₁ / λ₂ (SVD of Lie vectors)
- Angular dispersion σ
- Fraction of aligned edges

Classification:

- U(1)-like if ratio > 2.5 and σ < 0.4
- otherwise SO(3)-like

---

## Results

### Robustness (h = 0.1)

- U(1)-like in 95% of seeds
- OP ≈ 0.97 (vs 0.62 random)
- strong axis alignment

### Disorder scan (h_scale)

- stable U(1)-like regime for h ≲ 0.2
- degradation starts around h ≈ 0.3
- mostly lost for h ≥ 0.5

### Anisotropy scan (Jzz / Jxy)

- isotropic regime for Jzz/Jxy < 1
- transition near Jzz/Jxy ≈ 1
- robust U(1)-like for Jzz/Jxy ≥ 2

### Random control

- no U(1)-like detection
- OP ≈ 0.62
- confirms non-trivial structure

---

## Interpretation

The BuP effective vacuum is not isotropic.

It dynamically selects a preferred axis in a broad low-disorder regime,
leading to an emergent U(1)-like symmetry.

This effect:

- is absent in random states
- is controlled by XXZ anisotropy
- is degraded by strong disorder

---

## Implications for BuP

- The vacuum carries intrinsic structure (not neutral)
- Entanglement geometry is anisotropic
- A residual symmetry can emerge dynamically
- The vacuum behaves as a structured medium

This supports the idea that geometry and symmetry in BuP
are not imposed but emerge from the entanglement structure.

---

## Files

- `bup_vacuum_v3_standalone.py`
- `bup_vacuum_v3_results.csv`
- `bup_vacuum_v3_summary.json`

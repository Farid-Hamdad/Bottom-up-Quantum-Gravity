# BuP Cosmology — Emergent Spatial Dimension

**Paper:** Cosmological constraints on an emergent spatial dimension: variable dimension, dark energy and σ₈ in the Bottom-Up Quantum Gravity framework

**Author:** Farid Hamdad  
**HAL:** [hal-05590614v1](https://hal.science/hal-05590614v1) (April 2026)  
**GitHub:** [Farid-Hamdad/Bottom-up-Quantum-Gravity](https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity)

---

## Results

| Observable | BuP | Data | Tension |
|---|---:|---|---|
| r_drag (BAO) | unchanged | DESI 2024 | < 1σ ✓ |
| σ₈ | 0.772 | KiDS: 0.766 ± 0.020 | +0.3σ ✓ |
| d(z=0) | 2.40 | SPARC: 2.644 ± 0.295 | +0.8σ ✓ |
| ΔBIC | −3.71 | DESI + Pantheon+ | signal ✓ |
| CMB (θs) | unchanged | Planck 2018 | < 0.1σ ✓ |
| H₀ | 63.7 | Planck: 67.4 ± 0.5 | 3.7σ ⚠ |
| t₀ (age) | 14.3 Gyr | ΛCDM: 13.8 Gyr | +500 Myr |

---

## BuP Constants

```text
DC = 3.059842935509462   # critical dimension
α  = 1.78                # fixed by microscopic simulations (N = 6–16)
```

---

## Structure

```
bup_cosmology/
├── paper2_background-dimension/
│   ├── README.md  
│   ├── main.tex                    ← Full LaTeX article (EN)
│   ├── fig1_camb_scan.pdf          ← Figure 1: CAMB scan (σ₈, r_drag, θs)
│   ├── fig2_models_coherence.pdf   ← Figure 2: d(z), w(z), multi-scale
│   └── fig3_sigma8_geff.pdf        ← Figure 3: G_eff(z) and σ₈ scan
│
└── scripts/
    ├── bup_paper2_final.py         ← BAO+SN fit, Paper 2 model X(z)
    ├── bup_sigma8_final.py         ← σ₈ via coherent G_eff ODE
    ├── bup_camb_geff_patch.py      ← CAMB G_eff patch + instructions
    ├── bup_geff_module.f90         ← Fortran module for CAMB
    └── bup_propagator.py           ← Gravitational propagator on BuP graph
```

---

## Quickstart

```bash
conda create -n bottomup python=3.11
conda activate bottomup
pip install numpy scipy matplotlib camb

# Paper 2 — BAO + SN fit
python scripts/bup_paper2_final.py

# Paper 3 — σ₈ with coherent G_eff
python scripts/bup_sigma8_final.py

# CAMB G_eff patch (requires CAMB source)
python scripts/bup_camb_geff_patch.py --full
```

---

## Physical Mechanism

The key result is geometric. The critical dimension **DC = 3.0598 > 3** plays a triple role:

**1. CMB preserved**  
`d(z → ∞) → DC` ensures standard recombination physics and preserves the CMB angular scale.

**2. Modified gravitational coupling**  
`G_eff(z) = 2 / (d(z) − 1) · G`

- At high redshift: `d ≈ DC > 3` → `G_eff < G` → slight suppression of early structure growth  
- At low redshift: `d < 3` → `G_eff > G` → enhanced late-time clustering

**3. σ₈ resolved**  
Growth integrated from z ≈ 1000 to z = 0 gives:

```
σ₈ = 0.772
```

consistent with KiDS-1000 (0.766 ± 0.020) and DES Y3 (0.776 ± 0.017), without additional free parameters.

The slight over-dimensionality of the primordial Universe (**DC > 3**) provides a unified geometric explanation for late-time acceleration, suppressed structure growth, and preserved CMB physics — linking emergent quantum geometry to BAO, CMB, weak lensing, galaxy rotation curves, and early galaxy formation.

---

## Model (Paper 2)

```
X(z) = X₀ · (1 + z)^β              effective intermediate variable
d(z) = DC − Δd / (1 + X(z)^α)      emergent dimension profile
α    = 1.78                          fixed by BuP simulations (not a free parameter)
```

**Best-fit (DESI 2024 + Pantheon+):**

```
X₀ = 0.537    β  = 2.000    Δd = 0.878
H₀ = 63.67    Ωm = 0.3448
ΔAIC = −6.23  ΔBIC = −3.71
```

---

## Gravitational Propagator (Paper 3 — preliminary)

```bash
python scripts/bup_propagator.py --N 12 16 --lambda-n 6 --n-repeat 2
```

Tests whether `G(r) ~ r^{−(d_s−2)}` on the entanglement graph.

**Current result:** γ ≈ 0 for N ≤ 16 — insufficient system size.  
Reliable measurement likely requires N ≥ 100.

---

## JWST Prediction

BuP predicts a slightly older Universe:

```
t₀ ≈ 14.3 Gyr   (vs. ΛCDM: 13.8 Gyr)
```

This follows directly from H₀ ≈ 63.7 km/s/Mpc and Ωm = 0.3448, providing approximately **~500 Myr** of additional cosmic time for structure formation at z > 10.

Combined with the modified gravitational coupling `G_eff(z) = 2/(d(z)−1)·G`, BuP predicts:

- Slight suppression of very early growth (high z, d > 3, G_eff < G)
- Enhanced late-time clustering (low z, d < 3, G_eff > G)

The net effect allows both sufficient time for early massive galaxy formation and enhanced clustering without excessive σ₈.

This may help alleviate tensions raised by recent JWST observations (JADES, CEERS) concerning unexpectedly massive galaxies at very high redshift.

> **This should be understood as a qualitative and testable prediction**, not as a claimed explanation of all JWST observations. It arises from two independent mechanisms: (1) older Universe due to lower H₀, and (2) modified growth history from variable d(z).

---

## Citation

```bibtex
@misc{hamdad2026bup,
  author  = {Farid Hamdad},
  title   = {Cosmological constraints on an emergent spatial dimension:
             variable dimension, dark energy and $\sigma_8$
             in the Bottom-Up Quantum Gravity framework},
  year    = {2026},
  note    = {HAL preprint hal-05590614v1},
  url     = {https://hal.science/hal-05590614v1}
}
```

---

## Related

- **Main branch** — BuP microscopic simulations (N = 6–16 qubits, FSS, entanglement)
- **Paper** — Full LaTeX source in `paper/main.tex`
- **HAL** — https://hal.science/hal-05590614v1

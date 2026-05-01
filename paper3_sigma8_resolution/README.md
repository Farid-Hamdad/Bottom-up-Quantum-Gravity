# BuP Paper 3 — Geometric Resolution of the σ₈ Tension

**Title FR:** Résolution géométrique de la tension σ₈ par une dimension cosmologique émergente

**Title EN:** Geometric Resolution of the σ₈ Tension from an Emergent Cosmological Dimension

**Author:** Farid Hamdad
**ORCID:** [your ORCID]
**Year:** 2026

---

## Key Result

σ₈ = 0.772, consistent with KiDS-1000 at +0.3σ, without additional free parameters.

| Model | σ₈ | KiDS tension | Verdict |
|---|---|---|---|
| ΛCDM (Planck 2018) | 0.811 | +2.3σ | Tension |
| BuP — background only | 0.861 | +4.7σ | Artefact |
| BuP — coherent G_eff (this work) | **0.772** | +0.3σ | ✓ |

## Physical Mechanism

G_eff(z) = 2/(d(z)−1)·G applied coherently to both background and perturbations.

- High z: d → DC = 3.0598 > 3 → G_eff ≈ 0.97G → slight growth suppression (~99% of history)
- Low z: d < 3 → G_eff > G → enhanced clustering (~1% of history)
- Net: σ₈ = 0.772 ✓

## Structure

```
paper/
├── main_fr.tex          ← Article complet (FR)
├── main_en.tex          ← Article complet (EN)
└── fig3_sigma8_geff.pdf ← Figure : G_eff(z) et scan σ₈

scripts/
├── bup_sigma8_final.py      ← Calcul σ₈ via ODE croissance linéaire
├── bup_paper2_final.py      ← Ajustement BAO+SN, modèle X(z)
├── bup_camb_geff_patch.py   ← Patch CAMB avec G_eff(z)
└── bup_geff_module.f90      ← Module Fortran pour CAMB
```

## Quickstart

```bash
pip install numpy scipy matplotlib
python scripts/bup_sigma8_final.py
```

## BuP Constants

```
DC    = 3.059842935509462
alpha = 1.78  (fixed by microscopic simulations N=6–16)
X0    = 0.537,  beta = 2.000,  Delta_d = 0.878
H0    = 63.67 km/s/Mpc,  Omega_m = 0.3448
```

## Related

- **Paper 1/2** — HAL: [hal-05590614v1](https://hal.science/hal-05590614v1)
- **Paper 4** — Local dimensional reduction, JWST prediction
- **GitHub** — branch `cosmology`

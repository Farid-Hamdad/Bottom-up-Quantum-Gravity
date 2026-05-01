# BuP Paper 3 — Résolution géométrique de la tension σ₈

**Titre :** Résolution géométrique de la tension σ₈ par une dimension cosmologique émergente

**Auteur :** Farid Hamdad
**ORCID :** [votre ORCID]
**Année :** 2026

---

## Résultat central

σ₈ = 0.772, cohérent avec KiDS-1000 à +0.3σ, sans paramètre libre supplémentaire.

| Modèle | σ₈ | Tension KiDS | Verdict |
|---|---|---|---|
| ΛCDM (Planck 2018) | 0.811 | +2.3σ | Tension |
| BuP — fond seul | 0.861 | +4.7σ | Artefact |
| BuP — G_eff cohérent (ce travail) | **0.772** | +0.3σ | ✓ |

---

## Mécanisme physique

G_eff(z) = 2/(d(z)−1)·G appliqué de façon cohérente à la fois à l'expansion du fond et à la croissance linéaire des perturbations.

- Haut z : d → DC = 3.0598 > 3 → G_eff ≈ 0.97G → légère suppression de croissance (~99% de l'histoire cosmique)
- Bas z : d < 3 → G_eff > G → renforcement du regroupement tardif (~1% de l'histoire)
- Résultat net : σ₈ = 0.772 ✓

La valeur précédente σ₈ ≈ 0.86 était un artefact du traitement incohérent (fond seul modifié). Une fois G_eff(z) inclus dans les perturbations, la tension disparaît.

---

## Prédiction JWST

BuP prédit un Univers légèrement plus vieux :

```
t₀ = 14.24 Gyr  (vs 13.80 Gyr en ΛCDM)
Δt ≈ +430 Myr disponibles à z > 10
```

Une description quantitative de l'excès JWST nécessite l'introduction d'inhomogénéités locales d(z,δ), développées dans Paper 4.

---

## Structure

```
paper/
├── main_fr.tex              ← Article complet (version française)

figures
└── fig3_sigma8_geff.pdf     ← Figure : G_eff(z) et scan σ₈

scripts/
├── bup_sigma8_final.py      ← Calcul σ₈ via ODE de croissance linéaire
├── bup_paper2_final.py      ← Ajustement BAO+SN, modèle X(z)
├── bup_camb_geff_patch.py   ← Patch CAMB avec G_eff(z)
└── bup_geff_module.f90      ← Module Fortran pour CAMB
```

---

## Démarrage rapide

```bash
pip install numpy scipy matplotlib

# Calcul de σ₈ avec G_eff cohérent
python scripts/bup_sigma8_final.py

# Ajustement BAO+SN
python scripts/bup_paper2_final.py
```

---

## Constantes BuP

```
DC    = 3.059842935509462   # dimension critique
alpha = 1.78                # fixé par simulations microscopiques N=6–16
X0    = 0.537,  beta = 2.000,  Delta_d = 0.878
H0    = 63.67 km/s/Mpc,  Omega_m = 0.3448
sigma8 = 0.772
```

---

## Cohérence multi-échelle

| Échelle | d_eff | Source |
|---|---|---|
| Simulations BuP (N=6–16) | 2.61 ± 0.03 | Simulations microscopiques |
| Galactique (SPARC, 175 galaxies) | 2.644 ± 0.295 | Courbes de rotation |
| Cosmologique (DESI+Pantheon+) | 2.40–2.71 | BAO + supernovæ |

Cohérence sur douze ordres de grandeur en échelle — non imposée.

---

## Liens

- **Paper 2** — HAL : [hal-05590614v1](https://hal.science/hal-05590614v1)
- **Paper 4** — Réduction dimensionnelle locale, prédiction JWST
- **GitHub** — branche `cosmology`

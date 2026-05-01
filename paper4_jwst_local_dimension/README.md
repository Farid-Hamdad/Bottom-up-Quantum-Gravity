# BuP Paper 4 — Effondrement non-linéaire par réduction dimensionnelle locale

**Titre :** Effondrement non-linéaire par réduction dimensionnelle locale
dans la Gravité Quantique Bottom-Up : Origine microscopique de l'excès JWST de galaxies massives

**Auteur :** Farid Hamdad
**ORCID :** [votre ORCID]
**Année :** 2026

---

## Résultat central

| Résultat | Valeur |
|---|---|
| Modèle local | d(z,δ) = d_bg(z) − ε·δ |
| ε_JWST (calibré) | 0.0113 |
| ε_micro(N=20) | 0.0065  IC=[0.0024, 0.0103] |
| ε > 0 robuste pour | N ≥ 16 (IC bootstrap strictement positif) |
| σ₈ préservé | 0.772 au premier ordre (⟨δ⟩=0) ✓ |

---

## Tableau FSS

| N | ε_micro | IC 68% | ε>0? |
|---|---|---|---|
| 9 | −0.003 | [−0.037, +0.031] | — |
| 12 | +0.005 | [−0.013, +0.023] | ∼ |
| 16 | +0.016 | [+0.008, +0.024] | ✓ |
| 18 | +0.009 | [−0.003, +0.021] | ∼ (marginal) |
| 20 | +0.007 | [+0.002, +0.010] | ✓ |

---

## Mécanisme physique

```
Régions denses → intrication plus forte → connectivité accrue
→ spectre Laplacien décalé → d_local < d_bg
→ G_eff_local > G_eff_bg → collapse accéléré → excès JWST

Compatibilité σ₈ : ⟨δ⟩ = 0 → ⟨d_local⟩ = d_bg → σ₈ intact
```

---

## Structure

```
paper/
├── main_fr.tex                  ← Article complet (FR)
├── fig_epsilon_fss.pdf          ← FSS ε(N=9→20)
└── fig_jwst_halo_excess.pdf     ← δ_c, G_eff_local, excès halos vs JWST

scripts/
├── bup_paper4_collapse.py       ← Collapse + calibration JWST
├── bup_generate_mi_matrices.py  ← Génération matrices MI BuP
└── bup_kill_epsilon_micro.py    ← Mesure microscopique de ε
```

---

## Démarrage rapide

```bash
pip install numpy scipy matplotlib

# Générer les matrices MI (N=9 à 20)
python scripts/bup_generate_mi_matrices.py \
    --N 9 12 16 18 20 --lambda-n 8 \
    --output-dir results_mi

# Mesurer ε microscopiquement
python scripts/bup_kill_epsilon_micro.py \
    --mi-files "results_mi/MI_*.csv" \
    --k 5 --tau-min 0.01 --tau-max 50 \
    --output-dir results_epsilon

# Collapse + prédiction JWST
python scripts/bup_paper4_collapse.py
```

---

## Liens

- **Paper 2** — HAL : [hal-05590614v1](https://hal.science/hal-05590614v1)
- **Paper 3** — Résolution de σ₈ via G_eff cohérent
- **GitHub** — branche `cosmology`

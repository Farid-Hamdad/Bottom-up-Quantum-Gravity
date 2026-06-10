# Paper 11 — Le Flot Dimensionnel de la Gravité BuP

## Dimension Spectrale, Dimension de Marche et l'Exposant Gravitationnel Effectif

Paper 11 — The Dimensional Flow of BuP Gravity
Spectral Dimension, Walk Dimension and the Effective Gravitational Exponent

**Auteur :** Farid Hamdad  
**Projet :** Bottom-Up Quantum Gravity (BuP)  
**Année :** 2026  

---

## Vue d'ensemble

Ce dossier contient le matériel numérique et théorique associé au **Paper 11** du programme Gravité Quantique Bottom-Up (BuP).

Paper 11 introduit la **loi dimensionnelle** qui relie les propriétés spectrales d'un graphe d'intrication à un **exposant gravitationnel effectif** :

$$
\alpha_{\rm eff} = \frac{2d_s}{d_w} + d_w - 4
$$

où :

- $d_s$ est la **dimension spectrale** du graphe d'intrication ;
- $d_w$ est la **dimension de marche** ;
- $\alpha_{\rm eff}$ est l'**exposant gravitationnel effectif** qui contrôle la réponse émergente à grande échelle.

Ce papier constitue le **pont** entre le programme microscopique du graphe et les tests phénoménologiques réalisés plus tard dans la séquence BuP.

---

## Rôle dans le programme BuP

Paper 11 fournit la **colonne vertébrale dimensionnelle** pour les papiers ultérieurs :
Paper 11 :
W_ij → L_ent → d_s, d_w → α_eff

Paper 12 :
Σ(R) → W_ij → L_ent → α_eff → V(r)

Paper 14 :
Courbes de rotation des galaxies SPARC et régimes LOW/HIGH

Paper 21 :
Point fixe de lentillage fort SLACS et α_eff ≈ 1

text

Ainsi, Paper 11 **n'ajuste aucun catalogue de galaxies**. Il établit la **loi théorique et numérique** que les papiers suivants testent observationnellement.

---

## Équation centrale

Le résultat fondamental est :

$$
\boxed{\alpha_{\rm eff} = \frac{2d_s}{d_w} + d_w - 4}
$$

Cette formule combine deux quantités issues de la **diffusion sur le graphe** :

$$
P(t) \sim t^{-d_s/2}
$$

$$
\langle r^2(t) \rangle \sim t^{2/d_w}
$$

L'exposant gravitationnel BuP n'est donc **pas postulé** : il est **déduit** de la diffusion sur le graphe d'intrication.

---

## Point fixe newtonien

Le point fixe newtonien (ou baryonique) correspond à :

$$
\alpha_{\rm eff} = 1
$$

Par conséquent :

$$
\frac{2d_s}{d_w} + d_w - 4 = 1
$$

De manière équivalente :

$$
2d_s = d_w (5 - d_w)
$$

$$
d_s = \frac{d_w (5 - d_w)}{2}
$$

Pour la **diffusion brownienne** standard :

$$
d_w = 2
$$

on obtient :

$$
d_s = 3
$$

Ainsi, le **comportement newtonien tridimensionnel ordinaire** apparaît comme un **point fixe particulier** du flot dimensionnel BuP.

---

## Principaux résultats numériques

**Fichiers de résultats :**
results/
alpha_eff_table.csv
finite_size_summary.csv
alpha_predictions_vs_N.csv
paper11_summary.json

text

**Figures :**
figures/
fig1_pipeline_dimensional_flow.png
fig2_ds_vs_N.png
fig3_dw_vs_N.png
fig4_alpha_predictions_vs_N.png
fig5_alpha_fixed_point_curve.png
fig6_interpretation_regimes.png

text

---

## Interprétation

Paper 11 montre que le comportement gravitationnel effectif est contrôlé par le couple :

$$
(d_s,\; d_w)
$$

Différents régimes correspondent à différentes réponses gravitationnelles :

| Régime | Condition | Interprétation |
|--------|-----------|----------------|
| Point fixe newtonien | $\alpha_{\rm eff} = 1$ | Comportement baryonique / newtonien standard |
| Sous-newtonien | $\alpha_{\rm eff} < 1$ | Gravité **affaiblie** ou sous-couplée |
| Super-newtonien | $\alpha_{\rm eff} > 1$ | Gravité **renforcée** |
| Régime de transition | $\alpha_{\rm eff} \approx 1$ | Crossover entre phases du graphe |

---

## Reproductibilité

Exécuter :

```bash
cd papers/paper11_dimensional_flow
bash scripts/run_all.sh
Cela génère toutes les tables numériques et les figures.

Statut
Paper 11 doit être lu comme une dérivation théorique et numérique de la loi dimensionnelle de BuP.
Ses conséquences observationnelles sont testées plus tard dans les Papers 12, 14 et 21.

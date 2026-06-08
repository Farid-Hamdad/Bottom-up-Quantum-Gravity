# Paper 21 — Le Point Fixe de SLACS

## Validation par Lentille Gravitationnelle Forte du Potentiel Effectif BuP

**Auteur :** Farid Hamdad  
**Projet :** Gravité Quantique Bottom-Up  
**Année :** 2026

---

## Vue d'ensemble

Ce dossier contient les tests numériques et observationnels associés au **Paper 21** du programme de Gravité Quantique Bottom-Up.

Le Paper 21 teste une prédiction introduite dans le **Paper 9** : le potentiel gravitationnel effectif BuP,

$$
\mathcal{L}_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux},
$$

doit laisser une empreinte observable dans les données de lentille gravitationnelle.

Le test est réalisé sur l'échantillon de lentille forte SLACS. L'observable principal est le résidu de correction de lentille

$$
C_{\rm obs} = \log\left(
\frac{\theta_E^{\rm obs}}
{\theta_E^{\rm baryon}}
\right).
$$

Le résultat central est la découverte d'une transition par point fixe autour de

$$
\log M_\star \simeq 11,58-11,60.
$$

À cette échelle :

1. le résidu observé vérifie $C_{\rm obs}\simeq 0$ ;
2. le potentiel $\Phi_{\rm BuP}$ performe localement mieux que le proxy global $\log M_\star$ ;
3. le secteur dimensionnel atteint $\alpha_{\rm eff}\simeq 1$ à une échelle de diffusion intermédiaire.

Ceci fournit un test de consistance observationnelle non trivial du potentiel effectif BuP.

---

## Origine physique

Le Paper 9 prédit la chaîne

$$
S_{\rm flux}
\rightarrow
\Phi_{\rm BuP}
\rightarrow
\text{réponse gravitationnelle}.
$$

Le Paper 21 teste cette chaîne à l'aide de données de lentille forte. L'exposant effectif BuP est

$$
\alpha_{\rm eff}
= \frac{2d_s}{d_w}+d_w-4.
$$

Le point fixe newtonien ou baryonique correspond à

$$
\alpha_{\rm eff}=1.
$$

Pour la diffusion brownienne standard ($d_w=2$), cette condition donne

$$
d_s=3.
$$

Plus généralement, la condition de point fixe est

$$
2d_s = d_w(5-d_w).
$$

Dans le graphe SLACS, le point fixe n'est pas retrouvé dans le régime de diffusion le plus précoce ni le plus tardif, mais dans une fenêtre de diffusion intermédiaire.

---

## Résultats principaux

### 1. Signal global du potentiel BuP

La meilleure caractéristique dynamique BuP atteint une amélioration leave-one-out d'environ

$$
10,49\%
$$

sur l'ensemble de l'échantillon de travail SLACS.

Un contrôle dynamique par permutation, où l'amplitude de la masse stellaire est mélangée avant de construire le graphe, supprime fortement le signal. Ceci montre que le résultat dépend de l'association correcte entre l'amplitude stellaire et la structure du graphe.

---

### 2. Les indices de Sérsic mesurés renforcent le signal

Une comparaison contrôlée a été réalisée sur les mêmes 61 galaxies dont les indices de Sérsic sont mesurés.

Avec les indices $n_i$ mesurés :

$$
\text{Amélioration LOO} = 10,40\%.
$$

Pour les mêmes galaxies, mais avec $n=4$ imposé :

$$
\text{Amélioration LOO} = 9,12\%.
$$

Ainsi, la morphologie photométrique mesurée augmente le signal dynamique BuP.

---

### 3. Fenêtre de transition

La fenêtre de transition la plus forte est

$$
11,545 < \log M_\star < 11,645.
$$

Dans cette fenêtre, avec les indices de Sérsic mesurés,

$$
\Phi_{\rm BuP}
$$

surpasse $\log M_\star$ pour environ

$$
81,8\%
$$

des galaxies.

Cela indique que le potentiel BuP n'est pas simplement un proxy global de la masse stellaire. Il capture une correction dynamique localisée dans le régime de transition.

---

### 4. Point fixe observationnel

Le zéro du résidu de lentille observé,

$$
C_{\rm obs}=0,
$$

se trouve autour de

$$
\log M_\star \simeq 11,58-11,60.
$$

Cela coïncide avec la fenêtre de masse où $\Phi_{\rm BuP}$ surpasse localement $\log M_\star$.

---

### 5. Point fixe dimensionnel

Un balayage de la fenêtre de diffusion montre que le point fixe dimensionnel BuP est retrouvé à une échelle de diffusion intermédiaire.

Pour les indices de Sérsic mesurés, la fenêtre de diffusion optimale est

$$
t_{\min}=1,\qquad t_{\max}=21.
$$

Elle donne

$$
\langle \alpha_{\rm eff}\rangle = 1,014,
$$

et dans la fenêtre de masse fixée

$$
11,545 < \log M_\star < 11,645
$$

elle donne

$$
\langle \alpha_{\rm eff}\rangle = 1,014.
$$

Ainsi, le point fixe observationnel et le point fixe dimensionnel BuP coïncident.

---

## Interprétation

Le Paper 21 confirme l'énoncé suivant :

$$
C_{\rm obs}=0,
\qquad
\Phi_{\rm BuP} > \log M_\star,
\qquad
\alpha_{\rm eff}\simeq 1
$$

se produisent tous autour de la même échelle de transition,

$$
\log M_\star \simeq 11,6.
$$

Ceci est interprété comme le **point fixe SLACS** du potentiel gravitationnel effectif BuP.

---

## Structure du répertoire

```text
paper21_slacs_fixed_point/
  README.md
  paper21_slacs_fixed_point.tex
  references.bib
  data/
  scripts/
  results/
  figures/

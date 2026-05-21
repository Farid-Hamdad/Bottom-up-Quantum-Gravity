# Paper 12 — Test SPARC de la loi de propagateur effectif BuP

## Titre

**Courbes de rotation SPARC à partir des profils d'intrication locale en Gravité Quantique Bottom-Up**

## Idée centrale

Paper 12 applique la loi de propagateur effectif dérivée dans Paper 11 aux courbes de rotation galactiques de l'échantillon SPARC.

La chaîne théorique est :

$$
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s(r), d_w(r))
\rightarrow
\alpha_{\rm eff}(r)
\rightarrow
\Phi_{\rm BuP}(r)
\rightarrow
V(r).
$$

Paper 11 a montré que l'exposant du propagateur effectif n'est pas contrôlé uniquement par la dimension spectrale $d_s$, mais par la paire $(d_s, d_w)$ :

$$
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w} + d_w - 4
}
$$

Paper 12 teste cette relation sur les courbes de rotation SPARC en :

1. construisant un graphe d'intrication baryonique pour chaque galaxie,
2. mesurant les profils locaux $(d_s(r), d_w(r))$,
3. dérivant $\alpha_{\rm eff}(r)$,
4. projetant le résultat en une contribution à la vitesse circulaire.

## Interprétation physique

BuP n'assume pas de halo de matière noire universel. Au lieu de cela, chaque configuration baryonique induit une géométrie d'intrication effective. La réponse gravitationnelle dépend de la structure locale de ce graphe.

En ce sens, les galaxies ne sont pas forcées dans un régime universel unique. Elles se répartissent en **classes de profils d'intrication**.

Les précédentes analyses SPARC suggéraient déjà une telle séparation à travers les profils empiriques $d_{\rm emp}(r)$ et $d_{\rm fit}(r)$. Paper 12 reformule cette séparation en termes de dimensions d'intrication locale :

$$
d(r)
\quad\longrightarrow\quad
(d_s(r), d_w(r))
\quad\longrightarrow\quad
\alpha_{\rm eff}(r).
$$

## Résultat principal

Le modèle corrigé v3.4 combine :

1. une fenêtre d'extraction locale $w \in \{5, 8, 10\}$ kpc,
2. un critère de transition relatif
   $$
   \alpha_c = 0.97 \, \alpha_{\max},
   $$
3. une correction de compacité pour les profils très piqués :
   $$
   r_t^{\rm corr} = r_t^{(0)} (1 - 0.75\eta),
   $$
   où
   $$
   \eta = \frac{v_{b,\max} - v_{b,\rm med}}{v_{b,\rm med}}.
   $$

Sur un échantillon étendu de 18 galaxies, le scan v3.4 par meilleure fenêtre donne :

| Métrique | Valeur |
|----------|--------|
| Nombre de galaxies | 18 |
| Amélioration vs baryons seuls | 17 / 18 |
| $\chi^2_{\rm red} < 2$ | 4 / 18 |
| $\chi^2_{\rm red} < 5$ | 8 / 18 |
| $\chi^2_{\rm red} < 10$ | 11 / 18 |
| Médiane $\chi^2_{\rm red}$ | 6.80 |
| Moyenne $\chi^2_{\rm red}$ | 14.10 |

En ajoutant l'analyse ciblée de compacité pour NGC5055, on obtient un résumé effectif à 19 galaxies :

| Métrique | Valeur |
|----------|--------|
| Nombre de galaxies | 19 |
| Amélioration vs baryons seuls | 18 / 19 |
| $\chi^2_{\rm red} < 2$ | 4 / 19 |
| $\chi^2_{\rm red} < 5$ | 8 / 19 |
| $\chi^2_{\rm red} < 10$ | 12 / 19 |
| Médiane $\chi^2_{\rm red}$ | 6.61 |
| Moyenne $\chi^2_{\rm red}$ | 13.65 |

## Galaxies clés

| Galaxie | Interprétation | Correction | $\chi^2_{\rm red}$ |
|---------|----------------|------------|--------------------|
| NGC3198 | transition régulière canonique | $w=5$, $q=0.97$ | 0.994 |
| NGC2841 | disque étendu standard | $w=10$ | 4.586 |
| NGC5055 | système complexe (barre/bulbe) | $r_t$ corrigé par $\eta$ | 5.526 |
| NGC2403 | régime sous-newtonien précoce diffus | non résolu par sigmoïde simple | 115.8 |

## Statut scientifique

Paper 12 n'est **pas** présenté comme un modèle SPARC universel définitif. C'est une **preuve** que la loi de propagateur effectif issue de Paper 11 a un contenu observationnel et que les courbes de rotation galactiques peuvent être organisées par des profils d'intrication locale.

La conclusion la plus forte est :

$$
\boxed{
\text{Les galaxies SPARC se divisent en classes de profils d'intrication.}
}
$$

La prochaine étape consiste à remplacer la projection sigmoïde externe unique par des projections dépendantes de la classe ou multi-échelles, en particulier pour les systèmes diffus comme NGC2403.

---

# Résultats — Paper 12 : test SPARC de la loi $\alpha(d_s,d_w)$

Ce dossier contient les résultats numériques pour Paper 12.

## Fichiers de résultats principaux

| Fichier / Dossier | Description |
|---|---|
| `paper12_unified_sparc_classification.csv` | Table unifiée fusionnant les anciennes analyses SPARC $d(r)$ avec les résultats de profils alpha de Paper 12 |
| `paper12_unified_sparc_classification.md` | Version lisible en Markdown de la classification unifiée |
| `classe1_extended_window_scan_best_summary.csv` | Scan v3.4 par meilleure fenêtre sur 18 galaxies |
| `NGC5055_v3_4_window_scan_eta_summary.csv` | Scan de compacité v3.4 dédié pour NGC5055 |
| `NGC3198_v3_3_q097_dr2p5/` | Résultat canonique de NGC3198 avec déclenchement relatif |
| `NGC2841_v3_4_corrected/` | Correction par fenêtre adaptative pour NGC2841 |
| `NGC5055_v3_4_corrected/` | Résultat corrigé par compacité pour NGC5055 |

## Résumé du scan v3.4 par meilleure fenêtre

Le scan par meilleure fenêtre utilise :

$$
w \in \{5, 8, 10\}\ \text{kpc},\qquad q = 0.97,\qquad \Delta r = 2.5\ \text{kpc}.
$$

Résultats sur 18 galaxies :

| Métrique | Valeur |
|---|---|
| Galaxies | 18 |
| Amélioration vs baryons seuls | 17 |
| $\chi^2_{\rm red} < 2$ | 4 |
| $\chi^2_{\rm red} < 5$ | 8 |
| $\chi^2_{\rm red} < 10$ | 11 |
| Médiane $\chi^2_{\rm red}$ | 6.80 |
| Moyenne $\chi^2_{\rm red}$ | 14.10 |

## Résultat ciblé pour NGC5055

NGC5055 nécessite la correction de compacité :

$$
r_t^{\rm corr} = r_t^{(0)} (1 - 0.75\eta).
$$

Meilleur résultat obtenu :

| Fenêtre | $\alpha_{\max}$ | $\eta$ | $r_t$ (kpc) | $\chi^2_{\rm red}$ |
|---------|----------------|--------|-------------|--------------------|
| 5 kpc | 12.60 | 0.337 | 12.22 | 5.526 |

Les fenêtres plus larges effacent le pic de $\alpha_{\max}$ et suppriment la correction $\eta$, ce qui dégrade l'ajustement. Par conséquent, pour NGC5055, le pic élevé de $\alpha_{\max}$ est interprété comme un **signal physique de complexité interne**, et non comme du bruit numérique.

## Interprétation

Les résultats supportent quatre régimes :

| Classe | Signification | Exemple |
|--------|---------------|---------|
| I-A | transition externe régulière | NGC3198 |
| I-B | disque étendu standard, fenêtre locale plus large | NGC2841 |
| II | structure interne complexe, correction de compacité | NGC5055 |
| III | régime sous-newtonien précoce diffus | NGC2403 |

---

# Scripts — Paper 12

## Scripts principaux

| Script | Objectif |
|--------|----------|
| `bup_paper12_sparc_alpha_dsdw_v3_1_radial_outer.py` | Première projection réussie avec plateau externe |
| `bup_paper12_sparc_alpha_dsdw_v3_3_relative_trigger.py` | Déclenchement de transition relatif $\alpha_c = q \alpha_{\max}$ |
| `bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py` | Version avec fenêtre adaptative et correction de compacité |
| `paper12_merge_old_sparc_with_v33.py` | Fusionne les anciennes analyses SPARC $d(r)$ avec les profils Paper 12 |

## Modèle v3.4

Le modèle v3.4 utilise :

$$
\alpha_{\rm eff}(r) = \frac{2d_s(r)}{d_w(r)} + d_w(r) - 4.
$$

Le rayon de transition est dérivé de :

$$
\alpha_{\rm eff}(r_t) = q \alpha_{\max}, \qquad q = 0.97.
$$

Pour les galaxies à complexité interne, la correction de compacité est :

$$
r_t^{\rm corr} = r_t^{(0)} (1 - 0.75\eta),
$$

avec :

$$
\eta = \frac{v_{b,\max} - v_{b,\rm med}}{v_{b,\rm med}}.
$$

## Exemple d'exécution : NGC3198

```bash
python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/chemin/vers/SPARC/Rotmod_LTG 2/NGC3198_rotmod.dat" \
  --galaxy NGC3198 \
  --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
  --n-windows 8 --window-width 5 \
  --min-nodes-window 80 \
  --ups-disk 0.5 --ups-bulge 0.7 \
  --shape-mode flat_outer \
  --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
  --dr-fixed 2.5 --rt-mode alpha_crossing \
  --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
  --output-dir results/NGC3198_v3_4

  python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/chemin/vers/SPARC/Rotmod_LTG 2/NGC5055_rotmod.dat" \
  --galaxy NGC5055 \
  --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
  --n-windows 8 --window-width 5 \
  --min-nodes-window 80 \
  --ups-disk 0.5 --ups-bulge 0.7 \
  --shape-mode flat_outer \
  --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
  --dr-fixed 2.5 --rt-mode alpha_crossing \
  --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
  --output-dir results/NGC5055_v3_4_corrected

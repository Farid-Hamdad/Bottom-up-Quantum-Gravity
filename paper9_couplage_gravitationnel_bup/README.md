# BuP Paper 9 — Couplage gravitationnel de la matière émergente

**Titre FR :** Couplage gravitationnel de la matière émergente : de S_flux au potentiel effectif BuP  
**Titre EN :** Gravitational Coupling of Emergent Matter: from S_flux to the BuP Effective Potential  
**Auteur :** Farid Hamdad  
**Année :** 2026

---

## Idée centrale

Paper 8 a construit une source tensorielle matière candidate à partir d’une excitation locale du réseau d’intrication :

$$
\delta W^{loc} \longrightarrow S_{flux}.
$$

Paper 9 couple cette source à une équation de Poisson discrète sur le graphe d’intrication :

$$
L_{ent} \Phi_{BuP} = S_{flux}.
$$

La chaîne physique testée devient :

$$
\delta W^{loc} \longrightarrow S_{flux} \longrightarrow \Phi_{BuP} \longrightarrow |\delta R|.
$$

L’objectif est de vérifier si la source matière émergente produit un potentiel gravitationnel effectif cohérent avec la réponse de courbure.

---

## Source utilisée

La source reprise de Paper 8 est :

$$
S_{flux} = T_{00} - \frac12 T_{aa} + \frac12 T_{grad} + \|T_{0a}\|.
$$

Elle combine :

- densité d’énergie informationnelle \(T_{00}\)
- trace spatiale / pression effective \(T_{aa}\)
- contrainte interne \(T_{grad}\)
- flux informationnel \(T_{0a}\)

---

## Équation de Poisson discrète

Dans la limite faible effective, nous résolvons :

$$
L_{ent} \Phi_{BuP} = S_{flux}.
$$

Ici \(L_{ent}\) est le Laplacien du graphe d’intrication. Le potentiel \(\Phi_{BuP}\) est ensuite comparé à la réponse de courbure nodale :

$$
|\delta R_i| = |R_i[W^{(0)}+\delta W] - R_i[W^{(0)}]|.
$$

---

## Résultat principal

Pour :

$$
N = 20,\qquad \lambda = 0.57,\qquad k = 5,\qquad \sigma = 0.15,
$$

nous obtenons :

$$
\rho_{Spearman}(S_{flux},|\delta R|) = 0.741,\qquad p = 1.84 \times 10^{-4}.
$$

Après résolution de l’équation de Poisson discrète :

$$
L_{ent}\Phi_{BuP} = S_{flux},
$$

le potentiel vérifie :

$$
\rho_{Spearman}(\Phi_{BuP},|\delta R|) = 0.738,\qquad p = 2.01 \times 10^{-4}.
$$

Ainsi, le potentiel effectif conserve presque toute l’information géométrique portée par la source.

---

## Scan en largeur \(\sigma\)

| \(\sigma\) | \(S_{flux}\) vs \(|\delta R|\) | \(\Phi_{BuP}\) vs \(|\delta R|\) | Lecture |
|-----------|-------------------------------|----------------------------------|---------|
| 0.01 | 0.504 | 0.498 | signal modéré |
| 0.02 | 0.473 | 0.479 | signal modéré |
| 0.05 | 0.148 | 0.115 | non significatif |
| 0.08 | 0.176 | 0.157 | non significatif |
| 0.10 | 0.122 | 0.199 | non significatif |
| 0.15 | 0.741 | 0.738 | meilleur régime |

Le couplage gravitationnel est donc maximal pour une excitation suffisamment étendue :

$$
\sigma = 0.15.
$$

Cela suggère que la gravité émergente BuP est un phénomène collectif du réseau d’intrication.

---

## Stabilité au shift et à la normalisation

Une objection possible concerne le remplacement de la source brute par une source positive ou normalisée. Nous avons donc testé trois variantes :

$$
S_{\rm raw}=S_{\rm flux},
$$

$$
S_{\rm shift}
=
S_{\rm flux}
-
\min_i S_{\rm flux}(i),
$$

et

$$
S_{\rm norm}
=
\frac{S_{\rm shift}}{\sum_i S_{\rm shift}(i)}.
$$

Pour chaque variante, nous résolvons :

$$
\Phi
=
L_{\rm ent}^{+}S.
$$

Le test est effectué sur la matrice MI réelle `MI_N20_lam0.57.csv`, avec une réponse de courbure nodale issue de la courbure d’Ollivier-Ricci.

**Résultat :**

$$
\max_{\sigma} \Delta\rho_{\rm shift} = 0,
\qquad
\max_{\sigma} \Delta\rho_{\rm norm} = 0.
$$

Autrement dit :

$$
\rho(\Phi_{\rm raw},|\delta R|)
=
\rho(\Phi_{\rm shift},|\delta R|)
=
\rho(\Phi_{\rm norm},|\delta R|)
$$

pour toutes les valeurs de largeur testées.

Ce résultat confirme que le signal gravitationnel ne vient pas du shift ou de la normalisation de la source. Il vient de la projection de la partie non constante de \(S_{\rm flux}\) par la pseudo-inverse du Laplacien :

$$
\Phi_{\rm BuP}
=
L_{\rm ent}^{+}S_{\rm flux}.
$$

Cela renforce l’interprétation suivante :

$$
G_{ij}^{\rm ent}
=
(L_{\rm ent}^{+})_{ij}
$$

est un candidat naturel au propagateur gravitationnel discret du graphe d’intrication.

---

## Interprétation physique

Le résultat central est :

$$
\boxed{S_{flux} \text{ génère un potentiel effectif } \Phi_{BuP} \text{ fortement corrélé à la réponse de courbure.}}
$$

Cela établit un premier pont physique entre :

$$
\text{matière émergente}
\longrightarrow
\text{potentiel gravitationnel}
\longrightarrow
\text{courbure émergente}.
$$

Paper 9 montre que BuP ne produit pas seulement une corrélation locale entre source et courbure : la source peut être propagée par un opérateur de Poisson sur le graphe, comme dans une limite gravitationnelle effective.

---

## Prédictions

Paper 9 propose plusieurs prédictions internes au cadre BuP :

1. Le potentiel \(\Phi_{BuP}\) doit être corrélé à la réponse de courbure \(|\delta R|\).
2. La source gravitationnelle dépend de \(T_{00}\), \(T_{aa}\), \(T_{grad}\) et \(T_{0a}\), pas seulement de la densité.
3. Le couplage gravitationnel est maximal pour une excitation étendue.
4. Le gradient \(|\nabla\Phi_{BuP}|\) peut former une structure périphérique autour de la source.
5. À terme, \(\Phi_{BuP}\) peut être utilisé pour prédire des profils de rotation ou de lentille gravitationnelle.

---

## Structure

```text
paper9_couplage_gravitationnel_bup/
├── README.md
├── paper/
│   └── main_fr.tex
├── scripts/
│   ├── README.md
│   ├── bup_paper9_poisson_effective_potential_v1.py
│   └── bup_paper9_shift_stability_real_mi_v1.py
├── results/
│   ├── README.md
│   ├── paper9_sigma_scan_summary.csv
│   ├── paper9_poisson_summary.csv
│   ├── shift_stability_real_mi_summary.csv
│   ├── shift_stability_real_mi_summary.json
│   └── summary.json
└── figures/
    ├── fig_sigma_scan_correlations.png
    ├── fig_sigma_scan_phi.png
    ├── fig_sigma_scan_pvalues.png
    ├── fig_source_vs_curvature.png
    ├── fig_potential_vs_curvature.png
    ├── fig_acceleration_vs_curvature.png
    ├── fig_phi_profile.png
    ├── fig_acceleration_profile.png
    ├── fig_shift_stability_real_mi_correlations.png
    ├── fig_shift_stability_real_mi_delta_rho.png
    └── fig_shift_stability_real_mi_phi_difference.png

Citation
bibtex
@misc{hamdad2026paper9,
  author  = {Farid Hamdad},
  title   = {Couplage gravitationnel de la matière émergente :
             de S_flux au potentiel effectif BuP},
  year    = {2026},
  note    = {BuP Paper 9},
  url     = {https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity}
}
Liens
GitHub : Farid-Hamdad/Bottom-up-Quantum-Gravity

Paper 5 : paper5_effective_einstein/
Paper 6 : paper6_pure_bup_equation/
Paper 7 : paper7_variational_dynamics/
Paper 8 : paper8_matiere_emergente_bup/
Paper 9 (source) : paper9_couplage_gravitationnel_bup/paper/main_fr.tex




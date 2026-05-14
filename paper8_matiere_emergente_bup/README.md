# BuP Paper 8 — De l’excitation d’intrication à la matière émergente

**Titre FR :** De l’excitation d’intrication à la matière émergente : construction d’une source tensorielle BuP  
**Titre EN :** From Entanglement Excitations to Emergent Matter: Building a BuP Tensorial Source  
**Auteur :** Farid Hamdad  
**Année :** 2026

---

## Idée centrale

Paper 8 construit une source tensorielle candidate de matière à partir d’une excitation locale du réseau d’information mutuelle :

$$
W_{ij} = W_{ij}^{(0)} + \delta W_{ij}^{\rm loc}.
$$

La chaîne testée est :

$$
\delta W^{\rm loc} \longrightarrow T_{\mu\nu}^{\rm matter,cand} \longrightarrow \delta R.
$$

Le résultat principal est que la densité d’énergie $T_{00}$ seule est positive, normalisée et localisée, mais insuffisante pour prédire la réponse de courbure. Il faut ajouter les composantes spatiales, internes et de flux.

---

## Densité minimale

$$
T_{00}^{\rm matter,cand}(i) = \frac12 \sum_j (\delta W_{ij}^{\rm loc})^2.
$$

Elle vérifie :

$$
T_{00}^{\rm matter,cand}(i) \ge 0,
$$

et :

$$
\sum_i T_{00}^{\rm matter,cand}(i) = \sum_{i &lt; j} (\delta W_{ij}^{\rm loc})^2.
$$

---

## Source complète testée

La meilleure source candidate obtenue numériquement est :

$$
S_{\rm flux} = T_{00} - \frac12 T_{aa} + \frac12 T_{\rm grad} + \|T_{0a}\|.
$$

Les coefficients sont des poids effectifs dans un protocole normalisé. Ils ne doivent pas encore être interprétés comme des constantes fondamentales.

---

## Hiérarchie des résultats

| Source testée | Meilleurs paramètres | Spearman avec $|\delta R|$ | $p$-value | Lecture |
|---|---:|---:|---:|:---|
| $T_{00}$ | $\sigma=0.05$ | 0.195 | 0.409 | densité seule insuffisante |
| $T_{00} - \frac12 T_{aa}$ | $\sigma=0.15$ | 0.526 | 0.017 | source fluide effective |
| $T_{00} - \frac12 T_{aa} + \frac12 T_{\rm grad}$ | $\sigma=0.15$ | 0.612 | 0.00413 | contrainte interne utile |
| $T_{00} - \frac12 T_{aa} + \frac12 T_{\rm grad} + \|T_{0a}\|$ | $\sigma=0.15$ | 0.741 | $1.84 \times 10^{-4}$ | meilleure source complète |

---

## Test de conservation

La source matière candidate n’est pas conservée isolément. Le test spatial montre une non-conservation significative du secteur matière seul :

$$
\nabla_\mu T^{\mu\nu}_{\rm matter,cand} \neq 0.
$$

Dans BuP, cela n’est pas nécessairement un échec. La matière est une excitation du réseau d’intrication ; elle peut échanger énergie et impulsion informationnelles avec le fond géométrique.

La conservation attendue est plutôt :

$$
\nabla_\mu \left( T^{\mu\nu}_{\rm matter} + T^{\mu\nu}_{\rm ent} \right) \simeq 0.
$$

Un flux d’intrication construit avec :

$$
\phi_{\rm ent} = \delta d_s + \delta R + \frac12 \Delta\delta R
$$

réduit le résidu de divergence d’un facteur $2.63$ pour $\sigma = 0.01$ :

$$
\frac{ \|\mathrm{div}\,f_{\rm total}\|_2 }{ \|\mathrm{div}\,f_{\rm matter}\|_2 } = 0.380.
$$

Cela indique une compensation partielle matière‑intrication.

---

## Interprétation

Le résultat principal de Paper 8 est :

$$
\boxed{\text{la matière BuP est mieux décrite comme une source tensorielle candidate non fermée.}}
$$

Elle est fortement corrélée à la réponse de courbure locale, mais elle n’est pas conservée isolément. Sa non-conservation devient une signature d’échange avec le champ d’intrication.

---

## Structure

```text
paper8_matiere_emergente_bup/
├── README.md
├── paper/
│   └── main_fr.tex
├── scripts/
│   ├── README.md
│   ├── bup_paper8_Tmatter_positivity_v1.py
│   ├── bup_paper8_Tmatter_localization_v1.py
│   ├── bup_paper8_Tmatter_curvature_response_v1.py
│   ├── bup_paper8_Tmatter_fluid_source_v1.py
│   ├── bup_paper8_Tmatter_grad_source_v1.py
│   ├── bup_paper8_Tmatter_flux_source_v1.py
│   ├── bup_paper8_Tmatter_conservation_v1.py
│   ├── bup_paper8_total_conservation_test_v1.py
│   └── bup_paper8_total_conservation_test_v2.py
├── results/
│   ├── README.md
│   ├── paper8_summary.csv
│   ├── localization_sigma_scan.csv
│   ├── curvature_response_sigma_scan.csv
│   ├── fluid_source_scan.csv
│   ├── grad_source_chi_scan.csv
│   ├── flux_source_psi_scan.csv
│   ├── conservation_sigma_scan.csv
│   └── total_conservation_v2_sigma_summary.csv
└── figures/
    ├── fig_source_hierarchy.png
    ├── fig_flux_response.png
    ├── fig_total_conservation.png
    └── fig_best_source_vs_deltaR.png
Citation
bibtex
@misc{hamdad2026paper8,
  author  = {Farid Hamdad},
  title   = {De l'excitation d'intrication à la matière émergente :
             construction d'une source tensorielle BuP},
  year    = {2026},
  note    = {BuP Paper 8},
  url     = {https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity}
}
Liens
GitHub : Farid-Hamdad/Bottom-up-Quantum-Gravity

Paper 5 : paper5_effective_einstein/

Paper 6 : paper6_pure_bup_equation/

Paper 7 : paper7_variational_dynamics/

Paper 8 (source) : paper8_matiere_emergente_bup/paper/main_fr.tex

# BuP Paper 7 — Étude de la dynamique des excitations du réseau d’intrication BuP

**Titre FR :** Dynamique variationnelle des excitations BuP : métastabilité, phases internes et relaxation de la matière émergente  
**Titre EN :** Variational Dynamics of BuP Excitations: Metastability, Internal Phases and Relaxation of Emergent Matter  
**Auteur :** Farid Hamdad  
**Année :** 2026
---

## Idée centrale

Paper 6 a proposé l’équation pure BuP :

\[
\mathcal{F}_{ij}[W]=0,
\]

ou sous forme variationnelle :

\[
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0.
\]

Paper 7 étudie la dynamique induite par une action candidate \(S_{\rm BuP}[W]\).  
L’objectif est de comprendre comment une excitation localisée du réseau d’intrication évolue sous relaxation variationnelle.

Une excitation localisée est écrite :

\[
W_{ij}
=
W_{ij}^{(0)}
+
\delta W_{ij}^{\rm loc}.
\]

Le point central du papier est que le caractère “source de matière” d’une excitation n’est pas automatiquement conservé : il dépend de la dynamique de relaxation de \(S_{\rm BuP}\).

---

## Observable matière-courbure

Le proxy de matière est :

\[
T_{ij}^{\rm matter,proxy}
=
(\delta W_{ij})^2.
\]

Le proxy géométrique est :

\[
G_{ij}^{\rm proxy}
=
-\kappa_{ij}
+
\frac12 R_{\rm edge}.
\]

L’observable principale est :

\[
\mathcal{O}_{\rm matter}
=
\rho_{\rm Spearman}
\left(
T_{ij}^{\rm matter,proxy},
|\delta G_{ij}^{\rm proxy}|
\right).
\]

Interprétation :

- \( \mathcal{O}_{\rm matter} > 0 \) : l’excitation agit comme une source géométrique effective ;
- \( \mathcal{O}_{\rm matter} < 0 \) : la réponse devient anti-corrélée, donc l’excitation perd son caractère source-like.

---

## Résultat principal

Le scan en \(\lambda\) révèle une dynamique non triviale sous relaxation de \(S_{\rm BuP}\).

Trois régimes apparaissent :

\[
\lambda_{\rm init}
\simeq
0.265,
\]

\[
\lambda_{\rm meta}
\simeq
0.876,
\]

\[
\lambda_{\rm late}
\simeq
0.280.
\]

### Interprétation

| Régime | Valeur | Interprétation |
|---|---:|---|
| Initial | \( \lambda_{\rm init}\simeq0.265 \) | crossing juste après injection |
| Métastable | \( \lambda_{\rm meta}\simeq0.876 \) | plateau source-like observé entre 10 et 20 pas |
| Tardif | \( \lambda_{\rm late}\simeq0.280 \) | régime non-source stable de 40 à 320 pas |

Ce résultat suggère que les excitations BuP peuvent passer par un régime métastable source-like avant de relaxer vers un régime tardif non-source.

---

## Résultat temporel

La dépendance au nombre de pas de relaxation donne :

| Steps | Crossing final | Lecture |
|---:|---:|---|
| 5 | \(0.215\) | régime transitoire court |
| 10–20 | \(0.876\) | plateau métastable source-like |
| 40–320 | \(0.280\) | régime tardif non-source |

Le régime tardif :

\[
\lambda_{\rm late}
=
0.2800098754
\]

reste stable pour :

\[
40,\;80,\;120,\;160,\;320
\]

pas de relaxation.

---

## Interprétation physique prudente

Les résultats suggèrent que la matière émergente BuP pourrait être comprise comme une phase métastable du réseau d’intrication.

Une particule BuP ne serait pas simplement une perturbation locale de \(W_{ij}\), mais un défaut métastable du réseau :

\[
\text{particule BuP}
=
\text{excitation localisée source-like métastable}.
\]

Dans cette lecture :

- une excitation source-like correspond à une phase de matière effective ;
- une relaxation longue peut conduire vers un régime non-source ;
- la stabilité d’une particule nécessiterait une protection interne : topologique, spectrale ou symétrique.

Cette interprétation reste programmatique. Le temps de relaxation numérique n’est pas encore identifié au temps physique cosmologique.

---

## Action candidate

L’action utilisée dans ce papier est une action effective candidate :

\[
S_{\rm BuP}[W]
=
S_{\rm spec}
+
S_{\rm geom}
+
S_{\rm loc}
+
S_{\rm norm}
+
S_{\rm exc}.
\]

Elle contient :

- un terme spectral \(S_{\rm spec}\), contrôlant \(d_s\) ;
- un terme géométrique \(S_{\rm geom}\), basé sur \(G_{ij}^{\rm proxy}\) ;
- un terme de localité \(S_{\rm loc}\), basé sur \(d_{ij}^{\rm ent}\) ;
- un terme de normalisation \(S_{\rm norm}\) ;
- un terme d’excitation \(S_{\rm exc}\).

Cette action n’est pas encore l’action finale de BuP. Elle sert ici à explorer la dynamique variationnelle des excitations.

---

## Robustesse

Le plateau métastable autour de :

\[
\lambda_{\rm meta}
\simeq
0.876
\]

est robuste dans les scans effectués sur les coefficients :

\[
\gamma=0.005,\;0.01,\;0.02,\;0.05,
\]

\[
\beta=0.02,\;0.05,\;0.1,\;0.2,\;0.5,
\]

\[
m_{\rm exc}=0,\;0.01,\;0.05,\;0.1,\;0.2.
\]

Cette robustesse est observée dans le protocole numérique utilisé. Elle ne doit pas encore être interprétée comme une universalité fondamentale.

---

## Structure du dossier

```text
paper7_variational_dynamics/
├── README.md
├── paper/
│   └── main_fr.tex
├── figures/
│   ├── fig_phase_scan_lambda.png
│   ├── fig_delta_O_lambda.png
│   ├── fig_internal_observables_lambda.png
│   └── fig_crossing_vs_steps.png
├── scripts/
│   └── bup_paper7_phase_scan_v1.py
└── results/
    ├── paper7_phase_scan_summary.csv
    ├── paper7_steps_summary.csv
    └── summary.json

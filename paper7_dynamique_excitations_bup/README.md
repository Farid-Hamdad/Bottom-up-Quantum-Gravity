# BuP Paper 7 — Étude de la dynamique des excitations du réseau d’intrication BuP

**Titre FR :** Dynamique variationnelle des excitations BuP : métastabilité, phases internes et relaxation de la matière émergente  
**Titre EN :** Variational Dynamics of BuP Excitations: Metastability, Internal Phases and Relaxation of Emergent Matter  
**Auteur :** Farid Hamdad  
**Année :** 2026  

---

## Idée centrale

Paper 7 étudie la dynamique d'une excitation locale du réseau d'intrication sous relaxation d'une action candidate :

$$
S_{\rm BuP}[W].
$$

L'objectif est de comprendre si une excitation localisée $\delta W^{\rm loc}$ peut se comporter comme une source géométrique stable, métastable ou non-source.

L'observable principale est :

$$
\mathcal{O}_{\rm matter} = \rho_{\rm Spearman} \left( T_{ij}^{\rm matter,proxy}, |\delta G_{ij}^{\rm proxy}| \right),
$$

où :

$$
T_{ij}^{\rm matter,proxy} = (\delta W_{ij}^{\rm loc})^2.
$$

---

## Résultat principal

Pour $N = 20$, le système présente trois régimes :

$$
\lambda_{\rm init} \simeq 0.265,
$$

$$
\lambda_{\rm meta} \simeq 0.876,
$$

$$
\lambda_{\rm late} \simeq 0.280.
$$

Le régime autour de $\lambda_{\rm meta} \simeq 0.876$ apparaît pour des temps de relaxation intermédiaires, entre 10 et 20 pas. Le régime tardif autour de $\lambda_{\rm late} \simeq 0.280$ apparaît à partir de 40 pas et reste stable jusqu'à 320 pas.

---

## Interprétation

La relaxation de $S_{\rm BuP}$ révèle une structure de paysage variationnel :

- un état source-like **métastable** ;
- un régime tardif **non-source** ;
- une dépendance forte au temps de relaxation ;
- une dépendance non triviale à la taille $N$.

La matière émergente peut donc être interprétée comme une excitation source-like métastable du réseau d'intrication.

---

## Scan en temps de relaxation

| Steps | Crossing final | Interprétation |
|------:|---------------:|----------------|
| 5 | 0.215 | régime transitoire court |
| 10 | 0.876 | plateau métastable |
| 20 | 0.876 | plateau métastable confirmé |
| 40 | 0.280 | régime tardif |
| 80 | 0.280 | régime tardif stable |
| 120 | 0.280 | régime tardif stable |
| 160 | 0.280 | régime tardif stable |
| 320 | 0.280 | régime tardif stable |

---

## Scan en taille finie

| $N$ | Steps | $\lambda_{\rm init}$ | $\lambda_{\rm final}$ | Interprétation |
|----:|------:|---------------------:|----------------------:|----------------|
| 15 | 20 | 0.042 | 0.042 | pas de dynamique visible |
| 15 | 160 | 0.042 | 0.042 | pas de dynamique visible |
| 18 | 20 | 0.522 | 0.478 | régime précritique |
| 18 | 160 | 0.522 | 0.478 | régime précritique stable |
| 20 | 20 | 0.265 | 0.876 | plateau métastable source-like |
| 20 | 160 | 0.265 | 0.280 | régime tardif non-source |
| 22 | 20 | none | 0.359 | phase source-like initialement dominante |
| 22 | 160 | none | 0.360 | seuil dynamique stable |
| 24 | 20 | 0.427 | none | phase source-like globale après relaxation |
| 24 | 160 | 0.427 | 0.521 | sélection tardive |

---

## Conclusion courte

Paper 7 montre que les excitations BuP ne sont pas simplement des perturbations statiques de $W_{ij}$. Elles possèdent une dynamique interne, des phases métastables et des effets de taille finie. Une particule BuP peut être vue comme un défaut source-like métastable du réseau d'intrication.

---

## Structure

```text
paper7_variational_dynamics/
├── README.md
├── paper/
│   └── main_fr.tex
├── scripts/
│   ├── bup_paper7_phase_scan_v1.py
│   └── bup_generate_mi_matrices.py
├── results/
│   ├── paper7_phase_scan_summary.csv
│   ├── paper7_steps_summary.csv
│   ├── paper7_finite_size_summary.csv
│   └── summary.json
└── figures/
    ├── fig_phase_scan_lambda.png
    ├── fig_crossing_vs_steps.png
    └── fig_finite_size_crossings.png

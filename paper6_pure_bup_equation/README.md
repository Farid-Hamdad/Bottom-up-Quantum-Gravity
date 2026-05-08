# BuP Paper 6 — Équation pure BuP

**Titre FR :** Équation pure BuP : de l'intrication à l'équation d'Einstein effective  
**Titre EN :** Pure BuP Equation: from Entanglement to the Effective Einstein Equation  
**Auteur :** Farid Hamdad  
**Année :** 2026

---

## Idée centrale

Paper 5 a formulé l'équation d'Einstein effective de BuP :

$$
G_{\mu\nu}[g^{\rm ent}] = 8\pi G_{\rm eff}(d_s) T_{\mu\nu}^{\rm matter} + 8\pi G T_{\mu\nu}^{\rm ent}[d_s].
$$

Paper 6 propose l'équation pure BuP sous-jacente :

$$
\mathcal{F}_{ij}[W] = 0,
$$

où

$$
W_{ij} = I(i:j)
$$

est le réseau d'information mutuelle.

Dans cette formulation, la géométrie et la matière ne sont pas fondamentales. Elles émergent comme limites continues / grossières du réseau d'intrication.

---

## Formulation variationnelle

L'équation pure BuP peut s'écrire sous forme variationnelle :

$$
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}} = 0,
$$

où \( S_{\rm BuP}[W] \) est une action conjecturale (termes spectraux, de courbure, de localité, topologiques). Sa forme exacte n'est pas encore fixée et fait partie du programme futur.

---

## Chaîne d'émergence

$$
W_{ij}
\rightarrow
d_{ij}^{\rm ent}
\rightarrow
g_{\mu\nu}^{\rm ent}
\rightarrow
d_s(x)
\rightarrow
G_{\mu\nu}[g^{\rm ent}]
$$

et

$$
\delta W_{ij}^{\rm loc}
\rightarrow
T_{\mu\nu}^{\rm matter}.
$$

---

## Test numérique de cohérence

On injecte une excitation localisée :

$$
W_{ij} = W_{ij}^{(0)} + \delta W_{ij}^{\rm loc}.
$$

Le proxy de matière est :

$$
T_{ij}^{\rm matter,proxy} = (\delta W_{ij})^2.
$$

Le paramètre d'ordre est :

$$
\mathcal{O}_{\rm matter}(\lambda)
=
\rho_{\rm Spearman}
\left(
T_{ij}^{\rm matter,proxy},
|\delta G_{ij}^{\rm proxy}|
\right).
$$

Pour \(N = 20\), le scan fin donne :

$$
\lambda_c(N=20) \simeq 0.593.
$$

Ce seuil est interprété comme une estimation de taille finie d'un seuil critique de localité :

- Pour \( \lambda < \lambda_c \), l'excitation localisée agit comme une source positive de courbure, c'est-à-dire comme un proxy de matière effective.
- Pour \( \lambda > \lambda_c \), la réponse devient anti-corrélée, indiquant une perte de cohérence locale.

Ainsi, dans ce test, les excitations interprétables comme sources matérielles effectives n'apparaissent que dans la **phase locale** du réseau d'intrication.

---

## Aperçu exploratoire de la relaxation d'action

La relaxation préliminaire d'une action candidate \( S_{\rm BuP}^{(v2)} \) est présentée comme un prototype exploratoire, non comme une validation complète de l'action pure BuP. Cette relaxation préliminaire (qui conserve la masse et le rayon de l'excitation) montre que la cohérence géométrique \( \mathcal{O}_{\rm matter} \) diminue même lorsque la masse et le rayon sont fixés. Cela indique que la stabilité de la matière nécessite non seulement une localisation, mais aussi une **cohérence interne**, peut-être protégée par des invariants topologiques ou spectraux.

---

## Structure du dépôt

```text
paper/
└── main_fr.tex

figures/
├── fig_lambda_c_fine_positive.png
├── fig_lambda_c_fine_negative.png
└── fig_FSS_lambda_c.png

scripts/
├── bup_pure_matter_excitation_v2.py   # Mesure de O_matter(λ)
├── bup_lambda_critical_scan_v1.py     # Scan fin du seuil λ_c
└── bup_action_sim_v2.py               # Relaxation avec contraintes masse & rayon

results/
├── lambda_critical_scan_fine_positive/
├── lambda_critical_scan_fine_negative/
└── fss_lambda_c/

Citation
bibtex
@misc{hamdad2026paper6,
  author  = {Farid Hamdad},
  title   = {Équation pure BuP : de l'intrication à l'équation d'Einstein effective},
  year    = {2026},
  note    = {BuP Paper 6},
  url     = {https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity}
}
Liens
GitHub : Farid-Hamdad/Bottom-up-Quantum-Gravity

Paper 5 : paper5_effective_einstein/

Paper 6 (source) : paper6_pure_bup_equation/main_fr.tex

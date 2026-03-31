# BuP Vacuum — Sélection de Symétrie Effective

Cette expérience vise à tester si l’état fondamental d’un Hamiltonien inspiré de BuP
présente une structure de symétrie émergente.

On considère un Hamiltonien de type XXZ avec des champs locaux aléatoires faibles :

H = -Jzz Σ ZiZj - Jxy Σ (XiXj + YiYj) + Σ hi · σi

L’objectif est de déterminer si le vide reste isotrope (type SO(3))
ou s’il sélectionne dynamiquement un axe privilégié (type U(1)).

---

## Méthode

1. Calcul de l’état fondamental |ψ₀⟩ du Hamiltonien H
2. Extraction des matrices de densité réduites à 2 qubits
3. Construction des tenseurs de corrélation :
   T_ij = Tr(ρ_ij σ_a ⊗ σ_b)
4. Projection de T_ij vers SO(3) via décomposition polaire
5. Conversion des rotations en vecteurs de l’algèbre de Lie
6. Analyse de l’alignement de ces vecteurs :

* alignement fort → structure U(1)-like
* distribution isotrope → structure SO(3)-like

---

## Observables

* Paramètre d’ordre (OP) : alignement moyen avec l’axe dominant
* Ratio singulier : λ₁ / λ₂ (SVD des vecteurs de Lie)
* Dispersion angulaire σ
* Fraction d’arêtes alignées

Critère de classification :

* U(1)-like si ratio > 2.5 et σ < 0.4
* sinon SO(3)-like

---

## Résultats

### Robustesse (h = 0.1)

* U(1)-like dans 95 % des seeds
* OP ≈ 0.97 (contre ≈ 0.62 pour des états aléatoires)
* alignement fort et stable

### Scan en désordre (h_scale)

* régime U(1)-like stable pour h ≲ 0.2
* dégradation progressive à partir de h ≈ 0.3
* régime majoritairement perdu pour h ≥ 0.5

### Scan en anisotropie (Jzz / Jxy)

* régime isotrope pour Jzz/Jxy < 1
* transition autour de Jzz/Jxy ≈ 1
* régime U(1)-like robuste pour Jzz/Jxy ≥ 2

### Contrôle aléatoire

* aucune détection de structure U(1)-like
* OP ≈ 0.62
* confirme que l’effet n’est pas trivial

---

## Interprétation

Le vide effectif BuP n’est pas isotrope.

Dans un régime de faible désordre, il sélectionne dynamiquement un axe privilégié,
ce qui correspond à une symétrie résiduelle de type U(1).

Cet effet :

* est absent pour des états aléatoires
* dépend fortement de l’anisotropie XXZ
* est détruit par un désordre trop fort

---

## Implications pour BuP

* Le vide possède une structure intrinsèque (il n’est pas neutre)
* La géométrie d’intrication est anisotrope
* Une symétrie résiduelle peut émerger dynamiquement
* Le vide se comporte comme un milieu structuré

Ce résultat soutient l’idée que, dans BuP,
la géométrie et les symétries ne sont pas imposées,
mais émergent de la structure d’intrication.

---

## Fichiers

* `bup_vacuum_v3_standalone.py`
* `bup_vacuum_v3_results.csv`
* `bup_vacuum_v3_summary.json`

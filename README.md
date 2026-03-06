# Bottom-Up Quantum Gravity

**Émergence de l'espace, du temps et de la gravité à partir de l'intrication quantique**

Farid Hamdad — Février 2026

---

# Idée centrale

> **L'espace, le temps et la gravité ne sont pas fondamentaux.**  
> Ils émergent collectivement de la structure d'intrication d'un état quantique global fini.

Dans cette approche :

- l'intrication définit la connectivité
- la géométrie émerge de la structure informationnelle
- la gravité apparaît comme une thermodynamique de l'intrication

---

# Pourquoi ce projet ?

La physique moderne décrit avec une précision remarquable :

- la physique quantique (théories des champs, information quantique)
- la gravitation classique (relativité générale, thermodynamique des trous noirs)

mais laisse ouverte une question fondamentale :

**Pourquoi l'espace-temps existe-t-il ?**

et pourquoi la gravité possède-t-elle à la fois une structure :

- géométrique
- thermodynamique ?

Ce projet explore une hypothèse minimale :

> **L'espace-temps n'est pas le théâtre de la physique.  
> Il est reconstruit à partir de l'intrication quantique.**

---

# 1. Postulat minimal

Il existe un état quantique global pur

Ψ ∈ ⊗ᵢ Hᵢ

défini sur **N degrés de liberté élémentaires (qubits)**

sans :

- espace
- temps
- métrique préalable

Tout le reste doit émerger :

- espace
- temps
- dimension
- géométrie
- gravité effective

---

# 2. Méthodologie d'émergence

## 2.1 Émergence du temps — flot modulaire

Pour un sous-système A :

ρA = Tr_{Ā} |Ψ⟩⟨Ψ|

KA = − log(ρA)

Le flot modulaire :

O(τ) = e^{iKA τ} O e^{-iKA τ}

définit une **dynamique relationnelle intrinsèque**.

Lien conceptuel : **Page–Wootters mechanism**

→ le temps devient une propriété informationnelle interne.

---

## 2.2 Chaos modulaire (résultat)

On analyse le spectre du Hamiltonien modulaire KA.

Statistique des gaps :

⟨r⟩ = ⟨ min(Δn, Δn+1) / max(Δn, Δn+1) ⟩

Valeurs universelles :

| régime | valeur |
|------|------|
| Poisson | ≈ 0.386 |
| GOE | ≈ 0.536 |
| GUE | ≈ 0.603 |

Résultat obtenu :

⟨r⟩ ∈ [0.53 , 0.59]

→ **signature de chaos quantique (Random Matrix Theory)**

---

## 2.3 Spectral Form Factor

g₂(t) = (1 / d_A²) | Σ exp(-i t κ̃_n) |²

Structure observée :

dip → ramp → plateau

Scaling universel :

g₂_plateau ~ 1 / d_A

---

## 2.4 Constante modulaire topologique

On définit :

C = d_A × g₂_plateau

Résultats (d_A = 256)

| topologie | ⟨r⟩ | C |
|------|------|------|
| chaîne 1D | 0.594 | 1.21 |
| grille 3×6 | 0.575 | 1.42 |
| graphe ER | 0.528 | 1.63 |

➡ La géométrie d'intrication influence la dynamique modulaire.

---

# 3. Émergence de l'espace

La géométrie est reconstruite à partir de l'information mutuelle.

I(i:j) = S(ρi) + S(ρj) − S(ρij)

Distance informationnelle :

d_ij = − log( I(i:j) / I_max + ε )

Puis :

MDS → points xi ∈ ℝᵈ

La dimension émergente est la dimension minimale stabilisant l'erreur.

---

# 4. Résultats principaux

## Dimension émergente

| configuration | intrication | dimension |
|------|------|------|
| N=9 λ≈0 | locale | d≈2 |
| N=9 λ→1 | non-locale | d≈3 |
| N=16 λ≈0 | locale | d≈2 |
| N=16 λ→1 | non-locale | d≈3 |

---

## ER = EPR opérationnel

Des qubits éloignés topologiquement deviennent proches géométriquement lorsque l'intrication non-locale augmente.

→ signature wormhole-like discrète.

---

## Gravité thermodynamique

Test analogue à Jacobson :

δS ≈ β_eff δE

relation stable pour :

N = 9  
N = 16

---

# 5. Détection d'horizon bottom-up (N = 16)

On introduit un benchmark basé sur le graphe d'intrication.

Poids :

W_ij = I(i:j)

Seuil de densité :

ρ = 1/3

Score :

Score_BH = w_S z(S) − w_ϕ z(ϕ) + w_int z(internal)

où

- S = entropie région
- ϕ = conductance
- internal = intrication interne

---

## Région BH-like détectée

[0, 2, 3, 5, 6, 10, 11, 15]

| quantité | valeur |
|------|------|
| entropie | 7.1667 |
| cut | 3.2098 |
| internal | 5.0348 |
| conductance | 0.3951 |
| RMT ratio | 0.6041 |

---

## Signature

- intrication interne très élevée
- couplage extérieur faible
- bottleneck informationnel

➡ **un horizon peut émerger sans géométrie préalable**

---

# 6. Géométrie multipartite de l'intrication (nouveau résultat)

On étudie l'information mutuelle conditionnelle :

CMI(i:j|k)

Score triangulaire :

T(i,j|k) = I(i:j) − I(i:k) − I(j:k)

Résultats pour N=16 :

| λ | AUC | ρ |
|---|---|---|
| 0.0 | 0.69 | 0.31 |
| 0.2 | 0.73 | 0.42 |
| 0.4 | 0.71 | 0.37 |
| 0.5 | 0.72 | 0.39 |
| 0.7 | 0.73 | 0.42 |

Significativité :

p < 5×10⁻³

Conclusion :

> La CMI est prédite par une **géométrie triangulaire de l'intrication**.

La structure multipartite est donc encodée dans une géométrie d'ordre supérieur.

---

# 7. Ce que le projet établit

✓ une géométrie peut émerger d'un état quantique fini  
✓ la dimension dépend de l'intrication  
✓ ER = EPR devient mesurable  
✓ une thermodynamique de l'intrication apparaît  
✓ le flot modulaire est chaotique  
✓ le plateau SFF suit 1/d_A  
✓ la constante modulaire dépend de la topologie  
✓ un horizon informationnel peut émerger  
✓ la CMI possède une structure géométrique triangulaire

---

# 8. Limites

✗ limite continue N→∞ non démontrée  
✗ pas de dynamique relativiste complète  
✗ équations d'Einstein non dérivées  
✗ pas encore de prédictions observationnelles

---

# Organisation du dépôt
paper/
Bottom-up_Quantum_Gravity.tex
Bottom-up_Quantum_Gravity.pdf

theory/
theory_mathematics.tex

figures/

scripts/
bh_benchmark_louvain_N16.py

code/

data/


---

# Perspectives

- finite-size scaling N → 25+
- états critiques / topologiques
- connexion avec SYK
- dimension spectrale multi-échelle
- implémentation sur simulateurs quantiques

---

# Citation

@misc{hamdad2026bottomup,
author = {Hamdad, Farid},
title = {Bottom-Up Quantum Gravity: Emergence of Space, Time and Gravity from Quantum Entanglement},
year = {2026},
howpublished = {GitHub repository},
url = {https://github.com/Farid-Hamdad/Bottom-Up-Quantum-Gravity}
}

---

contact  
hamdadfarid54@gmail.com

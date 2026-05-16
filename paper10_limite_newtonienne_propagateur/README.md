# BuP Paper 10 — Limite quasi-newtonienne du propagateur d’intrication

**Titre FR :** Limite quasi-newtonienne du propagateur d’intrication : fonction de Green discrète, dimension spectrale et corrections BuP  
**Titre EN :** Quasi-Newtonian Limit of the Entanglement Propagator: Discrete Green Function, Spectral Dimension and BuP Corrections  
**Auteur :** Farid Hamdad  
**Année :** 2026

---

## Idée centrale

Paper 10 étudie le propagateur gravitationnel discret associé au graphe d’intrication BuP.

À partir du Laplacien d’intrication

\[
L_{\rm ent},
\]

on définit la réponse à une source ponctuelle par :

\[
L_{\rm ent}\Phi
=
\delta_{i_0}
-
\frac{1}{N}.
\]

La solution est :

\[
\Phi
=
L_{\rm ent}^{+}
\left(
\delta_{i_0}
-
\frac{1}{N}
\right),
\]

où \(L_{\rm ent}^{+}\) est la pseudo-inverse de Moore--Penrose du Laplacien.

L’objectif est de tester si, lorsque la dimension spectrale du graphe est proche de trois,

\[
d_s\simeq3,
\]

le propagateur possède une loi quasi-newtonienne :

\[
\Phi(r)
\simeq
C
+
A r^{-\alpha},
\qquad
\alpha\simeq1.
\]

---

## Résultat principal

Le résultat principal est obtenu sur une matrice MI géométrique :

\[
N=800,
\qquad
k=12,
\qquad
\lambda=0.3.
\]

Pour la réalisation `seed0`, la dimension spectrale mesurée par noyau de chaleur est :

\[
d_s=2.964.
\]

La relation spectrale attendue serait :

\[
\alpha_{\rm pred}=d_s-2=0.964.
\]

Dans la fenêtre intermédiaire :

\[
r/r_{\max}\in[0.15,0.50],
\]

le fit inter-centres du potentiel ponctuel donne, pour le filtre strict \(R^2>0.85\),

\[
\alpha_{\rm med}=1.1046,
\]

avec \(109/150\) centres valides.

Deux vérifications sur la même matrice donnent :

\[
\alpha_{\rm med}=1.1547
\]

avec \(300\) centres, et

\[
\alpha_{\rm med}=1.2022
\]

avec un autre seed de sélection des centres.

Une seconde réalisation indépendante de la matrice MI (`matrix seed1`) donne :

\[
d_s=2.970,
\qquad
\alpha_{\rm pred}=0.970,
\]

et, dans la même fenêtre :

\[
\alpha_{\rm med}=1.2597
\quad
(R^2>0.85).
\]

Le résultat robuste est donc :

\[
\boxed{
d_s\simeq3
\quad\Longrightarrow\quad
\alpha\simeq1.1-1.3.
}
\]

Cela établit une fenêtre quasi-newtonienne du propagateur d’intrication.

---

## Interprétation physique

Dans la gravité newtonienne continue, le potentiel gravitationnel vérifie :

\[
\Phi(r)\sim \frac{1}{r}.
\]

Dans BuP, cette loi n’est pas imposée. Elle doit émerger du spectre du graphe d’intrication.

Paper 10 montre que :

\[
L_{\rm ent}^{+}
\]

se comporte comme un propagateur de Green discret et qu’une loi proche de \(1/r\) apparaît lorsque le graphe atteint un régime spectral quasi tridimensionnel.

La loi newtonienne exacte

\[
\alpha=1
\]

doit être comprise comme une limite idéale continue. Sur un graphe fini, l’exposant effectif peut différer de \(1\), à cause :

- de la taille finie ;
- de la discrétisation ;
- de l’hétérogénéité du graphe ;
- du choix de la fenêtre de fit ;
- de la dimension spectrale effective \(d_s\).

---

## Résultat conceptuel

Paper 10 établit la chaîne :

\[
\text{matrice MI}
\longrightarrow
L_{\rm ent}
\longrightarrow
L_{\rm ent}^{+}
\longrightarrow
\Phi(r)
\longrightarrow
\alpha.
\]

Le résultat non trivial est que la limite newtonienne n’est pas postulée : elle apparaît comme une fenêtre effective du propagateur d’intrication.

\[
\boxed{
\text{Newton émerge comme approximation lorsque } d_s\simeq3.
}
\]

---

## Résultats numériques principaux

### 1. Calibration sur grilles contrôlées

Sur des grilles géométriques contrôlées, la réponse ponctuelle du Laplacien donne une loi de puissance en 3D.

Exemple :

\[
3D,\ L=11:
\qquad
\alpha=1.152,
\qquad
R^2=0.996.
\]

Ce test valide le protocole point-source.

### 2. Matrice MI quasi-3D — seed0

Pour la matrice MI \(N=800,\lambda=0.3,k=12\) :

\[
d_s=2.964,
\qquad
d_s-2=0.964.
\]

Fenêtre retenue :

\[
r/r_{\max}\in[0.15,0.50].
\]

Résultat principal :

\[
\alpha_{\rm med}=1.1046
\quad
(R^2>0.85).
\]

Vérifications :

| Test | Centres | Seed centres | \(R^2\) | \(\alpha_{\rm med}\) |
|---|---:|---:|---:|---:|
| principal | 150 | 0 | \(>0.85\) | 1.1046 |
| plus de centres | 300 | 0 | \(>0.85\) | 1.1547 |
| autre seed centres | 150 | 1 | \(>0.85\) | 1.2022 |

### 3. Matrice MI quasi-3D — matrix seed1

Une deuxième réalisation de la matrice MI donne :

\[
d_s=2.970,
\qquad
d_s-2=0.970.
\]

Dans la même fenêtre :

\[
r/r_{\max}\in[0.15,0.50],
\]

on obtient :

\[
\alpha_{\rm med}=1.2597
\quad
(R^2>0.85).
\]

Cela confirme que le régime quasi-newtonien ne dépend pas d’une unique réalisation de la matrice MI.

---

## Conclusion courte

Paper 10 ne démontre pas encore une limite newtonienne exacte. Il établit plus précisément :

\[
\boxed{
\text{une fenêtre quasi-newtonienne robuste du propagateur d’intrication.}
}
\]

La suite naturelle est Paper 11 :

- démontrer analytiquement la relation asymptotique \(\alpha\simeq d_s-2\) ;
- étudier la convergence grand \(N\) ;
- relier la correction \(\alpha-1\) à la dimension critique cosmologique BuP ;
- tester les conséquences sur les courbes de rotation et la cosmologie.

---

## Structure

```text
paper10_limite_newtonienne_propagateur/
├── README.md
├── paper/
│   └── main_fr.tex
├── scripts/
│   ├── README.md
│   ├── bup_paper10_point_source_green_v3.py
│   ├── bup_paper10_sparse_newton_compare_v4.py
│   ├── bup_paper10_generate_mi_ds_target_scan_v1.py
│   └── bup_paper10_mi_ds3_convergence.py
├── results/
│   ├── README.md
│   ├── mi_ds_check_N800_lam030_k12/
│   ├── mi_ds_check_N800_lam030_k12_seed1/
│   ├── mi_N800_alpha_window_015_050/
│   ├── mi_N800_alpha_window_015_050_n300_seed0/
│   ├── mi_N800_alpha_window_015_050_n150_seed1/
│   ├── mi_N800_alpha_window_015_050_matrix_seed1/
│   ├── mi_N800_alpha_window_010_040/
│   └── mi_N800_alpha_window_020_060/
└── figures/
    └── README.md

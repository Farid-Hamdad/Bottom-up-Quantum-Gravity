# Scripts — BuP Paper 10

Ce dossier contient les scripts utilisés pour le **BuP Paper 10 :**

*Limite quasi-newtonienne du propagateur d’intrication.*

L’objectif est de tester si la pseudo-inverse du Laplacien d’intrication,

$$
L_{\rm ent}^{+},
$$

se comporte comme un propagateur gravitationnel discret.

---

## 1. `bup_paper10_point_source_green_v3.py`

### Objectif

Tester le propagateur point-source sur des grilles contrôlées 1D, 2D et 3D.

Le problème résolu est :

$$
L \Phi = \delta_{i_0} - \frac{1}{N}.
$$

Le potentiel radial est ensuite ajusté par :

$$
\Phi(r) = C + A r^{-\alpha}.
$$

### Usage typique

```bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_point_source_green_v3.py \
  --L1 120 \
  --L2 31 \
  --L3 13 \
  --bins 18 \
  --tau-min 0.01 \
  --tau-max 100 \
  --tau-points 60 \
  --periodic \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/point_source_green_v3_L3_13_periodic_torus
Résultat attendu
En 3D :

α
≃
1.
α≃1.
Ce script valide le protocole point-source sur des graphes géométriques contrôlés.

2. bup_paper10_sparse_newton_compare_v4.py
Objectif
Comparer deux modèles sur grille 3D sparse.

Modèle newtonien fixé

Φ
(
r
)
=
C
+
A
r
.
Φ(r)=C+ 
r
A
​
 .
Modèle libre

Φ
(
r
)
=
C
+
A
r
−
α
.
Φ(r)=C+Ar 
−α
 .
La comparaison est faite avec :

R
2
R 
2
 

AIC

BIC

Usage typique
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_sparse_newton_compare_v4.py \
  --L-list 13 17 21 \
  --periodic \
  --bins 50 \
  --rmin-frac 0.10 \
  --rmax-frac 0.40 \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/sparse_newton_compare_v4_L13_17_21
Lecture
Si le modèle newtonien est comparable au modèle libre, la loi 
1
/
r
1/r est suffisante.

Si le modèle libre est préféré, le graphe présente une correction effective :

Φ
(
r
)
∼
r
−
α
,
α
≠
1.
Φ(r)∼r 
−α
 ,α

=1.
3. bup_paper10_generate_mi_ds_target_scan_v1.py
Objectif
Générer des matrices MI géométriques contrôlées et chercher des régimes où :

d
s
≃
3.
d 
s
​
 ≃3.
Les poids MI sont générés par :

W
i
j
=
exp
⁡
[
−
(
r
i
j
λ
)
p
]
.
W 
ij
​
 =exp[−( 
λ
r 
ij
​
 
​
 ) 
p
 ].
Le script construit ensuite un graphe 
k
k-NN et mesure la dimension spectrale.

Génération du candidat seed0
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_generate_mi_ds_target_scan_v1.py \
  --N-list 800 \
  --dim-list 3 \
  --lambda-list 0.30 \
  --k-list 12 \
  --seed-list 0 \
  --periodic \
  --target-ds 3.0 \
  --eig-k 600 \
  --tau-points 100 \
  --save-all-mi \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12
Génération du candidat matrix seed1
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_generate_mi_ds_target_scan_v1.py \
  --N-list 800 \
  --dim-list 3 \
  --lambda-list 0.30 \
  --k-list 12 \
  --seed-list 1 \
  --periodic \
  --target-ds 3.0 \
  --eig-k 600 \
  --tau-points 100 \
  --save-all-mi \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12_seed1
Résultats principaux
Seed0 : 
d
s
=
2.964
d 
s
​
 =2.964

Seed1 : 
d
s
=
2.970
d 
s
​
 =2.970

4. bup_paper10_mi_ds3_convergence.py
Objectif
Tester le propagateur point-source sur une matrice MI quasi-3D.

Le script :

charge une matrice MI ;

construit le graphe 
k
k-NN ;

calcule 
d
s
d 
s
​
  par noyau de chaleur ;

résout :

L
e
n
t
Φ
=
δ
i
0
−
1
N
L 
ent
​
 Φ=δ 
i 
0
​
 
​
 − 
N
1
​
 
pour plusieurs centres ;

ajuste :

Φ
(
r
)
=
C
+
A
r
−
α
Φ(r)=C+Ar 
−α
 
sur données brutes ;

filtre les centres par seuil 
R
2
R 
2
 .

Commande principale Paper 10 — seed0
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_mi_ds3_convergence.py \
  --mi-file papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12/mi_matrices/MI_N800_d3_lam0p3_k12_seed0.csv \
  --k 12 \
  --n-centers 150 \
  --seed-centers 0 \
  --r2-thresholds 0.70 0.85 \
  --r2-fit-min 0.15 \
  --r2-fit-max 0.50 \
  --tau-min 0.01 \
  --tau-max 100 \
  --tau-points 80 \
  --cg-rtol 1e-8 \
  --cg-maxiter 30000 \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_N800_alpha_window_015_050
Robustesse — 300 centres
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_mi_ds3_convergence.py \
  --mi-file papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12/mi_matrices/MI_N800_d3_lam0p3_k12_seed0.csv \
  --k 12 \
  --n-centers 300 \
  --seed-centers 0 \
  --r2-thresholds 0.70 0.85 \
  --r2-fit-min 0.15 \
  --r2-fit-max 0.50 \
  --tau-min 0.01 \
  --tau-max 100 \
  --tau-points 80 \
  --cg-rtol 1e-8 \
  --cg-maxiter 30000 \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_N800_alpha_window_015_050_n300_seed0
Robustesse — autre seed de centres
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_mi_ds3_convergence.py \
  --mi-file papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12/mi_matrices/MI_N800_d3_lam0p3_k12_seed0.csv \
  --k 12 \
  --n-centers 150 \
  --seed-centers 1 \
  --r2-thresholds 0.70 0.85 \
  --r2-fit-min 0.15 \
  --r2-fit-max 0.50 \
  --tau-min 0.01 \
  --tau-max 100 \
  --tau-points 80 \
  --cg-rtol 1e-8 \
  --cg-maxiter 30000 \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_N800_alpha_window_015_050_n150_seed1
Robustesse — matrix seed1
bash
python3 papers/paper10_limite_newtonienne_propagateur/scripts/bup_paper10_mi_ds3_convergence.py \
  --mi-file papers/paper10_limite_newtonienne_propagateur/results/mi_ds_check_N800_lam030_k12_seed1/mi_matrices/MI_N800_d3_lam0p3_k12_seed1.csv \
  --k 12 \
  --n-centers 150 \
  --seed-centers 0 \
  --r2-thresholds 0.70 0.85 \
  --r2-fit-min 0.15 \
  --r2-fit-max 0.50 \
  --tau-min 0.01 \
  --tau-max 100 \
  --tau-points 80 \
  --cg-rtol 1e-8 \
  --cg-maxiter 30000 \
  --output-dir papers/paper10_limite_newtonienne_propagateur/results/mi_N800_alpha_window_015_050_matrix_seed1

Dépendances Python
Les scripts utilisent :

numpy

pandas

scipy

matplotlib

Installation minimale :

bash
pip install numpy pandas scipy matplotlib

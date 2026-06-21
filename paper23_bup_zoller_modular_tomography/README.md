markdown
# Paper 23 — Tomographie modulaire BuP/Zoller

Ce dossier contient l'article scientifique connectant le principe tomographique BuP aux données d'intrication sur ions piégés de Zoller/Joshi 2023.

---

## Résultat central

La correspondance établie est :

$$ \boxed{ \rho_{\rm ent}^{\rm cut}(j) = \sum_{k \notin A} I(j:k) \propto T_{\rm mod}(j) = \frac{1}{\beta_j} } $$

où :
- $ \rho_{\rm ent}^{\rm cut}(j) $ est la **densité d'intrication de coupure** du graphe BuP,
- $ I(j:k) = S_j + S_k - S_{jk} $ est l'**information mutuelle** entre deux sites,
- $ \beta_j $ est le **coefficient modulaire** extrait du hamiltonien modulaire $ K_A = -\log \rho_A $,
- $ T_{\rm mod}(j) = 1/\beta_j $ est la **température modulaire locale**.

---

## Contenu du dossier

```text
main.tex
scripts/
  analyze_zoller_bup.py
  raw_pairwise_mi_test.py
  shadow_mi_test.py
  analyze_fig1_fig3_bup.py
results/
  modular_profiles/
  raw_mi/
  fig1_fig3/
figures/
Scripts
Script	Fonction
analyze_zoller_bup.py	Extrait les profils modulaires et teste la relation 
ρ
e
n
t
∝
1
/
β
j
ρ 
ent
​
 ∝1/β 
j
​
 
raw_pairwise_mi_test.py	Reconstruit l'information mutuelle paire-à-paire à partir des données Pauli brutes
shadow_mi_test.py	Compare les densités d'information mutuelle globale, interne et de coupure
analyze_fig1_fig3_bup.py	Analyse le scaling d'entropie, la décroissance de MI et les diagnostics de sous-systèmes disjoints
Résultats numériques principaux
1. Profils modulaires
État	Corrélation avec parabole CFT	
T
e
d
g
e
/
T
c
e
n
t
e
r
T 
edge
​
 /T 
center
​
 
Δ
=
1
Δ=1 GS	0.9986	4.42
Δ
=
1
Δ=1 ES	0.9948	1.66
Δ
=
1.7
Δ=1.7 GS	0.9903	7.17
Δ
=
1.7
Δ=1.7 ES	0.9795	1.76
Lecture BuP :

États fondamentaux (GS) : profil parabolique → température modulaire concentrée aux coupures → loi d'aire

États excités (ES) : profil plat → température modulaire uniforme → loi de volume

2. Scaling d'entropie
État	Scaling de 
S
A
S 
A
​
 	Lecture BuP
GS 
Δ
=
1
Δ=1	quasi constante, 
Δ
S
≈
0.115
ΔS≈0.115	loi d'aire
ES 
Δ
=
1
Δ=1	linéaire, 
R
2
=
0.995
R 
2
 =0.995, pente 0.412	loi de volume
GS 
Δ
=
1.7
Δ=1.7	quasi constante	phase localisée / gappée
ES 
Δ
=
1.7
Δ=1.7	linéaire, 
R
2
=
0.995
R 
2
 =0.995, pente 0.226	loi de volume plus faible
3. Comparaison des densités d'intrication
Densité testée	GS 
Δ
=
1
Δ=1 Pearson	GS 
Δ
=
1
Δ=1 Spearman
∑
k
I
(
j
:
k
)
∑ 
k
​
 I(j:k) globale	0.050	-0.088
∑
k
∈
A
I
(
j
:
k
)
∑ 
k∈A
​
 I(j:k) interne	-0.212	-0.286
∑
k
∉
A
I
(
j
:
k
)
∑ 
k∈
/
A
​
 I(j:k) coupure	0.588	0.777
Conclusion : La densité de coupure 
ρ
e
n
t
c
u
t
(
j
)
=
∑
k
∉
A
I
(
j
:
k
)
ρ 
ent
cut
​
 (j)=∑ 
k∈
/
A
​
 I(j:k) est significativement mieux corrélée à 
1
/
β
j
1/β 
j
​
  que la densité globale naïve.

4. Décroissance de l'information mutuelle avec la distance
M
I
(
d
=
0
)
≈
1.133
MI(d=0)≈1.133

M
I
(
d
=
max
⁡
)
≈
0.174
MI(d=max)≈0.174

Fit en loi de puissance : 
R
2
≈
0.993
R 
2
 ≈0.993

5. Liens inter-sous-systèmes et fidélité
Configuration	Fidélité
Avec liens	0.924
Sans liens	0.855
Gain moyen	0.069
Gain max	0.200
Conclusion : Les liens non-locaux du graphe d'intrication sont mesurables et améliorent la reconstruction.

Ce que ça implique pour BuP
Implication	Explication
BuP devient testable sur table	Pas seulement via galaxies ou ondes GW — la tomographie modulaire est mesurable en laboratoire.
La quasi-localité du hamiltonien modulaire est validée	Zoller montre que 
K
A
=
−
log
⁡
ρ
A
K 
A
​
 =−logρ 
A
​
  peut être appris et possède une structure locale.
La prescription BuP doit être raffinée	La bonne densité n'est pas 
∑
k
I
(
j
:
k
)
∑ 
k
​
 I(j:k) mais 
∑
k
∉
A
I
(
j
:
k
)
∑ 
k∈
/
A
​
 I(j:k).
Le proxy modulaire BuP est validé	
ρ
e
n
t
c
u
t
(
j
)
∝
1
/
β
j
ρ 
ent
cut
​
 (j)∝1/β 
j
​
 
Les phases BuP ont une signature expérimentale	GS → loi d'aire ; ES → loi de volume.
Les liens non-locaux sont mesurables	Les sous-systèmes disjoints montrent que les liens inter-régions améliorent la reconstruction.
Phrase finale
Les données de Zoller/Joshi 2023 fournissent une première validation expérimentale de la tomographie modulaire BuP : la température modulaire locale 
T
m
o
d
=
1
/
β
j
T 
mod
​
 =1/β 
j
​
  reconstruit une densité d'intrication de coupure, distingue les phases loi d'aire / loi de volume, et révèle des liens non-locaux compatibles avec un graphe d'intrication BuP mesurable.

Note sur les données
Les fichiers .mat bruts de Joshi et al. ne sont pas inclus par défaut dans ce dépôt. Les résumés dérivés (CSV/JSON) sont stockés dans results/, et la source des données originales est documentée dans l'article.

Prochaines étapes
Action	Description
Paper 23	Rédiger l'article court (ou section) sur cette validation expérimentale
Test systématique	Appliquer le même protocole à d'autres datasets de simulateurs quantiques
Prédiction BuP	La densité de coupure doit être proportionnelle à 
1
/
β
j
1/β 
j
​
  dans tout système critique
Lien avec Paper 13	La même tomographie s'applique aux données SPARC — ce n'est pas une analogie, c'est le même mécanisme

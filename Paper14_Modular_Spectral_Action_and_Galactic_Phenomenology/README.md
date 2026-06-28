markdown
# Paper 14 — Action modulaire/spectrale et phénoménologie galactique

**Bottom-Up Quantum Gravity — Pont modulaire-gravitationnel et prédiction des courbes de rotation SPARC**

Ce dossier contient les tests numériques et phénoménologiques pour le Paper 14 du programme de gravité quantique **Bottom-Up Quantum Gravity (BuP)**.

L'objectif du papier est de relier trois niveaux du cadre BuP :

1. l'hamiltonien modulaire 
```math
K_A=-\log\rho_A,
le laplacien d'intrication

math
L_{\rm ent},
la réponse gravitationnelle effective testée sur les courbes de rotation SPARC.

Le résultat central est que la même structure spectrale contrôle à la fois l'exposant modulaire

math
\beta_{\rm mod}
et l'exposant gravitationnel

math
\beta_{\rm grav},
et que la longueur de corrélation résultante

math
\lambda_{\rm corr}
peut être prédite pour les galaxies sans ajustement individuel.

Statut du dossier
Ce dossier contient les tests numériques et phénoménologiques pour le Paper 14.

La chaîne dynamique finale établie dans ce papier est :

math
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s,d_w,\lambda_2,\Delta_3)
\rightarrow
\beta_{\rm grav},\beta_{\rm mod}.
Ce résultat montre que BuP prédit les exposants modulaires et gravitationnels à partir des mêmes invariants spectraux du graphe d'intrication.

1. Motivation scientifique
La question de départ du Paper 14 est :

« La dynamique modulaire quantique et la gravité émergente peuvent-elles être dérivées du même spectre d'intrication ? »

Dans les articles précédents, BuP avait déjà introduit le laplacien d'intrication

math
L_{\rm ent}
comme le laplacien du graphe d'intrication, et l'exposant gravitationnel effectif :

math
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
Le Paper 14 demande si l'hamiltonien modulaire,

math
K_A=-\log\rho_A,
peut également être représenté comme une fonction spectrale du même laplacien d'intrication.

2. Hypothèse centrale
L'hypothèse de travail est :

math
K_A \sim f(L_{\rm ent}).
Plus précisément, les tests comparent le spectre de

math
K_A
avec les fonctions spectrales de trois laplaciens candidats :

math
L_A^{\rm induit},
\qquad
L_A^{\rm Schur},
\qquad
L_A^{\rm normalisé-MI}.
L'ajustement effectif prend la forme d'une loi de puissance :

math
K_A \sim (L_A)^{\beta_{\rm mod}}.
L'exposant

math
\beta_{\rm mod}
est ensuite comparé à l'exposant gravitationnel :

math
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
3. Test spectral modulaire
La première étape consistait à tester si

math
K_A
peut être représenté par une fonction spectrale du laplacien d'intrication.

Les meilleurs tests, réalisés avec des systèmes de taille

math
N=16,
ont donné un coefficient de détermination :

math
R^2 \simeq 0.99.
Cela montre que l'hamiltonien modulaire admet une représentation spectrale effective en termes du graphe d'intrication.

La famille spectrale positive fonctionne bien, tandis que la famille négative s'effondre avec un pouvoir explicatif proche de zéro.

4. Échec de la relation naïve
math
C\to\beta
Une première hypothèse naturelle était que la constante topologique modulaire

math
C_{\rm modular}
=
d_A g_2^{\rm plateau}
pourrait prédire directement

math
\beta_{\rm mod}.
Cette hypothèse a échoué.

Sur l'ensemble des jeux de données propres pour

math
N=9,
math
N=16
et le

math
N=16
optimal, la corrélation linéaire directe entre

math
C_{\rm modular}
et

math
\beta_{\rm mod}
reste très faible :

math
R^2 \approx 0.
Cet échec est important. Il montre que

math
\beta_{\rm mod}
n'est pas contrôlé par une seule constante globale de plateau. La structure pertinente est multivariée et spectrale.

5. Caractéristiques spectrales et diagnostics du facteur de forme spectral
Le Paper 14 a ensuite extrait des caractéristiques supplémentaires du facteur de forme spectral :

math
t_{\rm dip},
\qquad
t_{\rm ramp},
\qquad
{\rm pente}_{\rm ramp},
\qquad
\Delta_3.
La caractéristique

math
\Delta_3
est plus informative que

math
C_{\rm modular},
mais reste insuffisante par elle-même.

Ceci a conduit à un modèle multivarié utilisant neuf variables :

math
C_{\rm modular},
\quad
K_{\rm gap},
\quad
\Delta_3,
\quad
\lambda_2,
\quad
d_s,
\quad
d_w,
\quad
d_s/d_w,
\quad
\langle I\rangle,
\quad
I_{\max}.
6. Prédiction de
math
\beta_{\rm mod}
L'amélioration décisive vient de l'ajout des invariants de graphe et de diffusion :

math
d_s,
\qquad
d_w,
\qquad
\lambda_2,
\qquad
d_s/d_w.
Sur le sous-ensemble Schur/RMT (théorie des matrices aléatoires), les meilleurs modèles atteignent :

math
R^2 \simeq 0.96
en validation croisée k-fold, et restent solides dans les tests leave-one-regime-out.

Les caractéristiques les plus importantes sont :

math
\lambda_2^{\rm norm},
\qquad
d_s/d_w,
\qquad
\Delta_3,
\qquad
d_s.
Ceci établit que :

math
\beta_{\rm mod}
=
F(\lambda_2,d_s,d_w,\Delta_3,\ldots).
7. Pont modulaire-gravitationnel
En utilisant la relation gravitationnelle de BuP :

math
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4,
nous définissons :

math
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
La relation mesurée entre les deux exposants est :

math
\beta_{\rm mod}
\simeq
0.531
+
1.726\,\beta_{\rm grav}.
C'est le pont central du Paper 14.

Cela signifie que la dynamique modulaire et la gravité effective sont deux projections spectrales du même laplacien d'intrication.

8. Test du pont sur les données SPARC
Le pont a ensuite été testé sur des graphes d'intrication à l'échelle galactique reconstruits à partir des profils baryoniques SPARC.

La chaîne de traitement est :

math
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s,d_w,\lambda_2,\Delta_3)
\rightarrow
\beta_{\rm grav},\beta_{\rm mod}.
Pour chaque galaxie, le pont compare la valeur prédite

math
\beta_{\rm mod}^{\rm pont}
(issue de la relation linéaire) avec la valeur issue de l'apprentissage automatique

math
\beta_{\rm mod}^{\rm ML}.
L'erreur du pont est définie comme :

math
\epsilon_{\rm pont}
=
\frac{
|\beta_{\rm mod}^{\rm ML}-\beta_{\rm mod}^{\rm pont}|
}{
|\beta_{\rm mod}^{\rm ML}|
}.
9. Test prédictif de la longueur de corrélation
math
\lambda_{\rm corr}
Les premiers tests SPARC ont balayé plusieurs valeurs du rapport :

math
f_\lambda
=
\lambda_{\rm corr}/R_d
\in
\{0.5,1.0,1.5,2.0,3.0\},
où

math
R_d
est le rayon caractéristique du disque.

Le test décisif a ensuite été effectué sans ajuster

math
\lambda_{\rm corr}
galaxie par galaxie.

Un modèle de type leave-one-galaxy-out a été entraîné sur 174 galaxies et utilisé pour prédire :

math
\widehat f_\lambda
pour la galaxie exclue.

Le meilleur modèle était un classifieur Random Forest.

Résultat : sur les 175 galaxies,

math
173/175
passent le pont avec un

math
\lambda_{\rm corr}
prédit, sans ajustement individuel.

Les taux de réussite sont :

math
98.86\%
de pont fort ou modéré, et :

math
72.0\%
de pont fort uniquement.

Ainsi,

math
\lambda_{\rm corr}
n'est pas simplement une échelle ajustée. Elle est prédictible à partir des invariants du graphe.

10. Taxonomie unifiée des galaxies
Le Paper 14 construit également une taxonomie unifiée combinant trois critères :

la phase de dimension (basse ou haute),

la catégorie de qualité d'ajustement,

la classe de cohérence optimale correspondant à la valeur de

math
f_\lambda^{\rm opt}.
La répartition basse/haute est :

math
N_{\rm HAUTE}=88,
\qquad
N_{\rm BASSE}=87.
avec :

math
d_{\min}^{\rm HAUTE}=2.487107,
\qquad
d_{\min}^{\rm BASSE}=2.274448.
Les classes de cohérence du Paper 14 sont :

text
courte_0p5Rd          : 89 galaxies
standard_1Rd         : 38 galaxies
étendue_1p5_2Rd     : 39 galaxies
très_étendue_3Rd    : 9 galaxies
Résumé
Paper 14 établit un pont quantitatif entre la dynamique modulaire et la gravité effective à travers le spectre du laplacien d'intrication.

Le résultat final est :

math
\beta_{\rm mod}
\simeq
0.531
+
1.726\,\beta_{\rm grav}.
La longueur de corrélation

math
\lambda_{\rm corr}
est prédictible à partir des invariants spectraux, avec un taux de réussite de 98,86 % sur 175 galaxies SPARC.

La taxonomie unifiée montre une répartition cohérente entre phases de dimension et classes de cohérence, renforçant l'interprétation selon laquelle la structure du graphe d'intrication contrôle à la fois les propriétés modulaires et gravitationnelles à l'échelle galactique.

# Paper 14 — Action modulaire/spectrale et phénoménologie galactique

## Statut

Ce dossier contient les tests numériques et phénoménologiques pour le Paper 14 du programme de gravité quantique Bottom-Up (BuP).

Le Paper 14 relie trois niveaux du cadre BuP :

1. l'hamiltonien modulaire \(K_A=-\log\rho_A\),
2. le laplacien d'intrication \(L_{\rm ent}\),
3. la réponse gravitationnelle effective testée sur les courbes de rotation SPARC.

Le résultat central est que la même structure spectrale contrôle à la fois l'exposant modulaire \(\beta_{\rm mod}\) et l'exposant gravitationnel \(\beta_{\rm grav}\), et que la longueur de corrélation résultante \(\lambda_{\rm corr}\) peut être prédite pour les galaxies sans ajustement individuel.

---

## 1. Motivation scientifique

La question de départ du Paper 14 est :

\[
\text{La dynamique modulaire quantique et la gravité émergente peuvent-elles être dérivées du même spectre d'intrication ?}
\]

Dans les articles précédents, BuP avait déjà introduit :

\[
L_{\rm ent}
\]

comme le laplacien du graphe d'intrication, et l'exposant gravitationnel effectif :

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
\]

Le Paper 14 demande si l'hamiltonien modulaire,

\[
K_A=-\log\rho_A,
\]

peut également être représenté comme une fonction spectrale du même laplacien d'intrication.

---

## 2. Hypothèse centrale

L'hypothèse de travail est :

\[
K_A \sim f(L_{\rm ent}).
\]

Plus précisément, les tests comparent le spectre de \(K_A\) avec les fonctions spectrales de laplaciens candidats :

\[
L_A^{\rm induit},
\qquad
L_A^{\rm Schur},
\qquad
L_A^{\rm normalisé-MI}.
\]

L'ajustement effectif prend la forme :

\[
K_A \sim (L_A)^{\beta_{\rm mod}}.
\]

L'exposant \(\beta_{\rm mod}\) est ensuite comparé à l'exposant gravitationnel :

\[
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
\]

---

## 3. Test spectral modulaire

La première étape consistait à tester si \(K_A\) peut être représenté par une fonction spectrale du laplacien d'intrication.

Les meilleurs tests pour \(N=16\) ont donné :

\[
R^2 \simeq 0.99.
\]

Cela montre que l'hamiltonien modulaire admet une représentation spectrale effective en termes du graphe d'intrication.

La famille spectrale positive fonctionne bien, tandis que la famille négative s'effondre avec un pouvoir explicatif proche de zéro.

---

## 4. Échec de la relation naïve \(C\to\beta\)

Une première hypothèse naturelle était que la constante topologique modulaire

\[
C_{\rm modular}
=
d_A g_2^{\rm plateau}
\]

pourrait prédire directement \(\beta_{\rm mod}\).

Cette hypothèse a échoué.

Sur l'ensemble des jeux de données propres \(N=9\), \(N=16\) et \(N=16\) optimal, la corrélation linéaire directe entre \(C_{\rm modular}\) et \(\beta_{\rm mod}\) reste faible :

\[
R^2 \approx 0.
\]

Cet échec est important. Il montre que \(\beta_{\rm mod}\) n'est pas contrôlé par une seule constante globale de plateau. La structure pertinente est multivariée et spectrale.

---

## 5. Caractéristiques spectrales et diagnostics SFF

Le Paper 14 a ensuite extrait des caractéristiques supplémentaires du facteur de forme spectral :

\[
t_{\rm dip},
\qquad
t_{\rm ramp},
\qquad
{\rm pente}_{\rm ramp},
\qquad
\Delta_3.
\]

La caractéristique \(\Delta_3\) est plus informative que \(C_{\rm modular}\), mais reste insuffisante par elle-même.

Ceci a conduit à un modèle multivarié utilisant :

\[
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
\]

---

## 6. Prédiction de \(\beta_{\rm mod}\)

L'amélioration décisive vient de l'ajout des invariants de graphe et de diffusion :

\[
d_s,
\qquad
d_w,
\qquad
\lambda_2,
\qquad
d_s/d_w.
\]

Sur le sous-ensemble Schur/RMT, les meilleurs modèles atteignent :

\[
R^2 \simeq 0.96
\]

en validation croisée k-fold, et restent solides dans les tests leave-one-regime-out.

Les caractéristiques les plus importantes sont :

\[
\lambda_2^{\rm norm},
\qquad
d_s/d_w,
\qquad
\Delta_3,
\qquad
d_s.
\]

Ceci établit que :

\[
\beta_{\rm mod}
=
F(\lambda_2,d_s,d_w,\Delta_3,\ldots).
\]

---

## 7. Pont modulaire-gravitationnel

En utilisant la relation gravitationnelle de BuP :

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4,
\]

nous définissons :

\[
\beta_{\rm grav}
=
\frac{\alpha_{\rm eff}+1}{2}.
\]

La relation mesurée est :

\[
\beta_{\rm mod}
\simeq
0.531
+
1.726\,\beta_{\rm grav}.
\]

C'est le pont central du Paper 14.

Cela signifie que la dynamique modulaire et la gravité effective sont deux projections spectrales du même laplacien d'intrication.

---

## 8. Test du pont sur SPARC

Le pont a ensuite été testé sur des graphes d'intrication à l'échelle galactique reconstruits à partir des profils baryoniques SPARC.

La chaîne de traitement est :

\[
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s,d_w,\lambda_2,\Delta_3)
\rightarrow
\beta_{\rm grav},\beta_{\rm mod}.
\]

Pour chaque galaxie, le pont compare :

\[
\beta_{\rm mod}^{\rm pont}
\]

avec :

\[
\beta_{\rm mod}^{\rm ML}.
\]

L'erreur du pont est :

\[
\epsilon_{\rm pont}
=
\frac{
|\beta_{\rm mod}^{\rm ML}-\beta_{\rm mod}^{\rm pont}|
}{
|\beta_{\rm mod}^{\rm ML}|
}.
\]

---

## 9. Test prédictif de \(\lambda_{\rm corr}\)

Les premiers tests SPARC ont balayé :

\[
f_\lambda
=
\lambda_{\rm corr}/R_d
\in
\{0.5,1.0,1.5,2.0,3.0\}.
\]

Le test décisif a ensuite été effectué sans ajuster \(\lambda_{\rm corr}\) galaxie par galaxie.

Un modèle leave-one-galaxy-out a été entraîné sur 174 galaxies et utilisé pour prédire :

\[
\widehat f_\lambda
\]

pour la galaxie exclue.

Le meilleur modèle était un classifieur Random Forest.

Résultat :

\[
173/175
\]

galaxies passent le pont avec un \(\lambda_{\rm corr}\) prédit, sans ajustement individuel.

Les taux de réussite sont :

\[
98.86\%
\]

de pont fort ou modéré, et :

\[
72.0\%
\]

de pont fort.

Ainsi :

\[
\lambda_{\rm corr}
\]

n'est pas simplement une échelle ajustée. Elle est prédictible à partir des invariants du graphe.

---

## 10. Taxonomie unifiée des galaxies

Le Paper 14 construit également une taxonomie unifiée combinant :

1. la phase de dimension BASSE/HAUTE,
2. la catégorie de qualité d'ajustement,
3. la classe de cohérence optimale \(f_\lambda^{\rm opt}\).

La répartition BASSE/HAUTE est :

\[
N_{\rm HAUTE}=88,
\qquad
N_{\rm BASSE}=87.
\]

avec :

\[
d_{\min}^{\rm HAUTE}=2.487107,
\qquad
d_{\min}^{\rm BASSE}=2.274448.
\]

Les classes de cohérence du Paper 14 sont :

```text
courte_0p5Rd          : 89 galaxies
standard_1Rd         : 38 galaxies
étendue_1p5_2Rd     : 39 galaxies
très_étendue_3Rd    : 9 galaxies

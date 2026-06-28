# Paper 14 — Action modulaire/spectrale et phénoménologie galactique

## Statut

Ce dossier contient les tests numériques et phénoménologiques pour le Paper 14 du programme de gravité quantique Bottom-Up (BuP).

Le Paper 14 relie trois niveaux du cadre BuP :

1. l'hamiltonien modulaire *K indice A* égal à moins le logarithme de *rho indice A*,
2. le laplacien d'intrication *L indice ent*,
3. la réponse gravitationnelle effective testée sur les courbes de rotation SPARC.

Le résultat central est que la même structure spectrale contrôle à la fois l'exposant modulaire *bêta indice mod* et l'exposant gravitationnel *bêta indice grav*, et que la longueur de corrélation résultante *lambda indice corr* peut être prédite pour les galaxies sans ajustement individuel.

---

## 1. Motivation scientifique

La question de départ du Paper 14 est :

« La dynamique modulaire quantique et la gravité émergente peuvent-elles être dérivées du même spectre d'intrication ? »

Dans les articles précédents, BuP avait déjà introduit le laplacien d'intrication *L indice ent* comme le laplacien du graphe d'intrication, et l'exposant gravitationnel effectif *alpha indice eff* défini par la formule : alpha_eff = (2 * d_s / d_w) + d_w - 4, où *d_s* est la dimension spectrale et *d_w* la dimension de marche aléatoire.

Le Paper 14 demande si l'hamiltonien modulaire *K_A* (défini comme moins le logarithme de la matrice densité réduite *rho_A*) peut également être représenté comme une fonction spectrale du même laplacien d'intrication.

---

## 2. Hypothèse centrale

L'hypothèse de travail est que *K_A* est approximativement une fonction *f* du laplacien d'intrication *L_ent*.

Plus précisément, les tests comparent le spectre de *K_A* avec les fonctions spectrales de trois laplaciens candidats : le laplacien induit, le laplacien de Schur, et le laplacien normalisé par l'information mutuelle.

L'ajustement effectif prend la forme d'une loi de puissance : *K_A* est proportionnel à *L_A* élevé à la puissance *bêta indice mod*.

L'exposant *bêta mod* est ensuite comparé à l'exposant gravitationnel *bêta grav*, qui est défini par la relation : bêta_grav = (alpha_eff + 1) / 2.

---

## 3. Test spectral modulaire

La première étape consistait à tester si *K_A* peut être représenté par une fonction spectrale du laplacien d'intrication.

Les meilleurs tests, réalisés avec des systèmes de taille *N = 16*, ont donné un coefficient de détermination *R au carré* d'environ 0,99.

Cela montre que l'hamiltonien modulaire admet une représentation spectrale effective en termes du graphe d'intrication.

La famille spectrale positive fonctionne bien, tandis que la famille négative s'effondre avec un pouvoir explicatif proche de zéro.

---

## 4. Échec de la relation naïve entre la constante topologique C et bêta_mod

Une première hypothèse naturelle était que la constante topologique modulaire *C_modular*, définie comme le produit de la dimension *d_A* et de la hauteur de plateau *g_2* (notée g_2^plateau), pourrait prédire directement *bêta_mod*.

Cette hypothèse a échoué.

Sur l'ensemble des jeux de données propres (pour *N = 9*, *N = 16* et le *N = 16* optimal), la corrélation linéaire directe entre *C_modular* et *bêta_mod* reste très faible, avec un *R au carré* pratiquement nul.

Cet échec est important. Il montre que *bêta_mod* n'est pas contrôlé par une seule constante globale de plateau. La structure pertinente est multivariée et spectrale.

---

## 5. Caractéristiques spectrales et diagnostics du facteur de forme spectral (SFF)

Le Paper 14 a ensuite extrait des caractéristiques supplémentaires du facteur de forme spectral : le temps du dip (*t_dip*), le temps de la rampe (*t_ramp*), la pente de la rampe, et la statistique *Delta_3*.

La caractéristique *Delta_3* est plus informative que *C_modular*, mais reste insuffisante par elle-même.

Ceci a conduit à un modèle multivarié utilisant neuf variables : *C_modular*, l'écart spectral (*K_gap*), *Delta_3*, la seconde valeur propre du laplacien (*lambda_2*), la dimension spectrale *d_s*, la dimension de marche *d_w*, le rapport *d_s/d_w*, la moyenne de l'information mutuelle (*<I>*), et l'information mutuelle maximale (*I_max*).

---

## 6. Prédiction de *bêta_mod*

L'amélioration décisive vient de l'ajout des invariants de graphe et de diffusion que sont *d_s*, *d_w*, *lambda_2* et le rapport *d_s/d_w*.

Sur le sous-ensemble Schur/RMT (théorie des matrices aléatoires), les meilleurs modèles atteignent un *R au carré* d'environ 0,96 en validation croisée k-fold, et restent solides dans les tests où l'on retire un régime entier d'apprentissage.

Les caractéristiques les plus importantes sont : *lambda_2 normalisé*, le rapport *d_s/d_w*, *Delta_3*, et *d_s*.

Ceci établit que *bêta_mod* est une fonction de ces invariants : bêta_mod = F(lambda_2, d_s, d_w, Delta_3, ...).

---

## 7. Pont modulaire-gravitationnel

En utilisant la relation gravitationnelle de BuP : alpha_eff = (2*d_s/d_w) + d_w - 4, nous définissons bêta_grav = (alpha_eff + 1) / 2.

La relation mesurée entre les deux exposants est : bêta_mod ≈ 0,531 + 1,726 * bêta_grav.

C'est le pont central du Paper 14.

Cela signifie que la dynamique modulaire et la gravité effective sont deux projections spectrales du même laplacien d'intrication.

---

## 8. Test du pont sur les données SPARC

Le pont a ensuite été testé sur des graphes d'intrication à l'échelle galactique reconstruits à partir des profils baryoniques SPARC.

La chaîne de traitement est la suivante : à partir de la densité de surface *Sigma(R)*, on construit la matrice de poids *W_ij*, puis le laplacien d'intrication *L_ent*, on en extrait les invariants (*d_s, d_w, lambda_2, Delta_3*), et on calcule les exposants *bêta_grav* et *bêta_mod*.

Pour chaque galaxie, le pont compare la valeur prédite *bêta_mod^pont* (issue de la relation linéaire) avec la valeur issue de l'apprentissage automatique *bêta_mod^ML*.

L'erreur du pont est définie comme la valeur absolue de la différence entre *bêta_mod^ML* et *bêta_mod^pont*, divisée par la valeur absolue de *bêta_mod^ML*.

---

## 9. Test prédictif de la longueur de corrélation *lambda_corr*

Les premiers tests SPARC ont balayé plusieurs valeurs du rapport *f_lambda = lambda_corr / R_d* (où *R_d* est le rayon caractéristique du disque), à savoir 0,5, 1,0, 1,5, 2,0 et 3,0.

Le test décisif a ensuite été effectué sans ajuster *lambda_corr* galaxie par galaxie.

Un modèle de type « leave-one-galaxy-out » a été entraîné sur 174 galaxies et utilisé pour prédire la valeur de *f_lambda* pour la galaxie exclue.

Le meilleur modèle était un classifieur Random Forest.

Résultat : sur les 175 galaxies, 173 passent le pont avec un *lambda_corr* prédit, sans ajustement individuel.

Les taux de réussite sont de 98,86 % pour un pont fort ou modéré, et de 72,0 % pour un pont fort uniquement.

Ainsi, *lambda_corr* n'est pas simplement une échelle ajustée. Elle est prédictible à partir des invariants du graphe.

---

## 10. Taxonomie unifiée des galaxies

Le Paper 14 construit également une taxonomie unifiée combinant trois critères :

1. la phase de dimension (basse ou haute),
2. la catégorie de qualité d'ajustement,
3. la classe de cohérence optimale correspondant à la valeur de *f_lambda^opt*.

La répartition basse/haute est de 88 galaxies pour la phase haute et 87 pour la phase basse, avec des dimensions minimales respectives de 2,487 et 2,274.

Les classes de cohérence du Paper 14 sont :

```text
courte_0p5Rd          : 89 galaxies
standard_1Rd         : 38 galaxies
étendue_1p5_2Rd     : 39 galaxies
très_étendue_3Rd    : 9 galaxies

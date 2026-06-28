# Paper 14 — Action modulaire/spectrale et phénoménologie galactique

## Situation

Ce dossier rassemble l’ensemble des tests numériques et phénoménologiques menés pour l’article 14 du programme de gravité quantique *Bottom‑Up* (BuP).

L’article 14 établit un pont entre trois niveaux distincts du cadre BuP :

1. l’hamiltonien modulaire $K_A=-\log\rho_A$,
2. le laplacien d’intrication $L_{\rm ent}$,
3. la réponse gravitationnelle effective, testée sur les courbes de rotation SPARC.

Le résultat central est le suivant : une même structure spectrale gouverne à la fois l’exposant modulaire $\beta_{\rm mod}$ et l’exposant gravitationnel $\beta_{\rm grav}$, ce qui permet de prédire la longueur de corrélation $\lambda_{\rm corr}$ pour chaque galaxie sans avoir à l’ajuster individuellement.

---

## 1. Motivation scientifique

La question qui guide l’article 14 est la suivante :

$$
\text{La dynamique modulaire quantique et la gravité émergente peuvent-elles dériver d’un même spectre d’intrication ?}
$$

Dans les articles antérieurs, le programme BuP avait déjà introduit

$$
L_{\rm ent}
$$

comme le laplacien du graphe d’intrication, ainsi que l’exposant gravitationnel effectif

$$
\alpha_{\rm eff} = \frac{2d_s}{d_w}+d_w-4.
$$

L’article 14 se demande si l’hamiltonien modulaire,

$$
K_A=-\log\rho_A,
$$

peut, lui aussi, être représenté comme une fonction spectrale du même laplacien d’intrication.

---

## 2. Hypothèse centrale

L’hypothèse de travail est la suivante :

$$
K_A \sim f(L_{\rm ent}).
$$

Plus précisément, les tests comparent le spectre de $K_A$ aux fonctions spectrales de différents laplaciens candidats :

$$
L_A^{\rm induit}, \qquad L_A^{\rm Schur}, \qquad L_A^{\rm normalisé\ par\ l’IM}.
$$

L’ajustement effectif prend la forme :

$$
K_A \sim (L_A)^{\beta_{\rm mod}}.
$$

L’exposant $\beta_{\rm mod}$ est ensuite confronté à l’exposant gravitationnel :

$$
\beta_{\rm grav} = \frac{\alpha_{\rm eff}+1}{2}.
$$

---

## 3. Test spectral modulaire

La première étape consiste à vérifier si $K_A$ admet une représentation spectrale en termes du laplacien d’intrication.

Pour les meilleurs tests (avec $N=16$), on obtient :

$$
R^2 \simeq 0,99.
$$

Ce résultat démontre que l’hamiltonien modulaire possède effectivement une représentation spectrale sur le graphe d’intrication.

La famille spectrale positive fonctionne parfaitement, tandis que la famille négative s’effondre, avec un pouvoir explicatif quasi nul.

---

## 4. Échec de la relation naïve $C\to\beta$

Une première hypothèse, naturelle, consistait à penser que la constante topologique modulaire

$$
C_{\rm modulaire} = d_A g_2^{\rm plateau}
$$

permettait de prédire directement $\beta_{\rm mod}$.

Cette hypothèse est invalidée par les données. Sur les ensembles propres $N=9$, $N=16$ et $N=16$ optimal, la corrélation linéaire directe entre $C_{\rm modulaire}$ et $\beta_{\rm mod}$ reste très faible :

$$
R^2 \approx 0.
$$

Cet échec est important : il montre que $\beta_{\rm mod}$ n’est pas dicté par une seule constante de plateau globale. La structure pertinente est multidimensionnelle et de nature spectrale.

---

## 5. Caractéristiques spectrales et diagnostics SFF

L’article 14 extrait ensuite un ensemble de caractéristiques issues du facteur de forme spectral (SFF) :

$$
t_{\rm creux}, \qquad t_{\rm rampe}, \qquad {\rm pente}_{\rm rampe}, \qquad \Delta_3.
$$

La quantité $\Delta_3$ s’avère plus informative que $C_{\rm modulaire}$, mais elle reste insuffisante prise isolément.

On construit alors un modèle multivarié utilisant :

$$
C_{\rm modulaire}, \quad K_{\rm intervalle}, \quad \Delta_3, \quad \lambda_2, \quad d_s, \quad d_w, \quad d_s/d_w, \quad \langle I\rangle, \quad I_{\max}.
$$

---

## 6. Prédiction de $\beta_{\rm mod}$

L’amélioration décisive provient de l’ajout des invariants de graphe et de diffusion :

$$
d_s, \qquad d_w, \qquad \lambda_2, \qquad d_s/d_w.
$$

Sur le sous‑ensemble Schur/RMT, les meilleurs modèles atteignent :

$$
R^2 \simeq 0,96
$$

en validation croisée $k$-fold, et restent robustes lors des tests *leave‑one‑regime‑out*.

Les caractéristiques les plus importantes sont :

$$
\lambda_2^{\rm norm}, \qquad d_s/d_w, \qquad \Delta_3, \qquad d_s.
$$

On établit ainsi que :

$$
\beta_{\rm mod} = F(\lambda_2,d_s,d_w,\Delta_3,\ldots).
$$

---

## 7. Pont modulaire‑gravitationnel

En utilisant la relation gravitationnelle de BuP :

$$
\alpha_{\rm eff} = \frac{2d_s}{d_w}+d_w-4,
$$

on définit :

$$
\beta_{\rm grav} = \frac{\alpha_{\rm eff}+1}{2}.
$$

La relation mesurée s’écrit :

$$
\beta_{\rm mod} \simeq 0,531 + 1,726\,\beta_{\rm grav}.
$$

C’est le pont central de l’article 14.

Il signifie que la dynamique modulaire et la gravité effective ne sont que deux projections spectrales d’un même laplacien d’intrication.

---

## 8. Test du pont sur SPARC

Le pont est ensuite testé sur des graphes d’intrication à l’échelle galactique, reconstruits à partir des profils baryoniques SPARC.

Le pipeline s’articule comme suit :

$$
\Sigma(R) \rightarrow W_{ij} \rightarrow L_{\rm ent} \rightarrow (d_s,d_w,\lambda_2,\Delta_3) \rightarrow \beta_{\rm grav},\beta_{\rm mod}.
$$

Pour chaque galaxie, le pont confronte :

$$
\beta_{\rm mod}^{\rm pont}
$$

à

$$
\beta_{\rm mod}^{\rm ML}.
$$

L’erreur du pont est définie par :

$$
\epsilon_{\rm pont} = \frac{ |\beta_{\rm mod}^{\rm ML}-\beta_{\rm mod}^{\rm pont}| }{ |\beta_{\rm mod}^{\rm ML}| }.
$$

---

## 9. Test prédictif sur $\lambda_{\rm corr}$

Les premiers tests sur SPARC explorent les valeurs :

$$
f_\lambda = \lambda_{\rm corr}/R_d \in \{0,5\,;\ 1,0\,;\ 1,5\,;\ 2,0\,;\ 3,0\}.
$$

Le test décisif est ensuite réalisé sans ajuster $\lambda_{\rm corr}$ galaxie par galaxie.

Un modèle de type *Random Forest* est entraîné sur 174 galaxies et utilisé pour prédire

$$
\widehat f_\lambda
$$

pour la galaxie retenue.

Résultat :

$$
173/175
$$

galaxies satisfont au pont avec la valeur prédite de $\lambda_{\rm corr}$, sans aucun ajustement individuel.

Les taux de réussite sont de :

$$
98,86\%
$$

pour un pont fort ou modéré, et de :

$$
72,0\%
$$

pour un pont fort.

Ainsi,

$$
\lambda_{\rm corr}
$$

n’est pas une simple échelle ajustée *a posteriori* : elle est prédictible à partir des invariants du graphe.

---

## 10. Taxonomie unifiée des galaxies

L’article 14 propose enfin une taxonomie unifiée qui combine :

1. la phase de dimension LOW/HIGH,
2. la catégorie de qualité d’ajustement,
3. la classe de cohérence optimale définie par $f_\lambda^{\rm opt}$.

La répartition LOW/HIGH est la suivante :

$$
N_{\rm HIGH}=88, \qquad N_{\rm LOW}=87.
$$

avec :

$$
d_{\min}^{\rm HIGH}=2,487107, \qquad d_{\min}^{\rm LOW}=2,274448.
$$

Les classes de cohérence de l’article 14 sont :

```text
courte_0p5Rd          : 89 galaxies
standard_1Rd          : 38 galaxies
étendue_1p5_2Rd       : 39 galaxies
très_étendue_3Rd      : 9 galaxies

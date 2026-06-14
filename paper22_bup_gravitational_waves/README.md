# Paper 22 — Ondes gravitationnelles dans BuP

## Titre

**Ondes gravitationnelles dans BuP : modes informationnels du graphe d'intrication**

Version anglaise :

**Gravitational Waves in Bottom-Up Quantum Gravity: Information Modes of the Entanglement Graph**

---

## Objectif

Ce papier étudie les ondes gravitationnelles dans le cadre de la gravité quantique Bottom-Up (BuP).

Dans la relativité générale, une onde gravitationnelle est décrite comme une perturbation de la métrique :

\[
g_{\mu\nu} \;\rightarrow\; g_{\mu\nu} + h_{\mu\nu}.
\]

Dans BuP, la perturbation est plus profonde. Elle part du graphe d'intrication :

\[
W_{ij} = I(i : j),
\]

puis induit une perturbation de distance informationnelle, puis une perturbation métrique effective :

\[
\delta W_{ij}
\;\rightarrow\;
\delta d_{ij}^{\rm ent}
\;\rightarrow\;
\delta g_{\mu\nu}^{\rm eff}
\;\rightarrow\;
h_{\mu\nu}^{\rm eff}.
\]

Ainsi, les ondes gravitationnelles BuP sont interprétées comme des **modes collectifs informationnels** du réseau d'intrication.

---

## Résumé des prototypes numériques

### v2 — Polarisations \(+\) et \(\times\)

**Objectif :** vérifier que les deux polarisations gravitationnelles peuvent être reconstruites comme modes quadrupolaires de la perturbation du graphe.

**Résultat :** les deux modes sont séparés proprement dans une jauge fixe.

\[
\delta W_{ij}
\;\rightarrow\;
q_+(t),\; q_\times(t).
\]

---

### v3 — Propagation imposée

**Objectif :** imposer une perturbation de type onde plane :

\[
\delta W_{ij}(x,t) \sim A \cos(kx - \omega t),
\]

et vérifier que le mode quadrupolaire se propage avec le bon nombre d'onde.

**Résultat :** la pente de phase reconstruite vérifie :

\[
\frac{k_{\rm fit}}{k_{\rm input}} \simeq 1.
\]

---

### v4 — Source dynamique locale

**Objectif :** remplacer l'onde imposée par une source locale oscillante du graphe.

**Résultat pour la polarisation \(+\) :**

\[
v_+^{\rm arrivée} = 1,0036,
\qquad
R^2_{\rm arrivée} = 0,99993.
\]

**Résultat pour la polarisation \(\times\) :**

\[
v_\times^{\rm arrivée} = 1,0036,
\qquad
R^2_{\rm arrivée} = 0,99980.
\]

**Conclusion :**

\[
v_+ \simeq v_\times \simeq 1.
\]

Après calibration des unités du graphe :

\[
v_+ = v_\times = c.
\]

---

### v5 — Modes propres du Hessien de l'action BuP

**Objectif :** passer d'une simulation cinématique à une dynamique dérivée de l'action effective.

On linéarise autour d'un graphe d'équilibre :

\[
W_{ij}^{(0)}.
\]

On introduit une déformation nodale :

\[
W_{ij}(\phi) = W_{ij}^{(0)} \exp\left(\frac{\phi_i + \phi_j}{2}\right).
\]

Puis on calcule le Hessien :

\[
H_{ab} = \left. \frac{\partial^2 S_{\rm BuP}}{\partial \phi_a \partial \phi_b} \right|_{\phi = 0}.
\]

Les modes propres vérifient :

\[
H \phi_n = \omega_n^2 \phi_n.
\]

**Résultat :**

\[
\omega^2 \simeq c_{\rm graphe}^2 k^2 + m_{\rm eff}^2.
\]

Pour \(\eta_{\rm lisse} = 1,0\) :

\[
c_{\rm graphe} = 1,424,
\qquad
R^2 = 0,99962.
\]

Un scan de \(\eta_{\rm lisse}\) montre que le point naturel :

\[
\eta_{\rm lisse} \simeq 0,5
\]

donne :

\[
c_{\rm graphe} = 1,015,
\qquad
R^2 = 0,99855.
\]

**Ainsi, la vitesse de propagation est contrôlée par la rigidité informationnelle du graphe.**

---

## Interprétation physique

Dans BuP, une onde gravitationnelle n'est pas primitivement une ondulation d'un espace-temps déjà donné. C'est une **onde informationnelle du graphe d'intrication**, dont la projection géométrique est une onde de courbure.

\[
\delta W_{ij}
\;\rightarrow\;
\delta d_{ij}^{\rm ent}
\;\rightarrow\;
h_{\mu\nu}^{\rm eff}.
\]

Le graviton BuP n'est donc pas une particule fondamentale. Il est interprété comme le **quantum d'un mode collectif du réseau d'intrication**.

---

## Prédictions exploratoires

Dans la limite lisse actuelle, le graphe est au point fixe Einstein :

\[
c_{\rm graphe} = c.
\]

Des écarts éventuels :

\[
c_{\rm graphe} \neq c
\]

seraient des signatures hors limite lisse, associées à un secteur dynamique du tenseur d'émergence :

\[
\mathcal H_{\mu\nu}^{\rm dyn} \neq 0.
\]

Une hypothèse exploratoire est que la vitesse effective des ondes gravitationnelles pourrait sonder l'histoire informationnelle de l'univers :

\[
c_{\rm GW}^2(t) \sim \frac{\eta_{\rm lisse}(t)}{\rho_{\rm ent}(t)}.
\]

| Phase | Densité d'intrication | Vitesse \(c_{\rm GW}\) |
|-------|----------------------|----------------------|
| Univers primordial (fortement intriqué) | \(\rho_{\rm ent} \gg \rho_\ast\) | \(c_{\rm GW} < c\) |
| Univers actuel (lisse) | \(\rho_{\rm ent} = \rho_\ast\) | \(c_{\rm GW} = c\) |
| Univers futur (dilué) | \(\rho_{\rm ent} \ll \rho_\ast\) | \(c_{\rm GW} > c\) (apparent) |

Cette dernière possibilité ne doit pas être comprise comme une superluminalité ordinaire dans un espace-temps fixe, mais comme une **signature de perte progressive de cohérence géométrique**.

---

## Statut

Ce papier fournit une **preuve de concept numérique et variationnelle** :

- **v2** : extraction des polarisations ;
- **v3** : propagation imposée ;
- **v4** : source dynamique locale ;
- **v5** : modes propres du Hessien de l'action BuP effective.

Les résultats soutiennent l'idée que les ondes gravitationnelles BuP sont des **modes collectifs propagatifs du graphe d'intrication**.

---

## Limites

| Limitation | Explication |
|------------|-------------|
| v5 utilise une réduction nodale | Non encore le Hessien complet en variables d'arêtes \(W_{ij}\) |
| Graphes idéalisés | Les simulations sont effectuées sur des graphes simples |
| Calibration effective | \(c_{\rm graphe} = c\) est calibré, pas dérivé ab initio |
| Cosmologie exploratoire | Les prédictions sur \(c_{\rm GW}(t)\) sont exploratoires |
| Dérivation complète nécessaire | Il faudrait le Hessien edge-level : \(\displaystyle \frac{\delta^2 S_{\rm BuP}}{\delta W_{ij}\delta W_{kl}}\) |

---

## En une phrase

> **Paper 22 montre que les ondes gravitationnelles sont des modes collectifs du graphe d'intrication, avec deux polarisations, une vitesse \(c\), et une possible sonde de l'histoire cosmique de l'intrication.**

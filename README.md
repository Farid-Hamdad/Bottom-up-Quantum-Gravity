
# Bottom-Up Quantum Gravity

**Émergence de l'espace, du temps et de la gravité à partir de l'intrication quantique**  
**Farid Hamdad — Février 2026**

---

## Idée centrale

L'espace, le temps et la gravité ne sont pas fondamentaux.  
Ils émergent collectivement de la structure d'intrication d'un état quantique global fini.

Dans cette approche :

- l'intrication définit la connectivité,
- la géométrie émerge de la structure informationnelle,
- la gravité apparaît comme une thermodynamique — puis une courbure discrète — de l'intrication.

---

## Pourquoi ce projet ?

La physique moderne décrit avec une précision remarquable :

- la physique quantique,
- la gravitation classique,

mais laisse ouverte une question fondamentale :

**Pourquoi l'espace-temps existe-t-il ?**

et pourquoi la gravité possède-t-elle simultanément une structure :

- géométrique,
- thermodynamique ?

Ce projet explore une hypothèse minimale :

> **L'espace-temps n'est pas le théâtre de la physique.**  
> **Il est reconstruit à partir de l'intrication quantique.**

---

## 1. Postulat minimal

Il existe un état quantique global pur

\[
\Psi \in \bigotimes_i \mathcal{H}_i
\]

défini sur \(N\) degrés de liberté élémentaires (qubits),

sans :

- espace préalable,
- temps préalable,
- métrique préalable.

Tout le reste doit émerger :

- espace,
- temps,
- dimension,
- géométrie,
- gravité effective.

---

## 2. Méthodologie d'émergence

### 2.1 Émergence du temps — flot modulaire

Pour un sous-système \(A\) :

\[
\rho_A = \mathrm{Tr}_{\bar A} |\Psi\rangle\langle\Psi|
\]

\[
K_A = -\log(\rho_A)
\]

Le flot modulaire

\[
O(\tau) = e^{iK_A\tau} O e^{-iK_A\tau}
\]

définit une dynamique relationnelle intrinsèque.

**Lien conceptuel :** mécanisme de Page–Wootters.

> Le temps devient une propriété informationnelle interne.

---

### 2.2 Chaos modulaire

On analyse le spectre du Hamiltonien modulaire \(K_A\).

Statistique des gaps :

\[
\langle r \rangle
=
\left\langle
\frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})}
\right\rangle
\]

Valeurs universelles :

| Régime | Valeur |
|---|---:|
| Poisson | \(\approx 0.386\) |
| GOE | \(\approx 0.536\) |
| GUE | \(\approx 0.603\) |

**Résultat obtenu :**

\[
\langle r \rangle \in [0.53,\,0.59]
\]

> Signature compatible avec un chaos quantique de type Random Matrix Theory.

---

### 2.3 Spectral Form Factor

\[
g_2(t) = \frac{1}{d_A^2}\left|\sum_n e^{-it\tilde{\kappa}_n}\right|^2
\]

Structure observée :

- dip,
- ramp,
- plateau.

Scaling universel :

\[
g_{2,\mathrm{plateau}} \sim \frac{1}{d_A}
\]

---

### 2.4 Constante modulaire topologique

On définit :

\[
C = d_A \times g_{2,\mathrm{plateau}}
\]

Résultats typiques (\(d_A = 256\)) :

| Topologie | \(\langle r \rangle\) | \(C\) |
|---|---:|---:|
| chaîne 1D | 0.594 | 1.21 |
| grille 3×6 | 0.575 | 1.42 |
| graphe ER | 0.528 | 1.63 |

> La géométrie d'intrication influence la dynamique modulaire.

---

## 3. Émergence de l'espace

La géométrie est reconstruite à partir de l'information mutuelle :

\[
I(i:j)=S(\rho_i)+S(\rho_j)-S(\rho_{ij})
\]

Distance informationnelle :

\[
d_{ij}=-\log\left(\frac{I(i:j)}{I_{\max}}+\varepsilon\right)
\]

Puis :

- construction d'une matrice de distance,
- reconstruction par **MDS**,
- estimation de la dimension émergente.

La dimension effective est la plus petite dimension stabilisant l'erreur de reconstruction.

---

## 4. Résultats principaux

### 4.1 Dimension émergente

| Configuration | Intrication dominante | Dimension |
|---|---|---:|
| \(N=9,\ \lambda \approx 0\) | locale | \(d \approx 2\) |
| \(N=9,\ \lambda \to 1\) | non-locale | \(d \approx 3\) |
| \(N=16,\ \lambda \approx 0\) | locale | \(d \approx 2\) |
| \(N=16,\ \lambda \to 1\) | non-locale | \(d \approx 3\) |

### 4.2 ER = EPR opérationnel

Des qubits éloignés topologiquement deviennent proches géométriquement lorsque l'intrication non-locale augmente.

> Signature discrète de type wormhole-like.

### 4.3 Gravité thermodynamique

Test analogue à Jacobson :

\[
\delta S \approx \beta_{\mathrm{eff}}\delta E
\]

Relation stable observée pour :

- \(N=9\)
- \(N=16\)

---

## 5. Détection d'horizon bottom-up (\(N=16\))

On introduit un benchmark basé sur le graphe d'intrication.

Poids :

\[
W_{ij}=I(i:j)
\]

Seuil de densité :

\[
\rho = \frac{1}{3}
\]

Score :

\[
\mathrm{Score}_{BH}
=
w_S z(S)-w_\phi z(\phi)+w_{\mathrm{int}} z(\mathrm{internal})
\]

où :

- \(S\) = entropie de région,
- \(\phi\) = conductance,
- `internal` = intrication interne.

### Région BH-like détectée

\[
[0,2,3,5,6,10,11,15]
\]

| Quantité | Valeur |
|---|---:|
| entropie | 7.1667 |
| cut | 3.2098 |
| internal | 5.0348 |
| conductance | 0.3951 |
| RMT ratio | 0.6041 |

### Signature

- intrication interne très élevée,
- couplage extérieur faible,
- bottleneck informationnel.

> Un horizon peut émerger sans géométrie préalable.

---

## 6. Géométrie multipartite et limite du diagnostic triangulaire

Nous avons étudié la structure multipartite via l'information mutuelle conditionnelle :

\[
\mathrm{CMI}(i:j|k)
\]

et des descripteurs triangulaires de l'intrication.

### Résultat intermédiaire

Des corrélations significatives existent entre CMI et structures triangulaires locales de l'information mutuelle, ce qui indique qu'une part de la structure multipartite est bien encodée dans la géométrie relationnelle du graphe d'intrication.

### Mais

Les tests de **courbure triangulaire** et d'**holonomie discrète** (v1–v4) n'ont **pas** fourni de signal robuste de courbure locale.

> **Conclusion actuelle :**  
> la structure triangulaire est un bon descripteur géométrique local de la CMI,  
> mais **pas** encore une mesure robuste de courbure gravitationnelle émergente.

---

## 7. Nouveau résultat : courbure discrète d'Ollivier–Ricci

Pour dépasser les limites du test triangulaire, nous avons introduit une courbure discrète de type **Ollivier–Ricci** sur le graphe d'intrication.

### Objectif

Tester directement si une perturbation informationnelle locale agit comme une source de courbure :

\[
\text{défaut informationnel local}
\;\Longrightarrow\;
\text{courbure géométrique émergente}
\]

### Protocole

- graphe d'intrication reconstruit à partir de l'information mutuelle,
- perturbation locale contrôlée sur une région source,
- mesure de la variation de courbure :

\[
\Delta \kappa = \kappa_{\mathrm{pert}} - \kappa_{\mathrm{base}}
\]

- comparaison entre zone **proche** et **lointaine** :

\[
\Delta \kappa_{\mathrm{near-far}}
=
\Delta \kappa_{\mathrm{near}}
-
\Delta \kappa_{\mathrm{far}}
\]

### Résultat robuste

Le mode le plus propre est le mode **`frozen_baseline`**, où la connectivité du graphe de référence est maintenue fixe.

Batch multi-seeds (\(N=16\), 20 seeds, rayon \(r=1\), force \(s=0.6\)) :

\[
\Delta \kappa_{\mathrm{edge}} = 0.0758 \pm 0.0296
\]

avec :

- fraction positive = **1.00**
- \(p_{\mathrm{sign}} = 1.9\times10^{-6}\)

et

\[
\Delta \kappa_{\mathrm{near-far}} = 0.0927 \pm 0.0799
\]

avec :

- fraction positive = **0.95**
- \(p_{\mathrm{sign}} = 4.0\times10^{-5}\)

> Une perturbation informationnelle locale augmente la courbure moyenne,  
> et cette augmentation est plus forte près de la source qu'à distance.

### Loi de réponse en intensité

Pour un défaut compact (\(r=1\)), la courbure moyenne croît avec la force \(s\) :

| \(s\) | \(\Delta \kappa_{\mathrm{edge}}\) |
|---:|---:|
| 0.0 | 0.000 |
| 0.2 | 0.010 |
| 0.4 | 0.040 |
| 0.6 | 0.076 |
| 0.8 | 0.096 |
| 1.0 | 0.109 |

La localité suit aussi une croissance puis saturation :

| \(s\) | \(\Delta \kappa_{\mathrm{near-far}}\) |
|---:|---:|
| 0.0 | 0.000 |
| 0.2 | 0.031 |
| 0.4 | 0.071 |
| 0.6 | 0.093 |
| 0.8 | 0.099 |
| 1.0 | 0.095 |

### Effet du rayon

Le scan en rayon montre :

- **rayon petit** : réponse plus localisée,
- **rayon grand** : réponse plus diffuse.

Autrement dit :

> une source compacte courbe plus localement,  
> une source étendue diffuse la déformation géométrique.

### Forme effective

Le meilleur ajustement n'est pas strictement linéaire mais **quadratique concave** :

\[
\Delta \kappa(s)\approx as^2+bs+c,\qquad a<0
\]

Cela suggère une croissance initiale suivie d'un tassement ou d'une saturation.

### Conclusion physique

C'est à ce stade le signal le plus net obtenu dans le projet en faveur d'un couplage :

\[
\text{matière informationnelle locale}
\;\Longrightarrow\;
\text{courbure discrète émergente}
\]

---

## 8. Ce que le projet établit actuellement

- une géométrie peut émerger d'un état quantique fini,
- la dimension émergente dépend de la structure d'intrication,
- ER = EPR devient mesurable dans un cadre discret,
- une thermodynamique de l'intrication apparaît,
- le flot modulaire présente une signature chaotique,
- le plateau du SFF suit le scaling \(1/d_A\),
- la constante modulaire dépend de la topologie,
- un horizon informationnel peut émerger,
- la structure multipartite est partiellement encodée dans la géométrie relationnelle,
- une perturbation informationnelle locale peut induire une **courbure d'Ollivier–Ricci** croissante et localisée.

---

## 9. Limites actuelles

- limite continue \(N\to\infty\) non démontrée,
- pas de dynamique relativiste complète,
- pas d'équation d'Einstein discrète dérivée,
- pas encore de prédictions observationnelles directes,
- la courbure triangulaire / holonomie discrète n'a pas donné de signal robuste,
- le couplage matière-courbure est actuellement établi comme **preuve de concept numérique**, pas encore comme loi fondamentale dérivée.

---

## Organisation du dépôt

```text
paper/
  Bottom-up_Quantum_Gravity.tex
  Bottom-up_Quantum_Gravity.pdf

theory/
  theory_mathematics.tex

figures/

scripts/
  bh_benchmark_louvain_N16.py
  bup_ollivier_ricci_local_response_v1_2.py
  bup_ollivier_ricci_local_response_batch_v1_3.py
  bup_ollivier_ricci_local_response_scan_v1_4.py
  bup_ollivier_ricci_local_response_scan_radius_fit_v1_5.py

code/
data/

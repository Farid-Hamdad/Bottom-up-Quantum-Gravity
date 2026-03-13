
# Bottom-Up Quantum Gravity

**Émergence de l’espace, du temps et de la gravité à partir de l’intrication quantique**  
**Farid Hamdad — Février 2026**

---

## Vue d’ensemble

**Bottom-Up Quantum Gravity** explore une hypothèse minimale :

> l’espace, le temps et la gravité ne sont pas fondamentaux ;  
> ils émergent collectivement de la structure d’intrication d’un état quantique global fini.

Dans cette approche :

- l’intrication définit la connectivité,
- la géométrie émerge de la structure informationnelle,
- le temps apparaît via le flot modulaire,
- la gravité peut être approchée comme une thermodynamique — puis comme une courbure discrète — de l’intrication.

Le projet combine ainsi trois niveaux :

- une base conceptuelle issue de l’information quantique,
- des tests numériques sur systèmes finis,
- une première extension phénoménologique vers les courbes de rotation galactiques (**SPARC**).

---

## Idée centrale

On part d’un état quantique global pur

\[
\Psi \in \bigotimes_i \mathcal{H}_i
\]

défini sur \(N\) degrés de liberté élémentaires (qubits),

sans :

- espace préalable,
- temps préalable,
- métrique préalable.

L’objectif est alors de reconstruire, à partir des seules corrélations internes de l’état :

- une notion d’espace,
- une dynamique relationnelle,
- une dimension effective,
- une géométrie émergente,
- et une réponse gravitationnelle effective.

---

## Pourquoi ce projet ?

La physique moderne décrit avec succès :

- la mécanique quantique,
- la gravitation classique,

mais laisse ouverte une question plus profonde :

> **Pourquoi l’espace-temps existe-t-il ?**

et pourquoi la gravité semble-t-elle posséder à la fois une nature :

- géométrique,
- thermodynamique,
- informationnelle ?

Le programme bottom-up testé ici consiste à ne pas postuler l’espace-temps, mais à chercher à le **reconstruire** à partir de l’intrication.

---

## 1. Temps émergent : flot modulaire

Pour un sous-système \(A\),

\[
\rho_A = \mathrm{Tr}_{\bar A} |\Psi\rangle\langle\Psi|
\]

\[
K_A = -\log(\rho_A)
\]

Le flot modulaire

\[
O(\tau)=e^{iK_A\tau} O e^{-iK_A\tau}
\]

définit une dynamique relationnelle intrinsèque.

**Interprétation :**  
le temps n’est plus un paramètre de fond, mais une propriété informationnelle interne du système.

---

## 2. Diagnostics spectraux et chaos modulaire

Le spectre du Hamiltonien modulaire \(K_A\) est étudié à l’aide :

- des statistiques de gaps adjacents,
- des diagnostics de type Random Matrix Theory,
- du spectral form factor.

La quantité

\[
\langle r \rangle
=
\left\langle
\frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})}
\right\rangle
\]

permet de comparer le spectre observé aux classes universelles connues.

Valeurs de référence :

| Régime | Valeur |
|---|---:|
| Poisson | \(\approx 0.386\) |
| GOE | \(\approx 0.536\) |
| GUE | \(\approx 0.603\) |

Résultats typiques obtenus dans le projet :

\[
\langle r \rangle \in [0.53,\,0.59]
\]

ce qui est compatible avec un régime de chaos modulaire non trivial.

Le spectral form factor présente également la structure classique :

- dip,
- ramp,
- plateau,

avec un plateau suivant approximativement le scaling

\[
g_{2,\mathrm{plateau}} \sim \frac{1}{d_A}.
\]

---

## 3. Émergence de l’espace à partir de l’information mutuelle

La géométrie est reconstruite à partir de l’information mutuelle

\[
I(i:j)=S(\rho_i)+S(\rho_j)-S(\rho_{ij})
\]

et d’une distance informationnelle de type

\[
d_{ij}=-\log\left(\frac{I(i:j)}{I_{\max}}+\varepsilon\right).
\]

À partir de là, le projet construit :

- des matrices de distance,
- des plongements géométriques par **MDS**,
- des graphes d’intrication,
- des estimations de dimension effective.

L’idée centrale est la suivante :

> la géométrie n’est pas imposée ;  
> elle est reconstruite depuis la structure des corrélations.

---

## 4. Résultats internes principaux

### 4.1 Dimension émergente

Les expériences numériques montrent que la dimension effective dépend de la structure d’intrication.

| Configuration | Intrication dominante | Dimension émergente |
|---|---|---:|
| \(N=9,\ \lambda \approx 0\) | locale | \(d \approx 2\) |
| \(N=9,\ \lambda \to 1\) | non-locale | \(d \approx 3\) |
| \(N=16,\ \lambda \approx 0\) | locale | \(d \approx 2\) |
| \(N=16,\ \lambda \to 1\) | non-locale | \(d \approx 3\) |

Cela soutient l’idée que la dimension n’est pas un postulat, mais une grandeur émergente dépendant du motif d’intrication.

---

### 4.2 Signature opérationnelle ER = EPR

Lorsque l’intrication non-locale augmente, des qubits éloignés topologiquement peuvent devenir proches dans la géométrie reconstruite.

> C’est une signature discrète de type *wormhole-like*, compatible avec l’idée opérationnelle ER = EPR.

---

### 4.3 Gravité thermodynamique

Un test analogue à la relation de Jacobson

\[
\delta S \approx \beta_{\mathrm{eff}}\delta E
\]

présente une stabilité numérique sur des systèmes finis comme \(N=9\) et \(N=16\).

L’interprétation actuelle est qu’une thermodynamique de l’intrication peut jouer le rôle de précurseur d’une dynamique gravitationnelle effective.

---

## 5. Horizon bottom-up (\(N=16\))

Un benchmark de détection de région de type horizon a été introduit sur le graphe d’intrication.

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
| ratio RMT | 0.6041 |

Cette région combine :

- une intrication interne élevée,
- un couplage extérieur relativement faible,
- un bottleneck informationnel.

> Une structure de type horizon peut donc émerger sans géométrie préalable.

---

## 6. Géométrie multipartite : résultat utile et limite actuelle

Le projet a étudié la structure multipartite via :

- l’information mutuelle conditionnelle \(\mathrm{CMI}(i:j|k)\),
- des descripteurs triangulaires locaux,
- des tests de courbure triangulaire et d’holonomie discrète.

### Résultat intermédiaire

Une partie de la structure multipartite est bien corrélée à la géométrie relationnelle locale du graphe d’intrication.

### Limite actuelle

Les tests de courbure triangulaire et d’holonomie discrète (v1–v4) n’ont pas fourni de signal robuste de courbure gravitationnelle locale.

> Conclusion provisoire :  
> la structure triangulaire est un bon descripteur géométrique local de la CMI,  
> mais pas encore une mesure robuste de courbure émergente.

---

## 7. Nouveau résultat : courbure discrète d’Ollivier–Ricci

Pour dépasser les limites des diagnostics triangulaires, le projet a introduit une courbure discrète de type **Ollivier–Ricci** sur le graphe d’intrication.

### Objectif

Tester si une perturbation informationnelle locale agit comme une source de courbure émergente :

\[
\text{défaut informationnel local}
\;\Longrightarrow\;
\text{courbure géométrique émergente}
\]

### Protocole

- reconstruction du graphe à partir de l’information mutuelle,
- perturbation locale contrôlée sur une région source,
- mesure de

\[
\Delta \kappa = \kappa_{\mathrm{pert}} - \kappa_{\mathrm{base}}
\]

- comparaison entre zone proche et zone lointaine :

\[
\Delta \kappa_{\mathrm{near-far}}
=
\Delta \kappa_{\mathrm{near}}
-
\Delta \kappa_{\mathrm{far}}
\]

### Résultat robuste

Le mode le plus propre est **`frozen_baseline`**, où la connectivité de référence est maintenue fixe.

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

### Interprétation

Une perturbation informationnelle locale :

- augmente la courbure moyenne,
- et cet effet est plus marqué près de la source qu’à distance.

### Loi de réponse en intensité

Pour un défaut compact (\(r=1\)) :

| \(s\) | \(\Delta \kappa_{\mathrm{edge}}\) |
|---:|---:|
| 0.0 | 0.000 |
| 0.2 | 0.010 |
| 0.4 | 0.040 |
| 0.6 | 0.076 |
| 0.8 | 0.096 |
| 1.0 | 0.109 |

Indicateur de localité :

| \(s\) | \(\Delta \kappa_{\mathrm{near-far}}\) |
|---:|---:|
| 0.0 | 0.000 |
| 0.2 | 0.031 |
| 0.4 | 0.071 |
| 0.6 | 0.093 |
| 0.8 | 0.099 |
| 1.0 | 0.095 |

### Effet du rayon

Le scan en rayon indique que :

- un rayon petit produit une réponse plus localisée,
- un rayon grand diffuse davantage la déformation géométrique.

> Une source compacte courbe plus localement ;  
> une source étendue diffuse la réponse géométrique.

À ce stade, c’est le signal le plus net obtenu dans le projet en faveur d’un couplage :

\[
\text{matière informationnelle locale}
\;\Longrightarrow\;
\text{courbure discrète émergente}
\]

---

## 8. Extension phénoménologique : résultats SPARC

Le projet ne se limite plus à des reconstructions géométriques internes sur petits systèmes quantiques.  
Il commence aussi à tester si ses idées peuvent être prolongées vers des données astrophysiques réelles.

Dans cette optique, une première étude a été menée sur l’échantillon **SPARC** de courbes de rotation galactiques.

### Objectif

Tester si une géométrie émergente effective peut fournir une lecture phénoménologique pertinente de la dynamique galactique.

L’enjeu n’est pas encore de proposer une théorie observationnelle finale, mais de voir si le programme bottom-up peut relier :

- dimension émergente,
- réponse gravitationnelle effective,
- dynamique de rotation des galaxies.

---

### 8.1 Résultat principal : séparation de phases effectives

L’analyse actuelle suggère qu’une description unique n’est pas la plus naturelle.  
Le pipeline de diagramme de phase met en évidence **deux régimes effectifs** :

| Régime galactique | Dimension effective optimale |
|---|---:|
| phase massive | \(d_{\min}=2.4871 \pm 0.0294\) |
| phase naine | \(d_{\min}=2.2744 \pm 0.0341\) |

Diagnostics actuels issus du pipeline :

| Régime galactique | test \(\chi^2_{\mathrm{red}}\) |
|---|---:|
| phase massive | \(10.9353 \pm 12.6952\) |
| phase naine | \(25.2842 \pm 25.8802\) |

### 8.2 Lecture physique provisoire

Ces résultats suggèrent que la réponse gravitationnelle effective pourrait dépendre du **régime galactique**, et qu’une géométrie émergente peut conduire à des **dimensions effectives différentes** selon les classes de galaxies.

Dans l’état actuel du projet, cela se lit comme :

- une indication contre une lecture strictement monophasique,
- un soutien à une **géométrie émergente multi-régimes**,
- une première ouverture phénoménologique du programme bottom-up au-delà des systèmes quantiques finis.

### 8.3 Statut scientifique

Cette partie SPARC doit être lue comme :

- une **preuve de concept phénoménologique**,
- un **test de cohérence empirique exploratoire**,
- une étape vers une confrontation plus directe avec les observations.

Elle ne constitue pas encore :

- une théorie complète de la dynamique galactique,
- un modèle observationnel entièrement contraint,
- ni un remplacement des cadres standards d’ajustement astrophysique.

En revanche, elle marque une étape importante :

> le programme bottom-up ne reconstruit plus seulement une géométrie interne ;  
> il commence à être confronté à des données astrophysiques réelles.

---

## 9. Ce que le projet établit actuellement

À son stade actuel, le projet fournit des indices numériques en faveur des idées suivantes :

- une géométrie peut émerger d’un état quantique fini,
- la dimension émergente dépend de la structure d’intrication,
- ER = EPR devient testable dans un cadre discret,
- le flot modulaire fournit une dynamique interne non triviale,
- les spectres modulaires présentent une signature chaotique,
- le plateau du spectral form factor suit un scaling compatible avec \(1/d_A\),
- une structure de type horizon peut émerger,
- la structure multipartite est partiellement encodée dans la géométrie relationnelle,
- une perturbation informationnelle locale peut induire une courbure discrète d’Ollivier–Ricci,
- une première extension phénoménologique vers les galaxies SPARC suggère une description en **régimes effectifs distincts**.

---

## 10. Limites actuelles

Le dépôt présente actuellement un programme de **preuve de concept numérique**.  
Des limites importantes restent ouvertes :

- la limite continue \(N \to \infty\) n’est pas établie,
- il n’existe pas encore de dynamique relativiste complète,
- aucune équation d’Einstein discrète n’a encore été dérivée,
- la partie SPARC reste exploratoire,
- la courbure triangulaire / holonomie discrète n’a pas donné de signal robuste,
- le couplage matière–courbure est pour l’instant soutenu numériquement, pas encore dérivé de premiers principes.

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

experiments/
  sparc/
    results_phase_diagram/

code/
data/
Résultats SPARC actuellement générés

Le pipeline SPARC écrit notamment les fichiers suivants :

experiments/sparc/results_phase_diagram/phase_diagram_summary.json
experiments/sparc/results_phase_diagram/phase_table.csv
experiments/sparc/results_phase_diagram/phase_diagram.png
Ordre de lecture recommandé

Pour découvrir le projet :

README.md — vue d’ensemble

paper/Bottom-up_Quantum_Gravity.tex — présentation conceptuelle

theory/theory_mathematics.tex — base théorique plus formelle

scripts/ — protocoles numériques

experiments/ — résultats générés, en particulier SPARC

Statut du projet

Statut : recherche active / preuve de concept numérique

Les axes actuels sont :

consolider le couplage matière–courbure,

améliorer les analyses de scaling,

renforcer le lien entre dimension émergente et phénoménologie,

étendre les tests observationnels autour des résultats SPARC.

Auteur

Farid Hamdad
hamdadfarid54@gmail.com

Projet indépendant de recherche sur la géométrie émergente, la structure d’intrication et la gravité quantique bottom-up.

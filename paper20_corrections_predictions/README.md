# Paper 20 — Corrections, Prédictions et Déviation Observables

**Au-delà de la limite einsteinienne de Bottom-Up Quantum Gravity**

---

## Statut

Paper 19 a assemblé l'équation d'Einstein effective de Bottom-Up Quantum Gravity :

$$
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}.
$$

Paper 20 étudie le tenseur de correction :

$$
\mathcal{H}_{\mu\nu}.
$$

Ce tenseur contient les écarts de BuP par rapport à la limite einsteinienne lisse. Dans le régime continu de basse énergie,

$$
\mathcal{H}_{\mu\nu}\to0,
$$

et BuP se réduit à une équation d'Einstein effective. En dehors de ce régime, $\mathcal{H}_{\mu\nu}\neq0$ et BuP prédit des déviations observables.

---

## 1. Tenseur de correction

Paper 20 décompose :

$$
\mathcal{H}_{\mu\nu}
=
\mathcal{H}_{\mu\nu}^{\rm spec}
+
\mathcal{H}_{\mu\nu}^{\rm curv}
+
\mathcal{H}_{\mu\nu}^{\rm source}
+
\mathcal{H}_{\mu\nu}^{\rm dim}
+
\mathcal{H}_{\mu\nu}^{\rm nonlocal}
+
\mathcal{H}_{\mu\nu}^{\rm topo}
+
\mathcal{H}_{\mu\nu}^{\rm finite}.
$$

Chaque secteur a une signification physique distincte :

| Secteur | Signification | Cible principale |
|---|---|---|
| $\mathcal{H}^{\rm spec}_{\mu\nu}$ | résidus spectraux / du noyau de chaleur | corrections de courbure supérieure |
| $\mathcal{H}^{\rm curv}_{\mu\nu}$ | corrections d'étalonnage OR vers Ricci | étalonnage de la courbure |
| $\mathcal{H}^{\rm source}_{\mu\nu}$ | résidus de la source modulaire | reconstruction du tenseur énergie-impulsion |
| $\mathcal{H}^{\rm dim}_{\mu\nu}$ | corrections de dimension effective | SPARC / cosmologie |
| $\mathcal{H}^{\rm nonlocal}_{\mu\nu}$ | corrections d'intrication non locale / liens longs | lentille gravitationnelle / SLACS |
| $\mathcal{H}^{\rm topo}_{\mu\nu}$ | topologie globale / constantes modulaires | holonomies / topologie |
| $\mathcal{H}^{\rm finite}_{\mu\nu}$ | corrections à $N$ fini | mise à l'échelle numérique |

---

## 2. Secteurs de correction quantifiés

Le budget de correction v3 relie les résultats numériques déjà validés des Papers 16 à 19 aux branches phénoménologiques SPARC et lentille gravitationnelle.

| Secteur | Moyenne absolue de la grandeur proxy | Maximum absolu de la grandeur proxy | Statut |
|---|---:|---:|---|
| $H_{\rm spec}$ | 0.015850 | 0.026627 | borné numériquement |
| $H_{\rm source\_modular}$ | 0.010400 | 0.010519 | faible résidu |
| $H_{\rm source\_curvature}$ | 0.025260 | 0.032771 | faible résidu |
| $H_{\rm dim}$ | 8.508743 | 12.568829 | relié à la phénoménologie |
| $H_{\rm nonlocal}$ | 445.904384 | 889.087288 | relié à la phénoménologie |
| $H_{\rm localization}$ | 7.913268 | 9.533529 | localisation positive |

Les grandes valeurs dans $H_{\rm dim}$ et $H_{\rm nonlocal}$ ne sont pas des normes tensorielles. Ce sont des grandeurs proxy phénoménologiques : les résidus de SPARC pour $H_{\rm dim}$, et les amplitudes du pont de lentille gravitationnelle pour $H_{\rm nonlocal}$.

---

## 3. Corrections spectrales

Le secteur spectral provient du fait que

$$
\mathrm{Tr}\,L^{-\beta}
$$

ne se réduit pas uniquement au terme d'Einstein-Hilbert. Par le développement du noyau de chaleur,

$$
\mathrm{Tr}(e^{-t\Delta})
\sim
(4\pi t)^{-D/2}
\left(
a_0+a_1t+a_2t^2+\cdots
\right),
$$

les coefficients supérieurs génèrent des corrections du type :

$$
R^2,
\qquad
R_{\mu\nu}R^{\mu\nu},
\qquad
R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma},
\qquad
\nabla^2 R,
\qquad
\cdots.
$$

Dans le budget de correction, le résidu spectral est borné numériquement par Paper 16 :

$$
\langle |H_{\rm spec}|\rangle \simeq 0.015850,
$$

avec un maximum de :

$$
0.026627.
$$

Ainsi $H_{\rm spec}$ est faible sur les géométries contrôlées.

---

## 4. Corrections d'étalonnage de la courbure

Paper 17 a montré que la courbure d'Ollivier-Ricci porte un signal de Ricci, mais avec un étalonnage à échelle finie :

$$
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_N R_{\mu\nu}u^\mu u^\nu.
$$

À $N=512$,

$$
B_N=-0.301281,
$$

$$
C_N=0.391933,
$$

$$
A_N=\frac{1}{C_N}=2.551455.
$$

Ainsi $H_{\rm curv}$ n'est pas un échec du secteur de courbure. C'est la renormalisation à échelle finie requise pour convertir la courbure OR en une grandeur proxy pour le tenseur de Ricci.

---

## 5. Corrections de la source

Paper 18 a montré :

$$
\delta S_A^{\rm graphe}
\simeq
\delta\langle K_A^{\rm graphe}\rangle,
$$

avec des pentes de :

$$
0.989481
$$

sur le tore plat et de :

$$
0.989718
$$

sur la sphère.

Le résidu source-modulaire est donc approximativement :

$$
H_{\rm source\_modular}\sim 0.0104.
$$

Paper 18 a également montré que la source modulaire prédit la réponse de courbure près de la source, avec un résidu :

$$
H_{\rm source\_curvature}\sim0.0253.
$$

Cela fait du secteur source l'un des secteurs les mieux contrôlés dans le pipeline actuel de BuP.

---

## 6. Corrections dimensionnelles et SPARC

Le secteur dimensionnel est lié aux courbes de rotation SPARC.

Le budget de correction v3 donne :

$$
{\rm médiane}\ \chi^2_{\rm red}=12.568829,
$$

$$
{\rm fraction}(\chi^2_{\rm red}<2)=0.182857,
$$

$$
{\rm fraction}(\chi^2_{\rm red}<10)=0.485714.
$$

La taxonomie de qualité d'ajustement donne :

$$
f_{\rm excellent}=0.182857,
$$

$$
f_{\rm bon}=0.154286,
$$

$$
f_{\rm moyen}=0.274286,
$$

$$
f_{\rm faible}=0.388571.
$$

Par conséquent,

$$
f_{\rm excellent+bon+moyen}=0.611429.
$$

Environ $61.1\%$ des galaxies se situent au moins dans le régime moyen ou meilleur, tandis que $38.9\%$ restent dans un régime faible nécessitant des corrections morphologiques, de phase ou non locales supplémentaires.

La fenêtre dimensionnelle effective est :

$$
\langle d_{\rm utilisé,min}\rangle=2.023651,
$$

$$
\langle d_{\rm utilisé,max}\rangle=2.643653.
$$

Ainsi,

$$
\mathcal{H}_{\mu\nu}^{\rm dim}
=
\mathcal{H}_{\mu\nu}[d_s(x),d_w(x),\alpha_{\rm eff}(x)]
$$

est le secteur de correction qui encode le flux dimensionnel à l'échelle des galaxies.

---

## 7. Corrections non locales et lentille gravitationnelle

Le secteur non local est lié à la branche de lentille gravitationnelle de BuP.

Le budget de correction v3 actuel intègre le résumé du pont de lentille gravitationnelle et donne :

$$
\langle {\rm ratio}_{\rm pic}\rangle=2.721481,
$$

$$
{\rm médiane}({\rm ratio}_{\rm pic})=1.967831,
$$

$$
{\rm max}({\rm ratio}_{\rm pic})=13.807607.
$$

La réponse absolue au pic donne :

$$
\langle \Delta_{\rm pic}\rangle=889.087288,
$$

$$
{\rm médiane}(\Delta_{\rm pic})=704.916728,
$$

$$
{\rm max}(\Delta_{\rm pic})=3594.219422.
$$

La fraction de profils de pont de lentille gravitationnelle valides est :

$$
f_{\rm ok}=0.977143.
$$

Ainsi $H_{\rm nonlocal}$ est maintenant lié à la phénoménologie de la lentille gravitationnelle.

---

## 8. Branche historique de la lentille gravitationnelle / SLACS

Avant la formulation du budget de correction, BuP contenait déjà une branche de lentille gravitationnelle. Sa logique était :

$$
d(r)
\to
f_d(r)
\to
g_{\rm eff}(r)
\to
\alpha(R).
$$

Le facteur d'amplification dimensionnel était :

$$
f_d(r)=
\left(
\frac{r}{r_{\rm ref}}
\right)^{3-d(r)}.
$$

La réponse optique effective était modélisée comme :

$$
g_{\rm eff}(r)
=
g_{\rm bary}(r)
\left[
1+
\lambda w_{\rm tot}(r)(f_d(r)-1)
\right].
$$

Un régime optique centralisé de la forme

$$
(\lambda,\text{largeur},x_{\rm cap},r_{\rm cut})
=
(3.0,0.8,4.0,5.0)
$$

a produit des excès de déflexion robustes de l'ordre de $30\%$ à $50\%$ sur plusieurs spirales massives.

Cette branche historique est maintenant interprétée comme une première réalisation phénoménologique de :

$$
\mathcal{H}_{\mu\nu}^{\rm nonlocal/optique}.
$$

Le budget actuel de Paper 20 ne remplace pas les travaux antérieurs sur SLACS / lentille gravitationnelle. Il les reclassifie à l'intérieur du tenseur de correction :

$$
\boxed{
\text{branche lentille gravitationnelle / SLACS}
\subset
\mathcal{H}_{\mu\nu}^{\rm nonlocal}.
}
$$

---

## 9. Secteurs restants

Deux secteurs restent partiellement liés :

$$
H_{\rm topo}
$$

et

$$
H_{\rm finite}.
$$

Le secteur topologique devrait être connecté à :

- constantes modulaires ;
- holonomies ;
- boucles de type Wilson ;
- contraintes globales du graphe ;
- $S_{\rm topo}[W]$.

Le secteur fini devrait être connecté à :

- mise à l'échelle en $N$ ;
- mise à l'échelle de la densité ;
- mise à l'échelle de la résolution du graphe ;
- convergence de $d_s$, $\kappa^{OR}$ et des réponses modulaires.

---

## 10. Interprétation

Paper 20 montre que le tenseur de correction non-einsteinien n'est pas arbitraire. Il a une origine structurée :

$$
\mathcal{H}_{\mu\nu}^{\rm spec}
\quad
\text{provenant de la géométrie spectrale},
$$

$$
\mathcal{H}_{\mu\nu}^{\rm dim}
\quad
\text{provenant du flux dimensionnel effectif},
$$

$$
\mathcal{H}_{\mu\nu}^{\rm nonlocal}
\quad
\text{provenant de l'intrication à longue portée et de la réponse de lentille gravitationnelle},
$$

$$
\mathcal{H}_{\mu\nu}^{\rm source}
\quad
\text{provenant des résidus de la source modulaire},
$$

$$
\mathcal{H}_{\mu\nu}^{\rm topo}
\quad
\text{provenant de la topologie globale},
$$

$$
\mathcal{H}_{\mu\nu}^{\rm finite}
\quad
\text{provenant de la taille finie du graphe}.
$$

Les plus grands secteurs phénoménologiques actuellement liés sont $H_{\rm dim}$ et $H_{\rm nonlocal}$, qui correspondent respectivement aux courbes de rotation des galaxies et à la lentille gravitationnelle.

---

## 11. Résumé en une phrase

Paper 20 classifie le tenseur de correction $\mathcal{H}_{\mu\nu}$, relie le secteur dimensionnel à SPARC et le secteur non local à la lentille gravitationnelle, et identifie les déviations observables de BuP par rapport à la limite einsteinienne lisse.

---

*Dernière mise à jour : Juillet 2026*

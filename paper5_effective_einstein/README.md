# BuP Paper 5 — Équation d'Einstein effective émergente

**Auteur :** Farid Hamdad  
**Année :** 2026

---

## Résumé simple

En relativité générale, l’équation d’Einstein relie la courbure de l’espace-temps au contenu matériel de l’Univers.

Le programme **Bottom-Up Quantum Gravity (BuP)** propose une autre lecture : l’espace-temps et la gravité ne sont pas fondamentaux, mais **émergent d’un réseau d’intrication quantique**.

Les travaux précédents du programme ont exploré trois régimes :

- **Paper 2** : une dimension cosmologique variable \(d_{\rm bg}(z)\) permet de rendre compte de l’accélération tardive sans fluide d’énergie noire ajouté.
- **Paper 3** : le couplage gravitationnel effectif \(G_{\rm eff}(z)\) modifie la croissance linéaire et donne \(\sigma_8=0.772\).
- **Paper 4** : dans les régions denses, la dimension locale est plus faible, ce qui accélère l’effondrement des halos rares et fournit une piste pour l’excès JWST.

**Paper 5** rassemble ces résultats dans une équation effective unique :

- la métrique est une métrique émergente \(g_{\mu\nu}^{\rm ent}\) ;
- la force de gravité dépend de la dimension spectrale locale \(d_s\) ;
- un terme géométrique effectif dépend des gradients et variations de \(d_s\).

Nous avons également effectué un test numérique interne sur des réseaux d’intrication finis jusqu’à \(N=20\) :

- construction d’une courbure discrète ;
- construction d’un proxy de tenseur effectif à partir de \(d_s\) ;
- mesure de leur corrélation.

Le résultat montre une corrélation monotone robuste entre la courbure émergente et les variations locales de la dimension spectrale. C’est un premier test numérique interne de cohérence pour l’équation d’Einstein effective BuP.

---

## Équation centrale

$$
G_{\mu\nu}[g^{\rm ent}]
=
8\pi G_{\rm eff}(d_s) T_{\mu\nu}^{\rm matter}
+
8\pi G T_{\mu\nu}^{\rm ent}[d_s]
$$

avec

$$
G_{\rm eff}(d_s)
=
\frac{2}{d_s-1}G
$$

et

$$
T_{\mu\nu}^{\rm ent}[d_s]
=
\alpha\left(
\nabla_\mu d_s\nabla_\nu d_s
-
\frac12 g_{\mu\nu}^{\rm ent}(\nabla d_s)^2
\right)
+
\beta\left(
\nabla_\mu\nabla_\nu d_s
-
g_{\mu\nu}^{\rm ent}\Box d_s
\right).
$$

Aucune constante cosmologique \(\Lambda\) n’est postulée.  
Aucun fluide d’énergie sombre n’est ajouté.  
L’accélération tardive est interprétée comme une conséquence géométrique de la dimension effective variable \(d_{\rm bg}(z)<3\).

---

## Test numérique de cohérence interne

Nous comparons des proxys discrets de courbure et de tenseur effectif sur des graphes d’intrication BuP :

\[
N=9,12,16,20.
\]

La relation la plus stable pour \(N=20\), testée sur \(k=3,5,7\), est :

$$
G_{ij}^{\rm proxy}
=
-\kappa_{ij}
+
\frac12 R_{\rm edge}.
$$

Ce proxy de courbure est corrélé avec des proxys spectraux construits à partir de :

$$
(\nabla d_s)^2,
\qquad
\Delta d_s.
$$

| \(k\) | Proxy spectral dominant | \(\langle r \rangle\) | \(\langle\rho\rangle\) | Signe Spearman |
|---|---:|---:|---:|:---:|
| 3 | \((\nabla d_s)^2 + |\Delta d_s|\) | 0.707 | 0.598 | 100% |
| 5 | \((\nabla d_s)^2 + |\Delta d_s|\) | 0.558 | 0.566 | 100% |
| 7 | \((\nabla d_s)^2 - \Delta d_s\) | 0.995 | 0.417 | 100% |

Le signe de la corrélation de Spearman est **positif dans 100% des matrices MI testées** pour toutes les valeurs de \(k\).

Ce résultat ne constitue pas une preuve complète de l’équation tensorielle continue, mais il soutient numériquement l’idée centrale : la courbure émergente du graphe est liée aux variations locales de la dimension spectrale.

---

## Structure du dossier

```text
paper/
└── main_fr.tex                          ← Manuscrit complet en français avec section numérique

figures/
├── fig_best_G_vs_T_N20.png               ← Meilleure corrélation courbure / tenseur
└── fig_top10_stability_N20.png           ← Stabilité des proxys

scripts/
├── bup_einstein_tensor_correlation_v1.py
└── bup_einstein_tensor_correlation_v2.py

results/
├── results_einstein_corr_v2/
├── results_einstein_corr_N20_v2/
├── results_einstein_corr_N20_k3_v2/
├── results_einstein_corr_N20_k5_v2/
└── results_einstein_corr_N20_k7_v2/

Travaux liés
Paper 2 — Dimension cosmologique variable et contraintes BAO/SNe : HAL hal-05590614v1
Paper 3 — Croissance linéaire cohérente et σ8=0.772
Paper 4 — Réduction dimensionnelle locale d(z,δ) et galaxies massives JWST
Paper 6 — Formulation pure BuP : δS
BuP
	​

/δW=0

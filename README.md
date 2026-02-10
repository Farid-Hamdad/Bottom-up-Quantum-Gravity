# bottom-up-quantum-gravity

**Émergence de l'espace, du temps et de la gravité à partir de l'intrication quantique**

*Farid Hamdad — Février 2026*

---

## Idée centrale

&gt; **L'espace, le temps et la gravité ne sont pas fondamentaux.**
&gt; Ils émergent collectivement de la structure d'intrication d'un état quantique global fini.

---

## Pourquoi ce projet ?

La physique moderne décrit avec une précision remarquable :

- la physique quantique (théories des champs, information quantique)
- la gravitation classique (relativité générale, thermodynamique des trous noirs)

mais laisse ouverte une question fondamentale :

**Pourquoi l'espace-temps existe-t-il, et pourquoi la gravité possède-t-elle une structure à la fois géométrique et thermodynamique ?**

Ce projet explore une hypothèse minimale mais radicale : *l'espace-temps n'est pas le théâtre de la physique, mais un objet émergent reconstruit à partir de l'intrication quantique.*

---

## Architecture du projet

### 1. Fondements

#### Postulat minimal

Le projet repose sur un postulat unique :

&gt; Il existe un état quantique global pur $|\Psi\rangle \in \mathcal{H} = \bigotimes_{i=1}^N \mathcal{H}_i$,
&gt; défini sur $N$ degrés de liberté élémentaires (qubits),
&gt; **sans espace, sans temps, sans métrique préalable.**

Tout le reste — temps, espace, dimension, géométrie, gravité effective — doit émerger exclusivement de la structure interne de l'intrication de $|\Psi\rangle$.

*Note :* $N$ est fini mais arbitrairement extensible. Les simulations présentées portent sur $N=9$ et $N=16$, avec perspectives pour $N \geq 25$.

---

### 2. Méthodologie

#### Émergence du temps : flot modulaire

Le temps n'est pas un paramètre externe. Il est identifié au **flot modulaire** associé aux matrices de densité réduites :

$$K_A = -\log \rho_A \quad \text{où} \quad \rho_A = \mathrm{Tr}_{\bar{A}}|\Psi\rangle\langle\Psi|$$

Le flot modulaire définit une dynamique intrinsèque :
$$\mathcal{O}(\tau) = e^{iK_A\tau}\mathcal{O}e^{-iK_A\tau}$$

entièrement déterminée par les corrélations entre sous-systèmes. Le temps est interprété comme une notion **relationnelle** (Page-Wootters), issue de la structure informationnelle de l'état.

#### Émergence de l'espace : reconstruction géométrique

| Étape | Outil | Output |
|-------|-------|--------|
| 1. Corrélations | Information mutuelle $I(i:j) = S(\rho_i) + S(\rho_j) - S(\rho_{ij})$ | Matrice de proximité |
| 2. Distance | $d_{ij} = -\log\left(\frac{I(i:j) + \epsilon}{I_{\max}}\right)$ | Métrique informationnelle |
| 3. Géométrie | Multi-Dimensional Scaling (MDS) | Points $x_i \in \mathbb{R}^d$ |
| 4. Dimension | Erreur de reconstruction $\epsilon_d$ | Dimension effective |

*Propriétés de $d_{ij}$ :* symétrie ✓, positivité ✓, triangulaire approximative (violation &lt; 5% sur les données).

---

### 3. Résultats

#### Résultat clé n°1 — La dimension est émergente

| Configuration | Intrication | Dimension effective | Observation |
|-------------|-------------|---------------------|-------------|
| $N=9$, $\lambda \simeq 0$ | Locale (voisins) | $d \simeq 2$ | Reconstruction 2D stable |
| $N=9$, $\lambda \to 1$ | Non-locale (coins) | $d \simeq 3$ | Transition dimensionnelle |
| $N=16$, $\lambda \simeq 0$ | Locale | $d \simeq 2$ | Robustesse confirmée |
| $N=16$, $\lambda \to 1$ | Non-locale | $d \simeq 3$ | Transition plus nette |

**Observation clé :** La dimension spatiale n'est pas postulée, elle est **mesurée** et dépend de la structure de l'intrication (paramètre $\lambda$), pas seulement du nombre de qubits.

#### Résultat clé n°2 — ER = EPR devient mesurable

La correspondance ER = EPR (Maldacena-Susskind) est testée de manière opérationnelle :

| Paramètre | $I(\text{coins})$ | $d_{\text{topo}}$ | $d_{\text{bulk}}$ | Compression $C$ |
|-----------|-------------------|-------------------|-------------------|-----------------|
| $\lambda = 0$ | $\sim 0.01$ | 4 (ou 6 pour $N=16$) | $\sim 4$ | $\sim 1$ |
| $\lambda = 0.5$ | $\sim 0.45$ | 4 | $\sim 1.8$ | $\sim 2.2$ |
| $\lambda = 1.0$ | $\sim 1.28$ | 4 | $\sim 0.04$ | $\sim 100$ ($N=9$) / $\sim 10$ ($N=16$) |

Deux qubits fortement intriqués mais éloignés topologiquement deviennent **géométriquement proches** dans l'espace émergent. L'intrication maximale se manifeste comme un **raccourci géodésique**, signature *wormhole-like* mesurable.

*Note :* Dans ce cadre statique (pas d'évolution temporelle externe), il s'agit d'un **analogue géométrique** de ER=EPR, non d'une traversée dynamique de trou de ver.

#### Résultat clé n°3 — La gravité comme thermodynamique

Une relation de type **Jacobson** (1995) est testée numériquement :

$$\delta S_A \simeq \beta_{\text{eff}} \cdot \delta E$$

| Observable | Comportement | Interprétation |
|------------|--------------|----------------|
| $S_A(\lambda)$ vs $E(\lambda)$ | Relation monotone croissante | Loi d'état effective |
| $\beta_{\text{eff}} = dS/dE$ | Stable à $N=9$, robuste à $N=16$ | Température inverse de l'intrication |
| Test variationnel | $\delta S \approx \beta\,\delta E$ vérifié | Thermodynamique émergente |

La gravité apparaît comme une **loi d'état thermodynamique** de l'information, et non comme une interaction fondamentale.

---

### 4. Discussion : Portée et Limites

#### Ce que cette preuve de principe établit

✅ Un état quantique fini peut encoder une géométrie émergente  
✅ La dimension spatiale dépend quantitativement de la structure d'intrication  
✅ Des signatures analogues à ER=EPR et Jacobson sont observables numériquement  
✅ La robustesse des résultats augmente avec $N$ (9 → 16)

#### Ce qu'elle n'établit pas (limites explicites)

❌ **Limite continue :** $N=9,16$ sont finis ; la convergence $N \to \infty$ reste ouverte  
❌ **Dynamique temporelle :** Le "temps" est modulaire, pas une évolution hamiltonienne externe avec causalité macroscopique  
❌ **Gravité quantique complète :** Pas de dérivation des équations d'Einstein ni de principe holographique quantifié  
❌ **Prédiction expérimentale :** Aucun protocole de test sur système physique (ion piégé, supraconducteur, etc.) proposé

#### Nuance méthodologique

L'hamiltonien $H(\lambda)$ utilisé pour le test Jacobson est un **outil auxiliaire** pour définir une énergie effective. Il n'est pas dérivé fondamentalement de $|\Psi\rangle$ — cette derivation constitue un problème ouvert.

---

## Public visé et prérequis

| Profil | Prérequis suggérés | Ce qu'ils trouveront |
|--------|-------------------|---------------------|
| **M2/Doctorant physique** | MQ basique, entropie de von Neumann | Cadre conceptuel + Implémentation numérique |
| **Chercheur info-quantique** | Intrication, matrices de densité, MDS | Application nouvelle de techniques standards |
| **Physicien théoricien (gravité)** | Jacobson, ER=EPR, holographie | Perspective alternative bottom-up vs top-down (AdS/CFT) |
| **Scientifique généraliste** | Notions de base | Introduction accessible avec annexes explicatives |

*Aucune expertise préalable en gravité quantique n'est strictement requise pour comprendre les concepts principaux.*

---

## Organisation du dépôt
bottom-up-quantum-gravity/
├── paper/
│   ├── main.tex                    # Source LaTeX complet
│   ├── bottom_up_gravity.pdf       # PDF compilé (préprint)
│   └── macros.sty                  # Définitions et packages
├── figures/
│   ├── fig1_dimension_emergente_N9.png
│   ├── fig2_er_epr_N9.png
│   ├── fig3_jacobson_N9.png
│   ├── fig4_grid_N16.png
│   ├── fig5_mds_N16.png
│   ├── fig6_erepr_N16.png
│   └── fig7_jacobson_N16.png
├── code/
│   ├── src/
│   │   ├── states.py               # Génération des états Ψ(λ) 
│   │   ├── geometry.py             # MDS, distances, dimensions spectrales
│   │   ├── thermodynamics.py         # Tests Jacobson, entropies
│   │   └── analysis.py             # Pipelines de figures
│   ├── notebooks/
│   │   └── demo_reproductibilite.ipynb  # Reproduction interactive (N=9)
│   └── requirements.txt            # Python 3.9+, NumPy, SciPy, scikit-learn, matplotlib
├── data/
│   ├── N9_results.pkl              # Données numériques brutes (reproductibilité)
│   └── N16_results.pkl
└── README.md                       # Ce document

---

## Perspectives et extensions

### Court terme (en cours)
- [ ] **Monte Carlo tensoriel** pour $N \geq 25$ (PEPS, MERA) — test de convergence
- [ ] **Finite-size scaling** : la transition dimensionnelle devient-elle abrupte en limite thermodynamique ?
- [ ] **Autres familles d'états** : états aléatoires, états topologiques (toric code), états critiques

### Moyen terme (conceptuel)
- [ ] Connexion avec **Causal Set Theory** (rideaux de lumière discrets)
- [ ] Généralisation à **qudits** ($d>2$) et hamiltoniens non-triviaux
- [ ] Test de la **conjecture de Page-Wootters** pour le temps relationnel

### Long terme (spéculatif)
- [ ] **Principe holographique émergent** : volume $\leftrightarrow$ surface dans ce cadre ?
- [ ] Liens avec **modèles SYK** ou émulateurs de trous noirs en laboratoire
- [ ] **Prédiction expérimentale** : signatures mesurables sur simulateurs quantiques ?

---

## Citation

Si vous utilisez ou prolongez ce travail :

```bibtex
@misc{hamdad2026bottomup,
  author = {Hamdad, Farid},
  title = {Bottom-Up Quantum Gravity: \'Emergence de l'espace, du temps et 
           de la gravité à partir de l'intrication quantique},
  year = {2026},
  howpublished = {GitHub repository},
  url = {https://github.com/Farid-Hamdad/Bottom-Up-Quantum-Gravity}
}

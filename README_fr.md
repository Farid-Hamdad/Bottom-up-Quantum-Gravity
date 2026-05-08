# Bottom-Up Quantum Gravity (BuP)

**Émergence de l'espace, du temps et de la gravité depuis l'intrication quantique**  
Farid Hamdad — 2026

## Publication

F. Hamdad, "Contraintes cosmologiques sur une dimension spatiale émergente :
dimension variable, énergie noire et σ8 dans le cadre Bottom-Up Quantum Gravity",
Prépublication HAL hal-05590614v1 (2026).
https://hal.science/hal-05590614v1

---

## Vue d'ensemble

**Bottom-Up Quantum Gravity (BuP)** étudie une hypothèse minimale :

> L'espace, le temps et la gravité ne sont pas fondamentaux.  
> Ils émergent collectivement de la structure d'intrication d'un état quantique global fini.

Dans ce cadre :

- l'intrication définit la connectivité,
- la géométrie émerge de la structure informationnelle,
- le temps apparaît via le flot modulaire,
- la gravité est interprétée comme la courbure effective de cette structure.

Le projet combine une base conceptuelle en information quantique,
des simulations numériques exactes sur des systèmes finis, et une
extension phénoménologique aux courbes de rotation galactiques (données SPARC).

---

## Idée centrale

À partir d'un état quantique pur global \(|\Psi\rangle\) sans espace préexistant,
une distance informationnelle est définie :

$$
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max} + \epsilon}\right)
$$

où \(I(i:j)\) est l'information mutuelle entre les sous-systèmes \(i\) et \(j\).

À partir de cette distance, le cadre reconstruit :

- un graphe d'intrication et son plongement,
- un Laplacien discret,
- une dimension spectrale \(d_s(\tau)\),
- un profil radial effectif \(d(r)\).

---

## Résumé des résultats

| Résultat | Statut |
|----------|--------|
| Courbes de rotation émergentes (SPARC, 175 galaxies) | forte preuve numérique |
| Dimension émergente et régimes | robuste |
| Signal de courbure d'Ollivier–Ricci | statistiquement robuste |
| Diagnostics de chaos modulaire | cohérents avec un régime chaotique |
| Prédiction de lentillage | préliminaire |
| Symétrie effective du vide (type U(1)) | observée numériquement |
| Hamiltoniens modulaires comme sondes géométriques | établi numériquement |
| Unicité du vide et alignement d'axe | observé |
| Dimension émergente microscopique stable \(d_s \approx 2,61 \pm 0,03\) | robuste sur \(N = 6\)–\(16\) |

---

## 1. Temps émergent : flot modulaire

Pour un sous-système \(A\) de matrice densité réduite
\(\rho_A = \mathrm{Tr}_{\bar{A}} |\Psi\rangle\langle\Psi|\),
l'hamiltonien modulaire est :

$$
K_A = -\log \rho_A
$$

Le flot modulaire

$$
O(\tau) = e^{iK_A\tau}\, O\, e^{-iK_A\tau}
$$

définit une dynamique relationnelle intrinsèque. Le temps n'est pas un paramètre
d'arrière-plan mais une propriété informationnelle du système.

---

## 2. Diagnostics de chaos modulaire

Le spectre de \(K_A\) est analysé via :

- les statistiques du rapport des gaps adjacents,
- les diagnostics de la théorie des matrices aléatoires,
- le facteur de forme spectral.

$$
\langle r \rangle = \left\langle \frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})} \right\rangle
$$

| Régime | Valeur |
|--------|-------:|
| Poisson | 0,386 |
| GOE | 0,536 |
| GUE | 0,603 |

Les résultats numériques donnent \(\langle r \rangle \in [0,53 ; 0,59]\),
ce qui est cohérent avec un régime chaotique.

---

## 3. Espace émergent

La géométrie est reconstruite à partir de l'information mutuelle :

$$
I(i:j) = S(\rho_i) + S(\rho_j) - S(\rho_{ij})
$$

$$
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max}} + \varepsilon\right)
$$

La géométrie n'est pas postulée — elle est reconstruite à partir de la
structure d'intrication de \(|\Psi\rangle\).

---

## 4. Courbure discrète d'Ollivier–Ricci

Le signal de courbure discrète est le résultat numérique le plus
robuste du projet :

$$
\Delta\kappa_{\mathrm{edge}} \approx 0,0758
$$

- fraction positive : 1,00 sur toutes les graines
- \(p \sim 10^{-6}\)
- effet proche-lointain : \(\Delta\kappa_{\mathrm{near-far}} > 0\)

Le déficit informationnel produit un signal de courbure détectable
sans aucun postulat géométrique.

---

## 5. Horizon émergent (\(N = 16\))

Une région de type horizon est détectée, caractérisée par :

- une entropie d'intrication élevée,
- une faible conductance,
- un fort couplage interne.

Cette structure émerge sans géométrie pré-définie.

---

## 6. Gravité thermodynamique

La relation

$$
\delta S \approx \beta_{\mathrm{eff}}\,\delta E
$$

est numériquement stable sur l'ensemble des configurations, ce qui
est cohérent avec une interprétation thermodynamique émergente de la gravité.

---

## 7. Extension SPARC (175 galaxies)

Le modèle de dimension effective BuP est testé sur les courbes de rotation
SPARC.

| Régime | \(d_{\min}\) |
|--------|------------:|
| Galaxies massives (HIGH) | 2,49 |
| Galaxies naines (LOW) | 2,27 |

La structure en régimes est émergente, non imposée.

**Lien avec la structure microscopique.**

La dimension dépendante de l'échelle \(d(r)\) observée dans SPARC
peut être interprétée comme la manifestation à grande échelle
des structures d'interaction dépendantes de l'état observées dans
les hamiltoniens modulaires.

Cela suggère une chaîne continue :

$$
\text{intrication quantique}
\;\rightarrow\; \text{structure modulaire}
\;\rightarrow\; \text{graphe d'interaction}
\;\rightarrow\; \text{dimension effective } d(r).
$$

L'accord avec les données SPARC est donc cohérent avec les mécanismes
microscopiques identifiés dans ce travail.

---

## 8. Prédiction de lentillage (102 galaxies)

| Régime | \(x_{\mathrm{peak}}\) |
|--------|----------------------:|
| LOW | 0,365 |
| HIGH | 0,498 |

- Régime LOW : signal de lentillage centralisé
- Régime HIGH : signal de lentillage étendu

Cette prédiction constitue un discriminant potentiel entre
BuP et les modèles standards de matière noire.

---

## 9. Symétrie effective du vide

La structure de symétrie effective du vide BuP est étudiée
en utilisant un état fondamental d'un hamiltonien anisotrope
avec de faibles champs aléatoires locaux :

$$
H_{\mathrm{BuP}} =
-J_{zz}\sum_i Z_i Z_{i+1}
-J_{xy}\sum_i (X_iX_{i+1}+Y_iY_{i+1})
+\sum_i (h_i^x X_i + h_i^y Y_i + h_i^z Z_i)
$$

**Méthode.** À partir de l'état fondamental, des matrices de rotation effectives
\(R_{ij} \in SO(3)\) sont extraites par décomposition SVD polaire du tenseur
de corrélation de Pauli. L'algèbre de Lie engendrée par \(\{\log R_{ij}\}\) est
analysée pour classer le groupe de symétrie effectif.

**Résultat principal** (\(N=8\), \(J_{zz}=2,0\), \(J_{xy}=0,5\), \(h=0,1\)) :

- 95% des réalisations du désordre classées comme de type U(1)
- paramètre d'ordre moyen \(OP \approx 0,969\) contre \(OP_{\mathrm{aléatoire}} \approx 0,619\)

**Scan sur le désordre** \(h_{\mathrm{échelle}}\) :
régime U(1) robuste pour \(h \lesssim 0,2\),
transition près de \(h \sim 0,3\),
régime SO(3) retrouvé pour \(h \gtrsim 0,5\).

**Scan sur l'anisotropie** \(J_{zz}/J_{xy}\) :
transition du type SO(3) vers U(1) près de \(J_{zz}/J_{xy} \approx 1\),
régime U(1) robuste pour \(J_{zz}/J_{xy} \geq 2\).

**Interprétation.** Il ne s'agit pas d'une brisure spontanée de symétrie
au sens de la théorie des champs — les champs locaux brisent déjà
explicitement SO(3). Le résultat démontre une sélection robuste d'un axe effectif :
le vide BuP est un milieu structuré, non un état neutre isotrope.

**Stabilité kNN** (\(k = 1, \ldots, 7\)) : le signal de type U(1)
est stable pour toutes les valeurs de \(k\), confirmant qu'il s'agit
d'une propriété intrinsèque de l'état fondamental et non d'un artefact
de la reconstruction du graphe.

**Unicité du vide** : la dispersion angulaire entre les axes du vide
sur les réalisations du désordre est \(\Delta_{\mathrm{vide}} = 0,056\) rad
(normalisée : \(0,072\), où \(1\) correspond à une \(S^2\) isotrope).
L'axe du vide s'aligne sur la direction d'anisotropie \(J_{zz}\)
avec \(\langle|\cos z|\rangle = 0,999\).

---

## 10. Hamiltoniens modulaires comme sondes géométriques

**Article :** `bup_modular_paper_v3.tex`

La structure des hamiltoniens modulaires \(K_A\) est étudiée
systématiquement sur 180 configurations numériques
(\(N \in \{8,10,12\}\), \(|A| \in \{3,4,5,6\}\),
\(\eta_{\mathrm{micro}} \in \{1,0 ; 1,1 ; 1,25 ; 1,5 ; 2,0\}\),
trois classes d'états).

**Chaîne de dérivation :**

$$
|\Psi\rangle
\;\xrightarrow{\mathrm{Tr}_{\bar{A}}}\;
\rho_A
\;\xrightarrow{-\log}\;
K_A
\;\xrightarrow{\mathrm{ajustement}}\;
H_{\mathrm{eff}}^{1+2}
$$

**Résultat 1 — Quasi-localité :**

| Classe d'état | \(\langle R^2 \rangle\) | \(\sigma_{R^2}\) |
|---------------|------------------------:|-----------------:|
| Fondamental | 0,9948 | 0,0099 |
| Thermique | 0,9994 | 0,0011 |
| Aléatoire | 0,9205 | 0,0741 |

**Résultat 2 — Amplification anisotropique dépendante de l'état :**

$$
\mathcal{A} = \frac{\eta_{\mathrm{eff}}}{\eta_{\mathrm{micro}}}
$$

| Classe d'état | \(\langle\mathcal{A}\rangle\) | \(\sigma_{\mathcal{A}}\) |
|---------------|------------------------------:|-------------------------:|
| Fondamental | 1,048 | 0,114 |
| Thermique | 1,107 | 0,125 |
| Aléatoire | 0,805 | 0,392 |

Les états structurés (fondamental, thermique) amplifient l'anisotropie
microscopique ; les états aléatoires la suppriment. La hiérarchie

$$
\mathcal{A}_{\mathrm{thermique}} \gtrsim \mathcal{A}_{\mathrm{fondamental}} > 1 > \mathcal{A}_{\mathrm{aléatoire}}
$$

dépend de l'état : les états aléatoires construits à partir du même
hamiltonien montrent une suppression, confirmant que \(K_A\) agit comme
un filtre sensible à la cohérence.

**Résultat 3 — Amplification monotone avec \(\eta_{\mathrm{micro}}\) :**

Pour les états fondamentaux, \(\mathcal{A}\) augmente monotone de
\(1,000\) à \(\eta_{\mathrm{micro}} = 1,0\) jusqu'à \(1,144 \pm 0,190\)
à \(\eta_{\mathrm{micro}} = 2,0\).

**Implication pour BuP.** Cela établit un pont numérique concret
entre la structure d'intrication et la dynamique effective émergente,
remplaçant des hamiltoniens postulés par des opérateurs dérivés de \(|\Psi\rangle\).

**Lien avec le chaos modulaire.**

Les propriétés spectrales de \(K_A\) sont cohérentes avec un comportement
chaotique, avec \(\langle r \rangle \in [0,53 ; 0,59]\), correspondant aux
prédictions GOE/GUE.

Cela établit que les hamiltoniens modulaires sont non seulement
des opérateurs quasi-locaux, mais aussi qu'ils appartiennent à une
classe chaotique universelle, renforçant leur interprétation comme
générateurs dynamiques effectifs.

**Interprétation géométrique.**

Les couplages extraits \(J_{ij}\) définissent un graphe d'interaction
qui peut être interprété comme une géométrie émergente.

La localité de \(K_A\) correspond à une structure à courte portée,
tandis que l'anisotropie encode une déformation directionnelle.

Cela fournit une réalisation concrète du principe BuP :

$$
\text{intrication} \;\rightarrow\; \text{géométrie} \;\rightarrow\; \text{dynamique effective.}
$$

---

## 11. Structure du dépôt
experiments/
10_single_galaxy_ngc3198/ # Ajustement d'une galaxie (NGC 3198)
20_v1_3_scan_ngc3198/ # Scan de paramètres
30_multi_galaxy_regime_3_08_4_5/ # Analyse multi-galaxies par régime
40_micro_dimension_coherence/ # Dimension émergente stable sur N=6–16

derive_d_r/ # Pipeline de dimension effective
results_bup_hybrid_multi/ # Résultats multi-galaxies
results_bup_hybrid_multi_rmax_fixed/ # Variante à r_max fixé

jauge/ # Structure de jauge et symétrie du vide
bup_gauge_classification.py # Classification des groupes SO(3)
bup_symmetry_scan.py # Scan comparatif des mécanismes
bup_vacuum_v3_standalone.py # Symétrie effective du vide (v3)
bup_knn_scan.py # Scan de stabilité kNN
bup_vacuum_uniqueness.py # Analyse d'unicité du vide
bup_spontaneous_breaking.py # Test brisure explicite vs spontanée

modular/ # Étude des hamiltoniens modulaires
bup_modular_paper_v3.tex # Article (prêt à soumettre)
aggregate_rows.csv # Jeu de données numériques complet (180 configs)
fig_amplification_by_state.png
fig_amplification_vs_eta.png
fig_fit_r2_by_state.png
fig_locality_vs_eta.png

text

---

## 12. Interprétation physique

La matière modifie la structure informationnelle de \(|\Psi\rangle\),
ce qui modifie la géométrie émergente et produit un effet gravitationnel effectif.

Aucune particule de matière noire n'est postulée.
L'excès gravitationnel peut être interprété comme une conséquence
géométrique de la structure d'intrication.

---

## 13. Limitations

- Modèle effectif, pas une théorie relativiste complète
- Lentillage via observable proxy, pas une déflexion relativiste complète
- Dépendance au rapport masse/luminosité dans les ajustements SPARC
- États fondamentaux d'hamiltoniens XXZ utilisés comme proxys BuP —
  les états générés par la dynamique variationnelle BuP restent à étudier

---

## Image unifiée

Le cadre BuP peut être résumé en une chaîne unique :

- L'intrication définit la connectivité
- La connectivité définit la géométrie
- La géométrie définit la dynamique effective
- La dynamique effective reproduit la phénoménologie gravitationnelle

Ce travail fournit des preuves numériques à chaque étape de cette chaîne,
des systèmes quantiques ($N \leq 16$ qubits) aux échelles astrophysiques
(175 galaxies SPARC).

---

## 14. Perspectives

- Calcul complet de la lentille relativiste
- Confrontation avec les relevés de lentillage faible
- Extension cosmologique
- Dynamique BuP fondamentale à partir d'un principe variationnel
- Lien entre la structure des hamiltoniens modulaires et les champs de jauge

---

## 15. Dimension émergente stable aux échelles microscopique, galactique et cosmologique

Un scan multi-tailles dédié sur des états BuP structurés pour
\(N = 6, 8, 9, 12, 16\) montre que la dimension microscopique effective
extraite de l'estimateur de portée d'intrication reste remarquablement stable :

$$
d_s(\lambda, N) \approx 2,61 \pm 0,03.
$$

Deux faits robustes émergent des données brutes :

- **stabilité en couplage non-local** : \(d_s\) varie très faiblement avec
  le paramètre non-local \(\lambda \in [0,1]\) ;
- **stabilité en taille du système** : le même plateau persiste de
  \(N=6\) à \(N=16\), sans indication de dérive forte due à la taille
  finie dans les mesures brutes.

Ce plateau brut est directement cohérent avec l'inférence galactique SPARC
\(d_{\mathrm{eff}} = 2,644 \pm 0,295\), et reste compatible avec la
fenêtre cosmologique tardive préférée par DESI + Pantheon+.
Le recouvrement s'étend donc des simulations à l'échelle des qubits
jusqu'à la phénoménologie galactique et cosmologique sur de nombreux
ordres de grandeur.

L'expérience, les scripts, le scan brut, la sortie FSS exploratoire,
le tableau récapitulatif et la figure correspondants sont stockés dans :

- `experiments/40_micro_dimension_coherence/`

**Note importante.**
L'extrapolation par analyse de taille finie (FSS) stockée dans le fichier
exploratoire `fss_extrapolated.csv` n'est pas utilisée comme conclusion
principale. Dans le régime actuel des données, les mesures brutes sont
substantiellement plus stables et physiquement plus significatives que
l'extrapolation ajustée vers \(N \to \infty\).

---

## Auteur

**Farid Hamdad**  
Bottom-Up Quantum Gravity, 2026  
hamdadfarid54@gmail.com  
[github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity](https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity)
---


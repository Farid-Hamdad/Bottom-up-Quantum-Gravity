# Paper 15 — Des graphes d'intrication aux équations d'Einstein

**Dérivation effective de la gravité continue dans Bottom-Up Quantum Gravity**

---

## Statut

Paper 15 étudie si l'équation variationnelle discrète de BuP sur le réseau d'intrication,

$$
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0,
$$

admet une limite continue lisse de type Einstein :

$$
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}\,g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}\,T_{\mu\nu}^{\rm ent}.
$$

L'objectif n'est pas de postuler la relativité générale, mais de la retrouver comme limite continue stable d'une condition d'équilibre sur le graphe d'intrication.

---

## 1. Point de départ

L'objet fondamental est la matrice d'information mutuelle :

$$
W_{ij}=I(i:j).
$$

À partir d'elle, on construit le laplacien d'intrication :

$$
L_{\rm ent}=D-W,
\qquad
D_{ii}=\sum_j W_{ij}.
$$

Paper 15 part de l'action BuP minimale :

$$
S_{\rm BuP}[W]
=
\mathrm{Tr}\,L(W)^{-\beta[W]}
+
\sum_{ij}W_{ij}\,\kappa_{ij}[W]
+
\lambda\sum_{ij}W_{ij}\,d_{ij}^{2}
+
S_{\rm topo}[W].
$$

Les quatre termes sont :

| Terme | Signification | Rôle |
|-------|---------------|------|
| $\mathrm{Tr}\,L^{-\beta}$ | action spectrale | géométrie globale |
| $\sum W_{ij}\kappa_{ij}$ | courbure discrète | réponse de Ricci locale |
| $\lambda\sum W_{ij}d_{ij}^2$ | coût de localité | supprime les graphes complets non locaux |
| $S_{\rm topo}[W]$ | topologie | contraintes de connectivité globale |

Paper 14 a fixé l'exposant spectral :

$$
\beta
=
F(\lambda_2,\; d_s,\; d_w,\; \Delta_3),
$$

de sorte que le terme spectral n'est plus libre.

---

## 2. Objectif central

Paper 15 teste la chaîne :

$$
W_{ij}
\;\longrightarrow\;
L_\epsilon
\;\longrightarrow\;
\Delta_g,
$$

$$
\kappa_{ij}^{\rm OR}
\;\longrightarrow\;
R_{\mu\nu}\,u^\mu u^\nu,
$$

$$
\delta W_{\rm loc}
\;\longrightarrow\;
\delta\kappa(r)
\;\longrightarrow\;
T_{\mu\nu}^{\rm eff}.
$$

Si ces limites sont valides, l'action continue devrait prendre la forme schématique :

$$
S_{\rm cont}[g]
=
\int d^Dx\,\sqrt{|g|}
\left[
\frac{1}{16\pi G_{\rm eff}}\,R
+
\Lambda_{\rm ent}
+
\mathcal{L}_{\rm ent-source}
+
\mathcal{H}
\right],
$$

où $\mathcal{H}$ contient des corrections de courbure supérieure ou non locales.

Dans la limite lisse de basse énergie,

$$
\mathcal{H}_{\mu\nu}\to 0,
$$

et l'on s'attend à :

$$
G_{\mu\nu}
+
\Lambda_{\rm ent}\,g_{\mu\nu}
=
8\pi G_{\rm eff}\,T_{\mu\nu}^{\rm ent}.
$$

---

## 3. État numérique — trois premiers tests positifs

Paper 15 possède actuellement trois piliers numériques positifs.

---

### Étape B — Convergence spectrale

La première tâche est de tester :

$$
L_N\;\longrightarrow\;\Delta_g.
$$

Une première tentative utilisant le laplacien normalisé non mis à l'échelle,

$$
L_{\rm norm}=I-D^{-1/2}WD^{-1/2},
$$

a échoué sur les géométries 2D, la grille et la sphère s'effondrant toutes deux vers :

$$
d_s\simeq 1.2.
$$

La version corrigée utilise le laplacien remis à l'échelle continue :

$$
L_\epsilon=\frac{D-W}{\epsilon}.
$$

Avec $k=\sqrt{N}$, $\epsilon=0.5\,\epsilon_{\rm knn}$, et $N=1024$, les dimensions spectrales mesurées sont :

| Géométrie | Dimension cible $d_s$ | $d_s$ mesurée | Erreur |
|-----------|:---------------------:|:-------------:|:------:|
| cercle    | 1                     | 1.0021        | 0.0021 |
| intervalle | 1                    | 0.9445        | 0.0555 |
| grille 2D | 2                     | 1.9938        | 0.0062 |
| sphère    | 2                     | 1.8964        | 0.1036 |

Cela fournit la première preuve numérique positive pour :

$$
L_N\;\longrightarrow\;\Delta_g.
$$

---

### Étape C — Signal de courbure d'Ollivier–Ricci

La deuxième tâche est de tester si la courbure d'Ollivier–Ricci discrète détecte le signal de Ricci continu :

$$
\kappa_{ij}^{\rm OR}
\;\longrightarrow\;
R_{\mu\nu}\,u^\mu u^\nu.
$$

Une première version utilisant une grille plane avec bord était contaminée par des effets de bord. La version corrigée utilise un tore plat périodique comme référence de courbure nulle et des distances géodésiques pour la sphère.

À la résolution finale :

$$
\bar\kappa_{\rm flat}=-0.001671\;\simeq\;0,
$$

tandis que :

$$
\bar\kappa_{\rm sphere}=0.008854.
$$

L'excès de courbure relatif est :

$$
\Delta\bar\kappa_{\rm sphere-flat}=0.010525.
$$

Cela donne un premier signal numérique relatif pour :

$$
\kappa_{ij}^{\rm OR}
\;\to\;
R_{\mu\nu}\,u^\mu u^\nu.
$$

Le résultat reste *qualitatif* : il distingue la géométrie plate périodique de la courbure positive, mais ne prouve pas encore la convergence ponctuelle vers le tenseur de Ricci.

---

### Étape E — Test source‑réponse

La troisième tâche est de tester :

$$
\delta W_{\rm loc}
\;\longrightarrow\;
\delta\kappa(r).
$$

Le test source‑réponse v2 utilise une perturbation radiale lisse de l'intrication :

$$
\phi_i
=
\exp\left[
-\frac{d(i,\mathrm{source})^2}{2\sigma^2}
\right],
$$

et modifie les poids comme suit :

$$
W'_{ij}
=
W_{ij}
\left[
1+s\,\frac{\phi_i+\phi_j}{2}
\right].
$$

La réponse de courbure est :

$$
\Delta\kappa_{ij}
=
\kappa_{ij}^{\rm after}
-
\kappa_{ij}^{\rm before}.
$$

Pour $s=-0.30$, la réponse est fortement localisée :

| Géométrie     | rapport near/far | Spearman $(\phi_{\rm edge},|\Delta\kappa|)$ | p‑valeur |
|---------------|:----------------:|:-------------------------------------------:|:--------:|
| tore plat 2D  | 7.75             | 0.583                                       | $4.12\times10^{-83}$ |
| sphère        | 11.51            | 0.752                                       | $5.49\times10^{-165}$ |

Ainsi :

$$
\delta W_{\rm loc}
\;\longrightarrow\;
\delta\kappa(r)
$$

est numériquement soutenu sur des géométries contrôlées.

---

## 4. Lien avec Paper 7 et Paper 8

Le résultat de l'Étape E n'est pas isolé. Il confirme le même couplage déjà observé dans Paper 7 et Paper 8 :

$$
\delta W_{\rm loc}
\;\longrightarrow\;
\delta\kappa.
$$

| Paper                | Test                                | Cadre                         | Signal |
|----------------------|-------------------------------------|-------------------------------|--------|
| Paper 7              | réponse directe de courbure         | graphes MI quantiques, $N=16$ | $\Delta\kappa_{\rm edge}=0.076$, fraction positive $=100\%$ |
| Paper 8              | reconstruction de $T_{\mu\nu}^{\rm eff}$ | source depuis $\delta W_{\rm loc}$ | Spearman $\rho=0.741$ |
| Paper 15, Étape E    | réponse radiale source              | géométries contrôlées         | Spearman $\rho=0.752$, rapport near/far $=11.51$ sur sphère |

Le même signal apparaît sous trois angles indépendants :

$$
\text{réponse locale de courbure}
\;\leftrightarrow\;
\text{reconstruction du tenseur énergie‑impulsion effectif}
\;\leftrightarrow\;
\text{géométrie source‑réponse contrôlée}.
$$

Paper 15 est plus fort que les tests précédents dans un sens : la géométrie est contrôlée, et les cas de référence plats/courbés sont connus.

---

## 5. Relation avec la construction de source de Paper 8

Paper 8 a construit un tenseur énergie‑impulsion effectif à partir de perturbations locales d'intrication :

$$
T_{\mu\nu}^{\rm eff}
=
T_{\mu\nu}^{\rm matter}[\delta W^{\rm loc}]
+
T_{\mu\nu}^{\rm ent}[d_s].
$$

La meilleure source candidate était :

$$
S_{\rm flux}
=
T_{00}
-
\frac{1}{2}T_{aa}
+
\frac{1}{2}T_{\rm grad}
+
F.
$$

Paper 8 a trouvé :

$$
\rho_{\rm Spearman}=0.741,
\qquad
p=1.84\times10^{-4},
$$

entre la source reconstruite et les fluctuations de courbure.

Il a également montré que le secteur de type matière n'est pas conservé isolément :

$$
\nabla^\mu T_{\mu\nu}^{\rm matter}\neq 0.
$$

Dans BuP, cela est interprété comme un échange avec le fond d'intrication :

$$
J_\nu^{\rm exchange}
=
\nabla^\mu
\left[
G_{\rm eff}(d_s)\,T_{\mu\nu}^{\rm matter}
\right],
$$

avec compensation :

$$
\nabla^\mu
\left[
G\,T_{\mu\nu}^{\rm ent}
\right]
=
-
J_\nu^{\rm exchange}.
$$

Ainsi, l'Étape E de Paper 15 n'est pas une hypothèse ouverte : c'est la version sur géométries contrôlées du mécanisme de source de Paper 8.

---

## 6. Interprétation actuelle

Le résultat actuel peut être résumé par :

$$
\boxed{
\text{Paper 15 ne prouve pas encore Einstein, mais valide les trois flèches nécessaires à la dérivation effective.}
}
$$

Les trois flèches validées sont :

$$
W_{ij}\;\to\; L_\epsilon\;\to\;\Delta_g,
$$

$$
\kappa_{ij}^{\rm OR}\;\to\;\text{signal de Ricci},
$$

$$
\delta W_{\rm loc}\;\to\;\delta\kappa(r).
$$

Ensemble, elles soutiennent la cible continue :

$$
\frac{\delta S_{\rm BuP}}{\delta W_{ij}}=0
\quad
\xrightarrow[N\to\infty]{}
\quad
G_{\mu\nu}
+
\Lambda_{\rm ent}\,g_{\mu\nu}
=
8\pi G_{\rm eff}\,T_{\mu\nu}^{\rm ent}.
$$

---

## 7. Structure des dossiers

```text
papers/paper15_einstein_derivation/
  README.md
  paper15_einstein_derivation.tex

  scripts/
    paper15_spectral_convergence_v1.py
    paper15_spectral_convergence_v2.py
    paper15_ricci_convergence_v1.py
    paper15_ricci_convergence_v2.py
    paper15_source_response_v1.py
    paper15_source_response_v2.py
    paper15_build_numerical_summary_v1.py

  results/
    spectral_convergence_v1/
    spectral_convergence_v2/
    spectral_convergence_v2_best_unnormalized/
    ricci_convergence_v1/
    ricci_convergence_v2/
    source_response_v1/
    source_response_v2/
    paper15_numerical_summary_v1/

  figures/
    # figures finales copiées depuis les dossiers de résultats sélectionnés
      fig01_spectral_ds_grid2d.png
      fig02_spectral_ds_sphere.png
      fig03_spectral_error_vs_N_grid2d.png
      fig04_spectral_error_vs_N_sphere.png
      fig05_ricci_mean_kappa_vs_N.png
      fig06_ricci_final_mean_kappa_by_geometry.png
      fig07_source_response_correlation_sphere.png
      fig08_source_response_signed_sphere.png

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

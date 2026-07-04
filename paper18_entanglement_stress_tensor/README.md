# Paper 18 — Tenseur énergie-impulsion émergent à partir de l'intrication

**Des perturbations locales d'intrication à une source effective de tenseur énergie-impulsion**

---

## Statut

Paper 18 étudie la troisième flèche nécessaire à la limite einsteinienne effective de BuP :

$$
\delta W_{\rm loc}
\longrightarrow
T_{\mu\nu}^{\rm ent}.
$$

Paper 16 a étudié le côté géométrie spectrale :

$$
L_N\to\Delta_g.
$$

Paper 17 a étudié le côté courbure :

$$
\kappa_{ij}^{OR}
\to
R_{\mu\nu}u^\mu u^\nu.
$$

Paper 18 étudie le côté source : comment une perturbation locale du graphe d'intrication devient une source effective de tenseur énergie-impulsion.

---

## 1. Point de départ

L'objet fondamental est le graphe d'intrication :

$$
W_{ij}=I(i:j).
$$

Une source locale est modélisée comme une perturbation :

$$
W_{ij}
\to
W'_{ij}
=
W_{ij}
+
\delta W_{ij}^{\rm loc}.
$$

La question centrale est :

$$
\boxed{
\text{Peut-on interpréter } \delta W_{\rm loc} \text{ comme une source effective } T_{\mu\nu}^{\rm ent} \text{ ?}
}
$$

---

## 2. Preuves antérieures

### Paper 8

Paper 8 a reconstruit une grandeur proxy de source effective à partir de perturbations locales d'intrication :

$$
T_{\mu\nu}^{\rm eff}
=
T_{\mu\nu}^{\rm matter}[\delta W_{\rm loc}]
+
T_{\mu\nu}^{\rm ent}[d_s].
$$

Il a trouvé :

$$
\rho_{\rm Spearman}=0.741,
\qquad
p=1.84\times10^{-4}.
$$

Cela a soutenu :

$$
\delta W_{\rm loc}
\to
T_{\mu\nu}^{\rm eff}.
$$

### Paper 15 Étape E

Paper 15 a montré qu'un défaut radial d'intrication produit une réponse de courbure localisée :

$$
\delta W_{\rm loc}
\to
\delta\kappa(r).
$$

Sur la sphère :

$$
{\rm Spearman}(\phi_{\rm edge},|\Delta\kappa|)
=
0.752,
$$

avec un rapport near/far de :

$$
11.51.
$$

Ainsi Paper 18 part de deux signaux antérieurs positifs.

---

## 3. Identité centrale : première loi modulaire

L'entrée théorique clé est la première loi modulaire :

$$
\delta S_A
=
\delta\langle K_A\rangle.
$$

avec

$$
K_A=-\log\rho_A.
$$

Dans la TQFT continue, les hamiltoniens modulaires locaux relient les variations d'entropie aux variations du tenseur énergie-impulsion :

$$
\delta\langle K_A\rangle
\sim
\int_A \xi^\mu \delta T_{\mu\nu}d\Sigma^\nu.
$$

Dans BuP, l'analogue sur graphe est :

$$
\delta S_A[W]
\simeq
\delta\langle K_A[W]\rangle
\to
T_{\mu\nu}^{\rm ent}.
$$

---

## 4. Chaîne cible

Paper 18 vise à établir :

$$
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r)
\to
T_{\mu\nu}^{\rm ent}.
$$

---

## 5. Principaux résultats numériques

Paper 18 possède actuellement deux résultats positifs.

| Étape | Géométrie | Quantité | Valeur | Cible | Statut |
|---|---|---|---:|---|---|
| v1 première loi modulaire du graphe | tore plat 2D | pente $\delta S$ vs $\delta\langle K\rangle$ | 0.989481 | 1 | positif |
| v1 première loi modulaire du graphe | tore plat 2D | $R^2$ $\delta S$ vs $\delta\langle K\rangle$ | 0.996590 | proche de 1 | positif |
| v1 première loi modulaire du graphe | sphère | pente $\delta S$ vs $\delta\langle K\rangle$ | 0.989718 | 1 | positif |
| v1 première loi modulaire du graphe | sphère | $R^2$ $\delta S$ vs $\delta\langle K\rangle$ | 0.996592 | proche de 1 | positif |
| v2 source modulaire courbure | tore plat 2D | rapport de localisation near/far | 6.293007 | $>1$ | positif |
| v2 source modulaire courbure | tore plat 2D | $R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})$ | 0.967229 | proche de 1 | positif |
| v2 source modulaire courbure | tore plat 2D | Pearson signé $\delta K\to\Delta\kappa_{\rm near}$ | -0.999168 | $|r|$ proche de 1 | positif |
| v2 source modulaire courbure | sphère | rapport de localisation near/far | 9.533529 | $>1$ | positif |
| v2 source modulaire courbure | sphère | $R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})$ | 0.982252 | proche de 1 | positif |
| v2 source modulaire courbure | sphère | Pearson signé $\delta K\to\Delta\kappa_{\rm near}$ | -0.999428 | $|r|$ proche de 1 | positif |

---

## 6. Résultat v1 — Première loi modulaire du graphe

Une grandeur proxy de densité du graphe est définie sur une région $A$ par

$$
\rho_A
=
\frac{(L_A+\mu I)^{-1}}
{\mathrm{Tr}(L_A+\mu I)^{-1}}.
$$

Alors

$$
S_A=-\mathrm{Tr}(\rho_A\log\rho_A),
\qquad
K_A=-\log\rho_A.
$$

Après une perturbation radiale locale de $W_{ij}$, le test mesure

$$
\delta S_A=S_A(W')-S_A(W),
$$

et

$$
\delta\langle K_A\rangle
=
\mathrm{Tr}\left[(\rho'_A-\rho_A)K_A\right].
$$

Le résultat est fortement positif.

Sur le tore plat :

$$
\delta S_A
=
-0.000068
+
0.989481\,\delta\langle K_A\rangle,
$$

avec

$$
R^2=0.996590,
\qquad
{\rm Pearson}=0.998293.
$$

Sur la sphère :

$$
\delta S_A
=
-0.000081
+
0.989718\,\delta\langle K_A\rangle,
$$

avec

$$
R^2=0.996592,
\qquad
{\rm Pearson}=0.998295.
$$

Ainsi :

$$
\boxed{
\delta S_A^{\rm graphe}
\simeq
\delta\langle K_A^{\rm graphe}\rangle.
}
$$

---

## 7. Résultat v2 — La source modulaire prédit la réponse de courbure

Le second test relie la réponse modulaire à la courbure :

$$
\delta W_{\rm loc}
\to
\delta\langle K_A\rangle
\to
\delta\kappa(r).
$$

La réponse de courbure est localisée autour de la source :

| Géométrie | rapport near/far moyen | rapport near/far médian |
|---|---:|---:|
| tore plat 2D | 6.293 | 6.241 |
| sphère | 9.534 | 9.479 |

Plus important encore, l'amplitude de la source modulaire prédit la réponse de courbure près de la source :

| Géométrie | $R^2(|\delta\langle K_A\rangle|,\langle|\Delta\kappa|\rangle_{\rm near})$ | Pearson |
|---|---:|---:|
| tore plat 2D | 0.967229 | 0.983478 |
| sphère | 0.982252 | 0.991086 |

La relation signée est également presque parfaite, à la convention de signe près :

| Géométrie | $R^2$ signée | Pearson signé |
|---|---:|---:|
| tore plat 2D | 0.998336 | -0.999168 |
| sphère | 0.998857 | -0.999428 |

Par conséquent :

$$
\boxed{
\delta\langle K_A\rangle
\text{ se comporte comme une source effective pour la réponse de courbure.}
}
$$

---

## 8. Interprétation

Paper 18 possède deux piliers numériques positifs :

1. La première loi modulaire du graphe est vérifiée avec une pente d'environ $0.989$ et un $R^2\simeq0.9966$.
2. La réponse modulaire prédit une réponse de courbure localisée avec un $R^2$ de $0.967$ à $0.982$.

Ainsi :

$$
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r).
$$

Cela soutient la chaîne du côté source nécessaire pour la limite einsteinienne effective.

---

## 9. Limites

Les résultats actuels sont positifs mais pas définitifs.

1. La source est scalaire / modulaire, pas encore un tenseur complet $T_{\mu\nu}$.
2. La matrice densité $\rho_A$ est une grandeur proxy de graphe.
3. L'hamiltonien modulaire $K_A=-\log\rho_A$ n'est pas encore dérivé d'une véritable matrice densité de sous-système quantique dans ce test.
4. La relation est testée sur des géométries contrôlées, pas encore sur de véritables graphes d'information mutuelle quantique.
5. Le signe de la relation signée dépend de la convention de perturbation.
6. Une dérivation continue via la première loi modulaire reste à écrire.

---

## 10. Structure des dossiers

```text
papers/paper18_entanglement_stress_tensor/
  README.md
  paper18_entanglement_stress_tensor.tex

  scripts/
    paper18_graph_modular_first_law_v1.py
    paper18_modular_source_curvature_v2.py
    paper18_build_entanglement_stress_summary_v1.py

  results/
    graph_modular_first_law_v1/
    modular_source_curvature_v2/
    paper18_entanglement_stress_summary_v1/

  figures/
    # figures finales copiées depuis les dossiers de résultats sélectionnés

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

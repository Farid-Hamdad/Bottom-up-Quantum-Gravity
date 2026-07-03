# Paper 19 — Équations d'Einstein effectives à partir de l'équilibre d'intrication

**Assemblage de la limite continue de Bottom-Up Quantum Gravity**

---

## Statut

Paper 19 assemble les trois flèches continues établies dans les Papers 16 à 18 :

$$
L_N\to\Delta_g,
$$

$$
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
$$

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

L'équation effective cible est :

$$
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}
}
$$

où :

- \(g_{\mu\nu}^{\rm ent}\) est la métrique d'intrication émergente ;
- \(T_{\mu\nu}^{\rm ent}\) est la source modulaire de tenseur énergie-impulsion d'intrication ;
- \(\Lambda_{\rm ent}\) est un terme cosmologique d'intrication ;
- \(G_{\rm eff}\) est le couplage gravitationnel effectif ;
- \(\mathcal{H}_{\mu\nu}\) contient des corrections d'ordre supérieur, non locales et à échelle finie.

Dans la limite lisse de basse énergie :

$$
\mathcal{H}_{\mu\nu}\to0.
$$

Alors BuP se réduit à la forme effective d'Einstein :

$$
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
$$

---

## 1. Point de départ

L'équation BuP fondamentale est la condition d'équilibre variationnelle :

$$
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0.
$$

L'action discrète est :

$$
S_{\rm BuP}[W]
=
\mathrm{Tr}\,L(W)^{-\beta[W]}
+
\sum_{ij}W_{ij}\kappa_{ij}[W]
+
\lambda\sum_{ij}W_{ij}d_{ij}^{2}
+
S_{\rm topo}[W].
$$

Chaque terme a une interprétation continue :

| Terme discret | Rôle continu |
|---------------|--------------|
| \(\mathrm{Tr}\,L^{-\beta}\) | géométrie spectrale / secteur d'Einstein-Hilbert |
| \(\sum_{ij}W_{ij}\kappa_{ij}\) | secteur de courbure de Ricci |
| \(\lambda\sum_{ij}W_{ij}d_{ij}^2\) | secteur de localité et cosmologique |
| \(S_{\rm topo}[W]\) | topologie et contraintes globales |
| \(\delta\langle K_A\rangle\) | source / secteur du tenseur énergie-impulsion |

La limite continue de cet équilibre est l'objet central de Paper 19.

---

## 2. Pilier I — Géométrie spectrale de Paper 16

Paper 16 soutient la limite spectrale continue :

$$
c_NL_N\to-\Delta_g.
$$

Les tests du bas du spectre ont donné :

| Géométrie | Erreur spectrale relative moyenne | \(\lambda_1^{\rm mis\ à\ l'échelle}\) | Cible |
|-----------|---:|---:|---:|
| \(S^1\) | 0.005072 | 1.006212 | 1.000000 |
| \(T^2\) | 0.026627 | 40.961270 | 39.478418 |

Cela établit que le laplacien d'intrication reconstruit le bas du spectre de l'opérateur de Laplace-Beltrami.

Ainsi :

$$
L(W)\to-\Delta_g.
$$

---

## 3. Pilier II — Courbure de Ricci de Paper 17

Paper 17 soutient la limite de courbure de Ricci :

$$
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_N R_{\mu\nu}u^\mu u^\nu.
$$

Les résultats clés sont :

$$
\Delta\bar\kappa_{\rm sphere-flat}=0.013383,
$$

$$
\Delta\left\langle\frac{\kappa}{\epsilon}\right\rangle_{\rm sphere-flat}=0.197719,
$$

$$
\Delta\left\langle\frac{\kappa}{\ell^2}\right\rangle_{\rm sphere-flat}=0.557041.
$$

L'étalonnage affine à \(N=512\) donne :

$$
B_N=-0.301281,
$$

$$
C_N=0.391933,
$$

$$
A_N=\frac{1}{C_N}=2.551455.
$$

Par conséquent :

$$
\frac{\kappa^{OR}}{\epsilon}
\simeq
-0.301
+
0.392\,R_{\mu\nu}u^\mu u^\nu.
$$

Cela établit que la courbure d'Ollivier-Ricci discrète porte un signal de Ricci moyen étalonné.

---

## 4. Pilier III — Source modulaire de Paper 18

Paper 18 soutient la chaîne de source :

$$
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r).
$$

La première loi modulaire du graphe est validée :

| Géométrie | Pente \(\delta S_A\) vs \(\delta\langle K_A\rangle\) | \(R^2\) |
|-----------|---:|---:|
| tore plat | 0.989481 | 0.996590 |
| sphère | 0.989718 | 0.996592 |

La source modulaire prédit la réponse de courbure près de la source :

| Géométrie | Rapport de localisation near/far | \(R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})\) |
|-----------|---:|---:|
| tore plat | 6.293007 | 0.967229 |
| sphère | 9.533529 | 0.982252 |

Les corrélations de Pearson signées sont :

$$
r=-0.999168
$$

pour le tore plat, et

$$
r=-0.999428
$$

pour la sphère.

Ainsi, \(\delta\langle K_A\rangle\) se comporte comme une source modulaire effective pour la courbure.

---

## 5. Dictionnaire continu effectif

Paper 19 utilise le dictionnaire suivant :

| Objet BuP discret | Objet continu |
|-------------------|---------------|
| \(W_{ij}=I(i:j)\) | métrique d'intrication \(g_{\mu\nu}^{\rm ent}\) |
| \(L(W)\) | \(-\Delta_g\) |
| \(\mathrm{Tr}\,L^{-\beta}\) | action gravitationnelle spectrale |
| \(\kappa_{ij}^{OR}\) | \(R_{\mu\nu}u^\mu u^\nu\) |
| \(\delta\langle K_A\rangle\) | source modulaire / \(T_{\mu\nu}^{\rm ent}\) |
| pénalité de localité | secteur cosmologique / infrarouge |
| \(S_{\rm topo}[W]\) | contraintes topologiques et globales |
| corrections de graphe fini | \(\mathcal{H}_{\mu\nu}\) |

---

## 6. Équation continue

La combinaison des trois piliers donne :

$$
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}
}
$$

où :

$$
G_{\mu\nu}
=
R_{\mu\nu}
-
\frac{1}{2}Rg_{\mu\nu}.
$$

Le tenseur de correction \(\mathcal{H}_{\mu\nu}\) inclut :

1. les corrections spectrales d'ordre supérieur ;
2. les corrections d'intrication non locales ;
3. les corrections induites par la topologie ;
4. les corrections à \(N\) fini ;
5. les écarts par rapport au comportement de variété lisse.

Dans la limite lisse de basse énergie :

$$
\mathcal{H}_{\mu\nu}\to0.
$$

Alors :

$$
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
$$

---

## 7. Interprétation

Paper 19 ne prétend pas que la gravité d'Einstein a été entièrement dérivée à partir de premiers principes.

Il établit un assemblage contrôlé :

$$
L_N\to\Delta_g,
$$

$$
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
$$

$$
\delta W_{\rm loc}\to T_{\mu\nu}^{\rm ent}.
$$

Ensemble, ils soutiennent l'existence d'un régime effectif d'Einstein à l'intérieur de BuP.

La tâche restante est analytique :

1. dériver les coefficients \(G_{\rm eff}\) et \(\Lambda_{\rm ent}\) ;
2. contrôler le tenseur de correction \(\mathcal{H}_{\mu\nu}\) ;
3. prouver la limite continue au-delà des géométries numériques contrôlées ;
4. reconstruire un tenseur complet \(T_{\mu\nu}^{\rm ent}\), pas seulement une source scalaire modulaire.

---

## 8. Structure des dossiers

```text
papers/paper19_effective_einstein_equations/
  README.md
  paper19_effective_einstein_equations.tex

  scripts/
    paper19_build_einstein_limit_summary_v1.py

  results/
    einstein_limit_summary_v1/
      paper19_einstein_limit_key_results.csv
      paper19_einstein_limit_summary.json
      paper19_einstein_limit_summary.md

  figures/
    # figures de synthèse finales

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

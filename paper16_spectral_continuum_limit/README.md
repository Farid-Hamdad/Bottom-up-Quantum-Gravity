# Paper 16 — Limite spectrale continue

**Du laplacien d'intrication BuP à l'opérateur de Laplace-Beltrami**

---

## Statut

Paper 16 étudie la première flèche analytique nécessaire à la limite einsteinienne continue de BuP :

$$
L_N
=
\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
\Delta_g.
$$

Paper 15 a fourni une preuve numérique au niveau de la trace de chaleur :

$$
d_s(t)
=
-2\frac{d\log{\rm Tr}(e^{-tL})}{d\log t}.
$$

Paper 16 renforce ce résultat en testant la convergence des valeurs propres individuelles :

$$
\lambda_k(L_N)
\longrightarrow
\lambda_k(\Delta_g).
$$

Les résultats actuels soutiennent la convergence du bas du spectre sur deux géométries contrôlées :

$$
S^1,
\qquad
T^2.
$$

---

## 1. Point de départ

L'objet fondamental est la matrice d'information mutuelle :

$$
W_{ij}=I(i:j).
$$

À partir d'elle, on construit la matrice des degrés du graphe :

$$
D_{ii}=\sum_j W_{ij},
$$

et le laplacien d'intrication non mis à l'échelle :

$$
L_N^0=D_N-W_N.
$$

L'opérateur mis à l'échelle continue utilisé dans Paper 15 et Paper 16 est :

$$
L_N
=
\frac{D_N-W_N}{\epsilon_N}.
$$

L'opérateur continu cible est l'opérateur de Laplace-Beltrami :

$$
\Delta_g
=
\frac{1}{\sqrt{|g|}}
\partial_\mu
\left(
\sqrt{|g|}g^{\mu\nu}\partial_\nu
\right).
$$

Étant donné que le laplacien de graphe brut possède encore une normalisation multiplicative inconnue, les tests spectraux ajustent actuellement un facteur scalaire \(c_N\) :

$$
c_N L_N
\longrightarrow
-\Delta_g.
$$

L'une des principales tâches théoriques de Paper 16 est de dériver \(c_N\) analytiquement plutôt que de l'ajuster.

---

## 2. Lien avec Paper 15

Paper 15 a validé la dimension spectrale par trace de chaleur en utilisant le laplacien renormalisé :

$$
L_\epsilon=\frac{D-W}{\epsilon}.
$$

Les principaux résultats de Paper 15 étaient :

| Géométrie | Dimension cible \(d_s\) | \(d_s\) mesurée | Erreur |
|-----------|:---------------------:|:-------------:|:------:|
| cercle    | 1                     | 1.0021        | 0.0021 |
| intervalle | 1                    | 0.9445        | 0.0555 |
| grille 2D | 2                     | 1.9938        | 0.0062 |
| sphère    | 2                     | 1.8964        | 0.1036 |

Cela a soutenu :

$$
L_N\to\Delta_g
$$

au niveau de la trace de chaleur / dimension.

Paper 16 teste maintenant l'affirmation plus forte :

$$
\lambda_k(c_NL_N)
\simeq
\lambda_k(-\Delta_g).
$$

---

## 3. Principaux résultats numériques

Paper 16 possède actuellement deux tests positifs au niveau des valeurs propres.

| Test | Limite | \(N\) | Modes comparés | Erreur rel. moyenne | Erreur rel. médiane | Erreur rel. max | \(\lambda_1^{\rm mis\ à\ l'échelle}\) | \(\lambda_1^{\rm cible}\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Spectre du cercle | \(c_NL_N\to-\Delta_{S^1}\) | 1024 | 12 | 0.0051 | 0.0051 | 0.0107 | 1.0062 | 1.0000 |
| Spectre du tore plat | \(c_NL_N\to-\Delta_{T^2}\) | 1024 | 20 | 0.0266 | 0.0263 | 0.0522 | 40.9613 | 39.4784 |

Ainsi, Paper 16 renforce Paper 15 :

$$
d_s\to D
\quad
\Longrightarrow
\quad
\lambda_k(L_N)\to\lambda_k(\Delta_g).
$$

---

## 4. Test 1 — Spectre du cercle

Pour le cercle unité \(S^1\), le spectre positif analytique de l'opérateur de Laplace-Beltrami est :

$$
\lambda_m=m^2,
$$

avec dégénérescence sinus/cosinus :

$$
1,1,4,4,9,9,16,16,\ldots
$$

Le graphe est construit à partir de \(N\) points également espacés sur le cercle, en utilisant la distance géodésique intrinsèque :

$$
d(\theta_i,\theta_j)
=
\min(|\theta_i-\theta_j|,2\pi-|\theta_i-\theta_j|).
$$

Les poids sont :

$$
W_{ij}
=
\exp\left[
-\frac{d(\theta_i,\theta_j)^2}{4\epsilon_N}
\right]
$$

sur un graphe \(k\)-NN local, et l'opérateur est :

$$
L_N=\frac{D_N-W_N}{\epsilon_N}.
$$

Une normalisation scalaire \(c_N\) est ajustée de sorte que :

$$
c_N\lambda_k^{\rm graphe}
\simeq
\lambda_k^{S^1}.
$$

Les erreurs diminuent systématiquement avec \(N\) :

| \(N\) | Erreur rel. moyenne | Erreur rel. médiane | Erreur rel. max |
|---:|---:|---:|---:|
| 64 | 0.0937 | 0.0952 | 0.1852 |
| 128 | 0.0492 | 0.0480 | 0.1003 |
| 256 | 0.0214 | 0.0215 | 0.0446 |
| 512 | 0.0109 | 0.0094 | 0.0249 |
| 1024 | 0.0051 | 0.0051 | 0.0107 |

À \(N=1024\) :

$$
\lambda_1^{\rm mis\ à\ l'échelle}=1.0062,
$$

tandis que la valeur analytique est :

$$
\lambda_1^{S^1}=1.
$$

Cela fournit une preuve au niveau des valeurs propres pour :

$$
c_NL_N\longrightarrow-\Delta_{S^1}.
$$

---

## 5. Test 2 — Spectre du tore plat

Pour le tore plat unité :

$$
T^2=[0,1)^2,
$$

avec conditions aux limites périodiques, le spectre analytique est :

$$
\lambda_{m,n}
=
4\pi^2(m^2+n^2),
\qquad
(m,n)\in\mathbb{Z}^2\setminus\{(0,0)\}.
$$

Le graphe est construit à partir d'une grille périodique \(m\times m\) en utilisant la distance périodique intrinsèque.

Là encore, l'opérateur est :

$$
L_N=\frac{D_N-W_N}{\epsilon_N},
$$

et un scalaire \(c_N\) est ajusté :

$$
c_N\lambda_k^{\rm graphe}
\simeq
\lambda_k^{T^2}.
$$

Les erreurs diminuent avec la résolution :

| \(N\) | Erreur rel. moyenne | Erreur rel. médiane | Erreur rel. max |
|---:|---:|---:|---:|
| 64 | 0.0951 | 0.1340 | 0.1467 |
| 144 | 0.0757 | 0.0437 | 0.1367 |
| 256 | 0.0604 | 0.0574 | 0.1850 |
| 576 | 0.0334 | 0.0455 | 0.0504 |
| 1024 | 0.0266 | 0.0263 | 0.0522 |

À \(N=1024\) :

$$
\lambda_1^{\rm mis\ à\ l'échelle}=40.9613,
$$

tandis que la valeur analytique est :

$$
\lambda_1^{T^2}=4\pi^2=39.4784.
$$

Cela fournit une preuve au niveau des valeurs propres pour :

$$
c_NL_N\longrightarrow-\Delta_{T^2}.
$$

---

## 6. Interprétation

Les résultats de Paper 16 montrent que le laplacien de graphe renormalisé de BuP ne reproduit pas seulement la dimension globale par trace de chaleur. Il reconstruit également le bas du spectre de l'opérateur de Laplace-Beltrami correspondant.

Preuve actuelle :

$$
S^1:
\quad
{\rm erreur\ relative\ moyenne}=0.0051,
$$

$$
T^2:
\quad
{\rm erreur\ relative\ moyenne}=0.0266.
$$

Ainsi, la limite spectrale continue est soutenue en dimension 1 et en dimension 2 :

$$
c_N\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
-\Delta_g.
$$

Le problème théorique restant est d'identifier la normalisation correcte \(c_N\) et de prouver la convergence sous des hypothèses contrôlées.

---

## 7. Théorème cible

Un théorème cible possible pour Paper 16 est :

> Soit \(M\) une variété riemannienne compacte de dimension \(D\) échantillonnée par des points \(x_i\). Supposons que les poids BuP satisfassent localement :
>
> $$
> W_{ij}
> =
> \exp\left[
> -\frac{d_g(x_i,x_j)^2}{4\epsilon_N}
> \right]
> +
> o(1).
> $$
>
> Si :
>
> $$
> \epsilon_N\to0,
> \qquad
> N\epsilon_N^{D/2}\to\infty,
> $$
>
> alors il existe une normalisation \(c_N\) telle que :
>
> $$
> c_N\frac{D_N-W_N}{\epsilon_N}
> \longrightarrow
> -\Delta_g
> $$
>
> dans un sens spectral ou de noyau de chaleur contrôlé.

Le mode de convergence exact reste à choisir.

---

## 8. Modes de convergence possibles

Paper 16 peut utiliser un ou plusieurs des modes suivants :

1. convergence ponctuelle sur des fonctions tests lisses ;
2. convergence des formes quadratiques ;
3. convergence de Mosco ;
4. convergence des valeurs propres ;
5. convergence du noyau de chaleur ;
6. convergence de la trace de chaleur ;
7. convergence de la dimension spectrale.

---

## 9. Questions théoriques ouvertes

Les principales questions ouvertes sont :

1. Quelle est l'expression analytique de \(c_N\) ?
2. Comment \(c_N\) dépend-il de \(D\), \(\epsilon_N\), \(k_N\) et de la densité d'échantillonnage ?
3. Est-ce que \(W_{ij}=I(i:j)\) approxime naturellement un noyau de chaleur ?
4. Comment la distance d'intrication \(d_{\rm ent}\) est-elle reliée à \(d_g\) ?
5. Comment traiter l'échantillonnage non uniforme ?
6. Comment traiter les bords ?
7. La preuve peut-elle être étendue des variétés synthétiques aux véritables graphes d'information mutuelle quantique ?

---

## 10. Structure des dossiers

```text
papers/paper16_spectral_continuum_limit/
  README.md
  paper16_spectral_continuum_limit.tex

  scripts/
    paper16_circle_spectrum_convergence_v1.py
    paper16_flat_torus_spectrum_convergence_v1.py
    paper16_build_spectral_summary_v1.py

  results/
    circle_spectrum_convergence_v1/
    flat_torus_spectrum_convergence_v1/
    paper16_spectral_summary_v1/

  figures/
    # figures finales copiées depuis les dossiers de résultats sélectionnés

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

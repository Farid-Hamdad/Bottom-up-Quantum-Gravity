# Paper 13 — Fermeture Microscopique de la Densité Baryonique dans BuP

## Objectif

Le Paper 13 teste la fermeture microscopique de l'hypothèse BuP utilisée dans le Paper 12 :

$$
\Sigma(R) \simeq \rho_{\rm ent}(R).
$$

Le Paper 12 utilisait la densité surfacique baryonique $\Sigma(R)$ comme une entrée effective pour construire un graphe inspiré de l'intrication et en déduire les observables de rotation galactique.

Le Paper 13 demande si cette identification peut être reconstruite à partir d'un état quantique microscopique :

$$
|\Psi_{\rm gal}\rangle
\longrightarrow I_{ij}
\longrightarrow \rho_{\rm ent}^{\rm micro}(R)
\propto \Sigma(R).
$$

La question centrale est :

> La densité baryonique observée peut-elle être interprétée comme la projection macroscopique d'une densité d'intrication microscopique localisée ?

---

## Chaîne de fermeture microscopique

La chaîne de fermeture complète testée dans ce papier est :

$$
\Sigma(R)
\longrightarrow J_{ij}
\longrightarrow H
\longrightarrow |\Psi_{\rm gal}\rangle
\longrightarrow I_{ij}
\longrightarrow \rho_{\rm ent}^{\rm micro}(R)
\longrightarrow \Sigma(R).
$$

La matrice de couplage microscopique est définie par

$$
J_{ij}
=
J_0 \sqrt{\Sigma_i\Sigma_j}
\exp\left(-\frac{d_{ij}}{\xi}\right),
$$

où :

- $d_{ij}$ est la distance spatiale entre les cellules du disque,
- $\xi$ est la longueur de corrélation microscopique,
- $J_0$ est l'amplitude du couplage,
- $\Sigma_i$ est la densité baryonique cible au site $i$.

L'hamiltonien est un modèle corrélé de type transverse-Ising :

$$
H
=
-\sum_{i<j}J_{ij}X_iX_j
-
h_0\sum_i Z_i.
$$

À partir de l'état fondamental $|\Psi_{\rm gal}\rangle$, l'information mutuelle par paire est calculée comme

$$
I_{ij}
=
S_i+S_j-S_{ij}.
$$

La densité d'intrication microscopique locale est alors définie par

$$
\rho_{\rm ent}^{\rm micro}(i)
=
\sum_j I_{ij}.
$$

Après un lissage radial grossier, le test de fermeture est :

$$
\rho_{\rm ent}^{\rm micro}(R)
\propto
\Sigma(R).
$$

---

## Résultat principal du balayage

Un balayage affiné a été réalisé sur 110 points de paramètres dans le plan $(\xi, h_0)$.

Résultats :

| Quantité | Valeur |
|---|---:|
| Nombre total d'exécutions | 110 |
| Fermeture forte | 42 / 110 |
| Fermeture modérée | 25 / 110 |
| Fraction de fermeture forte | 38,2% |
| Fraction de fermeture modérée ou forte | 60,9% |

Le meilleur point est :

$$
\xi = 4,25,
\qquad
h_0 = 1,2.
$$

À ce point :

$$
\mathrm{corr}
(\rho_{\rm ent}^{\rm micro},\Sigma)
=
0,997608,
$$

$$
\mathrm{RMSE}
=
0,020802,
$$

$$
R_{\rm ent}
=
3,002641,
\qquad
R_d
=
3,000000.
$$

L'erreur relative sur l'échelle de longueur est donc :

$$
\frac{|R_{\rm ent}-R_d|}{R_d}
=
8,80 \times 10^{-4},
$$

soit inférieure à $0,1\%$.

---

## Crête de fermeture critique

Les points de fermeture forte forment une crête diagonale dans le plan $(\xi, h_0)$.

Pour le sous-ensemble à fermeture forte :

$$
\left\langle \frac{\xi}{h_0} \right\rangle
=
4,35,
$$

et

$$
\left\langle
\frac{J_{\rm eff}^{\rm spec}}{h_0}
\right\rangle
=
1,15.
$$

Ici $J_{\rm eff}^{\rm spec}$ est le rayon spectral de la matrice de couplage microscopique $J_{ij}$.

Cela suggère que la fermeture microscopique se produit près d'un équilibre critique entre la force de couplage effective et le champ transverse.

---

## Contrôles structurels

Plusieurs contrôles ont été réalisés.

### Absence de couplage : $J_0 = 0$

Lorsque $J_0 = 0$, l'information mutuelle s'annule :

$$
\langle I_{ij}\rangle = 0,
\qquad
\max(I_{ij}) = 0.
$$

Aucune densité d'intrication microscopique ne peut être reconstruite.

### Densité inversée

Le contrôle `inverted_sigma` inverse l'organisation radiale de la densité entrant dans la matrice de couplage.

Résultat :

$$
0/16
$$

points atteignent une fermeture modérée ou forte.

### Couplage aléatoire

Le contrôle `random_J` remplace la matrice de couplage structurée par une matrice aléatoire symétrique d'échelle comparable.

Résultat :

$$
0/16
$$

points atteignent une fermeture modérée ou forte.

### Couplage uniquement géométrique

Le contrôle `geometric_J` utilise

$$
J_{ij}=J_0 e^{-d_{ij}/\xi},
$$

sans le facteur de pondération baryonique $\sqrt{\Sigma_i\Sigma_j}$.

Résultat :

| Quantité | Valeur |
|---|---:|
| Fermeture forte | 1 / 16 |
| Fermeture modérée ou forte | 2 / 16 |

Cela montre que la géométrie seule peut occasionnellement produire une fermeture acceptable, mais ne génère pas de crête de fermeture robuste.

### Densité mélangée

Le contrôle `shuffled_sigma` permute aléatoirement les valeurs de densité entrant dans les couplages.

Sur 10 graines aléatoires :

| Quantité | Valeur |
|---|---:|
| Nombre total d'exécutions | 160 |
| Fermeture forte | 28 / 160 |
| Fermeture modérée ou forte | 42 / 160 |
| Fraction de fermeture forte | 17,5% |
| Fraction de fermeture modérée ou forte | 26,2% |

Cela montre que la fermeture est affaiblie par le désordre local, mais pas complètement détruite à petite échelle.

---

## Robustesse face à la taille finie

Le balayage a été répété pour $N = 8$, $N = 12$ et $N = 16$ qubits.

| N | Dimension de Hilbert | Fermeture forte | Modérée+forte | Meilleur $R_{\rm ent}$ | Erreur sur $R_d$ |
|---:|---:|---:|---:|---:|---:|
| 8 | 256 | 3 / 16 | 7 / 16 | 2,884991 | 3,83% |
| 12 | 4096 | 6 / 16 | 10 / 16 | 3,095187 | 3,17% |
| 16 | 65536 | 4 / 16 | 8 / 16 | 2,959983 | 1,33% |

Le phénomène de fermeture persiste jusqu'à $N = 16$, et la meilleure reconstruction de l'échelle radiale s'améliore à la plus grande taille testée.

---

## Résumé conceptuel

Les résultats numériques soutiennent l'énoncé de fermeture microscopique suivant :

$$
\boxed{
\text{La matière encode l'intrication.}
\qquad
\text{L'intrication reconstruit la matière.}
}
$$

Plus précisément :

$$
\boxed{
\text{La densité baryonique apparaît comme la trace macroscopique}
\quad
\text{d'une intrication critique spatialement organisée.}
}
$$

---

## Contenu du dépôt

```text
papers/paper13_micro_closure_bup/
  README.md
  paper13_micro_closure_bup.tex
  scripts/
    bup_paper13_micro_closure_scan_v1.py
    bup_paper13_micro_closure_scan_v2.py
    bup_paper13_micro_closure_scan_v2_geometric.py
  results/
    scan_v1/
    scan_v2_zoom/
    control_no_coupling/
    control_flat_sigma/
    control_shuffled_sigma_v2/
    control_inverted_sigma_v2/
    control_random_J_v2/
    control_geometric_J_v2/
    finite_size_N8/
    finite_size_N12/
    finite_size_N16/
  figures/
    fig_paper13_verdict_map.png
    fig_paper13_rmse_map.png
    fig_paper13_rd_error_map.png
    fig_paper13_jeff_over_h0.png

# Paper 13 — Fermeture Microscopique de la Densité Baryonique dans BuP

## Objectif

Le Paper 13 teste la fermeture microscopique de l'hypothèse BuP utilisée dans le Paper 12 :

$$
\boxed{\;\Sigma(R) \;\simeq\; \rho_{\text{ent}}(R)\;}
$$

Le Paper 12 utilisait la densité surfacique baryonique $\Sigma(R)$ comme une entrée effective pour construire un graphe inspiré de l'intrication et en déduire les observables de rotation galactique.

Le Paper 13 demande si cette identification peut être reconstruite à partir d'un état quantique microscopique :

$$
|\Psi_{\text{gal}}\rangle
\;\longrightarrow\;
I_{ij}
\;\longrightarrow\;
\rho_{\text{ent}}^{\text{micro}}(R)
\;\propto\;
\Sigma(R).
$$

**Question centrale :** La densité baryonique observée peut-elle être interprétée comme la projection macroscopique d'une densité d'intrication microscopique localisée ?

---

## Chaîne de fermeture microscopique

La chaîne de fermeture complète testée dans ce papier est :

$$
\Sigma(R)
\;\longrightarrow\;
J_{ij}
\;\longrightarrow\;
H
\;\longrightarrow\;
|\Psi_{\text{gal}}\rangle
\;\longrightarrow\;
I_{ij}
\;\longrightarrow\;
\rho_{\text{ent}}^{\text{micro}}(R)
\;\longrightarrow\;
\Sigma(R).
$$

### Matrice de couplage microscopique

$$
J_{ij} \;=\; J_0 \;\sqrt{\Sigma_i \,\Sigma_j} \;\exp\!\left(-\frac{d_{ij}}{\xi}\right)
$$

avec :
- $d_{ij}$ = distance spatiale entre les cellules du disque
- $\xi$ = longueur de corrélation microscopique
- $J_0$ = amplitude du couplage
- $\Sigma_i$ = densité baryonique cible au site $i$

### toto

$$
\mathcal{H} = - \sum_{i < j} J_{ij} \, \hat{X}_i \hat{X}_j - h_0 \sum_{i} \hat{Z}_i
$$


markdown$$
H = - \sum_{i<j} J_{ij} X_i X_j - h_0 \sum_i Z_i
$$

$$
H = - \sum_{i<j} J_{ij} X_i X_j - h_0 \sum_i Z_i
$$

$$
H = - \sum_{i<j} J_{ij} X_i X_j - h_0 \sum_i Z_i
$$

où :
- $X_i$, $Z_i$ sont les matrices de Pauli
- $J_{ij}$ est la matrice de couplage
- $h_0$ est le champ transverse

### Information mutuelle

À partir de l'état fondamental $|\Psi_{\text{gal}}\rangle$ :

$$
I_{ij} = S_i + S_j - S_{ij}
$$

### Densité d'intrication microscopique locale

$$
\rho_{\text{ent}}^{\text{micro}}(i) \;=\; \sum_{j} I_{ij}
$$

### Test de fermeture (après lissage radial)

$$
\boxed{\;\rho_{\text{ent}}^{\text{micro}}(R) \;\propto\; \Sigma(R)\;}
$$

---

## Résultat principal du balayage

Un balayage affiné a été réalisé sur **110 points** de paramètres dans le plan $(\xi, h_0)$.

| Quantité | Valeur |
|----------|-------:|
| Nombre total d'exécutions | 110 |
| Fermeture forte | 42 / 110 |
| Fermeture modérée | 25 / 110 |
| Fraction de fermeture forte | 38,2 % |
| Fraction de fermeture modérée ou forte | 60,9 % |

### Meilleur point

$$
\xi = 4{,}25 \qquad\text{et}\qquad h_0 = 1{,}2
$$

À ce point :

| Mesure | Valeur |
|--------|-------:|
| Corrélation $\text{corr}(\rho_{\text{ent}}^{\text{micro}},\Sigma)$ | $0,997608$ |
| RMSE | $0,020802$ |
| $R_{\text{ent}}$ | $3,002641$ |
| $R_d$ (vraie échelle) | $3,000000$ |

**Erreur relative sur l'échelle de longueur :**

$$
\frac{|R_{\text{ent}} - R_d|}{R_d} \;=\; 8{,}80 \times 10^{-4}
\;\;(\text{soit } 0{,}088\%)
$$

---

## Crête de fermeture critique

Les points de fermeture forte forment une **crête diagonale** dans le plan $(\xi, h_0)$.

Pour le sous-ensemble à fermeture forte :

$$
\left\langle \frac{\xi}{h_0} \right\rangle \;=\; 4{,}35
\qquad\text{et}\qquad
\left\langle \frac{J_{\text{eff}}^{\text{spec}}}{h_0} \right\rangle \;=\; 1{,}15
$$

où $J_{\text{eff}}^{\text{spec}}$ est le rayon spectral de la matrice $J_{ij}$.

> **Interprétation :** La fermeture microscopique se produit près d'un équilibre critique entre la force de couplage effective et le champ transverse.

---

## Contrôles structurels

Plusieurs contrôles ont été réalisés pour valider la robustesse du résultat.

### 1. Absence de couplage ($J_0 = 0$)

L'information mutuelle s'annule complètement :

$$
\langle I_{ij} \rangle = 0 \qquad\text{et}\qquad \max(I_{ij}) = 0
$$

✅ **Aucune reconstruction possible** — le couplage est nécessaire.

### 2. Densité inversée (`inverted_sigma`)

On inverse l'organisation radiale de la densité.

$$
\text{Résultat : } 0/16 \text{ points atteignent une fermeture modérée ou forte}
$$

✅ **L'organisation radiale est cruciale.**

### 3. Couplage aléatoire (`random_J`)

La matrice structurée est remplacée par une matrice aléatoire symétrique.

$$
\text{Résultat : } 0/16 \text{ points atteignent une fermeture modérée ou forte}
$$

✅ **La structure spécifique du couplage importe.**

### 4. Couplage uniquement géométrique (`geometric_J`)

$$
J_{ij} = J_0 \; e^{-d_{ij}/\xi} \quad\text{(sans le facteur } \sqrt{\Sigma_i\Sigma_j}\text{)}
$$

| Quantité | Valeur |
|----------|-------:|
| Fermeture forte | 1 / 16 |
| Fermeture modérée ou forte | 2 / 16 |

✅ **La géométrie seule ne suffit pas** — le couplage baryonique est nécessaire pour une crête robuste.

### 5. Densité mélangée (`shuffled_sigma`)

On permute aléatoirement les valeurs de densité (10 graines aléatoires).

| Quantité | Valeur |
|----------|-------:|
| Nombre total d'exécutions | 160 |
| Fermeture forte | 28 / 160 |
| Fermeture modérée ou forte | 42 / 160 |
| Fraction de fermeture forte | 17,5 % |
| Fraction de fermeture modérée ou forte | 26,2 % |

✅ **Le désordre local affaiblit la fermeture**, mais ne la détruit pas complètement à petite échelle.

---

## Robustesse face à la taille finie

Le balayage a été répété pour $N = 8$, $N = 12$ et $N = 16$ qubits.

| $N$ | Dimension de Hilbert | Fermeture forte | Modérée ou forte | Meilleur $R_{\text{ent}}$ | Erreur sur $R_d$ |
|:---:|---------------------:|----------------:|-----------------:|--------------------------:|-----------------:|
| 8 | 256 | 3 / 16 | 7 / 16 | 2,884991 | 3,83 % |
| 12 | 4 096 | 6 / 16 | 10 / 16 | 3,095187 | 3,17 % |
| 16 | 65 536 | 4 / 16 | 8 / 16 | 2,959983 | 1,33 % |

> **Conclusion :** Le phénomène de fermeture persiste jusqu'à $N = 16$, et la reconstruction de l'échelle radiale s'améliore à la plus grande taille testée.

---

## Résumé conceptuel

Les résultats numériques soutiennent l'énoncé suivant :

$$
\boxed{\text{La matière encode l'intrication.} \qquad \text{L'intrication reconstruit la matière.}}
$$

Plus précisément :

$$
\boxed{\text{La densité baryonique apparaît comme la trace macroscopique d'une intrication critique spatialement organisée.}}
$$

---

## Contenu du dépôt

# Paper 28 — Modular Time and Mutual-Information Geometry

**Frozen TFIM benchmarks from N = 6 to N = 10**

**Titre FR :** Temps modulaire et géométrie d'information mutuelle  
**Sous-titre FR :** Benchmarks TFIM gelés de N = 6 à N = 10

---

## Statut

Paper 28 est le papier BuP consacré au **temps modulaire** et à son lien opérationnel avec la géométrie informationnelle.

Le point de départ est le sous-système réduit

\[
\rho_A=\operatorname{Tr}_{\bar A}|\Psi\rangle\langle\Psi|,
\]

son Hamiltonien modulaire

\[
K_A=-\log\rho_A,
\]

et le flot modulaire intrinsèque

\[
O(\tau)=e^{iK_A\tau}Oe^{-iK_A\tau}.
\]

La question testée est précise : **le profil de propagation généré par \(K_A\) est-il descriptivement plus proche de la géométrie d'information mutuelle \(W\) que de l'adjacence microscopique nue de la chaîne TFIM ?**

La réponse, dans les systèmes finis et le protocole gelé étudiés ici, est oui.

---

## Résultat principal

Pour chaque paire \((i,j)\) du sous-système \(A\), on définit

\[
C_{ij}(\tau)=\frac{1}{9}\sum_{a,b=X,Y,Z}
\operatorname{Tr}\!\left[
\rho_A[\sigma_i^a(\tau),\sigma_j^b]^\dagger
[\sigma_i^a(\tau),\sigma_j^b]
\right].
\]

À court temps,

\[
C_{ij}(\tau)=A_{ij}\tau^2+O(\tau^3),
\qquad
v_{ij}^{\rm mod}=\sqrt{A_{ij}}.
\]

Le profil est comparé à deux candidats :

1. la géométrie d'information mutuelle
   \[
   W_{ij}=I(i:j),
   \]
2. l'adjacence nue de la chaîne TFIM, \(G_{\rm chain}\).

Après ajustement d'une échelle positive, on utilise

\[
R(g)=\frac{\|v-c(g)g\|}{\|v\|},
\qquad
c(g)=\max\!\left(0,\frac{v\cdot g}{g\cdot g}\right),
\]

et

\[
\boxed{\Delta_R^{\rm flow}=R(G_{\rm chain})-R(W)}.
\]

Une valeur positive signifie que \(W\) donne le plus faible résidu descriptif.

### Primaire — court temps

| Système | Sous-système | \(\Delta_R^{\rm flow}\) |
|---|---:|---:|
| N6 | A3 | +0.151958 |
| N8 | A4 | +0.122924 |
| N10 | A4 contrôle | +0.124199 |
| N10 | A5 primaire | +0.220862 |

### Secondaire — grille gelée \(\tau\in[0,120]\)

La grille contient 1024 points, dont 1023 avec \(\tau>0\).

| Système | \(\langle\Delta_R\rangle\) | Fraction \(R_W<R_{chain}\) | \(\Delta_R^{\min}\) | \(\Delta_R^{\max}\) |
|---|---:|---:|---:|---:|
| N6/A3 | +0.132052 | 0.993157 | -0.030579 | +0.191262 |
| N8/A4 | +0.068709 | 0.991202 | -0.026039 | +0.225529 |
| N10/A4 | +0.069911 | 0.992180 | -0.024641 | +0.222309 |
| N10/A5 | **+0.147941** | **1.000000** | **+0.126910** | **+0.283654** |

Le bras N10/A5 est favorable à \(W\) sur **1023/1023** points de temps modulaire de la grille gelée.

---

## Contrôle fixe N8/A4 → N10/A4

Ce contrôle sépare partiellement l'effet de la taille globale \(N\) de celui de la taille du sous-système \(|A|\).

\[
\Delta_R^{\rm primary}:
0.1229235652\rightarrow0.1241991943,
\]

et

\[
\langle\Delta_R\rangle:
0.0687091263\rightarrow0.0699112761.
\]

La fraction de temps modulaire favorable à \(W\) passe de

\[
0.9912023460\rightarrow0.9921798631.
\]

Ainsi, à \(|A|=4\) fixé, le comportement change très peu lorsque l'environnement global passe de N=8 à N=10.

Cette stabilité est descriptive. Elle ne constitue pas une preuve de limite thermodynamique ou d'universalité.

---

## Interprétation dans BuP

Paper 28 clarifie la séparation entre les deux secteurs :

### Secteur spatial / informationnel

\[
W_{ij}=I(i:j)
\rightarrow
d_{ij}^{\rm ent}
\rightarrow
L_{\rm ent}=D-W
\rightarrow
\text{géométrie émergente}.
\]

### Secteur temporel / modulaire

\[
\rho_A
\rightarrow
K_A=-\log\rho_A
\rightarrow
O(\tau)=e^{iK_A\tau}Oe^{-iK_A\tau}.
\]

Le résultat de Paper 28 n'identifie pas ces deux objets. Il montre plutôt une **compatibilité dynamique** : dans les benchmarks étudiés, le profil de propagation produit par \(K_A\) s'aligne mieux sur \(W\) que sur la connectivité microscopique nue.

C'est le pont opérationnel recherché entre **temps modulaire** et **géométrie informationnelle**.

---

## Ce que Paper 28 n'établit pas

Paper 28 ne démontre pas :

- que le temps physique macroscopique est entièrement identifié au paramètre modulaire \(\tau\) ;
- une causalité relativiste ou un cône de lumière émergent ;
- les équations d'Einstein ;
- une limite continue ;
- une universalité au-delà du protocole TFIM étudié ;
- une signification statistique des fractions de grille secondaire.

Les résultats secondaires sont explicitement **descriptifs post-unblind**.

---

## Reproductibilité

Le dossier contient :

```text
paper28_modular_time_geometry/
├── README.md
├── main.tex
├── figures/
│   ├── fig01_primary_deltaR.tex
│   ├── fig02_secondary_mean_deltaR.tex
│   ├── fig03_fixed_A4_control.tex
│   └── fig04_W_win_fraction.tex
├── notes/
│   ├── claim_boundary.md
│   ├── provenance.md
│   └── reproducibility.md
├── results/paper28_modular_time_v1/
│   ├── frozen_reference_metrics.csv
│   ├── primary_metrics.csv
│   ├── secondary_metrics.csv
│   └── summary.json
└── scripts/
    ├── paper28_reproduce_modular_time_v1.py
    ├── paper28_verify_reproduction_v1.py
    ├── paper28_make_figures_v1.py
    └── run_all.sh
```

Reproduction :

```bash
cd paper28_modular_time_geometry
./scripts/run_all.sh
```

Le script reconstruit exactement les états TFIM finis, \(\rho_A\), \(W\), \(K_A\), les coefficients à court temps et les courbes secondaires. La vérification publiée compare ensuite les sorties aux valeurs gelées de la campagne d'origine avec une tolérance numérique de \(10^{-10}\).

---

## Auteur

**Farid Hamdad**  
Bottom-Up Quantum Gravity Program — 2026

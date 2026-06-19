# Paper 22 — Ondes gravitationnelles BuP

**Bottom-Up Quantum Gravity — Modes informationnels du graphe d’intrication**

Ce dossier contient le secteur « ondes gravitationnelles » de la théorie **Bottom-Up Quantum Gravity (BuP)**.

L’objectif du papier est d’étudier comment des modes tensoriels de type ondes gravitationnelles émergent à partir de perturbations du graphe d’intrication

```math
W_{ij}=I(i:j),
```

plutôt qu’à partir d’une métrique fondamentale de l’espace-temps.

Dans BuP, la métrique est reconstruite à partir de la structure informationnelle de l’état quantique. Les ondes gravitationnelles sont donc interprétées comme des excitations collectives du graphe d’intrication.

La chaîne dynamique finale établie dans ce papier est :

```math
S_{\rm flux}
\rightarrow
J_e(t)
\rightarrow
J_n(t)
\rightarrow
q_n(t)
\rightarrow
\delta W_{ij}(t)
\rightarrow
h_{\mu\nu}^{\rm eff}.
```

Ce résultat montre que BuP ne se contente pas de supporter des modes de type ondes gravitationnelles : la théorie les génère dynamiquement à partir de la source émergente de flux de matière.

---

## Statut du dossier

Ce dossier contient deux étapes du développement de Paper 22 :

1. **Prototype historique avec hessien nodal limité**
2. **Secteur final avec vrai hessien d’arêtes**

Le secteur historique est conservé pour reproductibilité. Les résultats finaux de Paper 22 sont ceux du vrai hessien d’arêtes, correspondant aux scripts v10–v17.

---

## 1. Prototype historique : hessien nodal limité

La première implémentation utilisait une réduction en champ nodal :

```math
W_{ij}(\phi)=W_{ij}^{(0)}
\exp\left(\frac{\phi_i+\phi_j}{2}\right),
```

avec une perturbation scalaire (\phi_i) par nœud.

Cette réduction donne un hessien nodal :

```math
H_{ij}
=
\frac{\partial^2 S_{\rm BuP}}{\partial \phi_i \partial \phi_j}.
```

Elle a permis de montrer que BuP supporte des modes hessiens de type onde, avec une relation de dispersion :

```math
\omega^2\simeq c_{\rm graph}^2 k^2+m_{\rm eff}^2.
```

Cette étape est conservée comme test préliminaire de compatibilité.

Cependant, la réduction nodale possède deux limites :

* elle tronque les vrais degrés de liberté d’arêtes du graphe d’intrication ;
* elle ne ferme pas la chaîne dynamique matière (\rightarrow) onde.

Emplacement recommandé :

```text
scripts/legacy/
results/legacy_v5_node_hessian/
figures/legacy_v5/
```

---

## 2. Secteur final : vrai hessien d’arêtes

L’implémentation finale utilise des variables d’arêtes normalisées :

```math
x_e=\frac{\delta W_e}{\sqrt{W_e^{(0)}}},
```

où (e=(i,j)) parcourt les arêtes du graphe d’intrication.

L’équation dynamique est :

```math
\ddot{x}_e+\gamma\dot{x}_e+\sum_f K_{ef}x_f=J_e(t),
```

où (K_{ef}) est le hessien dans le secteur des variables d’arêtes de l’action effective BuP, et (J_e(t)) est la projection de la source émergente (S_{\rm flux}) sur les arêtes du graphe.

En décomposant sur les modes propres du hessien :

```math
\ddot{q}_n+\gamma\dot{q}_n+\omega_n^2 q_n=J_n(t).
```

L’observable numérique centrale est la corrélation modale :

```math
\mathrm{corr}(|J_n|,q_n^{\rm peak}).
```

---

## Résultats numériques principaux

### v10 — Génération dynamique par vrai hessien d’arêtes

Meilleur run de génération haut-(k) :

```math
\eta_{\rm edge}=0.75,
\qquad
\text{source-width}=0.12,
\qquad
\text{pulse-}\sigma=0.25.
```

Résultats :

```math
\overline{\mathrm{corr}}(|J_n|,q_n^{\rm peak})
=
0.999718,
```

```math
\overline{E}_{\rm low25}=0.006216,
```

```math
\overline{k}_{\rm weighted}=4.1556.
```

Ce résultat valide la génération dynamique des modes du hessien d’arêtes par la source de flux BuP.

---

### v12 — Scaling en taille finie à (\eta_{\rm edge}=1)

Le point Einstein calibré a été testé pour :

```math
N=121,\ 256,\ 400,\ 625.
```

| (N) | (c_{\rm edge}) |   (v_g) | (R^2_{\rm disp}) | corrélation moyenne |
| --: | -------------: | ------: | ---------------: | ------------------: |
| 121 |       1.000020 | 0.99995 |              1.0 |            0.999806 |
| 256 |       1.000008 | 0.99998 |              1.0 |            0.999454 |
| 400 |       1.000003 | 0.99998 |              1.0 |            0.994903 |
| 625 |       0.999999 | 0.99998 |              1.0 |            0.975093 |

Le mécanisme de génération reste robuste jusqu’à (N=625).

---

### v13–v14 — Analyse du gap de masse

L’offset résiduel

```math
m_{\rm eff}^2\simeq2\times10^{-3}
```

a été testé contre des effets de discrétisation et de volume fini.

Le scan à (N) fixé et densité variable montre que le gap ne disparaît pas lorsque le pas effectif (a) varie.

Le scan à densité fixée et volume variable donne :

```math
m_{\rm eff}^2(L)
=
m_0^2+\frac{B}{L^2},
```

avec :

```math
m_0^2\simeq2.08\times10^{-3}.
```

Cela indique que le gap n’est pas un simple artefact de densité locale ni un pur effet infrarouge de volume fini dans la plage testée.

---

### v15–v16 — Origine spectrale du gap résiduel

Le gap résiduel corrèle fortement avec la première valeur propre non nulle du laplacien normalisé d’intrication :

```math
\lambda_1(\mathcal L_{\rm norm}).
```

La loi ajustée est :

```math
m_{\rm eff}^2
=
m_0^2
\left(
1-\frac{\lambda_1(\mathcal L_{\rm norm})}{\lambda_*}
\right),
```

avec :

```math
m_0^2\simeq2.03\times10^{-3},
\qquad
\lambda_*\simeq0.723.
```

Ce résultat suggère que le gap résiduel du hessien d’arêtes est contrôlé par un déficit spectral infrarouge relatif du graphe d’intrication.

---

### v17 — Calibration LVK

La vitesse des modes vérifie :

```math
c_{\rm edge}^2\simeq\eta_{\rm edge}.
```

Un scan local autour de (\eta_{\rm edge}=1) donne :

```math
\frac{dc_{\rm edge}}{d\eta_{\rm edge}}
\bigg|_{\eta=1}
\simeq0.5.
```

Donc :

```math
\frac{\delta c_{\rm GW}}{c}
\simeq
\frac12
\frac{\delta\eta_{\rm edge}}{\eta_*}.
```

La contrainte LVK :

```math
\left|\frac{c_{\rm GW}}{c}-1\right|
<
5\times10^{-16}
```

se traduit alors par :

```math
\left|
\frac{\eta_{\rm edge}}{\eta_*}-1
\right|
\lesssim10^{-15}.
```

Dans BuP, cette condition est interprétée comme une calibration du point fixe géométrique, et non comme un ajustement libre.

---

## Structure finale du point fixe

Paper 22 identifie une double condition critique :

```math
\eta_{\rm edge}=\eta_*,
\qquad
\lambda_1(\mathcal L_{\rm norm})=\lambda_*.
```

À ce point :

```math
c_{\rm GW}=c,
\qquad
m_{\rm eff}^2=0.
```

Avec :

```math
\eta_*=1,
\qquad
\lambda_*\simeq0.723.
```

Ce point définit le **point fixe tensoriel BuP–Einstein**.

---

## Lien avec Paper 21 — Point fixe SLACS

Paper 21 correspond au papier **SLACS fixed point**.

Paper 22 correspond au papier **ondes gravitationnelles BuP**.

Les deux papiers sont distincts, mais connectés.

Paper 21 a identifié un point fixe observationnel candidat autour de :

```math
\log M_\star\simeq11.6,
```

où plusieurs signaux SLACS convergent :

```math
C_{\rm obs}=0,
```

```math
\Phi_{\rm BuP}>\log M_\star,
```

```math
\alpha_{\rm eff}\simeq1.
```

Paper 22 ajoute le secteur tensoriel dynamique :

```math
c_{\rm GW}=c,
\qquad
m_{\rm eff}^2=0.
```

Ensemble, ces résultats suggèrent une convergence quadruple vers un point fixe BuP–Einstein.

Cette liaison est interprétée comme une cross-validation candidate, pas comme une preuve d’universalité.

---

## Organisation recommandée du dossier

```text
papers/paper22_bup_gravitational_waves/
├── README.md
├── main.tex
├── scripts/
│   ├── legacy/
│   │   └── bup_gw_hessian_modes_v5.py
│   └── true_edge/
│       ├── bup_gw_true_edge_hessian_v10.py
│       ├── scan_true_edge_hessian_v10.py
│       ├── scan_finite_size_true_edge_hessian_v12_eta1_fix.py
│       ├── scan_mass_gap_v13.py
│       ├── scan_fixed_density_mass_gap_v14.py
│       ├── correlate_mass_gap_laplacian_v15.py
│       ├── test_relative_spectral_gap_v16.py
│       └── scan_lvk_eta_calibration_v17.py
├── results/
│   ├── legacy_v5_node_hessian/
│   ├── true_edge_hessian_v10/
│   ├── scan_true_edge_hessian_v10/
│   ├── finite_size_true_edge_hessian_v12_eta1/
│   ├── mass_gap_v13_fixedN_density/
│   ├── mass_gap_v14_fixed_density/
│   ├── mass_gap_v15_laplacian_corr_combined/
│   ├── mass_gap_v16_relative_spectral/
│   └── lvk_eta_calibration_v17_fine/
├── figures/
│   ├── legacy_v5/
│   ├── true_edge_hessian/
│   ├── finite_size/
│   ├── mass_gap/
│   └── lvk_calibration/
└── archive/
```

---

## Commandes de reproduction

### v17 — Calibration LVK

```bash
cd ~/bottomup

python3 papers/paper22_bup_gravitational_waves/scripts/true_edge/scan_lvk_eta_calibration_v17.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/true_edge/bup_gw_true_edge_hessian_v10.py \
  --eta-edge-list 0.999,0.9995,0.9999,1.0,1.0001,1.0005,1.001 \
  --N-side 11 \
  --extent 1.0 \
  --source-width 0.12 \
  --pulse-sigma 0.25 \
  --output-dir papers/paper22_bup_gravitational_waves/results/lvk_eta_calibration_v17_fine
```

### v15 — Corrélation avec le gap spectral

```bash
cd ~/bottomup

python3 papers/paper22_bup_gravitational_waves/scripts/true_edge/correlate_mass_gap_laplacian_v15.py \
  --summary-csv papers/paper22_bup_gravitational_waves/results/mass_gap_v13_fixedN_density/results/mass_gap_v13_summary.csv \
  --summary-csv papers/paper22_bup_gravitational_waves/results/mass_gap_v14_fixed_density/results/mass_gap_v14_summary.csv \
  --output-dir papers/paper22_bup_gravitational_waves/results/mass_gap_v15_laplacian_corr_combined
```

### v12 — Scan en taille finie à (\eta_{\rm edge}=1)

```bash
cd ~/bottomup

python3 papers/paper22_bup_gravitational_waves/scripts/true_edge/scan_finite_size_true_edge_hessian_v12_eta1_fix.py \
  --v10-script papers/paper22_bup_gravitational_waves/scripts/true_edge/bup_gw_true_edge_hessian_v10.py \
  --N-side-list 11,16,20,25 \
  --eta-edge 1.0 \
  --source-width 0.12 \
  --pulse-sigma 0.25 \
  --skip-existing \
  --output-dir papers/paper22_bup_gravitational_waves/results/finite_size_true_edge_hessian_v12_eta1
```

---

## Limitations

Le vrai hessien d’arêtes utilisé ici est le hessien du secteur effectif en variables d’arêtes de l’action BuP. Un hessien complet par différences finies de toute l’action spectrale microscopique sur toutes les variables d’arêtes reste une extension future.

Les scans en taille finie couvrent (N=121) à (N=625). Des graphes plus grands seront nécessaires pour tester plus fortement la limite continue.

La loi spectrale :

```math
m_{\rm eff}^2
=
m_0^2
\left(
1-\frac{\lambda_1}{\lambda_*}
\right)
```

a été validée numériquement sur les points v13+v14, mais l’universalité de (\lambda_*\simeq0.723) sur d’autres familles de graphes reste à tester.

Le lien SLACS doit être interprété comme une cross-validation, pas comme une dérivation directe du point fixe tensoriel à partir des données de lentilles.

---

## Résumé

Paper 22 transforme le prototype historique avec hessien nodal limité en un mécanisme dynamique avec vrai hessien d’arêtes.

Le résultat final est :

```math
S_{\rm flux}
\rightarrow
J_e(t)
\rightarrow
J_n(t)
\rightarrow
q_n(t)
\rightarrow
\delta W_{ij}(t)
\rightarrow
h_{\mu\nu}^{\rm eff}.
```

Le régime tensoriel Einstein est atteint lorsque :

```math
\eta_{\rm edge}=\eta_*,
\qquad
\lambda_1(\mathcal L_{\rm norm})=\lambda_*,
```

avec :

```math
\eta_*=1,
\qquad
\lambda_*\simeq0.723.
```

À ce point fixe :

```math
c_{\rm GW}=c,
\qquad
m_{\rm eff}^2=0.
```

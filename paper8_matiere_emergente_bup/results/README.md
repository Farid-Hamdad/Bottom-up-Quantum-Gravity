# Résultats numériques — Paper 8

Ce dossier contient les résultats complets des simulations et analyses effectuées pour la construction et la validation du tenseur matière émergent dans BuP.

---

## Organisation des fichiers

| Fichier | Description |
|--------|-------------|
| `paper8_summary.csv` | Synthèse globale : meilleures sources, corrélations, paramètres optimaux |
| `localization_sigma_scan.csv` | Scan de la largeur $\sigma$ : rayon de localisation $R_{\rm loc}$ en fonction de $\sigma$ |
| `curvature_response_sigma_scan.csv` | Corrélation (Spearman/Pearson) entre $T_{00}$ et $\vert\delta R\vert$ pour différents $\sigma$ |
| `fluid_source_scan.csv` | Scan du coefficient $\omega$ dans $T_{00} + \omega T_{aa}$ – meilleur $\omega = -0.5$ |
| `grad_source_chi_scan.csv` | Scan du coefficient $\chi$ dans $T_{\rm grad}$ – meilleur $\chi = 0.5$ |
| `flux_source_psi_scan.csv` | Scan du coefficient du flux $\Vert T_{0a}\Vert$ – meilleur $\psi = 1.0$ |
| `conservation_sigma_scan.csv` | Divergence de $T_{\mu\nu}^{\rm matter}$ seule pour différents $\sigma$ |
| `total_conservation_v2_sigma_summary.csv` | Résumé de la compensation matière-intrication : réduction de divergence pour $\sigma$ et coefficients optimaux |

---

## Description détaillée des colonnes

### `paper8_summary.csv`

| Colonne | Signification |
|---------|---------------|
| `source_type` | Type de source testée ($T_{00}$, fluide, avec gradient, avec flux) |
| `best_sigma` | Largeur optimale de l’excitation |
| `spearman_abs` | Corrélation de Spearman avec $\vert\delta R\vert$ |
| `p_value` | Significativité statistique |
| `interpretation` | Brève interprétation |

### `fluid_source_scan.csv`

| Colonne | Signification |
|---------|---------------|
| `sigma` | Largeur de l’excitation |
| `omega` | Coefficient de $T_{aa}$ (pression isotrope) |
| `spearman_abs` | Corrélation avec $\vert\delta R\vert$ |
| `p_value` | Significativité |

### `total_conservation_v2_sigma_summary.csv`

| Colonne | Signification |
|---------|---------------|
| `sigma` | Largeur de l’excitation |
| `div_matter` | Norme L2 de la divergence de $T_{\mu\nu}^{\rm matter}$ seule |
| `div_total` | Norme L2 de la divergence après ajout du flux d’intrication |
| `reduction_factor` | Rapport `div_matter / div_total` |
| `beta` | Coefficient de $\delta R$ dans $\phi_{\rm ent}$ |
| `gamma` | Coefficient de $\Delta\delta R$ dans $\phi_{\rm ent}$ |
| `optimal` | Indique si c’est le meilleur point pour ce $\sigma$ |

---

## Comment reproduire ces résultats

Les scripts correspondants se trouvent dans `../scripts/` :

- `bup_paper8_Tmatter_positivity_v1.py` → vérification de la positivité de $T_{00}$
- `bup_paper8_Tmatter_localization_v1.py` → scan de la localisation ($R_{\rm loc}$ vs $\sigma$)
- `bup_paper8_Tmatter_curvature_response_v1.py` → corrélation $T_{00}$ vs $\vert\delta R\vert$
- `bup_paper8_Tmatter_fluid_source_v1.py` → scan de la source fluide ($\omega$)
- `bup_paper8_Tmatter_grad_source_v1.py` → scan du gradient interne ($\chi$)
- `bup_paper8_Tmatter_flux_source_v1.py` → scan du flux ($\psi$)
- `bup_paper8_Tmatter_conservation_v1.py` → divergence de $T_{\mu\nu}^{\rm matter}$ seule
- `bup_paper8_total_conservation_test_v2.py` → test de conservation totale avec flux d’intrication

Exemple de commande pour lancer le test de conservation totale :

```bash
python3 ../scripts/bup_paper8_total_conservation_test_v2.py \
  --mi-file ../data/MI_N20_lam0.57.csv \
  --k 5 \
  --sigma 0.01 \
  --beta -0.5 \
  --gamma -1.0 \
  --output total_conservation_v2_sigma_0.01.csv

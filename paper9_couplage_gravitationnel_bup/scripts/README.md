# Scripts — BuP Paper 9

Ce dossier contient le script utilisé dans le **BuP Paper 9 — Couplage gravitationnel de la matière émergente**.

L’objectif est de tester le passage physique suivant :

$$
\delta W^{loc} \longrightarrow S_{flux} \longrightarrow \Phi_{BuP} \longrightarrow |\delta R|.
$$

Paper 8 a construit la source matière candidate \( S_{flux} \). Paper 9 teste si cette source peut générer un potentiel gravitationnel effectif \( \Phi_{BuP} \) via une équation de Poisson discrète sur le graphe d’intrication.

---

## Script principal

### `bup_paper9_poisson_effective_potential_v1.py`

Ce script :

1. charge une matrice d’information mutuelle \( W_{ij} = I(i:j) \) ;
2. construit le graphe d’intrication \( k \)-NN ;
3. injecte une excitation locale \( \delta W^{loc} \) ;
4. reconstruit la source matière candidate de Paper 8 :

$$
S_{flux} = T_{00} - \frac{1}{2} T_{aa} + \frac{1}{2} T_{grad} + \|T_{0a}\| ;
$$

5. résout l’équation de Poisson discrète :

$$
L_{ent} \Phi_{BuP} = S_{flux} ;
$$

6. calcule le gradient effectif :

$$
|\nabla \Phi_{BuP}| ;
$$

7. compare \( S_{flux} \), \( \Phi_{BuP} \) et \( |\nabla \Phi_{BuP}| \) à la réponse de courbure nodale :

$$
|\delta R_i| = \left| R_i[W^{(0)} + \delta W] - R_i[W^{(0)}] \right|.
$$

---

## Commande principale

Run principal utilisé dans le papier :

```bash
python3 papers/paper9_couplage_gravitationnel_bup/scripts/bup_paper9_poisson_effective_potential_v1.py \
  --mi-file results_mi_N20_full/MI_N20_lam0.57.csv \
  --k 5 \
  --amp 0.15 \
  --sigma 0.15 \
  --omega -0.5 \
  --chi 0.5 \
  --psi 1.0 \
  --normalize-components \
  --output-dir papers/paper9_couplage_gravitationnel_bup/results/run_N20_lam057_sigma_0.15

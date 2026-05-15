# Results — BuP Paper 9

Ce dossier contient les résultats numériques réduits du Paper 9.

---

## Fichiers

| Fichier | Description |
|---|---|
| `paper9_sigma_scan_summary.csv` | Résumé du scan en largeur \( \sigma \) |
| `paper9_poisson_summary.csv` | Résumé du run principal \( \sigma = 0.15 \) |
| `summary.json` | Configuration et métriques détaillées du run principal |

---

## Résultat principal

Le run principal utilise :

$$
N = 20,\qquad \lambda = 0.57,\qquad k = 5,\qquad \sigma = 0.15.
$$

La source de Paper 8 est :

$$
S_{flux} = T_{00} - \frac12 T_{aa} + \frac12 T_{grad} + \|T_{0a}\|.
$$

Elle vérifie :

$$
\rho_{Spearman}(S_{flux}, |\delta R|) = 0.741,\qquad p = 1.84 \times 10^{-4}.
$$

Après résolution :

$$
L_{ent} \Phi_{BuP} = S_{flux},
$$

le potentiel obtenu vérifie :

$$
\rho_{Spearman}(\Phi_{BuP}, |\delta R|) = 0.738,\qquad p = 2.01 \times 10^{-4}.
$$

---

## Scan en \( \sigma \)

Le scan montre que le meilleur couplage gravitationnel est obtenu pour :

$$
\sigma = 0.15.
$$

Les excitations compactes \( \sigma = 0.01 \) et \( \sigma = 0.02 \) donnent un signal modéré, tandis que les valeurs intermédiaires \( \sigma = 0.05, 0.08, 0.10 \) ne donnent pas de corrélation significative.

---

## Interprétation

Ces résultats suggèrent que la source tensorielle émergente de Paper 8 peut générer un potentiel gravitationnel effectif cohérent avec la réponse de courbure. Cela constitue le premier test physique du couplage entre matière émergente et géométrie effective dans BuP.

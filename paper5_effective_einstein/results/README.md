# Résultats numériques — Paper 5

Ce dossier contient les résultats compacts des tests internes de cohérence au niveau du graphe, entre les proxys de courbure émergente et les proxys de contrainte spectrale (tenseur effectif).

Pour chaque session de calcul, seules les sorties compactes et reproductibles sont conservées :

- `summary.json` — configuration de la simulation et meilleures combinaisons de proxys (globales / stables)
- `stability_summary.csv` — classement de stabilité des proxys sur l’ensemble des matrices MI
- `scan_global.csv` — scan global des proxys
- `scan_by_N.csv` — scan des proxys regroupé par taille du graphe (quand disponible)
- `fig_*.png` — synthèses visuelles sélectionnées

Les grandes tables de diagnostics au niveau des arêtes ou des nœuds ne sont pas versionnées par défaut. Elles peuvent être régénérées avec la commande :

```bash
python3 scripts/bup_einstein_tensor_correlation_v2.py \
  --mi-files "<chemin-vers-fichiers-MI>/*.csv" \
  --k 5 \
  --tau-min 0.01 \
  --tau-max 50 \
  --output-dir results/<nom_de_la_session>

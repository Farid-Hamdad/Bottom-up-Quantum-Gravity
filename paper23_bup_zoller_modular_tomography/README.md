# Article 23 — Tomographie modulaire BuP/Zoller

Ce dossier contient l'article scientifique établissant le lien entre le principe tomographique BuP et les données d'intrication d'ions piégés de Zoller/Joshi (2023).

Résultat principal :

```math
\rho_{\rm ent}^{\rm cut}(j)
=
\sum_{k\notin A} I(j:k)
\propto
T_{\rm mod}(j)
=
\frac{1}{\beta_j}.

Contenu
text
main.tex
main.pdf
scripts/
results/
  modular_profiles/
  raw_mi/
  fig1_fig3/
figures/

Scripts
analyze_zoller_bup.py : extrait les profils modulaires et teste la relation 1/beta_j.

raw_pairwise_mi_test.py : reconstruit l'information mutuelle par paires à partir des données Pauli brutes.

shadow_mi_test.py : compare les densités d'information mutuelle globale, interne et tronquée.

analyze_fig1_fig3_bup.py : analyse la mise à l'échelle de l'entropie, la décroissance de l'IM et les diagnostics de sous-systèmes disjoints.

Principales affirmations numériques
États fondamentaux : profils modulaires paraboliques compatibles avec la forme CFT j(L-j).

États excités : entropie en loi de volume avec ajustements linéaires R^2 ≈ 0,995.

La densité globale naïve sum_k I(j:k) ne correspond pas à 1/beta_j.

La densité tronquée sum_{k notin A} I(j:k) est mieux corrélée à 1/beta_j.

Les liens entre sous-systèmes disjoints améliorent la fidélité de la reconstruction.

Note sur les données
Les fichiers .mat bruts de Joshi et al. ne sont pas inclus par défaut. Conservez les résumés CSV/JSON dérivés dans results/ et documentez la source des données originales dans l'article.

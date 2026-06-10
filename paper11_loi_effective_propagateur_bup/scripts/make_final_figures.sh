#!/usr/bin/env bash
set -e

BASE="papers/paper11_loi_effective_propagateur_bup"
RUN="$BASE/results/dw_scaling_N3000_k24_seed0_recheck"
FIG="$BASE/figures"

mkdir -p "$FIG"

cp -av "$RUN/figures/fig_ds_vs_N.png" "$FIG/"
cp -av "$RUN/figures/fig_dw_vs_N.png" "$FIG/"
cp -av "$RUN/figures/fig_alpha_predictions_vs_N.png" "$FIG/"

echo "Figures finales copiées dans:"
echo "$FIG"
ls -lh "$FIG"

#!/usr/bin/env bash
set -e

BASE="papers/paper11_loi_effective_propagateur_bup"
RUN="$BASE/results/dw_scaling_N3000_k24_seed0_recheck"
FINAL="$BASE/results/final"
FIG="$BASE/figures"

echo "============================================================"
echo "Paper 11 — Loi effective du propagateur BuP"
echo "============================================================"

mkdir -p "$FINAL" "$FIG"

echo "[1] Copy final numerical results"

cp -av "$RUN/dw_scaling_summary.csv" "$FINAL/"
cp -av "$RUN/spectral_dimension_all.csv" "$FINAL/"
cp -av "$RUN/walk_msd_all.csv" "$FINAL/"
cp -av "$RUN/summary.json" "$FINAL/"

echo "[2] Copy final figures"

cp -av "$RUN/figures/fig_ds_vs_N.png" "$FIG/"
cp -av "$RUN/figures/fig_dw_vs_N.png" "$FIG/"
cp -av "$RUN/figures/fig_alpha_predictions_vs_N.png" "$FIG/"

echo "[3] Done"

echo
echo "Final results:"
ls -lh "$FINAL"

echo
echo "Final figures:"
ls -lh "$FIG"

echo "============================================================"
echo "Paper 11 GitHub package ready."
echo "============================================================"

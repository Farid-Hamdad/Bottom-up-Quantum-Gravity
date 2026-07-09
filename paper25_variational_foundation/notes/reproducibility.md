# Paper 25 — Reproducibility Notes

## Repository branch

Paper 25 is built on the branch:

```text
bup_cosmology
Farid-Hamdad/Bottom-up-Quantum-Gravity`
paper25_variational_foundation/
  README.md
  paper25_variational_foundation.tex

  scripts/
    paper25_build_foundation_tables_v1.py
    paper25_make_fig01_variational_chain_v1.py
    paper25_make_fig02_dictionary_v1.py
    paper25_make_fig03_correction_hierarchy_v1.py
    paper25_make_fig04_cross_scale_fixed_points_v1.py

  results/
    paper25_foundation_summary_v1/
      paper25_chain_table.csv
      paper25_canonical_dictionary.csv
      paper25_evidence_status_table.csv
      paper25_open_debts_table.csv
      paper25_summary.json

  figures/
    fig01_variational_chain.png
    fig01_variational_chain.pdf
    fig02_discrete_to_continuum_dictionary.png
    fig02_discrete_to_continuum_dictionary.pdf
    fig03_correction_tensor_hierarchy.png
    fig03_correction_tensor_hierarchy.pdf
    fig04_cross_scale_fixed_points.png
    fig04_cross_scale_fixed_points.pdf

  notes/
    roadmap.md
    open_problems.md
    referee_notes.md
    reproducibility.md
python paper25_variational_foundation/scripts/paper25_make_fig01_variational_chain_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig02_dictionary_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig03_correction_hierarchy_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig04_cross_scale_fixed_points_v1.py
Paper 25 folder

Expected structure:

paper25_variational_foundation/
  README.md
  paper25_variational_foundation.tex

  scripts/
    paper25_build_foundation_tables_v1.py
    paper25_make_fig01_variational_chain_v1.py
    paper25_make_fig02_dictionary_v1.py
    paper25_make_fig03_correction_hierarchy_v1.py
    paper25_make_fig04_cross_scale_fixed_points_v1.py

  results/
    paper25_foundation_summary_v1/
      paper25_chain_table.csv
      paper25_canonical_dictionary.csv
      paper25_evidence_status_table.csv
      paper25_open_debts_table.csv
      paper25_summary.json

  figures/
    fig01_variational_chain.png
    fig01_variational_chain.pdf
    fig02_discrete_to_continuum_dictionary.png
    fig02_discrete_to_continuum_dictionary.pdf
    fig03_correction_tensor_hierarchy.png
    fig03_correction_tensor_hierarchy.pdf
    fig04_cross_scale_fixed_points.png
    fig04_cross_scale_fixed_points.pdf

  notes/
    roadmap.md
    open_problems.md
    referee_notes.md
    reproducibility.md
Generate summary tables

From the repository root:

python paper25_variational_foundation/scripts/paper25_build_foundation_tables_v1.py
Generate figures

From the repository root:

python paper25_variational_foundation/scripts/paper25_make_fig01_variational_chain_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig02_dictionary_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig03_correction_hierarchy_v1.py
python paper25_variational_foundation/scripts/paper25_make_fig04_cross_scale_fixed_points_v1.py
Compile LaTeX

From the Paper 25 folder:

cd paper25_variational_foundation
pdflatex paper25_variational_foundation.tex
pdflatex paper25_variational_foundation.tex

Expected output:

paper25_variational_foundation.pdf
Dependencies

Required Python package:

matplotlib

Required LaTeX packages:

amsmath
amssymb
amsfonts
mathtools
graphicx
booktabs
longtable
array
hyperref
xcolor
Source papers synthesized

Paper 25 synthesizes:

Papers 2--4     cosmological dimensional sector
Papers 8--11    emergent source and weak-field propagator
Papers 12--14   SPARC and micro-closure
Papers 15--19   Einstein-limit foundation
Paper 20        correction tensor
Paper 21        SLACS fixed point
Paper 22        gravitational-wave sector
Paper 23        modular tomography
Reproducibility status

Paper 25 is a synthesis and foundation paper.

It does not introduce new numerical experiments, but organizes existing results into:

canonical chain table,
dictionary table,
evidence-status table,
open-debts table,
conceptual figures.

All generated artifacts are deterministic.

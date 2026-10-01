#!/bin/sh
set -eu
python scripts/paper28_reproduce_modular_time_v1.py
python scripts/paper28_verify_reproduction_v1.py
python scripts/paper28_make_figures_v1.py

#!/usr/bin/env python3
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results'/'paper28_modular_time_v1'

def rows(path):
    with path.open(newline='',encoding='utf-8') as f: return {r['arm']:r for r in csv.DictReader(f)}
ref=rows(RES/'frozen_reference_metrics.csv')
pri=rows(RES/'primary_metrics.csv'); sec=rows(RES/'secondary_metrics.csv')
checks=[]
for arm,r in ref.items():
    checks += [
        (arm,'primary_Delta_R',float(pri[arm]['Delta_R_primary']),float(r['primary_Delta_R'])),
        (arm,'secondary_mean_Delta_R',float(sec[arm]['mean_Delta_R']),float(r['secondary_mean_Delta_R'])),
        (arm,'secondary_fraction_W_lower',float(sec[arm]['fraction_W_lower_residual']),float(r['secondary_fraction_W_lower'])),
        (arm,'secondary_min_Delta_R',float(sec[arm]['min_Delta_R']),float(r['secondary_min_Delta_R'])),
        (arm,'secondary_max_Delta_R',float(sec[arm]['max_Delta_R']),float(r['secondary_max_Delta_R'])),
        (arm,'order_transition_count',float(sec[arm]['order_transition_count']),float(r['order_transition_count'])),
        (arm,'observed_order_count',float(sec[arm]['observed_order_count']),float(r['observed_order_count'])),
    ]
maxerr=max(abs(a-b) for _,_,a,b in checks)
for arm,key,a,b in checks:
    tol=0.0 if key in ('order_transition_count','observed_order_count') else 1e-10
    if abs(a-b)>tol: raise SystemExit(f'REPRODUCTION_MISMATCH {arm} {key}: got {a}, ref {b}, abs_err {abs(a-b)}')
print('max_abs_error =',repr(maxerr))
print('PAPER28_REFERENCE_REPRODUCTION_PASS')

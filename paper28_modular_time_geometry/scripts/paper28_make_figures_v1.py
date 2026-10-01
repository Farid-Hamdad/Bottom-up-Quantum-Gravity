#!/usr/bin/env python3
import csv
from pathlib import Path
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results'/'paper28_modular_time_v1'
FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)

def read(path):
    with path.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
pri=read(RES/'primary_metrics.csv'); sec=read(RES/'secondary_metrics.csv')
labels=[r['arm'].replace('_',' / ') for r in pri]
fig,ax=plt.subplots(figsize=(7.4,4.5)); vals=[float(r['Delta_R_primary']) for r in pri]
ax.bar(labels,vals); ax.axhline(0,linewidth=1); ax.set_ylabel(r'$\Delta_R^{\rm flow}=R_{chain}-R_W$'); ax.set_title('Primary short-time modular-flow contrast'); ax.tick_params(axis='x',rotation=18); fig.tight_layout(); fig.savefig(FIG/'fig01_primary_deltaR.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(7.4,4.5)); vals=[float(r['mean_Delta_R']) for r in sec]
ax.bar(labels,vals); ax.axhline(0,linewidth=1); ax.set_ylabel(r'$\langle\Delta_R(\tau)\rangle_{\tau>0}$'); ax.set_title('Frozen-grid secondary modular-flow contrast'); ax.tick_params(axis='x',rotation=18); fig.tight_layout(); fig.savefig(FIG/'fig02_secondary_mean_deltaR.png',dpi=180); plt.close(fig)
arms=['N8_A4','N10_A4']; d={r['arm']:r for r in pri}; s={r['arm']:r for r in sec}; x=[0,1]; width=.34
fig,ax=plt.subplots(figsize=(7.0,4.5)); ax.bar([i-width/2 for i in x],[float(d[a]['Delta_R_primary']) for a in arms],width,label='Primary'); ax.bar([i+width/2 for i in x],[float(s[a]['mean_Delta_R']) for a in arms],width,label='Secondary mean'); ax.set_xticks(x,['N=8, |A|=4','N=10, |A|=4']); ax.set_ylabel(r'$\Delta_R$'); ax.set_title('Fixed-subsystem-size control'); ax.legend(); fig.tight_layout(); fig.savefig(FIG/'fig03_fixed_A4_control.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(7.4,4.5)); vals=[float(r['fraction_W_lower_residual']) for r in sec]; ax.bar(labels,vals); ax.set_ylim(0.96,1.002); ax.set_ylabel('Fraction with lower residual for W'); ax.set_title('Fraction of frozen modular times favoring W'); ax.tick_params(axis='x',rotation=18); fig.tight_layout(); fig.savefig(FIG/'fig04_W_win_fraction.png',dpi=180); plt.close(fig)
print('PAPER28_FIGURES_PASS')

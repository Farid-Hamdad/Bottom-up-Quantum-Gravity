#!/usr/bin/env python3
import argparse, os, pandas as pd, numpy as np, matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-csv', required=True)
    ap.add_argument('--output-dir', default='results')
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.raw_csv)

    summary = df.groupby('lambda').agg(
        d_port_mean=('d_port','mean'),
        d_port_std=('d_port','std')
    ).reset_index()
    summary.to_csv(os.path.join(args.output_dir,'d_lambda_summary.csv'), index=False)

    overall = df['d_port'].mean()

    plt.figure(figsize=(6,4))
    plt.errorbar(summary['lambda'], summary['d_port_mean'], yerr=summary['d_port_std'], fmt='o')
    plt.axhline(overall, ls='--', label=f'd ~ {overall:.2f}')
    plt.xlabel('lambda'); plt.ylabel('d_port')
    plt.legend(); plt.grid(alpha=0.3)
    out = os.path.join(args.output_dir,'fig_micro_coherence.svg')
    plt.savefig(out)
    print('saved', out)

if __name__=='__main__':
    main()

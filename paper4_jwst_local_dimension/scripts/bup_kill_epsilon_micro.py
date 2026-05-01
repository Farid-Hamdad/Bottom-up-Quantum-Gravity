#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bup_kill_epsilon_micro.py — BuP Paper 4"""
import argparse, glob, json, os, warnings
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
from scipy.stats import linregress
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

EPS_PAPER4 = 0.0113

def load_matrix_csv(path):
    arr = np.genfromtxt(path, delimiter=",", dtype=float)
    if arr.ndim==2 and arr.shape[0]==arr.shape[1]: return clean_matrix(arr)
    raise ValueError(f"Pas une matrice carrée : {path}")

def clean_matrix(W):
    W=np.array(W,dtype=float); W[~np.isfinite(W)]=0.0
    W=np.maximum(W,0.0); W=0.5*(W+W.T); np.fill_diagonal(W,0.0)
    mx=np.max(W);
    if mx>0: W=W/mx
    return W

def laplacian(W):
    return np.diag(W.sum(1))-W

def spectral_dimension(W, tau_min=0.05, tau_max=20.0, n_tau=80, min_points=8):
    n=W.shape[0]
    if n<4 or np.max(W)<=0: return np.nan,np.nan
    L=laplacian(W); evals=np.maximum(eigvalsh(L),0.0)
    tau=np.logspace(np.log10(tau_min),np.log10(tau_max),n_tau)
    P=np.array([np.mean(np.exp(-t*evals)) for t in tau])
    mask=np.isfinite(P)&(P>1e-12)&(P<0.98)
    if mask.sum()<min_points: mask=np.isfinite(P)&(P>1e-12)
    if mask.sum()<min_points: return np.nan,np.nan
    res=linregress(np.log(tau[mask]),np.log(P[mask]))
    return float(-2.0*res.slope),float(res.rvalue**2)

def local_subgraph(W, i, k=5):
    n=W.shape[0]; k=min(k,n-1)
    neigh=np.argsort(W[i])[::-1][:k]
    idx=np.unique(np.concatenate([[i],neigh]))
    return W[np.ix_(idx,idx)],idx

def analyze_matrix(W, label, k=5, tau_min=0.05, tau_max=20.0):
    n=W.shape[0]; strength=W.sum(1)
    mean_s=float(np.mean(strength)); std_s=float(np.std(strength)+1e-12)
    d_bg,r2_bg=spectral_dimension(W,tau_min=tau_min,tau_max=tau_max)
    rows=[]
    for i in range(n):
        Wloc,idx=local_subgraph(W,i,k=k)
        dloc,r2loc=spectral_dimension(Wloc,tau_min=tau_min,tau_max=tau_max)
        rows.append({
            "label":label,"N":int(n),"node":int(i),"k":int(k),
            "strength":float(strength[i]),
            "delta_rel":float(strength[i]/mean_s-1.0),
            "delta_z":float((strength[i]-mean_s)/std_s),
            "d_bg":float(d_bg),"d_local":float(dloc) if np.isfinite(dloc) else np.nan,
            "delta_d":float(d_bg-dloc) if np.isfinite(dloc) and np.isfinite(d_bg) else np.nan,
            "r2_bg":float(r2_bg),"r2_local":float(r2loc) if np.isfinite(r2loc) else np.nan,
            "subgraph_size":int(Wloc.shape[0]),
        })
    return rows

def fit_epsilon(rows, delta_key="delta_rel"):
    xs=[r[delta_key] for r in rows if np.isfinite(r.get(delta_key,np.nan)) and np.isfinite(r.get("d_local",np.nan))]
    ys=[r["d_local"]  for r in rows if np.isfinite(r.get(delta_key,np.nan)) and np.isfinite(r.get("d_local",np.nan))]
    xs=np.array(xs,float); ys=np.array(ys,float)
    if len(xs)<4 or np.std(xs)<1e-12:
        return {"epsilon":np.nan,"intercept":np.nan,"slope":np.nan,"r2":np.nan,"pvalue":np.nan,"stderr":np.nan,"n":int(len(xs))}
    res=linregress(xs,ys)
    return {"epsilon":float(-res.slope),"intercept":float(res.intercept),"slope":float(res.slope),
            "r2":float(res.rvalue**2),"pvalue":float(res.pvalue),"stderr":float(res.stderr),"n":int(len(xs))}

def bootstrap_epsilon(rows, delta_key="delta_rel", n_boot=2000, seed=42):
    rng=np.random.default_rng(seed)
    valid=[r for r in rows if np.isfinite(r.get(delta_key,np.nan)) and np.isfinite(r.get("d_local",np.nan))]
    if len(valid)<5: return [np.nan,np.nan,np.nan]
    eps=[]
    for _ in range(n_boot):
        s=[valid[j] for j in rng.integers(0,len(valid),len(valid))]
        f=fit_epsilon(s,delta_key)
        if np.isfinite(f["epsilon"]): eps.append(f["epsilon"])
    if not eps: return [np.nan,np.nan,np.nan]
    eps=np.array(eps)
    return [float(np.percentile(eps,16)),float(np.median(eps)),float(np.percentile(eps,84))]

def make_figure(rows, output_dir):
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    fig.subplots_adjust(wspace=0.35)
    colors={'N09':'#1a3a7a','N12':'#7a1a3a','N16':'#1a7a3a'}
    for ax,key,title in [
        (axes[0],"delta_rel",r"$\delta_{\rm rel}=s_i/\langle s\rangle-1$"),
        (axes[1],"delta_z",  r"$\delta_z=(s_i-\langle s\rangle)/\sigma_s$"),
    ]:
        # Par groupe N
        for tag,col in colors.items():
            sub=[r for r in rows if tag in r['label'] and np.isfinite(r.get(key,np.nan)) and np.isfinite(r.get('d_local',np.nan))]
            if not sub: continue
            xs=np.array([r[key] for r in sub]); ys=np.array([r['d_local'] for r in sub])
            N_val=sub[0]['N']
            ax.scatter(xs,ys,s=25,alpha=0.5,color=col,label=f'N={N_val}')
        # Fit global
        fit=fit_epsilon(rows,delta_key=key)
        valid=[r for r in rows if np.isfinite(r.get(key,np.nan)) and np.isfinite(r.get('d_local',np.nan))]
        if valid and np.isfinite(fit['epsilon']):
            xs_all=np.array([r[key] for r in valid])
            xx=np.linspace(xs_all.min(),xs_all.max(),100)
            yy=fit['intercept']+fit['slope']*xx
            ax.plot(xx,yy,'k-',lw=2,label=fr"$\varepsilon={fit['epsilon']:.4f}$ (fit global)")
        ax.axvline(0,color='gray',ls=':',lw=1,alpha=0.5)
        ax.set_xlabel(title,fontsize=11); ax.set_ylabel(r"$d_{s,\rm local}$",fontsize=12)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax.set_title(f'ε = {fit["epsilon"]:.5f}  (R²={fit["r2"]:.3f})',fontsize=10)
    fig.suptitle("BuP Paper 4 — Test microscopique : d_local vs connectivité locale",fontsize=11,y=1.01)
    path=os.path.join(output_dir,"fig_epsilon_micro.png")
    plt.savefig(path,dpi=180,bbox_inches='tight'); plt.close()
    return path

def print_verdict(name, fit, ci):
    lo,med,hi=ci; eps=fit['epsilon']
    print(f"\n{name}")
    print("─"*65)
    print(f"  ε (fit)    = {eps:.6f}")
    print(f"  bootstrap  = [{lo:.6f},  {med:.6f},  {hi:.6f}]")
    print(f"  R²         = {fit['r2']:.4f}")
    print(f"  p-value    = {fit['pvalue']:.3g}")
    print(f"  N points   = {fit['n']}")
    if not np.isfinite(eps):
        verdict="impossible à estimer"
    elif np.isfinite(lo) and hi<0:
        verdict="signe opposé à Paper 4 — hypothèse fragilisée"
    elif np.isfinite(lo) and lo<=0<=hi:
        verdict="compatible avec zéro — ε non détecté à cette résolution"
    else:
        verdict=f"ε positif détecté ✓  (distance Paper4={abs(eps-EPS_PAPER4):.5f})"
    print(f"  VERDICT    = {verdict}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mi-files',nargs='*',default=[])
    ap.add_argument('--simulate',action='store_true')
    ap.add_argument('--n-sim',type=int,default=30)
    ap.add_argument('--N',type=int,default=16)
    ap.add_argument('--k',type=int,default=5)
    ap.add_argument('--tau-min',type=float,default=0.05)
    ap.add_argument('--tau-max',type=float,default=20.0)
    ap.add_argument('--output-dir',default='results_epsilon_micro')
    args=ap.parse_args()
    os.makedirs(args.output_dir,exist_ok=True)

    print("="*65)
    print("BuP Paper 4 — Test microscopique de ε")
    print(f"k={args.k}  τ∈[{args.tau_min},{args.tau_max}]")
    print("="*65)

    all_rows=[]
    files=[]
    for pat in args.mi_files: files.extend(glob.glob(pat))

    if files:
        print(f"\nAnalyse de {len(files)} matrice(s) MI BuP")
        for path in sorted(files):
            try:
                W=load_matrix_csv(path)
                rows=analyze_matrix(W,label=Path(path).stem,k=args.k,
                                    tau_min=args.tau_min,tau_max=args.tau_max)
                all_rows.extend(rows)
                d_bg=rows[0]['d_bg']
                s_std=np.std([r['strength'] for r in rows])
                print(f"  ✓ {Path(path).name}  N={W.shape[0]}  d_bg={d_bg:.4f}  s_std={s_std:.4f}")
            except Exception as e:
                print(f"  ✗ {path}: {e}")

    if args.simulate or not files:
        from scipy.linalg import eigvalsh as _ev
        print(f"\nMode simulation : {args.n_sim} matrices synthétiques N={args.N}")
        rng=np.random.default_rng(42)
        for s in range(args.n_sim):
            n=args.N
            x=np.linspace(0,1,n); dist=np.abs(x[:,None]-x[None,:])
            W=np.exp(-dist/0.25)*np.exp(0.15*rng.normal(size=(n,n)))
            W=0.5*(W+W.T); np.fill_diagonal(W,0)
            start=rng.integers(0,max(1,n-4))
            idx=np.arange(start,min(start+4,n))
            for a in idx:
                for b in idx:
                    if a!=b: W[a,b]*=1.5+rng.random()*0.5
            W=clean_matrix(W)
            rows=analyze_matrix(W,label=f"sim_{s:04d}",k=args.k,
                                tau_min=args.tau_min,tau_max=args.tau_max)
            all_rows.extend(rows)

    if not all_rows: raise RuntimeError("Aucune donnée.")

    # Sauvegardes
    csv_path=os.path.join(args.output_dir,"node_table.csv")
    keys=list(all_rows[0].keys())
    with open(csv_path,'w') as f:
        f.write(",".join(keys)+"\n")
        for r in all_rows:
            f.write(",".join(str(r.get(k,'')) for k in keys)+"\n")

    # Fits et bootstrap
    fit_rel=fit_epsilon(all_rows,"delta_rel")
    fit_z  =fit_epsilon(all_rows,"delta_z")
    ci_rel =bootstrap_epsilon(all_rows,"delta_rel")
    ci_z   =bootstrap_epsilon(all_rows,"delta_z")

    # Résumé JSON
    summary={"n_rows":len(all_rows),"k":args.k,"epsilon_paper4":EPS_PAPER4,
             "fit_delta_rel":fit_rel,"ci_rel_16_50_84":ci_rel,
             "fit_delta_z":fit_z,"ci_z_16_50_84":ci_z}
    with open(os.path.join(args.output_dir,"summary.json"),'w') as f:
        json.dump(summary,f,indent=2)

    fig_path=make_figure(all_rows,args.output_dir)

    print(f"\n{'='*65}")
    print("RÉSULTATS")
    print_verdict("Fit δ_rel = s_i/<s> − 1  (sur-connectivité relative)", fit_rel, ci_rel)
    print_verdict("Fit δ_z   = z-score connectivité", fit_z, ci_z)

    print(f"\n{'='*65}")
    print("INTERPRÉTATION BuP Paper 4")
    print("─"*65)
    eps_best=fit_rel['epsilon']
    if np.isfinite(eps_best):
        if 0.002 < eps_best < 0.05:
            print(f"  ε_micro = {eps_best:.5f}")
            print(f"  ε_JWST  = {EPS_PAPER4:.5f}")
            print(f"  Ratio   = {eps_best/EPS_PAPER4:.2f}x")
            print(f"  → Signal dans la plage attendue")
        elif eps_best <= 0:
            print(f"  ε ≤ 0 : pas de réduction dimensionnelle locale détectée")
            print(f"  → Limitation de taille système (N ≤ 16)")
            print(f"  → Nécessite N ≥ 30 pour signal robuste")
        else:
            print(f"  ε = {eps_best:.5f} (hors plage attendue)")
    print(f"\n  Fichiers : {csv_path}, {fig_path}")
    print("\nDONE")

if __name__=='__main__':
    main()

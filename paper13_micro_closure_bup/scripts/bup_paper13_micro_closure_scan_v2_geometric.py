#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuP Paper 13 — microscopic closure scan v2

Adds to v1:
- J_eff diagnostics: mean/sum/spectral coupling scales and ratios to h0.
- Control modes: none, shuffled_sigma, inverted_sigma, random_J.

Core test:
|Psi_gal> -> I_ij -> rho_ent^micro(R) ? proportional to Sigma(R)
"""
from __future__ import annotations

import argparse, json, time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import eigsh
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

EPS = 1e-12


def parse_args():
    p = argparse.ArgumentParser(description="BuP Paper 13 scan v2: J_eff + controls")
    p.add_argument("--n-rings", type=int, default=3)
    p.add_argument("--n-theta", type=int, default=4)
    p.add_argument("--r-max", type=float, default=10.0)
    p.add_argument("--rd", type=float, default=3.0)
    p.add_argument("--sigma0", type=float, default=1.0)
    p.add_argument("--xi-list", type=float, nargs="+", default=[0.8,1.2,1.6,2.0,2.5,3.0,4.0])
    p.add_argument("--h0-list", type=float, nargs="+", default=[0.1,0.2,0.3,0.5,0.8,1.0])
    p.add_argument("--j0", type=float, default=1.0)
    p.add_argument("--control-mode", choices=["none","shuffled_sigma","inverted_sigma","random_J","geometric_J"], default="none")
    p.add_argument("--seed", type=int, default=13)
    p.add_argument("--max-qubits", type=int, default=14)
    p.add_argument("--output-dir", type=str, default="results_paper13_scan_v2")
    p.add_argument("--save-all-radial", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


def make_disk(n_rings, n_theta, r_max, rd, sigma0, control_mode, seed):
    xs=[]; ys=[]; rs=[]; ths=[]; rings=[]; sig=[]
    dr = r_max / n_rings
    for k in range(n_rings):
        R = (k + 0.5) * dr
        for m in range(n_theta):
            th = 2*np.pi*m/n_theta
            xs.append(R*np.cos(th)); ys.append(R*np.sin(th)); rs.append(R); ths.append(th); rings.append(k)
            sig.append(sigma0*np.exp(-R/rd))
    x=np.array(xs,float); y=np.array(ys,float); r=np.array(rs,float); th=np.array(ths,float)
    ring=np.array(rings,int); sigma_target=np.array(sig,float)
    sigma_coupling=sigma_target.copy()
    rng=np.random.default_rng(seed)
    if control_mode == "shuffled_sigma":
        sigma_coupling = rng.permutation(sigma_target)
    elif control_mode == "inverted_sigma":
        order=np.argsort(r)
        sigma_coupling[order] = sigma_target[order][::-1]
    coords=np.c_[x,y]
    dist=np.sqrt(((coords[:,None,:]-coords[None,:,:])**2).sum(axis=-1))
    return dict(x=x,y=y,r=r,theta=th,ring_id=ring,sigma_target=sigma_target,
                sigma_coupling=sigma_coupling,dist=dist)


def build_J(disk, j0, xi, control_mode, seed):
    s=disk["sigma_coupling"]; dist=disk["dist"]
    J = j0*np.sqrt(np.outer(s,s))*np.exp(-dist/xi)
    np.fill_diagonal(J,0.0)

    if control_mode == "geometric_J":
        J = j0*np.exp(-dist/xi)
        np.fill_diagonal(J,0.0)

    if control_mode == "random_J":
        rng=np.random.default_rng(seed + int(round(1000*xi)))
        upper=J[np.triu_indices_from(J,k=1)]
        scale=float(np.mean(np.abs(upper))) if upper.size else abs(j0)
        scale=max(scale,EPS)
        R=rng.lognormal(mean=np.log(scale), sigma=0.5, size=J.shape)
        J=0.5*(R+R.T); np.fill_diagonal(J,0.0)
    return J


def coupling_diag(J, xi, h0):
    n=J.shape[0]
    upper=np.abs(J[np.triu_indices(n,k=1)])
    mean=float(upper.mean()) if upper.size else 0.0
    sumn=float(upper.sum()/upper.size) if upper.size else 0.0
    maxj=float(upper.max()) if upper.size else 0.0
    spec=float(np.max(np.abs(np.linalg.eigvalsh(J))))
    return dict(
        J_eff_mean=mean,
        J_eff_sum_norm=sumn,
        J_eff_max=maxj,
        J_eff_spectral=spec,
        J_eff_mean_over_h0=float(mean/h0) if abs(h0)>EPS else np.nan,
        J_eff_spectral_over_h0=float(spec/h0) if abs(h0)>EPS else np.nan,
        xi_over_h0=float(xi/h0) if abs(h0)>EPS else np.nan,
        h0_over_xi=float(h0/xi) if abs(xi)>EPS else np.nan,
    )


def paulis():
    X=csr_matrix(np.array([[0.,1.],[1.,0.]]))
    Z=csr_matrix(np.array([[1.,0.],[0.,-1.]]))
    I=identity(2,format="csr",dtype=float)
    return X,Z,I


def one_op(op, site, n):
    _,_,I=paulis(); out=None
    for k in range(n):
        f=op if k==site else I
        out=f if out is None else kron(out,f,format="csr")
    return out


def xx_op(i,j,n):
    X,_,I=paulis(); out=None
    for k in range(n):
        f=X if k in (i,j) else I
        out=f if out is None else kron(out,f,format="csr")
    return out


def precompute_ops(n):
    _,Z,_=paulis()
    z=[one_op(Z,i,n) for i in range(n)]
    xx={(i,j):xx_op(i,j,n) for i in range(n) for j in range(i+1,n)}
    return z,xx


def build_H(J,h0,z_ops,xx_ops):
    n=J.shape[0]; H=csr_matrix((2**n,2**n),dtype=float)
    for i in range(n): H += -h0*z_ops[i]
    for i in range(n):
        for j in range(i+1,n):
            if abs(J[i,j])>EPS: H += -J[i,j]*xx_ops[(i,j)]
    return H


def ground(H):
    vals,vecs=eigsh(H,k=1,which="SA")
    psi=np.asarray(vecs[:,0],complex); psi/=np.linalg.norm(psi)
    return float(vals[0]),psi


def rdm(psi, keep, n):
    keep=list(keep); trace=[i for i in range(n) if i not in keep]
    T=psi.reshape([2]*n).transpose(keep+trace)
    M=T.reshape(2**len(keep), 2**(n-len(keep)))
    rho=M@M.conj().T
    return 0.5*(rho+rho.conj().T)


def entropy(rho):
    w=np.linalg.eigvalsh(rho).real
    w=w[w>EPS]
    return float(-np.sum(w*np.log2(w))) if len(w) else 0.0


def MI_matrix(psi,n):
    S1=np.array([entropy(rdm(psi,[i],n)) for i in range(n)])
    MI=np.zeros((n,n),float)
    for i in range(n):
        for j in range(i+1,n):
            mij=S1[i]+S1[j]-entropy(rdm(psi,[i,j],n))
            if mij<0 and abs(mij)<1e-10: mij=0.0
            MI[i,j]=MI[j,i]=mij
    return MI,S1


def radial_df(values,disk,n_rings):
    rows=[]
    for k in range(n_rings):
        m=disk["ring_id"]==k
        rows.append(dict(ring_id=k,R_mean=float(disk["r"][m].mean()),n_cells=int(m.sum()),
                         sigma_target_mean=float(disk["sigma_target"][m].mean()),
                         sigma_coupling_mean=float(disk["sigma_coupling"][m].mean()),
                         rho_ent_mean=float(values[m].mean()),rho_ent_sum=float(values[m].sum())))
    return pd.DataFrame(rows)


def norm(y):
    y=np.asarray(y,float); s=y.sum()
    return np.zeros_like(y) if s<=EPS else y/s


def kl(p,q):
    p=norm(p)+EPS; q=norm(q)+EPS; p=p/p.sum(); q=q/q.sum()
    return float(np.sum(p*np.log(p/q)))


def exp_prof(R,A,Rscale): return A*np.exp(-R/Rscale)


def fit_scale(R,y,guess):
    y=np.maximum(np.asarray(y,float),EPS)
    try:
        popt,_=curve_fit(exp_prof,np.asarray(R,float),y,p0=(float(y.max()),float(guess)),
                         bounds=([0.,EPS],[np.inf,np.inf]),maxfev=10000)
        return float(popt[0]),float(popt[1])
    except Exception:
        return np.nan,np.nan


def diagnostics(rad,rd):
    sig=rad["sigma_target_mean"].to_numpy(float); rho=rad["rho_ent_mean"].to_numpy(float); R=rad["R_mean"].to_numpy(float)
    sn=norm(sig); rn=norm(rho)
    corr=float(pearsonr(sn,rn).statistic) if len(R)>=2 and sn.std()>EPS and rn.std()>EPS else np.nan
    rmse=float(np.sqrt(np.mean((rn-sn)**2)))
    _,Rent=fit_scale(R,rho,rd); _,Rs=fit_scale(R,sig,rd)
    rderr=float(abs(Rent-rd)/rd) if np.isfinite(Rent) else np.nan
    return dict(corr_rho_sigma=corr,rmse_normalized=rmse,kl_rho_ent_to_sigma=kl(rho,sig),
                kl_sigma_to_rho_ent=kl(sig,rho),R_ent_fit=Rent,R_sigma_fit=Rs,
                R_d_target=float(rd),R_d_relative_error=rderr)


def verdict(d):
    c=d["corr_rho_sigma"]; r=d["rmse_normalized"]; e=d["R_d_relative_error"]
    if np.isfinite(c) and np.isfinite(r) and np.isfinite(e):
        if c>=0.90 and r<=0.10 and e<=0.20: return "strong_closure"
        if c>=0.80 and r<=0.15 and e<=0.35: return "moderate_closure"
        return "weak_or_failed_closure"
    return "insufficient_diagnostics"


def pivot(df,metric): return df.pivot(index="h0",columns="xi",values=metric).sort_index().sort_index(axis=1)


def heatmap(df,metric,path,title):
    P=pivot(df,metric); xv=P.columns.to_numpy(float); yv=P.index.to_numpy(float); Z=P.to_numpy(float)
    plt.figure(figsize=(8,5.5)); im=plt.imshow(Z,origin="lower",aspect="auto",extent=[xv.min(),xv.max(),yv.min(),yv.max()])
    plt.colorbar(im,label=metric); plt.xlabel("xi"); plt.ylabel("h0"); plt.title(title); plt.tight_layout(); plt.savefig(path,dpi=200); plt.close()


def verdict_map(df,path):
    mp={"weak_or_failed_closure":0,"insufficient_diagnostics":0,"failed":0,"moderate_closure":1,"strong_closure":2}
    D=df.copy(); D["verdict_code"]=D["verdict"].map(mp).fillna(0)
    P=pivot(D,"verdict_code"); xv=P.columns.to_numpy(float); yv=P.index.to_numpy(float); Z=P.to_numpy(float)
    plt.figure(figsize=(8,5.5)); im=plt.imshow(Z,origin="lower",aspect="auto",extent=[xv.min(),xv.max(),yv.min(),yv.max()],vmin=0,vmax=2)
    c=plt.colorbar(im,ticks=[0,1,2]); c.ax.set_yticklabels(["weak/failed","moderate","strong"])
    plt.xlabel("xi"); plt.ylabel("h0"); plt.title("Microscopic closure verdict map"); plt.tight_layout(); plt.savefig(path,dpi=200); plt.close()


def scatter(df,xcol,ycol,path,title):
    plt.figure(figsize=(7,5.5))
    for v,label in [("weak_or_failed_closure","weak/failed"),("insufficient_diagnostics","insufficient"),("moderate_closure","moderate"),("strong_closure","strong")]:
        sub=df[df["verdict"]==v]
        if len(sub): plt.scatter(sub[xcol],sub[ycol],s=55,label=label)
    plt.xlabel(xcol); plt.ylabel(ycol); plt.title(title); plt.legend(); plt.tight_layout(); plt.savefig(path,dpi=200); plt.close()


def summarize_ratios(df):
    out={}
    for name,sub in [("strong",df[df.verdict=="strong_closure"]),("moderate_or_strong",df[df.verdict.isin(["moderate_closure","strong_closure"])])]:
        for col in ["xi_over_h0","h0_over_xi","J_eff_spectral_over_h0","J_eff_mean_over_h0"]:
            vals=pd.to_numeric(sub.get(col,pd.Series(dtype=float)),errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
            out[f"{name}_{col}_mean"]=float(vals.mean()) if len(vals) else None
            out[f"{name}_{col}_std"]=float(vals.std(ddof=0)) if len(vals) else None
            out[f"{name}_{col}_min"]=float(vals.min()) if len(vals) else None
            out[f"{name}_{col}_max"]=float(vals.max()) if len(vals) else None
    return out


def main():
    a=parse_args(); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    n=a.n_rings*a.n_theta
    if n>a.max_qubits: raise ValueError(f"N={n} > max-qubits={a.max_qubits}")
    disk=make_disk(a.n_rings,a.n_theta,a.r_max,a.rd,a.sigma0,a.control_mode,a.seed)
    pd.DataFrame(dict(site=np.arange(n),x=disk["x"],y=disk["y"],R=disk["r"],theta=disk["theta"],ring_id=disk["ring_id"],
                      sigma_target=disk["sigma_target"],sigma_coupling=disk["sigma_coupling"])).to_csv(out/"disk_sites.csv",index=False)
    print("\n=== BuP Paper 13 — microscopic closure scan v2 ===")
    print(f"N qubits/cells        : {n}")
    print(f"Hilbert dimension     : {2**n}")
    print(f"R_d target            : {a.rd:g}")
    print(f"j0                    : {a.j0:g}")
    print(f"control_mode          : {a.control_mode}")
    print(f"xi values             : {a.xi_list}")
    print(f"h0 values             : {a.h0_list}")
    print(f"Total runs            : {len(a.xi_list)*len(a.h0_list)}")
    print("\nPrecomputing operators...")
    z_ops,xx_ops=precompute_ops(n)
    radial_dir=out/"radial_profiles_by_point"
    if a.save_all_radial: radial_dir.mkdir(exist_ok=True)
    rows=[]; best=None; best_payload=None; t0=time.time(); idx=0; total=len(a.xi_list)*len(a.h0_list)
    for xi in a.xi_list:
        for h0 in a.h0_list:
            idx+=1; print(f"\n[{idx}/{total}] Running xi={xi:g}, h0={h0:g}...")
            try:
                J=build_J(disk,a.j0,xi,a.control_mode,a.seed); jd=coupling_diag(J,xi,h0)
                H=build_H(J,h0,z_ops,xx_ops); E,psi=ground(H); MI,S1=MI_matrix(psi,n)
                rho=MI.sum(axis=1); rad=radial_df(rho,disk,a.n_rings); dg=diagnostics(rad,a.rd); vd=verdict(dg)
                upper=MI[np.triu_indices(n,k=1)]
                row=dict(xi=float(xi),h0=float(h0),j0=float(a.j0),control_mode=a.control_mode,ground_state_energy=E,
                         mean_single_site_entropy=float(S1.mean()),max_single_site_entropy=float(S1.max()),
                         mean_mutual_information=float(upper.mean()),max_mutual_information=float(MI.max()),verdict=vd)
                row.update(dg); row.update(jd); rows.append(row)
                print(f"  corr={row['corr_rho_sigma']:.6f} | rmse={row['rmse_normalized']:.6f} | R_ent={row['R_ent_fit']:.6f} | rd_err={row['R_d_relative_error']:.6f} | Jspec/h0={row['J_eff_spectral_over_h0']:.6f} | xi/h0={row['xi_over_h0']:.6f} | {vd}")
                if a.save_all_radial: rad.to_csv(radial_dir/(f"radial_xi_{xi:g}_h0_{h0:g}.csv".replace('.','p')),index=False)
                if best is None or row["rmse_normalized"] < best["rmse_normalized"]:
                    best=row.copy(); best_payload=(rad.copy(),J.copy(),MI.copy())
            except Exception as e:
                print(f"  FAILED: {e}"); rows.append(dict(xi=float(xi),h0=float(h0),j0=float(a.j0),control_mode=a.control_mode,verdict="failed",error=str(e)))
    df=pd.DataFrame(rows); df.to_csv(out/"scan_results.csv",index=False)
    if best_payload:
        best_payload[0].to_csv(out/"best_radial_profile.csv",index=False)
        pd.DataFrame(best_payload[1]).to_csv(out/"best_coupling_matrix.csv",index=False)
        pd.DataFrame(best_payload[2]).to_csv(out/"best_mutual_information_matrix.csv",index=False)
    nstrong=int((df.verdict=="strong_closure").sum()); nmod=int((df.verdict=="moderate_closure").sum()); ntot=len(df)
    summary=dict(experiment="BuP Paper 13 microscopic closure scan v2",n_qubits=n,hilbert_dimension=2**n,n_rings=a.n_rings,n_theta=a.n_theta,
                 r_max=a.r_max,rd_target=a.rd,sigma0=a.sigma0,j0=a.j0,control_mode=a.control_mode,seed=a.seed,
                 xi_list=[float(x) for x in a.xi_list],h0_list=[float(x) for x in a.h0_list],n_total_runs=ntot,
                 n_successful_runs=int((df.verdict!="failed").sum()),n_strong_closure=nstrong,n_moderate_closure=nmod,
                 n_weak_failed_or_insufficient=int(df.verdict.isin(["weak_or_failed_closure","failed","insufficient_diagnostics"]).sum()),
                 fraction_strong_closure=float(nstrong/max(ntot,1)),fraction_moderate_or_strong=float((nstrong+nmod)/max(ntot,1)),
                 best_by_rmse={k:(None if isinstance(v,float) and not np.isfinite(v) else (float(v) if isinstance(v,(int,float,np.floating,np.integer)) else str(v))) for k,v in best.items()} if best else None,
                 ratio_summary=summarize_ratios(df),runtime_seconds=float(time.time()-t0))
    with open(out/"scan_summary.json","w") as f: json.dump(summary,f,indent=2)
    if not a.no_plots and len(df):
        D=df.replace([np.inf,-np.inf],np.nan)
        for metric,title in [("corr_rho_sigma",r"Closure correlation"),("rmse_normalized","Normalized RMSE"),("R_d_relative_error","Relative scale error"),("J_eff_spectral_over_h0",r"J_eff spectral / h0")]:
            try: heatmap(D,metric,out/f"fig_{metric}_xi_h0.png",title)
            except Exception as e: print(f"Could not save heatmap {metric}: {e}")
        try: verdict_map(D,out/"fig_verdict_map.png")
        except Exception as e: print(f"Could not save verdict map: {e}")
        try: scatter(D,"xi_over_h0","rmse_normalized",out/"fig_scatter_xi_over_h0_vs_rmse.png","RMSE versus xi/h0")
        except Exception as e: print(f"Could not save xi/h0 scatter: {e}")
        try: scatter(D,"J_eff_spectral_over_h0","rmse_normalized",out/"fig_scatter_jeff_over_h0_vs_rmse.png","RMSE versus J_eff_spectral/h0")
        except Exception as e: print(f"Could not save J_eff/h0 scatter: {e}")
    print("\n=== Scan summary ===")
    print(f"Control mode              : {a.control_mode}")
    print(f"Total runs                : {ntot}")
    print(f"Strong closure            : {nstrong}")
    print(f"Moderate closure          : {nmod}")
    print(f"Fraction strong           : {summary['fraction_strong_closure']:.3f}")
    print(f"Fraction moderate+strong  : {summary['fraction_moderate_or_strong']:.3f}")
    if best:
        print("\nBest point by RMSE:")
        for k in ["xi","h0","corr_rho_sigma","rmse_normalized","R_ent_fit","R_d_relative_error","xi_over_h0","J_eff_spectral_over_h0","verdict"]:
            print(f"  {k:24s}: {best[k]}")
    print("\nFiles written:")
    for name in ["scan_results.csv","scan_summary.json","disk_sites.csv","fig_verdict_map.png","fig_rmse_normalized_xi_h0.png","fig_J_eff_spectral_over_h0_xi_h0.png","fig_scatter_jeff_over_h0_vs_rmse.png"]:
        if (out/name).exists(): print(f"  - {out/name}")
    print("\nDone.")

if __name__ == "__main__":
    main()

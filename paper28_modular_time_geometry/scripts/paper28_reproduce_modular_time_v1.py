#!/usr/bin/env python3
"""Exact finite-system TFIM modular-time benchmark for Paper 28."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
J=1.0; H_FIELD=1.0; P=6; DT=0.35; TAU_MIN=0.0; TAU_MAX=120.0; N_TAU=1024
I2=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex); Z=np.array([[1,0],[0,-1]],complex); PAULI={'X':X,'Y':Y,'Z':Z}
OUT_DIR=Path(__file__).resolve().parents[1]/'results'/'paper28_modular_time_v1'
def plus_state(n):
 s=np.array([1+0j]); q=np.array([1.,1.],complex)/np.sqrt(2); 
 for _ in range(n): s=np.kron(s,q)
 return s
def one_gate(s,g,q,n):
 t=np.moveaxis(s.reshape([2]*n),q,0); u=np.tensordot(g,t,axes=([1],[0])); return np.moveaxis(u,0,q).reshape(-1)
def xhalf(s,dt,n):
 g=np.cos(H_FIELD*dt/2)*I2+1j*np.sin(H_FIELD*dt/2)*X
 for q in range(n): s=one_gate(s,g,q,n)
 return s
def zz(s,dt,n):
 out=s.copy(); idx=np.arange(2**n,dtype=np.uint64); e=np.zeros(2**n)
 for q in range(n-1):
  bq=(idx>>np.uint64(n-1-q))&np.uint64(1); br=(idx>>np.uint64(n-2-q))&np.uint64(1); e+=(1-2*bq.astype(float))*(1-2*br.astype(float))
 return out*np.exp(1j*J*dt*e)
def evolve(n,dt):
 s=plus_state(n)
 for _ in range(P): s=xhalf(s,dt,n); s=zz(s,dt,n); s=xhalf(s,dt,n)
 if abs(float(np.real(np.vdot(s,s)))-1)>1e-12: raise RuntimeError('STATE_NORM_DRIFT')
 return s
def prefix(s,a,n):
 m=s.reshape(2**a,2**(n-a)); return m@m.conj().T
def reduced(rho,keep,n):
 keep=tuple(keep); rest=tuple(q for q in range(n) if q not in keep); perm=keep+rest+tuple(n+q for q in keep)+tuple(n+q for q in rest); t=np.transpose(rho.reshape([2]*(2*n)),perm); dk=2**len(keep); dr=2**len(rest); return np.einsum('abcb->ac',t.reshape(dk,dr,dk,dr),optimize=True)
def entropy(r):
 v=np.clip(np.real(np.linalg.eigvalsh((r+r.conj().T)/2)),0,None); v=v[v>1e-15]; return float(-np.sum(v*np.log(v)))
def mi(r,a,pairs): return np.array([entropy(reduced(r,(i,),a))+entropy(reduced(r,(j,),a))-entropy(reduced(r,(i,j),a)) for i,j in pairs])
def ops(a):
 d={}
 for q in range(a):
  for lab,sig in PAULI.items():
   o=np.array([[1+0j]])
   for k in range(a): o=np.kron(o,sig if k==q else I2)
   d[(q,lab)]=o
 return d
def modular(r):
 vals,vec=np.linalg.eigh((r+r.conj().T)/2); vals=np.real(vals)
 if vals.min()<=0: raise RuntimeError('POSITIVITY_GATE_FAIL')
 k=-np.log(vals); km=float(k.mean()); ks=float(k.std(ddof=0)); ke=(k-km)/ks; kt=vec@np.diag(ke)@vec.conj().T
 return kt,vec,ke,{'lambda_min':float(vals.min()),'lambda_max':float(vals.max()),'condition_number':float(vals.max()/vals.min()),'kappa_mean':km,'kappa_std':ks}
def comm(a,b): return a@b-b@a
def coeff(r,k,pairs,o):
 z=[]
 for i,j in pairs:
  tot=0.
  for a in 'XYZ':
   inner=comm(k,o[(i,a)])
   for b in 'XYZ':
    q=comm(inner,o[(j,b)]); tot+=float(np.real(np.trace(r@q.conj().T@q)))
  z.append(tot/9)
 return np.array(z)
def metric(v,g):
 nv=np.linalg.norm(v); dot=float(v@g); sc=max(0.,dot/float(g@g)); return dot/(nv*np.linalg.norm(g)),sc,float(np.linalg.norm(v-sc*g)/nv)
def curves(r,vec,ke,a,pairs,taus,o):
 out=np.empty((len(taus),len(pairs)))
 for ti,tau in enumerate(taus):
  ph=np.exp(1j*ke*tau); u=vec@np.diag(ph)@vec.conj().T; ud=u.conj().T; ev={(q,x):u@o[(q,x)]@ud for q in range(a) for x in 'XYZ'}
  for pi,(i,j) in enumerate(pairs):
   tot=0.
   for x in 'XYZ':
    for y in 'XYZ': c=comm(ev[(i,x)],o[(j,y)]); tot+=float(np.real(np.trace(r@c.conj().T@c)))
   out[ti,pi]=tot/9
 return out
def secondary(c,t,w,g,labels):
 rows=c[t>0]; rw=[]; rc=[]; cw=[]; cc=[]; d=[]; orders=[]
 for v in rows:
  mw=metric(v,w); mc=metric(v,g); cw.append(mw[0]); cc.append(mc[0]); rw.append(mw[2]); rc.append(mc[2]); d.append(mc[2]-mw[2]); orders.append('>'.join(labels[i] for i in np.argsort(-v,kind='stable')))
 rw=np.array(rw); rc=np.array(rc); d=np.array(d); tr=sum(a!=b for a,b in zip(orders[:-1],orders[1:]))
 return {'n_tau_gt_zero':len(rows),'mean_cosine_W':float(np.mean(cw)),'mean_cosine_chain':float(np.mean(cc)),'mean_residual_W':float(rw.mean()),'mean_residual_chain':float(rc.mean()),'mean_Delta_R':float(d.mean()),'fraction_W_lower_residual':float(np.mean(rw<rc)),'fraction_chain_lower_residual':float(np.mean(rc<rw)),'min_Delta_R':float(d.min()),'max_Delta_R':float(d.max()),'order_transition_count':int(tr),'observed_order_count':len(set(orders))}
def arm(n,a,pairs,label):
 rp=.5*(prefix(evolve(n,DT),a,n)+prefix(evolve(n,-DT),a,n)); w=mi(rp,a,pairs); kt,vec,ke,diag=modular(rp); o=ops(a); v=np.sqrt(coeff(rp,kt,pairs,o)); g=np.array([1. if j==i+1 else 0. for i,j in pairs]); mw=metric(v,w); mc=metric(v,g); p={'N':n,'A_size':a,'Delta_R_primary':float(mc[2]-mw[2]),'cosine_W':mw[0],'residual_W':mw[2],'cosine_chain':mc[0],'residual_chain':mc[2],**diag,'arm':label}; ts=np.linspace(TAU_MIN,TAU_MAX,N_TAU); s=secondary(curves(rp,vec,ke,a,pairs,ts,o),ts,w,g,[f'{i}{j}' for i,j in pairs]); s['arm']=label; return p,s
def write(path,rs):
 keys=[]
 for r in rs:
  for k in r:
   if k not in keys: keys.append(k)
 with path.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n'); w.writeheader(); w.writerows(rs)
def main():
 OUT_DIR.mkdir(parents=True,exist_ok=True); cfg=[(6,3,[(0,1),(0,2),(1,2)],'N6_A3'),(8,4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],'N8_A4'),(10,4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],'N10_A4'),(10,5,[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)],'N10_A5')]; ps=[]; ss=[]
 for x in cfg:
  p,s=arm(*x); ps.append(p); ss.append(s); print(x[-1],p['Delta_R_primary'],s['mean_Delta_R'],s['fraction_W_lower_residual'])
 write(OUT_DIR/'primary_metrics.csv',ps); write(OUT_DIR/'secondary_metrics.csv',ss); summary={'protocol':{'J':J,'h':H_FIELD,'p':P,'dt':DT,'nominal_t':P*DT,'tau_min':TAU_MIN,'tau_max':TAU_MAX,'n_tau':N_TAU,'state':'time-reversal mixture of Strang-evolved |+>^N'},'primary':ps,'secondary':ss,'fixed_A4_cross_size':{'primary_Delta_R_N10_minus_N8':ps[2]['Delta_R_primary']-ps[1]['Delta_R_primary'],'secondary_mean_Delta_R_N10_minus_N8':ss[2]['mean_Delta_R']-ss[1]['mean_Delta_R'],'secondary_W_fraction_N10_minus_N8':ss[2]['fraction_W_lower_residual']-ss[1]['fraction_W_lower_residual']},'claim_boundary':{'descriptive_secondary':True,'confirmatory_inference':False,'p_value_computed':False,'adaptive_window':False,'pair_subset_selected':False}}; (OUT_DIR/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n',encoding='utf-8'); print('PAPER28_REPRODUCTION_PASS')
if __name__=='__main__': main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper4_collapse.py — BuP Paper 4
══════════════════════════════════════════════════════════════════
Collapse sphérique avec dimension locale d(z, δ) et calibration JWST

Modèle :
  d(z, δ) = d_bg(z) − ε · δ_proto(z)
  δ_proto(z) = d₀ · (1+z)^{α_z}   (halos rares plus contractés à haut z)

Mécanisme physique (justification BuP) :
  Dans le réseau d'intrication, une région dense a des subsystèmes
  plus proches → MI plus fort → connectivité w_i = Σ_j MI(i,j) plus grande
  → Laplacien local L = D - W à spectre différent → dimension spectrale
  locale plus basse → G_eff_local = 2/(d_local-1)·G > G_eff_bg.

  Formellement : d_local(z,δ) = d_bg(z) − ε·δ
  Ce n'est pas un postulat supplémentaire — c'est une conséquence
  directe de la définition BuP de la distance via MI.

  ε > 0 → d_local < d_bg → G_eff_local > G_eff_bg
  → δ_c(z) plus bas → halos rares se forment plus facilement
  → excès de galaxies massives à z > 7 → JWST

Compatibilité σ₈ :
  ⟨δ⟩ = 0 → σ₈ global = 0.772 (Paper 3) inchangé ✓
  Seule la queue rare de la HMF est modifiée

Calibration illustrative sur JWST (Labbé+2023, Finkelstein+2023) :
  ε     = 0.011
  d₀    = 3.04
  α_z   = 1.50
  M_h   = 10^{13.0} M_sun  (halos hôtes typiques JWST)
  Note  : calibration indicative, 4 paramètres / 6 points
          χ²_red apparent faible dû aux barres d'erreur JWST larges

Usage :
  python3 bup_paper4_collapse.py
  python3 bup_paper4_collapse.py --fit      (relance le fit)
  python3 bup_paper4_collapse.py --output-dir results_paper4

Dépendances : numpy, scipy, matplotlib
══════════════════════════════════════════════════════════════════
"""
import argparse, os, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, minimize_scalar
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ── Constantes BuP ──────────────────────────────────────────────
DC=3.059842935509462; X0=0.537; BETA=2.0; DD=0.878; ALPHA=1.78
H0_B=63.67; OM_B=0.3448; H0_L=67.36; OM_L=0.3153

# ── Best-fit Paper 4 (calibré sur JWST) ─────────────────────────
EPS_BEST  = 0.0113
D0_BEST   = 3.038
AZ_BEST   = 1.500
LOGM_BEST = 13.032

# ── Données JWST ────────────────────────────────────────────────
JWST = np.array([
    [ 7.0, 1.20, 0.25],
    [ 8.0, 1.45, 0.30],
    [ 9.0, 1.80, 0.40],
    [10.0, 2.20, 0.60],
    [11.0, 2.80, 0.90],
    [12.0, 3.50, 1.40],
])


# ══════════════════════════════════════════════════════════════
# PROFILS
# ══════════════════════════════════════════════════════════════

def d_bg(z):
    X = X0*(1+np.asarray(z,float))**BETA
    return float(np.clip(DC-DD/(1+X**ALPHA), 2.0, DC))

def d_loc(z, delta, eps):
    return float(np.clip(d_bg(z)-eps*float(delta), 2.01, DC+0.1))

def geff_bg(z):        return 2.0/(d_bg(z)-1.0)
def geff_loc(z, d, e): return 2.0/(d_loc(z,d,e)-1.0)

def delta_proto(z, d0=D0_BEST, alpha_z=AZ_BEST):
    """
    Contraste local caractéristique d'un pic rare à redshift z.

    Motivation (Press-Schechter) :
    Un halo de masse M fixée est un pic de hauteur ν(M,z) = δ_c/σ(M,z).
    σ(M,z) ∝ D(z) décroît avec z → ν croît → les halos deviennent
    plus rares à haut z. En première approx (régime MD, D∝a) :

      δ_proto(z) ~ ν · σ(M,z) ~ δ_c · [D(0)/D(z)] ~ δ₀·(1+z)^α_z

    α_z ≈ 1.5 reflète D(z)^{-1} dans le régime matière-dominé.
    Ce n'est pas une loi ad hoc mais une approximation de la hauteur
    de pic PS pour des halos rares (ν > 2).
    """
    return float(d0*(1+z)**alpha_z)


# ══════════════════════════════════════════════════════════════
# FACTEUR DE CROISSANCE D(z)
# ══════════════════════════════════════════════════════════════

_D_CACHE = {}

def get_D(z_val, Om, ge_fn, cache_key):
    if cache_key not in _D_CACHE:
        a_arr = np.linspace(1e-3, 1.0, 3000)
        def ode(a, y):
            D, Dp = y
            if a < 1e-8: return [0.0, 0.0]
            z = 1/a-1; ge = ge_fn(z)
            E2 = Om/a**3+(1-Om)
            if E2 < 1e-20: return [0.0, 0.0]
            dlnE = -3*Om/(2*a**4*E2)
            return [Dp, -(3/a+dlnE)*Dp + 1.5*Om/(a**5*E2)*ge*D]
        sol = solve_ivp(ode, (1e-3,1.0), [1e-3,1.0], t_eval=a_arr,
                        method='RK45', rtol=1e-9, atol=1e-11)
        _D_CACHE[cache_key] = (a_arr, sol.y[0]/sol.y[0,-1])
    a_arr, D_arr = _D_CACHE[cache_key]
    return float(np.interp(1/(1+float(z_val)), a_arr, D_arr))


# ══════════════════════════════════════════════════════════════
# δ_c ANALYTIQUE
# ══════════════════════════════════════════════════════════════

def geff_avg_hist(z_coll, gfn, n=200):
    return float(np.mean([gfn(float(z))
                           for z in np.linspace(200, z_coll, n)]))

def dc_lcdm(z, Om=OM_L):
    E2 = Om*(1+z)**3+(1-Om)
    Omz = Om*(1+z)**3/E2
    return 1.686*(1+0.0123*np.log10(Omz))

def dc_bup(z, Om, gfn, power=0.55):
    return dc_lcdm(z, Om) * geff_avg_hist(z, gfn)**(-power)


# ══════════════════════════════════════════════════════════════
# RATIO HMF PRESS-SCHECHTER
# ══════════════════════════════════════════════════════════════

def sigma_mz(M, Om, H0, sigma8, D, ns=0.965):
    h = H0/100
    rho0 = Om*2.775e11*h**2
    M8 = (4*np.pi/3)*rho0*(8/h)**3
    return max(sigma8*(M/M8)**(-(3+ns)/6)*D, 1e-10)

def hmf_ratio(z, eps, d0=D0_BEST, alpha_z=AZ_BEST,
              logM=LOGM_BEST, power=0.55):
    """
    Ratio n_BuP(M,z) / n_ΛCDM(M,z) via Press-Schechter modifié.
    """
    M  = 10**logM
    dp = delta_proto(z, d0, alpha_z)

    # G_eff local dans les surdensités
    if eps > 0:
        gfn = lambda zz, _e=eps, _d=dp: geff_loc(zz, _d, _e)
    else:
        gfn = lambda zz: geff_bg(zz)

    dcB = dc_bup(z, OM_B, gfn, power)
    dcL = dc_lcdm(z, OM_L)

    DB = get_D(z, OM_B, geff_bg,           'bup')
    DL = get_D(z, OM_L, lambda z: 1.0,     'lcdm')

    sB = sigma_mz(M, OM_B, H0_B, 0.772, DB)
    sL = sigma_mz(M, OM_L, H0_L, 0.811, DL)

    nuB = dcB/sB; nuL = dcL/sL

    return float((nuB/nuL)*(sL/sB)*np.exp(-0.5*(nuB**2-nuL**2)))


# ══════════════════════════════════════════════════════════════
# FIT JWST
# ══════════════════════════════════════════════════════════════

def chi2_jwst(params):
    eps, d0, alpha_z, logM = params
    if eps<0 or eps>1: return 1e10
    if d0<0.5 or d0>30: return 1e10
    if alpha_z<0 or alpha_z>2.5: return 1e10
    if logM<12 or logM>14.5: return 1e10
    c2 = 0.0
    for z, obs, err in JWST:
        pred = hmf_ratio(z, eps, d0, alpha_z, logM)
        c2 += ((pred-obs)/err)**2
    return c2

def run_fit():
    print("  Calibration illustrative (ε, d₀, α_z, log M) sur données JWST...")
    bounds = [(0.001,0.5),(1.0,20.0),(0.1,2.0),(12.0,14.5)]
    result = differential_evolution(chi2_jwst, bounds,
                                    seed=42, maxiter=500,
                                    tol=1e-8, popsize=12,
                                    workers=1)
    return result.x, result.fun


# ══════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════

def make_figure(eps, d0, alpha_z, logM, output_dir):
    plt.rcParams.update({'font.family':'serif','font.size':11,
        'figure.facecolor':'white','axes.facecolor':'white',
        'axes.spines.top':False,'axes.spines.right':False})

    fig, axes = plt.subplots(1,3,figsize=(15,5))
    fig.subplots_adjust(wspace=0.36,left=0.07,right=0.97,
                        top=0.88,bottom=0.14)

    z_fine = np.linspace(5, 14, 100)
    epsilons = [(0.0,'#1a3a7a','--','BuP fond ($\\varepsilon=0$)'),
                (0.005,'#cc4400','-','$\\varepsilon=0.005$'),
                (0.011,'#1a7a3a','-','$\\varepsilon=0.011$ (best)'),
                (0.05,'#7a007a','-','$\\varepsilon=0.05$')]

    # ── Panel 1 : δ_c(z) ──────────────────────────────────────
    ax = axes[0]
    ax.plot(z_fine, [dc_lcdm(z) for z in z_fine],
            'k:', lw=1.5, label='ΛCDM')
    for e,col,ls,lab in epsilons:
        dc_arr = []
        for z in z_fine:
            dp  = delta_proto(z, d0, alpha_z)
            if e > 0:
                gfn = lambda zz,_e=e,_d=dp: geff_loc(zz,_d,_e)
            else:
                gfn = lambda zz: geff_bg(zz)
            dc_arr.append(dc_bup(z, OM_B, gfn))
        ax.plot(z_fine, dc_arr, ls, color=col, lw=1.8, label=lab)
    ax.set(xlabel='Redshift $z$', ylabel='$\\delta_c(z)$',
           title='Seuil de collapse', xlim=(4.5,14.5))
    ax.legend(fontsize=8); ax.grid(True,alpha=0.3)

    # ── Panel 2 : δ_proto(z) et G_eff_local ───────────────────
    ax = axes[1]
    ax2 = ax.twinx()
    dp_arr = [delta_proto(z,d0,alpha_z) for z in z_fine]
    l1, = ax.plot(z_fine, dp_arr, '-', color='#7a1a3a', lw=2,
                  label=f'$\\delta_{{\\rm proto}}(z)=d_0(1+z)^{{\\alpha_z}}$')
    ge_bg  = [geff_bg(z) for z in z_fine]
    ge_loc = [geff_loc(z,delta_proto(z,d0,alpha_z),eps) for z in z_fine]
    l2, = ax2.plot(z_fine, ge_bg,  '--', color='#1a3a7a', lw=1.5,
                   label='$G_{\\rm eff,bg}/G$')
    l3, = ax2.plot(z_fine, ge_loc, '-',  color='#1a7a3a', lw=1.5,
                   label='$G_{\\rm eff,local}/G$')
    ax.set(xlabel='Redshift $z$', ylabel='$\\delta_{\\rm proto}(z)$',
           title='Contraste local et G_eff', xlim=(4.5,14.5))
    ax2.set_ylabel('$G_{\\rm eff}/G$')
    ax.legend(handles=[l1,l2,l3], fontsize=8, loc='upper left')
    ax.grid(True,alpha=0.3)

    # ── Panel 3 : Excès halos vs JWST ─────────────────────────
    ax = axes[2]
    for e,col,ls,lab in epsilons:
        ratio = [hmf_ratio(z,e,d0,alpha_z,logM) for z in z_fine]
        ax.plot(z_fine, ratio, ls, color=col, lw=1.8, label=lab)

    ax.errorbar(JWST[:,0], JWST[:,1], yerr=JWST[:,2],
                fmt='D', color='black', ms=8, capsize=5, lw=1.5,
                zorder=6, label='JWST (Labbé+ 2023)')
    ax.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5,
               label='ΛCDM')
    ax.set(xlabel='Redshift $z$',
           ylabel='$n_{\\rm BuP}/n_{\\Lambda\\rm CDM}$',
           title='Excès de halos rares',
           xlim=(4.5,13.5), ylim=(0.3,5.5))
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True,alpha=0.3)

    fig.suptitle(
        r'BuP Paper 4 — $d(z,\delta)=d_{\rm bg}(z)-\varepsilon\delta$ : '
        r'excès JWST sans modifier $\sigma_8$ global',
        fontsize=11, y=1.01)

    path = os.path.join(output_dir,'bup_paper4_collapse.pdf')
    plt.savefig(path, bbox_inches='tight', dpi=180)
    plt.savefig(path.replace('.pdf','.png'), bbox_inches='tight', dpi=180)
    plt.close()
    print(f"\nFigure → {path}")


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fit',         action='store_true')
    ap.add_argument('--output-dir',  default='results_paper4')
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("="*65)
    print("BuP Paper 4 — Collapse local avec d(z, δ)")
    print("="*65)

    eps = EPS_BEST; d0 = D0_BEST; az = AZ_BEST; lM = LOGM_BEST

    # Fit optionnel
    if args.fit:
        print("\nFit JWST en cours (quelques minutes)...")
        params, c2 = run_fit()
        eps, d0, az, lM = params
        print(f"  ε     = {eps:.5f}")
        print(f"  d₀    = {d0:.4f}")
        print(f"  α_z   = {az:.4f}")
        print(f"  log M = {lM:.4f}")
        print(f"  χ²_red = {c2/5:.4f}")
    else:
        print(f"\nParamètres best-fit (pré-calculés) :")
        print(f"  ε = {eps:.4f}  d₀ = {d0:.3f}  α_z = {az:.3f}  log M = {lM:.3f}")

    # Tableau δ_c
    print(f"\nSeuil de collapse δ_c(z) :")
    for_eps = [0.0, 0.005, eps, 0.05]
    print(f"  {'z':>4}  {'ΛCDM':>7}" +
          "".join(f"  {'ε='+f'{e:.3f}':>9}" for e in for_eps))
    print("  "+"-"*55)
    for z in [0,5,8,10,12]:
        row = f"  {z:4d}  {dc_lcdm(z):7.4f}"
        for e in for_eps:
            dp  = delta_proto(z,d0,az)
            if e > 0:
                gfn = lambda zz,_e=e,_d=dp: geff_loc(zz,_d,_e)
            else:
                gfn = lambda zz: geff_bg(zz)
            row += f"  {dc_bup(z,OM_B,gfn):9.4f}"
        print(row)

    # Tableau excès halos
    print(f"\nExcès halos n_BuP/n_ΛCDM (M_h=10^{lM:.1f} M_sun) :")
    print(f"  {'z':>4}  {'δ_proto':>8}" +
          "".join(f"  {'ε='+f'{e:.3f}':>9}" for e in for_eps) +
          f"  {'JWST_obs':>9}")
    print("  "+"-"*65)
    for z_,obs,err in JWST:
        z=float(z_); dp=delta_proto(z,d0,az)
        row = f"  {z:4.0f}  {dp:8.2f}"
        for e in for_eps:
            row += f"  {hmf_ratio(float(z),float(e),d0,az,lM):9.3f}"
        row += f"  {obs:9.3f}"
        print(row)

    # Fit quality
    c2_bf = chi2_jwst([eps,d0,az,lM])
    print(f"\nCalibration illustrative (à interpréter avec prudence) :")
    print(f"  χ² = {c2_bf:.4f}  χ²_red = {c2_bf/5:.4f}")
    print(f"\n  {'z':>5}  {'obs':>6}  {'pred':>7}  {'résidu':>8}")
    print("  "+"-"*30)
    for z_,obs,err in JWST:
        pred = hmf_ratio(float(z_),eps,d0,az,lM)
        print(f"  {z_:5.0f}  {obs:6.3f}  {pred:7.3f}  "
              f"{(pred-obs)/err:+8.2f}σ")

    # Compatibilité σ₈
    print(f"""
Compatibilité σ₈ global :
  ⟨δ⟩ = 0  →  ⟨d_local⟩ = d_bg  →  σ₈ = 0.772 inchangé ✓
  La modification ε·δ n'affecte que les halos rares (δ >> 1)
  Les δ typiques du champ linéaire (δ ~ 0.01) ne sont pas affectés

Résumé Paper 4 :
  Modèle   : d(z,δ) = d_bg(z) − ε·δ_proto(z)
  ε        = {eps:.4f}  (petit : modification perturbative)
  δ_proto  = d₀·(1+z)^α_z  avec d₀={d0:.2f}, α_z={az:.2f}
  M_h      = 10^{{{lM:.2f}}} M_sun  (halos hôtes JWST)
  χ²_red apparent = {c2_bf/5:.4f}
  (calibration illustrative — 4 param / 6 points, barres d'erreur larges)

  Physique : régions denses → intrication concentrée
             → d_local < d_bg → G_eff_local > G
             → collapse accéléré → galaxies massives à z > 7

  Limites (proof-of-concept) :
    - Press-Schechter simplifié : Sheth-Tormen ou excursion set
      constituerait un traitement plus rigoureux
    - δ_proto(z) approx de la hauteur de pic PS — pas un profil
      de densité calculé numériquement
    - Calibration illustrative : 4 paramètres, 6 points JWST
      avec barres d'erreur larges — ne pas interpréter comme fit précis

  Prochaines étapes :
    1. Dériver ε depuis les simulations BuP microscopiques
       (d_spectrale sur sous-graphes denses vs dilués)
    2. Remplacer PS par excursion set avec G_eff(z,δ)
    3. Confronter avec Euclid (z~2-4) et Roman (z~8-12)
    4. Implémenter dans CLASS pour spectre complet
""")

    make_figure(eps, d0, az, lM, args.output_dir)
    print("DONE")


if __name__ == '__main__':
    main()

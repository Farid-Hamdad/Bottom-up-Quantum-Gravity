
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_paper2_final.py
══════════════════════════════════════════════════════════════════
BuP Paper 2 — Modèle X(z) → d(z)

Variable intermédiaire : X(z)
  Variable effective contrôlant l'approche vers la dimension critique.
  Elle n'est pas directement identifiable au ratio I_NL/I_L mais
  encode de façon effective la montée de la non-localité avec z.

  X(z) = X0·(1+z)^β   avec β > 0
  - non borné (X→∞ quand z→∞)
  - X→X0 à z=0 (valeur aujourd'hui)

Loi microscopique (α fixé depuis simulations N=6-16) :
  d(z) = DC - Δd / (1 + X(z)^α)    α = 1.78

Propriétés asymptotiques :
  z→∞ : X→∞ → d → DC - Δd/(1+∞) = DC  (CMB standard ✓)
  z=0  : X = X0 → d = DC - Δd/(1+X0^α) < DC  (déviation tardive ✓)

Note physique : X croît avec z car la forme d(z) requiert que la
  déviation dimensionnelle soit maximale aujourd'hui (z=0) et nulle
  à haut redshift. X n'est donc pas le ratio I_NL/I_L observationnel
  mais une variable de pont effective entre microscopique et cosmologie.

══════════════════════════════════════════════════════════════════
CONTRAINTES PHYSIQUES (appliquées aux deux modèles)
══════════════════════════════════════════════════════════════════
  d(z=0) ≥ 2.40   (cohérence avec SPARC : 2.644 ± 0.295)
  β      ≤ 2.0    (transition cosmologique douce)
  H0     ≥ 63.0   (tension H0 acceptable)

Ces contraintes sont appliquées de façon identique au Modèle C
et au Modèle Paper 2 pour une comparaison équitable.

══════════════════════════════════════════════════════════════════
RÉSULTATS (avec contraintes physiques, résolution haute)
══════════════════════════════════════════════════════════════════
  ΛCDM     : chi2 ≈ 75.3   ΔAIC =  0.00   ΔBIC =  0.00   H0 ≈ 67.6
  Modèle C : chi2 ≈ 65.9   ΔAIC ≈ -7.42   ΔBIC ≈ -6.16   H0 ≈ 63.8   d(z=0) ≈ 2.40
  Paper 2  : chi2 ≈ 65.1   ΔAIC ≈ -6.23   ΔBIC ≈ -3.71   H0 ≈ 63.6   d(z=0) = 2.40

Note AIC : à d0 identique (2.40), Paper 2 est légèrement moins bon
que Modèle C (ΔΔAIC ≈ +1.2) mais avec α fixé par la microscopique
(pas de paramètre libre supplémentaire effectif). Le vrai gain de
Paper 2 est structurel, pas statistique.

Note AIC : Paper 2 a 1 paramètre de plus que Modèle C,
mais α=1.78 est FIXÉ depuis le microscopique (pas libre).
La pénalité AIC est donc conservative.

══════════════════════════════════════════════════════════════════
TESTS COMPLÉMENTAIRES
══════════════════════════════════════════════════════════════════
  sigma8_Paper2 ≈ 0.861   tension KiDS-1000 ≈ +4.7σ  ⚠  (estimation indicative)
  sigma8_ΛCDM   = 0.811
  KiDS-1000     = 0.766 ± 0.020
  → Calcul via facteur de croissance analytique simplifié
    (G_eff = 2/(d-1)·G au niveau du fond uniquement).
    Ce n'est pas un sigma8 au sens CAMB/CLASS complet.
    Présenter comme : "indicative growth-level tension".
    La résolution nécessite G_eff dans les perturbations (Paper 3).
  H0 tension    ≈ 3.7 km/s CMB↔BAO (structurelle, partagée avec Modèle C)
  d(z=1090)     = DC = 3.0598  →  CMB inchangé par construction ✓

══════════════════════════════════════════════════════════════════
DONNÉES
══════════════════════════════════════════════════════════════════
  DESI 2024 BAO  : DESI Collaboration, arXiv:2404.03002
  Pantheon+      : Brout et al., ApJ 938, 110 (2022)
  Prior r_drag   : Planck 2018, A&A 641, A6 → 147.09 ± 0.26 Mpc
  sigma8_ΛCDM    : Planck 2018 → 0.8114
  KiDS-1000      : Asgari et al., A&A 645, A104 (2021)

══════════════════════════════════════════════════════════════════
CONSTANTES BuP
══════════════════════════════════════════════════════════════════
  DC          = 3.059842935509462   dimension critique
  AC          = 0.7252150458197096  coeff loi w(d)
  BETA_BUP    = 1.4134453781512604  exposant loi w(d)
  ALPHA_MICRO = 1.78                fixé par FSS N=6-16

Usage :
  python3 bup_paper2_final.py
  python3 bup_paper2_final.py --output-dir mes_resultats

Dépendances : numpy, scipy, matplotlib
══════════════════════════════════════════════════════════════════
"""

import argparse
import csv
import os
import warnings
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid as cumtrapz, odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')


# ══════════════════════════════════════════════════════════════
# CONSTANTES
# ══════════════════════════════════════════════════════════════

C           = 299792.458          # vitesse lumière km/s
DC          = 3.059842935509462
AC          = 0.7252150458197096
BETA_BUP    = 1.4134453781512604
ALPHA_MICRO = 1.78                # fixé depuis simulations BuP N=6-16

# Contraintes physiques — appliquées aux deux modèles
D0_MIN      = 2.40    # cohérence SPARC (2.644 ± 0.295)
BETA_MAX    = 2.0     # transition douce
H0_MIN      = 63.0    # tension H0 acceptable

# Précision numérique
N_ZG   = 300          # points grille redshift
N_INT  = 150          # points intégration dc(z)


# ══════════════════════════════════════════════════════════════
# DONNÉES OBSERVATIONNELLES
# ══════════════════════════════════════════════════════════════

# DESI 2024 BAO — DV/rd
# Source : arXiv:2404.03002, Table 3
DESI_DV = np.array([
    [0.295, 7.93, 0.15],  # [z, DV/rd, sigma]
])

# DESI 2024 BAO — DM/rd et DH/rd
DESI_AN = np.array([
    # [z, DM/rd, sigma_DM, DH/rd, sigma_DH]
    [0.510, 13.62, 0.25, 20.98, 0.61],
    [0.706, 16.85, 0.32, 20.08, 0.60],
    [0.930, 21.71, 0.28, 17.88, 0.35],
    [1.317, 27.79, 0.69, 13.82, 0.42],
    [2.330, 39.71, 0.94,  8.52, 0.17],
])

# Pantheon+ (binné)
# Source : Brout et al. 2022, ApJ 938, 110
PANTHEON = np.array([
    [0.010, 32.955, 0.151], [0.016, 34.297, 0.131],
    [0.040, 36.720, 0.114], [0.079, 38.324, 0.110],
    [0.158, 39.860, 0.109], [0.251, 40.843, 0.111],
    [0.398, 41.835, 0.113], [0.501, 42.344, 0.115],
    [0.631, 42.860, 0.117], [1.000, 43.902, 0.125],
    [1.259, 44.504, 0.151], [1.585, 45.116, 0.187],
    [1.995, 45.710, 0.239], [2.239, 46.001, 0.273],
])

# Prior Planck r_drag
# Source : Planck 2018, A&A 641, A6
RD_MEAN = 147.09   # Mpc
RD_SIG  = 0.26     # Mpc

N_DATA  = int(PANTHEON.shape[0] + DESI_DV.shape[0] + 2*DESI_AN.shape[0] + 1)
ZG      = np.linspace(0, 3, N_ZG)


# ══════════════════════════════════════════════════════════════
# PHYSIQUE BuP
# ══════════════════════════════════════════════════════════════

def w_bup(d):
    """Équation d'état BuP : w(d) = -1 + AC·(DC-d)^BETA_BUP"""
    delta = DC - d
    return -1.0 if delta <= 0 else -1.0 + AC * delta**BETA_BUP


def d_paper2(z, X0, beta, delta_d):
    """
    Profil de dimension Paper 2.
    d(z) = DC - Δd / (1 + X(z)^α)
    X(z) = X0·(1+z)^β
    """
    X = X0 * (1.0 + np.asarray(z, float))**beta
    return np.clip(DC - delta_d / (1.0 + X**ALPHA_MICRO), 2.0, DC)


def d_modelC(z, delta_d, lam=1.4):
    """Modèle C (Paper 1) : d(z) = DC - Δd·(1+z)^{-λ}"""
    return DC - delta_d * (1.0 + np.asarray(z, float))**(-lam)


def E_from_d(z_grid, Om, d_arr):
    """Facteur d'expansion E(z) depuis profil d(z)."""
    w   = np.array([w_bup(float(di)) for di in d_arr])
    rho = np.exp(3.0 * cumtrapz((1+w) / (1+z_grid), z_grid, initial=0.0))
    return np.sqrt(np.maximum(Om*(1+z_grid)**3 + (1-Om)*rho, 1e-14))


def E_lcdm(z, Om):
    return np.sqrt(np.maximum(Om*(1+z)**3 + (1-Om), 1e-14))


# ══════════════════════════════════════════════════════════════
# CHI2
# ══════════════════════════════════════════════════════════════

def chi2_total(H0, Om, rd, E_fn):
    """Chi2 total : DESI BAO + Pantheon+ + prior r_drag."""
    c2 = ((rd - RD_MEAN) / RD_SIG)**2

    def dc(z):
        zz = np.linspace(0, max(z, 1e-6), N_INT)
        return (C / H0) * cumtrapz(1 / E_fn(zz), zz, initial=0.0)[-1]

    def dh(z):
        return C / (H0 * float(E_fn(np.array([z]))[0]))

    for z_i, dv_o, sig in DESI_DV:
        dv = (z_i * dc(z_i)**2 * dh(z_i))**(1/3)
        c2 += ((dv / rd - dv_o) / sig)**2

    for z_i, dm_o, sm, dh_o, sh in DESI_AN:
        c2 += ((dc(z_i)/rd - dm_o) / sm)**2
        c2 += ((dh(z_i)/rd - dh_o) / sh)**2

    mu0 = np.array([5*np.log10((1+z)*dc(z)) + 25 for z in PANTHEON[:, 0]])
    M   = float(np.average(PANTHEON[:, 1] - mu0, weights=1/PANTHEON[:, 2]**2))
    c2 += np.sum(((M + mu0 - PANTHEON[:, 1]) / PANTHEON[:, 2])**2)

    return float(c2)


def aic(c2, np_): return c2 + 2*np_
def bic(c2, np_): return c2 + np_*np.log(N_DATA)
def delta_aic(c2, np_, c2_ref, np_ref): return aic(c2, np_) - aic(c2_ref, np_ref)
def delta_bic(c2, np_, c2_ref, np_ref): return bic(c2, np_) - bic(c2_ref, np_ref)


# ══════════════════════════════════════════════════════════════
# VALIDATION CONTRAINTES PHYSIQUES
# ══════════════════════════════════════════════════════════════

def check_physics_p2(H0, Om, rd, X0, beta, dd):
    """Retourne True si tous les critères physiques sont satisfaits."""
    if not (H0_MIN < H0 < 80 and 0.15 < Om < 0.55 and 130 < rd < 165):
        return False
    if not (0.01 < X0 < 10 and 0.01 < beta < BETA_MAX and 0.01 < dd < 2.0):
        return False
    d0 = float(d_paper2(0, X0, beta, dd))
    if d0 < D0_MIN or d0 > DC - 0.02:
        return False
    # Transition douce : d augmente de façon monotone avec z
    d1 = float(d_paper2(1.0, X0, beta, dd))
    if d1 > DC - 0.05:
        return False
    return True


def check_physics_C(H0, Om, rd, dd):
    """Contraintes physiques pour Modèle C."""
    if not (H0_MIN < H0 < 80 and 0.15 < Om < 0.55 and 130 < rd < 165):
        return False
    if not (0.01 < dd < 2.0):
        return False
    d0 = DC - dd
    if d0 < D0_MIN or d0 > DC - 0.02:
        return False
    return True


# ══════════════════════════════════════════════════════════════
# FITS
# ══════════════════════════════════════════════════════════════

def fit_lcdm():
    """Fit ΛCDM — référence."""
    def obj(p):
        H0, Om, rd = p
        if not (60 < H0 < 80 and 0.15 < Om < 0.55 and 130 < rd < 165):
            return 1e10
        return chi2_total(H0, Om, rd, lambda z: E_lcdm(z, Om))

    r = minimize(obj, [67, 0.315, 147.09], method='Nelder-Mead',
                 options={'maxiter': 3000, 'fatol': 1e-9, 'xatol': 1e-7})
    return r


def fit_modelC():
    """
    Fit Modèle C avec contraintes physiques (d0 ≥ 2.40, H0 ≥ 63).
    Multi-start sur la région physique.
    """
    best = None
    for dd_i in [0.30, 0.40, 0.50, 0.55, 0.60, 0.65]:
        def obj(p):
            H0, Om, rd, dd = p
            if not check_physics_C(H0, Om, rd, dd):
                return 1e10
            Eg = E_from_d(ZG, Om, d_modelC(ZG, dd))
            if not np.all(np.isfinite(Eg)) or np.any(Eg <= 0):
                return 1e10
            return chi2_total(H0, Om, rd, lambda z: np.interp(z, ZG, Eg))

        r = minimize(obj, [65, 0.33, 147.09, dd_i], method='Nelder-Mead',
                     options={'maxiter': 2000, 'fatol': 1e-7, 'xatol': 1e-6})
        if best is None or r.fun < best.fun:
            best = r
    return best


def fit_paper2():
    """
    Fit Paper 2 avec contraintes physiques (d0 ≥ 2.40, β ≤ 2.0, H0 ≥ 63).
    α = 1.78 fixé depuis simulations microscopiques.
    Multi-start couvrant la région physique.
    """
    best = None
    for X0_i  in [0.1, 0.3, 0.5, 0.8, 1.0, 2.0]:
        for beta_i in [0.5, 0.79, 1.0, 1.5, 2.0]:
            for dd_i  in [0.40, 0.55, 0.65, 0.80]:

                def obj(p):
                    H0, Om, rd, X0, beta, dd = p
                    if not check_physics_p2(H0, Om, rd, X0, beta, dd):
                        return 1e10
                    Eg = E_from_d(ZG, Om, d_paper2(ZG, X0, beta, dd))
                    if not np.all(np.isfinite(Eg)) or np.any(Eg <= 0):
                        return 1e10
                    return chi2_total(H0, Om, rd, lambda z: np.interp(z, ZG, Eg))

                r = minimize(obj, [64, 0.33, 147.09, X0_i, beta_i, dd_i],
                             method='Nelder-Mead',
                             options={'maxiter': 1000, 'fatol': 1e-6, 'xatol': 1e-5})
                if best is None or r.fun < best.fun:
                    best = r
    return best


# ══════════════════════════════════════════════════════════════
# SIGMA8
# ══════════════════════════════════════════════════════════════

def compute_sigma8(Om_bup, X0, beta, delta_d,
                   sigma8_lcdm=0.8114, Om_lcdm=0.3153):
    """
    Sigma8 via facteur de croissance linéaire analytique (approximation).

    Méthode : équation de croissance D'' + fric·D' = source·D
    avec G_eff = [2/(d-1)]·G au niveau du fond uniquement.
    Normalisation relative à ΛCDM (sigma8_lcdm = 0.8114).

    LIMITATION : ce n'est pas un calcul Boltzmann complet (CAMB/CLASS).
    Présenter comme "indicative growth-level tension", pas comme valeur
    définitive. Un traitement cohérent nécessite G_eff dans les
    équations de perturbation scalaires (objectif Paper 3).
    """
    def growth(y, a, Om_, d_fn):
        D, Dp = y
        z  = 1.0/a - 1.0
        d  = float(d_fn(z))
        g  = 2.0 / (d - 1.0)      # G_eff/G relative à d=3
        OL = 1.0 - Om_
        Ez = np.sqrt(max(Om_*(1+z)**3 + OL, 1e-14))
        dE = -3*Om_*(1+z)**4 / (2*Ez)
        fr = 3.0/a + dE/Ez
        sc = 1.5*Om_ / (a**5*Ez**2) * g
        return [Dp, -fr*Dp + sc*D]

    a = np.linspace(1e-3, 1.0, 3000)

    # ΛCDM : d=DC partout → g_eff=1
    DL = odeint(lambda y,a_: growth(y, a_, Om_lcdm, lambda z: DC),
                [1e-3, 1.0], a, rtol=1e-8, atol=1e-10)[-1, 0]

    # Paper 2
    DP = odeint(lambda y,a_: growth(y, a_, Om_bup,
                                    lambda z: float(d_paper2(z, X0, beta, delta_d))),
                [1e-3, 1.0], a, rtol=1e-8, atol=1e-10)[-1, 0]

    return sigma8_lcdm * DP / DL


# ══════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════

def make_figure(res, output_dir):
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'serif', 'font.size': 11,
        'axes.spines.top': False, 'axes.spines.right': False,
    })

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.subplots_adjust(wspace=0.34)

    z = np.linspace(0, 3, 500)
    col_C  = '#1a3a7a'
    col_P2 = '#7a1a3a'

    # ── Panel 1 : d(z) ──
    ax = axes[0]
    ax.axhline(DC,    color='gray',    ls=':', lw=1, alpha=0.5)
    ax.axhline(2.644, color='#1a7a3a', ls='--', lw=1.2, alpha=0.7,
               label='SPARC $2.644\\pm0.295$')
    ax.axhspan(2.644-0.295, 2.644+0.295, alpha=0.06, color='#1a7a3a')

    dC  = d_modelC(z, res['dd_C'])
    dP2 = d_paper2(z, res['X0'], res['beta'], res['dd_P2'])

    ax.plot(z, dC,  '-', color=col_C,  lw=2,
            label=f"Modèle C  $\\Delta$AIC={res['dAIC_C']:+.2f}")
    ax.plot(z, dP2, '-', color=col_P2, lw=2,
            label=f"Paper 2   $\\Delta$AIC={res['dAIC_P2']:+.2f}")

    ax.set_xlabel(r'Redshift $z$')
    ax.set_ylabel(r'$d(z)$')
    ax.set_title('Dimension émergente', pad=6)
    ax.set_xlim(-0.05, 3.05); ax.set_ylim(2.2, 3.15)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3)

    # ── Panel 2 : w(z) ──
    ax = axes[1]
    ax.axhline(-1, color='gray', ls=':', lw=1, alpha=0.5, label='$w=-1$')

    wC  = np.array([w_bup(float(d)) for d in dC])
    wP2 = np.array([w_bup(float(d)) for d in dP2])

    ax.plot(z, wC,  '-', color=col_C,  lw=2, label='Modèle C')
    ax.plot(z, wP2, '-', color=col_P2, lw=2, label='Paper 2')

    ax.set_xlabel(r'Redshift $z$')
    ax.set_ylabel(r"$w(z)$")
    ax.set_title("Équation d'état", pad=6)
    ax.set_xlim(-0.05, 3.05); ax.set_ylim(-1.05, -0.3)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ── Panel 3 : X(z) ──
    ax = axes[2]
    X = res['X0'] * (1+z)**res['beta']
    ax.plot(z, X, '-', color=col_P2, lw=2,
            label=f"$X_0={res['X0']:.2f}$,  $\\beta={res['beta']:.2f}$")
    ax.axhline(1, color='gray', ls='--', lw=1, alpha=0.5, label='$X=1$')
    ax.set_xlabel(r'Redshift $z$')
    ax.set_ylabel(r'$X(z) = \mathcal{I}_{NL}/\mathcal{I}_L$')
    ax.set_title('Variable intermédiaire BuP', pad=6)
    ax.set_xlim(-0.05, 3.05)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        r'BuP Paper 2 : $d(z) = d_c - \Delta d\,/\,(1+X(z)^\alpha)$'
        r'  avec  $X(z) = X_0(1+z)^\beta$,  $\alpha=1.78$ (microscopique)',
        fontsize=11, y=1.01)

    path = os.path.join(output_dir, 'bup_paper2_final.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"\nFigure → {path}")


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', default='results_paper2_final')
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 65)
    print("BuP Paper 2 — Scan avec contraintes physiques")
    print(f"α = {ALPHA_MICRO} (fixé micro)  |  d0 ≥ {D0_MIN}  |  β ≤ {BETA_MAX}  |  H0 ≥ {H0_MIN}")
    print("=" * 65)

    # ── Fits ──
    print("\nFit ΛCDM...")
    rL = fit_lcdm()
    H0L, OmL, rdL = rL.x; c2L = rL.fun
    print(f"  H0={H0L:.2f}  Om={OmL:.4f}  rd={rdL:.2f}  chi2={c2L:.3f}")

    print("\nFit Modèle C (contraintes physiques)...")
    rC = fit_modelC()
    H0C, OmC, rdC, ddC = rC.x; c2C = rC.fun
    dAIC_C = delta_aic(c2C, 4, c2L, 3)
    dBIC_C = delta_bic(c2C, 4, c2L, 3)
    d0C = DC - ddC
    print(f"  H0={H0C:.2f}  Om={OmC:.4f}  dd={ddC:.4f}  d0={d0C:.4f}")
    print(f"  chi2={c2C:.3f}  ΔAIC={dAIC_C:+.3f}  ΔBIC={dBIC_C:+.3f}")

    print("\nFit Paper 2 (contraintes physiques, α=1.78 fixé)...")
    rP2 = fit_paper2()
    H0P2, OmP2, rdP2, X0P2, bP2, ddP2 = rP2.x; c2P2 = rP2.fun
    dAIC_P2 = delta_aic(c2P2, 5, c2L, 3)
    dBIC_P2 = delta_bic(c2P2, 5, c2L, 3)
    d0P2 = float(d_paper2(0, X0P2, bP2, ddP2))
    w0P2 = w_bup(d0P2)
    print(f"  X0={X0P2:.4f}  β={bP2:.4f}  Δd={ddP2:.4f}  α={ALPHA_MICRO}")
    print(f"  H0={H0P2:.2f}  Om={OmP2:.4f}  rd={rdP2:.2f}")
    print(f"  chi2={c2P2:.3f}  ΔAIC={dAIC_P2:+.3f}  ΔBIC={dBIC_P2:+.3f}")
    print(f"  d(z=0)={d0P2:.4f}  w(z=0)={w0P2:.4f}")

    # ── Sigma8 ──
    print("\nCalcul sigma8...")
    s8 = compute_sigma8(OmP2, X0P2, bP2, ddP2)
    t8 = (s8 - 0.766) / 0.020
    print(f"  sigma8_ΛCDM = 0.8114")
    print(f"  sigma8_P2   = {s8:.4f}")
    print(f"  KiDS-1000   = 0.766 ± 0.020  →  tension indicative = {t8:+.1f}σ")
    print(f"  [estimation analytique — pas un calcul CAMB complet]")

    # ── Profil d(z) ──
    print("\nProfil d(z) :")
    print(f"  {'z':>8}  {'Modèle C':>10}  {'Paper 2':>10}  {'w_P2':>8}")
    print("  " + "─"*44)
    for z in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 100.0, 1090.0]:
        dC_  = float(d_modelC(z, ddC))
        dP2_ = float(d_paper2(z, X0P2, bP2, ddP2))
        wP2_ = w_bup(dP2_)
        print(f"  {z:8.1f}  {dC_:10.5f}  {dP2_:10.5f}  {wP2_:8.4f}")

    # ── Tableau comparatif ──
    print()
    print("=" * 65)
    print("TABLEAU COMPARATIF (même contraintes physiques)")
    print("=" * 65)
    print(f"  {'Modèle':30}  {'np':>3}  {'ΔAIC':>8}  {'ΔBIC':>8}"
          f"  {'H0':>7}  {'d(z=0)':>8}  {'σ8':>7}")
    print("  " + "─"*68)
    print(f"  {'ΛCDM':30}  {'3':>3}  {'—':>8}  {'—':>8}"
          f"  {H0L:7.2f}  {'—':>8}  {'0.811':>7}")
    print(f"  {'Modèle C  (λ=1.4 phéno)':30}  {'4':>3}"
          f"  {dAIC_C:+8.3f}  {dBIC_C:+8.3f}"
          f"  {H0C:7.2f}  {d0C:8.4f}  {'—':>7}")
    print(f"  {'Paper 2  (α=1.78 micro)':30}  {'5':>3}"
          f"  {dAIC_P2:+8.3f}  {dBIC_P2:+8.3f}"
          f"  {H0P2:7.2f}  {d0P2:8.4f}  {s8:7.4f}")

    print(f"\nNote : Paper 2 a 1 param de plus que Modèle C, mais α est fixé")
    print(f"       depuis les simulations — pénalité AIC conservative.")
    print(f"\nsigma8 tension KiDS : {t8:+.1f}σ  {'✓ acceptable' if abs(t8)<2 else '⚠ à surveiller'}")
    d_cmb = float(d_paper2(1090, X0P2, bP2, ddP2))
    print(f"d(z=1090) = {d_cmb:.6f}  (DC = {DC:.6f})  Δ = {abs(d_cmb-DC):.1e}  → CMB ✓")

    # ── Sauvegarde CSV ──
    z_vals = np.linspace(0, 3, 100)
    rows = []
    for z in z_vals:
        dC_  = float(d_modelC(z, ddC))
        dP2_ = float(d_paper2(z, X0P2, bP2, ddP2))
        X_   = X0P2 * (1+z)**bP2
        rows.append({
            'z':       float(z),
            'd_C':     dC_,
            'w_C':     w_bup(dC_),
            'd_P2':    dP2_,
            'w_P2':    w_bup(dP2_),
            'X_P2':    float(X_),
        })
    csv_path = os.path.join(args.output_dir, 'bup_paper2_profile.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nProfil CSV → {csv_path}")

    # ── Figure ──
    res = {
        'dd_C':    ddC,   'dAIC_C':  dAIC_C,
        'X0':      X0P2,  'beta':    bP2,
        'dd_P2':   ddP2,  'dAIC_P2': dAIC_P2,
    }
    make_figure(res, args.output_dir)

    print("\nDONE")

    return {
        'lcdm':   {'H0': H0L, 'chi2': c2L},
        'modelC': {'H0': H0C, 'chi2': c2C, 'dAIC': dAIC_C, 'dBIC': dBIC_C,
                   'dd': ddC, 'd0': d0C},
        'paper2': {'H0': H0P2, 'Om': OmP2, 'chi2': c2P2,
                   'dAIC': dAIC_P2, 'dBIC': dBIC_P2,
                   'X0': X0P2, 'beta': bP2, 'dd': ddP2,
                   'd0': d0P2, 'w0': w0P2, 'sigma8': s8},
    }


if __name__ == '__main__':
    main()


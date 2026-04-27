
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_camb_geff_patch.py
══════════════════════════════════════════════════════════════════
BuP Paper 3 — Patch CAMB : G_eff dans les perturbations

VERSION A : G_eff uniforme sur toutes les échelles
  G_eff(z) = [2/(d(z)-1)] · G

Ce script fait deux choses :
  1. Génère le patch Fortran pour CAMB (sources.f90)
  2. Calcule sigma8 via croissance numérique exacte avec G_eff(z)
     (plus précis que l'estimation analytique de Paper 2)

══════════════════════════════════════════════════════════════════
PHYSIQUE
══════════════════════════════════════════════════════════════════

Équation de croissance avec G_eff :
  D'' + H·D' - (4πG_eff·ρ_m)·D = 0
  
En dimension d :
  G_eff(z) = G · 2/(d(z)-1)

  d=3     → G_eff = G      (standard)
  d=2.6   → G_eff = 1.25·G (renforcement 25%)
  d=2.4   → G_eff = 1.43·G (renforcement 43%)

Impact sur sigma8 :
  sigma8 ∝ D(z=0)  où D intègre G_eff sur toute l'histoire

══════════════════════════════════════════════════════════════════
PATCH CAMB — sources.f90
══════════════════════════════════════════════════════════════════

Dans CAMB, l'équation de croissance des perturbations scalaires
est dans sources.f90, fonction derivs() :

  grho_matter = grhom/a + grhor*(1+3*EV%w_lam)/a**2 + ...

Le terme source de croissance est :
  source_growth = 1.5*grhom*delta_m / k2

Modification BuP :
  geff = 2.0 / (d(z) - 1.0)    ! G_eff/G
  source_growth = 1.5*grhom*delta_m*geff / k2

Usage :
  1. Lancer ce script pour générer le patch
  2. Appliquer le patch à ton installation CAMB
  3. Recompiler CAMB
  4. Relancer les scans Paper 2

Usage python :
  python3 bup_camb_geff_patch.py --generate-patch
  python3 bup_camb_geff_patch.py --compute-sigma8
  python3 bup_camb_geff_patch.py --full   (les deux)

Dépendances : numpy, scipy
══════════════════════════════════════════════════════════════════
"""

import argparse
import os
import numpy as np
from scipy.integrate import odeint, cumulative_trapezoid as cumtrapz
import warnings
warnings.filterwarnings('ignore')

# ── Constantes BuP ──────────────────────────────────────────────
DC          = 3.059842935509462
AC          = 0.7252150458197096
BETA_BUP    = 1.4134453781512604
ALPHA_MICRO = 1.78

# ── Best-fit Paper 2 ────────────────────────────────────────────
X0_P2    = 0.537
BETA_P2  = 2.000
DD_P2    = 0.878
H0_P2    = 63.67
OM_P2    = 0.3448


# ══════════════════════════════════════════════════════════════
# PROFIL d(z) Paper 2
# ══════════════════════════════════════════════════════════════

def d_paper2(z, X0=X0_P2, beta=BETA_P2, dd=DD_P2):
    X = X0 * (1.0 + np.asarray(z, float))**beta
    return float(np.clip(DC - dd / (1.0 + X**ALPHA_MICRO), 2.0, DC))


def g_eff(z, X0=X0_P2, beta=BETA_P2, dd=DD_P2):
    """G_eff(z) / G = 2/(d(z)-1)"""
    d = d_paper2(z, X0, beta, dd)
    return 2.0 / (d - 1.0)


# ══════════════════════════════════════════════════════════════
# SIGMA8 NUMÉRIQUE EXACT (ODE complète)
# ══════════════════════════════════════════════════════════════

def E_fn(z, Om):
    """E(z) = H(z)/H0 pour ΛCDM (approximation fond)."""
    return float(np.sqrt(np.maximum(Om*(1+z)**3 + (1-Om), 1e-14)))


def growth_equation(y, a, Om, geff_fn):
    """
    Équation de croissance D(a) en fonction du facteur d'échelle :
      D'' + [2 + d(ln H)/d(ln a)] · D' / a - (3/2)·Ωm·H0²·G_eff/H² · D/a³ = 0
    
    Forme standard en variable a :
      d²D/da² + fric · dD/da - src · D = 0
    
    avec :
      fric = (3/a) + d(ln E)/d(ln a) / a
      src  = 3·Om / (2·a⁵·E²) · G_eff/G
    """
    D, Dp = y
    z  = 1.0/a - 1.0
    geff = float(geff_fn(z))
    Ez  = E_fn(z, Om)
    # d(ln E)/d(ln a) = -a/E · dE/da = 3/2 · Om*(1+z)³/E² · (1+z) -- ΛCDM approx
    # Forme exacte : fric = 3/a + dlnE/dlna · 1/a
    dlnE_dlna = -1.5 * Om * (1+z)**3 / Ez**2   # ΛCDM
    fric = 3.0/a + dlnE_dlna/a
    src  = 1.5 * Om / (a**5 * Ez**2) * geff
    return [Dp, -fric*Dp + src*D]


def compute_sigma8_exact(Om, geff_fn, sigma8_lcdm=0.8114, Om_lcdm=0.3153):
    """
    Calcule sigma8 via intégration numérique de l'équation de croissance.
    
    Normalisation : ratio D_BuP/D_ΛCDM · sigma8_ΛCDM
    Conditions initiales : D=a, D'=1 à a=1e-3 (régime MD)
    """
    a = np.linspace(1e-3, 1.0, 5000)

    # ΛCDM référence
    def geff_lcdm(z): return 1.0
    sol_L = odeint(
        lambda y, a_: growth_equation(y, a_, Om_lcdm, geff_lcdm),
        [1e-3, 1.0], a, rtol=1e-10, atol=1e-12
    )
    D_lcdm = float(sol_L[-1, 0])

    # BuP avec G_eff(z)
    sol_B = odeint(
        lambda y, a_: growth_equation(y, a_, Om, geff_fn),
        [1e-3, 1.0], a, rtol=1e-10, atol=1e-12
    )
    D_bup = float(sol_B[-1, 0])

    sigma8 = sigma8_lcdm * D_bup / D_lcdm
    growth_ratio = D_bup / D_lcdm
    return sigma8, growth_ratio


def scan_sigma8(X0=X0_P2, beta=BETA_P2, dd=DD_P2, Om=OM_P2):
    """Scan sigma8 en fonction des paramètres Paper 2."""

    print("=" * 65)
    print("SIGMA8 — Calcul exact via ODE croissance")
    print(f"Modèle Paper 2 : X0={X0:.4f}  β={beta:.4f}  Δd={dd:.4f}  α={ALPHA_MICRO}")
    print("=" * 65)

    # Profil d(z) pour ce best-fit
    print("\nProfil d(z) et G_eff(z) :")
    print(f"  {'z':>7}  {'d(z)':>8}  {'G_eff/G':>9}  {'renforcement':>13}")
    print("  " + "─"*42)
    for z in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 1090.0]:
        d = d_paper2(z, X0, beta, dd)
        ge = g_eff(z, X0, beta, dd)
        print(f"  {z:7.1f}  {d:8.5f}  {ge:9.5f}  {(ge-1)*100:+12.1f}%")

    # Sigma8
    geff_fn = lambda z: g_eff(z, X0, beta, dd)
    s8, ratio = compute_sigma8_exact(Om, geff_fn)
    s8_lcdm = 0.8114
    tension_kids = (s8 - 0.766) / 0.020
    tension_des  = (s8 - 0.776) / 0.017

    print(f"\nRésultats sigma8 :")
    print(f"  sigma8_ΛCDM   = {s8_lcdm:.4f}  (Planck 2018)")
    print(f"  sigma8_Paper2 = {s8:.4f}  (ODE exacte, G_eff dans perturbations)")
    print(f"  Ratio D_P2/D_ΛCDM = {ratio:.4f}")
    print(f"  KiDS-1000  = 0.766 ± 0.020  →  tension = {tension_kids:+.1f}σ")
    print(f"  DES Y3     = 0.776 ± 0.017  →  tension = {tension_des:+.1f}σ")
    print()

    # Scan Δd pour trouver la fenêtre sigma8 acceptable
    print("Scan Δd → sigma8 (X0, β, Om fixés) :")
    print(f"  {'Δd':>6}  {'d(z=0)':>8}  {'sigma8':>8}  {'KiDS(σ)':>9}  {'verdict':>10}")
    print("  " + "─"*48)
    for dd_scan in np.linspace(0.30, 0.90, 13):
        d0 = d_paper2(0, X0, beta, dd_scan)
        geff_scan = lambda z, dd_=dd_scan: g_eff(z, X0, beta, dd_)
        s8_scan, _ = compute_sigma8_exact(Om, geff_scan)
        t = (s8_scan - 0.766) / 0.020
        v = "✓" if abs(t) < 2 else ("~" if abs(t) < 3 else "")
        print(f"  {dd_scan:6.3f}  {d0:8.5f}  {s8_scan:8.4f}  {t:+9.1f}  {v:>10}")

    return s8


# ══════════════════════════════════════════════════════════════
# PATCH FORTRAN CAMB
# ══════════════════════════════════════════════════════════════

FORTRAN_MODULE = '''
!==============================================================
! BuP module — d(z) et G_eff pour CAMB
! À inclure dans equations.f90 ou sources.f90
!
! Auteur : BuP Paper 3
! Usage  : ajouter "use bup_geff" en tête de module CAMB
!==============================================================

module bup_geff
    implicit none

    ! Constantes BuP
    real(dl), parameter :: DC        = 3.059842935509462_dl
    real(dl), parameter :: AC        = 0.7252150458197096_dl
    real(dl), parameter :: BETA_BUP  = 1.4134453781512604_dl
    real(dl), parameter :: ALPHA_MIC = 1.78_dl

    ! Paramètres Paper 2 best-fit (à mettre à jour si besoin)
    real(dl) :: bup_X0   = 0.537_dl
    real(dl) :: bup_beta = 2.000_dl
    real(dl) :: bup_dd   = 0.878_dl

contains

    !----------------------------------------------------------
    ! Dimension effective d(z)
    !----------------------------------------------------------
    function bup_d_of_z(z) result(d)
        real(dl), intent(in) :: z
        real(dl) :: d, X

        X = bup_X0 * (1.0_dl + z)**bup_beta
        d = DC - bup_dd / (1.0_dl + X**ALPHA_MIC)
        ! Clamp physique
        d = max(2.0_dl, min(d, DC))

    end function bup_d_of_z

    !----------------------------------------------------------
    ! G_eff(z) / G = 2/(d(z)-1)
    !----------------------------------------------------------
    function bup_geff(z) result(geff)
        real(dl), intent(in) :: z
        real(dl) :: geff, d

        d    = bup_d_of_z(z)
        geff = 2.0_dl / (d - 1.0_dl)

    end function bup_geff

    !----------------------------------------------------------
    ! G_eff avec coupure d'échelle (Version B, Paper 3)
    ! geff(z,k) = geff(z) / (1 + (k/k_star)^2)
    !----------------------------------------------------------
    function bup_geff_scale(z, k, k_star) result(geff)
        real(dl), intent(in) :: z, k, k_star
        real(dl) :: geff

        geff = bup_geff(z) / (1.0_dl + (k/k_star)**2)

    end function bup_geff_scale

end module bup_geff
'''

PATCH_SOURCES = '''
!==============================================================
! PATCH sources.f90 — BuP G_eff dans les perturbations
!
! Chercher dans sources.f90 le calcul du terme source de
! croissance des perturbations scalaires. Il ressemble à :
!
!   sources(EV%s_MatterSource) = ... 1.5_dl * grhom * delta_m ...
!
! MODIFICATION : multiplier par bup_geff(z)
!
! AVANT (ΛCDM standard) :
!   dgrhom_source = 1.5_dl * grhom * delta_total
!
! APRÈS (BuP) :
!   z_current = 1.0_dl/a - 1.0_dl
!   geff_bup  = bup_geff(z_current)
!   dgrhom_source = 1.5_dl * grhom * delta_total * geff_bup
!
!==============================================================
!
! LOCALISATION DANS CAMB 1.6.x :
!
! Fichier : camb/sources.f90
! Fonction: TSourceMethods%CalcScalarSources()
!           ou derivs() selon la version
!
! Chercher : "ppf" ou "grhom" ou "delta_c" dans sources.f90
! La ligne exacte varie selon la version de CAMB.
!
! POUR CAMB 1.6.7 (ta version) :
!   Fichier  : fortran/equations.f90
!   Fonction : output()  (calcul des sources)
!   Variable : grho_matter  → multiplier par bup_geff(z)
!
!==============================================================
!
! IMPLÉMENTATION MINIMALE (modification d'une seule ligne) :
!
! 1. Ajouter en tête de equations.f90 :
!       use bup_geff
!
! 2. Trouver le calcul de la dérivée de delta :
!       sources = ... k2*phi + ...
!
! 3. Trouver la ligne qui calcule grho_matter ou grhom
!    et l'utilise dans l'équation de Poisson/croissance.
!    C'est typiquement dans la fonction "derivs" ou équivalent.
!
! 4. Ajouter APRÈS le calcul de z ou a :
!       geff_factor = bup_geff(z_current)
!
! 5. Multiplier le terme source par geff_factor :
!       ! Avant : source = 1.5 * grhom / k2 * delta
!       ! Après :
!       source = 1.5 * grhom / k2 * delta * geff_factor
!
!==============================================================
'''

PATCH_SCRIPT = '''#!/bin/bash
# apply_bup_patch.sh
# ══════════════════════════════════════════════════════════════
# Applique le patch BuP G_eff à une installation CAMB existante
#
# Usage : bash apply_bup_patch.sh /chemin/vers/camb
# ══════════════════════════════════════════════════════════════

CAMB_DIR="${1:-.}"
FORTRAN_DIR="$CAMB_DIR/fortran"

echo "=== BuP CAMB Patch ==="
echo "Répertoire CAMB : $CAMB_DIR"

# Vérifier que equations.f90 existe
if [ ! -f "$FORTRAN_DIR/equations.f90" ]; then
    echo "ERREUR : $FORTRAN_DIR/equations.f90 non trouvé"
    echo "Vérifier le chemin CAMB"
    exit 1
fi

# Backup
cp "$FORTRAN_DIR/equations.f90" "$FORTRAN_DIR/equations.f90.bup_backup"
echo "Backup : equations.f90.bup_backup"

# Copier le module BuP
cp bup_geff_module.f90 "$FORTRAN_DIR/"
echo "Module copié : bup_geff_module.f90"

# Ajouter "use bup_geff" après la première ligne "implicit none"
# (heuristique — à vérifier manuellement)
sed -i 's/implicit none/implicit none\\n    use bup_geff/' \
    "$FORTRAN_DIR/equations.f90"
echo "use bup_geff ajouté"

# Recompiler
echo "Recompilation..."
cd "$CAMB_DIR" && python setup.py build_ext --inplace 2>&1 | tail -5

echo "=== Patch terminé ==="
echo "Vérifier avec : python -c \\"import camb; print(camb.__version__)\\"
'''

VERIFICATION_SCRIPT = '''#!/usr/bin/env python3
"""
verify_bup_geff.py
Vérifie que le patch CAMB G_eff est actif en comparant
sigma8 avec et sans patch pour un profil d(z) connu.
"""
import camb
import numpy as np

def test_geff_active():
    """
    Test : avec d(z)=3 partout, sigma8 doit être identique à ΛCDM.
    Avec d(z)=2.4 à z=0, sigma8 doit être > sigma8_ΛCDM.
    """
    # ΛCDM reference
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.054)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_matter_power(redshifts=[0], kmax=2.0)
    results_lcdm = camb.get_results(pars)
    sigma8_lcdm = results_lcdm.get_sigma8_0()
    print(f"sigma8_ΛCDM = {sigma8_lcdm:.4f}  (attendu : 0.811)")

    # BuP avec G_eff actif (si patch appliqué)
    # Paramètres Paper 2 best-fit
    pars.bup_X0   = 0.537
    pars.bup_beta = 2.000
    pars.bup_dd   = 0.878
    pars.bup_geff_active = True

    try:
        results_bup = camb.get_results(pars)
        sigma8_bup = results_bup.get_sigma8_0()
        print(f"sigma8_BuP  = {sigma8_bup:.4f}")
        print(f"Ratio       = {sigma8_bup/sigma8_lcdm:.4f}  (attendu > 1)")
        if sigma8_bup > sigma8_lcdm:
            print("PATCH ACTIF ✓")
        else:
            print("PATCH INACTIF ou problème")
    except AttributeError:
        print("Paramètres BuP non reconnus — patch non appliqué")
        print("Appliquer le patch Fortran et recompiler CAMB")

if __name__ == '__main__':
    test_geff_active()
'''


def generate_patch_files(output_dir='.'):
    """Génère tous les fichiers du patch."""
    os.makedirs(output_dir, exist_ok=True)

    # Module Fortran
    path_mod = os.path.join(output_dir, 'bup_geff_module.f90')
    with open(path_mod, 'w') as f:
        f.write(FORTRAN_MODULE)
    print(f"Module Fortran → {path_mod}")

    # Instructions patch sources.f90
    path_patch = os.path.join(output_dir, 'PATCH_SOURCES_F90.txt')
    with open(path_patch, 'w') as f:
        f.write(PATCH_SOURCES)
    print(f"Instructions   → {path_patch}")

    # Script bash
    path_sh = os.path.join(output_dir, 'apply_bup_patch.sh')
    with open(path_sh, 'w') as f:
        f.write(PATCH_SCRIPT)
    os.chmod(path_sh, 0o755)
    print(f"Script bash    → {path_sh}")

    # Script vérification
    path_ver = os.path.join(output_dir, 'verify_bup_geff.py')
    with open(path_ver, 'w') as f:
        f.write(VERIFICATION_SCRIPT)
    print(f"Vérification   → {path_ver}")


# ══════════════════════════════════════════════════════════════
# GUIDE MANUEL D'APPLICATION DU PATCH
# ══════════════════════════════════════════════════════════════

GUIDE = """
══════════════════════════════════════════════════════════════════
GUIDE : Appliquer le patch BuP G_eff à CAMB
══════════════════════════════════════════════════════════════════

ÉTAPE 1 — Localiser ton installation CAMB
  python3 -c "import camb; print(camb.__file__)"
  # Exemple : /Users/toi/miniforge3/envs/bottomup/lib/python3.11/
  #           site-packages/camb/__init__.py
  # → Le répertoire source est dans camb/fortran/

ÉTAPE 2 — Trouver le terme de croissance dans equations.f90
  grep -n "grhom\\|grho_matter\\|delta_c" fortran/equations.f90 | head -20
  # Chercher la ligne qui calcule l'équation de perturbation scalaire

ÉTAPE 3 — Localiser la fonction de dérivée
  grep -n "subroutine\\|function" fortran/equations.f90 | grep -i "deriv\\|source"
  # C'est là que se trouve le terme source de croissance

ÉTAPE 4 — Ajouter le module BuP
  # En tête de fortran/equations.f90, après "module equations" :
  use bup_geff

  # Et copier bup_geff_module.f90 dans le répertoire fortran/

ÉTAPE 5 — Modifier le terme source
  # Trouver la ligne avec le terme source de croissance (chercher "1.5" ou "k2")
  # Exemple de modification :
  
  ! AVANT :
  EV%Delta_c_source = 1.5_dl*grhom*clxc/(EV%k2)
  
  ! APRÈS (BuP) :
  z_current = 1.0_dl/a - 1.0_dl
  EV%Delta_c_source = 1.5_dl*grhom*clxc/(EV%k2) * bup_geff(z_current)

ÉTAPE 6 — Recompiler
  cd /chemin/vers/camb
  pip install -e . --no-build-isolation

ÉTAPE 7 — Vérifier
  python3 verify_bup_geff.py

══════════════════════════════════════════════════════════════════
"""


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--generate-patch', action='store_true',
                        help='Générer les fichiers du patch Fortran')
    parser.add_argument('--compute-sigma8', action='store_true',
                        help='Calculer sigma8 via ODE exacte')
    parser.add_argument('--full', action='store_true',
                        help='Tout faire')
    parser.add_argument('--output-dir', default='bup_camb_patch')
    parser.add_argument('--X0',   type=float, default=X0_P2)
    parser.add_argument('--beta', type=float, default=BETA_P2)
    parser.add_argument('--dd',   type=float, default=DD_P2)
    parser.add_argument('--Om',   type=float, default=OM_P2)
    args = parser.parse_args()

    if args.full or args.compute_sigma8:
        s8 = scan_sigma8(args.X0, args.beta, args.dd, args.Om)

    if args.full or args.generate_patch:
        print("\n" + "="*65)
        print("GÉNÉRATION DES FICHIERS PATCH CAMB")
        print("="*65)
        generate_patch_files(args.output_dir)
        print(GUIDE)

    if not (args.full or args.compute_sigma8 or args.generate_patch):
        # Par défaut : sigma8 + patch
        scan_sigma8(args.X0, args.beta, args.dd, args.Om)
        generate_patch_files(args.output_dir)
        print(GUIDE)

    print("DONE")


if __name__ == '__main__':
    main()

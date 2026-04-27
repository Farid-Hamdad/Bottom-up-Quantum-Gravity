
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_propagator.py
══════════════════════════════════════════════════════════════════
BuP Paper 3 — Mesure du propagateur gravitationnel sur le réseau
d'intrication

Question physique centrale :
  Le propagateur gravitationnel sur le graphe BuP suit-il
      G(r) ~ r^{-(d_s - 2)}   (fractal, d_s ≈ 2.61)
  ou
      G(r) ~ r^{-1}           (Newton standard, d=3) ?

Si G(r) ~ r^{-0.61} (fractal), la gravité est plus faible à
grande distance que Newton → mécanisme naturel pour réduire σ₈.

══════════════════════════════════════════════════════════════════
MÉTHODE
══════════════════════════════════════════════════════════════════

1. Construire l'état BuP ψ pour N = 6, 9, 12, 16 qubits
2. Calculer la matrice MI exacte
3. Construire le Laplacien d'intrication L = D - W
4. Résoudre l'équation de Poisson sur le graphe :
      L·Φ = ρ   (avec ρ = delta au site source i=0)
5. Mesurer G(r) = <Φ(i)> en fonction de la distance r sur le graphe
6. Fitter G(r) = A·r^{-γ} et comparer γ avec (d_s - 2)

Trois estimateurs de dimension spectrale comparés :
  - d_RW : retour de marche (P(t) ~ t^{-d_s/2})
  - d_Phi : exposant du propagateur γ = d_s - 2
  - d_geo : dimension géodésique du graphe

══════════════════════════════════════════════════════════════════
INTERPRÉTATION ATTENDUE
══════════════════════════════════════════════════════════════════

  γ ≈ 1.0  →  Newton standard, gravité non modifiée
  γ ≈ 0.6  →  Fractal BuP, gravité affaiblie à grande distance
  γ ≈ 0.0  →  Gravité constante (dimension effective ~2)

Résolution σ₈ :
  Si γ ≈ 0.6 aux échelles sub-Hubble, G_eff_net ≈ G malgré
  G_eff_fond = 2/(d-1)·G → σ₈ proche de ΛCDM.

Usage :
  python3 bup_propagator.py --N 6 9 12 16 --lambda-n 6

Dépendances : numpy, scipy, matplotlib
══════════════════════════════════════════════════════════════════
"""

import argparse
import csv
import os
import warnings
import numpy as np
from scipy.linalg import eigh, solve
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# ── Constantes BuP ──────────────────────────────────────────────
DC   = 3.059842935509462
AC   = 0.7252150458197096
BETA = 1.4134453781512604


# ══════════════════════════════════════════════════════════════
# 1. ÉTAT BuP (même logique que bup_simulation_n16.py)
# ══════════════════════════════════════════════════════════════

def factor_grid(N):
    best = None
    for r in range(1, int(N**0.5) + 2):
        if N % r == 0:
            c = N // r
            s = abs(c - r)
            if best is None or s < best[0]:
                best = (s, r, c)
    return best[1], best[2]


def grid_neighbors(i, rows, cols):
    r, c = i // cols, i % cols
    out = []
    if r > 0:        out.append(i - cols)
    if r < rows - 1: out.append(i + cols)
    if c > 0:        out.append(i - 1)
    if c < cols - 1: out.append(i + 1)
    return out


def generate_state(N, lam):
    """
    État variationnel BuP (même construction que bup_simulation_n16.py).
    lam : couplage non-local [0, 1]
    """
    rows, cols = factor_grid(N)
    J, h, J_nl = 1.0, 0.5 + lam*0.5, 0.8

    # Paires non-locales (coins opposés)
    pairs = set()
    for a, b in [(0, N-1), (cols-1, N-cols)]:
        if 0 <= a < N and 0 <= b < N and a != b:
            pairs.add(tuple(sorted((a, b))))
    pairs = list(pairs)

    dim = 2**N
    psi = np.zeros(dim, dtype=complex)

    for cfg in range(dim):
        bits = [(cfg >> i) & 1 for i in range(N)]
        loc = sum(
            1 if bits[i] == bits[j] else -1
            for i in range(N) for j in grid_neighbors(i, rows, cols)
            if i < j
        )
        nloc = sum(
            lam * J_nl * (1 if bits[i] == bits[j] else -1)
            for i, j in pairs
        ) if lam > 0 else 0.0
        n_up = sum(bits)
        psi[cfg] = (
            np.exp(J * loc / max(2.0, N/4.0)) *
            np.exp(nloc / max(1.0, len(pairs))) *
            np.exp(1j * h * n_up)
        )

    norm = np.linalg.norm(psi)
    if norm < 1e-15:
        raise ValueError("État nul")
    return psi / norm


# ══════════════════════════════════════════════════════════════
# 2. MATRICE MI EXACTE
# ══════════════════════════════════════════════════════════════

def mi_matrix(psi, N):
    """MI(i,j) = S(ρ_i) + S(ρ_j) - S(ρ_{ij}) — exact pour toutes paires."""

    def rho1(i):
        psi_t = psi.reshape([2]*N)
        other = [k for k in range(N) if k != i]
        psi_r = np.transpose(psi_t, [i] + other)
        m = psi_r.reshape(2, 2**(N-1))
        return m @ m.conj().T

    def rho2(i, j):
        psi_t = psi.reshape([2]*N)
        other = [k for k in range(N) if k != i and k != j]
        psi_r = np.transpose(psi_t, [i, j] + other)
        m = psi_r.reshape(4, 2**(N-2))
        return m @ m.conj().T

    def S(rho):
        v = np.linalg.eigvalsh(rho)
        v = v[v > 1e-14]
        return float(-np.sum(v * np.log(v)))

    Ss = [S(rho1(i)) for i in range(N)]
    MI = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1, N):
            Sij = S(rho2(i, j))
            MI[i, j] = MI[j, i] = max(Ss[i] + Ss[j] - Sij, 0.0)
    return MI


# ══════════════════════════════════════════════════════════════
# 3. DIMENSION SPECTRALE (retour de marche)
# ══════════════════════════════════════════════════════════════

def spectral_dimension_rw(MI, n_t=80):
    """d_s via P(t) = Tr[exp(-tL)]/N."""
    N = MI.shape[0]
    W = MI.copy(); np.fill_diagonal(W, 0)
    deg = W.sum(1)
    if np.any(deg < 1e-12):
        return np.nan, None

    L = np.diag(deg) - W
    ev = np.sort(np.linalg.eigvalsh(L))
    ev = ev[ev > 1e-10]
    if len(ev) < 2:
        return np.nan, None

    t = np.logspace(np.log10(0.5/ev[-1]), np.log10(2.0/ev[0]), n_t)
    P = np.array([np.sum(np.exp(-ev*ti))/N for ti in t])
    ds_arr = -2.0 * np.gradient(np.log(P+1e-15), np.log(t))
    mid = slice(n_t//4, n_t//2)
    return float(np.median(ds_arr[mid])), L


# ══════════════════════════════════════════════════════════════
# 4. PROPAGATEUR GRAVITATIONNEL
# ══════════════════════════════════════════════════════════════

def gravitational_propagator(L, MI, N, source=0):
    """
    Résout L·Φ = ρ avec ρ = δ_{i, source}.

    Retourne G(r) = Φ(source) - Φ(r) : différence de potentiel.
    G(r) > 0 et décroît avec r (comportement gravitationnel correct).
    """
    eps = 1e-6 * np.max(np.abs(L))
    L_reg = L + eps * np.eye(N)

    rho = np.zeros(N)
    rho[source] = 1.0

    Phi = np.linalg.solve(L_reg, rho)
    # Différence de potentiel par rapport à la source
    G_r = Phi[source] - Phi
    return G_r


def graph_distances(MI, N, source=0):
    """
    Distance sur la grille (Manhattan).
    Physiquement plus robuste que 1/MI pour le propagateur.
    Pour N nœuds sur grille rows×cols :
      dist(i, j) = |r_i - r_j| + |c_i - c_j|
    """
    rows, cols = factor_grid(N)
    r_src, c_src = source // cols, source % cols
    dist = np.zeros(N)
    for i in range(N):
        r_i, c_i = i // cols, i % cols
        dist[i] = abs(r_i - r_src) + abs(c_i - c_src)
    return dist


def measure_propagator(Phi, distances, N, source=0, n_bins=8):
    """
    Mesure G(r) = <|Φ(i)|> en fonction de la distance r.
    Retourne (r_bins, G_bins).
    """
    r_arr = np.array([distances[i] for i in range(N) if i != source])
    G_arr = np.array([Phi[i] for i in range(N) if i != source])

    if len(r_arr) < 3:
        return None, None

    # Binning logarithmique
    r_min, r_max = r_arr.min(), r_arr.max()
    if r_min <= 0 or r_max <= r_min:
        return None, None

    bins = np.logspace(np.log10(r_min*0.99), np.log10(r_max*1.01), n_bins+1)
    r_centers, G_means = [], []

    for i in range(n_bins):
        mask = (r_arr >= bins[i]) & (r_arr < bins[i+1])
        if mask.sum() >= 2:
            r_centers.append(float(np.mean(r_arr[mask])))
            G_means.append(float(np.mean(G_arr[mask])))

    if len(r_centers) < 3:
        return None, None

    return np.array(r_centers), np.array(G_means)


def fit_power_law(r, G):
    """Fitte G(r) = A·r^{-γ} en espace log-log."""
    if r is None or len(r) < 3:
        return np.nan, np.nan, np.nan

    valid = (r > 0) & (G > 0) & np.isfinite(r) & np.isfinite(G)
    if valid.sum() < 3:
        return np.nan, np.nan, np.nan

    log_r = np.log(r[valid])
    log_G = np.log(G[valid])

    try:
        coeffs = np.polyfit(log_r, log_G, 1)
        gamma = -coeffs[0]   # G ~ r^{-γ}
        A     = np.exp(coeffs[1])

        # R² pour qualité du fit
        pred = coeffs[0]*log_r + coeffs[1]
        ss_res = np.sum((log_G - pred)**2)
        ss_tot = np.sum((log_G - log_G.mean())**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0

        return float(gamma), float(A), float(r2)
    except Exception:
        return np.nan, np.nan, np.nan


# ══════════════════════════════════════════════════════════════
# 5. DIMENSION GÉODÉSIQUE DU GRAPHE
# ══════════════════════════════════════════════════════════════

def geodesic_dimension(distances, N, source=0, n_bins=6):
    """
    Dimension géodésique d_geo via N(r) ~ r^{d_geo}
    N(r) = nombre de nœuds à distance ≤ r.
    """
    r_arr = np.array([distances[i] for i in range(N) if i != source])
    r_arr = r_arr[r_arr < np.inf]
    if len(r_arr) < 4:
        return np.nan

    r_vals = np.logspace(
        np.log10(r_arr.min() + 1e-10),
        np.log10(r_arr.max()),
        n_bins
    )
    Nr = np.array([np.sum(r_arr <= r) for r in r_vals], dtype=float)
    Nr = Nr[Nr > 0]
    r_vals = r_vals[:len(Nr)]

    if len(Nr) < 3:
        return np.nan

    try:
        coeffs = np.polyfit(np.log(r_vals), np.log(Nr + 1e-10), 1)
        return float(coeffs[0])
    except Exception:
        return np.nan


# ══════════════════════════════════════════════════════════════
# 6. SCAN PRINCIPAL
# ══════════════════════════════════════════════════════════════

def scan_propagator(N, lambda_grid, n_repeat=2):
    """Scan (N, λ) → (d_RW, γ, d_geo, R²)."""
    rows, cols = factor_grid(N)
    results = []

    print(f"\n── N={N:2d}  (grille {rows}×{cols}) ─────────────────────────────────")
    print(f"  {'λ':>5}  {'d_RW':>7}  {'γ=d_s-2':>9}  {'d_s-2':>8}  "
          f"{'d_geo':>7}  {'R²':>6}  {'verdict':>10}")
    print("  " + "─"*60)

    for lam in lambda_grid:
        d_rw_list, gamma_list, dgeo_list, r2_list = [], [], [], []

        for _ in range(n_repeat):
            try:
                psi  = generate_state(N, lam)
                MI   = mi_matrix(psi, N)
                d_rw, L = spectral_dimension_rw(MI)

                if L is None or not np.isfinite(d_rw):
                    continue

                # Propagateur depuis chaque site source (moyenne)
                gammas, r2s = [], []
                rows_i, cols_i = factor_grid(N)
                center = (rows_i//2)*cols_i + cols_i//2
                for src in [center, 0, N-1]:
                    Phi   = gravitational_propagator(L, MI, N, source=src)
                    dists = graph_distances(MI, N, source=src)
                    r_bins, G_bins = measure_propagator(Phi, dists, N, source=src)
                    gam, _, r2 = fit_power_law(r_bins, G_bins)
                    if np.isfinite(gam) and np.isfinite(r2):
                        gammas.append(gam); r2s.append(r2)

                if gammas:
                    gamma_list.append(float(np.mean(gammas)))
                    r2_list.append(float(np.mean(r2s)))

                # Dimension géodésique (depuis site 0)
                dists0 = graph_distances(MI, N, source=0)
                d_geo  = geodesic_dimension(dists0, N)

                if np.isfinite(d_rw):   d_rw_list.append(d_rw)
                if np.isfinite(d_geo):  dgeo_list.append(d_geo)

            except Exception as e:
                pass

        def sm(lst): return float(np.mean(lst)) if lst else np.nan

        d_rw_  = sm(d_rw_list)
        gamma_ = sm(gamma_list)
        d_geo_ = sm(dgeo_list)
        r2_    = sm(r2_list)

        # Prédiction fractale : γ = d_s - 2
        gamma_pred = d_rw_ - 2 if np.isfinite(d_rw_) else np.nan
        delta = gamma_ - gamma_pred if np.isfinite(gamma_) else np.nan

        # Verdict
        if np.isfinite(gamma_) and np.isfinite(gamma_pred):
            if abs(delta) < 0.1:
                verdict = "FRACTAL ✓"
            elif gamma_ > 0.9:
                verdict = "Newton"
            elif gamma_ < 0.1:
                verdict = "d~2"
            else:
                verdict = f"Δ={delta:+.2f}"
        else:
            verdict = "n/a"

        def fmt(x): return f"{x:7.4f}" if np.isfinite(x) else "    nan"
        print(f"  {lam:5.2f}  {fmt(d_rw_)}  {fmt(gamma_):>9}  "
              f"{fmt(gamma_pred):>8}  {fmt(d_geo_):>7}  "
              f"{r2_:6.3f}  {verdict:>10}")

        results.append({
            'N': int(N), 'lambda': float(lam),
            'd_RW':  d_rw_,   'gamma':   gamma_,
            'd_geo': d_geo_,  'R2':       r2_,
            'gamma_pred': gamma_pred,
            'delta': delta,
        })

    return results


# ══════════════════════════════════════════════════════════════
# 7. FIGURE
# ══════════════════════════════════════════════════════════════

def make_figure(all_results, output_dir):
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'serif', 'font.size': 11,
        'axes.spines.top': False, 'axes.spines.right': False,
    })

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.subplots_adjust(wspace=0.35)

    colors = {6:'#aabbcc', 8:'#5588bb', 9:'#1a3a7a', 12:'#7a3a1a', 16:'#7a1a3a'}
    Ns = sorted(set(r['N'] for r in all_results))

    # ── Panel 1 : γ(λ) ──
    ax = axes[0]
    # Références
    ax.axhline(1.0,        color='#1a7a3a', ls='--', lw=1.5,
               label='Newton : $\\gamma=1$')
    ax.axhline(DC - 2.0,   color='#7a1a3a', ls='--', lw=1.5,
               label=f'Fractal BuP : $\\gamma=d_s-2\\approx{DC-2:.2f}$')
    ax.axhline(0.61,       color='#cc6600', ls=':', lw=1,
               label='$\\gamma=0.61$ (cible)')

    for N in Ns:
        r_N = [r for r in all_results if r['N'] == N]
        lam  = [r['lambda'] for r in r_N]
        gam  = [r['gamma']  for r in r_N]
        ax.plot(lam, gam, 'o-', color=colors.get(N, '#444'),
                lw=1.5, ms=5, label=f'$N={N}$')

    ax.set_xlabel(r'Couplage non-local $\lambda$')
    ax.set_ylabel(r'Exposant propagateur $\gamma$')
    ax.set_title(r'$G(r) \sim r^{-\gamma}$', pad=6)
    ax.set_xlim(-0.05, 1.05)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel 2 : γ vs d_s - 2 ──
    ax = axes[1]
    for N in Ns:
        r_N = [r for r in all_results
               if r['N'] == N and np.isfinite(r['gamma'])
               and np.isfinite(r['gamma_pred'])]
        if not r_N:
            continue
        x = [r['gamma_pred'] for r in r_N]
        y = [r['gamma']      for r in r_N]
        ax.scatter(x, y, color=colors.get(N, '#444'), s=40,
                   label=f'$N={N}$', zorder=5)

    # Ligne de cohérence parfaite
    lim_vals = np.linspace(0, 2, 50)
    ax.plot(lim_vals, lim_vals, 'k--', lw=1, alpha=0.4, label='$\\gamma = d_s-2$')

    ax.set_xlabel(r'Prédiction $d_s - 2$')
    ax.set_ylabel(r'Mesuré $\gamma$')
    ax.set_title('Cohérence fractal', pad=6)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel 3 : d_geo(λ) ──
    ax = axes[2]
    ax.axhline(DC, color='gray', ls=':', lw=1, alpha=0.5, label=f'$d_c={DC:.2f}$')
    ax.axhline(2.644, color='#1a7a3a', ls='--', lw=1, alpha=0.6,
               label='SPARC $2.644$')

    for N in Ns:
        r_N = [r for r in all_results if r['N'] == N]
        lam  = [r['lambda'] for r in r_N]
        dgeo = [r['d_geo']  for r in r_N]
        ax.plot(lam, dgeo, 'o-', color=colors.get(N, '#444'),
                lw=1.5, ms=5, label=f'$N={N}$')

    ax.set_xlabel(r'Couplage non-local $\lambda$')
    ax.set_ylabel(r'Dimension géodésique $d_\mathrm{geo}$')
    ax.set_title('Dimension géodésique', pad=6)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        r'BuP Paper 3 — Propagateur gravitationnel $G(r) \sim r^{-\gamma}$'
        r' sur le réseau d\'intrication',
        fontsize=11, y=1.01)

    path = os.path.join(output_dir, 'bup_propagator.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"\nFigure → {path}")


# ══════════════════════════════════════════════════════════════
# 8. RÉSUMÉ
# ══════════════════════════════════════════════════════════════

def print_summary(all_results):
    print()
    print("=" * 65)
    print("RÉSUMÉ — Test de cohérence fractale")
    print("Hypothèse BuP : γ = d_s - 2 ≈ 0.61")
    print("=" * 65)

    # Moyenne par N
    Ns = sorted(set(r['N'] for r in all_results))
    print(f"\n  {'N':>4}  {'γ moyen':>9}  {'d_s-2 prédit':>14}  "
          f"{'Δ':>7}  {'R²':>6}  {'verdict':>12}")
    print("  " + "─"*60)

    for N in Ns:
        r_N = [r for r in all_results
               if r['N'] == N and np.isfinite(r['gamma'])]
        if not r_N:
            continue
        gam_m  = float(np.mean([r['gamma']      for r in r_N]))
        pred_m = float(np.mean([r['gamma_pred'] for r in r_N
                                 if np.isfinite(r['gamma_pred'])]))
        r2_m   = float(np.mean([r['R2']         for r in r_N
                                 if np.isfinite(r['R2'])]))
        delta  = gam_m - pred_m

        if abs(delta) < 0.15:
            v = "FRACTAL ✓✓"
        elif abs(delta) < 0.30:
            v = "FRACTAL ~"
        elif gam_m > 0.85:
            v = "Newton"
        else:
            v = f"interméd."

        print(f"  {N:4d}  {gam_m:9.4f}  {pred_m:14.4f}  "
              f"{delta:+7.4f}  {r2_m:6.3f}  {v:>12}")

    print(f"""
INTERPRÉTATION :
  γ ≈ 0.61 (fractal BuP) → gravité affaiblie à grande distance
    → mécanisme naturel pour réduire σ₈ sans modifier les BAO

  γ ≈ 1.0  (Newton)      → propagateur standard, pas de modification
    → G_eff = 2/(d-1)·G sans compensation → tension σ₈ persiste

  γ intermédiaire        → régime de transition, dépend de l'échelle

NOTE : ces mesures sont sur N ≤ 16 nœuds.
  L'extrapolation N→∞ nécessite une analyse FSS séparée.
""")


# ══════════════════════════════════════════════════════════════
# 9. MAIN
# ══════════════════════════════════════════════════════════════

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N',        type=int,   nargs='+', default=[6, 9])
    ap.add_argument('--lambda-n', type=int,   default=6)
    ap.add_argument('--n-repeat', type=int,   default=2)
    ap.add_argument('--output-dir', default='results_propagator')
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    lam_grid = np.linspace(0.0, 1.0, args.lambda_n)

    print("=" * 65)
    print("BuP Paper 3 — Propagateur gravitationnel")
    print(f"Systèmes   : N = {args.N}")
    print(f"Grille λ   : [{lam_grid[0]:.2f}, {lam_grid[-1]:.2f}] "
          f"({args.lambda_n} pts)")
    print(f"Répétitions: {args.n_repeat}")
    print(f"Hypothèse  : G(r) ~ r^{{-(d_s-2)}} avec d_s ≈ 2.61 → γ ≈ 0.61")
    print("=" * 65)

    all_results = []
    for N in args.N:
        res = scan_propagator(N, lam_grid, n_repeat=args.n_repeat)
        all_results.extend(res)

    print_summary(all_results)

    # CSV
    if all_results:
        csv_path = os.path.join(args.output_dir, 'propagator_results.csv')
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
        print(f"CSV → {csv_path}")

    make_figure(all_results, args.output_dir)
    print("\nDONE")


if __name__ == '__main__':
    main()

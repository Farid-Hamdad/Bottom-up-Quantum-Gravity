#!/usr/bin/env python3
"""
BuP Paper 10 — Convergence inter-centres sur graphe MI avec d_s → 3

But : vérifier si α_médian converge vers d_s − 2 ≈ 1 quand d_s → 3.

Méthode :
  - Charger la matrice MI réelle (CSV)
  - Construire le graphe kNN sparse
  - Résoudre L·Φ = δ_{i0} − 1/N par CG pour un sous-échantillon de centres
  - Fitter Φ(r) = C + A·r^{−α} sur données brutes
  - Rapporter α_médian, σ_α, CV, |Δα| vs α_pred = d_s − 2
  - Deux seuils R² : >0.70 et >0.85

Usage :
    python bup_paper10_mi_ds3_convergence.py \
        --mi-file path/to/MI_matrix.csv \
        --k 5 \
        --n-centers 50 \
        --r2-thresholds 0.70 0.85 \
        --output-dir results/mi_ds3_convergence

    # Pour scanner plusieurs k :
    for k in 5 10 15; do
        python bup_paper10_mi_ds3_convergence.py \
            --mi-file MI_ds3.csv --k $k \
            --output-dir results/mi_ds3_k${k}
    done
"""

import os
import re
import json
import argparse
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import cg
from scipy.sparse.csgraph import shortest_path
from scipy.optimize import curve_fit
from scipy.linalg import expm


# ─── Utilitaires ──────────────────────────────────────────

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def parse_lambda(path):
    """Reconnaît 'lam0.3', 'lam0_3', 'lam0p3' → 0.3."""
    m = re.search(r"lam([0-9]+(?:[._p][0-9]+)?)", os.path.basename(path))
    if not m:
        return None
    return float(m.group(1).replace('_', '.').replace('p', '.'))

def load_mi(path):
    """Charge une matrice MI carrée depuis un CSV."""
    try:
        W = np.loadtxt(path, delimiter=",")
        if W.ndim == 2 and W.shape[0] == W.shape[1]:
            pass
        else:
            raise ValueError
    except Exception:
        df = pd.read_csv(path, header=None)
        W = df.values.astype(float)
        if not np.issubdtype(W.dtype, np.number):
            df = pd.read_csv(path, index_col=0)
            W = df.values.astype(float)
    W = np.nan_to_num(W, nan=0., posinf=0., neginf=0.)
    W = 0.5*(W + W.T); np.fill_diagonal(W, 0.); W[W<0] = 0.
    return W

def build_sparse_laplacian(W, k):
    """Construit L sparse kNN depuis la matrice MI pleine."""
    N = W.shape[0]
    mask = np.zeros((N,N), dtype=bool)
    for i in range(N):
        idx = np.argsort(W[i])[::-1]
        idx = [j for j in idx if j != i and W[i,j] > 0]
        for j in idx[:k]: mask[i,j] = True
    mask = np.logical_or(mask, mask.T); np.fill_diagonal(mask, False)
    A = np.where(mask, W, 0.); A = 0.5*(A+A.T); np.fill_diagonal(A, 0.)
    # Sparse Laplacian
    deg = A.sum(1)
    rr, cc, vv = list(range(N)), list(range(N)), list(deg)
    for i in range(N):
        for j in np.where(A[i]>0)[0]:
            rr.append(i); cc.append(int(j)); vv.append(-A[i,j])
    Lsp = sparse.csr_matrix((vv,(rr,cc)), shape=(N,N))
    return Lsp, A, mask

def entanglement_distance(A, mask, eps=1e-12):
    """Distance d_ent = −log(W/W_max) sur les arêtes kNN."""
    Wmax = A[mask].max()
    L = np.full_like(A, np.inf)
    L[mask] = -np.log((A[mask]+eps)/(Wmax+eps))
    L[mask] = np.maximum(L[mask], 1e-9)
    np.fill_diagonal(L, 0.)
    D = shortest_path(L, method="FW", directed=False)
    D[np.isinf(D)] = np.nan
    return D

def spectral_dim_eigenvalue(A, N):
    """d_s via ratio λ_max / λ_2 du Laplacien."""
    deg = A.sum(1); L = np.diag(deg) - A
    ev = np.linalg.eigvalsh(L); ep = ev[ev > 1e-10]
    if len(ep) < 2: return np.nan
    return 2*np.log(N) / (np.log(ep[-1]/ep[0]) + 1e-12)

def spectral_dim_heat(A, N, tau_min=0.01, tau_max=50, tau_pts=40):
    """d_s via noyau de chaleur (méthode Farid)."""
    deg = A.sum(1); L = np.diag(deg) - A
    taus = np.logspace(np.log10(tau_min), np.log10(tau_max), tau_pts)
    Z = np.array([np.trace(expm(-t*L)) for t in taus])
    logt = np.log(taus); logZ = np.log(np.maximum(Z, 1e-300))
    slope = np.gradient(logZ, logt); ds = -2.*slope
    n = len(ds); lo = max(1,n//4); hi = min(n-1,3*n//4)
    return taus, ds, float(np.nanmedian(ds[lo:hi]))

def solve_point_source(Lsp, center, N, rtol=1e-8, atol=1e-10, maxiter=30000):
    """L·Φ = δ_{center} − 1/N par CG sparse."""
    b = np.full(N, -1.0/N); b[center] += 1.0
    try:
        phi, info = cg(Lsp, b, rtol=rtol, atol=atol, maxiter=maxiter)
    except TypeError:
        phi, info = cg(Lsp, b, tol=rtol, maxiter=maxiter)
    phi = np.asarray(phi, dtype=float)
    return phi - phi.mean(), info

def fit_alpha_raw(D_or_r, phi, center):
    """
    Fit Φ(r) = C + A·r^{−α} sur toutes les paires (center, j).

    D_or_r peut être :
      - soit la matrice complète D, auquel cas on prend D[center]
      - soit directement le vecteur de distances r_use
    """
    if np.ndim(D_or_r) == 1:
        r = D_or_r
    else:
        r = D_or_r[center]

    y = phi

    msk = np.isfinite(r) & np.isfinite(y) & (r > 0)
    r_m, y_m = r[msk], y[msk]
    if len(r_m) < 5: return np.nan, np.nan, np.nan, np.nan

    def m(r, C, A, a): return C + A*r**(-a)

    # Newton fixé (α=1) pour référence
    X = np.column_stack([np.ones_like(r_m), 1./r_m])
    coef, *_ = np.linalg.lstsq(X, y_m, rcond=None)
    C_N, A_N = coef; yh_N = C_N + A_N/r_m
    R2_N = 1 - np.sum((y_m-yh_N)**2)/(np.sum((y_m-y_m.mean())**2)+1e-15)

    # Libre
    best_a, best_R2 = np.nan, -np.inf
    for a0 in [-1.5, -0.9, -0.5, 0., 0.5, 1.0, 1.5]:
        try:
            C0 = float(np.median(y_m)); A0 = y_m[0] - C0
            p, _ = curve_fit(m, r_m, y_m, p0=[C0, A0, a0],
                             bounds=([-np.inf,-np.inf,-5],[np.inf,np.inf,5]),
                             maxfev=5000)
            yh = m(r_m, *p)
            R2 = 1 - np.sum((y_m-yh)**2)/(np.sum((y_m-y_m.mean())**2)+1e-15)
            if np.isfinite(R2) and R2 > best_R2:
                best_a, best_R2 = float(p[2]), R2
        except: pass

    return best_a, best_R2, float(R2_N), int(len(r_m))


# ─── Main ─────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="BuP Paper 10 — Convergence inter-centres sur graphe MI ds→3")
    parser.add_argument("--mi-file",   required=True,
                        help="Matrice MI carrée (CSV)")
    parser.add_argument("--k",         type=int,   default=5,
                        help="Nombre de voisins kNN")
    parser.add_argument("--n-centers", type=int,   default=50,
                        help="Nombre de centres à tester (sous-échantillon)")
    parser.add_argument("--seed-centers", type=int, default=0,
                        help="Seed pour le sous-échantillon de centres")
    parser.add_argument("--r2-thresholds", type=float, nargs="+",
                        default=[0.70, 0.85],
                        help="Seuils R² pour filtrer les centres")
    parser.add_argument("--r2-fit-min",    type=float, default=None,
                        help="Fenêtre fit : fraction min de r_max (ex. 0.1)")
    parser.add_argument("--r2-fit-max",    type=float, default=None,
                        help="Fenêtre fit : fraction max de r_max (ex. 0.4)")
    parser.add_argument("--output-dir",    default="results/mi_ds3_convergence")
    parser.add_argument("--tau-min",       type=float, default=0.01)
    parser.add_argument("--tau-max",       type=float, default=50.)
    parser.add_argument("--tau-points",    type=int,   default=40)
    parser.add_argument("--cg-rtol",       type=float, default=1e-8)
    parser.add_argument("--cg-maxiter",    type=int,   default=30000)
    args = parser.parse_args()
    ensure_dir(args.output_dir)

    t0 = time.time()
    lam = parse_lambda(args.mi_file)

    print("="*70)
    print("BuP Paper 10 — MI point-source convergence inter-centres")
    print("="*70)
    print(f"Fichier MI : {args.mi_file}")
    print(f"λ parsé    : {lam}")

    # ─── Chargement + graphe ─────────────────
    W = load_mi(args.mi_file)
    N = W.shape[0]
    print(f"N          : {N}")
    print(f"k          : {args.k}")

    Lsp, A, mask = build_sparse_laplacian(W, args.k)
    D             = entanglement_distance(A, mask)

    # ─── Dimension spectrale (deux méthodes) ─
    d_s_eig = spectral_dim_eigenvalue(A, N)
    taus, ds_curve, d_s_heat = spectral_dim_heat(
        A, N, args.tau_min, args.tau_max, args.tau_points)
    alpha_pred = d_s_heat - 2.

    print(f"d_s (eigenvalue ratio) : {d_s_eig:.4f}")
    print(f"d_s (heat kernel)      : {d_s_heat:.4f}")
    print(f"α_pred = d_s − 2       : {alpha_pred:.4f}")
    print(f"Cible Newton (α=1)     : d_s = 3 → α_pred = 1")

    # ─── Sous-échantillon de centres ─────────
    rng   = np.random.default_rng(args.seed_centers)
    n_ctr = min(args.n_centers, N)
    ctrs  = rng.choice(N, size=n_ctr, replace=False)
    print(f"\nCentres testés : {n_ctr} / {N}")

    # ─── Boucle centres ──────────────────────
    records = []
    n_cg_warn = 0

    for c in ctrs:
        phi, info = solve_point_source(
            Lsp, c, N, rtol=args.cg_rtol, maxiter=args.cg_maxiter)
        if info != 0: n_cg_warn += 1

        # Fenêtre optionnelle
        if args.r2_fit_min is not None or args.r2_fit_max is not None:
            r_full = D[c].copy()
            r_finite = r_full[np.isfinite(r_full) & (r_full > 0)]
            r_max_val = r_finite.max() if len(r_finite) else 1.
            lo = (args.r2_fit_min or 0.) * r_max_val
            hi = (args.r2_fit_max or 1.) * r_max_val
            msk_w = (D[c] >= lo) & (D[c] <= hi)
            r_use = D[c].copy(); r_use[~msk_w] = np.nan
        else:
            r_use = D[c]

        # Note : on utilise le même phi pour tous les appels à fit_alpha_raw
        # phi a la même longueur que N, on utilise r_use comme premier argument
        alpha, R2_f, R2_N, n_pairs = fit_alpha_raw(r_use, phi, c)
        records.append({'center':int(c), 'alpha':alpha,
                        'R2_free':R2_f, 'R2_newton':R2_N, 'n_pairs':n_pairs})

    if n_cg_warn:
        print(f"[WARN] CG non convergé pour {n_cg_warn}/{n_ctr} centres")

    # ─── Statistiques par seuil R² ───────────
    results_by_thresh = {}
    for thresh in args.r2_thresholds:
        sel = [r for r in records
               if np.isfinite(r['alpha']) and np.isfinite(r['R2_free'])
               and r['R2_free'] >= thresh]
        if len(sel) < 3:
            print(f"[WARN] Trop peu de centres pour R²>{thresh}: n={len(sel)}")
            results_by_thresh[thresh] = None
            continue
        alphas = np.array([s['alpha'] for s in sel])
        results_by_thresh[thresh] = {
            'n':       int(len(sel)),
            'median':  float(np.median(alphas)),
            'mean':    float(np.mean(alphas)),
            'std':     float(np.std(alphas)),
            'iqr':     float(np.percentile(alphas,75) - np.percentile(alphas,25)),
            'cv':      float(np.std(alphas)/(np.abs(np.mean(alphas))+1e-12)),
            'delta':   float(np.median(alphas) - alpha_pred),
            'q05':     float(np.percentile(alphas, 5)),
            'q25':     float(np.percentile(alphas,25)),
            'q75':     float(np.percentile(alphas,75)),
            'q95':     float(np.percentile(alphas,95)),
        }

    # ─── Affichage ───────────────────────────
    print(f"\n{'─'*70}")
    print(f"RÉSULTATS — N={N}, k={args.k}, d_s={d_s_heat:.3f}, α_pred={alpha_pred:.3f}")
    print(f"{'─'*70}")
    print(f"{'Seuil R²':>10} {'n':>5} {'α méd':>8} {'α moy':>8} "
          f"{'σ':>7} {'IQR':>7} {'CV':>7} {'|Δα|':>8}")
    print("-"*65)
    for thresh in args.r2_thresholds:
        s = results_by_thresh.get(thresh)
        if s is None:
            print(f"{thresh:>10.2f}  trop peu de points")
            continue
        print(f"{thresh:>10.2f} {s['n']:5d} {s['median']:8.4f} "
              f"{s['mean']:8.4f} {s['std']:7.4f} {s['iqr']:7.4f} "
              f"{s['cv']:7.3f} {abs(s['delta']):8.4f}")

    # ─── Outputs ─────────────────────────────
    # CSV brut centres
    df_rec = pd.DataFrame(records)
    df_rec.to_csv(os.path.join(args.output_dir, "centers_alpha.csv"), index=False)

    # CSV spectral dimension
    pd.DataFrame({'tau':taus,'d_s':ds_curve}).to_csv(
        os.path.join(args.output_dir, "spectral_dim_curve.csv"), index=False)

    # JSON summary
    summary = {
        'mi_file':    args.mi_file,
        'lambda':     lam,
        'N':          N,
        'k':          args.k,
        'n_centers':  n_ctr,
        'd_s_eigenvalue': float(d_s_eig),
        'd_s_heat':       float(d_s_heat),
        'alpha_pred':     float(alpha_pred),
        'results_by_threshold': {
            str(t): results_by_thresh[t] for t in args.r2_thresholds
        },
        'elapsed_s': round(time.time()-t0, 1)
    }
    with open(os.path.join(args.output_dir, "summary.json"), 'w') as f:
        json.dump(summary, f, indent=2)

    # CSV résumé
    rows_sum = []
    for thresh, s in results_by_thresh.items():
        if s is None: continue
        rows_sum.append({'mi_file': os.path.basename(args.mi_file),
                         'N':N, 'k':args.k, 'lambda':lam,
                         'd_s_heat':d_s_heat, 'alpha_pred':alpha_pred,
                         'R2_threshold':thresh, **s})
    pd.DataFrame(rows_sum).to_csv(
        os.path.join(args.output_dir, "summary.csv"), index=False)

    # ─── Figures ─────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle(f"BuP Paper 10 — MI ds→3  N={N}, k={args.k}, "
                 f"d_s={d_s_heat:.3f}, α_pred={alpha_pred:.3f}", fontsize=11)

    # (0) Distribution α par seuil R²
    ax = axes[0]
    colors = ['steelblue', 'tomato', 'mediumseagreen']
    for i, thresh in enumerate(args.r2_thresholds):
        sel = [r['alpha'] for r in records
               if np.isfinite(r['alpha']) and r.get('R2_free',0) >= thresh]
        if sel:
            ax.hist(sel, bins=20, alpha=0.6, color=colors[i%3],
                    label=f'R²>{thresh} (n={len(sel)})', density=True)
    ax.axvline(alpha_pred, color='r', ls='--', lw=2, label=f'α_pred={alpha_pred:.3f}')
    ax.axvline(1.0, color='k', ls=':', lw=2, label='Newton α=1')
    ax.set_xlabel('α_libre'); ax.set_ylabel('densité')
    ax.set_title('Distribution α (tous centres)'); ax.legend(fontsize=7); ax.grid(alpha=0.3)

    # (1) α vs R²_free — scatter avec correction robuste
    ax = axes[1]
    pairs = [
        (r["R2_free"], r["alpha"])
        for r in records
        if np.isfinite(r["R2_free"]) and np.isfinite(r["alpha"])
    ]
    if pairs:
        R2s = np.array([p[0] for p in pairs])
        alphs = np.array([p[1] for p in pairs])
        ax.scatter(R2s, alphs, s=20, alpha=0.5, color="steelblue")
    ax.axhline(alpha_pred, color='r', ls='--', lw=2, label=f'α_pred={alpha_pred:.3f}')
    ax.axhline(1.0, color='k', ls=':', lw=2, label='Newton')
    for thresh in args.r2_thresholds:
        ax.axvline(thresh, color='gray', ls=':', lw=1)
    ax.set_xlabel('R²_libre'); ax.set_ylabel('α')
    ax.set_title('α vs R²  (sélection par seuil)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (2) Courbe d_s(τ)
    ax = axes[2]
    ax.semilogx(taus, ds_curve, 'o-', ms=3, lw=1.5, color='mediumseagreen')
    ax.axhline(d_s_heat, color='r', ls='--', lw=2,
               label=f'd_s effectif={d_s_heat:.3f}')
    ax.axhline(3.0, color='k', ls=':', lw=1.5, label='d_s=3 (Newton 3D)')
    ax.axhline(2.0, color='gray', ls='--', lw=1)
    ax.set_xlabel('τ'); ax.set_ylabel('d_s(τ)')
    ax.set_title('Dimension spectrale (noyau chaleur)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "figures.pdf"),
                dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\n[OK] Résultats → {args.output_dir}/")
    print(f"     centers_alpha.csv | summary.csv | summary.json | figures.pdf")
    print(f"     Temps total : {time.time()-t0:.1f}s")
    print("DONE")

if __name__ == '__main__':
    main()

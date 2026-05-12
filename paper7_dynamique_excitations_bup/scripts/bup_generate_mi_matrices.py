#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bup_generate_mi_matrices.py
══════════════════════════════════════════════════════════════════
Génère les matrices MI depuis les simulations BuP et les sauvegarde
en CSV pour bup_kill_epsilon_micro.py

Usage :
  python3 bup_generate_mi_matrices.py
  python3 bup_generate_mi_matrices.py --N 9 12 16 --lambda-n 8

Sortie :
  results_mi_matrices/MI_N{N}_lam{lam:.2f}.csv
  pour chaque (N, λ)
══════════════════════════════════════════════════════════════════
"""
import argparse, os
import numpy as np

# ── Constantes BuP ──────────────────────────────────────────────
DC = 3.059842935509462


# ══════════════════════════════════════════════════════════════
# Génération état BuP
# ══════════════════════════════════════════════════════════════

def factor_grid(N):
    best = None
    for r in range(1, int(N**0.5)+2):
        if N % r == 0:
            c = N // r; s = abs(c-r)
            if best is None or s < best[0]: best = (s, r, c)
    return best[1], best[2]


def generate_state(N, lam):
    rows, cols = factor_grid(N)
    J, h, J_nl = 1.0, 0.5+lam*0.5, 0.8
    pairs = set()
    for a, b in [(0, N-1), (cols-1, N-cols)]:
        if 0 <= a < N and 0 <= b < N and a != b:
            pairs.add(tuple(sorted((a,b))))
    pairs = list(pairs)
    dim = 2**N
    psi = np.zeros(dim, dtype=complex)
    for cfg in range(dim):
        bits = [(cfg>>i)&1 for i in range(N)]
        neighbors = []
        for i in range(N):
            r, c = i//cols, i%cols
            if r > 0:        neighbors.append((i, i-cols))
            if r < rows-1:   neighbors.append((i, i+cols))
            if c > 0:        neighbors.append((i, i-1))
            if c < cols-1:   neighbors.append((i, i+1))
        loc  = sum(1 if bits[i]==bits[j] else -1 for i,j in neighbors if i<j)
        nloc = sum(lam*J_nl*(1 if bits[i]==bits[j] else -1)
                   for i,j in pairs) if lam>0 else 0.0
        psi[cfg] = (np.exp(J*loc/max(2.0,N/4.0))
                    * np.exp(nloc/max(1.0,len(pairs)))
                    * np.exp(1j*h*sum(bits)))
    norm = np.linalg.norm(psi)
    return psi/norm if norm > 1e-15 else psi


# ══════════════════════════════════════════════════════════════
# Matrice MI exacte
# ══════════════════════════════════════════════════════════════

def mi_matrix(psi, N):
    def rho1(i):
        t = psi.reshape([2]*N)
        other = [k for k in range(N) if k != i]
        m = np.transpose(t, [i]+other).reshape(2, 2**(N-1))
        return m @ m.conj().T

    def rho2(i, j):
        t = psi.reshape([2]*N)
        other = [k for k in range(N) if k != i and k != j]
        m = np.transpose(t, [i,j]+other).reshape(4, 2**(N-2))
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
            MI[i,j] = MI[j,i] = max(Ss[i] + Ss[j] - Sij, 0.0)
    return MI


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N',         type=int, nargs='+', default=[9, 12, 16])
    ap.add_argument('--lambda-n',  type=int,            default=8)
    ap.add_argument('--output-dir',                     default='results_mi_matrices')
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    lam_grid = np.linspace(0.0, 1.0, args.lambda_n)

    print("="*60)
    print("BuP — Génération des matrices MI")
    print(f"N ∈ {args.N}  λ ∈ [{lam_grid[0]:.2f}, {lam_grid[-1]:.2f}]"
          f" ({args.lambda_n} pts)")
    print("="*60)

    total = len(args.N) * len(lam_grid)
    done  = 0

    for N in args.N:
        rows, cols = factor_grid(N)
        print(f"\nN={N:2d}  grille {rows}×{cols}")
        for lam in lam_grid:
            try:
                psi = generate_state(N, lam)
                MI  = mi_matrix(psi, N)

                # Sauvegarde CSV
                fname = f"MI_N{N:02d}_lam{lam:.2f}.csv"
                path  = os.path.join(args.output_dir, fname)
                np.savetxt(path, MI, delimiter=",", fmt="%.8f")

                # Statistiques rapides
                W = MI.copy(); np.fill_diagonal(W, 0)
                strength = W.sum(1)
                print(f"  λ={lam:.2f}  s_mean={strength.mean():.4f}"
                      f"  s_std={strength.std():.4f}"
                      f"  MI_max={MI.max():.4f}  → {fname}")
                done += 1

            except Exception as e:
                print(f"  λ={lam:.2f}  ERREUR: {e}")

    print(f"\n{done}/{total} matrices générées → {args.output_dir}/")
    print("\nCommande pour tester epsilon :")
    print(f"  python3 bup_kill_epsilon_micro.py \\")
    print(f"    --mi-files \"{args.output_dir}/MI_*.csv\" \\")
    print(f"    --k 5 --output-dir results_epsilon_micro")
    print("\nDONE")


if __name__ == '__main__':
    main()

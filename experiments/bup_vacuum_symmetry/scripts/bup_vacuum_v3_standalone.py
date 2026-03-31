
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bup_vacuum_v3_standalone.py

Standalone version of the BuP vacuum effective symmetry test.
No imports from local helper modules are required.
"""

import argparse
import csv
import json
import math
import warnings
from dataclasses import dataclass

import numpy as np

warnings.filterwarnings("ignore")

PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
PAULIS = [PAULI_X, PAULI_Y, PAULI_Z]


def kron_op(op, i, n):
    ops = [I2] * n
    ops[i] = op
    out = ops[0]
    for k in range(1, n):
        out = np.kron(out, ops[k])
    return out


@dataclass
class VacuumConfig:
    n_qubits: int = 8
    n_source: int = 2
    J_zz: float = 2.0
    J_xy: float = 0.5
    h_scale: float = 0.1
    lie_tol: float = 1e-3
    n_seeds: int = 20
    output_csv: str = "bup_vacuum_v3_results.csv"
    output_json: str = "bup_vacuum_v3_summary.json"
    knn_k: int = 3


def build_H_BUP(n, J_zz, J_xy, h_scale, seed):
    rng = np.random.default_rng(seed)
    dim = 2 ** n
    H = np.zeros((dim, dim), dtype=complex)

    for i in range(n):
        j = (i + 1) % n

        H -= J_zz * kron_op(PAULI_Z, i, n) @ kron_op(PAULI_Z, j, n)
        H -= J_xy * (
            kron_op(PAULI_X, i, n) @ kron_op(PAULI_X, j, n)
            + kron_op(PAULI_Y, i, n) @ kron_op(PAULI_Y, j, n)
        )

        hx, hy, hz = rng.standard_normal(3) * h_scale
        H += hx * kron_op(PAULI_X, i, n)
        H += hy * kron_op(PAULI_Y, i, n)
        H += hz * kron_op(PAULI_Z, i, n)

    return H


def ground_state(H):
    eigvals, eigvecs = np.linalg.eigh(H)
    return eigvecs[:, 0], float(np.real(eigvals[0]))


def spectral_gap(H):
    eigvals = np.linalg.eigvalsh(H)
    return float(np.real(eigvals[1] - eigvals[0]))


def reduced_density_matrix_pure(state, keep, n):
    keep = sorted(keep)
    rest = [i for i in range(n) if i not in keep]
    psi = state.reshape((2,) * n)
    psi = np.transpose(psi, axes=keep + rest)
    dim_keep = 2 ** len(keep)
    dim_rest = 2 ** len(rest)
    psi = psi.reshape(dim_keep, dim_rest)
    rho = psi @ psi.conj().T
    rho = 0.5 * (rho + rho.conj().T)
    return rho


def pauli_tensor(rho_ij):
    T = np.zeros((3, 3), dtype=np.float64)
    for a, A in enumerate(PAULIS):
        for b, B in enumerate(PAULIS):
            T[a, b] = float(np.real_if_close(np.trace(rho_ij @ np.kron(A, B))))
    return T


def polar_so3(T, eps=1e-12):
    if np.linalg.norm(T) < eps:
        return np.eye(3)
    U, _, Vt = np.linalg.svd(T, full_matrices=False)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1.0
        R = U @ Vt
    return R


def rotation_to_lie_vector(R, tol=1e-12):
    tr = (np.trace(R) - 1.0) / 2.0
    tr = min(max(float(np.real_if_close(tr)), -1.0), 1.0)
    theta = math.acos(tr)

    if theta < tol:
        return np.zeros(3, dtype=float)

    denom = 2.0 * math.sin(theta)
    if abs(denom) < tol:
        return np.zeros(3, dtype=float)

    axis = np.array(
        [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]],
        dtype=float,
    ) / denom

    norm = np.linalg.norm(axis)
    if norm < tol:
        return np.zeros(3, dtype=float)

    axis = axis / norm
    return theta * axis


def build_rotation_graph(psi, cfg):
    n = cfg.n_qubits
    weights = {}
    tensors = {}

    for i in range(n):
        for j in range(i + 1, n):
            rho_ij = reduced_density_matrix_pure(psi, [i, j], n)
            T = pauli_tensor(rho_ij)
            w = float(np.linalg.norm(T, ord="fro"))
            weights[(i, j)] = w
            tensors[(i, j)] = T

    neigh = {i: [] for i in range(n)}
    for (i, j), w in weights.items():
        neigh[i].append((w, j))
        neigh[j].append((w, i))

    edges = set()
    for i in range(n):
        neigh[i].sort(reverse=True, key=lambda x: x[0])
        for _, j in neigh[i][:cfg.knn_k]:
            edges.add(tuple(sorted((i, j))))

    R = {}
    S = {}
    for i, j in edges:
        T = tensors[(i, j)]
        R[(i, j)] = polar_so3(T)
        S[(i, j)] = weights[(i, j)]
    return R, S


def lie_vectors(R, lie_tol):
    vecs = []
    for Rij in R.values():
        v = rotation_to_lie_vector(Rij, tol=lie_tol)
        if np.linalg.norm(v) > lie_tol:
            vecs.append(v)
    if not vecs:
        return np.zeros((0, 3), dtype=float)
    return np.array(vecs, dtype=float)


def lie_effective_dim(R, lie_tol):
    vecs = lie_vectors(R, lie_tol)
    if len(vecs) == 0:
        return 0
    _, s, _ = np.linalg.svd(vecs, full_matrices=False)
    if len(s) == 0:
        return 0
    smax = s[0]
    return int(np.sum(s > max(lie_tol, 0.1 * smax)))


def mean_commutator(R, lie_tol):
    vecs = lie_vectors(R, lie_tol)
    m = len(vecs)
    if m < 2:
        return 0.0
    vals = []
    for i in range(m):
        for j in range(i + 1, m):
            vals.append(float(np.linalg.norm(np.cross(vecs[i], vecs[j]))))
    return float(np.mean(vals)) if vals else 0.0


def _empty_result():
    return {
        "u1_like": False,
        "order_parameter": 0.5,
        "sigma_angles": 1.0,
        "ratio_singular": 1.0,
        "u1_axis_fraction": 0.0,
        "dim_lie": 0,
        "mean_commutator": 0.0,
        "n_lie_vecs": 0,
        "best_axis": None,
    }


def analyze_effective_symmetry(psi, cfg):
    R, _ = build_rotation_graph(psi, cfg)
    vecs = lie_vectors(R, cfg.lie_tol)

    if len(vecs) < 3:
        return _empty_result()

    norms = np.linalg.norm(vecs, axis=1, keepdims=True)
    vn = vecs / (norms + 1e-12)

    _, sv, Vt = np.linalg.svd(vn, full_matrices=False)
    n_hat = Vt[0]
    ratio_s = float(sv[0] / (sv[1] + 1e-10))

    dots = np.abs(vn @ n_hat)
    dots = np.clip(dots, 0.0, 1.0)
    angles = np.arccos(dots)

    op = float(np.mean(dots))
    sigma = float(np.std(angles))
    u1_axis_fraction = float(np.mean(angles < np.pi / 6))

    u1_like = (ratio_s > 2.5) and (sigma < 0.4)

    dim_lie = lie_effective_dim(R, cfg.lie_tol)
    comm = mean_commutator(R, cfg.lie_tol)

    return {
        "u1_like": u1_like,
        "order_parameter": op,
        "sigma_angles": sigma,
        "ratio_singular": ratio_s,
        "u1_axis_fraction": u1_axis_fraction,
        "dim_lie": dim_lie,
        "mean_commutator": comm,
        "n_lie_vecs": int(len(vecs)),
        "best_axis": n_hat.tolist(),
    }


def experiment_seeds(cfg):
    results = []
    print(f"\n  ── Robustesse sur {cfg.n_seeds} seeds (h={cfg.h_scale}) ──")
    for seed in range(cfg.n_seeds):
        H = build_H_BUP(cfg.n_qubits, cfg.J_zz, cfg.J_xy, cfg.h_scale, seed)
        psi, E0 = ground_state(H)
        gap = spectral_gap(H)
        an = analyze_effective_symmetry(psi, cfg)

        row = {
            "seed": seed,
            "E0": E0,
            "gap": gap,
            "h_scale": cfg.h_scale,
            "J_zz": cfg.J_zz,
            "J_xy": cfg.J_xy,
            **an,
        }
        results.append(row)

        if seed % 5 == 0:
            verdict = "U(1)-like ✓" if an["u1_like"] else "SO(3)-like"
            print(
                f"  seed={seed:3d}: ratio_s={an['ratio_singular']:6.2f}  "
                f"OP={an['order_parameter']:.3f}  "
                f"σ={an['sigma_angles']:.3f}  → {verdict}"
            )
    return results


def experiment_h_scan(cfg):
    h_values = [0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 1.00, 1.50, 2.00]
    results = []
    print(f"\n  ── Scan sur h_scale (sélection effective d'axe) ──")
    print(f"  {'h':>6}  {'OP':>7}  {'ratio_s':>8}  {'σ':>7}  {'U(1)%':>7}  verdict")
    print(f"  {'─'*55}")

    for h in h_values:
        seed_results = []
        for seed in range(10):
            H = build_H_BUP(cfg.n_qubits, cfg.J_zz, cfg.J_xy, h, seed)
            psi, _ = ground_state(H)
            an = analyze_effective_symmetry(psi, cfg)
            seed_results.append(an)

        op_mean = np.mean([r["order_parameter"] for r in seed_results])
        ratio_mean = np.mean([r["ratio_singular"] for r in seed_results])
        sigma_mean = np.mean([r["sigma_angles"] for r in seed_results])
        u1_frac = np.mean([r["u1_like"] for r in seed_results])

        row = {
            "h_scale": h,
            "op_mean": float(op_mean),
            "ratio_s_mean": float(ratio_mean),
            "sigma_mean": float(sigma_mean),
            "u1_fraction": float(u1_frac),
            "verdict": "U(1)-like" if u1_frac > 0.5 else "SO(3)-like",
        }
        results.append(row)

        print(
            f"  {h:6.3f}  {op_mean:7.3f}  {ratio_mean:8.2f}  "
            f"{sigma_mean:7.3f}  {u1_frac*100:6.1f}%  {row['verdict']}"
        )
    return results


def experiment_J_scan(cfg):
    J_ratios = [0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]
    results = []
    print(f"\n  ── Scan sur J_zz/J_xy (anisotropie) ──")

    for ratio in J_ratios:
        J_zz = ratio
        J_xy = 0.5
        seed_results = []

        for seed in range(10):
            H = build_H_BUP(cfg.n_qubits, J_zz, J_xy, cfg.h_scale, seed)
            psi, _ = ground_state(H)
            an = analyze_effective_symmetry(psi, cfg)
            seed_results.append(an)

        op_mean = np.mean([r["order_parameter"] for r in seed_results])
        ratio_mean = np.mean([r["ratio_singular"] for r in seed_results])
        u1_frac = np.mean([r["u1_like"] for r in seed_results])

        row = {
            "J_ratio": ratio,
            "J_zz": J_zz,
            "J_xy": J_xy,
            "op_mean": float(op_mean),
            "ratio_s_mean": float(ratio_mean),
            "u1_fraction": float(u1_frac),
            "verdict": "U(1)-like" if u1_frac > 0.5 else "SO(3)-like",
        }
        results.append(row)

        print(
            f"  J_zz/J_xy={ratio:4.1f}: OP={op_mean:.3f}  "
            f"ratio_s={ratio_mean:.2f}  U(1)={u1_frac*100:.0f}%  [{row['verdict']}]"
        )
    return results


def experiment_random_control(cfg, n_random=50):
    rng = np.random.default_rng(999)
    dim = 2 ** cfg.n_qubits
    results = []

    for _ in range(n_random):
        psi = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        psi /= np.linalg.norm(psi)
        an = analyze_effective_symmetry(psi, cfg)
        results.append(an)

    return {
        "op_mean": float(np.mean([r["order_parameter"] for r in results])),
        "op_std": float(np.std([r["order_parameter"] for r in results])),
        "ratio_s_mean": float(np.mean([r["ratio_singular"] for r in results])),
        "sigma_mean": float(np.mean([r["sigma_angles"] for r in results])),
        "u1_fraction": float(np.mean([r["u1_like"] for r in results])),
    }


def run_all(cfg):
    print(f"\n{'='*65}")
    print(f"  BUP Vacuum Effective Symmetry — v3.0")
    print(f"  N={cfg.n_qubits} qubits | h={cfg.h_scale} | J_zz={cfg.J_zz} | J_xy={cfg.J_xy}")
    print(f"  H_BUP = XXZ + champs aléatoires locaux")
    print(f"{'='*65}")

    seeds_res = experiment_seeds(cfg)
    h_scan = experiment_h_scan(cfg)
    J_scan = experiment_J_scan(cfg)
    control = experiment_random_control(cfg)

    print(f"\n  ── Contrôle aléatoire (50 états) ──")
    print(f"  OP moyen      : {control['op_mean']:.3f} ± {control['op_std']:.3f}")
    print(f"  ratio_s moyen : {control['ratio_s_mean']:.2f}")
    print(f"  U(1)-like %   : {control['u1_fraction']*100:.1f}%")

    return {
        "seeds": seeds_res,
        "h_scan": h_scan,
        "J_scan": J_scan,
        "control": control,
    }


def print_final_summary(results, cfg):
    seeds = results["seeds"]
    u1_rate = np.mean([r["u1_like"] for r in seeds])
    op_mean = np.mean([r["order_parameter"] for r in seeds])
    ctrl = results["control"]

    print(f"\n{'='*65}")
    print(f"  CONCLUSION")
    print(f"{'='*65}")
    print(f"  État fondamental H_BUP (h={cfg.h_scale}) :")
    print(f"    U(1)-like émergent : {u1_rate*100:.0f}% des seeds")
    print(f"    OP moyen           : {op_mean:.3f}  (vs {ctrl['op_mean']:.3f} aléatoire)")
    print(f"    Δ(OP)              : {op_mean - ctrl['op_mean']:+.3f}")

    if u1_rate > 0.8:
        print(f"\n  → le vide effectif sélectionne robustement un axe U(1)-like ✓")
    elif u1_rate > 0.4:
        print(f"\n  → sélection d’axe partielle / dépendante du désordre")
    else:
        print(f"\n  → pas de sélection d’axe robuste dans ce régime")

    h_scan = results["h_scan"]
    onset = next((r["h_scale"] for r in h_scan if r["u1_fraction"] > 0.5), None)
    if onset is not None:
        print(f"  → sélection U(1)-like déjà visible pour h ≲ {onset:.3f} dans ce scan")
    print(f"{'='*65}\n")


def save_csv(results, path):
    rows = []
    for r in results.get("seeds", []):
        r2 = {k: v for k, v in r.items() if k != "best_axis"}
        rows.append({"category": "seed", **r2})
    for r in results.get("h_scan", []):
        rows.append({"category": "h_scan", **r})
    for r in results.get("J_scan", []):
        rows.append({"category": "J_scan", **r})

    if not rows:
        return

    keys = sorted(set().union(*[row.keys() for row in rows]))
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  ✓ CSV → {path}")


def save_json(results, path):
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  ✓ JSON → {path}")


def parse_args():
    p = argparse.ArgumentParser(description="BUP Vacuum Effective Symmetry v3 standalone")
    p.add_argument("--n_qubits", type=int, default=8)
    p.add_argument("--n_source", type=int, default=2)
    p.add_argument("--h_scale", type=float, default=0.1)
    p.add_argument("--J_zz", type=float, default=2.0)
    p.add_argument("--J_xy", type=float, default=0.5)
    p.add_argument("--n_seeds", type=int, default=20)
    p.add_argument("--lie_tol", type=float, default=1e-3)
    p.add_argument("--knn_k", type=int, default=3)
    p.add_argument("--output_csv", type=str, default="bup_vacuum_v3_results.csv")
    p.add_argument("--output_json", type=str, default="bup_vacuum_v3_summary.json")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg = VacuumConfig(
        n_qubits=args.n_qubits,
        n_source=args.n_source,
        h_scale=args.h_scale,
        J_zz=args.J_zz,
        J_xy=args.J_xy,
        n_seeds=args.n_seeds,
        lie_tol=args.lie_tol,
        output_csv=args.output_csv,
        output_json=args.output_json,
        knn_k=args.knn_k,
    )
    results = run_all(cfg)
    print_final_summary(results, cfg)
    save_csv(results, cfg.output_csv)
    save_json(results, cfg.output_json)

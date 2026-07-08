#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper 25 — Variational Foundation of Bottom-Up Quantum Gravity
Builds canonical summary tables for the foundation paper.

Outputs:
  - paper25_chain_table.csv
  - paper25_canonical_dictionary.csv
  - paper25_evidence_status_table.csv
  - paper25_open_debts_table.csv
  - paper25_summary.json
"""

from pathlib import Path
import csv
import json
from datetime import datetime

OUTDIR = Path("paper25_variational_foundation/results/paper25_foundation_summary_v1")
OUTDIR.mkdir(parents=True, exist_ok=True)


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


chain_rows = [
    {
        "block": "Cosmological dimensional sector",
        "papers": "Papers 2--4",
        "core_chain": "d(z) -> G_eff(z) -> sigma8 / JWST local dimension",
        "main_equation": "G_eff(z)=2G/(d(z)-1)",
        "status": "phenomenological cosmological branch",
        "role_in_paper25": "Feeds H_dim and cross-scale dimensional flow",
    },
    {
        "block": "Emergent matter and weak-field source",
        "papers": "Papers 8--9",
        "core_chain": "delta W_loc -> S_flux -> Phi_BuP -> |delta R|",
        "main_equation": "L_ent Phi_BuP = S_flux",
        "status": "finite graph weak-field precursor",
        "role_in_paper25": "Defines entanglement Green function and source propagation",
    },
    {
        "block": "Newtonian propagator limit",
        "papers": "Papers 10--11",
        "core_chain": "L_ent^+ -> Green function -> alpha_eff",
        "main_equation": "alpha_eff = 2 d_s/d_w + d_w - 4",
        "status": "effective propagator law",
        "role_in_paper25": "Connects graph dimensions to weak-field radial behavior",
    },
    {
        "block": "Galactic sector",
        "papers": "Papers 12--14",
        "core_chain": "Sigma(R) -> d_s(r), d_w(r) -> alpha_eff(r) -> V(r)",
        "main_equation": "V_BuP from alpha_eff(r)",
        "status": "phenomenological validation on SPARC",
        "role_in_paper25": "Macroscopic expression of dimensional propagator law",
    },
    {
        "block": "Spectral continuum limit",
        "papers": "Papers 15--16",
        "core_chain": "W -> L_N -> -Delta_g",
        "main_equation": "c_N L_N -> -Delta_g",
        "status": "controlled numerical convergence",
        "role_in_paper25": "Geometric foundation of continuum limit",
    },
    {
        "block": "Discrete Ricci limit",
        "papers": "Paper 17",
        "core_chain": "kappa_OR -> R_mu nu u^mu u^nu",
        "main_equation": "kappa_OR/epsilon ~= B_N + C_N R_mu nu u^mu u^nu",
        "status": "mean-level affine calibration",
        "role_in_paper25": "Curvature side of Einstein limit",
    },
    {
        "block": "Modular source sector",
        "papers": "Paper 18",
        "core_chain": "delta W_loc -> delta S_A ~= delta<K_A> -> delta kappa",
        "main_equation": "delta S_A ~= delta<K_A>",
        "status": "controlled modular-source precursor",
        "role_in_paper25": "Source side of Einstein limit",
    },
    {
        "block": "Effective Einstein assembly",
        "papers": "Paper 19",
        "core_chain": "spectral + Ricci + source arrows -> Einstein effective equation",
        "main_equation": "G_mu nu + Lambda g_mu nu = 8 pi G_eff T_mu nu^ent + H_mu nu",
        "status": "controlled assembly / theorem target",
        "role_in_paper25": "Main continuum target",
    },
    {
        "block": "Correction tensor",
        "papers": "Paper 20",
        "core_chain": "Einstein limit deviations -> H_mu nu sectors",
        "main_equation": "H_static = H_spec + H_curv + H_source + H_dim + H_nonlocal + H_topo + H_finite",
        "status": "classification plus proxies",
        "role_in_paper25": "Organizes deviations from smooth Einstein regime",
    },
    {
        "block": "Strong-lensing fixed point",
        "papers": "Paper 21",
        "core_chain": "SLACS -> alpha_eff ~= 1",
        "main_equation": "<alpha_eff> ~= 1.014 near log M_star ~= 11.58--11.60",
        "status": "phenomenological fixed point",
        "role_in_paper25": "Nonlocal / optical correction sector",
    },
    {
        "block": "Gravitational-wave sector",
        "papers": "Paper 22",
        "core_chain": "W -> edge Hessian -> modes q_n(t)",
        "main_equation": "H_full = H_static + H_dyn",
        "status": "dynamic numerical sector",
        "role_in_paper25": "Adds H_dyn and wave-speed constraints",
    },
    {
        "block": "Modular tomography",
        "papers": "Paper 23",
        "core_chain": "I(i:j) tomography -> rho_ent^cut(j) -> 1/beta_j",
        "main_equation": "rho_ent^cut(j) = sum_{k notin A} I(j:k)",
        "status": "experimental testability proposal",
        "role_in_paper25": "Microscopic observability of W_ij",
    },
]

dictionary_rows = [
    {
        "bup_object": "W_ij = I(i:j)",
        "meaning": "mutual-information graph",
        "continuum_target": "entanglement geometry seed",
        "status": "fundamental postulate",
    },
    {
        "bup_object": "L_ent = D - W",
        "meaning": "entanglement Laplacian",
        "continuum_target": "-Delta_g",
        "status": "supported by Papers 15--16",
    },
    {
        "bup_object": "L_ent^+",
        "meaning": "Moore--Penrose pseudoinverse",
        "continuum_target": "(-Delta_g)^(-1)",
        "status": "weak-field Green function candidate",
    },
    {
        "bup_object": "kappa_OR",
        "meaning": "Ollivier--Ricci curvature on W",
        "continuum_target": "R_mu nu u^mu u^nu",
        "status": "Paper 17 affine calibration",
    },
    {
        "bup_object": "delta<K_A>",
        "meaning": "modular energy variation",
        "continuum_target": "source precursor for T_mu nu^ent",
        "status": "Paper 18 modular first law",
    },
    {
        "bup_object": "alpha_eff",
        "meaning": "effective propagator exponent",
        "continuum_target": "weak-field radial behavior",
        "status": "Papers 11, 12, 21",
    },
    {
        "bup_object": "H_mu nu^BuP",
        "meaning": "deviation from effective Einstein limit",
        "continuum_target": "higher-order / finite / nonlocal corrections",
        "status": "Paper 20 classification",
    },
]

evidence_rows = [
    {
        "claim": "Entanglement Laplacian approximates Laplace-Beltrami",
        "main_result": "circle error 0.0051; flat torus error 0.0266",
        "source_papers": "Papers 15--16",
        "status": "controlled numerical convergence",
        "limitation": "normalization c_N still analytical debt",
    },
    {
        "claim": "Ollivier--Ricci curvature carries Ricci signal",
        "main_result": "B_N=-0.301281, C_N=0.391933 at N=512",
        "source_papers": "Paper 17",
        "status": "mean-level affine calibration",
        "limitation": "not yet pointwise tensor theorem",
    },
    {
        "claim": "Graph modular first law holds",
        "main_result": "slope ~=0.989; R2 ~=0.9966",
        "source_papers": "Paper 18",
        "status": "controlled numerical evidence",
        "limitation": "rho_A is graph proxy; T_mu nu not yet fully reconstructed",
    },
    {
        "claim": "Modular source predicts curvature response",
        "main_result": "R2 0.967--0.982 for |delta K| vs near curvature response",
        "source_papers": "Paper 18",
        "status": "controlled numerical evidence",
        "limitation": "scalar source precursor, not full tensor",
    },
    {
        "claim": "Discrete Poisson propagation preserves source-curvature signal",
        "main_result": "rho(S_flux,|delta R|)=0.741; rho(Phi,|delta R|)=0.738",
        "source_papers": "Paper 9",
        "status": "finite graph weak-field precursor",
        "limitation": "small N and candidate source",
    },
    {
        "claim": "SPARC branch performs well phenomenologically",
        "main_result": "median chi2_red BuP=0.467 vs NFW2p=1.332; 83.4 percent better",
        "source_papers": "Paper 14",
        "status": "phenomenological validation",
        "limitation": "not a derivation of Einstein sector",
    },
    {
        "claim": "Correction tensor has structured sectors",
        "main_result": "7 static sectors plus dynamic sector from Paper 22",
        "source_papers": "Papers 20--22",
        "status": "classification",
        "limitation": "dim/nonlocal are currently proxy diagnostics",
    },
]

debt_rows = [
    {
        "debt": "Derive c_N",
        "why_it_matters": "Needed for rigorous L_N -> -Delta_g",
        "target_paper": "Paper 26",
        "priority": "critical",
    },
    {
        "debt": "Reconstruct full Ricci tensor",
        "why_it_matters": "Need R_mu nu, not only R_mu nu u^mu u^nu",
        "target_paper": "Paper 27",
        "priority": "high",
    },
    {
        "debt": "Construct T_mu nu^ent",
        "why_it_matters": "Paper 18 gives scalar modular precursor only",
        "target_paper": "Paper 28",
        "priority": "high",
    },
    {
        "debt": "Define common norm for H_mu nu",
        "why_it_matters": "Needed to compare correction sectors consistently",
        "target_paper": "Paper 29",
        "priority": "high",
    },
    {
        "debt": "Unify cross-scale alpha_eff constraints",
        "why_it_matters": "Connect Mercury, SPARC, SLACS, GW and cosmology",
        "target_paper": "Paper 30",
        "priority": "medium-high",
    },
    {
        "debt": "General micro-to-macro theorem",
        "why_it_matters": "Need theorem for |Psi> -> I_ij -> matter profiles",
        "target_paper": "future",
        "priority": "medium",
    },
]

write_csv(
    OUTDIR / "paper25_chain_table.csv",
    chain_rows,
    ["block", "papers", "core_chain", "main_equation", "status", "role_in_paper25"],
)

write_csv(
    OUTDIR / "paper25_canonical_dictionary.csv",
    dictionary_rows,
    ["bup_object", "meaning", "continuum_target", "status"],
)

write_csv(
    OUTDIR / "paper25_evidence_status_table.csv",
    evidence_rows,
    ["claim", "main_result", "source_papers", "status", "limitation"],
)

write_csv(
    OUTDIR / "paper25_open_debts_table.csv",
    debt_rows,
    ["debt", "why_it_matters", "target_paper", "priority"],
)

summary = {
    "paper": "Paper 25 — Variational Foundation of Bottom-Up Quantum Gravity",
    "generated_at": datetime.utcnow().isoformat() + "Z",
    "central_variable": "W_ij = I(i:j)",
    "fundamental_equation": "delta S_BuP[W,rho] / delta W_ij = 0",
    "continuum_target": "G_mu nu + Lambda_ent g_mu nu = 8 pi G_eff T_mu nu^ent + H_mu nu^BuP",
    "num_chain_blocks": len(chain_rows),
    "num_dictionary_entries": len(dictionary_rows),
    "num_evidence_entries": len(evidence_rows),
    "num_open_debts": len(debt_rows),
    "next_priority": "Paper 26 — derive normalization c_N for L_N -> -Delta_g",
}

with open(OUTDIR / "paper25_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"[OK] Wrote Paper 25 summary tables to {OUTDIR}")

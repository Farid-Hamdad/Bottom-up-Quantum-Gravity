# Paper 18 — Emergent Stress-Energy from Entanglement

**From local entanglement perturbations to an effective stress-energy source**

---

## Status

Paper 18 studies the third arrow required for the effective Einstein limit of BuP:

\[
\delta W_{\rm loc}
\longrightarrow
T_{\mu\nu}^{\rm ent}.
\]

Paper 16 studied the spectral geometry side:

\[
L_N\to\Delta_g.
\]

Paper 17 studied the curvature side:

\[
\kappa_{ij}^{OR}
\to
R_{\mu\nu}u^\mu u^\nu.
\]

Paper 18 studies the source side: how a local perturbation of the entanglement graph becomes an effective stress-energy source.

---

## 1. Starting point

The fundamental object is the entanglement graph:

\[
W_{ij}=I(i:j).
\]

A local source is modeled as a perturbation:

\[
W_{ij}
\to
W'_{ij}
=
W_{ij}
+
\delta W_{ij}^{\rm loc}.
\]

The central question is:

\[
\boxed{
\text{Can }\delta W_{\rm loc}\text{ be interpreted as an effective source }T_{\mu\nu}^{\rm ent}?
}
\]

---

## 2. Previous evidence

### Paper 8

Paper 8 reconstructed an effective source proxy from local entanglement perturbations:

\[
T_{\mu\nu}^{\rm eff}
=
T_{\mu\nu}^{\rm matter}[\delta W_{\rm loc}]
+
T_{\mu\nu}^{\rm ent}[d_s].
\]

It found:

\[
\rho_{\rm Spearman}=0.741,
\qquad
p=1.84\times10^{-4}.
\]

This supported:

\[
\delta W_{\rm loc}
\to
T_{\mu\nu}^{\rm eff}.
\]

### Paper 15 Step E

Paper 15 showed that a radial entanglement defect produces a localized curvature response:

\[
\delta W_{\rm loc}
\to
\delta\kappa(r).
\]

On the sphere:

\[
{\rm Spearman}(\phi_{\rm edge},|\Delta\kappa|)
=
0.752,
\]

with near/far response:

\[
11.51.
\]

Thus Paper 18 starts from two positive prior signals.

---

## 3. Core identity: modular first law

The key theoretical input is the modular first law:

\[
\delta S_A
=
\delta\langle K_A\rangle.
\]

with

\[
K_A=-\log\rho_A.
\]

In continuum QFT, local modular Hamiltonians relate entropy variations to stress-energy variations:

\[
\delta\langle K_A\rangle
\sim
\int_A \xi^\mu \delta T_{\mu\nu}d\Sigma^\nu.
\]

In BuP, the graph analogue is:

\[
\delta S_A[W]
\simeq
\delta\langle K_A[W]\rangle
\to
T_{\mu\nu}^{\rm ent}.
\]

---

## 4. Target chain

Paper 18 aims to establish:

\[
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r)
\to
T_{\mu\nu}^{\rm ent}.
\]

---

## 5. Main numerical results

Paper 18 currently has two positive results.

| Step | Geometry | Quantity | Value | Target | Status |
|---|---|---|---:|---|---|
| v1 graph modular first law | flat_torus2d | slope \(\delta S\) vs \(\delta\langle K\rangle\) | 0.989481 | 1 | positive |
| v1 graph modular first law | flat_torus2d | \(R^2\) \(\delta S\) vs \(\delta\langle K\rangle\) | 0.996590 | near 1 | positive |
| v1 graph modular first law | sphere | slope \(\delta S\) vs \(\delta\langle K\rangle\) | 0.989718 | 1 | positive |
| v1 graph modular first law | sphere | \(R^2\) \(\delta S\) vs \(\delta\langle K\rangle\) | 0.996592 | near 1 | positive |
| v2 modular source curvature | flat_torus2d | near/far localization ratio | 6.293007 | \(>1\) | positive |
| v2 modular source curvature | flat_torus2d | \(R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})\) | 0.967229 | near 1 | positive |
| v2 modular source curvature | flat_torus2d | signed Pearson \(\delta K\to\Delta\kappa_{\rm near}\) | -0.999168 | \(|r|\) near 1 | positive |
| v2 modular source curvature | sphere | near/far localization ratio | 9.533529 | \(>1\) | positive |
| v2 modular source curvature | sphere | \(R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})\) | 0.982252 | near 1 | positive |
| v2 modular source curvature | sphere | signed Pearson \(\delta K\to\Delta\kappa_{\rm near}\) | -0.999428 | \(|r|\) near 1 | positive |

---

## 6. Result v1 — Graph modular first law

A graph density proxy is defined on a region \(A\) by

\[
\rho_A
=
\frac{(L_A+\mu I)^{-1}}
{\mathrm{Tr}(L_A+\mu I)^{-1}}.
\]

Then

\[
S_A=-\mathrm{Tr}(\rho_A\log\rho_A),
\qquad
K_A=-\log\rho_A.
\]

After a local radial perturbation of \(W_{ij}\), the test measures

\[
\delta S_A=S_A(W')-S_A(W),
\]

and

\[
\delta\langle K_A\rangle
=
\mathrm{Tr}\left[(\rho'_A-\rho_A)K_A\right].
\]

The result is strongly positive.

On the flat torus:

\[
\delta S_A
=
-0.000068
+
0.989481\,\delta\langle K_A\rangle,
\]

with

\[
R^2=0.996590,
\qquad
{\rm Pearson}=0.998293.
\]

On the sphere:

\[
\delta S_A
=
-0.000081
+
0.989718\,\delta\langle K_A\rangle,
\]

with

\[
R^2=0.996592,
\qquad
{\rm Pearson}=0.998295.
\]

Thus:

\[
\boxed{
\delta S_A^{\rm graph}
\simeq
\delta\langle K_A^{\rm graph}\rangle.
}
\]

---

## 7. Result v2 — Modular source predicts curvature response

The second test connects the modular response to curvature:

\[
\delta W_{\rm loc}
\to
\delta\langle K_A\rangle
\to
\delta\kappa(r).
\]

The curvature response is localized around the source:

| Geometry | mean near/far ratio | median near/far ratio |
|---|---:|---:|
| flat_torus2d | 6.293 | 6.241 |
| sphere | 9.534 | 9.479 |

Most importantly, the modular source amplitude predicts the near-source curvature response:

| Geometry | \(R^2(|\delta\langle K_A\rangle|,\langle|\Delta\kappa|\rangle_{\rm near})\) | Pearson |
|---|---:|---:|
| flat_torus2d | 0.967229 | 0.983478 |
| sphere | 0.982252 | 0.991086 |

The signed relation is also nearly perfect, up to the sign convention:

| Geometry | signed \(R^2\) | signed Pearson |
|---|---:|---:|
| flat_torus2d | 0.998336 | -0.999168 |
| sphere | 0.998857 | -0.999428 |

Therefore:

\[
\boxed{
\delta\langle K_A\rangle
\text{ behaves as an effective source for the curvature response.}
}
\]

---

## 8. Interpretation

Paper 18 has two positive numerical pillars:

1. The graph modular first law holds with slope approximately \(0.989\) and \(R^2\simeq0.9966\).
2. The modular response predicts localized curvature response with \(R^2\simeq0.967\) to \(0.982\).

Thus:

\[
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r).
\]

This supports the source-side chain needed for the effective Einstein limit.

---

## 9. Limits

The current results are positive but not final.

1. The source is scalar/modular, not yet a full tensor \(T_{\mu\nu}\).
2. The density matrix \(\rho_A\) is a graph proxy.
3. The modular Hamiltonian \(K_A=-\log\rho_A\) is not yet derived from a true quantum subsystem density matrix in this test.
4. The relation is tested on controlled geometries, not yet on full quantum MI graphs.
5. The sign of the signed relation depends on the perturbation convention.
6. A continuum derivation via the modular first law remains to be written.

---

## 10. Folder structure

```text
papers/paper18_entanglement_stress_tensor/
  README.md
  paper18_entanglement_stress_tensor.tex

  scripts/
    paper18_graph_modular_first_law_v1.py
    paper18_modular_source_curvature_v2.py
    paper18_build_entanglement_stress_summary_v1.py

  results/
    graph_modular_first_law_v1/
    modular_source_curvature_v2/
    paper18_entanglement_stress_summary_v1/

  figures/
    # final figures copied from selected result folders

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

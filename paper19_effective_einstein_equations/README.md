# Paper 19 — Effective Einstein Equations from Entanglement Equilibrium

**Assembling the continuum limit of Bottom-Up Quantum Gravity**

---

## Status

Paper 19 assembles the three continuum arrows established in Papers 16--18:

\[
L_N\to\Delta_g,
\]

\[
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
\]

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

The target effective equation is:

\[
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}
}
\]

where:

- \(g_{\mu\nu}^{\rm ent}\) is the emergent entanglement metric;
- \(T_{\mu\nu}^{\rm ent}\) is the modular entanglement stress-energy source;
- \(\Lambda_{\rm ent}\) is an entanglement cosmological term;
- \(G_{\rm eff}\) is the effective gravitational coupling;
- \(\mathcal{H}_{\mu\nu}\) contains higher-order, nonlocal and finite-scale corrections.

In the smooth low-energy limit:

\[
\mathcal{H}_{\mu\nu}\to0.
\]

Then BuP reduces to the effective Einstein form:

\[
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

---

## 1. Starting point

The fundamental BuP equation is the variational equilibrium condition:

\[
\frac{\delta S_{\rm BuP}[W]}{\delta W_{ij}}=0.
\]

The discrete action is:

\[
S_{\rm BuP}[W]
=
\mathrm{Tr}\,L(W)^{-\beta[W]}
+
\sum_{ij}W_{ij}\kappa_{ij}[W]
+
\lambda\sum_{ij}W_{ij}d_{ij}^{2}
+
S_{\rm topo}[W].
\]

Each term has a continuum interpretation:

| Discrete term | Continuum role |
|---|---|
| \(\mathrm{Tr}\,L^{-\beta}\) | spectral geometry / Einstein--Hilbert sector |
| \(\sum_{ij}W_{ij}\kappa_{ij}\) | Ricci curvature sector |
| \(\lambda\sum_{ij}W_{ij}d_{ij}^2\) | locality and cosmological sector |
| \(S_{\rm topo}[W]\) | topology and global constraints |
| \(\delta\langle K_A\rangle\) | source / stress-energy sector |

The continuum limit of this equilibrium is the central object of Paper 19.

---

## 2. Pillar I — Spectral geometry from Paper 16

Paper 16 supports the continuum spectral limit:

\[
c_NL_N\to-\Delta_g.
\]

The low-spectrum tests gave:

| Geometry | Mean relative spectral error | \(\lambda_1^{\rm scaled}\) | Target |
|---|---:|---:|---:|
| \(S^1\) | 0.005072 | 1.006212 | 1.000000 |
| \(T^2\) | 0.026627 | 40.961270 | 39.478418 |

This establishes that the entanglement Laplacian reconstructs the low spectrum of the Laplace--Beltrami operator.

Thus:

\[
L(W)\to-\Delta_g.
\]

---

## 3. Pillar II — Ricci curvature from Paper 17

Paper 17 supports the Ricci curvature limit:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_N R_{\mu\nu}u^\mu u^\nu.
\]

The key results are:

\[
\Delta\bar\kappa_{\rm sphere-flat}=0.013383,
\]

\[
\Delta\left\langle\frac{\kappa}{\epsilon}\right\rangle_{\rm sphere-flat}=0.197719,
\]

\[
\Delta\left\langle\frac{\kappa}{\ell^2}\right\rangle_{\rm sphere-flat}=0.557041.
\]

The affine calibration at \(N=512\) gives:

\[
B_N=-0.301281,
\]

\[
C_N=0.391933,
\]

\[
A_N=\frac{1}{C_N}=2.551455.
\]

Therefore:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
-0.301
+
0.392\,R_{\mu\nu}u^\mu u^\nu.
\]

This establishes that the discrete Ollivier--Ricci curvature carries a calibrated mean Ricci signal.

---

## 4. Pillar III — Modular source from Paper 18

Paper 18 supports the source chain:

\[
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
\delta\kappa(r).
\]

The graph modular first law is validated:

| Geometry | Slope \(\delta S_A\) vs \(\delta\langle K_A\rangle\) | \(R^2\) |
|---|---:|---:|
| flat torus | 0.989481 | 0.996590 |
| sphere | 0.989718 | 0.996592 |

The modular source predicts the near-source curvature response:

| Geometry | near/far localization ratio | \(R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})\) |
|---|---:|---:|
| flat torus | 6.293007 | 0.967229 |
| sphere | 9.533529 | 0.982252 |

The signed Pearson correlations are:

\[
r=-0.999168
\]

for the flat torus, and

\[
r=-0.999428
\]

for the sphere.

Thus, \(\delta\langle K_A\rangle\) behaves as an effective modular source for curvature.

---

## 5. Effective continuum dictionary

Paper 19 uses the following dictionary:

| Discrete BuP object | Continuum object |
|---|---|
| \(W_{ij}=I(i:j)\) | entanglement metric \(g_{\mu\nu}^{\rm ent}\) |
| \(L(W)\) | \(-\Delta_g\) |
| \(\mathrm{Tr}\,L^{-\beta}\) | spectral gravitational action |
| \(\kappa_{ij}^{OR}\) | \(R_{\mu\nu}u^\mu u^\nu\) |
| \(\delta\langle K_A\rangle\) | modular source / \(T_{\mu\nu}^{\rm ent}\) |
| locality penalty | cosmological / infrared sector |
| \(S_{\rm topo}[W]\) | topological and global constraints |
| finite graph corrections | \(\mathcal{H}_{\mu\nu}\) |

---

## 6. Continuum equation

Combining the three pillars gives:

\[
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal{H}_{\mu\nu}
}
\]

where:

\[
G_{\mu\nu}
=
R_{\mu\nu}
-
\frac{1}{2}Rg_{\mu\nu}.
\]

The correction tensor \(\mathcal{H}_{\mu\nu}\) includes:

1. spectral higher-order corrections;
2. nonlocal entanglement corrections;
3. topology-induced corrections;
4. finite-\(N\) corrections;
5. deviations from smooth manifold behavior.

In the smooth low-energy limit:

\[
\mathcal{H}_{\mu\nu}\to0.
\]

Then:

\[
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

---

## 7. Interpretation

Paper 19 does not claim that Einstein gravity has been fully derived from first principles.

It establishes a controlled assembly:

\[
L_N\to\Delta_g,
\]

\[
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
\]

\[
\delta W_{\rm loc}\to T_{\mu\nu}^{\rm ent}.
\]

Together, these support the existence of an effective Einstein regime inside BuP.

The remaining task is analytic:

1. derive the coefficients \(G_{\rm eff}\) and \(\Lambda_{\rm ent}\);
2. control the correction tensor \(\mathcal{H}_{\mu\nu}\);
3. prove the continuum limit beyond controlled numerical geometries;
4. reconstruct a full tensor \(T_{\mu\nu}^{\rm ent}\), not only a modular scalar source.

---

## 8. Folder structure

```text
papers/paper19_effective_einstein_equations/
  README.md
  paper19_effective_einstein_equations.tex

  scripts/
    paper19_build_einstein_limit_summary_v1.py

  results/
    einstein_limit_summary_v1/
      paper19_einstein_limit_key_results.csv
      paper19_einstein_limit_summary.json
      paper19_einstein_limit_summary.md

  figures/
    # final synthesis figures

  notes/
    roadmap.md
    open_problems.md
    reproducibility.md

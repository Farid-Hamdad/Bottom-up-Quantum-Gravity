# Bottom-Up Quantum Gravity (BuP)

**Emergence of space, time, and gravity from quantum entanglement**  
Farid Hamdad — 2026

publication
F. Hamdad, "Contraintes cosmologiques sur une dimension spatiale émergente :
dimension variable, énergie noire et σ8 dans le cadre Bottom-Up Quantum Gravity",
HAL preprint hal-05590614v1 (2026).
https://hal.science/hal-05590614v1

## Overview

**Bottom-Up Quantum Gravity (BuP)** investigates a minimal hypothesis:

> Space, time, and gravity are not fundamental.  
> They emerge collectively from the entanglement structure of a finite global quantum state.

In this framework:

- entanglement defines connectivity,
- geometry emerges from informational structure,
- time appears via modular flow,
- gravity is interpreted as effective curvature of that structure.

The project combines a conceptual foundation in quantum information,
exact numerical simulations on finite systems, and phenomenological
extension to galactic rotation curves (SPARC dataset).

---

## Core Idea

Starting from a global pure quantum state $|\Psi\rangle$ with no
pre-existing space, an informational distance is defined:

$$
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max} + \epsilon}\right)
$$

where $I(i:j)$ is the mutual information between subsystems $i$ and $j$.

From this distance, the framework reconstructs:

- an entanglement graph and embedding,
- a discrete Laplacian,
- a spectral dimension $d_s(\tau)$,
- an effective radial profile $d(r)$.

---

## Results Summary

| Result | Status |
|--------|--------|
| Emergent rotation curves (SPARC, 175 galaxies) | strong numerical evidence |
| Emergent dimension and regimes | robust |
| Ollivier–Ricci curvature signal | statistically robust |
| Modular chaos diagnostics | consistent with chaotic regime |
| Lensing prediction | preliminary |
| Effective vacuum symmetry (U(1)-like) | numerically observed |
| Modular Hamiltonians as geometric probes | established numerically |
| Vacuum uniqueness and axis alignment | observed |
| Stable microscopic emergent dimension $d_s \approx 2.61 \pm 0.03$ | robust on $N=6$–$16$ |

---

## 1. Emergent Time: Modular Flow

For a subsystem $A$ with reduced density matrix
$\rho_A = \mathrm{Tr}_{\bar{A}} |\Psi\rangle\langle\Psi|$,
the modular Hamiltonian is:

$$
K_A = -\log \rho_A
$$

The modular flow

$$
O(\tau) = e^{iK_A\tau}\, O\, e^{-iK_A\tau}
$$

defines an intrinsic relational dynamics. Time is not a background
parameter but an informational property of the system.

---

## 2. Modular Chaos Diagnostics

The spectrum of $K_A$ is analyzed via:

- adjacent gap ratio statistics,
- Random Matrix Theory diagnostics,
- spectral form factor.

$$
\langle r \rangle = \left\langle \frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})} \right\rangle
$$

| Regime | Value |
|--------|------:|
| Poisson | 0.386 |
| GOE | 0.536 |
| GUE | 0.603 |

Numerical results give $\langle r \rangle \in [0.53, 0.59]$,
consistent with a chaotic regime.

---

## 3. Emergent Space

Geometry is reconstructed from mutual information:

$$
I(i:j) = S(\rho_i) + S(\rho_j) - S(\rho_{ij})
$$

$$
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max}} + \varepsilon\right)
$$

The geometry is not assumed — it is reconstructed from the
entanglement structure of $|\Psi\rangle$.

---

## 4. Ollivier–Ricci Discrete Curvature

The discrete curvature signal is the most statistically robust
numerical result of the project:

$$
\Delta\kappa_{\mathrm{edge}} \approx 0.0758
$$

- positive fraction: 1.00 across all seeds
- $p \sim 10^{-6}$
- near-far effect: $\Delta\kappa_{\mathrm{near-far}} > 0$

The informational deficit produces a detectable curvature signal
without any geometric postulate.

---

## 5. Emergent Horizon ($N = 16$)

A horizon-like region is detected characterized by:

- high entanglement entropy,
- low conductance,
- strong internal coupling.

This structure emerges without a pre-defined geometry.

---

## 6. Thermodynamic Gravity

The relation

$$
\delta S \approx \beta_{\mathrm{eff}}\,\delta E
$$

is numerically stable across configurations, consistent with
an emergent thermodynamic interpretation of gravity.

---

## 7. SPARC Extension (175 galaxies)

The BuP effective dimension model is tested against the SPARC
rotation curve dataset.

| Regime | $d_{\min}$ |
|--------|----------:|
| Massive galaxies (HIGH) | 2.49 |
| Dwarf galaxies (LOW) | 2.27 |

The regime structure is emergent, not imposed.

**Connection to microscopic structure.**

The scale-dependent dimension $d(r)$ observed in SPARC
can be interpreted as the large-scale manifestation of
state-dependent interaction structures observed in modular Hamiltonians.

This suggests a continuous chain:

$$
\text{quantum entanglement}
\;\rightarrow\; \text{modular structure}
\;\rightarrow\; \text{interaction graph}
\;\rightarrow\; \text{effective dimension } d(r).
$$

The agreement with SPARC data therefore is consistent with 
the microscopic mechanisms identified in this work.

---

## 8. Lensing Prediction (102 galaxies)

| Regime | $x_{\mathrm{peak}}$ |
|--------|--------------------:|
| LOW | 0.365 |
| HIGH | 0.498 |

- LOW regime: centralized lensing signal
- HIGH regime: extended lensing signal

This prediction constitutes a potential discriminant between
BuP and standard dark matter models.

---

## 9. Effective Vacuum Symmetry

The effective symmetry structure of the BuP vacuum is investigated
using a ground state of an anisotropic Hamiltonian with weak
local random fields:

$$
H_{\mathrm{BuP}} =
-J_{zz}\sum_i Z_i Z_{i+1}
-J_{xy}\sum_i (X_iX_{i+1}+Y_iY_{i+1})
+\sum_i (h_i^x X_i + h_i^y Y_i + h_i^z Z_i)
$$

**Method.** From the ground state, effective rotation matrices $R_{ij} \in SO(3)$
are extracted via polar SVD decomposition of the Pauli correlation tensor.
The Lie algebra generated by $\{\log R_{ij}\}$ is analyzed to classify
the effective symmetry group.

**Main result** ($N=8$, $J_{zz}=2.0$, $J_{xy}=0.5$, $h=0.1$):

- 95% of disorder realizations classified as U(1)-like
- mean order parameter $OP \approx 0.969$ vs $OP_{\mathrm{random}} \approx 0.619$

**Scan over disorder** $h_{\mathrm{scale}}$:
U(1)-like regime robust for $h \lesssim 0.2$,
transition near $h \sim 0.3$,
SO(3) regime recovered for $h \gtrsim 0.5$.

**Scan over anisotropy** $J_{zz}/J_{xy}$:
transition from SO(3)-like to U(1)-like near $J_{zz}/J_{xy} \approx 1$,
U(1)-like robust for $J_{zz}/J_{xy} \geq 2$.

**Interpretation.** This is not spontaneous symmetry breaking in the
field-theoretic sense — the local fields already break full SO(3)
explicitly. The result demonstrates robust effective axis selection:
the BuP vacuum is a structured medium, not an isotropic neutral state.

**kNN stability** ($k = 1, \ldots, 7$): U(1)-like signal
is stable for all values of $k$, confirming it is an intrinsic
property of the ground state rather than an artifact of the
graph reconstruction.

**Vacuum uniqueness**: the angular dispersion between vacuum axes
across disorder realizations is $\Delta_{\mathrm{vide}} = 0.056$ rad
(normalized: $0.072$, where $1$ corresponds to isotropic $S^2$).
The vacuum axis aligns with the $J_{zz}$ anisotropy direction
with $\langle|\cos z|\rangle = 0.999$.

---

## 10. Modular Hamiltonians as Geometric Probes

**Paper:** `bup_modular_paper_v3.tex`

The structure of modular Hamiltonians $K_A$ is investigated
systematically across 180 numerical configurations
($N \in \{8,10,12\}$, $|A| \in \{3,4,5,6\}$,
$\eta_{\mathrm{micro}} \in \{1.0, 1.1, 1.25, 1.5, 2.0\}$,
three state classes).

**Derivation chain:**

$$
|\Psi\rangle
\;\xrightarrow{\mathrm{Tr}_{\bar{A}}}\;
\rho_A
\;\xrightarrow{-\log}\;
K_A
\;\xrightarrow{\mathrm{fit}}\;
H_{\mathrm{eff}}^{1+2}
$$

**Result 1 — Quasi-locality:**

| State class | $\langle R^2 \rangle$ | $\sigma_{R^2}$ |
|-------------|----------------------:|---------------:|
| Ground | 0.9948 | 0.0099 |
| Thermal | 0.9994 | 0.0011 |
| Random | 0.9205 | 0.0741 |

**Result 2 — State-dependent anisotropy amplification:**

$$
\mathcal{A} = \frac{\eta_{\mathrm{eff}}}{\eta_{\mathrm{micro}}}
$$

| State class | $\langle\mathcal{A}\rangle$ | $\sigma_{\mathcal{A}}$ |
|-------------|----------------------------:|-----------------------:|
| Ground | 1.048 | 0.114 |
| Thermal | 1.107 | 0.125 |
| Random | 0.805 | 0.392 |

Structured states (ground, thermal) amplify microscopic anisotropy;
random states suppress it. The hierarchy

$$
\mathcal{A}_{\mathrm{thermal}} \gtrsim \mathcal{A}_{\mathrm{ground}} > 1 > \mathcal{A}_{\mathrm{random}}
$$

is state-dependent: random states constructed from the same
Hamiltonian show suppression, confirming that $K_A$ acts as a
coherence-sensitive filter.

**Result 3 — Monotonic amplification with $\eta_{\mathrm{micro}}$:***

For ground states, $\mathcal{A}$ increases monotonically from
$1.000$ at $\eta_{\mathrm{micro}} = 1.0$ to $1.144 \pm 0.190$
at $\eta_{\mathrm{micro}} = 2.0$.

**Implication for BuP.** This establishes a concrete numerical
bridge between entanglement structure and emergent effective dynamics,
replacing postulated Hamiltonians with operators derived from $|\Psi\rangle$.

**Link to modular chaos.**

The spectral properties of $K_A$ are consistent with chaotic behavior,
with $\langle r \rangle \in [0.53, 0.59]$, matching GOE/GUE predictions.

This establishes that modular Hamiltonians are not only quasi-local
operators, but also are consistent with a universal chaotic class,
reinforcing their interpretation as effective dynamical generators.

**Geometric interpretation.**

The extracted couplings $J_{ij}$ define an interaction graph,
which can be interpreted as an emergent geometry.

Locality of $K_A$ corresponds to short-range structure,
while anisotropy encodes directional deformation.

This provides a concrete realization of the BuP principle:

$$
\text{entanglement} \;\rightarrow\; \text{geometry} \;\rightarrow\; \text{effective dynamics.}
$$

---

## 11. Repository Structure

```
experiments/
    10_single_galaxy_ngc3198/         # Single galaxy fit (NGC 3198)
    20_v1_3_scan_ngc3198/             # Parameter scan
    30_multi_galaxy_regime_3_08_4_5/  # Multi-galaxy regime analysis
    40_micro_dimension_coherence/     # Stable emergent dimension on N=6–16

derive_d_r/                           # Effective dimension pipeline
results_bup_hybrid_multi/             # Multi-galaxy results
results_bup_hybrid_multi_rmax_fixed/  # Fixed r_max variant

jauge/                                # Gauge structure and vacuum symmetry
    bup_gauge_classification.py       # SO(3) group classification
    bup_symmetry_scan.py              # Comparative mechanism scan
    bup_vacuum_v3_standalone.py       # Vacuum effective symmetry (v3)
    bup_knn_scan.py                   # kNN stability scan
    bup_vacuum_uniqueness.py          # Vacuum uniqueness analysis
    bup_spontaneous_breaking.py       # Explicit vs spontaneous breaking test

modular/                              # Modular Hamiltonian study
    bup_modular_paper_v3.tex          # Paper (ready for submission)
    aggregate_rows.csv                # Full numerical dataset (180 configs)
    fig_amplification_by_state.png
    fig_amplification_vs_eta.png
    fig_fit_r2_by_state.png
    fig_locality_vs_eta.png
```

---

## 12. Physical Interpretation

Matter modifies the informational structure of $|\Psi\rangle$,
which modifies the emergent geometry,
which produces an effective gravitational effect.

No dark matter particle is postulated.
The gravitational excess may be interpreted as a geometric consequence of entanglement structure.

---

## 13. Limitations

- Effective model, not a complete relativistic theory
- Lensing via proxy observable, not full relativistic deflection
- Mass-to-light ratio dependence in SPARC fits
- Ground states of XXZ Hamiltonian used as BuP proxies —
  states generated by BuP variational dynamics remain to be studied

---

## Unified Picture

The BuP framework can be summarized as a single pipeline:

- Entanglement defines connectivity
- Connectivity defines geometry
- Geometry defines effective dynamics
- Effective dynamics reproduces gravitational phenomenology

This work provides numerical evidence at each step of this chain,
from quantum systems ($N \leq 16$ qubits) to astrophysical scales
(175 SPARC galaxies).

---

## 14. Perspectives

- Full relativistic lensing computation
- Confrontation with weak lensing surveys
- Cosmological extension
- Fundamental BuP dynamics from variational principle
- Connection between modular Hamiltonian structure and gauge fields

---

## 15. Stable emergent dimension across microscopic, galactic, and cosmological scales

A dedicated multi-size scan on structured BuP states for
$N = 6, 8, 9, 12, 16$ shows that the effective microscopic
dimension extracted from the entanglement-range estimator remains
remarkably stable:

$$
d_s(\lambda, N) \approx 2.61 \pm 0.03.
$$

Two robust facts emerge from the raw data:

- **stability in non-local coupling**: $d_s$ varies only weakly with the
  non-local parameter $\lambda \in [0,1]$;
- **stability in system size**: the same plateau persists from
  $N=6$ to $N=16$, with no evidence for a strong finite-size drift
  in the raw measurements.

This raw plateau is directly consistent with the galactic SPARC
inference $d_{\mathrm{eff}} = 2.644 \pm 0.295$, and remains compatible
with the cosmological late-time window preferred by DESI + Pantheon+.
The overlap therefore extends from qubit-scale simulations to galactic
and cosmological phenomenology over many orders of magnitude.

The corresponding experiment, scripts, raw scan, exploratory FSS output,
summary table, and figure are stored in:

- `experiments/40_micro_dimension_coherence/`

**Important note.**
The finite-size-scaling extrapolation stored in the exploratory
`fss_extrapolated.csv` file is not used as the primary conclusion.
In the current data regime, the raw measurements are substantially
more stable and physically meaningful than the fitted $N \to \infty$
extrapolation.

---

## Author

**Farid Hamdad**  
Bottom-Up Quantum Gravity, 2026  
hamdadfarid54@gmail.com  
[github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity](https://github.com/Farid-Hamdad/Bottom-up-Quantum-Gravity)

# Localized Modular Curvature Response — Results Synthesis

## 1. Status and scope

This document summarizes the frozen confirmatory results of:

`experiments/localized_modular_curvature_response`

The experiment was designed as a fresh N=12 follow-up to the previous N=10 modular-curvature response benchmark. Its purpose was to test whether the failure of the previous global Ollivier–Ricci curvature observable could be explained by a locality mismatch: a modular observable tied to a subsystem A was previously compared to a curvature mean over the full graph.

The new experiment therefore preregistered a curvature observable localized structurally around the same subsystem A=(0,1,2,3), without selecting edges from response data and without evaluating the new statistic on the previous N=10 raw curvature results before the fresh N=12 confirmatory run.

Frozen preregistration commit:

`00aceb88318f42fbc3c2a59f99d483a74b8bfabb`

Frozen preregistration SHA-256:

`a96703fee17d0bdd53bbe910aca687672c8144a3cbaf836e87be10547031ae8d`

Frozen implementation commit:

`8e43201bfb1dd4f56b20b5ef0151f937603e2e7b`

Frozen baseline commit:

`02a85a78875f210740eb1a7da53eef0adc86edbe`

Frozen confirmatory-results commit:

`9b1f8d105bf28a95906d2c0a8f954bdeb43520ca`

---

## 2. Frozen confirmatory design

The experiment used the same open-chain TFIM construction as the previous modular-time benchmark:

\[
H=-J\sum_i Z_iZ_{i+1}-h\sum_i X_i,
\]

with:

- N=12,
- J=h=1,
- P=6,
- Δt=0.35,
- nominal t=2.1,
- initial state |+>^12,
- subsystem A=(0,1,2,3).

The analyzed mixed state was the time-reversal mixture

\[
\rho_{TR}=\frac12\left(|\psi(+t)\rangle\langle\psi(+t)|+|\psi(-t)\rangle\langle\psi(-t)|\right).
\]

A local perturbation R_y^(q)(epsilon) was applied before both +t and -t evolutions.

The fresh confirmatory grid was

\[
q=0,\ldots,11,\qquad \epsilon\in\{0.075,0.15,0.30\},
\]

for a total of 36 conditions.

Only positive epsilon values were preregistered because the preceding N=10 experiment had shown numerical evenness under epsilon <-> -epsilon.

### 2.1 Local information observable

For the six internal pairs of A,

\[
A_{pairs}=\{(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)\},
\]

the local information observable was

\[
M_W=\frac{1}{6}\sum_{(i,j)\in A_{pairs}}W_{ij},
\]

with mutual information computed in natural-log units.

### 2.2 Local modular observable

From

\[
K_A=-\ln\rho_A,
\]

the modular response coefficients were converted to

\[
v_{ij}=\sqrt{A_{ij}},
\]

and the local modular observable was

\[
M_K=\frac{1}{6}\sum_{(i,j)\in A_{pairs}}v_{ij}.
\]

Strict positivity of rho_A was required; no spectral floor was introduced.

### 2.3 Frozen Ollivier–Ricci graph

The baseline mutual-information graph used density 0.333, corresponding to 22 support edges for the N=12 complete graph.

The global support E0 was frozen from the baseline and retained for every condition.

Ollivier–Ricci curvature was computed with alpha=0.5, using the historical internal backend and the same per-condition similarity normalization rule as the earlier BuP Ollivier–Ricci implementation.

The localized edge set was defined before unblinding as

\[
E_A=\{(u,v)\in E_0:u\in A\ \text{or}\ v\in A\}.
\]

No response magnitude, graph distance, fitted radius, or post hoc selection entered this definition.

For the frozen N=12 baseline,

\[
|E_0|=22,\qquad |E_A|=9.
\]

The nine localized support edges were

\[
(0,1),(0,2),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(3,5).
\]

The primary localized curvature observable was

\[
\kappa_A^{OR}=\frac{1}{|E_A|}\sum_{(u,v)\in E_A}\kappa_{uv}^{OR}.
\]

The full-graph curvature mean was retained only as a descriptive control.

---

## 3. Baseline pre-unblind validation

The baseline-only run passed before any confirmatory perturbation was executed.

Baseline values:

\[
\lambda_{min}(\rho_A)=3.1424600940622966\times10^{-4},
\]

\[
\mathrm{cond}(\rho_A)=1218.6512813308993,
\]

\[
M_W=0.11220107989331128,
\]

\[
M_K=0.89347671123521966,
\]

\[
\kappa_A^{OR}=0.40068672493142266,
\]

\[
\kappa_{global}^{OR}=0.33818424594450647.
\]

Structural gates:

- global support edge count: 22,
- local support edge count: 9,
- every local edge incident to A: true,
- global graph connected: true,
- number of density-matrix diagnostics: 79.

Numerical diagnostics:

\[
\max |\mathrm{Tr}\rho-1|=1.0436096431476471\times10^{-14},
\]

\[
\max \|\rho-\rho^\dagger\|=0,
\]

and the minimum reduced-state eigenvalue observed in the baseline diagnostics was

\[
3.142460094062331\times10^{-4}.
\]

No confirmatory-results directory existed at the time the baseline was frozen.

---

## 4. Confirmatory results

All 36 preregistered conditions completed with `status=ok`.

Across the full grid:

- all global supports matched the frozen baseline support,
- all local supports matched the frozen E_A,
- all graphs remained connected,
- global support count was always 22,
- local support count was always 9.

The minimum modular reduced-state eigenvalue across the confirmatory grid was

\[
1.295593932639744\times10^{-4},
\]

and the largest condition number was

\[
2943.094994866602.
\]

The modular rank remained defined for the full confirmatory grid.

### 4.1 H1 — modular / information alignment

The first frozen hypothesis tested ΔM_K against ΔM_W.

Observed Spearman correlation:

\[
\rho_S=-0.6120978120978121.
\]

The preregistered stratified permutation test used 100,000 permutations, permuting the site label q within each fixed epsilon stratum.

Result:

\[
p_{perm}=0.0002599974000259997.
\]

The permutation extreme count was 25 / 100,000.

Therefore H1 passed the frozen alpha=0.05 criterion.

### 4.2 H2 — modular / localized-curvature alignment

H2 was tested only because H1 passed, as required by the frozen hierarchical protocol.

The second hypothesis tested ΔM_K against Δkappa_A^OR.

Observed Spearman correlation:

\[
\rho_S=0.47078507078507076.
\]

Permutation result:

\[
p_{perm}=0.013579864201357986.
\]

The permutation extreme count was 1357 / 100,000.

Therefore H2 also passed the frozen alpha=0.05 criterion.

### 4.3 Frozen classification

Because both H1 and H2 passed, the preregistered classification is

`LOCALIZED_MODULAR_CURVATURE_RESPONSE_SUPPORTED`

This classification is fixed by the confirmatory protocol and is not altered by any diagnostic analysis reported below.

---

## 5. Preregistered robustness

The protocol preregistered leave-one-amplitude-level-out robustness checks as diagnostics.

### 5.1 H1 robustness

Full-grid value:

\[
\rho_S=-0.6120978120978121.
\]

Leaving out each amplitude level in turn gave:

| Excluded epsilon | Spearman rho_S |
|---:|---:|
| 0.075 | -0.6652173913043479 |
| 0.15 | -0.6669565217391304 |
| 0.30 | -0.5739130434782609 |

All three values have the same sign as the full result.

Thus the H1 association is not carried by a single perturbation amplitude.

### 5.2 H2 robustness

Full-grid value:

\[
\rho_S=0.47078507078507076.
\]

Leaving out each amplitude level gave:

| Excluded epsilon | Spearman rho_S |
|---:|---:|
| 0.075 | 0.47130434782608693 |
| 0.15 | 0.43217391304347824 |
| 0.30 | 0.4782608695652174 |

Again, all three values have the same sign as the full result.

Thus the localized modular-curvature association is also not carried by one amplitude level.

These leave-one-level-out calculations are robustness diagnostics, not additional confirmatory hypothesis tests.

---

## 6. Post-confirmatory diagnostics

Everything in this section was examined only after the frozen confirmatory result had been committed.

These diagnostics are descriptive and must not be used to redefine the confirmatory claim.

### 6.1 Spatial localization relative to subsystem A

The Spearman correlation between distance to A and absolute response amplitude was:

\[
\rho_S(d_A,|\Delta M_K|)=-0.8212277046551728,
\]

\[
\rho_S(d_A,|\Delta M_W|)=-0.8661724758172898,
\]

\[
\rho_S(d_A,|\Delta\kappa_A|)=-0.6954611958769682.
\]

By contrast, for the global curvature mean,

\[
\rho_S(d_A,|\Delta\kappa_{global}|)=0.04731028543380736.
\]

Thus the modular response, local information response, and localized curvature response all show strong attenuation with distance from A, whereas the full-graph curvature mean does not.

This is the principal post-confirmatory diagnostic supporting the locality interpretation.

### 6.2 Mean absolute response by distance

| Distance d_A | n | |ΔM_K| | |ΔM_W| | |Δkappa_A| | |Δkappa_global| |
|---:|---:|---:|---:|---:|---:|
| 0 | 12 | 0.0230948584078 | 8.93258254782e-4 | 7.92169854157e-4 | 3.63740845648e-4 |
| 1 | 3 | 0.0104899924471 | 5.06204012241e-4 | 7.26418575056e-4 | 6.28864461317e-4 |
| 2 | 3 | 0.00924390420896 | 1.05114427688e-4 | 7.02445779720e-4 | 4.28011210700e-5 |
| 3 | 3 | 0.00776503697243 | 1.36245429647e-6 | 4.37765024090e-4 | 2.04762461119e-5 |
| 4 | 3 | 0.00424466456174 | 1.93645877221e-5 | 3.25396260589e-4 | 5.64755503964e-4 |
| 5 | 3 | 2.12664725263e-4 | 8.94455793853e-7 | 1.32392707349e-4 | 3.38119102535e-4 |
| 6 | 3 | 1.00648087975e-6 | 3.76428459699e-9 | 2.13947628844e-5 | 2.50793255381e-4 |
| 7 | 3 | 1.40628249786e-15 | 5.08852219620e-17 | 1.81243487192e-6 | 2.27769851259e-4 |
| 8 | 3 | 8.88178419700e-16 | 2.03540887848e-16 | 1.66147355688e-8 | 9.00354040977e-4 |

The modular and information responses become numerically negligible at the most distant sites. The localized curvature response also decreases strongly, although it remains nonzero farther from A because Ollivier–Ricci curvature depends on transport geometry on the full frozen support graph.

### 6.3 Local curvature versus global curvature on the same N=12 conditions

Using the same 36 conditions:

\[
\rho_S(\Delta M_K,\Delta\kappa_A)=0.47078507078507076,
\]

whereas

\[
\rho_S(\Delta M_K,\Delta\kappa_{global})=0.20592020592020593.
\]

The local and global curvature responses themselves satisfy

\[
\rho_S(\Delta\kappa_A,\Delta\kappa_{global})=0.584041184041184.
\]

Therefore the localized curvature observable tracks the local modular response substantially better than the full-graph mean within the same N=12 dataset.

This comparison is descriptive because global curvature was not the confirmatory H2 target in this experiment.

### 6.4 Dependence on perturbation amplitude

| epsilon | |ΔM_K| | |ΔM_W| | |Δkappa_A| | |Δkappa_global| |
|---:|---:|---:|---:|---:|
| 0.075 | 0.00169763923401 | 5.09876921858e-5 | 6.60157481390e-5 | 5.21465961576e-5 |
| 0.15 | 0.00650128443329 | 2.02987202293e-4 | 2.63695203679e-4 | 2.09166688131e-4 |
| 0.30 | 0.0228852520895 | 7.97519285810e-4 | 1.04936944216e-3 | 8.45910957013e-4 |

Doubling epsilon produces approximately fourfold growth in several response amplitudes.

This is compatible with a leading even, approximately quadratic response over the tested range,

\[
|\Delta O|\sim O(\epsilon^2),
\]

and is consistent with the near-exact epsilon <-> -epsilon symmetry observed in the preceding N=10 experiment.

This scaling observation is post-confirmatory and was not itself a preregistered hypothesis test.

### 6.5 Pairwise modular and information profiles

The six internal pairs of A do not show a universal pairwise relationship between Δv_ij and ΔW_ij.

| Pair | rho_S(Δv_ij, ΔW_ij) |
|---|---:|
| (0,1) | 0.41931475786977124 |
| (0,2) | -0.1236078152942135 |
| (0,3) | -0.11445124127151067 |
| (1,2) | -0.7602973226922074 |
| (1,3) | -0.49043732880882757 |
| (2,3) | 0.8324108636890204 |

The response is therefore heterogeneous at the individual-pair level.

The confirmatory H1 result should be interpreted as a relationship between the preregistered aggregate local observables M_K and M_W, not as evidence for a universal microscopic identity v_ij proportional to W_ij.

### 6.6 Local curvature-edge profile

| Edge | Mean |Δkappa_uv| |
|---|---:|
| (0,1) | 8.68036357085e-4 |
| (0,2) | 2.40795449946e-3 |
| (1,2) | 8.22937949818e-4 |
| (1,3) | 1.49539584161e-3 |
| (1,4) | 2.44888488850e-3 |
| (2,3) | 1.31900717734e-3 |
| (2,4) | 1.30863385460e-3 |
| (3,4) | 8.08832996031e-4 |
| (3,5) | 1.17068188564e-3 |

The strongest average responses occur on both an internal edge, (0,2), and a boundary-crossing edge, (1,4).

Therefore the relevant localized geometric response is not restricted to edges entirely inside A. It includes the structural interface between A and its environment, exactly as intended by the preregistered definition of E_A.

---

## 7. Comparison with the previous N=10 experiment

The previous experiment used N=10, subsystem A=(0,1,2,3), and a global mean Ollivier–Ricci curvature response.

Its frozen H1 result was

\[
\rho_S(\Delta M_K,\Delta M_W)=-0.5780420716404548,
\]

with

\[
p_{perm}=9.99990000099999\times10^{-6}.
\]

Thus the modular-information association was strongly supported.

However, the previous global-curvature H2 result was

\[
\rho_S(\Delta M_K,\Delta\kappa_{global})=0.5087529178356222,
\]

with

\[
p_{perm}=0.06737932620673794,
\]

which failed the preregistered alpha=0.05 criterion.

Post-unblind analysis of the N=10 experiment further showed strong spatial localization of the modular and information responses but essentially no corresponding localization of the global curvature mean:

\[
\rho_S(d_A,|\Delta M_K|)\approx-0.692,
\]

\[
\rho_S(d_A,|\Delta M_W|)\approx-0.768,
\]

\[
\rho_S(d_A,|\Delta\kappa_{global}|)\approx-0.073.
\]

This motivated, but did not determine from the old data, the fresh localized-curvature experiment.

The N=12 experiment then preregistered E_A structurally and tested it on fresh conditions.

The resulting H2 passed:

\[
\rho_S(\Delta M_K,\Delta\kappa_A)=0.47078507078507076,
\]

\[
p_{perm}=0.013579864201357986.
\]

Furthermore, within the same N=12 dataset, the local curvature response shows strong spatial attenuation,

\[
\rho_S(d_A,|\Delta\kappa_A|)=-0.6954611958769682,
\]

whereas the global curvature response does not,

\[
\rho_S(d_A,|\Delta\kappa_{global}|)=0.04731028543380736.
\]

### 7.1 Correct interpretation of the N=10 / N=12 comparison

The N=10 and N=12 experiments are not a controlled single-variable ablation.

Between them, both the system size and perturbation grid changed, in addition to the primary curvature observable.

Therefore one must not claim that replacing global curvature by localized curvature alone caused the change from H2 failure to H2 success.

The stronger evidence for the locality interpretation comes from the internal N=12 comparison:

- same states,
- same graph construction,
- same 36 perturbations,
- local curvature versus global curvature.

On these identical conditions,

\[
\rho_S(\Delta M_K,\Delta\kappa_A)=0.471
\]

while

\[
\rho_S(\Delta M_K,\Delta\kappa_{global})=0.206,
\]

and only the localized curvature exhibits strong distance attenuation relative to A.

Thus the combined experiments support the hypothesis that matching the spatial domain of the modular and geometric observables is important.

---

## 8. Interpretation for BuP

The combined result supports the following microscopic response structure:

\[
|\Psi\rangle\longrightarrow\rho_A\longrightarrow K_A,
\]

together with a local information response

\[
\delta K_A\leftrightarrow\delta W_A,
\]

and a corresponding response of the localized discrete geometry constructed from the information graph,

\[
\delta W\longrightarrow\delta\kappa_A^{OR}[W].
\]

The experimentally supported response chain is therefore most accurately summarized as

\[
\boxed{\delta K_A\leftrightarrow\delta W_A\longrightarrow\delta\kappa_A^{OR}[W]}
\]

with a common local spatial structure relative to subsystem A.

This is relevant to the Bottom-Up Quantum Gravity program because it supplies evidence, within the present finite microscopic model, that:

1. modular response is associated with local information redistribution;
2. the information-derived geometric response is spatially compatible with the same subsystem when curvature is localized structurally;
3. a full-graph curvature average can erase or dilute this locality;
4. the relevant geometric object may therefore need to be attached to a region and its information-theoretic boundary rather than averaged indiscriminately over the full emergent graph.

The N=12 result consequently strengthens the microscopic information-to-geometry component of BuP.

It does not yet establish a continuum gravitational dynamics.

---

## 9. Claim boundaries

### 9.1 Ollivier–Ricci curvature is not independent of W

The graph geometry is constructed from the information weights:

\[
\kappa^{OR}=F(W).
\]

Therefore the experiment does not establish three statistically independent observables K, W, and kappa.

The positive H2 result cannot be interpreted as showing that K_A independently predicts curvature after controlling for the complete information graph.

The correct claim is that the modular response is associated with an information redistribution whose induced localized Ollivier–Ricci geometry responds coherently.

### 9.2 No pairwise law has been established

The heterogeneous pairwise correlations rule out interpreting the result as evidence for a simple universal relation such as

\[
v_{ij}=f(W_{ij})
\]

for every pair.

The supported relation is at the level of the frozen aggregate observables.

### 9.3 No continuum spacetime has been derived

This experiment does not derive:

- a continuum metric g_mu_nu,
- Einstein's field equations,
- Lorentzian spacetime,
- a Newtonian limit,
- a gravitational coupling constant,
- or a phenomenological prediction for a physical gravitational system.

The result concerns a finite quantum-information model and an emergent discrete graph geometry.

### 9.4 No causal direction between K_A and W_A is established

The observed Spearman association does not by itself establish K_A -> W_A or W_A -> K_A.

Both arise from the perturbed quantum state and may reflect a common underlying state-dependent structure.

The notation

\[
\delta K_A\leftrightarrow\delta W_A
\]

is therefore deliberately non-causal.

### 9.5 The local-curvature definition was not optimized on the N=12 response data

The edge set

\[
E_A=\{(u,v)\in E_0:u\in A\ \text{or}\ v\in A\}
\]

was fixed before unblinding.

No response threshold, fitted distance radius, selected subset of favorable edges, or post hoc tuning was used to define the confirmatory H2 observable.

This is central to the evidential value of the fresh N=12 result.

---

## 10. Main conclusion

The fresh N=12 preregistered experiment independently reproduced the modular-information response alignment and, unlike the preceding global-curvature test, also found a statistically significant association between the local modular response and a structurally preregistered localized Ollivier–Ricci curvature response.

The confirmatory results are

\[
\rho_S(\Delta M_K,\Delta M_W)=-0.61210,
\qquad
p_{perm}=2.60\times10^{-4},
\]

and

\[
\rho_S(\Delta M_K,\Delta\kappa_A^{OR})=0.47079,
\qquad
p_{perm}=0.01358.
\]

The spatial diagnostics show a common decay away from A for modularity, mutual information, and localized curvature, while the global curvature mean has essentially no corresponding distance dependence.

The strongest defensible conclusion is therefore:

\[
\boxed{\text{local modular response}\leftrightarrow\text{local information response}\rightarrow\text{localized information-derived geometric response}}
\]

within the tested finite TFIM model.

For BuP, this provides evidence that preserving locality across the modular, informational, and geometric levels is a necessary structural ingredient of the proposed bottom-up emergence mechanism.

The next scientific problem is no longer merely whether these local responses are associated, but whether a quantitative microscopic response law can be derived that predicts the induced geometric change from the modular/information dynamics without fitting the observed curvature response.

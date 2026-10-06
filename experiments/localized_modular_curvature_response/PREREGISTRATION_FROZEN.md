# BuP — Preregistration

## Localized Modular Response and Emergent Ollivier–Ricci Curvature

**Status:** PRE-EXECUTION — FROZEN BEFORE N=12 BASELINE
**Date:** 2026-10-06
**Repository:** `Bottom-up-Quantum-Gravity-bup-cosmology`
**Branch at freeze:** `publish_p29_p30_clean`
**HEAD before this preregistration:** `3dccfbc7c342c1b55428a1e3f4cb4139f8a51e2c`

## 0. Relationship to the previous experiment

This is a new experiment.

It is motivated by the completed experiment:

`experiments/modular_curvature_response`

whose frozen confirmatory classification was:

`MODULAR_INFORMATION_ALIGNMENT_SUPPORTED_CURVATURE_ALIGNMENT_NOT_SUPPORTED`.

The previous experiment found:

\[
\rho_S(\Delta M_K,\Delta M_W)
=
-0.5780420716404548,
\]

with

\[
p_{\rm perm}
=
9.99990000099999\times10^{-6},
\]

while the global mean Ollivier–Ricci curvature response gave

\[
\rho_S(\Delta M_K,\Delta\kappa_{\rm global})
=
0.5087529178356222,
\]

with

\[
p_{\rm perm}
=
0.06737932620673794.
\]

Post-unblind diagnostics showed a structural mismatch:

\[
\rho_S(d(q,A),|\Delta M_K|)
=
-0.6916855063469646,
\]

\[
\rho_S(d(q,A),|\Delta M_W|)
=
-0.7684625975514775,
\]

but

\[
\rho_S(d(q,A),|\Delta\kappa_{\rm global}|)
=
-0.07262697816643127.
\]

The global curvature scalar was also reflection symmetric to numerical precision under the reflection of the N=10 chain.

This motivates, but does not predetermine, the present hypothesis:

> a curvature observable localized relative to the same modular subsystem \(A\) may retain spatial information discarded by the global curvature mean.

The present experiment does not alter, repair, reclassify, or replace the previous experiment.

---

# 1. Contamination boundary

The previous experiment contains a raw file:

`experiments/modular_curvature_response/results/confirmatory_v1/curvature_edges.csv`

which would permit retrospective construction of many alternative curvature summaries.

To prevent post hoc selection:

**the newly defined localized curvature observable in this preregistration must not be evaluated on the previous N=10 confirmatory dataset before the new N=12 confirmatory run is frozen.**

In particular, before the new confirmatory run it is prohibited to compute the new localized curvature statistic from the previous:

- `curvature_edges.csv`;
- `per_condition.csv`;
- or any reconstructed N=10 result.

The definition below is selected structurally from the subsystem \(A\), not by optimizing its performance on prior raw data.

---

# 2. Scientific question

The experiment tests whether the microscopic response chain

\[
|\Psi\rangle
\rightarrow
\{\rho_A,\ W\}
\rightarrow
\{K_A,\ \kappa_A^{OR}[W]\}
\]

exhibits aligned perturbative response when the curvature observable is localized relative to the same subsystem used to define \(K_A\).

The empirical chain is

\[
\boxed{
\delta K_A
\leftrightarrow
\delta W_A
\rightarrow
\delta\kappa_A^{OR}[W].
}
\]

The arrow

\[
W\rightarrow\kappa_A^{OR}
\]

remains an analysis construction.

Therefore the experiment does not test whether curvature is statistically independent of \(W\).

It tests whether the local modular response is aligned with:

1. the local information-geometric response;
2. the response of a curvature statistic defined on the same spatial sector.

---

# 3. Freshness of the new benchmark

No N=12 result from the protocol defined below has been inspected at preregistration time.

The new confirmatory conditions differ from the previous experiment in two independent ways:

1. system size changes from \(N=10\) to \(N=12\);
2. perturbation amplitudes are new.

The frozen amplitudes are:

\[
\epsilon
\in
\{
0.075,\ 0.15,\ 0.30
\}.
\]

None of these amplitudes was used in the previous confirmatory grid.

Only positive amplitudes are used because the previous experiment established, post-unblind, that the measured N=10 responses were even under

\[
\epsilon\rightarrow-\epsilon
\]

to numerical precision.

This choice is made before any N=12 result is generated.

---

# 4. Frozen quantum system

The Hamiltonian remains the open-chain transverse-field Ising model:

\[
H
=
-J\sum_{i=0}^{N-2}Z_iZ_{i+1}
-h\sum_{i=0}^{N-1}X_i.
\]

Frozen parameters:

\[
N=12,
\qquad
J=1,
\qquad
h=1.
\]

Strang evolution parameters:

\[
P=6,
\qquad
\Delta t=0.35,
\qquad
t=P\Delta t=2.1.
\]

The microscopic topology is an open chain.

The unperturbed initial state is

\[
|\Psi_0\rangle
=
|+\rangle^{\otimes12}.
\]

The modular subsystem remains fixed to

\[
A=(0,1,2,3).
\]

No subsystem scan is allowed.

The choice of \(A\) is inherited unchanged from the previous experiment so that the new test probes cross-size generalization rather than selecting a new favorable subsystem.

---

# 5. Analysis state

For each initial condition, compute the two Strang-evolved branches:

\[
|\Psi(+t)\rangle,
\qquad
|\Psi(-t)\rangle.
\]

The analysis state is

\[
\rho_{\rm TR}
=
\frac12
\left[
|\Psi(+t)\rangle\langle\Psi(+t)|
+
|\Psi(-t)\rangle\langle\Psi(-t)|
\right].
\]

For the modular sector:

\[
\rho_A
=
\operatorname{Tr}_{\bar A}\rho_{\rm TR}.
\]

All reduced density matrices must be obtained from the mixed state, or equivalently by linear reduction of its two pure-state branches.

---

# 6. Fresh perturbation grid

For every site

\[
q\in\{0,1,\ldots,11\},
\]

apply

\[
R_y^{(q)}(\epsilon)
=
e^{-i\epsilon Y_q/2}
\]

to the initial state.

The perturbed initial state is

\[
|\Psi_{q,\epsilon}(0)\rangle
=
R_y^{(q)}(\epsilon)
|+\rangle^{\otimes12}.
\]

Frozen amplitudes:

\[
\epsilon
\in
\{
0.075,\ 0.15,\ 0.30
\}.
\]

Thus the confirmatory dataset contains:

\[
12\times3=36
\]

fresh perturbation conditions.

The unique baseline is

\[
\epsilon=0.
\]

The same perturbed initial state is used for the \(+t\) and \(-t\) branches.

No perturbation position or amplitude may be removed after unblinding.

---

# 7. Mutual-information geometry

All von Neumann entropies use the natural logarithm:

\[
S(\rho)
=
-\operatorname{Tr}(\rho\ln\rho).
\]

Thus mutual information is measured in nats:

\[
W_{ij}
=
S(\rho_i)
+
S(\rho_j)
-
S(\rho_{ij}).
\]

A full

\[
12\times12
\]

mutual-information matrix is computed for every condition.

## 7.1 Local information observable

Inside

\[
A=(0,1,2,3),
\]

the frozen six pairs are

\[
(0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
\]

Define

\[
M_W
=
\frac16
\sum_{i<j,\ i,j\in A}
W_{ij}.
\]

The response is

\[
\Delta M_W(c)
=
M_W(c)-M_W(0).
\]

This definition is unchanged from the previous experiment.

---

# 8. Modular sector

For the fixed subsystem \(A\),

\[
K_A
=
-\log\rho_A.
\]

The strict Paper-28 positivity gate is retained.

If

\[
\lambda_{\min}(\rho_A)\le0,
\]

then the modular observable is undefined.

No:

- eigenvalue floor;
- pseudolog;
- pseudoinverse;
- clipping;
- regularized logarithm

is permitted.

Let the eigenvalues of \(\rho_A\) be \(\lambda_n\), and define

\[
\kappa_n
=
-\log\lambda_n.
\]

As in Paper 28, standardize the modular spectrum:

\[
\widetilde K_A
=
\frac{
K_A-\langle\kappa\rangle I
}{
\operatorname{std}(\kappa)
}.
\]

For every pair \(i<j\) in \(A\),

\[
Q_{ij}^{ab}
=
[
[\widetilde K_A,\sigma_i^a],
\sigma_j^b
],
\]

and

\[
A_{ij}
=
\frac19
\sum_{a,b=X,Y,Z}
\operatorname{Tr}
\left[
\rho_A
Q_{ij}^{ab\dagger}
Q_{ij}^{ab}
\right].
\]

Define

\[
v_{ij}
=
\sqrt{A_{ij}}.
\]

The frozen scalar modular observable is

\[
M_K
=
\frac16
\sum_{i<j,\ i,j\in A}
v_{ij}.
\]

Its response is

\[
\Delta M_K(c)
=
M_K(c)-M_K(0).
\]

---

# 9. Global OR graph construction

The Ollivier–Ricci graph construction remains fixed to the historically audited internal pipeline.

Frozen parameters:

\[
\text{density}=0.333,
\]

\[
\alpha=0.5,
\]

\[
\text{backend}=\texttt{internal}.
\]

For

\[
N=12,
\]

the complete graph contains

\[
\binom{12}{2}=66
\]

possible edges.

The target frozen support size is therefore

\[
\operatorname{round}(0.333\times66)
=
22.
\]

The unperturbed baseline matrix \(W_0\) is used to construct the connected threshold graph:

1. maximum spanning tree;
2. strongest remaining MI edges;
3. stop at 22 edges.

Call this baseline support

\[
E_0.
\]

For every perturbation condition, the graph support remains exactly \(E_0\).

Thus the graph mode remains

\[
\texttt{frozen\_baseline}.
\]

For each condition, the historical per-condition normalization is retained:

\[
W_{ij}^{\rm norm}
=
\frac{
W_{ij}
}{
\max_{(u,v)\in E_0}W_{uv}
},
\]

and

\[
\ell_{ij}
=
\frac1{
\max(W_{ij}^{\rm norm},10^{-9})
}.
\]

The support is frozen but edge metric weights respond.

The internal OR implementation uses:

- weighted local measures;
- Dijkstra shortest paths using `length`;
- linear-programming Wasserstein distance;
- \(\alpha=0.5\);
- and

\[
\kappa_{uv}^{OR}
=
1-\frac{W_1(m_u,m_v)}{d(u,v)}.
\]

No `GraphRicciCurvature` fallback is allowed.

---

# 10. Frozen localized-curvature observable

This is the new central definition.

After the baseline support

\[
E_0
\]

has been constructed, define the local edge set:

\[
\boxed{
E_A
=
\{
(u,v)\in E_0:
u\in A
\ \text{or}\
v\in A
\}.
}
\]

Thus an edge belongs to \(E_A\) if and only if at least one endpoint lies in the same subsystem used to define \(K_A\).

This rule is purely structural.

It contains no observed curvature value, no perturbation response, no distance optimization, and no fitted radius.

The primary local curvature scalar is

\[
\boxed{
\kappa_A^{OR}
=
\frac1{|E_A|}
\sum_{(u,v)\in E_A}
\kappa_{uv}^{OR}.
}
\]

For every condition \(c\),

\[
\Delta\kappa_A(c)
=
\kappa_A^{OR}(c)
-
\kappa_A^{OR}(0).
\]

The edge set \(E_A\) is frozen once from the baseline support and must remain identical for all 36 perturbation conditions.

No local radius scan is allowed.

No alternative shell is allowed in the confirmatory analysis.

No edge may be selected based on its observed response.

---

# 11. Baseline validity gates

Before any of the 36 confirmatory perturbations are executed, a baseline-only run is allowed.

The baseline must satisfy all of the following:

1. state normalization residual

\[
\le10^{-12};
\]

2. modular positivity:

\[
\lambda_{\min}(\rho_A)>0;
\]

3. global OR graph is connected;

4. global support size is exactly

\[
|E_0|=22;
\]

5. local edge set is nonempty:

\[
|E_A|>0;
\]

6. every edge in \(E_A\) belongs to \(E_0\).

If one of these structural gates fails, the confirmatory grid must not be executed under this preregistration.

A new preregistration would then be required.

No change to \(A\), density, or the definition of \(E_A\) is permitted merely to make the baseline pass.

---

# 12. Global curvature control

For comparison only, retain the previous global curvature scalar:

\[
\kappa_{\rm global}^{OR}
=
\frac1{|E_0|}
\sum_{e\in E_0}
\kappa_e^{OR}.
\]

Its response is

\[
\Delta\kappa_{\rm global}.
\]

This quantity is **not** part of the primary confirmatory decision.

It is a descriptive control intended to determine whether localization changes the spatial behavior of the curvature observable.

The confirmatory status of the experiment must not depend on whether the local statistic outperforms the global statistic.

---

# 13. Confirmatory H1 — cross-size replication of modular / information alignment

The first confirmatory question is whether the modular / information relationship survives the move from \(N=10\) to fresh \(N=12\) data.

Frozen statistic:

\[
T_{KW}
=
\rho_S
\left(
\Delta M_K,
\Delta M_W
\right)
\]

across all 36 conditions.

The sign is not preregistered.

The test is two-sided.

H1 is necessary before any local-curvature claim is interpreted confirmatorily.

---

# 14. Confirmatory H2 — modular / localized-curvature alignment

Only if H1 passes the frozen criterion, test

\[
T_{K\kappa_A}
=
\rho_S
\left(
\Delta M_K,
\Delta\kappa_A
\right).
\]

The sign is not preregistered.

The test is two-sided.

H2 is the central new hypothesis.

No confirmatory test using a different curvature localization may substitute for H2.

---

# 15. Frozen permutation test

The 36 conditions are not treated as independent identically distributed replicates across perturbation strengths.

For each fixed amplitude

\[
\epsilon\in\{0.075,0.15,0.30\},
\]

there are 12 perturbation positions.

The null distribution is generated by independently permuting the 12 position labels of the modular response within each fixed-\(\epsilon\) stratum:

\[
\Delta M_K(q,\epsilon)
\rightarrow
\Delta M_K(\pi_\epsilon(q),\epsilon).
\]

Each amplitude stratum receives an independent position permutation.

This preserves perturbation strength exactly.

Frozen number of Monte-Carlo permutations:

\[
B=100000.
\]

Frozen random seed:

`20261006`

For observed statistic \(T_{\rm obs}\),

\[
p_{\rm perm}
=
\frac{
1+
\#\{
|T_b|\ge|T_{\rm obs}|
\}
}{
B+1
}.
\]

Frozen significance threshold:

\[
\alpha=0.05.
\]

The same permutation method is used for H1 and, if reached, H2.

---

# 16. Confirmatory decision rule

All 36 perturbation conditions must pass the modular positivity gate.

If any condition is modular-rank undefined:

`LOCALIZED_CURVATURE_CONFIRMATORY_TEST_UNDEFINED_DUE_TO_MODULAR_RANK`

and no failed condition may be deleted.

If H1 fails:

`CROSS_SIZE_MODULAR_INFORMATION_ALIGNMENT_NOT_SUPPORTED`

and H2 is not performed confirmatorily.

If H1 passes but H2 fails:

`MODULAR_INFORMATION_ALIGNMENT_SUPPORTED_LOCAL_CURVATURE_ALIGNMENT_NOT_SUPPORTED`

If H1 passes and H2 passes:

`LOCALIZED_MODULAR_CURVATURE_RESPONSE_SUPPORTED`

The successful local-curvature label therefore requires both:

\[
p_{KW}\le0.05
\]

and

\[
p_{K\kappa_A}\le0.05.
\]

---

# 17. Frozen robustness diagnostics

For every confirmatory relationship that passes, perform three leave-one-amplitude-level-out analyses.

Each removes one complete amplitude stratum:

- remove \(0.075\);
- remove \(0.15\);
- remove \(0.30\).

Each retained dataset contains 24 conditions.

Report:

- full Spearman coefficient;
- the three leave-one-level-out coefficients;
- whether all three retain the sign of the full coefficient;
- minimum coefficient;
- maximum coefficient.

These are robustness diagnostics only.

They do not alter the frozen confirmatory classification.

---

# 18. Secondary spatial diagnostics

The following are preregistered as descriptive only.

Define the open-chain distance from a perturbed site to \(A\):

\[
d(q,A)
=
\begin{cases}
0,&q\in A,\\
q-3,&q>3.
\end{cases}
\]

Report Spearman correlations between distance and:

\[
|\Delta M_K|,
\]

\[
|\Delta M_W|,
\]

\[
|\Delta\kappa_A|,
\]

and

\[
|\Delta\kappa_{\rm global}|.
\]

This directly checks whether localized curvature retains spatial dependence relative to \(A\).

No p-value from these secondary diagnostics changes the confirmatory result.

---

# 19. Secondary reflection diagnostic

For the N=12 open chain, reflection is

\[
q\rightarrow11-q.
\]

Report the maximum absolute reflection mismatch for:

\[
\Delta M_K,
\]

\[
\Delta M_W,
\]

\[
\Delta\kappa_A,
\]

and

\[
\Delta\kappa_{\rm global}.
\]

The expectation that localization may break the reflection symmetry of the global scalar is theoretical motivation only.

It is not a confirmatory success condition.

---

# 20. Secondary perturbation-strength diagnostic

For each amplitude

\[
\epsilon\in\{0.075,0.15,0.30\},
\]

report the mean absolute values of:

\[
|\Delta M_K|,
\quad
|\Delta M_W|,
\quad
|\Delta\kappa_A|,
\quad
|\Delta\kappa_{\rm global}|.
\]

This is descriptive.

No monotonicity requirement is frozen.

---

# 21. Secondary pairwise modular / information profile

Inside \(A\), preserve all six pair-level values:

\[
\Delta v_{ij}
\]

and

\[
\Delta W_{ij}.
\]

Pairwise profile correlations may be reported descriptively.

No pair subset may be promoted post hoc to the primary analysis.

---

# 22. Numerical integrity diagnostics

For every pure branch, record the state norm.

For every reduced density matrix used in the analysis, record:

- trace;
- Hermiticity residual;
- minimum eigenvalue;
- maximum eigenvalue;
- condition number over the positive spectrum.

For \(\rho_A\), strict positivity is mandatory.

For every graph condition, record:

- number of nodes;
- number of global support edges;
- graph connectedness;
- global support identity;
- local support identity;
- local support size;
- minimum and maximum similarity;
- minimum and maximum normalized similarity;
- minimum and maximum length.

The global support \(E_0\) and local support \(E_A\) must be identical across all 36 perturbation conditions.

---

# 23. Required outputs

At minimum preserve:

`config.json`

`provenance.json`

`baseline.json`

`per_condition.csv`

`modular_pairwise.csv`

`information_pairwise.csv`

`curvature_edges.csv`

`density_diagnostics.csv`

`permutation_results.json`

`summary.json`

The per-condition table must include at least:

- `q`
- `epsilon`
- `distance_to_A`
- `status`
- `lambda_min`
- `condition_number`
- `M_K`
- `delta_M_K`
- `M_W`
- `delta_M_W`
- `kappa_A`
- `delta_kappa_A`
- `kappa_global`
- `delta_kappa_global`
- global support count
- local support count
- global support identity status
- local support identity status
- graph connectedness

Failed or undefined conditions must remain in raw outputs.

---

# 24. Provenance requirements

Record:

- preregistration commit;
- preregistration SHA-256;
- implementation script SHA-256;
- source commit from the completed previous experiment;
- Paper 28 audited source blob;
- historical OR audited source blob;
- Python version;
- NumPy version;
- SciPy version;
- NetworkX version.

The implementation must not silently import scientific constants from mutable external files.

Frozen scientific constants must be explicit in the implementation.

---

# 25. Prohibited adaptive actions

After this preregistration is frozen, do not:

- inspect an N=12 perturbation response before implementation freeze;
- inspect the new localized curvature statistic on the previous N=10 raw confirmatory data before the N=12 run;
- change \(N\);
- change \(A\);
- change perturbation axis;
- change perturbation amplitudes;
- add negative amplitudes to increase sample count;
- remove an amplitude;
- remove a perturbation site;
- change graph density;
- change OR \(\alpha\);
- change OR backend;
- change graph mode;
- redefine \(E_A\);
- scan a local radius;
- select edges by observed curvature response;
- regularize a failed modular logarithm;
- change the MI logarithm base;
- change the permutation scheme;
- change the random seed;
- change the significance threshold;
- convert a descriptive global-curvature control into a replacement primary endpoint after seeing H2;
- reinterpret OR curvature as continuum Ricci curvature.

Any new localization definition after unblinding requires a separate future experiment.

---

# 26. Claim boundary

If H1 and H2 both pass, the allowed claim is:

> In a fresh N=12 finite TFIM benchmark, controlled microscopic perturbations produce aligned response between a state-derived modular observable, local mutual-information geometry, and an Ollivier–Ricci curvature statistic localized on baseline graph edges incident to the same subsystem \(A\).

This would support a localized microscopic bridge:

\[
|\Psi\rangle
\rightarrow
\rho_A
\rightarrow
K_A
\leftrightarrow
W_A
\rightarrow
\kappa_A^{OR}.
\]

It would **not** establish:

- \(K_A=T_{\mu\nu}\);
- \(\kappa_A^{OR}=R_{\mu\nu}\);
- an Einstein equation;
- continuum gravity;
- Lorentz invariance;
- relativistic causality;
- a conservation law;
- independence of curvature from \(W\);
- universality outside the tested finite TFIM systems.

If H2 fails, the correct conclusion is that the incident-edge localization defined here is not supported as the missing modular-curvature bridge.

No alternative localization may be substituted post hoc.

---

# 27. Execution order

The frozen execution sequence is:

1. save and commit this preregistration;
2. implement the N=12 benchmark without running any N=12 perturbation condition;
3. perform static/code audit;
4. execute baseline-only;
5. verify the baseline structural and numerical gates;
6. if any baseline gate fails, stop and create a new preregistration rather than altering this one;
7. freeze and commit the valid baseline;
8. execute the 36-condition confirmatory grid exactly once;
9. preserve all raw outputs immediately;
10. classify using the frozen hierarchical rule;
11. commit the confirmatory outputs;
12. only then perform secondary spatial, reflection, strength, and global-control analyses.

No scientific threshold or endpoint may change after step 8.

The previous N=10 confirmatory raw data must remain untouched throughout this sequence.

# BuP — Preregistration

## Microscopic Modular Response and Emergent Ollivier–Ricci Curvature

**Status:** PRE-EXECUTION — FROZEN BEFORE IMPLEMENTATION
**Date:** 2026-10-05
**Repository:** `Bottom-up-Quantum-Gravity-bup-cosmology`
**Branch at freeze:** `publish_p29_p30_clean`
**HEAD at freeze:** `e32105f4df3335f47f957e09f9aeae04f23ca838`

### Audited source provenance

Paper 28 modular benchmark:

`paper28_modular_time_geometry/scripts/paper28_reproduce_modular_time_v1.py`

Git blob:

`5fc9e845115a68135199e31f65324290ae0f210e`

Ollivier–Ricci reference:

`experiments/ollivier_ricci/scripts/bup_ollivier_ricci_local_response_scan_radius_fit_v1_5.py`

Git blob:

`3fd6eceb9c107075c80abf9eefd8fa73db1b69cd`

No result from the experiment described below has been generated at the time of freezing this preregistration.

---

# 1. Scientific question

The experiment tests whether a controlled microscopic perturbation produces an aligned response in:

\[
\rho_A
\rightarrow
K_A=-\log\rho_A,
\]

the mutual-information geometry

\[
W_{ij}=I(i:j),
\]

and the Ollivier–Ricci curvature of the graph constructed from \(W\).

The empirical chain tested is

\[
\boxed{
\delta K_A
\leftrightarrow
\delta W
\rightarrow
\delta\kappa^{OR}[W].
}
\]

The arrow

\[
W\rightarrow\kappa^{OR}
\]

is a **construction of the analysis pipeline**, not an independently inferred physical relationship.

Therefore this experiment will **not** test whether \(K_A\) provides information about \(\kappa^{OR}\) beyond the complete matrix \(W\).

The empirical questions are instead:

1. Does the modular response covary with the information-geometric response?
2. Does the modular response covary with the curvature response induced by that information geometry?
3. Is this alignment robust across perturbation positions and perturbation strengths?

---

# 2. Frozen quantum system

The benchmark inherits the finite TFIM protocol from Paper 28.

\[
H=
-J\sum_{i=0}^{N-2} Z_iZ_{i+1}
-h\sum_{i=0}^{N-1}X_i.
\]

Frozen parameters:

\[
N=10,
\qquad
J=1,
\qquad
h=1,
\]

\[
P=6,
\qquad
\Delta t=0.35,
\qquad
t=P\Delta t=2.1.
\]

The microscopic topology is an **open chain**, not a ring.

Initial unperturbed state:

\[
|\Psi_0\rangle=|+\rangle^{\otimes 10}.
\]

The modular subsystem is fixed to

\[
A=(0,1,2,3).
\]

No subsystem scan or adaptive selection is allowed.

---

# 3. Analysis state

For each initial condition, two Strang evolutions are performed:

\[
|\Psi(+t)\rangle,
\qquad
|\Psi(-t)\rangle.
\]

The analysis state is the Paper 28 time-reversal mixture

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

The modular reduced density matrix is

\[
\rho_A
=
\operatorname{Tr}_{\bar A}
\rho_{\rm TR}.
\]

This is a genuine reduced density matrix of a mixed quantum state.

It is not replaced by a graph proxy.

---

# 4. Controlled perturbation

The perturbation is applied to the initial state before TFIM evolution.

For site

\[
q\in\{0,\ldots,9\},
\]

define

\[
|\Psi_{q,\epsilon}(0)\rangle
=
R_y^{(q)}(\epsilon)
|+\rangle^{\otimes10},
\]

with

\[
R_y(\epsilon)
=
e^{-i\epsilon Y/2}.
\]

Frozen signed perturbation amplitudes, in radians:

\[
\epsilon\in
\{
-0.20,
-0.10,
-0.05,
+0.05,
+0.10,
+0.20
\}.
\]

Thus the confirmatory dataset contains

\[
10\times6=60
\]

perturbed conditions.

The baseline is

\[
\epsilon=0.
\]

The same perturbed initial state is used for the \(+t\) and \(-t\) branches before forming the time-reversal mixture.

No perturbation strength or position may be removed because of the observed result.

---

# 5. Mutual-information geometry

All von Neumann entropies use the natural logarithm:

\[
S(\rho)
=
-\operatorname{Tr}(\rho\ln\rho).
\]

Therefore all mutual information values are expressed in **nats**:

\[
W_{ij}
=
S(\rho_i)
+
S(\rho_j)
-
S(\rho_{ij}).
\]

The pairwise reduced density matrices must be obtained from the mixed analysis state \(\rho_{\rm TR}\), or equivalently by linear reduction of its two pure-state branches.

The pure-state SVD entropy implementation from the historical Ollivier–Ricci script must **not** be used for this benchmark.

Two uses of \(W\) are distinguished.

### Local Paper-28 sector

Inside \(A=(0,1,2,3)\), the six pairs are

\[
(0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
\]

Define

\[
M_W
=
\frac{1}{6}
\sum_{i<j,\ i,j\in A}
W_{ij}.
\]

### Global curvature sector

The full

\[
10\times10
\]

mutual-information matrix is used to construct the Ollivier–Ricci graph.

---

# 6. Modular sector

For the fixed subsystem \(A\),

\[
K_A=-\log\rho_A.
\]

The Paper 28 positivity gate is retained.

If the eigenvalues of \(\rho_A\) are

\[
\lambda_n,
\]

then the modular Hamiltonian is defined only when

\[
\lambda_{\min}>0.
\]

No eigenvalue floor, pseudolog, clipping, regularized logarithm, or pseudoinverse is allowed.

Define

\[
\kappa_n=-\log\lambda_n.
\]

As in Paper 28, the eigenvalues are centered and standardized:

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
A_{ij}
=
\frac19
\sum_{a,b=X,Y,Z}
\operatorname{Tr}
\left[
\rho_A
Q_{ij}^{ab\dagger}
Q_{ij}^{ab}
\right],
\]

where

\[
Q_{ij}^{ab}
=
[
[\widetilde K_A,\sigma_i^a],
\sigma_j^b
].
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

No pair subset may be selected after viewing the results.

---

# 7. Ollivier–Ricci graph

The graph construction is inherited from the audited Ollivier–Ricci pipeline, except that \(W\) is computed from the mixed TFIM state defined above.

Frozen graph parameters:

\[
\text{density}=0.333,
\]

\[
\alpha=0.5,
\]

\[
\text{backend}=\texttt{internal}.
\]

`backend=auto` is prohibited for the confirmatory run.

For \(N=10\), density \(0.333\) gives a target of 15 graph edges.

The unperturbed baseline matrix \(W_0\) is used to construct the connected threshold graph:

1. maximum spanning tree;
2. strongest remaining mutual-information edges until the frozen target density is reached.

This baseline edge set is

\[
E_0.
\]

For every perturbed condition, the graph support remains exactly \(E_0\).

Thus the confirmatory graph mode is

\[
\boxed{\texttt{frozen\_baseline}}.
\]

No free, union, intersection, or perturbed-support graph is used in the primary test.

For each condition, edge similarities are normalized according to the historical pipeline:

\[
W_{ij}^{\rm norm}
=
\frac{W_{ij}}
{\max_{(u,v)\in E_0}W_{uv}},
\]

and edge length is

\[
\ell_{ij}
=
\frac{1}
{\max(W_{ij}^{\rm norm},10^{-9})}.
\]

Thus the support is frozen, while the metric weights respond to the perturbed state.

The internal Ollivier–Ricci implementation uses weighted shortest paths and linear-programming Wasserstein distance.

No GraphRicciCurvature fallback is permitted.

---

# 8. Frozen response observables

For every condition

\[
c=(q,\epsilon),
\]

define differences relative to the unique unperturbed baseline.

### Modular response

\[
\Delta M_K(c)
=
M_K(c)-M_K(0).
\]

### Information response

\[
\Delta M_W(c)
=
M_W(c)-M_W(0).
\]

### Curvature response

Let

\[
\bar\kappa_{E_0}
=
\frac1{|E_0|}
\sum_{e\in E_0}
\kappa_e^{OR}.
\]

Then

\[
\Delta\kappa_{\rm edge}(c)
=
\bar\kappa_{E_0}(c)
-
\bar\kappa_{E_0}(0).
\]

These definitions are frozen before observing any response values.

---

# 9. Confirmatory hypothesis hierarchy

Two relationships are tested hierarchically.

## H1 — Modular / information alignment

Primary statistic:

\[
T_{KW}
=
\rho_S
\left(
\Delta M_K,
\Delta M_W
\right),
\]

where \(\rho_S\) is Spearman rank correlation across the 60 frozen perturbation conditions.

The sign is **not preregistered**.

The test is two-sided.

## H2 — Modular / curvature alignment

Only if H1 passes the confirmatory criterion, test

\[
T_{K\kappa}
=
\rho_S
\left(
\Delta M_K,
\Delta\kappa_{\rm edge}
\right).
\]

The sign is again not preregistered.

The test is two-sided.

This hierarchical order avoids treating the deterministic transformation

\[
W\rightarrow\kappa^{OR}[W]
\]

as an independent third source of evidence.

---

# 10. Permutation test

The 60 conditions must not be treated as 60 independent identically distributed experimental replicates.

For each signed perturbation amplitude \(\epsilon\), there are ten positions \(q\).

The null distribution is generated by independently permuting the ten position labels of the modular response **within each fixed signed \(\epsilon\) stratum**.

Therefore perturbation strength is preserved exactly.

For each permutation:

\[
\Delta M_K(q,\epsilon)
\rightarrow
\Delta M_K(\pi_\epsilon(q),\epsilon),
\]

with an independent permutation \(\pi_\epsilon\) for every signed value of \(\epsilon\).

Frozen number of Monte-Carlo permutations:

\[
B=100000.
\]

Frozen random seed:

`20261005`

For an observed statistic \(T_{\rm obs}\),

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

No alternative permutation scheme may replace this test after unblinding.

---

# 11. Confirmatory decision rule

The complete microscopic bridge receives the label

`MODULAR_CURVATURE_RESPONSE_SUPPORTED`

only if:

1. all 60 modular conditions pass the positivity gate;
2. H1 gives
   \[
   p_{KW}\le0.05;
   \]
3. H2 subsequently gives
   \[
   p_{K\kappa}\le0.05.
   \]

If H1 fails:

`MODULAR_INFORMATION_ALIGNMENT_NOT_SUPPORTED`

and H2 is reported descriptively but not treated as confirmatory evidence.

If H1 passes but H2 fails:

`MODULAR_INFORMATION_ALIGNMENT_SUPPORTED_CURVATURE_ALIGNMENT_NOT_SUPPORTED`.

If any of the 60 conditions has undefined \(K_A\) because of failure of the positivity gate:

`CONFIRMATORY_TEST_UNDEFINED_DUE_TO_MODULAR_RANK`.

No failed condition may simply be deleted from the confirmatory dataset.

---

# 12. Frozen robustness control

For each successful confirmatory relationship, perform six leave-one-\(\epsilon\)-level-out analyses.

Each analysis removes one complete signed perturbation stratum and retains the other five.

These are robustness diagnostics only.

They do not alter the primary decision rule.

Report:

- the six leave-one-level-out Spearman coefficients;
- whether all six retain the sign of the full statistic;
- minimum and maximum coefficient.

No threshold from these diagnostics affects the frozen confirmatory classification.

---

# 13. Secondary descriptive analyses

The following are explicitly secondary.

### Pairwise profile response

Compare the six-dimensional vectors

\[
\Delta v_{ij}
\]

and

\[
\Delta W_{ij}
\]

inside subsystem \(A\).

### Distance from modular subsystem

Define open-chain distance from perturbed site \(q\) to \(A\):

\[
d(q,A)=
\begin{cases}
0,&q\in A,\\
q-3,&q>3.
\end{cases}
\]

Study the dependence of

\[
|\Delta M_K|,
\quad
|\Delta M_W|,
\quad
|\Delta\kappa_{\rm edge}|
\]

on this distance.

No ring distance is allowed.

### Signed perturbation symmetry

Compare responses for

\[
+\epsilon
\quad\text{and}\quad
-\epsilon.
\]

### Strength response

Inspect dependence on

\[
|\epsilon|.
\]

These analyses are descriptive unless a separate preregistration is frozen before they are promoted to confirmatory status.

---

# 14. Numerical integrity gates

State normalization must satisfy

\[
\left|
\langle\Psi|\Psi\rangle-1
\right|
\le10^{-12}.
\]

For every reduced density matrix, record:

- trace;
- Hermiticity residual;
- minimum eigenvalue;
- maximum eigenvalue;
- condition number.

For the modular subsystem, no regularization is allowed if positivity fails.

For every graph condition, record:

- number of nodes;
- number of edges;
- connectedness;
- frozen support identity;
- minimum and maximum edge weight;
- minimum and maximum edge length.

The frozen graph support must be identical across all 60 perturbed conditions.

---

# 15. Required raw outputs

The implementation must preserve enough information to reconstruct every aggregate statistic.

At minimum:

`config.json`

`provenance.json`

`baseline.json`

`per_condition.csv`

`modular_pairwise.csv`

`information_pairwise.csv`

`curvature_edges.csv`

`permutation_results.json`

`summary.json`

The per-condition table must include at least:

- `q`
- `epsilon`
- `lambda_min`
- `condition_number`
- `M_K`
- `delta_M_K`
- `M_W`
- `delta_M_W`
- `kappa_edge_mean`
- `delta_kappa_edge`
- status fields

No raw failed or null condition may be silently omitted.

---

# 16. Prohibited adaptive actions

After execution begins, do not:

- change \(N\);
- change subsystem \(A\);
- select a different perturbation site subset;
- remove inconvenient \(\epsilon\) values;
- choose only positive or negative rotations;
- change graph density;
- change OR \(\alpha\);
- switch OR backend;
- choose a graph mode based on the result;
- alter the MI logarithm base;
- regularize rank-deficient \(\rho_A\);
- choose favorable pairs inside \(A\);
- change the permutation scheme;
- change the significance threshold;
- redefine the primary scalar observables;
- reinterpret \(\kappa^{OR}\) as a continuum Ricci tensor.

Any extension requires an explicit amendment recorded before examining the corresponding new result.

---

# 17. Claim boundary

A positive result would support the statement:

> In the frozen finite TFIM benchmark, controlled microscopic perturbations produce an aligned response between the state-derived modular sector and the mutual-information geometry, with the modular response also tracking the Ollivier–Ricci curvature response induced by that geometry.

It would **not** establish:

- \(K_A=T_{\mu\nu}\);
- \(\kappa^{OR}=R_{\mu\nu}\);
- an Einstein equation;
- a continuum limit;
- a relativistic causal structure;
- a fundamental conservation law;
- that modular information independently determines curvature beyond \(W\);
- universality outside the tested finite TFIM protocol.

The logical interpretation remains

\[
\boxed{
\text{microscopic state}
\rightarrow
\{W,K_A\},
\qquad
W\rightarrow\kappa^{OR},
}
\]

with the experiment testing whether the two state-derived sectors respond coherently under the same controlled perturbation.

---

# 18. Execution order

The frozen sequence is:

1. save and commit this preregistration;
2. implement the benchmark without running the 60-condition analysis;
3. perform static/code audit;
4. perform baseline-only integrity checks;
5. freeze any necessary pre-unblinding amendment;
6. execute the full confirmatory grid exactly once;
7. preserve all outputs;
8. classify using the frozen rule;
9. only afterward perform secondary interpretation.

No scientific threshold may be changed after step 6.

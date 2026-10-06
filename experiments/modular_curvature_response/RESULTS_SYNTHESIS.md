# Modular Curvature Response — Results Synthesis

**Experiment:** `experiments/modular_curvature_response`
**Preregistration commit:** `c90bcc35df4f303149e60eaf9a4bb98c025128ff`
**Implementation commit:** `4cc3e9aa3e36aef0cbce966c280557d75614628a`
**Pre-unblind baseline commit:** `2de8145854a19cee90b3f69d501291b645be022f`
**Frozen confirmatory-results commit:** `2aab1eb557010e4d6edc4305ffd992f40d41159a`

## 1. Question

The experiment tested the frozen response chain

\[
\delta K_A \leftrightarrow \delta W
ightarrow \delta \kappa^{OR}[W],
\]

with the important interpretation that the arrow

\[
W
ightarrow\kappa^{OR}
\]

is part of the analysis construction, not independent evidence.

The fixed modular subsystem was

\[
A=(0,1,2,3).
\]

The confirmatory perturbation grid contained 60 conditions:

\[
q\in\{0,\ldots,9\},
\]

and

\[
\epsilon\in\{-0.20,-0.10,-0.05,+0.05,+0.10,+0.20\}.
\]

## 2. Baseline integrity

The pre-unblind baseline passed all frozen integrity checks.

Key baseline values:

- `lambda_min = 0.00031424600940622896`
- `condition_number = 1218.6512813309037`
- `M_W = 0.11220107989331`
- `M_K = 0.89347671123521843`
- `kappa_edge_mean = 0.28248142014698069`
- frozen OR support: 15 edges
- graph connected: true

The reduced-density-matrix audit contained 56 diagnostics:

- 10 one-site reduced matrices;
- 45 two-site reduced matrices;
- the modular reduced matrix \(
ho_A\).

The maximum trace deviation was

\[
8.88	imes 10^{-15},
\]

and the maximum Hermiticity residual was

\[
5.56	imes10^{-17}.
\]

No pre-unblind amendment was required.

## 3. Confirmatory execution

All 60 frozen perturbation conditions completed successfully.

- `N_CONDITIONS = 60`
- `N_OK = 60`
- frozen support matched in every condition
- all OR graphs remained connected
- minimum modular eigenvalue across the grid:
  `0.00021234469988031744`
- maximum modular condition number:
  `1800.0197598228376`

No condition was removed.

## 4. Confirmatory H1 — modular / information alignment

Frozen statistic:

\[
T_{KW}
=

ho_S(\Delta M_K,\Delta M_W).
\]

Observed result:

\[

ho_S=-0.5780420716404548.
\]

Frozen stratified permutation test:

- permutations: 100000
- seed: 20261005
- extreme count: 0

Therefore

\[
p_{
m perm}
=
9.99990000099999	imes10^{-6}.
\]

H1 passed the preregistered threshold.

The confirmatory interpretation is:

\[
oxed{\delta K_A\leftrightarrow\delta W_A\ 	ext{supported in the frozen benchmark}.}
\]

The negative sign was not preregistered and is therefore interpreted as the observed direction of association, not as a post hoc success criterion.

## 5. Confirmatory H2 — modular / global OR-curvature alignment

Because H1 passed, the hierarchical H2 test was performed:

\[
T_{K\kappa}
=

ho_S(\Delta M_K,\Delta\kappa_{
m edge}).
\]

Observed result:

\[

ho_S=0.5087529178356222.
\]

Frozen stratified permutation test:

- permutations: 100000
- seed: 20261005
- extreme count: 6737

Therefore

\[
p_{
m perm}
=
0.06737932620673794.
\]

This does not pass the frozen threshold

\[
lpha=0.05.
\]

The frozen classification is therefore:

`MODULAR_INFORMATION_ALIGNMENT_SUPPORTED_CURVATURE_ALIGNMENT_NOT_SUPPORTED`

This classification must not be changed by the post-unblind diagnostics below.

## 6. Frozen H1 robustness control

The preregistered leave-one-signed-epsilon-level-out analysis retained the H1 sign in all six cases.

Observed Spearman coefficients:

- exclude `-0.20`: `-0.5798194244401575`
- exclude `-0.10`: `-0.5843738001359612`
- exclude `-0.05`: `-0.5892042459851088`
- exclude `+0.05`: `-0.5892042459851088`
- exclude `+0.10`: `-0.5843738001359612`
- exclude `+0.20`: `-0.5798194244401575`

Range:

\[
-0.5892042459851088
\le

ho_S
\le
-0.5798194244401575.
\]

Thus the H1 association is not carried by a single signed perturbation level.

---

# Post-unblind descriptive diagnostics

Everything from this section onward is explicitly post-unblind and does not alter the frozen confirmatory classification.

## 7. Exact response symmetry under \(\epsilon
ightarrow-\epsilon\)

The maximum differences between the responses at \(+\epsilon\) and \(-\epsilon\) were:

\[
\max|\Delta M_K(+\epsilon)-\Delta M_K(-\epsilon)|
=
3.6637359812630166	imes10^{-15},
\]

\[
\max|\Delta M_W(+\epsilon)-\Delta M_W(-\epsilon)|
=
3.7470027081099033	imes10^{-16},
\]

\[
\max|\Delta\kappa(+\epsilon)-\Delta\kappa(-\epsilon)|
=
9.43689570931383	imes10^{-16}.
\]

Thus, to numerical precision, the measured response is even in the perturbation amplitude:

\[
R(+\epsilon)=R(-\epsilon).
\]

The six signed perturbation levels therefore reduce to three distinct response-strength levels in practice.

## 8. Symmetry-aware post hoc permutation diagnostic

The \(+\epsilon\) and \(-\epsilon\) conditions were averaged pairwise, producing 30 post hoc conditions:

\[
10\ q	ext{-positions}	imes3\ |\epsilon|	ext{-levels}.
\]

Using the same 100000-permutation procedure and seed for comparability:

### H1 symmetry-aware diagnostic

\[

ho_S=-0.5786429365962181,
\]

\[
p_{
m perm}=0.0008199918000819992.
\]

Thus the modular / information alignment remains statistically strong after removal of the signed duplication.

### H2 symmetry-aware diagnostic

\[

ho_S=0.51412680756396,
\]

\[
p_{
m perm}=0.14419855801441986.
\]

Thus the global-curvature alignment remains unsupported and becomes less compelling under the symmetry-aware diagnostic.

These are post hoc robustness diagnostics, not replacements for the frozen confirmatory statistics.

## 9. Spatial localization relative to \(A\)

Define the open-chain distance to the fixed subsystem \(A=(0,1,2,3)\):

\[
d(q,A)=
egin{cases}
0,&q\in A,\
q-3,&q>3.
\end{cases}
\]

Using the 30 symmetry-aware conditions, the Spearman correlations between distance and absolute response are:

\[

ho_S(d,|\Delta M_K|)
=
-0.6916855063469646,
\]

\[

ho_S(d,|\Delta M_W|)
=
-0.7684625975514775,
\]

whereas

\[

ho_S(d,|\Delta\kappa_{
m edge}|)
=
-0.07262697816643127.
\]

Therefore the modular and local mutual-information responses are strongly localized around \(A\), while the global mean OR-curvature response is essentially not localized with respect to \(A\).

## 10. Mean absolute response by distance

Symmetry-aware mean absolute responses:

| distance | mean \(|\Delta M_K|\) | mean \(|\Delta M_W|\) | mean \(|\Delta\kappa_{
m edge}|\) |
|---:|---:|---:|---:|
| 0 | 0.0112192480307 | 0.000401819144117 | 0.000151664056007 |
| 1 | 0.00479739447223 | 0.000226443769598 | 0.000237795099639 |
| 2 | 0.00420954882427 | 0.0000457550726836 | 0.000237795099639 |
| 3 | 0.00354900887648 | 0.000000652285958032 | 0.00024723915985 |
| 4 | 0.00192424143567 | 0.00000879202407295 | 0.000128071403511 |
| 5 | 0.0000958610993379 | 0.000000403620902080 | 0.00000334132433393 |
| 6 | 0.000000453431616846 | 0.00000000169584928739 | 0.000228004336333 |

The modular response decreases strongly with distance.

The local MI response shows the same broad localization.

The global OR-curvature mean does not exhibit the same monotonic spatial structure.

## 11. Strength dependence

Mean absolute responses by perturbation magnitude:

| \(|\epsilon|\) | mean \(|\Delta M_K|\) | mean \(|\Delta M_W|\) | mean \(|\Delta\kappa_{
m edge}|\) |
|---:|---:|---:|---:|
| 0.05 | 0.000913611778959 | 0.0000272174482140 | 0.0000238751587475 |
| 0.10 | 0.00357822556563 | 0.000108639743360 | 0.0000957563213939 |
| 0.20 | 0.0133442127342 | 0.000430940322086 | 0.000387039314059 |

All three response magnitudes increase strongly with perturbation strength.

## 12. Reflection diagnostic

Compare the response generated at \(q\) with the reflected perturbation at \(9-q\).

Maximum reflection differences:

\[
\max_q
|
\Delta M_K(q)-\Delta M_K(9-q)
|
=
0.05848793366656263,
\]

\[
\max_q
|
\Delta M_W(q)-\Delta M_W(9-q)
|
=
0.001993514435335818,
\]

but

\[
\max_q
|
\Delta\kappa_{
m edge}(q)
-
\Delta\kappa_{
m edge}(9-q)
|
=
9.71445146547012	imes10^{-16}.
\]

Thus the global curvature scalar is reflection symmetric to numerical precision, whereas the local modular and MI observables are not.

This is consistent with their definitions:

- \(K_A\) is explicitly tied to the fixed left-edge subsystem \(A=(0,1,2,3)\);
- \(M_W\) is the MI average inside the same \(A\);
- \(\Delta\kappa_{
m edge}\) is an average over the entire frozen graph support.

The global curvature scalar therefore cannot distinguish a perturbation close to \(A\) from its mirror image far from \(A\).

## 13. Scientific interpretation

The strongest supported microscopic statement from this experiment is

\[
oxed{
|\Psi
angle

ightarrow

ho_A

ightarrow
K_A
\leftrightarrow
W_A.
}
\]

The controlled perturbation test shows that the state-derived modular sector and the local mutual-information geometry respond coherently.

The experiment does **not** establish

\[
K_A
\leftrightarrow
ar\kappa_{
m global}^{OR}.
\]

The negative H2 result should not be generalized to all local notions of emergent curvature.

The post-unblind diagnostics identify a specific structural mismatch:

\[
oxed{
	ext{local modular observable}
\quad	ext{vs}\quad
	ext{global reflection-symmetric curvature scalar}.
}
\]

This mismatch provides the scientific motivation for a new, separately preregistered experiment in which the curvature observable is localized relative to the same subsystem \(A\).

## 14. Next hypothesis

The next experiment should test a new quantity defined before any new data are generated:

\[
\kappa_A^{OR},
\]

with the curvature localization tied explicitly to the same fixed subsystem \(A\).

The new test must use fresh, previously unobserved perturbation conditions and a new preregistration.

The current experiment is complete and must not be retrospectively modified.

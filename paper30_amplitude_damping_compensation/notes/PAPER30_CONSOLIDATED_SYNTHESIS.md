# Paper 30 — Consolidated Scientific Synthesis

## Amplitude-damping compensation, causal exposure, and local causal-geometry landscape

### 1. Scope and status

This document consolidates the full Paper 30 analysis chain.

It distinguishes explicitly between:

1. frozen preregistered results;
2. post-unblinding implementation amendments;
3. failed preregistered controls;
4. post-run numerical diagnostics;
5. descriptive exploratory structure;
6. scientific statements that remain defensible.

No threshold, frozen result, classification, or preregistered decision rule is modified by this synthesis.

The strict preregistered reflection control for the final causal-compensation landscape remains failed.

Accordingly, the final landscape is not promoted to a fully preregistered positive result.

Post-run structural observations are reported descriptively.

---

## 2. Original Paper 30 amplitude-damping experiment

### 2.1 Protocol

The original Paper 30 experiment applied local amplitude damping to the finite-depth TFIM state.

Primary system:

- N = 12;
- subsystem A = (0,1,2,3);
- circuit depth P = 6.

Replication systems:

- N = 8;
- N = 10.

The damping parameter was:

gamma = 1 - exp(-s),

with:

s in [0,3]

sampled at step:

0.1.

Five frozen local-linear detrending bandwidth fractions were used:

0.05,
0.10,
0.20,
0.30,
0.40.

For each bandwidth, standardized residuals of the information-sector observable and modular-sector observable were compared.

The historical frozen support conditions were:

P1:
residual Pearson correlation < 0;

P2:
0.8 <= alpha_optimal <= 1.2;

P3:
R_C <= 0.20.

The cell was classified PAPER30_SUPPORTED only if every frozen bandwidth satisfied the support conditions.

### 2.2 Original result

For the primary N=12 configuration, all five bandwidths were supported.

Representative values were:

bandwidth 0.05:
rho = -0.9925486404537318,
alpha = 0.9925486404537315,
R_C = 0.01484719633344938;

bandwidth 0.10:
rho = -0.9814036850954925,
alpha = 0.981403685095492,
R_C = 0.03684680688098779;

bandwidth 0.20:
rho = -0.9804103192479913,
alpha = 0.9804103192479905,
R_C = 0.038795605912052035;

bandwidth 0.30:
rho = -0.986933312639062,
alpha = 0.9869333126390623,
R_C = 0.02596263640328707;

bandwidth 0.40:
rho = -0.9902599419001374,
alpha = 0.9902599419001378,
R_C = 0.01938524746793603.

The original Paper 30 classification was therefore:

PAPER30_SUPPORTED.

### 2.3 Limitation discovered after the original run

The N=10 and N=12 local observables were found to agree to approximately machine precision.

This raised a causal-boundary concern:

the original positive result might be a consequence of the finite circuit depth preventing the local subsystem from sensing the additional degrees of freedom.

The original positive result therefore required a causal-size robustness test.

---

## 3. Causal-size extension

### 3.1 Frozen configurations

The causal-size extension compared N=10 and N=12 under four configurations:

C0:
P=6, edge subsystem;

C1:
P=6, centered subsystem;

C2:
P=8, edge subsystem;

C3:
P=8, centered subsystem.

The primary size-sensitivity observable was:

D_A^(10,12)
=
(1/2) ||rho_A^(10) - rho_A^(12)||_1.

Frozen interpretation:

- equivalent if D_A <= 1e-10;
- size-sensitive if D_A > 1e-10.

### 3.2 Results

C0:

D_A_max
approximately
1.11998e-15;

size-equivalent;

compensation supported.

C1:

D_A_max
approximately
0.08348386;

size-sensitive;

compensation not supported.

C2:

D_A_max
approximately
0.0001680998;

size-sensitive;

compensation not supported.

C3:

D_A_max
approximately
0.40511198;

size-sensitive;

compensation not supported.

### 3.3 Frozen causal-size classification

The frozen extension classification was:

COMPENSATION_NOT_ROBUST_AFTER_CAUSAL_SIZE_EXPOSURE.

The original Paper 30 result therefore remains valid only for the tested finite-depth edge geometry.

It does not establish universal compensation.

An important nuance is that C2 is size-sensitive but still loses compensation.

Size sensitivity itself is therefore not a sufficient explanation of the compensation failure.

---

## 4. Causal-compensation landscape preregistration

To map the dependence on local causal geometry, a larger frozen landscape was preregistered.

Grid:

N in {8,10,12,14,16};

P in {1,...,12};

all contiguous four-qubit subsystems:

A = (a,a+1,a+2,a+3),

with:

a = 0,...,N-4.

Total preregistered cells:

540.

For each cell the following geometric descriptors were recorded:

d_L = a;

d_R = N-(a+4);

d_b = min(d_L,d_R);

finite-depth causal footprint:

[max(0,a-P), min(N-1,a+3+P)];

L_causal;

f_causal.

These are operational finite-depth circuit descriptors.

They are not identified with relativistic causal cones.

---

## 5. Attempt 001: rank-deficient modular sector

The first official landscape unblinding passed all C0-C3 recovery anchors.

It then failed on the first grid cell:

N=8,
P=1,
A=(0,1,2,3),
gamma=0.

The reduced state had structural rank deficiency.

Representative spectrum:

lambda_min
approximately
-3.325e-17;

lambda_max
approximately
0.740469;

effective rank:

4 / 16.

The modular Hamiltonian:

K_A = -log(rho_A)

was therefore not defined on the full Hilbert space under the frozen positivity requirement.

This was not a compression implementation bug.

The compressed and original full-global routes reproduced the same reduced state and spectrum.

No eigenvalue clipping, floor, pseudolog, or modified positivity tolerance was introduced.

Attempt 001 was therefore archived as an aborted official attempt.

---

## 6. Post-unblinding rank-deficiency amendment

A frozen amendment was introduced only to classify structurally undefined modular sectors explicitly.

The amendment allowed two specific nonfatal statuses:

MODULAR_SECTOR_UNDEFINED_RANK_DEFICIENT

and

MODULAR_SECTOR_UNDEFINED_STD_ZERO.

All 540 preregistered cells remained in the landscape.

Undefined modular cells were not treated as ordinary compensation failures.

All other unexpected RuntimeErrors remained fatal.

No support threshold or modular definition was altered.

---

## 7. Attempt 002: complete landscape

Attempt 002 completed the full preregistered grid.

Counts:

total cells:
540;

valid compensation-analysis cells:
485;

modular-sector undefined cells:
55;

valid fraction:
485 / 540
=
0.8981481481481481;

PAPER30_SUPPORTED valid cells:
68;

supported / all:
68 / 540
=
0.1259259259259259;

supported / valid:
68 / 485
=
0.1402061855670103.

### 7.1 Modular domain-of-definition structure

At P=1:

all subsystem positions are modular-sector undefined at gamma=0.

At P=2:

only two edge positions remain rank-deficient.

At P>=3:

all tested subsystem positions are valid for the frozen modular analysis.

Thus modular definability itself has a finite-depth activation structure.

However, modular definability is not sufficient for compensation.

For example:

P=3 contains valid cells but no supported causal classes.

---

## 8. Preregistered reflection control

The frozen reflection control compared mirrored subsystem configurations at tolerance:

1e-10.

Observed:

reflection checks:
300;

passes:
242;

failures:
58;

maximum absolute discrepancy:
approximately
1.3766901396650155e-08.

Therefore:

reflection_control_passed = False.

The preregistered interpretation gate was not satisfied.

This fact is retained unchanged.

---

## 9. Post-run reflection diagnostic

All 58 failed reflection pairs had:

identical Paper 30 support classification.

No self-mirror configuration failed.

Two numerical amplification mechanisms were identified.

### 9.1 Spectral-log amplification

Five extreme failures occurred for P=3 edge classes.

At gamma=0 the reduced state had condition number of order:

1.44e9.

Machine-scale perturbations in rho_A were amplified by:

K_A = -log(rho_A),

producing discrepancies of order:

1e-08

in downstream modular observables.

### 9.2 Detrending and z-score amplification

The remaining 53 failed pairs were not spectrally extreme.

Their raw information-sector discrepancies remained approximately:

1e-15.

After normalization, detrending, and z-scoring, these differences were amplified to approximately:

1e-09 to 1e-08.

The amplification magnitude was consistent with inverse residual standard deviations.

### 9.3 Safe interpretation

The strict frozen reflection control failed numerically.

The post-run audit found no evidence that the failures represent genuine physical reflection-symmetry breaking.

However, the frozen reflection criterion is not retroactively relaxed.

---

## 10. Causal-class collapse

For each valid cell define finite-depth left and right reaches:

ell = min(d_L,P);

r = min(d_R,P).

After quotienting left-right reflection define:

C
=
(P, min(ell,r), max(ell,r)).

The 485 valid finite-size cells collapse onto:

180 distinct local causal-geometry classes.

Observed:

pure supported classes:
27;

pure unsupported classes:
153;

mixed support classes:
0.

Therefore Paper 30 support is exactly constant inside every tested causal class.

All 68 supported finite-size cells belong to the 27 supported causal classes.

### 10.1 Bandwidth-resolved collapse

Five bandwidths across 180 causal classes produce:

900 causal-class x bandwidth groups.

Observed:

mixed bandwidth-support groups:
0.

Thus every frozen bandwidth-support decision is also constant inside each causal class.

### 10.2 Continuous metric collapse

Maximum within-class spans were of order:

few x 1e-09.

Examples:

pearson_residuals:
6.7423042149350465e-09;

alpha_optimal:
6.742303881868139e-09;

R_C:
6.034392430187552e-09;

rho_weak:
4.751270821223841e-09.

Therefore:

binary support classification collapses exactly;

continuous observables collapse numerically to approximately 1e-09.

---

## 11. Reduced geometric descriptors are insufficient

The complete descriptor:

(P, ell_near, ell_far)

is required to obtain zero mixed support classes.

Examples of weaker descriptors:

P alone:

11 groups,
8 mixed;

P + ell_near:

66 groups,
15 mixed,
70 percent of classes in pure groups;

P + ell_far:

66 groups,
20 mixed;

P + L_causal:

65 groups,
16 mixed;

P + causal asymmetry:

66 groups,
17 mixed.

Thus no tested simple lower-dimensional geometric descriptor determines Paper 30 support.

Total causal-footprint size alone is insufficient.

Global causal fraction alone is also insufficient.

The local left-right causal geometry matters.

---

## 12. Non-monotonic depth structure

The 27 supported classes do not form a single connected region in unit-Manhattan causal-coordinate space.

Observed:

connected components:
19;

isolated supported classes:
12;

maximum supported-graph degree:
2.

Several fixed causal geometries repeatedly enter and leave support as P increases.

Examples:

(ell_near,ell_far) = (0,4)

supported at:

P = 4,6,9;

and:

(ell_near,ell_far) = (1,3)

supported at:

P = 4,7,8,11.

Thus increasing circuit depth does not produce a monotonic transition toward or away from compensation.

---

## 13. Scalar reduction of the historical support criteria

After residual z-scoring:

alpha_optimal = -rho

and:

R_C = 1-rho^2.

The historical frozen conditions were:

rho < 0;

0.8 <= alpha_optimal <= 1.2;

R_C <= 0.20.

On the negative Pearson branch these are algebraically equivalent to:

rho <= -sqrt(0.8).

Numerically:

-sqrt(0.8)
=
-0.8944271909999159.

This is not a new threshold.

It is an algebraic reduction of the already frozen criteria.

### 13.1 Complete numerical verification

Valid bandwidth rows:

2425.

Maximum deviations from the algebraic identities:

max |alpha + rho|
=
1.3322676295501879e-15;

max |R_C - (1-rho^2)|
=
1.9567680809018384e-15.

Bandwidth-support mismatches against the scalar threshold:

0 / 2425.

### 13.2 Cell-level scalar field

For each valid cell:

rho_weak
=
maximum residual Pearson correlation across the five bandwidths.

Observed:

rho_weak
=
max_bandwidth_rho

exactly in the frozen output.

Cell-support mismatches against:

rho_weak <= -sqrt(0.8)

were:

0 / 485.

Causal-class mismatches:

0 / 180.

Thus the 27 supported causal classes are exactly the sub-threshold region of the scalar field:

rho_weak(P,ell_near,ell_far).

---

## 14. Robustness of the scalar support boundary

The closest supported causal class lies above the required support margin by:

0.0030594363237904654.

The closest unsupported causal class lies from the threshold by:

0.00044193137920456316.

The maximum observed within-class numerical spread of rho_weak is:

4.751270821223841e-09.

Therefore the nearest unsupported class is approximately:

9.3e4

times farther from the threshold than the largest observed within-class numerical spread.

The known reflection numerical discrepancies therefore do not threaten any frozen support classification.

---

## 15. Bottleneck bandwidth

For each causal class define the bottleneck bandwidth as the bandwidth that realizes:

rho_weak.

Across all 180 causal classes:

bandwidth 0.05:
54 classes;

bandwidth 0.10:
18 classes;

bandwidth 0.20:
17 classes;

bandwidth 0.30:
6 classes;

bandwidth 0.40:
85 classes.

Among the 27 supported classes:

0.05:
1;

0.10:
6;

0.20:
11;

0.30:
2;

0.40:
7.

The bottleneck identity is perfectly stable across finite-size repetitions inside each causal class:

classes_with_mixed_cell_bottleneck = 0.

### 15.1 Descriptive enrichment

Conditional support fractions are:

0.05:
1 / 54
=
0.018519;

0.10:
6 / 18
=
0.333333;

0.20:
11 / 17
=
0.647059;

0.30:
2 / 6
=
0.333333;

0.40:
7 / 85
=
0.082353.

The global causal-class support fraction is:

27 / 180
=
0.15.

Bandwidth 0.20 is therefore descriptively enriched among supported classes.

This is not used as a new support criterion.

---

## 16. Multiscale depth transitions

For fixed causal geometry there are:

135 consecutive-depth transitions.

Observed joint counts:

support unchanged,
bottleneck unchanged:
53;

support unchanged,
bottleneck changed:
46;

support changed,
bottleneck unchanged:
3;

support changed,
bottleneck changed:
33.

Thus:

36 transitions change support status.

Of these:

33 / 36

also change bottleneck bandwidth.

However:

46 bottleneck changes occur without a support transition.

Therefore bottleneck change is neither necessary nor sufficient for support change.

The scalar rho_weak remains the direct support variable.

### 16.1 Descriptive association

Observed:

P(support change | bottleneck change)
=
0.417722;

P(support change | no bottleneck change)
=
0.053571.

Descriptive risk ratio:

approximately 7.797.

Descriptive odds ratio:

approximately 12.674.

These are not inferential statistics because the transitions are not independent samples.

### 16.2 Directional scale reorganization

Using the explicitly post-run descriptive grouping:

intermediate bandwidths:
{0.10,0.20,0.30};

extreme bandwidths:
{0.05,0.40};

there are:

18 support entries;

18 support exits.

For support entries:

13 / 18 begin with an extreme bottleneck;

13 / 18 end with an intermediate bottleneck.

For support exits:

13 / 18 begin with an intermediate bottleneck;

16 / 18 end with an extreme bottleneck.

Thus support entry is descriptively associated with migration toward intermediate limiting scales, while support exit is descriptively associated with migration toward the extreme tested scales.

This pattern is post-run and is not claimed as a universal law.

---

## 17. What Paper 30 establishes

Within the tested finite-depth TFIM amplitude-damping protocol, the analysis supports the following statements.

### 17.1 Direct frozen outcomes

1. The original finite-depth edge configuration satisfies the frozen Paper 30 support criteria across all five preregistered bandwidths.

2. The causal-size extension has the frozen classification:

   COMPENSATION_NOT_ROBUST_AFTER_CAUSAL_SIZE_EXPOSURE.

3. The modular observable has a nontrivial domain-of-definition structure at shallow circuit depth, including structurally rank-deficient cells that remain explicitly undefined rather than being reclassified as failures.

4. Attempt 002 contains:

   - 540 total preregistered cells;
   - 485 valid compensation-analysis cells;
   - 55 modular-sector undefined cells;
   - 68 PAPER30_SUPPORTED valid cells.

5. The strict preregistered reflection control fails:

   - 242 / 300 checks pass;
   - 58 / 300 checks fail;
   - reflection_control_passed = False.

These frozen outcomes are retained without reinterpretation or threshold modification.

### 17.2 Post-run structural conclusions derived from frozen outputs

The following statements are post-run analyses of the frozen results rather than separately preregistered outcome criteria:

1. The 485 valid finite-size cells collapse onto 180 reflection-quotiented local causal-geometry classes with:

   - 27 pure supported classes;
   - 153 pure unsupported classes;
   - 0 mixed-support classes.

2. Every frozen bandwidth-support decision is also constant within these causal classes.

3. After residual z-scoring, the historical three-condition support rule reduces algebraically to:

   rho_weak <= -sqrt(0.8).

4. This scalar reduction has:

   - 0 / 2425 bandwidth-level mismatches;
   - 0 / 485 cell-level mismatches;
   - 0 / 180 causal-class mismatches.

5. The 27 supported causal classes are exactly the sub-threshold region of the post-run scalar field:

   rho_weak(P, ell_near, ell_far).

6. This scalar field is strongly non-monotonic with circuit depth over the tested grid.

7. The observed causal-class collapse is consistent with finite-depth circuit locality and shows that global system size alone does not determine the frozen support classification.

These conclusions are strongly supported by the frozen numerical outputs, but their formulation is post-run.

### 17.3 Post-run exploratory multiscale structure

The following observations are additionally descriptive and exploratory:

- the supported causal classes form a sparse disconnected set;
- bottleneck bandwidth varies systematically over causal geometry;
- intermediate bottleneck scales, especially 0.20, are enriched among supported classes;
- support transitions are frequently accompanied by a change in bottleneck scale;
- entry and exit transitions show an asymmetric extreme/intermediate bottleneck flow.

These observations require independent preregistered confirmation before being treated as new laws or predictive rules.

---

## 18. What Paper 30 does not establish

Paper 30 does not establish:

- a universal conservation law between information and modular sectors;
- literal transfer of information into modular energy;
- universal compensation under open-system dynamics;
- independence from finite-depth circuit geometry;
- independence from boundary conditions;
- a relativistic causal structure;
- emergent spacetime from this experiment alone;
- a thermodynamic arrow of time;
- genuine physical reflection-symmetry breaking;
- a universal preferred detrending bandwidth.

The failed strict reflection control must remain visible in any final scientific presentation.

---

## 19. Recommended final interpretation

The safest consolidated interpretation is:

The original amplitude-damping compensation signal is numerically strong within its tested finite-depth local geometry but is not universal.

Once causal exposure is varied systematically, compensation occupies a sparse and strongly non-monotonic subset of local causal geometries.

Post-run analysis shows that the finite-size landscape collapses onto reflection-quotiented local causal classes, consistent with finite-depth local dependence and showing that global system size alone does not determine the observed support classification.

After z-scoring, the historical support conditions reduce to a single scalar field rho_weak with a frozen threshold at -sqrt(0.8).

The resulting support structure is therefore best viewed as a structured threshold landscape over local finite-depth causal geometry.

A secondary multiscale organization is visible in which detrending bandwidth realizes rho_weak, but this bottleneck identity does not replace the scalar correlation field as the primary descriptor.

The strict preregistered reflection control failed numerically, so the causal-landscape interpretation remains formally descriptive rather than a fully preregistered positive result.

---

## 20. Recommended manuscript-level summary

A concise manuscript-safe summary is:

"Under local amplitude damping, the original finite-depth edge configuration exhibits a strong quasi-compensated anti-correlation between information-sector and modular-sector residuals, but this behavior is not robust once the subsystem is causally exposed to additional degrees of freedom. A preregistered 540-cell causal landscape yields 485 valid modular-sector cells and 68 supported finite-size cells. Post-run analysis of these frozen outputs collapses the valid cells onto 180 reflection-quotiented local causal geometries, of which 27 satisfy the frozen compensation criterion. After residual z-scoring, the historical three-condition support rule is algebraically equivalent to rho_weak <= -sqrt(0.8), so the supported region is exactly the sub-threshold portion of a continuous correlation field over local causal geometry. The field varies strongly and non-monotonically with circuit depth. Post-run diagnostics further reveal a multiscale organization of the limiting detrending bandwidth. However, the strict preregistered reflection control failed numerically, and these landscape-level structural conclusions are therefore retained as descriptive rather than promoted to an unqualified preregistered positive result."


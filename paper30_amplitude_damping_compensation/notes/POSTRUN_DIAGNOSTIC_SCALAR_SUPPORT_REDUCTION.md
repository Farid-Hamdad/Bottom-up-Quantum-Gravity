# Paper 30 — Causal Compensation Landscape

## Post-run diagnostic: scalar reduction of the frozen support criterion

### 1. Status

This is a post-run algebraic and numerical diagnostic of the frozen Attempt 002 landscape.

Relevant frozen commits:

- 23acc83 Freeze Paper 30 landscape attempt 002 results
- 0198946 Document Paper 30 landscape reflection diagnostics
- 61acee0 Document Paper 30 causal-class collapse

No preregistered threshold, observable, bandwidth, support criterion, landscape cell, or frozen result is modified here.

The strict preregistered reflection control remains:

reflection_control_passed = False

### 2. Frozen bandwidth-level criteria

For each frozen bandwidth, Paper 30 support requires simultaneously:

P1: pearson_residuals < 0

P2: 0.8 <= alpha_optimal <= 1.2

P3: R_C <= 0.20

After residual z-scoring, the frozen definitions imply:

alpha_optimal = -rho

and

R_C = 1 - rho^2

where rho is the residual Pearson correlation.

### 3. Algebraic reduction

Because a Pearson correlation satisfies

-1 <= rho <= 1,

P1 restricts the relevant branch to

rho < 0.

On that branch:

alpha_optimal = -rho

lies between 0 and 1.

Therefore P2 reduces to:

-rho >= 0.8

or equivalently:

rho <= -0.8.

P3 gives:

1 - rho^2 <= 0.20,

hence:

rho^2 >= 0.8.

Together with rho < 0:

rho <= -sqrt(0.8).

Numerically:

-sqrt(0.8) = -0.8944271909999159.

Because

-sqrt(0.8) < -0.8,

P3 is stricter than the lower-bound part of P2 once P1 is satisfied.

Thus the complete frozen bandwidth-level support criterion is algebraically equivalent to the single scalar condition:

bandwidth_supported
if and only if
rho <= -sqrt(0.8).

This is a reduction of the already frozen definitions, not a new threshold.

### 4. Numerical verification over all valid bandwidth rows

The frozen Attempt 002 landscape contains:

2425 valid bandwidth rows.

Across all of them:

max_abs(alpha_optimal + rho)
= 1.3322676295501879e-15

and

max_abs(R_C - (1-rho^2))
= 1.9567680809018384e-15.

Observed mismatches between the original three-criterion bandwidth classification and the scalar threshold:

0 / 2425.

Therefore the algebraic reduction is numerically verified across the complete valid bandwidth dataset.

### 5. Cell-level reduction

For every valid cell, define:

rho_weak = maximum residual Pearson correlation over the five frozen bandwidths.

The frozen cell-level PAPER30_SUPPORTED classification requires all five bandwidths to be supported.

Therefore:

PAPER30_SUPPORTED
if and only if
rho_weak <= -sqrt(0.8).

The frozen data verify:

max_abs(
    rho_weak
    - max_bandwidth_pearson_residuals
)
= 0.

Observed cell-level mismatches against the scalar criterion:

0 / 485.

Thus the complete Paper 30 support decision for every valid landscape cell is encoded by one continuous scalar observable:

rho_weak.

### 6. Causal-class reduction

The 485 valid cells collapse onto:

180 reflection-quotiented local causal-geometry classes.

Observed mismatches between causal-class support and the scalar criterion:

0 / 180.

The maximum observed within-class spread of rho_weak is:

4.751270821223841e-09.

Thus the causal compensation landscape may be represented as a scalar field:

rho_weak(P, ell_near, ell_far)

over the discrete local causal-geometry space, together with the frozen support threshold:

rho_weak <= -sqrt(0.8).

The 27 supported causal classes are exactly the sub-threshold region of this scalar field.

### 7. Robustness of the support boundary

The supported class closest to the threshold has positive support margin:

0.0030594363237904654.

The unsupported class closest to the threshold has distance:

0.00044193137920456316.

The maximum observed within-class numerical spread of rho_weak is:

4.751270821223841e-09.

Therefore the nearest unsupported class lies approximately:

9.30e4

times farther from the support boundary than the largest observed within-class rho_weak numerical spread.

The closest supported class lies approximately:

6.44e5

times farther from the boundary than that same spread.

Accordingly, the known reflection/numerical discrepancies are far too small to alter any causal-class support decision in the frozen dataset.

### 8. Non-monotonic depth crossings

The sparse topology of the 27 supported causal classes can therefore be interpreted as threshold crossings of the continuous field rho_weak rather than as an additional independent binary structure.

Examples include the fixed causal geometry:

(ell_near, ell_far) = (0,4),

which is supported at:

P = 4, 6, 9,

but unsupported at the other tested depths.

Its rho_weak values alternate strongly across the fixed threshold.

Likewise:

(ell_near, ell_far) = (1,3)

is supported at:

P = 4, 7, 8, 11,

and unsupported at the other tested depths.

Thus increasing circuit depth does not move the system monotonically toward or away from compensation.

Depth changes the internal finite-depth dynamics, and rho_weak can cross the same frozen support threshold multiple times.

### 9. Scientific interpretation

The historical Paper 30 support criteria should not be counted as three statistically independent signatures after z-scoring.

They are algebraically linked.

The independent descriptive object is the residual correlation field, in particular its least favorable bandwidth value:

rho_weak.

The causal landscape can therefore be summarized as:

1. a domain-of-definition structure for the modular observable;
2. a scalar field rho_weak over valid local causal geometries;
3. a frozen threshold at -sqrt(0.8);
4. a sparse, non-monotonic sub-threshold region containing 27 of 180 valid causal classes.

The nontrivial result is not the algebraic threshold reduction itself.

The nontrivial result is the structured and non-monotonic dependence of rho_weak on circuit depth and local causal geometry.

### 10. Safe manuscript statement

A concise statement suitable for later manuscript use is:

"After residual z-scoring, the three historical Paper 30 bandwidth-support conditions are algebraically equivalent to a single threshold on the residual Pearson correlation, rho <= -sqrt(0.8). At cell level this becomes rho_weak <= -sqrt(0.8), where rho_weak is the least negative correlation across the five frozen bandwidths. This equivalence holds for all 2425 valid bandwidth rows, all 485 valid cells, and all 180 causal classes. The 27 supported causal classes are therefore exactly the sub-threshold region of a continuous rho_weak field on local causal-geometry space. Their non-monotonic recurrence with circuit depth reflects repeated crossings of the same frozen threshold rather than independent binary phases."

No new threshold is introduced and no frozen result is reclassified.

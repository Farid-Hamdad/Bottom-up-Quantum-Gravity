# Paper 30 — Causal Compensation Landscape

## Post-run diagnostic: collapse onto local causal-geometry classes

### 1. Status

This note is a post-run diagnostic of the frozen Attempt 002 landscape.

Frozen results commit:

23acc83 Freeze Paper 30 landscape attempt 002 results

Frozen reflection-diagnostic commit:

0198946 Document Paper 30 landscape reflection diagnostics

No preregistered threshold, observable, support criterion, landscape cell, or frozen result is modified here.

The strict preregistered reflection control remains:

reflection_control_passed = False

### 2. Landscape counts

The frozen landscape contains:

- 540 preregistered cells;
- 485 valid compensation-analysis cells;
- 55 modular-sector undefined cells;
- 68 Paper 30 supported valid cells.

Among valid cells, the supported fraction is:

68 / 485 = 0.1402061855670103.

### 3. Local causal-geometry signature

For a contiguous four-qubit subsystem with distances to the physical chain boundaries

d_L and d_R,

define the finite-depth causal reaches

ell = min(d_L, P)

and

r = min(d_R, P).

After quotienting left-right reflection, define the local causal-geometry signature

C = (P, min(ell,r), max(ell,r)).

This descriptor retains:

- circuit depth P;
- accessible distance to the nearer boundary;
- accessible distance to the farther boundary;

while removing absolute system size and left-right orientation.

### 4. Exact collapse of Paper 30 support classification

The 485 valid cells collapse onto:

180 distinct local causal-geometry classes.

Observed class counts:

- pure supported classes: 27;
- pure unsupported classes: 153;
- mixed support classes: 0.

Therefore every valid cell belonging to the same frozen causal-geometry class has the same PAPER30_SUPPORTED classification.

All 68 supported cells are contained in the 27 pure supported classes.

There are no supported cells outside those classes.

### 5. Bandwidth-resolved classification collapse

The five frozen bandwidth fractions generate:

900 causal-class x bandwidth groups.

Observed:

mixed_bandwidth_support_groups = 0.

Thus bandwidth-level support is also exactly constant within every local causal-geometry class.

At cell level:

f_supported

has zero within-class span over all 180 classes.

The classification collapse therefore holds both for:

- the final PAPER30_SUPPORTED cell classification; and
- every individual frozen bandwidth-support decision.

### 6. Continuous observables

The continuous compensation metrics also approximately collapse within each causal class, but not identically at machine precision.

Maximum observed within-class spans are:

- pearson_residuals:
  6.7423042149350465e-09;

- alpha_optimal:
  6.742303881868139e-09;

- R_C:
  6.034392430187552e-09;

- rho_strong:
  6.7423042149350465e-09;

- rho_weak:
  4.751270821223841e-09;

- R_C_worst:
  2.3733656151492255e-09.

Spearman residual correlation has zero observed within-class span.

These approximately 1e-9 discrepancies are compatible in scale with the numerical amplification mechanisms identified independently in the frozen reflection-control audit.

Accordingly, this note does not claim exact equality of continuous metrics.

The robust statement is:

- continuous metrics collapse to approximately few-parts-in-1e9 numerical agreement;
- all frozen threshold-based support classifications collapse exactly.

### 7. Total causal size is insufficient

A weaker descriptor using only

(P, L_causal)

produces:

65 classes,

of which:

16 are mixed with respect to Paper 30 support.

Therefore total causal-footprint size alone is insufficient to determine the compensation classification.

The left-right distribution of causal reach relative to physical boundaries carries additional information.

### 8. Global causal fraction is not the determining variable

Within a fixed local causal-geometry class:

L_causal has zero observed span,

whereas

f_causal = L_causal / N

can vary substantially as N changes.

The maximum observed within-class span of f_causal is:

0.5.

Nevertheless, the Paper 30 support classification remains exactly invariant within those classes.

This shows that the observed classification is not organized simply by the fraction of the global chain covered by the finite-depth causal footprint.

### 9. Depth-dependent supported causal classes

Support is highly nonuniform in circuit depth.

Number of supported causal classes by depth:

- P=2: 0;
- P=3: 0;
- P=4: 3;
- P=5: 0;
- P=6: 3;
- P=7: 2;
- P=8: 1;
- P=9: 8;
- P=10: 2;
- P=11: 6;
- P=12: 2.

This is consistent with the previously observed support landscape:

- a strong outer-layer band at P=4;
- predominantly edge support at P=6;
- migration toward interior causal geometries at larger depths, especially P=9 to P=12.

The depth dependence is therefore structured and non-monotonic.

### 10. Interpretation

The causal-class collapse should not by itself be interpreted as an independent dynamical law.

For a finite-depth local circuit acting on the same product-state preparation, observables inside identical local causal neighborhoods are expected to be insensitive to degrees of freedom outside the finite-depth causal domain, up to numerical effects and reflection.

The collapse is therefore primarily:

1. a strong structural validation that the frozen numerical landscape respects finite-depth circuit locality; and
2. a dimensional reduction of the Paper 30 landscape from 485 valid finite-size cells to 180 physically distinct local causal geometries.

The nontrivial empirical content is not that equivalent causal neighborhoods agree.

The nontrivial content is which causal geometries support compensation.

In the tested landscape:

27 of 180 valid local causal-geometry classes support Paper 30 compensation.

Thus the compensation regime is not universal and is not determined by system size alone, distance to a single boundary alone, total causal-footprint size alone, or global causal fraction alone.

Instead, it is organized by a structured combination of circuit depth and local left/right causal geometry.

### 11. Safe scientific statement

A concise statement suitable for later manuscript use is:

"Within the tested finite-depth landscape, the Paper 30 compensation classification collapses exactly onto reflection-quotiented local causal-geometry classes. Continuous compensation metrics agree within each class to approximately 1e-9, while all frozen bandwidth and cell-level support decisions are identical. The total causal-footprint size alone is insufficient: the left-right geometry of the finite-depth causal neighborhood is required. This collapse is consistent with finite-depth circuit locality; the substantive result is the sparse, structured subset of causal geometries that supports compensation."

No claim of relativistic causality, emergent spacetime, conservation law, or universal compensation is made by this diagnostic.

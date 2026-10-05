# Paper 30 — Causal Compensation Landscape

## Post-run diagnostic: bottleneck bandwidth and multiscale support transitions

### 1. Status

This note is a post-run descriptive diagnostic of the frozen Attempt 002 landscape.

Relevant frozen commits:

- 23acc83 Freeze Paper 30 landscape attempt 002 results
- 0198946 Document Paper 30 landscape reflection diagnostics
- 61acee0 Document Paper 30 causal-class collapse
- 2618561 Document Paper 30 scalar support reduction

No preregistered threshold, bandwidth, support criterion, observable, landscape cell, or frozen result is modified here.

The strict preregistered reflection control remains:

reflection_control_passed = False

### 2. Definition

For each valid cell or causal class, define the bottleneck bandwidth as the frozen bandwidth whose residual Pearson correlation is maximal:

rho_weak = max_b rho_b.

The bottleneck bandwidth is therefore the bandwidth that realizes rho_weak and hence provides the least favorable residual correlation for the frozen all-bandwidth support requirement.

The support decision itself remains:

rho_weak <= -sqrt(0.8).

The bottleneck identity is not an additional support criterion.

### 3. Bottleneck stability inside causal classes

Across the 180 valid local causal-geometry classes:

classes_with_mixed_cell_bottleneck = 0.

Thus every finite-size repetition belonging to the same causal class identifies the same bottleneck bandwidth.

The bottleneck identity therefore collapses onto the same local causal-geometry classes as the Paper 30 support classification.

### 4. Bottleneck distribution over causal classes

Across all 180 causal classes:

- bandwidth 0.05: 54 classes;
- bandwidth 0.10: 18 classes;
- bandwidth 0.20: 17 classes;
- bandwidth 0.30: 6 classes;
- bandwidth 0.40: 85 classes.

Among the 27 supported causal classes:

- bandwidth 0.05: 1 class;
- bandwidth 0.10: 6 classes;
- bandwidth 0.20: 11 classes;
- bandwidth 0.30: 2 classes;
- bandwidth 0.40: 7 classes.

Among the 153 unsupported causal classes:

- bandwidth 0.05: 53 classes;
- bandwidth 0.10: 12 classes;
- bandwidth 0.20: 6 classes;
- bandwidth 0.30: 4 classes;
- bandwidth 0.40: 78 classes.

The support fraction conditional on bottleneck bandwidth is:

- 0.05: 1 / 54 = 0.018519;
- 0.10: 6 / 18 = 0.333333;
- 0.20: 11 / 17 = 0.647059;
- 0.30: 2 / 6 = 0.333333;
- 0.40: 7 / 85 = 0.082353.

The global causal-class support fraction is:

27 / 180 = 0.15.

Therefore bandwidth 0.20 is descriptively enriched among supported classes, whereas the two extreme bandwidths 0.05 and 0.40 are comparatively depleted.

This is a post-run descriptive association and is not introduced as a new classification rule.

### 5. Robustness of bottleneck identity

The gap between the largest and second-largest residual Pearson correlation was evaluated for every causal class.

Observed:

- minimum bottleneck gap:
  2.8327210223055843e-05;

- median bottleneck gap:
  1.4450308589428415e-01;

- maximum bottleneck gap:
  4.0712209024788704e-01.

No causal class has a bottleneck gap below 1e-06.

Only:

- 2 / 180 classes have gap <= 1e-04;
- 12 / 180 classes have gap <= 1e-03;
- 32 / 180 classes have gap <= 1e-02.

Thus the identity of the bottleneck is not determined by the approximately 1e-09 numerical variations identified in the reflection audit.

### 6. Consecutive-depth transitions

For fixed reflection-quotiented geometry (ell_near, ell_far), there are:

135 consecutive-depth transitions.

Their joint behavior is:

- support unchanged, bottleneck unchanged: 53;
- support unchanged, bottleneck changed: 46;
- support changed, bottleneck unchanged: 3;
- support changed, bottleneck changed: 33.

Thus:

36 / 135 transitions change Paper 30 support status.

Among those 36 support transitions:

33 / 36 = 0.916667

also change bottleneck bandwidth.

However:

46 bottleneck changes occur without any support transition.

Therefore a bottleneck change is strongly associated with a support transition, but is neither necessary nor sufficient for one.

The scalar rho_weak remains the variable that directly determines the frozen support boundary.

### 7. Descriptive association strength

Among transitions with a bottleneck change:

P(support change | bottleneck change)
= 33 / 79
= 0.417722.

Among transitions without a bottleneck change:

P(support change | no bottleneck change)
= 3 / 56
= 0.053571.

The corresponding descriptive ratios are approximately:

- risk ratio: 7.797;
- odds ratio: 12.674.

These values are descriptive only.

The transitions are not independent statistical samples, because multiple transitions belong to the same causal geometries and depth sequences.

No inferential p-value or causal interpretation is assigned.

### 8. Direction of scale reorganization

For descriptive purposes only, define:

intermediate bandwidths:
{0.10, 0.20, 0.30}

and extreme bandwidths:
{0.05, 0.40}.

This grouping is post-run and was not preregistered.

There are:

18 entries into Paper 30 support

and

18 exits from Paper 30 support.

For entries into support:

- 13 / 18 start from an extreme bottleneck;
- 5 / 18 start from an intermediate bottleneck;

while:

- 13 / 18 terminate at an intermediate bottleneck;
- 5 / 18 terminate at an extreme bottleneck.

The entry flow counts are:

- extreme -> intermediate: 9;
- extreme -> extreme: 4;
- intermediate -> intermediate: 4;
- intermediate -> extreme: 1.

For exits from support:

- 13 / 18 start from an intermediate bottleneck;
- 5 / 18 start from an extreme bottleneck;

while:

- 16 / 18 terminate at an extreme bottleneck;
- 2 / 18 terminate at an intermediate bottleneck.

The exit flow counts are:

- intermediate -> extreme: 12;
- extreme -> extreme: 4;
- extreme -> intermediate: 1;
- intermediate -> intermediate: 1.

Thus the observed transitions show a directional descriptive tendency:

entry into the supported regime is frequently accompanied by migration of the limiting scale from an extreme bandwidth toward an intermediate bandwidth,

whereas exit from the supported regime is frequently accompanied by migration toward an extreme bandwidth.

This pattern is not used as a classifier and is not claimed as a universal law.

### 9. Same-bottleneck threshold crossings

Three support transitions occur without any bottleneck change:

1. geometry (1,3), P=8 -> 9:
   bottleneck = 0.20,
   rho_weak = -0.975582511 -> -0.864696572,
   support True -> False;

2. geometry (1,3), P=11 -> 12:
   bottleneck = 0.40,
   rho_weak = -0.897486627 -> 0.119146647,
   support True -> False;

3. geometry (0,4), P=4 -> 5:
   bottleneck = 0.05,
   rho_weak = -0.964434726 -> 0.971747723,
   support True -> False.

These cases demonstrate directly that a support transition does not require a change of limiting bandwidth.

The value of rho_weak itself can cross the frozen threshold while the same bandwidth remains the bottleneck.

### 10. Changes in rho_weak

Median absolute changes in rho_weak across consecutive depths are:

- support change + bottleneck change:
  1.526920;

- support change + same bottleneck:
  1.016633;

- same support + bottleneck change:
  0.490250;

- same support + same bottleneck:
  0.127204.

Support transitions therefore coincide with substantially larger changes in the continuous rho_weak field than transitions that preserve support.

### 11. Interpretation

The bottleneck bandwidth should be interpreted as a secondary multiscale descriptor of the scalar rho_weak field.

The frozen classification remains completely determined by:

rho_weak <= -sqrt(0.8).

The bottleneck identifies which detrending scale realizes the least favorable correlation at each causal geometry.

The observed landscape is therefore multiscale in the following limited sense:

- different causal geometries and depths are limited by different frozen bandwidths;
- the limiting bandwidth is stable inside each causal class;
- support transitions are usually accompanied by a change in limiting scale;
- but limiting-scale changes also frequently occur without crossing the support boundary;
- and threshold crossings can occur while the limiting bandwidth remains unchanged.

Consequently, bottleneck identity does not replace rho_weak as the primary continuous descriptor.

### 12. Safe manuscript statement

A concise statement suitable for later manuscript use is:

"The bandwidth that realizes rho_weak is itself organized by local causal geometry and varies across the landscape. Support transitions are frequently accompanied by a change in this limiting detrending scale (33 of 36 consecutive-depth support transitions), although bottleneck changes are neither necessary nor sufficient for crossing the frozen support threshold. Descriptively, entries into the supported regime tend to shift the limiting scale toward intermediate bandwidths, whereas exits tend to shift it toward the extreme tested bandwidths. The support boundary nevertheless remains determined solely by rho_weak <= -sqrt(0.8)."

No new support criterion, threshold, or physical law is introduced by this post-run diagnostic.

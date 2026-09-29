# Paper 23 — Modular Cut-Density Audit v2

## Status

This directory contains a retrospective reproducibility and adversarial-control audit of the modular cut-density test introduced in Paper 23.

The historical Paper 23 manuscript is not modified by this audit.

The reference computation was executed twice with identical source files, parameters and random seeds. The generated JSON and CSV artifacts were bitwise identical across both runs.

## Scientific question

Paper 23 investigated whether a local entanglement quantity crossing the boundary of a subsystem,

\[
\rho_{\rm cut}(j) = \sum_{k \notin A} I(j:k),
\]

tracks the local modular-temperature profile

\[
T_{\rm mod}(j) = \frac{1}{\beta_j}.
\]

The original motivation was that, in the BuP construction, information crossing the subsystem boundary and modular thermodynamics might represent two descriptions of the same emergent structure.

The purpose of this audit is not to test BuP as a whole. It is to determine what is actually supported by the specific rho_cut versus 1/beta comparison in the Joshi/Kokail/Zoller trapped-ion data.

## Source data

Official archived data:

DOI: 10.5281/zenodo.8279583

Files used:

- RawData/dataRawParam_GS_Delta_1.mat
  SHA256: ee4737a8b111ebfebd4502b4bb540ccddc2f049d06d32d01d4009edea863b135

- RawData/dataRawParam_ES_Delta_1.mat
  SHA256: b2f33c879d14c1a4c516df86ac8d2ae458db27b6b433a8688a5fd548fd5c2c18

- AnalyzedData/Figure2/GroundStateBetas.mat
  SHA256: 74d359adc97d0a474171790020933f6beea5a011502229fabb6e352be7349a2a

- AnalyzedData/Figure2/ExcitedStateBetas.mat
  SHA256: 02abda667711fe971b37c639830baad6d0361ec49525ccfb3615213049fc80ea

- Readme.docx
  SHA256: ae7ee4c1ac14bf34a03f6277f398dfff81f3f9f9e354b6c6feb46b0fdb3a8510

## Tomographic measurement design

The official README defines tomography-history codes as:

1 = X
2 = -X
3 = Y
4 = -Y
5 = Z
6 = -Z

For the Delta = 1 GS and ES files:

- 51 ions;
- 243 measurement settings;
- 200 outcomes per setting;
- GS and ES use identical measurement histories;
- every ion is measured 81 times along X, 81 times along Y and 81 times along Z.

Single-ion axis sampling is therefore exactly uniform.

Pairwise sampling is not uniformly independent.

Of the 1275 possible ion pairs:

- 1040 have all nine joint axis combinations exactly 27 times each;
- 235 have incomplete joint-axis coverage;
- the 235 incomplete pairs are exactly those satisfying i = j modulo 5.

For these pairs, only XX, YY and ZZ occur, each 81 times.

This explains why an estimator using a universal factor of 9 for two-qubit Pauli inclusion probabilities is not valid for the complete measurement design.

The audit therefore uses the conditional estimator for the reconstructible pairwise Pauli coefficients and treats incompletely tomographed pairs as unavailable.

## Coverage control

Because 235 pairwise reduced states cannot be completely reconstructed, the historical cut density is more precisely an observed cut density over reconstructible links.

The audit therefore compares the historical sum

\[
\rho_{\rm cut,sum}(j)
\]

with a coverage-normalized quantity

\[
\rho_{\rm cut,mean}(j)
=
\rho_{\rm cut,sum}(j) / C_cut(j),
\]

where C_cut is the number of reconstructible links crossing the subsystem boundary.

For the five experimental GS profiles:

Historical cut sum:

- mean Pearson = 0.5880493209434857
- mean Spearman = 0.7774891774891776

Coverage-normalized cut mean:

- mean Pearson = 0.5861597170037138
- mean Spearman = 0.8008225108225109

The cut-density association therefore survives normalization by the number of reconstructible links.

The simple measurement-coverage pattern is not sufficient to explain the historical correlation.

## Cut specificity

After coverage normalization, the five-profile GS averages are:

### Global density

- Pearson = 0.049485933278633695
- Spearman = -0.08792207792207789

### Internal density

- Pearson = -0.17485865169644688
- Spearman = -0.2737229437229437

### Cut density

- Pearson = 0.5861597170037138
- Spearman = 0.8008225108225109

Thus the descriptive association with the modular profile is specific to information crossing the subsystem boundary. It is not reproduced by the global or internal pairwise-MI density.

This cut specificity survives the measurement-coverage correction.

## Edge-geometry control

Both the modular-temperature profile and the cut density can naturally acquire an edge-enhanced spatial profile.

A geometric edge predictor proportional to

\[
\frac{1}{j(L+1-j)}
\]

was therefore introduced as a confound control.

For GS experimental profiles with L >= 5, after residualizing against this edge predictor:

- mean partial Pearson = -0.029928251590637167
- mean partial Spearman = 0.17655949451036312

The original Pearson association is therefore not retained after removing the common edge geometry.

## Empirical distance-decay control

A stronger null was constructed directly from the experimental mutual-information matrix.

For each physical separation d, a reference MI scale I(d) was estimated from experimentally reconstructed pairs. A cut-density profile expected solely from:

- one-dimensional geometry,
- MI decay with physical separation,
- the actual tomography mask,

was then generated without using beta.

For GS profiles with L >= 5:

- mean partial Pearson after controlling the distance null = 0.0060309584346819045
- mean partial Spearman = 0.16343629806233373

No stable residual Pearson association remains.

## Distance-preserving permutation test

A permutation null was next constructed by shuffling experimental MI values only among pairs with the same physical separation.

This preserves the empirical distance dependence while destroying ion-specific structure.

No GS or ES experimental profile produced a one-sided empirical p value below 0.05 for either Pearson or Spearman correlation.

Therefore the observed rho_cut versus 1/beta correlations are not exceptional relative to a null that preserves the measured locality structure.

## Leave-A-out distance-matched null

The strictest control avoids using links touching the tested subsystem A when constructing the reference distribution.

For every observed link from A to its complement, surrogate MI values are sampled only from fully external pairs with the same physical separation.

### Ground state

Pearson empirical p values across L = 3, 5, 7, 9, 11:

0.5311, 0.7732, 0.8092, 0.9578, 0.6931

Spearman empirical p values:

0.3243, 0.3589, 0.1958, 0.4429, 0.3949

### Excited state

Pearson empirical p values:

0.3981, 0.3245, 0.6877, 0.9870, 0.5319

Spearman empirical p values:

0.7319, 0.6413, 0.8842, 0.9302, 0.8410

None of the ten experimental profiles exceeds the distance-matched leave-A-out null at p < 0.05.

## What survives the audit

The following empirical result survives:

> The pairwise mutual information crossing the subsystem boundary reproduces the spatial boundary structure of the modular-temperature profile substantially better than either the global or internal pairwise-MI density.

This result:

- is reproducible from the official archived data;
- is specific to the cut observable;
- survives correction for incomplete pairwise tomography;
- occurs in a physically meaningful boundary-localized quantity.

## What does not survive the adversarial controls

The stronger interpretation is not established:

> rho_cut contains modular information independent of ordinary spatial locality.

After explicitly controlling for boundary geometry and experimental MI decay with distance, no robust residual association is detected.

Distance-preserving and leave-A-out permutation nulls reproduce correlations of comparable or greater magnitude.

Therefore the historical rho_cut versus 1/beta correlation cannot, on these data alone, be interpreted as an independent validation of modular thermodynamics or of BuP.

## Revised interpretation

The result should be described as a boundary-structure correspondence:

\[
\rho_{\rm cut}
\quad \text{and} \quad
\frac{1}{\beta}
\]

select compatible spatial boundary structure, whereas rho_all and rho_in do not.

The present data do not distinguish this correspondence from an explanation based on locality plus the geometry of the subsystem boundary.

A stronger test requires a benchmark in which modular structure can vary independently of simple physical distance or boundary geometry.

This motivates a prospective benchmark on a programmable quantum device, with states, partitions, null hypotheses and decision criteria fixed before examining the final data.

## Claims not supported by this audit

This audit does not show that:

- BuP is false;
- modular thermodynamics is false;
- rho_cut is physically irrelevant;
- the Zoller/Joshi experiment validates BuP;
- the Zoller/Joshi experiment refutes BuP;
- the cut-density correlation contains an independently identified modular component.

It establishes the evidential scope of this particular benchmark.

## Reproducibility

Reference script:

run_audit.py

SHA256:

e1326f49b4fe71e57decaf79b15e3cbde106125538da9b968e9104e59f403e6b

Reference artifacts:

results/audit_summary.json

SHA256:

f9959014b9679d096b4d84546edf5db36d4a67fa65449687537fa9b97eec369d

results/profile_correlations.csv

SHA256:

b028a3fa591a6b68e0fd9e0517d288a32bc05a1fe3a39cf778370fb3f5cf1556

The reference audit was executed twice with 5000 permutations per permutation family and fixed seeds.

Both executions generated bitwise-identical JSON and CSV artifacts.

## Audit status

RETROSPECTIVE_AUDIT_REPRODUCED

The numerical audit is frozen at this stage.

Additional hypotheses or controls should be treated as a subsequent audit version rather than silently added to this reference analysis.

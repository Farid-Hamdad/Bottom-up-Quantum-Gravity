
# SPARC

This folder contains the SPARC phenomenology module of the **Bottom-Up Quantum Gravity** project.

Its purpose is to test whether the Bottom-Up framework, initially developed on finite quantum systems, can be extended toward **real galaxy rotation-curve data** through an **effective emergent dimension**.

---

## Scientific objective

The central idea explored here is that galaxy dynamics may be described, at least phenomenologically, through a variable effective dimension \( d(x) \), rather than through an explicit dark-matter halo.

The main tested model family is based on a dimension profile of the form

\[
d(x)=3-\Delta d \exp\!\left[-\frac{\left(\ln(x/x_0)\right)^2}{2\sigma^2}\right]
\]

with different choices for the variable \(x\), including:

- \(x=r\)
- \(x=r/R_d\)
- \(x=g_{\mathrm{bar}}\)

The goal is not to claim a final galactic theory at this stage, but to build a **transparent and falsifiable experimental chain**:
- historical claim re-evaluation,
- real-data baselines,
- subgroup analysis,
- robustness tests,
- out-of-sample prediction,
- empirical inversion of \(d(x)\),
- phase separation between massive and dwarf galaxies.

---

## Main conclusions

### 1. Historical claim not reproduced
The old internal claim

\[
\chi^2/\mathrm{DOF}\approx 1.01,
\qquad
d_{\min}\approx 2.63
\]

was **not reproduced on the real SPARC tables**.

It is therefore archived here as:
- an old internal proof of concept,
- non-reproducible in its current form,
- not usable as a main scientific result.

---

### 2. Best current baseline
Among the tested model families, the most successful baseline is:

\[
d = d(r/R_d)
\]

This scaled-radius variable performs significantly better than:
- absolute radius \(r\),
- baryonic acceleration \(g_{\mathrm{bar}}\),
- simple additive corrections on \(V^2\).

---

### 3. Massive-galaxy phase
The strongest result of the session is the identification of a robust **massive spiral phase**.

Main leave-one-out result on the 5 massive spirals:

\[
d_{\min}^{\mathrm{train}} = 2.4871 \pm 0.0294
\]

with typical out-of-sample performance

\[
\chi_{\mathrm{red,test}}^2 \approx 10.94 \pm 12.70
\]

This indicates that the geometric minimum near

\[
d_{\min}\approx 2.49\text{--}2.50
\]

is stable and predictive within the massive class.

---

### 4. Dwarf-galaxy phase
A separate dwarf sample suggests a distinct lower-dimensional regime:

\[
d_{\min}^{\mathrm{train}} = 2.2744 \pm 0.0341
\]

with broader dispersion and weaker predictive quality than the massive class, but still substantially better than baryons-only fits.

This supports the existence of a **distinct dwarf phase**.

---

### 5. No single universal dimension
A major conceptual result is that the data do **not** point toward a unique universal dimension.

Instead, they suggest:
- a local fit minimum around \(d_{\min}\approx 2.50\) for massive spirals,
- but a broader empirical effective regime around

\[
d_{\mathrm{eff}}\sim 2.6\text{--}2.8
\]

for much of the massive-galaxy data.

This supports the idea of an **interval of effective dimensions**, rather than a single fixed value.

---

### 6. Hierarchical compactness model
A baryonic-compactness variable improves **in-sample** fits, including difficult cases such as NGC5055.

However, it does **not** improve out-of-sample transferability.

So at the current stage, compactness should be interpreted as:
- a useful descriptive internal parameter,
- not yet a validated universal hierarchical variable.

---

## Folder structure

```text
sparc/
  README.md
  scripts/
    sparc_historical_nonreproduced_v1_0.py
    sparc_radius_baseline_v2_0.py
    sparc_scaled_radius_bump_v3_0.py
    sparc_scaled_radius_additive_v3_1.py
    sparc_gbar_bump_v4_0.py
    sparc_scaled_radius_subgroups_v3_2.py
    sparc_massive_leave_one_out_v3_3.py
    sparc_massive_oos_v3_4.py
    sparc_empirical_dimension_inversion_v3_35.py
    sparc_empirical_dimension_binned_v3_36.py
    sparc_hierarchical_compacity_v3_5.py
    sparc_hierarchical_compacity_oos_v3_5.py
    sparc_dwarf_phase_v34.py
    sparc_phase_diagram_builder_v1_0.py

  results/
    00_historical_nonreproduced/
    10_real_sparc_baselines/
    20_massive_phase/
    30_empirical_dimension/
    40_hierarchical_models/
    50_dwarf_phase/
    60_phase_diagram/

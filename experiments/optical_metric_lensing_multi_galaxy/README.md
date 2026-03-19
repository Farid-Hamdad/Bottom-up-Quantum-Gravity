# Optical Metric Lensing — Multi-Galaxy BuP Test

**Bottom-Up Quantum Gravity**  
**Farid Hamdad — 2026**

---

## Overview

This experiment explores a new phenomenological branch of the Bottom-Up program:

> instead of correcting lensing through a projected local density alone,  
> light is propagated through an **effective BuP optical response** derived from the emergent dimension profile calibrated on SPARC.

The goal is to test whether a **centralized optical BuP regime** can generate a controlled excess of gravitational deflection across several massive spiral galaxies.

This branch should be read as:

- a proof of concept for an optical extension of BuP,
- an intermediate step between SPARC phenomenology and real lensing observables,
- not yet a full first-principles derivation of an emergent optical metric.

---

## Minimal model

We start from a radial BuP dimension profile

\[
d(r)=3-\delta_d \exp\!\left[-\frac{\ln^2(r/r_0)}{2\sigma^2}\right]
\]

and define the associated geometric factor

\[
f_d(r)=\left(\frac{r}{r_{\mathrm{ref}}}\right)^{3-d(r)}.
\]

A first optical response is built through an effective radial acceleration:

\[
g_{\mathrm{eff}}(r)=g_{\mathrm{bary}}(r)\,f_d(r).
\]

However, this raw version tends to generate an excess deflection that is too external.

To control that behavior, a **centralized localized optical response** is introduced:

\[
g_{\mathrm{eff}}(r)=g_{\mathrm{bary}}(r)\Bigl[1+\lambda\,w_{\mathrm{tot}}(r)\,\bigl(f_d(r)-1\bigr)\Bigr],
\]

where \(w_{\mathrm{tot}}(r)\) combines:

- a central localization term,
- an outer exponential cutoff.

The light deflection is then estimated through

\[
\alpha(R)\simeq \frac{2}{c^2}\int_{-\infty}^{+\infty} g_\perp(R,z)\,dz,
\qquad
g_\perp(R,z)=g_{\mathrm{eff}}(r)\frac{R}{r},
\qquad
r=\sqrt{R^2+z^2}.
\]

---

## Version history

### v1.0
First optical proof of concept:
- direct effective acceleration lensing,
- real BuP excess,
- but too cumulative at large radius.

### v1.1
Improved diagnostics:
- robust ratio,
- cleaner radial interpretation,
- confirms external-drift problem.

### v1.2
Localized optical response:
- partial control of the outer tail,
- effect still too extended.

### v1.3
Centralized optical response:
- central cap on the active BuP region,
- external exponential cutoff,
- produces a much cleaner radial excess.

### v1.3 scan
Phase scan on **NGC3198** to identify a regime where the effect is:

- centralized,
- non-negligible,
- and still physically interpretable.

### multi-galaxy v1
Application of the most promising optical regime to several spiral galaxies.

---

## Main current result

The most promising centralized optical regime found in the current exploration is:

\[
(\lambda,\mathrm{width},x_{\mathrm{cap}},r_{\mathrm{cut}})=(3.0,\ 0.8,\ 4.0,\ 5.0)
\]

This regime was first identified from the v1.3 scan on **NGC3198**, then applied without re-adjustment to several SPARC massive spirals.

It generates a **robust excess of lensing deflection** typically of order:

- **+30% to +50%** in robust mean enhancement,
- with non-trivial but still reasonably centralized radial excess profiles.

---

## Multi-galaxy results

### Best centralized optical BuP regime tested

\[
(\lambda,\mathrm{width},x_{\mathrm{cap}},r_{\mathrm{cut}})=(3.0,\ 0.8,\ 4.0,\ 5.0)
\]

### Best galaxies obtained so far

| Galaxy   | Robust mean enhancement | Robust median enhancement | \(\Delta \alpha_{\max}\) [arcsec] | \(R_{\Delta\alpha,\max}\) [kpc] | \(\Delta \theta_E\) [arcsec] |
|----------|-------------------------|----------------------------|-----------------------------------|----------------------------------|------------------------------|
| NGC3198  | 1.358 | 1.331 | 0.1224 | 11.04 | 0.0055 |
| NGC5055  | 1.325 | 1.282 | 0.2243 | 8.62  | 0.0566 |
| NGC2841  | 1.390 | 1.387 | 0.5205 | 11.27 | 0.2434 |
| NGC7331  | 1.509 | 1.541 | 0.3000 | 13.94 | 0.0878 |
| UGC02885 | 1.493 | 1.543 | 0.4420 | 31.20 | 0.0620 |

### Preliminary interpretation

These results suggest that a **centralized BuP optical regime** may be quasi-stable across a subclass of **massive spiral galaxies**.

The most plausible current reading is:

1. local projected-density corrections alone are insufficient,
2. an effective optical BuP response is more promising,
3. a centralized optical regime can produce a real and robust excess of deflection,
4. this regime appears more naturally adapted to **massive spirals** than to all galaxy classes uniformly.

---

## Current interpretation

At this stage, the optical-metric branch should be interpreted as a **phenomenological optical extension** of the SPARC BuP program.

It does **not** yet establish:

- a first-principles emergent metric tensor,
- a full GR-equivalent lensing theory,
- or a direct explanation of strong-lens samples such as SLACS.

What it does establish is more modest but important:

> a Bottom-Up dimension profile can be translated into a controlled optical response that produces a centralized and significant excess of gravitational deflection on a non-trivial set of galaxies.

---

## Limitations

This branch remains exploratory.

Important limitations include:

- the optical response is not yet derived from entanglement first principles,
- the radial weight \(w_{\mathrm{tot}}(r)\) is phenomenological,
- only a small number of galaxies has been tested,
- the link to real strong-lensing systems remains open,
- one must still determine whether the same optical regime survives across broader galaxy classes.

---

## Directory structure

Suggested structure:

```text
optical_metric_lensing_multi_galaxy/
  README.md
  scripts/
    bup_optical_metric_lensing_v1.py
    bup_optical_metric_lensing_v1_1.py
    bup_optical_metric_lensing_v1_2.py
    bup_optical_metric_lensing_v1_3.py
    bup_optical_metric_lensing_v1_3_scan.py
    bup_optical_metric_lensing_multi_galaxy_v1.py
  results/
    10_single_galaxy_ngc3198/
    20_v1_3_scan_ngc3198/
    30_multi_galaxy_regime_3_08_4_5/

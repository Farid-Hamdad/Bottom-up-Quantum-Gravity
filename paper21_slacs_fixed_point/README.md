# Paper 21 — The SLACS Fixed Point

## Strong-Lensing Validation of the BuP Effective Potential

**Author:** Farid Hamdad  
**Project:** Bottom-Up Quantum Gravity  
**Year:** 2026

---

## Overview

This folder contains the numerical and observational tests associated with **Paper 21** of the Bottom-Up Quantum Gravity program.

Paper 21 tests a prediction introduced in **Paper 9**: the effective BuP gravitational potential,

$$
\mathcal{L}_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux},
$$

should leave an observable imprint in gravitational lensing data.

The test is performed on the SLACS strong-lensing sample. The main observable is the residual lensing correction

$$
C_{\rm obs} = \log\left(
\frac{\theta_E^{\rm obs}}
{\theta_E^{\rm baryon}}
\right).
$$

The central result is the discovery of a fixed-point transition around

$$
\log M_\star \simeq 11.58-11.60.
$$

At this scale:

1. the observed residual satisfies $C_{\rm obs}\simeq 0$;
2. the BuP potential $\Phi_{\rm BuP}$ locally outperforms the global proxy $\log M_\star$;
3. the dimensional sector reaches $\alpha_{\rm eff}\simeq 1$ at an intermediate diffusion scale.

This provides a non-trivial observational consistency test of the BuP effective potential.

---

## Physical origin

Paper 9 predicts the chain

$$
S_{\rm flux}
\rightarrow
\Phi_{\rm BuP}
\rightarrow
\text{gravitational response}.
$$

Paper 21 tests this chain using strong-lensing data. The BuP effective exponent is

$$
\alpha_{\rm eff}
= \frac{2d_s}{d_w}+d_w-4.
$$

The Newtonian or baryonic fixed point corresponds to

$$
\alpha_{\rm eff}=1.
$$

For standard Brownian diffusion ($d_w=2$), this condition gives

$$
d_s=3.
$$

More generally, the fixed-point condition is

$$
2d_s = d_w(5-d_w).
$$

In the SLACS graph, the fixed point is recovered not in the earliest or latest diffusion regime, but in an intermediate diffusion window.

---

## Main results

### 1. Global BuP potential signal

The best dynamical BuP feature reaches a leave-one-out improvement of approximately

$$
10.49\%
$$

on the full SLACS working sample.

A dynamic shuffle control, where the stellar-mass amplitude is permuted before constructing the graph, strongly suppresses the signal. This shows that the result depends on the correct association between stellar amplitude and graph structure.

---

### 2. Measured Sérsic indices strengthen the signal

A controlled comparison was performed on the same 61 galaxies with measured Sérsic indices.

Using measured $n_i$:

$$
\text{LOO improvement} = 10.40\%.
$$

For the exact same galaxies forced to $n=4$:

$$
\text{LOO improvement} = 9.12\%.
$$

Thus, the measured photometric morphology increases the BuP dynamical signal.

---

### 3. Transition window

The strongest transition window is

$$
11.545 < \log M_\star < 11.645.
$$

In this window, with measured Sérsic indices,

$$
\Phi_{\rm BuP}
$$

outperforms $\log M_\star$ for approximately

$$
81.8\%
$$

of galaxies.

This indicates that the BuP potential is not merely a global stellar-mass proxy. It captures a localized dynamical correction in the transition regime.

---

### 4. Observational fixed point

The zero of the observed lensing residual,

$$
C_{\rm obs}=0,
$$

is found around

$$
\log M_\star\simeq 11.58-11.60.
$$

This coincides with the mass window where $\Phi_{\rm BuP}$ locally outperforms $\log M_\star$.

---

### 5. Dimensional fixed point

A diffusion-window scan shows that the BuP dimensional fixed point is recovered at an intermediate diffusion scale.

For measured Sérsic indices, the optimal diffusion window is

$$
t_{\min}=1,\qquad t_{\max}=21.
$$

It gives

$$
\langle \alpha_{\rm eff}\rangle = 1.014,
$$

and in the fixed mass window

$$
11.545 < \log M_\star < 11.645
$$

it gives

$$
\langle \alpha_{\rm eff}\rangle = 1.014.
$$

Thus, the observational fixed point and the BuP dimensional fixed point coincide.

---

## Interpretation

Paper 21 supports the following statement:

$$
C_{\rm obs}=0,
\qquad
\Phi_{\rm BuP} > \log M_\star,
\qquad
\alpha_{\rm eff}\simeq 1
$$

all occur around the same transition scale,

$$
\log M_\star\simeq 11.6.
$$

This is interpreted as the **SLACS fixed point** of the BuP effective gravitational potential.

---

## Directory structure

```text
paper21_slacs_fixed_point/
  README.md
  paper21_slacs_fixed_point.tex
  references.bib
  data/
  scripts/
  results/
  figures/

markdown
# Paper 11 — The Dimensional Flow of BuP Gravity

## Spectral Dimension, Walk Dimension and the Effective Gravitational Exponent

**Author:** Farid Hamdad  
**Project:** Bottom-Up Quantum Gravity (BuP)  
**Year:** 2026  

---

## Overview

This folder contains the numerical and theoretical material associated with **Paper 11** of the Bottom-Up Quantum Gravity program.

Paper 11 introduces the **dimensional law** connecting the spectral properties of an entanglement graph to an **effective gravitational exponent**:

\[
\alpha_{\rm eff} = \frac{2d_s}{d_w} + d_w - 4
\]

where:

- \(d_s\) is the **spectral dimension** of the entanglement graph;
- \(d_w\) is the **walk dimension**;
- \(\alpha_{\rm eff}\) is the **effective gravitational exponent** controlling the emergent large-scale response.

This paper is the **bridge** between the microscopic graph program and the phenomenological tests performed later in the BuP sequence.

---

## Role in the BuP Program

Paper 11 provides the **dimensional backbone** for later papers:

```text
Paper 11 :
    W_ij → L_ent → d_s, d_w → α_eff

Paper 12 :
    Σ(R) → W_ij → L_ent → α_eff → V(r)

Paper 14 :
    SPARC galaxy rotation curves and LOW/HIGH regimes

Paper 21 :
    SLACS strong-lensing fixed point and α_eff ≈ 1
Thus, Paper 11 does not primarily fit a galaxy catalogue. It establishes the theoretical and numerical law that later papers test observationally.

Central Equation
The central result is:

α
e
f
f
=
2
d
s
d
w
+
d
w
−
4
α 
eff
​
 = 
d 
w
​
 
2d 
s
​
 
​
 +d 
w
​
 −4
​
 
This formula combines two graph-diffusion quantities:

P
(
t
)
∼
t
−
d
s
/
2
P(t)∼t 
−d 
s
​
 /2
 
⟨
r
2
(
t
)
⟩
∼
t
2
/
d
w
⟨r 
2
 (t)⟩∼t 
2/d 
w
​
 
 
The BuP gravitational exponent is therefore not inserted by hand. It is inferred from diffusion on the entanglement graph.

Newtonian Fixed Point
The Newtonian (or baryonic) fixed point corresponds to:

α
e
f
f
=
1
α 
eff
​
 =1
Therefore:

2
d
s
d
w
+
d
w
−
4
=
1
d 
w
​
 
2d 
s
​
 
​
 +d 
w
​
 −4=1
Equivalently:

2
d
s
=
d
w
(
5
−
d
w
)
2d 
s
​
 =d 
w
​
 (5−d 
w
​
 )
d
s
=
d
w
(
5
−
d
w
)
2
d 
s
​
 = 
2
d 
w
​
 (5−d 
w
​
 )
​
 
For Brownian diffusion:

d
w
=
2
d 
w
​
 =2
this gives:

d
s
=
3
d 
s
​
 =3
Thus, ordinary three-dimensional Newtonian behavior appears as a special fixed point of the BuP dimensional flow.

Main Numerical Outputs
Results:

text
results/
  alpha_eff_table.csv
  finite_size_summary.csv
  alpha_predictions_vs_N.csv
  paper11_summary.json
Figures:

text
figures/
  fig1_pipeline_dimensional_flow.png
  fig2_ds_vs_N.png
  fig3_dw_vs_N.png
  fig4_alpha_predictions_vs_N.png
  fig5_alpha_fixed_point_curve.png
  fig6_interpretation_regimes.png
Interpretation
Paper 11 shows that the effective gravitational behavior is controlled by the pair:

(
d
s
,
  
d
w
)
(d 
s
​
 ,d 
w
​
 )
Different regimes correspond to different gravitational responses:

Regime	Condition	Interpretation
Newtonian fixed point	
α
e
f
f
=
1
α 
eff
​
 =1	ordinary baryonic/Newtonian scaling
sub-Newtonian	
α
e
f
f
<
1
α 
eff
​
 <1	softened or under-coupled regime
super-Newtonian	
α
e
f
f
>
1
α 
eff
​
 >1	enhanced effective response
transition regime	
α
e
f
f
≈
1
α 
eff
​
 ≈1	crossover between graph phases
Reproducibility
Run:

bash
cd papers/paper11_dimensional_flow
bash scripts/run_all.sh
This generates all numerical tables and figures.

Status
Paper 11 should be read as a theoretical and numerical derivation of the BuP dimensional law.
Its observational consequences are tested later in Paper 12, Paper 14, and Paper 21.

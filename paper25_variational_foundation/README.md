# Paper 25 — Variational Foundation of Bottom-Up Quantum Gravity

**From Entanglement Equilibrium to Effective Einstein Equations**

**Titre FR :** Fondation variationnelle de Bottom-Up Quantum Gravity  
**Sous-titre FR :** De l’équilibre d’intrication aux équations d’Einstein effectives

---

## Statut

Paper 25 est le papier de synthèse théorique du programme Bottom-Up Quantum Gravity.

Son objectif est de reformuler BuP non plus comme une collection de fits ou de tests numériques séparés, mais comme une théorie variationnelle discrète de l’équilibre d’intrication.

La thèse centrale est :

\[
\boxed{
\text{Bottom-Up Quantum Gravity is a variational theory of entanglement equilibrium whose smooth local limit is an effective Einstein equation.}
}
\]

En français :

\[
\boxed{
\text{BuP est une théorie variationnelle de l'équilibre d'intrication dont la limite locale lisse est une équation d'Einstein effective.}
}
\]

La variable fondamentale est le graphe d’intrication :

\[
W_{ij}=I(i:j),
\]

où \(I(i:j)\) est l’information mutuelle entre les degrés de liberté \(i\) et \(j\).

La métrique, la courbure, la source matière, le potentiel gravitationnel, les corrections galactiques, les corrections cosmologiques et les modes gravitationnels sont ensuite reconstruits comme des structures émergentes issues de \(W\).

---

## 1. Position in the literature

Bottom-Up Quantum Gravity belongs to the broad family of approaches in which geometry is not taken as fundamental, but is reconstructed from quantum information.

Holographic entanglement entropy established a quantitative relation between boundary entanglement and bulk geometry through the Ryu--Takayanagi formula. Van Raamsdonk later argued that spacetime connectivity itself is controlled by entanglement, while tensor-network approaches showed how multi-scale entanglement structures can encode emergent geometries.

BuP shares this information-geometric motivation, but differs in its starting point. It does not assume an AdS/CFT duality, a pre-existing bulk geometry, or a fixed tensor-network architecture. Instead, its fundamental variable is the mutual-information graph

\[
W_{ij}=I(i:j),
\]

from which one constructs the entanglement Laplacian

\[
L_{\rm ent}=D-W,
\qquad
D_{ii}=\sum_j W_{ij}.
\]

The metric, curvature, matter source and weak-field propagator are treated as emergent structures derived from \(W\).

The closest conceptual antecedent is Jacobson's entanglement-equilibrium derivation of the Einstein equation. In BuP, however, the equilibrium principle is formulated directly as a discrete variational equation,

\[
\frac{\delta S_{\rm BuP}[W,\rho]}{\delta W_{ij}}=0.
\]

The smooth local limit of this equation is conjectured and numerically supported to take the form

\[
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal H_{\mu\nu}^{\rm BuP}.
\]

BuP also intersects spectral geometry and discrete Ricci curvature. Its spectral sector uses an action of the form

\[
S_{\rm spec}[W]=\mathrm{Tr}\,L_{\rm ent}^{-\beta[W]},
\]

while its curvature sector uses Ollivier--Ricci curvature on the entanglement graph,

\[
S_{\rm curv}[W]
=
\sum_{ij}W_{ij}\kappa^{OR}_{ij}[W].
\]

The continuum targets are

\[
c_NL_N\to-\Delta_g,
\]

and

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_NR_{\mu\nu}u^\mu u^\nu.
\]

Finally, BuP differs from standard modified-gravity phenomenology by deriving its weak-field, galactic and cosmological corrections from graph spectral data. The discrete Poisson equation

\[
L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}
\]

identifies

\[
G_{ij}^{\rm ent}=(L_{\rm ent}^{+})_{ij}
\]

as the entanglement Green function. The effective radial law is then controlled by

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}
+
d_w
-
4,
\]

rather than by an externally imposed modification of Newtonian gravity.

---

## 2. Fundamental postulates

Paper 25 uses the following canonical postulates.

### Postulate 1 — Entanglement graph

The fundamental object is a finite weighted graph of quantum correlations:

\[
W_{ij}=I(i:j),
\]

with

\[
W_{ij}\ge 0,
\qquad
W_{ij}=W_{ji},
\qquad
W_{ii}=0.
\]

The nodes are elementary quantum degrees of freedom. The edge weights encode mutual information.

---

### Postulate 2 — Entanglement Laplacian

From \(W\), one defines

\[
D_{ii}=\sum_j W_{ij},
\]

and

\[
L_{\rm ent}=D-W.
\]

This operator is the discrete geometric operator of BuP.

In the continuum limit, the target is

\[
c_NL_N\to-\Delta_g.
\]

---

### Postulate 3 — Emergent metric

The metric is not fundamental. It is reconstructed from the entanglement structure:

\[
g_{\mu\nu}^{\rm ent}
=
\mathcal G[W,\rho].
\]

The notation \(\mathcal G\) denotes the reconstruction map from the entanglement graph and state data to a smooth effective geometry.

---

### Postulate 4 — BuP action

The canonical BuP action is written as

\[
S_{\rm BuP}[W,\rho]
=
S_{\rm spec}[W]
+
S_{\rm curv}[W]
+
S_{\rm loc}[W]
+
S_{\rm topo}[W]
+
S_{\rm source}[W,\rho].
\]

A minimal representative is

\[
S_{\rm BuP}[W,\rho]
=
\alpha\,\mathrm{Tr}\,L_{\rm ent}^{-\beta[W]}
+
\gamma\sum_{ij}W_{ij}\kappa_{ij}^{OR}[W]
+
\lambda\sum_{ij}W_{ij}d_{ij}^{2}
+
\zeta S_{\rm topo}[W]
+
\eta S_{\rm source}[W,\rho].
\]

The terms have the following roles:

| Term | Meaning |
|---|---|
| \(S_{\rm spec}\) | spectral geometry / Einstein-Hilbert sector |
| \(S_{\rm curv}\) | discrete Ricci curvature sector |
| \(S_{\rm loc}\) | locality and infrared sector |
| \(S_{\rm topo}\) | topology and global constraints |
| \(S_{\rm source}\) | modular source / emergent stress-energy sector |

---

### Postulate 5 — Entanglement equilibrium

The fundamental dynamical condition is

\[
\boxed{
\frac{\delta S_{\rm BuP}[W,\rho]}{\delta W_{ij}}=0.
}
\]

This equation is the discrete equilibrium equation of Bottom-Up Quantum Gravity.

The continuum target is

\[
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal H_{\mu\nu}^{\rm BuP}.
}
\]

---

## 3. Discrete weak-field limit and entanglement Green function

Papers 8, 9 and 10 develop the weak-field sector of BuP.

Paper 8 constructs a candidate emergent matter source from local entanglement excitations:

\[
\delta W^{\rm loc}
\to
S_{\rm flux}.
\]

The source used in Paper 9 is

\[
S_{\rm flux}
=
T_{00}
-
\frac12 T_{aa}
+
\frac12 T_{\rm grad}
+
\|T_{0a}\|.
\]

Paper 9 then solves the discrete Poisson equation

\[
L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}.
\]

Equivalently,

\[
\Phi_{\rm BuP}
=
L_{\rm ent}^{+}S_{\rm flux},
\]

where \(L_{\rm ent}^{+}\) is the Moore--Penrose pseudoinverse of the entanglement Laplacian.

This identifies

\[
\boxed{
G_{ij}^{\rm ent}
=
(L_{\rm ent}^{+})_{ij}
}
\]

as a natural candidate for the discrete gravitational Green function of BuP.

For the best tested regime,

\[
N=20,
\qquad
\lambda=0.57,
\qquad
k=5,
\qquad
\sigma=0.15,
\]

Paper 9 found

\[
\rho_{\rm Spearman}(S_{\rm flux},|\delta R|)
=
0.741,
\qquad
p=1.84\times10^{-4},
\]

and, after solving the discrete Poisson equation,

\[
\rho_{\rm Spearman}(\Phi_{\rm BuP},|\delta R|)
=
0.738,
\qquad
p=2.01\times10^{-4}.
\]

Thus the propagated potential retains almost all of the geometric information contained in the source.

Together with the spectral convergence results of Paper 16,

\[
L_{\rm ent}\to-\Delta_g,
\]

this supports the weak-field limit

\[
L_{\rm ent}^{+}
\to
(-\Delta_g)^{-1}.
\]

Therefore,

\[
\Phi_{\rm BuP}
=
L_{\rm ent}^{+}S_{\rm flux}
\]

is the discrete BuP precursor of the continuum gravitational potential.

---

## 4. Effective propagator law

Paper 11 identifies the effective radial law of the BuP propagator.

The central relation is

\[
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}
+
d_w
-
4.
}
\]

Here:

- \(d_s\) is the spectral dimension;
- \(d_w\) is the walk dimension;
- \(\alpha_{\rm eff}\) controls the effective power law of the gravitational propagator.

The Newtonian fixed point is

\[
\alpha_{\rm eff}=1.
\]

For the diffusive value

\[
d_w=2,
\]

this gives

\[
d_s=3.
\]

Therefore, the ordinary Newtonian regime corresponds to the fixed point

\[
\boxed{
(d_s,d_w)=(3,2).
}
\]

This relation is one of the central bridges between microscopic entanglement geometry and macroscopic gravity in BuP.

---

## 5. Cosmological dimensional sector

Papers 2, 3 and 4 develop the cosmological branch of BuP.

The key idea is that the effective background dimension can vary with redshift:

\[
d(z)
=
D_C
-
\frac{\Delta d}{1+X(z)^\alpha}.
\]

The gravitational coupling is then modified through

\[
\boxed{
G_{\rm eff}(z)
=
\frac{2}{d(z)-1}G.
}
\]

Paper 2 obtained a cosmological fit with

\[
\sigma_8=0.772,
\]

\[
D_C=3.0598,
\]

\[
H_0=63.67,
\]

\[
\Omega_m=0.3448,
\]

and

\[
\Delta AIC=-6.23,
\qquad
\Delta BIC=-3.71.
\]

Paper 3 confirmed that the same dimensional mechanism reduces the \(\sigma_8\) tension relative to \(\Lambda\)CDM:

\[
\sigma_8^{\rm BuP}=0.772,
\qquad
\sigma_8^{\Lambda{\rm CDM}}\simeq0.811.
\]

Paper 4 adds a local density-dependent dimensional correction:

\[
d(z,\delta)
=
d_{\rm bg}(z)-\varepsilon\delta.
\]

The numerical comparison gives

\[
\varepsilon_{\rm JWST}=0.0113,
\]

and

\[
\varepsilon_{\rm micro}(N=20)=0.0065.
\]

Thus the cosmological branch is naturally assigned to the dimensional correction sector

\[
\mathcal H_{\mu\nu}^{\rm dim}.
\]

---

## 6. Galactic sector and SPARC

Papers 12, 13 and 14 apply the propagator law to galaxies.

Paper 12 uses the chain

\[
\Sigma(R)
\to
(d_s(r),d_w(r))
\to
\alpha_{\rm eff}(r)
\to
V(r).
\]

This connects baryonic structure to an effective radial gravitational response.

The core law remains

\[
\alpha_{\rm eff}(r)
=
\frac{2d_s(r)}{d_w(r)}
+
d_w(r)
-
4.
\]

Paper 12 showed that this mechanism can reproduce individual SPARC rotation curves. A representative result is

\[
\chi^2_{\rm red}({\rm NGC3198})=0.994.
\]

Paper 13 closes the microscopic bridge:

\[
|\Psi\rangle
\to
I_{ij}
\to
\rho_{\rm ent}(R)
\to
\Sigma(R).
\]

Its central message is:

\[
\boxed{
\text{Matter encodes entanglement; entanglement reconstructs matter.}
}
\]

Paper 14 provides the publication-grade SPARC comparison. Its constrained full-sample run over 175 galaxies gives

\[
{\rm median}\ \chi^2_{\rm red,BuP}=0.467,
\]

\[
{\rm median}\ \chi^2_{\rm red,NFW2p}=1.332,
\]

with BuP better than NFW2p on

\[
146/175
\]

galaxies, i.e.

\[
83.4\%.
\]

Thus the galactic branch is interpreted as a macroscopic expression of the dimensional propagator law.

---

## 7. Spectral continuum limit

Papers 15 and 16 establish the spectral geometry side of the Einstein limit.

The target is

\[
c_NL_N\to-\Delta_g.
\]

Paper 16 tests the low spectrum of the graph Laplacian against the Laplace-Beltrami spectrum on controlled geometries.

For the circle,

\[
\text{mean relative spectral error}=0.0051.
\]

For the flat torus,

\[
\text{mean relative spectral error}=0.0266.
\]

This supports the identification

\[
L_{\rm ent}
\to
-\Delta_g.
\]

The spectral action

\[
\mathrm{Tr}\,L_{\rm ent}^{-\beta}
\]

then admits a heat-kernel interpretation:

\[
\mathrm{Tr}(e^{-t\Delta_g})
\sim
(4\pi t)^{-D/2}
\left(
a_0+a_1t+a_2t^2+\cdots
\right).
\]

The \(a_1\) term contains the scalar curvature contribution, while higher coefficients generate higher-curvature corrections.

---

## 8. Discrete Ricci curvature limit

Paper 17 studies the Ricci curvature side of the Einstein limit.

The central relation is

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_NR_{\mu\nu}u^\mu u^\nu.
\]

At \(N=512\), Paper 17 obtained

\[
B_N=-0.301281,
\]

\[
C_N=0.391933,
\]

and

\[
A_N=\frac{1}{C_N}=2.551455.
\]

Thus the calibrated Ollivier--Ricci estimator is

\[
\widehat R_{OR}
=
2.551455
\left(
\frac{\kappa^{OR}}{\epsilon}
+
0.301281
\right).
\]

This is currently a mean-level calibrated Ricci proxy rather than a complete pointwise convergence theorem.

Its role in Paper 25 is to support the chain

\[
\kappa^{OR}
\to
R_{\mu\nu}u^\mu u^\nu.
\]

---

## 9. Modular source sector

Paper 18 develops the source side of the Einstein limit.

The key identity is the graph modular first law:

\[
\delta S_A^{\rm graph}
\simeq
\delta\langle K_A^{\rm graph}\rangle.
\]

The graph density matrix proxy is

\[
\rho_A
=
\frac{(L_A+\mu I)^{-1}}
{\mathrm{Tr}(L_A+\mu I)^{-1}},
\]

with

\[
S_A=-\mathrm{Tr}(\rho_A\log\rho_A),
\qquad
K_A=-\log\rho_A.
\]

Paper 18 found, on the flat torus,

\[
\delta S_A
=
-0.000068
+
0.989481\,\delta\langle K_A\rangle,
\]

with

\[
R^2=0.996590.
\]

On the sphere,

\[
\delta S_A
=
-0.000081
+
0.989718\,\delta\langle K_A\rangle,
\]

with

\[
R^2=0.996592.
\]

The modular source also predicts localized curvature response:

\[
R^2(|\delta K|,\langle|\Delta\kappa|\rangle_{\rm near})
=
0.967229
\]

on the flat torus, and

\[
0.982252
\]

on the sphere.

Thus

\[
\delta\langle K_A\rangle
\]

acts as a scalar modular precursor of the emergent stress-energy tensor.

The current limitation is that Paper 18 does not yet construct a full tensor

\[
T_{\mu\nu}^{\rm ent}.
\]

It establishes the modular-source precursor.

---

## 10. Effective Einstein equation

Paper 19 assembles the three continuum arrows:

\[
L_N\to\Delta_g,
\]

\[
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu,
\]

and

\[
\delta W_{\rm loc}
\to
\delta S_A
\simeq
\delta\langle K_A\rangle
\to
T_{\mu\nu}^{\rm ent}.
\]

The target equation is

\[
\boxed{
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal H_{\mu\nu}^{\rm BuP}.
}
\]

Here,

\[
G_{\mu\nu}
=
R_{\mu\nu}
-
\frac12Rg_{\mu\nu}.
\]

In the smooth, local, low-energy limit,

\[
\mathcal H_{\mu\nu}^{\rm BuP}\to0,
\]

and BuP reduces to the effective Einstein form

\[
G_{\mu\nu}
+
\Lambda_{\rm ent}g_{\mu\nu}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}.
\]

Paper 19 should be interpreted as a controlled assembly, not yet a complete theorem.

---

## 11. Correction tensor

Paper 20 studies the deviation tensor

\[
\mathcal H_{\mu\nu}^{\rm BuP}.
\]

The static correction hierarchy is

\[
\mathcal H_{\mu\nu}^{\rm static}
=
\mathcal H_{\mu\nu}^{\rm spec}
+
\mathcal H_{\mu\nu}^{\rm curv}
+
\mathcal H_{\mu\nu}^{\rm source}
+
\mathcal H_{\mu\nu}^{\rm dim}
+
\mathcal H_{\mu\nu}^{\rm nonlocal}
+
\mathcal H_{\mu\nu}^{\rm topo}
+
\mathcal H_{\mu\nu}^{\rm finite}.
\]

Paper 22 later adds the dynamical sector,

\[
\mathcal H_{\mu\nu}^{\rm dyn}.
\]

Therefore the full hierarchy is

\[
\boxed{
\mathcal H_{\mu\nu}^{\rm full}
=
\mathcal H_{\mu\nu}^{\rm static}
+
\mathcal H_{\mu\nu}^{\rm dyn}.
}
\]

Paper 20 also clarifies that the current \(H_{\rm dim}\) and \(H_{\rm nonlocal}\) quantities are phenomenological proxies, not homogeneous tensor norms.

A future goal is therefore to define a common norm such as

\[
\|\mathcal H^{(a)}\|
=
\frac{
\left(
\int d^Dx\sqrt g\,
\mathcal H^{(a)}_{\mu\nu}
\mathcal H^{(a)\mu\nu}
\right)^{1/2}
}{
\left(
\int d^Dx\sqrt g\,
G_{\mu\nu}G^{\mu\nu}
\right)^{1/2}
}.
\]

---

## 12. SLACS and strong-lensing fixed point

Paper 21 connects BuP to the SLACS strong-lensing branch.

The central result is the existence of a fixed point around

\[
\log M_\star\simeq11.58-11.60,
\]

where

\[
\langle\alpha_{\rm eff}\rangle\simeq1.014.
\]

This connects the lensing branch to the same Newtonian fixed point identified by the propagator law:

\[
\alpha_{\rm eff}=1.
\]

Thus Paper 21 is naturally assigned to the nonlocal/optical sector:

\[
\mathcal H_{\mu\nu}^{\rm nonlocal}.
\]

---

## 13. Gravitational-wave sector

Paper 22 develops the dynamical tensor sector of BuP.

The main chain is

\[
W_{ij}
\to
H_{\rm edge}
\to
q_n(t)
\to
\text{gravitational-wave-like modes}.
\]

In the true-edge Hessian formulation, Paper 22 found extremely strong source-mode coupling:

\[
\overline{\mathrm{corr}}(|J_n|,q_n^{\rm peak})
=
0.999718.
\]

Finite-size scans show a stable wave-speed sector with

\[
c_{\rm edge}\simeq1,
\]

and group velocities compatible with the relativistic propagation regime.

The dynamical correction sector is therefore

\[
\mathcal H_{\mu\nu}^{\rm dyn}.
\]

Together, Papers 21 and 22 suggest a cross-scale fixed-point structure:

\[
\alpha_{\rm eff}\simeq1,
\qquad
m_{\rm eff}^2\simeq0,
\qquad
v_g\simeq1.
\]

This is the natural place to connect future tests involving solar-system constraints, GW170817-type speed constraints, and lensing fixed points.

---

## 14. Experimental modular tomography

Paper 23 connects BuP to modular tomography protocols inspired by cold-atom and quantum-simulator experiments.

The key test is whether the entanglement cut density

\[
\rho_{\rm ent}^{\rm cut}(j)
=
\sum_{k\notin A}I(j:k)
\]

tracks the inverse modular temperature profile

\[
T_{\rm mod}(j)=\frac{1}{\beta_j}.
\]

This supports the possibility that \(W_{ij}\) is not only a theoretical object, but a reconstructible experimental observable through mutual-information tomography.

Thus Paper 23 provides an experimental direction for testing the microscopic premise of BuP:

\[
W_{ij}=I(i:j).
\]

---

## 15. Target theorem

The mathematical target of Paper 25 can be stated as follows.

Let \((W_N,\rho_N)\) be a sequence of finite entanglement graphs satisfying

\[
W_{ij}\ge0,
\qquad
W_{ij}=W_{ji},
\qquad
W_{ii}=0.
\]

Let

\[
L_N=D_N-W_N.
\]

Assume that there exists a compact smooth manifold \((M,g)\), a scale \(\epsilon_N\to0\), and a normalization \(c_N\) such that

\[
c_N\frac{D_N-W_N}{\epsilon_N}
\to
-\Delta_g
\]

in a spectral, heat-kernel, or quadratic-form sense.

Assume also that the calibrated Ollivier--Ricci curvature satisfies

\[
\frac{\kappa_N^{OR}}{\epsilon_N}
=
B_N
+
C_NR_{\mu\nu}u^\mu u^\nu
+
o(1),
\]

and that local entanglement perturbations satisfy a modular first law,

\[
\delta S_A^{(N)}
=
\delta\langle K_A^{(N)}\rangle
+
o(1).
\]

Then the target statement is that the discrete equilibrium equation

\[
\frac{\delta S_{\rm BuP}[W_N,\rho_N]}{\delta W_{ij}}=0
\]

admits a smooth local continuum limit of the form

\[
G_{\mu\nu}[g^{\rm ent}]
+
\Lambda_{\rm ent}g_{\mu\nu}^{\rm ent}
=
8\pi G_{\rm eff}T_{\mu\nu}^{\rm ent}
+
\mathcal H_{\mu\nu}^{\rm BuP},
\]

with

\[
\mathcal H_{\mu\nu}^{\rm BuP}\to0
\]

in the smooth, local, low-energy and large-\(N\) limit.

This is not yet claimed as a proven theorem. It is the theorem target supported by the numerical and phenomenological chain of Papers 2--23.

---

## 16. Open mathematical debts

Paper 25 identifies the following theoretical debts.

### 1. Derivation of \(c_N\)

The normalization in

\[
c_NL_N\to-\Delta_g
\]

must be derived analytically as a function of dimension, density, graph scale and sampling.

---

### 2. Tensor reconstruction from directional Ricci

Paper 17 currently supports

\[
\kappa^{OR}\to R_{\mu\nu}u^\mu u^\nu.
\]

A full derivation requires reconstructing

\[
R_{\mu\nu}
\]

and then

\[
G_{\mu\nu}
=
R_{\mu\nu}
-
\frac12Rg_{\mu\nu}.
\]

---

### 3. Full stress-energy tensor

Paper 18 gives a modular scalar source,

\[
\delta\langle K_A\rangle,
\]

but a complete construction of

\[
T_{\mu\nu}^{\rm ent}
\]

requires inverting relations of the form

\[
\delta\langle K_A\rangle
=
\int_A
\xi^\mu
\delta T_{\mu\nu}^{\rm ent}
d\Sigma^\nu.
\]

---

### 4. Homogeneous correction norm

Paper 20 separates controlled residuals and phenomenological proxies.

A common tensor norm for all sectors of

\[
\mathcal H_{\mu\nu}
\]

is still required.

---

### 5. Micro-to-macro closure

Paper 13 supports

\[
|\Psi\rangle
\to
I_{ij}
\to
\rho_{\rm ent}(R)
\to
\Sigma(R),
\]

but a general theorem connecting quantum states, mutual-information graphs and astrophysical matter profiles remains open.

---

### 6. Cross-scale constraints

Future work must unify the fixed-point constraints from:

\[
\text{solar system},
\qquad
\text{SPARC},
\qquad
\text{SLACS},
\qquad
\text{GW speed},
\qquad
\text{cosmology}.
\]

The relevant unifying quantity is expected to be

\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}
+
d_w
-
4.
\]

---

## 17. Summary in one sentence

\[
\boxed{
\text{Paper 25 reformulates Bottom-Up Quantum Gravity as a variational theory of entanglement equilibrium, whose weak-field, galactic, cosmological, lensing, gravitational-wave and Einstein limits all descend from the mutual-information graph }W_{ij}=I(i:j).
}
\]

---

## Repository structure

```text
paper25_variational_foundation/
  README.md
  paper25_variational_foundation.tex

  scripts/
    paper25_build_foundation_tables_v1.py

  results/
    paper25_foundation_summary_v1/
      paper25_chain_table.csv
      paper25_open_debts_table.csv
      paper25_canonical_dictionary.csv
      paper25_summary.json

  figures/
    fig01_variational_chain.png
    fig02_discrete_to_continuum_dictionary.png
    fig03_correction_tensor_hierarchy.png
    fig04_cross_scale_fixed_points.png

  notes/
    roadmap.md
    open_problems.md
    referee_notes.md
    reproducibility.md
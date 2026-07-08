# Paper 25 — Open Problems

## 1. Normalization of the entanglement Laplacian

Current target:

\[
c_N L_N \to -\Delta_g.
\]

Open problem:

Derive \(c_N\) analytically as a function of graph scale, sampling density, dimension, kernel bandwidth and normalization convention.

This is the highest-priority mathematical debt.

---

## 2. Mutual-information graph convergence

BuP assumes:

\[
W_{ij}=I(i:j).
\]

Open problem:

Show under what physical and information-theoretic assumptions the mutual-information graph approximates a heat-kernel or diffusion-kernel graph.

Target structure:

\[
W_{ij}
\sim
\exp\left[-\frac{d_g(x_i,x_j)^2}{4\epsilon_N}\right].
\]

---

## 3. Ricci tensor reconstruction

Paper 17 supports:

\[
\frac{\kappa^{OR}}{\epsilon}
\simeq
B_N+C_NR_{\mu\nu}u^\mu u^\nu.
\]

Open problem:

Upgrade the mean directional Ricci proxy into a local tensor reconstruction:

\[
R_{\mu\nu}u^\mu u^\nu
\to
R_{\mu\nu}
\to
G_{\mu\nu}.
\]

---

## 4. Full modular stress-energy tensor

Paper 18 supports:

\[
\delta S_A
\simeq
\delta\langle K_A\rangle.
\]

Open problem:

Construct a full tensor:

\[
T_{\mu\nu}^{\rm ent}.
\]

The expected continuum relation is:

\[
\delta\langle K_A\rangle
=
\int_A
\xi^\mu
\delta T_{\mu\nu}^{\rm ent}
d\Sigma^\nu.
\]

The inversion of this relation remains open.

---

## 5. Tensor norm for the correction sector

Paper 20 defines:

\[
\mathcal H_{\mu\nu}^{\rm static}
=
\mathcal H^{\rm spec}
+
\mathcal H^{\rm curv}
+
\mathcal H^{\rm source}
+
\mathcal H^{\rm dim}
+
\mathcal H^{\rm nonlocal}
+
\mathcal H^{\rm topo}
+
\mathcal H^{\rm finite}.
\]

Paper 22 adds:

\[
\mathcal H_{\mu\nu}^{\rm dyn}.
\]

Open problem:

Define a homogeneous norm for all sectors:

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

## 6. Micro-to-macro closure

Paper 13 supports:

\[
|\Psi\rangle
\to
I_{ij}
\to
\rho_{\rm ent}(R)
\to
\Sigma(R).
\]

Open problem:

Prove a general theorem connecting quantum state structure, mutual-information graphs and astrophysical matter profiles.

---

## 7. Cross-scale fixed-point unification

BuP contains several fixed-point signals:

\[
\alpha_{\rm eff}\simeq1,
\]

\[
m_{\rm eff}^2\simeq0,
\]

\[
v_g\simeq1.
\]

Open problem:

Unify solar-system, SPARC, SLACS, gravitational-wave and cosmological constraints in one correction-tensor framework.

---

## 8. Experimental reconstructability

Paper 23 suggests that \(W_{ij}=I(i:j)\) may be reconstructible through modular tomography.

Open problem:

Define a reproducible experimental protocol connecting measured mutual information to the BuP graph \(W\).

---

## Priority ranking

| Priority | Problem |
|---|---|
| Critical | derive \(c_N\) |
| Critical | construct \(T_{\mu\nu}^{\rm ent}\) |
| High | Ricci tensor reconstruction |
| High | correction tensor norm |
| Medium-high | cross-scale fixed-point unification |
| Medium | micro-to-macro closure |
| Medium | experimental tomography |

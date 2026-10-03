# Paper 30 — Preregistration Frozen

## Title

**Robustness of Information-Sector Compensation under Amplitude Damping**

## Status

**PREREGISTRATION — TO BE FROZEN BEFORE ANY PAPER 30 RESULT IS COMPUTED**

This document defines the prospective protocol for Paper 30.

Paper 29 is the discovery study. Paper 30 is an independent prospective validation under a qualitatively different irreversible channel and includes a system size not used to discover the compensation structure.

No Paper 30 scientific result may be inspected before this preregistration and the implementation script are frozen in version control.

---

## 1. Scientific question

Paper 29 found, under local dephasing, a robust residual complementarity between:

- the mutual-information geometry \(W\),
- and the modular sector derived from \(K_A=-\log\rho_A\).

After removing the smooth common dependence on the irreversible trajectory and standardizing the residuals, Paper 29 found approximately

\[
z(\epsilon_{\rm mod})
\simeq
-z(\epsilon_W),
\]

with a low-variance quasi-invariant

\[
\mathcal C_\alpha
=
z(\epsilon_{\rm mod})
+
\alpha z(\epsilon_W),
\qquad
\alpha \approx 1.
\]

Paper 30 asks:

> **Does this residual compensation structure survive under local amplitude damping, and does it remain present at a previously untested system size?**

The goal is prospective validation, not post-hoc rediscovery.

---

## 2. Primary hypothesis

The primary Paper 30 prediction is

\[
z(\epsilon_{\rm mod})
+
\alpha z(\epsilon_W)
\approx 0,
\]

with a compensation coefficient \(\alpha\) close to unity.

The preregistered primary acceptance window is

\[
0.8
\le
\alpha
\le
1.2.
\]

This interval is frozen before execution and must not be changed after unblinding.

---

## 3. Systems and validation hierarchy

### 3.1 Primary prospective validation system

The primary Paper 30 validation system is

\[
N=12,
\qquad
A=\{0,1,2,3\}.
\]

N=12 was not used to discover the Paper 29 compensation pattern and is therefore the primary prospective size validation.

The six fixed unordered pairs inside \(A\) are

\[
(0,1),\,
(0,2),\,
(0,3),\,
(1,2),\,
(1,3),\,
(2,3).
\]

No pair selection may be changed after execution.

### 3.2 Direct replication systems

The following systems are retained as direct replication controls:

\[
N=8,
\qquad
N=10,
\]

with the same subsystem

\[
A=\{0,1,2,3\}.
\]

These sizes were used in Paper 29 and therefore do not constitute new size validation.

Their role in Paper 30 is to test channel robustness under amplitude damping.

### 3.3 Hierarchy

Paper 30 therefore separates:

- **new-size prospective validation:** N=12;
- **same-size channel replication:** N=8 and N=10.

The primary Paper 30 classification is based on N=12.

N8 and N10 are secondary replication systems.

---

## 4. Initial state and Hamiltonian

Paper 30 must use the same initial-state preparation protocol as Paper 29 unless a pre-execution implementation audit reveals a literal incompatibility.

The transverse-field Ising model is

\[
H
=
-J\sum_i Z_i Z_{i+1}
-
h\sum_i X_i,
\]

with

\[
J=h=1.
\]

The initial product state is

\[
|+\rangle^{\otimes N}.
\]

The interacting state preparation must reproduce the same frozen Paper 29 / Paper 28 convention:

- Strang splitting;
- \(p=6\);
- \(\Delta t=0.35\);
- total preparation time \(t=2.1\);
- same boundary-condition convention as Paper 29;
- same time-reversal mixture convention as Paper 29.

Recovery checks must be performed before Paper 30 observables are accepted.

---

## 5. Irreversible dynamics

Paper 30 uses independent local amplitude damping.

For one qubit,

\[
\Phi_\gamma(\rho)
=
E_0\rho E_0^\dagger
+
E_1\rho E_1^\dagger,
\]

with

\[
E_0=
\begin{pmatrix}
1&0\\
0&\sqrt{1-\gamma}
\end{pmatrix},
\qquad
E_1=
\begin{pmatrix}
0&\sqrt{\gamma}\\
0&0
\end{pmatrix}.
\]

The global channel is

\[
\Phi_\gamma^{\otimes N}.
\]

The trajectory parameter is

\[
s_n=0.1n,
\qquad
n=0,\ldots,30,
\]

with

\[
\gamma(s)=1-e^{-s}.
\]

Thus there are 31 states and 30 transitions.

Each trajectory point is evaluated directly at the preregistered \(\gamma(s_n)\) from the same prepared initial density matrix.

---

## 6. Primary progress variable

The primary progress variable is

\[
\boxed{
\gamma(s)=1-e^{-s}
}
\]

rather than \(Q=N\ln2-S(\rho)\).

This choice is frozen because amplitude damping is non-unital and the maximally mixed state is not its natural fixed point.

The primary detrending analysis is therefore performed against \(\gamma\).

The auxiliary parameter \(s\) is retained only as the generator coordinate used to define the grid.

No primary Paper 30 conclusion may depend on replacing \(\gamma\) by another trajectory coordinate after unblinding.

---

## 7. Secondary regularized relative-entropy diagnostic

Paper 30 does not use quantum relative entropy as a primary progress variable.

A secondary diagnostic is nevertheless preregistered using a full-rank regularized reference state:

\[
\sigma_\varepsilon
=
(1-\varepsilon)
|0\cdots0\rangle\langle0\cdots0|
+
\varepsilon\frac{I}{2^N}.
\]

Define

\[
D_\varepsilon(\rho\Vert\sigma_\varepsilon)
=
\operatorname{Tr}
\left[
\rho
\left(
\log\rho-\log\sigma_\varepsilon
\right)
\right].
\]

The regularization parameter is frozen as

\[
\varepsilon=10^{-6}.
\]

This quantity is:

- secondary only;
- not used in the primary success rule;
- not used to define or tune bandwidths;
- not used to replace \(\gamma\) after unblinding.

Its purpose is only to test whether the qualitative trajectory ordering is compatible with a regularized distance toward the amplitude-damping fixed state.

Any sensitivity analysis in \(\varepsilon\) must be labeled post hoc unless separately preregistered before execution.

---

## 8. Descriptive physical observables

The following are tracked descriptively.

### 8.1 Von Neumann entropy

\[
S(\rho)
=
-\operatorname{Tr}(\rho\log\rho).
\]

### 8.2 Excitation content

\[
N_{\rm exc}
=
\sum_i
\frac{I-Z_i}{2},
\]

\[
E_{\rm exc}(\gamma)
=
\operatorname{Tr}
\left[
\rho(\gamma)N_{\rm exc}
\right].
\]

These are descriptive controls and are not primary success criteria.

---

## 9. Mutual-information sector

For every trajectory point, compute

\[
W_{ij}(\gamma)
=
I(i:j)
=
S(\rho_i)+S(\rho_j)-S(\rho_{ij})
\]

for the six fixed pairs in \(A\).

Construct

\[
w(\gamma)
=
\left(
W_{01},
W_{02},
W_{03},
W_{12},
W_{13},
W_{23}
\right).
\]

Retain

\[
M_W(\gamma)=\|w(\gamma)\|_2
\]

and

\[
\hat w(\gamma)
=
\frac{w(\gamma)}
{\|w(\gamma)\|_2}.
\]

If

\[
\|w\|_2
\le
10^{-12},
\]

the normalized direction is undefined and must be explicitly recorded.

---

## 10. Modular sector

For the fixed subsystem \(A\),

\[
\rho_A(\gamma)
=
\operatorname{Tr}_{\bar A}\rho(\gamma),
\]

\[
K_A(\gamma)
=
-\log\rho_A(\gamma).
\]

Use the same modular normalization convention as Paper 29:

\[
\widetilde K_A
=
\frac{
K_A-\langle\kappa\rangle I
}{
\operatorname{std}(\kappa)
}.
\]

For each fixed pair,

\[
Q_{ij}^{ab}
=
\left[
\left[
\widetilde K_A,
\sigma_i^a
\right],
\sigma_j^b
\right],
\]

\[
A_{ij}
=
\frac{1}{9}
\sum_{a,b}
\operatorname{Tr}
\left[
\rho_A
(Q_{ij}^{ab})^\dagger
Q_{ij}^{ab}
\right],
\]

\[
v_{ij}
=
\sqrt{A_{ij}}.
\]

Construct

\[
v(\gamma)
=
(v_{01},v_{02},v_{03},v_{12},v_{13},v_{23}),
\]

with

\[
M_{\rm mod}(\gamma)
=
\|v(\gamma)\|_2
\]

and

\[
\hat v(\gamma)
=
\frac{v(\gamma)}
{\|v(\gamma)\|_2}.
\]

If the norm is \(\le10^{-12}\), the normalized direction is undefined and must be retained as such.

---

## 11. Transition magnitudes

For every adjacent trajectory step,

\[
d_W(n)
=
\left\|
\hat w(\gamma_{n+1})
-
\hat w(\gamma_n)
\right\|_2,
\]

and

\[
d_{\rm mod}(n)
=
\left\|
\hat v(\gamma_{n+1})
-
\hat v(\gamma_n)
\right\|_2.
\]

These 30 transition magnitudes per system are the inputs to the compensation analysis.

Undefined transitions must not be silently removed.

---

## 12. Primary detrending procedure

The common smooth dependence on the irreversible trajectory is removed separately from \(d_W\) and \(d_{\rm mod}\).

The primary method is frozen as leave-one-out local-linear Gaussian kernel regression against \(\gamma\).

For each

\[
y\in\{d_W,d_{\rm mod}\},
\]

fit

\[
y_n
=
f(\gamma_n)
+
\epsilon_n.
\]

The preregistered bandwidth fractions of the total \(\gamma\)-range are

\[
0.05,\,
0.10,\,
0.20,\,
0.30,\,
0.40.
\]

All five bandwidths must be reported.

No preferred bandwidth may be selected after unblinding.

---

## 13. Residual standardization

For each system and bandwidth separately,

\[
z(\epsilon_W)
=
\frac{
\epsilon_W-\langle\epsilon_W\rangle
}{
\operatorname{std}(\epsilon_W)
},
\]

\[
z(\epsilon_{\rm mod})
=
\frac{
\epsilon_{\rm mod}
-
\langle\epsilon_{\rm mod}\rangle
}{
\operatorname{std}(\epsilon_{\rm mod})
}.
\]

No clipping, trimming, robust rescaling, transition deletion, or hidden filtering is allowed in the primary analysis.

---

## 14. Compensation coefficient

Define

\[
\mathcal C_\alpha
=
z(\epsilon_{\rm mod})
+
\alpha z(\epsilon_W).
\]

For each system and bandwidth,

\[
\alpha^\star
=
\arg\min_\alpha
\operatorname{Var}
(\mathcal C_\alpha).
\]

Compute:

- \(\alpha^\star_{N12}\);
- \(\alpha^\star_{N8}\);
- \(\alpha^\star_{N10}\).

Also compute pooled coefficients for descriptive cross-size comparisons, but the primary classification is based on N=12.

---

## 15. Primary endpoints — N=12

Paper 30 has three primary endpoints, all evaluated prospectively on N=12.

### P1 — residual sign

For every preregistered bandwidth,

\[
\operatorname{Pearson}
(\epsilon_W,\epsilon_{\rm mod})
<
0.
\]

### P2 — compensation coefficient

For every preregistered bandwidth,

\[
0.8
\le
\alpha^\star_{N12}
\le
1.2.
\]

### P3 — quasi-invariant variance suppression

Define

\[
R_C
=
\frac{
\operatorname{Var}
(\mathcal C_{\alpha^\star_{N12}})
}{
\operatorname{Var}
[z(\epsilon_{\rm mod})]
}.
\]

The preregistered target is

\[
R_C
\le
0.20
\]

for every frozen bandwidth.

This corresponds to at least 80% suppression of standardized modular residual variance.

---

## 16. Primary success rule

The Paper 30 prospective validation is classified as **SUPPORTED** only if all of the following hold for N=12:

1. residual Pearson correlation is negative at all five frozen bandwidths;
2. \(\alpha^\star_{N12}\in[0.8,1.2]\) at all five frozen bandwidths;
3. \(R_C\le0.20\) at all five frozen bandwidths;
4. no transition is removed or bandwidth selected after unblinding.

If any condition fails, the exact pattern must be reported.

A failed prediction remains a valid scientific result.

---

## 17. N8/N10 replication endpoints

For N=8 and N=10, repeat the full Paper 30 pipeline under amplitude damping.

Report for all bandwidths:

- residual Pearson coefficient;
- \(\alpha^\star\);
- quasi-invariant variance ratio;
- parameter-free \(\alpha=1\) control;
- cross-size transfer.

These are secondary replication endpoints and do not override the N=12 primary classification.

---

## 18. Cross-size transfer

For every bandwidth:

- learn \(\alpha^\star_{N12}\) and apply it unchanged to N8 and N10;
- learn \(\alpha^\star_{N8}\) and apply it to N12;
- learn \(\alpha^\star_{N10}\) and apply it to N12.

Report the resulting quasi-invariant variance ratios.

No hard threshold is assigned to cross-size transfer for primary classification.

---

## 19. Parameter-free alpha = 1 control

Evaluate

\[
\mathcal C_1
=
z(\epsilon_{\rm mod})
+
z(\epsilon_W)
\]

for all systems and all bandwidths.

Report

\[
R_{C,1}
=
\frac{
\operatorname{Var}(\mathcal C_1)
}{
\operatorname{Var}[z(\epsilon_{\rm mod})]
}.
\]

This is a preregistered secondary control.

---

## 20. Raw-sector controls

Retain and report:

- \(M_W(\gamma)\);
- \(M_{\rm mod}(\gamma)\);
- \(d_W(n)\);
- \(d_{\rm mod}(n)\);
- \(S(\rho(\gamma))\);
- \(E_{\rm exc}(\gamma)\);
- \(D_\varepsilon(\rho\Vert\sigma_\varepsilon)\);
- all six \(W_{ij}(\gamma)\);
- all six \(v_{ij}(\gamma)\).

The compensation claim must not rely solely on normalized quantities.

---

## 21. Numerical safeguards

Frozen before execution:

- state normalization tolerance: \(10^{-12}\);
- zero-vector threshold: \(10^{-12}\);
- entropy eigenvalue cutoff: \(10^{-15}\);
- strict positivity checks before modular logarithms;
- no hidden clipping beyond documented numerical safeguards;
- no silent NaN deletion;
- no silent transition deletion;
- all undefined cases counted and reported.

Any required pre-unblinding implementation correction must be recorded in a numbered amendment.

---

## 22. Recovery controls

Before scientific results are accepted:

1. initial prepared states for N8 and N10 must reproduce the frozen Paper 29 / Paper 28 recovery references within the established tolerance;
2. N12 preparation must pass internal normalization, Hermiticity, and positivity checks;
3. at \(\gamma=0\), amplitude damping must reduce to the identity channel within floating-point precision;
4. \(\operatorname{Tr}\rho(\gamma)=1\) at every step;
5. \(\rho(\gamma)\) must remain Hermitian and positive semidefinite within numerical tolerance;
6. excitation content must not increase under the implemented amplitude-damping trajectory, up to numerical tolerance.

A failed recovery gate invalidates the run before scientific interpretation.

---

## 23. Interpretation boundary

If the primary N12 prediction is supported, Paper 30 may conclude only that:

> the residual compensation structure discovered under dephasing generalizes prospectively to local amplitude damping at a previously untested finite system size.

Paper 30 must not claim from this test alone that:

- information is fundamentally conserved between \(W\) and \(K_A\);
- \(W\) and \(K_A\) are the same physical quantity;
- a literal information-transfer law has been established;
- the thermodynamic arrow emerges from BuP;
- a black-hole horizon has been reproduced;
- the relation holds in the thermodynamic limit;
- \(\alpha\) is a universal constant.

Those claims require separate tests.

---

## 24. Relation to Paper 29

Paper 29 is the discovery dataset:

- irreversible channel: local dephasing;
- sizes used for discovery: N8 and N10.

Paper 30 is the prospective validation dataset:

- irreversible channel: local amplitude damping;
- new primary validation size: N12;
- direct same-size channel replications: N8 and N10.

This separation must be maintained in reporting.

---

## 25. Post-unblinding rule

After execution:

- preregistered endpoints must be reported first;
- any additional analysis must be labeled `POST_HOC_EXPLORATORY`;
- exploratory analyses may not replace failed primary endpoints;
- negative, null, undefined, and unexpected outcomes must be retained.

---

## 26. Freeze requirement

Before any Paper 30 scientific result is generated:

1. this preregistration must be committed;
2. the implementation script must be committed;
3. SHA256 hashes must be recorded;
4. the Git working tree must be clean;
5. execution may begin only after these conditions are satisfied.

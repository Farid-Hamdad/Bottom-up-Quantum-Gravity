# Paper 29 — Preregistration FROZEN

## Title
**Information Loss and the Arrow of Modular Time**

## Status
**PREREGISTRATION FROZEN — NO PAPER 29 RESULT USED**

No primary metric, grid, normalization, control, exclusion rule, or interpretation criterion may be modified after unblinding.

## Scientific question
Under a prospectively fixed irreversible CPTP evolution, do informational geometry and modular propagation evolve coherently along the independently defined direction of information loss?

This paper does **not** test whether irreversibility itself emerges from BuP.

## Distinct parameters
- \(s\): irreversible-channel parameter.
- \(\tau\): modular-flow parameter at fixed \(s\).

\[
s \neq \tau.
\]

## Primary system
\[
N=8,\qquad A=\{0,1,2,3\}.
\]

Frozen pair ordering:

\[
(0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
\]

## Initial state
Inherited exactly from Paper 28:

\[
H_{\rm TFIM}=-J\sum_i Z_iZ_{i+1}-h\sum_i X_i,
\]

with

\[
J=h=1,\qquad p=6,\qquad dt=0.35,\qquad t=2.1.
\]

Initial pure state:

\[
|+\rangle^{\otimes N}.
\]

Time-reversal mixture:

\[
\rho_0=
\frac12\left[
|\Psi(+t)\rangle\langle\Psi(+t)|
+
|\Psi(-t)\rangle\langle\Psi(-t)|
\right].
\]

## Primary irreversible channel
Independent local dephasing:

\[
\Phi_s(\rho)
=
\frac{1+e^{-s}}2\rho
+
\frac{1-e^{-s}}2 Z\rho Z.
\]

Global evolution:

\[
\rho(s)=\Phi_s^{\otimes N}[\rho_0].
\]

## Frozen irreversible grid
\[
s_n=0.1n,\qquad n=0,\ldots,30.
\]

Thus:

\[
0\le s\le3,
\]

with 31 states and 30 transitions.

## Reference state
\[
\sigma_N=\frac{I}{2^N}.
\]

Information coordinate:

\[
Q(s)=D[\rho(s)\|\sigma_N]
=N\log2-S[\rho(s)].
\]

## Irreversible orientation
\[
\Sigma_n=Q(s_n)-Q(s_{n+1}).
\]

Expected monotonicity:

\[
\Sigma_n\ge0.
\]

This direction is defined independently of \(W\) and the modular sector.

## Informational geometry
For each \(s\),

\[
W_{ij}(s)=I(i:j).
\]

Frozen six-component representation:

\[
w(s)=
\begin{pmatrix}
I_{01}(s)\\
I_{02}(s)\\
I_{03}(s)\\
I_{12}(s)\\
I_{13}(s)\\
I_{23}(s)
\end{pmatrix}.
\]

Normalize by:

\[
\widehat w(s)=\frac{w(s)}{\|w(s)\|_2}.
\]

Raw magnitude retained:

\[
M_W(s)=\|w(s)\|_2.
\]

## Modular sector
\[
\rho_A(s)=\operatorname{Tr}_{\bar A}\rho(s),
\]

\[
K_A(s)=-\log\rho_A(s).
\]

Paper 28 normalization retained exactly:

\[
\widetilde K_A(s)
=
\frac{
K_A(s)-\langle\kappa(s)\rangle I
}{
\operatorname{std}[\kappa(s)]
}.
\]

## Modular propagation observable
For each pair:

\[
C_{ij}(\tau;s)
=
\frac19
\sum_{a,b=X,Y,Z}
\operatorname{Tr}
\left[
\rho_A(s)
[\sigma_i^a(\tau;s),\sigma_j^b]^\dagger
[\sigma_i^a(\tau;s),\sigma_j^b]
\right].
\]

The short-time coefficient is computed analytically, not by fitting:

\[
Q_{ij}^{ab}(s)
=
[
[
\widetilde K_A(s),\sigma_i^a
],\sigma_j^b
],
\]

\[
A_{ij}(s)
=
\frac19
\sum_{a,b}
\operatorname{Tr}
[
\rho_A(s)
Q_{ij}^{ab\dagger}(s)
Q_{ij}^{ab}(s)
].
\]

Then

\[
v_{ij}^{\rm mod}(s)=\sqrt{A_{ij}(s)}.
\]

Frozen modular vector:

\[
v(s)=
\begin{pmatrix}
v_{01}^{\rm mod}(s)\\
v_{02}^{\rm mod}(s)\\
v_{03}^{\rm mod}(s)\\
v_{12}^{\rm mod}(s)\\
v_{13}^{\rm mod}(s)\\
v_{23}^{\rm mod}(s)
\end{pmatrix}.
\]

Normalize by:

\[
\widehat v(s)=\frac{v(s)}{\|v(s)\|_2}.
\]

Raw magnitude retained:

\[
M_{\rm mod}(s)=\|v(s)\|_2.
\]

## Incremental displacements
\[
\delta\widehat w_n
=
\widehat w(s_{n+1})-\widehat w(s_n),
\]

\[
\delta\widehat v_n
=
\widehat v(s_{n+1})-\widehat v(s_n).
\]

Magnitudes:

\[
d_W(n)=\|\delta\widehat w_n\|_2,
\]

\[
d_{\rm mod}(n)=\|\delta\widehat v_n\|_2.
\]

## Zero-step rule
Frozen threshold:

\[
\varepsilon_{\rm step}=10^{-12}.
\]

A transition is directionally valid only if both step norms exceed this threshold.

Otherwise:

\[
\chi_n=\mathrm{NaN},
\]

with status:

`DIRECTION_UNDEFINED_ZERO_STEP`.

No undefined transition may be silently discarded.

## Primary directional statistic
For valid transitions:

\[
\chi_n
=
\frac{
\delta\widehat w_n\cdot\delta\widehat v_n
}{
\|\delta\widehat w_n\|_2
\|\delta\widehat v_n\|_2
}.
\]

Primary endpoint:

\[
\overline\chi
=
\frac1{N_{\rm valid}}
\sum_{n\in\mathcal V}\chi_n.
\]

Mandatory directional descriptor:

\[
f_+
=
\frac{
\#\{n\in\mathcal V:\chi_n>0\}
}{
N_{\rm valid}
}.
\]

Also report:

\[
N_{\rm valid},
\qquad
N_{\rm undefined}.
\]

## Secondary analyses
Frozen secondary quantities:

\[
\rho_S(d_W,d_{\rm mod}),
\]

\[
\rho_S(\Sigma,d_W),
\]

\[
\rho_S(\Sigma,d_{\rm mod}),
\]

\[
D_W^{(0)}(s)
=
\|\widehat w(s)-\widehat w(0)\|_2,
\]

\[
D_{\rm mod}^{(0)}(s)
=
\|\widehat v(s)-\widehat v(0)\|_2,
\]

plus raw magnitudes \(M_W(s)\) and \(M_{\rm mod}(s)\).

## Controls

### Control 1 — exact Paper 28 recovery
At \(s=0\), Paper 29 must recover the corresponding Paper 28 N8/A4 values within the frozen implementation tolerance.

### Control 2 — reversible TFIM trajectory
\[
U(s)=e^{-isH_{\rm TFIM}},
\]

\[
\rho_U(s)=U(s)\rho_0U^\dagger(s).
\]

Using the same grid, the maximally mixed reference implies:

\[
D[\rho_U(s)\|I/2^N]=\text{constant},
\]

so

\[
\Sigma_n^{(U)}\simeq0.
\]

### Control 3 — raw vs normalized quantities
Both raw and normalized \(w\) and \(v\) must be retained.

### Control 4 — N10/A4 replication
Replication system:

\[
N=10,\qquad A=\{0,1,2,3\}.
\]

The entire analysis pipeline is reused unchanged.

## Numerical safeguards
State-vector norm tolerance:

\[
|\langle\psi|\psi\rangle-1|\le10^{-12}.
\]

Entropy eigenvalue cutoff:

\[
\lambda>10^{-15}.
\]

Strict modular positivity gate:

\[
\lambda_{\min}(\rho_A)>0.
\]

No undocumented eigenvalue clipping or regularization is permitted.

## Interpretation boundary
A positive result may support only the statement that, in the preregistered finite-system benchmark, informational geometry and modular propagation reorganize coherently along an independently imposed irreversible information-loss direction.

The following conclusions are not permitted from this experiment alone:

- “The arrow of time emerges from modular flow.”
- “The arrow of time emerges from \(W\).”
- “BuP explains the thermodynamic arrow of time.”

The irreversible direction is externally introduced through the CPTP channel.

## Freeze statement
The hypothesis, preparation, channel, reference state, grid, pair ordering, normalizations, primary statistic, zero-step rule, controls, replication, numerical safeguards, and interpretation boundaries are frozen before any Paper 29 result is inspected.

\[
\boxed{\text{NO PAPER 29 RESULT USED TO DEFINE THIS PROTOCOL}}
\]

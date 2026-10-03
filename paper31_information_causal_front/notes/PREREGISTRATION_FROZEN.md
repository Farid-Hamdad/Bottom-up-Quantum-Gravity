# Paper 31 — Information Causal Front and Emergent Propagation Speed

## PREREGISTRATION_FROZEN — DRAFT v1

**Status:** prospective preregistration draft.  
**No Paper 31 scientific result may be inspected before this protocol and its implementation are frozen.**

---

## 1. Scientific question

Paper 30 showed that, for the finite-depth local TFIM preparation used in the BuP information-sector benchmark, enlarging the system from N=10 to N=12 leaves the local observables of the subsystem

A={0,1,2,3}

unchanged to floating-point precision.

This establishes the presence of a finite local domain of influence for the chosen circuit architecture.

Paper 31 asks a different and more stringent question:

> Does the physical propagation of a local information perturbation define an emergent information-front velocity v_info, distinct from the trivial circuit-support velocity, and is that velocity compatible with the calibrated BuP edge-wave speed c_edge from Paper 22?

The central comparison is therefore

v_info / c_edge ?≈ 1.

Paper 22 uses calibrated Einstein-fixed-point units in which

c_edge = 1

at eta_edge = 1.

Paper 31 does not calibrate v_info to this value.

---

## 2. Core methodological distinction

Paper 31 distinguishes two different velocities.

### 2.1 Circuit-support velocity

v_circuit is the maximum propagation speed implied purely by the finite-depth nearest-neighbor circuit architecture.

This quantity is kinematic and is not a primary scientific result.

### 2.2 Information-front velocity

v_info is extracted from the measured propagation of a finite perturbation signal in reduced-state space.

Only v_info is compared with c_edge.

A result that merely reproduces the exact finite support imposed by the circuit is not sufficient to establish the Paper 31 hypothesis.

---

## 3. Reference dynamics

The reference state preparation follows the same local TFIM protocol used in Papers 28–30:

H = -J sum_i Z_i Z_{i+1} - h sum_i X_i,

with J=1 and h=1.

The initial state is |+>^N.

Time evolution is approximated by the same second-order Strang splitting convention used previously.

The time increment is Delta t = 0.35.

For an evolution depth P,

t = P Delta t.

The system uses open nearest-neighbor boundaries.

---

## 4. System size and observed subsystem

The observed subsystem is frozen as

A={0,1,2,3}.

The primary system size is

N=16.

The choice N=16 is intended to provide sufficient spatial separation between A and perturbation sites while remaining computationally tractable for reduced-state diagnostics.

If the frozen implementation demonstrates before scientific unblinding that N=16 is computationally infeasible, any change in N requires a documented pre-unblinding amendment.

---

## 5. Perturbation protocol

A reference evolution and a perturbed evolution are compared.

The perturbation is local and applied at a single site r outside A.

The perturbation operator is frozen as a local Pauli-Z operation:

U_r = Z_r.

The perturbed initial state is

|psi_pert(0)> = Z_r |+>^N.

The reference initial state is

|psi_ref(0)> = |+>^N.

Both states are then evolved under the identical frozen TFIM evolution.

No perturbation strength or operator may be tuned after unblinding.

---

## 6. Perturbation distances

The perturbation site r is varied over all sites outside A:

r in {4,5,...,N-1}.

The graph distance from the boundary of A is

d = r - 3.

Thus d=1,2,...,N-4.

The analysis is performed as a function of d, not absolute site label.

---

## 7. Evolution-time grid

The preregistered depth grid is

P in {1,2,3,4,5,6,7,8}.

With Delta t=0.35, the corresponding evolution times are

t in {0.35,0.70,1.05,1.40,1.75,2.10,2.45,2.80}.

Every (d,t) pair is evaluated independently from the same frozen initial-state definitions.

---

## 8. Primary signal observable

For each distance d and time t, compute the reduced states

rho_A^ref(t)

and

rho_A^pert(d,t).

The primary perturbation signal is the trace distance

Delta_A(d,t) = 1/2 ||rho_A^pert(d,t)-rho_A^ref(t)||_1.

No alternative norm may replace the trace distance in the primary analysis after unblinding.

---

## 9. Secondary signal observables

Secondary diagnostics only:

1. Frobenius distance;
2. change in the six mutual-information coordinates Delta W(d,t);
3. change in the six modular coefficients Delta v(d,t).

These diagnostics may support interpretation but do not override the primary trace-distance result.

---

## 10. Removal of exact-zero circuit-support region

The exact finite-depth circuit architecture can force Delta_A(d,t)=0 outside its algebraic support.

Paper 31 does not identify this exact support boundary with the physical information front.

For each time t, distinguish:

- support boundary: furthest d with numerically nonzero signal;
- signal front: location where a fixed fraction of the observable signal amplitude is reached.

Only the signal front enters the primary v_info estimate.

---

## 11. Globally normalized signal profile

Define one global signal scale over the complete preregistered distance-time grid:

Delta_A_global_max = max_{d,t} Delta_A(d,t).

The normalized signal profile is

F(d,t) = Delta_A(d,t) / Delta_A_global_max.

The same denominator is therefore used for every distance and every evolution time.

This prevents an apparent front displacement from being generated solely by a time-dependent rescaling of the near-field signal.

If Delta_A_global_max <= 1e-14, the complete Paper 31 trajectory is classified as signal-undefined.

No time-dependent renormalization is permitted in the primary analysis.

---

## 12. Preregistered front levels

Three fixed global signal levels are used:

q in {0.10,0.25,0.50}.

For each q and time t, the front position r_q(t) is defined by the outward spatial crossing

F(d,t)=q.

Equivalently, each front tracks the fixed absolute signal level

Delta_A(d,t) = q Delta_A_global_max

throughout the complete trajectory.

Linear interpolation between the two nearest adjacent distances is used when the crossing lies between integer lattice sites.

If multiple crossings occur because of spatial oscillations, the furthest outward crossing is used.

If no crossing exists at a given time, that (q,t) point is marked undefined.

No threshold level or normalization convention may be selected after unblinding.

---

## 13. Front-velocity extraction

For each threshold q, fit

r_q(t) = v_q t + b_q

by ordinary least squares over all valid preregistered time points.

Record v_q, b_q, and R_q^2.

The primary information-front velocity is

v_info = (v_0.10 + v_0.25 + v_0.50)/3.

The spread is

sigma_v = std(v_0.10, v_0.25, v_0.50).

---

## 14. Internal front-consistency endpoint

The front is internally coherent only if all of the following hold:

1. each of the three threshold fits has at least 4 valid time points;
2. each threshold fit satisfies R_q^2 >= 0.95;
3. all three velocities are positive;
4. sigma_v / v_info <= 0.15.

If any condition fails, Paper 31 is classified as NOT INTERNALLY COHERENT.

---

## 15. Comparison with Paper 22

Paper 22 defines the calibrated true-edge propagation scale

c_edge = 1

at the Einstein fixed point eta_edge=1.

Paper 31 does not use c_edge during front extraction.

Only after v_info has been fully computed is

R_c = v_info / c_edge

evaluated.

---

## 16. Unit-compatibility gate

A direct physical comparison v_info ?≈ c_edge is permitted only if the implementation documents a common mapping between:

- Paper 31 lattice distance d;
- Paper 31 evolution time t=P Delta t;
- Paper 22 graph spatial units;
- Paper 22 propagation-time units.

If this mapping is not established independently of the observed Paper 31 velocity, the result must be reported only as a dimensionless circuit velocity and no claim of equality with c_edge is allowed.

No scale factor may be fitted using the requirement v_info=c_edge.

---

## 17. Cross-sector compatibility endpoint

If and only if the unit-compatibility gate is passed, define the preregistered compatibility condition

|v_info/c_edge - 1| <= 0.10.

This 10% tolerance is a finite-size numerical compatibility window, not an observational gravitational-wave bound.

Paper 31 is classified as CROSS-SECTOR COMPATIBLE only if:

1. the information front is internally coherent;
2. the unit-compatibility gate passes;
3. the 10% condition is satisfied.

---

## 18. Circuit-support diagnostic

For each depth P, record the largest distance d satisfying

Delta_A(d,t) > 1e-12.

This diagnostic estimates the maximum finite-depth support velocity v_circuit.

It is descriptive and must not be substituted for v_info.

---

## 19. Numerical safeguards

Frozen controls:

- state normalization tolerance: 1e-12;
- reduced-density trace tolerance: 1e-12;
- Hermiticity tolerance: 1e-12;
- positivity tolerance: 1e-12;
- trace-distance numerical zero diagnostic: 1e-14;
- circuit-support diagnostic threshold: 1e-12.

No data point may be silently deleted.

Undefined crossings, failed density checks, zero-signal time points, or failed fits must be retained and reported.

---

## 20. Recovery controls

Before Paper 31 results are accepted:

1. the unperturbed forward P=6 evolution must reproduce the forward Strang branch generated by the frozen Papers 28–30 state-preparation implementation; equality with the Paper 30 time-reversal mixture is not required;
2. for perturbations outside the exact finite-depth support, the observed signal must be compatible with floating-point zero;
3. all reduced density matrices must have unit trace within tolerance;
4. all reduced density matrices must be Hermitian and positive semidefinite within tolerance;
5. the trace distance must satisfy 0 <= Delta_A(d,t) <= 1 within numerical tolerance.

---

## 21. Primary interpretation boundary

If the Paper 31 front is internally coherent, the allowed conclusion is:

> A local perturbation in the frozen BuP information-sector dynamics exhibits a reproducible finite-speed information front distinct from the formal circuit-support boundary.

If the unit-compatibility gate also passes and the cross-sector compatibility criterion is satisfied, the additional allowed conclusion is:

> The independently extracted information-front velocity is compatible, within the preregistered finite-size tolerance, with the calibrated true-edge propagation speed of Paper 22.

Not allowed solely from Paper 31:

- derivation of the physical speed of light from first principles;
- proof of Lorentz invariance;
- proof of continuum causality;
- proof of a black-hole horizon;
- proof of a universal Lieb-Robinson velocity;
- proof that information and gravitational waves are literally the same excitation;
- thermodynamic-limit universality.

---

## 22. Failure / null-result interpretation

A negative or null result is retained as scientific information.

Examples:

- if the front is not linear in time, the ballistic-front hypothesis is not supported;
- if the three threshold velocities disagree, a unique front velocity is not supported;
- if the measured front saturates v_circuit, the signal may primarily reflect the circuit architecture;
- if v_info differs from c_edge, cross-sector velocity universality is not supported;
- if the unit mapping cannot be established, no direct v_info/c_edge claim is made.

No endpoint is redefined after unblinding.

---

## 23. Required frozen outputs

The implementation must preserve:

- full Delta_A(d,t) matrix;
- normalized profiles F(d,t);
- Frobenius-distance matrix;
- secondary W-sector differences;
- secondary modular-sector differences;
- all threshold crossings;
- all v_q, b_q, R_q^2;
- v_info;
- sigma_v;
- v_circuit;
- recovery diagnostics;
- provenance SHA256;
- frozen manifest;
- final classification.

---

## 24. Freeze requirements before first scientific run

Before any Paper 31 result is inspected:

1. this preregistration must be finalized;
2. the preregistration must be committed;
3. the implementation must be syntax-checked;
4. the implementation SHA256 must be recorded;
5. the implementation must be committed;
6. the Git working tree must be clean.

Any protocol correction made before scientific unblinding must be documented as a pre-unblinding amendment.

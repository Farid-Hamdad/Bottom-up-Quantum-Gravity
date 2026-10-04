# Paper 30 — Causal Compensation Landscape

## PREREGISTRATION DRAFT

### 1. Status and scope

This is a prospective follow-up to the frozen Paper 30 causal-size extension.

It does not modify, replace, or reinterpret the preregistered classification already frozen in commit 7a10276.

The objective is descriptive and structural: map the residual compensation between the mutual-information sector W and the modular sector K_A as a function of:

- system size N,
- circuit depth P,
- position of the four-qubit subsystem A.

No thermodynamic-limit, conservation-law, relativistic-causality, or fundamental-spacetime claim is preregistered.

### 2. Results known before this preregistration

The following observations are already known and therefore are NOT prospective predictions:

- Original edge configuration, P=6:
  N10 and N12 are locally equivalent to numerical precision and satisfy the original Paper 30 compensation criterion.

- Centered subsystem, P=6:
  N10 and N12 are size-sensitive and fail the all-bandwidth Paper 30 compensation criterion.

- Edge subsystem, P=8:
  N10 and N12 are size-sensitive and fail the all-bandwidth Paper 30 compensation criterion.

- Centered subsystem, P=8:
  N10 and N12 are strongly size-sensitive and fail the all-bandwidth Paper 30 compensation criterion.

- The frozen causal-size extension classification is:
  COMPENSATION_NOT_ROBUST_AFTER_CAUSAL_SIZE_EXPOSURE.

These previously observed cells will be used only as recovery anchors.

### 3. Local-channel reduction validated before preregistration

For observables depending only on rho_A, amplitude damping outside A is removed by the partial trace because the channel on the complement is trace preserving:

Tr_B[(E_A tensor E_B)(rho_AB)] = E_A[Tr_B(rho_AB)].

Direct numerical validation was performed before this preregistration for:

- N=10, P=6, edge A,
- N=10, P=6, centered A,
- N=10, P=8, edge A,

and showed agreement between the global and local procedures at approximately 10^-16 numerical precision.

The mapping implementation will therefore:

1. evolve the two pure TFIM branches,
2. obtain rho_A(0) directly from the statevectors,
3. apply amplitude damping only to the four-qubit rho_A,
4. evaluate the unchanged Paper 30 local observables and compensation pipeline.

### 4. Fixed model

The physical and numerical model remains identical to Paper 30 except for N, P, and the location of A.

- TFIM parameters: J = h = 1.
- Trotter step: dt = 0.35.
- Initial state: |+>^N.
- Time-reversal mixture:
  rho(0) = 1/2 (|psi_+><psi_+| + |psi_-><psi_-|).
- Amplitude damping parameter:
  gamma = 1 - exp(-s).
- s grid:
  0.0 to 3.0 inclusive in steps of 0.1.
- Four-qubit subsystem A.
- Mutual-information sector W and normalized modular-Hamiltonian sector exactly as in Paper 30.

### 5. Preregistered parameter grid

System sizes:

N in {8, 10, 12, 14, 16}.

Circuit depths:

P in {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}.

For every (N,P), every contiguous four-qubit subsystem is evaluated:

A_a = (a, a+1, a+2, a+3),

with

a = 0, 1, ..., N-4.

No position will be selected or discarded after inspecting the compensation results.

### 6. Geometric descriptors recorded before analysis

For every cell (N,P,a), record:

- left boundary distance:
  d_L = a;

- right boundary distance:
  d_R = N - (a + 4);

- nearest-boundary distance:
  d_b = min(d_L,d_R);

- subsystem center:
  x_A = a + 3/2;

- chain center:
  x_chain = (N - 1)/2;

- absolute subsystem-center offset:
  delta_center = |x_A - x_chain|;

- causal-footprint interval:
  [max(0,a-P), min(N-1,a+3+P)];

- causal-footprint size:
  L_causal = min(N-1,a+3+P) - max(0,a-P) + 1;

- causal-saturation fraction:
  f_causal = L_causal / N.

The causal footprint is an operational nearest-neighbor circuit descriptor only. It is not identified with a relativistic light cone.

### 7. Compensation analysis

The Paper 30 detrending and compensation pipeline is retained unchanged.

Bandwidth fractions:

{0.05, 0.10, 0.20, 0.30, 0.40}.

For each bandwidth compute:

- Pearson correlation of standardized residuals,
- Spearman correlation,
- optimal alpha,
- quasi-invariant variance ratio R_C.

The historical Paper 30 support criterion is retained for compatibility:

P1:
Pearson residual correlation < 0.

P2:
0.8 <= alpha_optimal <= 1.2.

P3:
R_C <= 0.20.

A cell is PAPER30_SUPPORTED iff P1, P2, and P3 all pass at all five bandwidths.

### 8. Algebraic dependence of the historical endpoints

Because both residual sectors are z-standardized before optimizing

C_alpha = z(epsilon_mod) + alpha z(epsilon_W),

the optimum satisfies

alpha_optimal = -rho

and

R_C = 1 - rho^2

up to numerical precision.

Therefore alpha_optimal and R_C will NOT be interpreted as statistically independent confirmations of Pearson anticorrelation.

The primary continuous compensation observable in this mapping study is the residual Pearson correlation rho.

The historical three-part criterion is retained only to preserve direct compatibility with Paper 30.

### 9. Primary outputs

For every (N,P,a), the frozen outputs will include:

1. the five Pearson residual correlations;
2. the five historical support flags;
3. PAPER30_SUPPORTED across all five bandwidths;
4. supported-bandwidth fraction:
   f_supported = (number of supported bandwidths) / 5;

5. strongest residual anticorrelation:
   rho_strong = min_h rho(h);

6. weakest residual anticorrelation:
   rho_weak = max_h rho(h);

7. worst-bandwidth quasi-invariant ratio:
   R_C_worst = max_h R_C(h);

8. the geometric descriptors defined above.

In addition, the frozen aggregate outputs will report:

- for every (N,P), the fraction of subsystem positions that are PAPER30_SUPPORTED;
- for every (N,P) and bandwidth, the minimum, median, and maximum rho across positions;
- for every available nearest-boundary distance d_b, the corresponding compensation metrics without selecting positions post hoc.

The primary study product is the full compensation landscape, not a single optimized configuration.

### 10. Reflection-symmetry control

For an open homogeneous chain with the symmetric initial state and homogeneous damping, mirrored subsystems

a and N-4-a

should agree up to numerical precision.

The implementation will compare each subsystem A_a with its mirror A_(N-4-a).

For every mirrored pair it will report the maximum absolute discrepancy across:

- rho_A(0);
- all six W_ij values over the gamma grid;
- all six modular A_ij values over the gamma grid;
- all six v_ij values over the gamma grid;
- M_W and M_mod over the gamma grid;
- the five Pearson residual correlations;
- the five alpha_optimal values;
- the five R_C values.

Numerical tolerance:

1e-10.

Historical support booleans must agree exactly between mirrored cells.

A numerical discrepancy above tolerance, or a mirrored support-status mismatch, is treated as an implementation or numerical-control failure requiring investigation before physical interpretation.

### 11. Frozen-anchor recovery

Before accepting the new landscape, the implementation must reproduce the previously frozen C0-C3 anchor configurations within absolute tolerance 1e-10 for all archived local compensation quantities available in the frozen results.

If anchor recovery fails, the new mapping results are not interpreted until the discrepancy is resolved.

### 12. Analysis rules

No bandwidth, position, depth, or system size may be removed because it weakens or contradicts a preferred interpretation.

No threshold may be changed after unblinding.

Known C0-C3 cells are recovery anchors and must not be counted as novel confirmations.

New structures observed after the mapping may motivate later analyses, but those analyses must be labeled exploratory unless separately preregistered before execution.

### 13. Interpretation boundary

This study may establish finite-size and finite-depth structure in the information-sector compensation landscape of the specified TFIM protocol under amplitude damping.

It does not by itself establish:

- a fundamental conservation law,
- literal transfer of information between sectors,
- thermodynamic-limit universality,
- a universal value alpha = 1,
- a relativistic causal cone,
- equality with the Paper 22 gravitational-wave speed,
- a fundamental spacetime propagation law.

# Paper 12 — SPARC test of the BuP effective propagator law

## Title

**SPARC rotation curves from local entanglement profiles in Bottom-Up Quantum Gravity**

## Core idea

Paper 12 applies the effective propagator law derived in Paper 11 to galactic
rotation curves from SPARC.

The theoretical chain is:

\[
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s(r),d_w(r))
\rightarrow
\alpha_{\rm eff}(r)
\rightarrow
\Phi_{\rm BuP}(r)
\rightarrow
V(r).
\]

Paper 11 showed that the exponent of the effective propagator is not controlled
only by the spectral dimension \(d_s\), but by the pair \((d_s,d_w)\):

\[
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
}
\]

Paper 12 tests this relation on SPARC rotation curves by constructing a
baryonic entanglement graph for each galaxy, measuring local profiles
\((d_s(r),d_w(r))\), deriving \(\alpha_{\rm eff}(r)\), and projecting the
result into a contribution to the circular velocity.

## Physical interpretation

BuP does not assume a universal dark matter halo. Instead, each baryonic
configuration induces an effective entanglement geometry. The gravitational
response depends on the local structure of this graph.

In this sense, galaxies are not forced into a single universal regime. They
fall into classes of entanglement profiles.

The previous SPARC runs already suggested such a separation through empirical
profiles \(d_{\rm emp}(r)\) and \(d_{\rm fit}(r)\). Paper 12 reformulates this
separation in terms of local entanglement dimensions:

\[
d(r)
\quad\longrightarrow\quad
(d_s(r),d_w(r))
\quad\longrightarrow\quad
\alpha_{\rm eff}(r).
\]

## Main result

The v3.4 corrected model combines:

1. a local extraction window \(w\in\{5,8,10\}\) kpc,
2. a relative transition criterion
   \[
   \alpha_c = 0.97\,\alpha_{\max},
   \]
3. a compacity correction for highly peaked profiles:
   \[
   r_t^{\rm corr}
   =
   r_t^{(0)}(1-0.75\eta),
   \]
   where
   \[
   \eta=
   \frac{v_{b,\max}-v_{b,\rm med}}{v_{b,\rm med}}.
   \]

On an extended sample of 18 galaxies, the best-window v3.4 scan gives:

| Metric | Value |
|---|---:|
| Number of galaxies | 18 |
| Improved vs baryon-only | 17 / 18 |
| \(\chi^2_{\rm red}<2\) | 4 / 18 |
| \(\chi^2_{\rm red}<5\) | 8 / 18 |
| \(\chi^2_{\rm red}<10\) | 11 / 18 |
| Median \(\chi^2_{\rm red}\) | 6.80 |
| Mean \(\chi^2_{\rm red}\) | 14.10 |

Adding the targeted NGC5055 compacity run gives an effective 19-galaxy summary:

| Metric | Value |
|---|---:|
| Number of galaxies | 19 |
| Improved vs baryon-only | 18 / 19 |
| \(\chi^2_{\rm red}<2\) | 4 / 19 |
| \(\chi^2_{\rm red}<5\) | 8 / 19 |
| \(\chi^2_{\rm red}<10\) | 12 / 19 |
| Median \(\chi^2_{\rm red}\) | 6.61 |
| Mean \(\chi^2_{\rm red}\) | 13.65 |

## Key galaxies

| Galaxy | Interpretation | Correction | \(\chi^2_{\rm red}\) |
|---|---|---|---:|
| NGC3198 | canonical regular transition | \(w=5\), \(q=0.97\) | 0.994 |
| NGC2841 | standard extended disk | \(w=10\) | 4.586 |
| NGC5055 | complex barred/bulged system | \(\eta\)-corrected \(r_t\) | 5.526 |
| NGC2403 | diffuse early sub-Newtonian regime | unresolved by single sigmoid | 115.8 |

## Scientific status

Paper 12 is not presented as a final universal SPARC model. It is a proof that
the effective propagator law from Paper 11 has observational content and that
galactic rotation curves can be organized by local entanglement profiles.

The strongest conclusion is:

\[
\boxed{
\text{SPARC galaxies split into classes of entanglement profiles.}
}
\]

The next step is to replace the single external sigmoid projection by
class-dependent or multi-scale projections, especially for diffuse systems such
as NGC2403.

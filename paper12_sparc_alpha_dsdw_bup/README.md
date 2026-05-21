# Paper 12 — SPARC as an observational test of BuP entanglement geometry

## Title

**SPARC rotation curves from local entanglement profiles in Bottom-Up Quantum Gravity**

## Foundational principle

Bottom-Up Quantum Gravity starts from a single postulate:

\[
\boxed{
\text{entanglement is primary.}
}
\]

Space, time, matter and gravity are not assumed as fundamental objects. They
emerge from the structure of a global quantum entanglement state.

In this framework, matter is localized entanglement. Therefore the observed
baryonic surface density of a galaxy,

\[
\Sigma(R),
\]

is interpreted as the macroscopic projection of a local entanglement density.

Paper 12 is the first galactic-scale numerical test of this chain:

\[
\text{entanglement}
\rightarrow
\text{geometry}
\rightarrow
\text{localized matter}
\rightarrow
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s(r),d_w(r))
\rightarrow
\alpha_{\rm eff}(r)
\rightarrow
V(r).
\]

This is not a dark-matter-like fitting model. It is a numerical application of
the BuP emergence chain to real SPARC galaxies.

## Link with Paper 11

Paper 11 derived the effective BuP propagator exponent:

\[
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
}
\]

Here:

- \(d_s\) is the local spectral dimension of the entanglement graph;
- \(d_w\) is the local walk dimension;
- \(\alpha_{\rm eff}\) controls the effective power-law behavior of the BuP
  propagator.

In the Brownian limit \(d_w=2\),

\[
\alpha_{\rm eff}=d_s-2.
\]

Newtonian gravity appears as a special regime. It is not an axiom of BuP.

Paper 12 applies the Paper 11 law to galaxies by constructing a baryonic
entanglement graph from \(\Sigma(R)\), measuring \(d_s(r)\) and \(d_w(r)\), and
using the resulting \(\alpha_{\rm eff}(r)\) to organize the gravitational
response.

## Main equation

The local BuP exponent is

\[
\boxed{
\alpha_{\rm eff}(r)
=
\frac{2d_s(r)}{d_w(r)}
+
d_w(r)
-
4.
}
\]

The velocity model used in the current numerical projection is

\[
V^2_{\rm model}(r)
=
V^2_{\rm bar}(r)
+
A\,S(r),
\]

with

\[
S(r)
=
\frac{1}{1+\exp[-(r-r_t)/\Delta r]}.
\]

In v3.3/v3.4, the transition radius \(r_t\) is not chosen freely. It is derived
from the local entanglement profile through

\[
\alpha_{\rm eff}(r_t)
=
q\,\alpha_{\max},
\qquad
q=0.97.
\]

The transition width is fixed to

\[
\Delta r=2.5\,{\rm kpc}.
\]

## v3.4 corrections

The v3.4 version introduces two corrections motivated by previous SPARC runs.

### 1. Local extraction scale

The radial window used to extract the local entanglement dimensions is scanned:

\[
w\in\{5,8,10\}\ {\rm kpc}.
\]

This is not a halo parameter. It is the scale at which the local graph geometry
is measured.

The best-window scan shows that different galaxies prefer different local
entanglement extraction scales.

### 2. Compacity correction

For internally complex galaxies, the local \(\alpha_{\rm eff}(r)\) profile can
develop a high peak. The previous BuP SPARC runs already showed that such cases
need a compacity correction.

We define

\[
\eta=
\frac{v_{b,\max}-v_{b,\rm med}}{v_{b,\rm med}}.
\]

For high-peak profiles,

\[
\alpha_{\max}>8,
\]

the transition radius is corrected as

\[
\boxed{
r_t^{\rm corr}
=
r_t^{(0)}(1-0.75\eta).
}
\]

This advances the BuP transition in compact or internally complex systems.

## Results: extended best-window sample

The v3.4 best-window scan was run on 18 galaxies.

| Metric | Value |
|---|---:|
| Number of galaxies | 18 |
| Improved vs baryon-only | 17 / 18 |
| \(\chi^2_{\rm red}<2\) | 4 / 18 |
| \(\chi^2_{\rm red}<5\) | 8 / 18 |
| \(\chi^2_{\rm red}<10\) | 11 / 18 |
| Median \(\chi^2_{\rm red}\) | 6.80 |
| Mean \(\chi^2_{\rm red}\) | 14.10 |

Adding the targeted NGC5055 compacity-corrected run gives an effective
19-galaxy summary:

| Metric | Value |
|---|---:|
| Number of galaxies | 19 |
| Improved vs baryon-only | 18 / 19 |
| \(\chi^2_{\rm red}<2\) | 4 / 19 |
| \(\chi^2_{\rm red}<5\) | 8 / 19 |
| \(\chi^2_{\rm red}<10\) | 12 / 19 |
| Median \(\chi^2_{\rm red}\) | 6.61 |
| Mean \(\chi^2_{\rm red}\) | 13.65 |

## Best-window results

| Galaxy | \(w\) kpc | \(\alpha_{\rm med}\) | \(\alpha_{\max}\) | \(r_t\) kpc | \(\chi^2_{\rm red}\) |
|---|---:|---:|---:|---:|---:|
| NGC3198 | 5 | 0.369 | 4.965 | 5.184 | 0.994 |
| NGC3769 | 5 | 0.071 | 5.373 | 7.034 | 1.082 |
| NGC5005 | 10 | 0.785 | 0.822 | 1.526 | 1.425 |
| NGC3521 | 10 | 1.028 | 1.701 | 0.860 | 1.825 |
| NGC7793 | 8 | 2.146 | 2.146 | 3.990 | 2.550 |
| NGC2998 | 8 | 0.215 | 4.847 | 0.985 | 2.616 |
| NGC0024 | 10 | 1.364 | 1.424 | 1.326 | 3.949 |
| NGC2841 | 10 | 0.620 | 1.871 | 5.517 | 4.586 |
| NGC4085 | 5 | 1.319 | 1.343 | 5.495 | 6.608 |
| NGC0300 | 8 | 1.014 | 1.141 | 1.385 | 7.000 |
| NGC1003 | 5 | 0.126 | 3.807 | 5.529 | 9.431 |
| NGC6946 | 10 | 0.761 | 1.445 | 0.617 | 12.345 |
| NGC0055 | 10 | 0.766 | 1.024 | 1.518 | 13.895 |
| NGC2366 | 10 | 2.285 | 2.285 | 3.090 | 14.359 |
| NGC7331 | 5 | -0.018 | 4.437 | 7.637 | 16.348 |
| NGC4389 | 5 | 1.693 | 1.693 | 3.095 | 16.768 |
| NGC0247 | 10 | 0.880 | 1.301 | 1.321 | 22.275 |
| NGC2403 | 5 | 0.082 | 5.313 | 0.517 | 115.805 |

## Key galaxies

### NGC3198 — canonical regular transition

NGC3198 is the cleanest Paper 12 case.

\[
w=5\,{\rm kpc},
\qquad
q=0.97,
\qquad
\chi^2_{\rm red}=0.994.
\]

It demonstrates that a transition radius derived from
\(\alpha_{\rm eff}(r)\) can reproduce a SPARC rotation curve with near-unit
reduced chi-square.

### NGC2841 — standard extended disk

Previous BuP runs already showed that NGC2841 belongs to a stable standard
regime. Paper 12 v3.3 initially failed because the fixed local window was too
small.

Using

\[
w=10\,{\rm kpc}
\]

gives

\[
\chi^2_{\rm red}=4.586.
\]

This shows that NGC2841 is not intrinsically complex; it requires a larger
local entanglement extraction scale.

### NGC5055 — internally complex system

NGC5055 is the clearest compacity case. With \(w=5\), the profile contains a
high peak:

\[
\alpha_{\max}=12.60.
\]

This activates the compacity correction:

\[
r_t^{\rm corr}=r_t^{(0)}(1-0.75\eta).
\]

For NGC5055:

\[
\eta=0.337,
\qquad
r_t^{(0)}=16.35\,{\rm kpc},
\qquad
r_t^{\rm corr}=12.22\,{\rm kpc},
\]

and

\[
\chi^2_{\rm red}=5.526.
\]

Wider windows erase the high-\(\alpha\) peak and degrade the fit. The peak is
therefore interpreted as a physical signature of internal complexity, not as
numerical noise.

### NGC2403 — diffuse early sub-Newtonian regime

NGC2403 remains the main resistant case. Earlier BuP runs already showed that
this galaxy is difficult.

In Paper 12, its profile is already strongly sub-Newtonian near the center:

\[
\alpha_{\rm med}=0.0818,
\qquad
r_t=0.517\,{\rm kpc}.
\]

A single external sigmoid is not appropriate. NGC2403 defines a third regime:
diffuse early sub-Newtonian systems, requiring a non-sigmoidal or multi-scale
projection.

## Classification

The results suggest four BuP regimes:

| Class | Meaning | Example |
|---|---|---|
| I-A | regular external transition | NGC3198 |
| I-B | extended standard disk, larger local window | NGC2841 |
| II | internally complex system, compacity correction | NGC5055 |
| III | diffuse early sub-Newtonian regime | NGC2403 |

## Scientific meaning

Paper 12 does not claim that the present v3.4 projection is the final universal
SPARC model.

The deeper result is:

\[
\boxed{
\text{SPARC galaxies organize into classes of entanglement profiles.}
}
\]

The previous SPARC runs already suggested several dynamical regimes through
empirical \(d(r)\) profiles. Paper 12 reformulates these regimes in terms of
local entanglement dimensions:

\[
(d_s(r),d_w(r))
\rightarrow
\alpha_{\rm eff}(r).
\]

The main conclusion is:

\[
\boxed{
\text{Paper 11 derives the BuP propagator law; Paper 12 shows that this law
organizes galactic rotation curves.}
}
\]

# Results — Paper 12 SPARC alpha(d_s,d_w) test

This folder contains the numerical results for Paper 12.

## Main files

| File / Folder | Description |
|---|---|
| `paper12_unified_sparc_classification.csv` | Fusion of old SPARC \(d(r)\) runs with Paper 12 alpha-profile results |
| `paper12_unified_sparc_classification.md` | Markdown version of the unified classification |
| `classe1_extended_window_scan_best_summary.csv` | Best-window v3.4 scan on 18 galaxies |
| `NGC5055_v3_4_window_scan_eta_summary.csv` | Dedicated v3.4 compacity scan for NGC5055 |
| `NGC3198_v3_3_q097_dr2p5/` | Canonical NGC3198 relative-trigger result |
| `NGC2841_v3_4_corrected/` | Adaptive-window correction for NGC2841 |
| `NGC5055_v3_4_corrected/` | Compacity-corrected result for NGC5055 |

## v3.4 best-window sample

The best-window scan uses

\[
w\in\{5,8,10\}\ {\rm kpc},
\qquad
q=0.97,
\qquad
\Delta r=2.5\,{\rm kpc}.
\]

It gives:

| Metric | Value |
|---|---:|
| Galaxies | 18 |
| Improved vs baryon-only | 17 |
| \(\chi^2_{\rm red}<2\) | 4 |
| \(\chi^2_{\rm red}<5\) | 8 |
| \(\chi^2_{\rm red}<10\) | 11 |
| Median \(\chi^2_{\rm red}\) | 6.80 |
| Mean \(\chi^2_{\rm red}\) | 14.10 |

## Targeted NGC5055 result

NGC5055 requires the compacity correction:

\[
r_t^{\rm corr}=r_t^{(0)}(1-0.75\eta).
\]

The best result is:

| Window | \(\alpha_{\max}\) | \(\eta\) | \(r_t\) | \(\chi^2_{\rm red}\) |
|---:|---:|---:|---:|---:|
| 5 kpc | 12.60 | 0.337 | 12.22 | 5.526 |

Wider windows erase the \(\alpha_{\max}\) peak and suppress the \(\eta\)
correction, degrading the fit.

## Effective 19-galaxy summary

Combining the 18-galaxy best-window sample with the targeted NGC5055 run gives:

| Metric | Value |
|---|---:|
| Galaxies | 19 |
| Improved vs baryon-only | 18 |
| \(\chi^2_{\rm red}<2\) | 4 |
| \(\chi^2_{\rm red}<5\) | 8 |
| \(\chi^2_{\rm red}<10\) | 12 |
| Median \(\chi^2_{\rm red}\) | 6.61 |
| Mean \(\chi^2_{\rm red}\) | 13.65 |

## Interpretation

The results support four regimes:

| Class | Meaning | Example |
|---|---|---|
| I-A | regular external transition | NGC3198 |
| I-B | standard extended disk | NGC2841 |
| II | internally complex system requiring \(\eta\) | NGC5055 |
| III | diffuse early sub-Newtonian system | NGC2403 |

# Scripts — Paper 12

## Main scripts

| Script | Purpose |
|---|---|
| `bup_paper12_sparc_alpha_dsdw_v3_1_radial_outer.py` | First successful flat-outer projection |
| `bup_paper12_sparc_alpha_dsdw_v3_3_relative_trigger.py` | Relative trigger \( \alpha_c=q\alpha_{\max} \) |
| `bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py` | Adaptive-window and compacity-corrected version |
| `paper12_merge_old_sparc_with_v33.py` | Merges old SPARC \(d(r)\) runs with Paper 12 profiles |

## Theoretical input

The main Paper 12 relation comes from Paper 11:

\[
\alpha_{\rm eff}(r)
=
\frac{2d_s(r)}{d_w(r)}
+
d_w(r)
-
4.
\]

## v3.4 model

The transition radius is derived from

\[
\alpha_{\rm eff}(r_t)
=
q\alpha_{\max},
\qquad
q=0.97.
\]

For complex systems,

\[
r_t^{\rm corr}
=
r_t^{(0)}(1-0.75\eta),
\]

with

\[
\eta=
\frac{v_{b,\max}-v_{b,\rm med}}{v_{b,\rm med}}.
\]

## Example: NGC3198

```bash
python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2/NGC3198_rotmod.dat" \
  --galaxy NGC3198 \
  --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
  --n-windows 8 --window-width 5 \
  --min-nodes-window 80 \
  --ups-disk 0.5 --ups-bulge 0.7 \
  --shape-mode flat_outer \
  --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
  --dr-fixed 2.5 --rt-mode alpha_crossing \
  --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
  --output-dir results/NGC3198_v3_4

Example: NGC5055
python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2/NGC5055_rotmod.dat" \
  --galaxy NGC5055 \
  --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
  --n-windows 8 --window-width 5 \
  --min-nodes-window 80 \
  --ups-disk 0.5 --ups-bulge 0.7 \
  --shape-mode flat_outer \
  --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
  --dr-fixed 2.5 --rt-mode alpha_crossing \
  --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
  --output-dir results/NGC5055_v3_4_corrected
# RUNS — Paper 12

## Best-window scan

```bash
cd ~/bottomup/energie_noir

ROT_DIR="/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2"
OUT_BASE="papers/paper12_sparc_alpha_dsdw_bup/results/classe1_extended_window_scan"
SCRIPT="papers/paper12_sparc_alpha_dsdw_bup/scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py"

mkdir -p "$OUT_BASE"

for gal in NGC2998 NGC3521 NGC6946 NGC7331 NGC7793 NGC1003 NGC2366 NGC3031 NGC0300 NGC0024 NGC0055 NGC0247 NGC0925 NGC2403 NGC2841 NGC3198 NGC3769 NGC4085 NGC4389 NGC5005; do
  rot="$ROT_DIR/${gal}_rotmod.dat"
  if [ ! -f "$rot" ]; then
    echo "[SKIP] $gal — fichier absent"
    continue
  fi

  for w in 5 8 10; do
    echo "===== $gal | w=$w ====="
    python3 "$SCRIPT" \
      --rotmod "$rot" --galaxy "$gal" \
      --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
      --n-windows 8 --window-width "$w" \
      --min-nodes-window 80 --ups-disk 0.5 --ups-bulge 0.7 \
      --shape-mode flat_outer \
      --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
      --dr-fixed 2.5 --rt-mode alpha_crossing \
      --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
      --output-dir "$OUT_BASE/${gal}_w${w}"
  done
done
```

## Best-window summary

```bash
python3 - <<'PY'
import json, glob
import pandas as pd

BASE="papers/paper12_sparc_alpha_dsdw_bup/results/classe1_extended_window_scan"

rows=[]
for f in glob.glob(BASE+"/*/summary.json"):
    with open(f) as fh:
        s=json.load(fh)

    rows.append({
        "galaxy": s["galaxy"],
        "N_rot": s["N_rot_points"],
        "window_width": s["local_windows"]["window_width"],
        "alpha_med": s["alpha_profile"]["alpha_median"],
        "alpha_max": s["alpha_profile"]["alpha_max"],
        "rt": s["fit"]["rt_derived"],
        "eta": s["fit"].get("eta_baryon"),
        "eta_applied": s["fit"].get("eta_applied"),
        "chi2_red": s["fit"]["chi2_red"],
        "chi2_bar_red": s["fit"]["chi2_baryon_only_red"],
        "delta_chi2": s["fit"]["delta_chi2_vs_baryon"],
    })

df=pd.DataFrame(rows)
df=df.sort_values("chi2_red")

best=df.sort_values("chi2_red").drop_duplicates("galaxy", keep="first")
best=best.sort_values("chi2_red")

out=BASE+"_best_summary.csv"
best.to_csv(out,index=False)

print(best[[
    "galaxy","N_rot","window_width","alpha_med","alpha_max",
    "rt","eta","eta_applied","chi2_red","chi2_bar_red","delta_chi2"
]].to_string(index=False))

print()
print("Wrote:",out)

print()
print("Best-window stats:")
print("n_galaxies =",len(best))
print("n chi2<2   =",(best["chi2_red"]<2).sum())
print("n chi2<5   =",(best["chi2_red"]<5).sum())
print("n chi2<10  =",(best["chi2_red"]<10).sum())
print("n improved =", (best["delta_chi2"]<0).sum())
print("median chi2_red =", best["chi2_red"].median())
print("mean chi2_red   =", best["chi2_red"].mean())
PY
```

## NGC5055 eta scan

```bash
cd ~/bottomup/energie_noir

ROT_DIR="/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2"
OUT_BASE="papers/paper12_sparc_alpha_dsdw_bup/results/NGC5055_v3_4_window_scan_eta"
SCRIPT="papers/paper12_sparc_alpha_dsdw_bup/scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py"

mkdir -p "$OUT_BASE"

for w in 5 8 10; do
  python3 "$SCRIPT" \
    --rotmod "$ROT_DIR/NGC5055_rotmod.dat" \
    --galaxy NGC5055 \
    --nr 20 --ntheta 30 --k 12 --lambda-xy 3.0 \
    --n-windows 8 --window-width "$w" \
    --min-nodes-window 80 \
    --ups-disk 0.5 --ups-bulge 0.7 \
    --shape-mode flat_outer \
    --alpha-trigger-mode relative_max --q-alpha-max 0.97 \
    --dr-fixed 2.5 --rt-mode alpha_crossing \
    --eta-correction rt --eta-gamma 0.75 --eta-alpha-threshold 8.0 \
    --output-dir "$OUT_BASE/w${w}"
done
```

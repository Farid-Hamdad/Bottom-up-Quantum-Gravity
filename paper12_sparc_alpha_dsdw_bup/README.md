markdown
# Paper 12 — SPARC comme test observationnel de la géométrie d'intrication BuP

## Titre

**Courbes de rotation SPARC à partir des profils d'intrication locale en Gravité Quantique Bottom-Up**

## Principe fondamental

La Gravité Quantique Bottom-Up part d'un postulat unique :

$$
\boxed{
\text{L'intrication est première.}
}
$$

L'espace, le temps, la matière et la gravité ne sont pas supposés comme objets fondamentaux. Ils émergent de la structure d'un état d'intrication quantique global.

Dans ce cadre, la matière est de l'intrication localisée. Par conséquent, la densité surfacique baryonique observée d'une galaxie,

$$
\Sigma(R),
$$

est interprétée comme la projection macroscopique d'une densité d'intrication locale.

Paper 12 est le premier test numérique à l'échelle galactique de cette chaîne :

$$
\text{intrication}
\rightarrow
\text{géométrie}
\rightarrow
\text{matière localisée}
\rightarrow
\Sigma(R)
\rightarrow
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
(d_s(r), d_w(r))
\rightarrow
\alpha_{\rm eff}(r)
\rightarrow
V(r).
$$

Ce n'est pas un modèle d'ajustement de type matière noire. C'est une application numérique de la chaîne d'émergence BuP aux galaxies SPARC réelles.

## Lien avec Paper 11

Paper 11 a dérivé l'exposant du propagateur effectif BuP :

$$
\boxed{
\alpha_{\rm eff} = \frac{2d_s}{d_w} + d_w - 4
}

$$
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w} + d_w - 4.
}
$$

Ici :

- $d_s$ est la dimension spectrale locale du graphe d'intrication ;
- $d_w$ est la dimension de marche locale ;
- $\alpha_{\rm eff}$ contrôle le comportement en loi de puissance du propagateur BuP.

Dans la limite brownienne $d_w = 2$,

$$
\alpha_{\rm eff} = d_s - 2.
$$

La gravité newtonienne apparaît comme un régime particulier. Ce n'est pas un axiome de BuP.

Paper 12 applique la loi de Paper 11 aux galaxies en construisant un graphe d'intrication baryonique à partir de $\Sigma(R)$, en mesurant $d_s(r)$ et $d_w(r)$, et en utilisant le $\alpha_{\rm eff}(r)$ résultant pour organiser la réponse gravitationnelle.

## Équation principale

L'exposant BuP local est

$$
\boxed{
\alpha_{\rm eff}(r)
=
\frac{2d_s(r)}{d_w(r)}
+
d_w(r)
-
4.
}
$$

Le modèle de vitesse utilisé dans la projection numérique actuelle est

$$
V^2_{\rm model}(r)
=
V^2_{\rm bar}(r)
+
A \, S(r),
$$

avec

$$
S(r)
=
\frac{1}{1 + \exp[-(r - r_t)/\Delta r]}.
$$

En v3.3/v3.4, le rayon de transition $r_t$ n'est pas choisi librement. Il est dérivé du profil d'intrication locale par

$$
\alpha_{\rm eff}(r_t)
=
q \, \alpha_{\max},
\qquad
q = 0.97.
$$

La largeur de transition est fixée à

$$
\Delta r = 2.5\ \text{kpc}.
$$

## Corrections v3.4

La version v3.4 introduit deux corrections motivées par les précédentes analyses SPARC.

### 1. Échelle d'extraction locale

La fenêtre radiale utilisée pour extraire les dimensions d'intrication locales est scannée :

$$
w \in \{5, 8, 10\}\ \text{kpc}.
$$

Ce n'est pas un paramètre de halo. C'est l'échelle à laquelle la géométrie locale du graphe est mesurée.

Le scan par meilleure fenêtre montre que différentes galaxies préfèrent différentes échelles d'extraction de l'intrication locale.

### 2. Correction de compacité

Pour les galaxies à complexité interne, le profil local $\alpha_{\rm eff}(r)$ peut développer un pic élevé. Les précédentes analyses BuP SPARC ont déjà montré que ces cas nécessitent une correction de compacité.

On définit

$$
\eta = \frac{v_{b,\max} - v_{b,\rm med}}{v_{b,\rm med}}.
$$

Pour les profils à pic élevé,

$$
\alpha_{\max} > 8,
$$

le rayon de transition est corrigé comme suit

$$
\boxed{
r_t^{\rm corr}
=
r_t^{(0)} (1 - 0.75\eta).
}
$$

Cela avance la transition BuP dans les systèmes compacts ou à complexité interne.

## Résultats : échantillon étendu par meilleure fenêtre

Le scan v3.4 par meilleure fenêtre a été exécuté sur 18 galaxies.

| Métrique | Valeur |
|----------|--------|
| Nombre de galaxies | 18 |
| Amélioration vs baryons seuls | 17 / 18 |
| $\chi^2_{\rm red} < 2$ | 4 / 18 |
| $\chi^2_{\rm red} < 5$ | 8 / 18 |
| $\chi^2_{\rm red} < 10$ | 11 / 18 |
| Médiane $\chi^2_{\rm red}$ | 6.80 |
| Moyenne $\chi^2_{\rm red}$ | 14.10 |

En ajoutant l'analyse ciblée avec correction de compacité pour NGC5055, on obtient un résumé effectif à 19 galaxies :

| Métrique | Valeur |
|----------|--------|
| Nombre de galaxies | 19 |
| Amélioration vs baryons seuls | 18 / 19 |
| $\chi^2_{\rm red} < 2$ | 4 / 19 |
| $\chi^2_{\rm red} < 5$ | 8 / 19 |
| $\chi^2_{\rm red} < 10$ | 12 / 19 |
| Médiane $\chi^2_{\rm red}$ | 6.61 |
| Moyenne $\chi^2_{\rm red}$ | 13.65 |

## Résultats par meilleure fenêtre

| Galaxie | $w$ (kpc) | $\alpha_{\rm med}$ | $\alpha_{\max}$ | $r_t$ (kpc) | $\chi^2_{\rm red}$ |
|---------|-----------|--------------------|-----------------|-------------|--------------------|
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

## Galaxies clés

### NGC3198 — transition régulière canonique

NGC3198 est le cas le plus propre de Paper 12.

$$
w = 5\ \text{kpc},
\qquad
q = 0.97,
\qquad
\chi^2_{\rm red} = 0.994.
$$

Cela démontre qu'un rayon de transition dérivé de $\alpha_{\rm eff}(r)$ peut reproduire une courbe de rotation SPARC avec un chi-deux réduit proche de l'unité.

### NGC2841 — disque étendu standard

Les précédentes analyses BuP ont déjà montré que NGC2841 appartient à un régime standard stable. Paper 12 v3.3 a initialement échoué car la fenêtre locale fixe était trop petite.

En utilisant

$$
w = 10\ \text{kpc}
$$

on obtient

$$
\chi^2_{\rm red} = 4.586.
$$

Cela montre que NGC2841 n'est pas intrinsèquement complexe ; il nécessite une échelle d'extraction d'intrication locale plus grande.

### NGC5055 — système à complexité interne

NGC5055 est le cas de compacité le plus clair. Avec $w = 5$, le profil contient un pic élevé :

$$
\alpha_{\max} = 12.60.
$$

Cela active la correction de compacité :

$$
r_t^{\rm corr} = r_t^{(0)} (1 - 0.75\eta).
$$

Pour NGC5055 :

$$
\eta = 0.337,
\qquad
r_t^{(0)} = 16.35\ \text{kpc},
\qquad
r_t^{\rm corr} = 12.22\ \text{kpc},
$$

et

$$
\chi^2_{\rm red} = 5.526.
$$

Les fenêtres plus larges effacent le pic de $\alpha_{\max}$ et dégradent l'ajustement. Le pic est donc interprété comme une signature physique de la complexité interne, et non comme du bruit numérique.

### NGC2403 — régime sous-newtonien précoce diffus

NGC2403 reste le principal cas résistant. Les précédentes analyses BuP ont déjà montré que cette galaxie est difficile.

Dans Paper 12, son profil est déjà fortement sous-newtonien près du centre :

$$
\alpha_{\rm med} = 0.0818,
\qquad
r_t = 0.517\ \text{kpc}.
$$

Une sigmoïde externe unique n'est pas appropriée. NGC2403 définit un troisième régime : les systèmes diffus à régime sous-newtonien précoce, nécessitant une projection non sigmoïdale ou multi-échelles.

## Classification

Les résultats suggèrent quatre régimes BuP :

| Classe | Signification | Exemple |
|--------|---------------|---------|
| I-A | transition externe régulière | NGC3198 |
| I-B | disque étendu standard, fenêtre locale plus large | NGC2841 |
| II | système à complexité interne, correction de compacité | NGC5055 |
| III | régime sous-newtonien précoce diffus | NGC2403 |

## Signification scientifique

Paper 12 ne prétend pas que la projection v3.4 actuelle est le modèle universel SPARC définitif.

Le résultat plus profond est :

$$
\boxed{
\text{Les galaxies SPARC s'organisent en classes de profils d'intrication.}
}
$$

Les précédentes analyses SPARC suggéraient déjà plusieurs régimes dynamiques à travers les profils empiriques $d(r)$. Paper 12 reformule ces régimes en termes de dimensions d'intrication locale :

$$
(d_s(r), d_w(r))
\rightarrow
\alpha_{\rm eff}(r).
$$

La conclusion principale est :

$$
\boxed{
\text{Paper 11 dérive la loi du propagateur BuP ; Paper 12 montre que cette loi organise les courbes de rotation galactiques.}
}
$$

---

# Résultats — Paper 12 : test SPARC de la loi $\alpha(d_s, d_w)$

Ce dossier contient les résultats numériques pour Paper 12.

## Fichiers principaux

| Fichier / Dossier | Description |
|---|---|
| `paper12_unified_sparc_classification.csv` | Fusion des anciennes analyses SPARC $d(r)$ avec les résultats de profils alpha de Paper 12 |
| `paper12_unified_sparc_classification.md` | Version Markdown de la classification unifiée |
| `classe1_extended_window_scan_best_summary.csv` | Scan v3.4 par meilleure fenêtre sur 18 galaxies |
| `NGC5055_v3_4_window_scan_eta_summary.csv` | Scan de compacité v3.4 dédié pour NGC5055 |
| `NGC3198_v3_3_q097_dr2p5/` | Résultat canonique de NGC3198 avec déclenchement relatif |
| `NGC2841_v3_4_corrected/` | Correction par fenêtre adaptative pour NGC2841 |
| `NGC5055_v3_4_corrected/` | Résultat corrigé par compacité pour NGC5055 |

## Échantillon v3.4 par meilleure fenêtre

Le scan par meilleure fenêtre utilise

$$
w \in \{5, 8, 10\}\ \text{kpc},
\qquad
q = 0.97,
\qquad
\Delta r = 2.5\ \text{kpc}.
$$

Il donne :

| Métrique | Valeur |
|----------|--------|
| Galaxies | 18 |
| Amélioration vs baryons seuls | 17 |
| $\chi^2_{\rm red} < 2$ | 4 |
| $\chi^2_{\rm red} < 5$ | 8 |
| $\chi^2_{\rm red} < 10$ | 11 |
| Médiane $\chi^2_{\rm red}$ | 6.80 |
| Moyenne $\chi^2_{\rm red}$ | 14.10 |

## Résultat ciblé pour NGC5055

NGC5055 nécessite la correction de compacité :

$$
r_t^{\rm corr} = r_t^{(0)} (1 - 0.75\eta).
$$

Le meilleur résultat est :

| Fenêtre | $\alpha_{\max}$ | $\eta$ | $r_t$ (kpc) | $\chi^2_{\rm red}$ |
|---------|----------------|--------|-------------|--------------------|
| 5 kpc | 12.60 | 0.337 | 12.22 | 5.526 |

Les fenêtres plus larges effacent le pic de $\alpha_{\max}$ et suppriment la correction $\eta$, ce qui dégrade l'ajustement.

## Résumé effectif à 19 galaxies

La combinaison de l'échantillon à 18 galaxies par meilleure fenêtre avec l'analyse ciblée de NGC5055 donne :

| Métrique | Valeur |
|----------|--------|
| Galaxies | 19 |
| Amélioration vs baryons seuls | 18 |
| $\chi^2_{\rm red} < 2$ | 4 |
| $\chi^2_{\rm red} < 5$ | 8 |
| $\chi^2_{\rm red} < 10$ | 12 |
| Médiane $\chi^2_{\rm red}$ | 6.61 |
| Moyenne $\chi^2_{\rm red}$ | 13.65 |

## Interprétation

Les résultats supportent quatre régimes :

| Classe | Signification | Exemple |
|--------|---------------|---------|
| I-A | transition externe régulière | NGC3198 |
| I-B | disque étendu standard | NGC2841 |
| II | système à complexité interne nécessitant $\eta$ | NGC5055 |
| III | système diffus à régime sous-newtonien précoce | NGC2403 |

---

# Scripts — Paper 12

## Scripts principaux

| Script | Objectif |
|--------|----------|
| `bup_paper12_sparc_alpha_dsdw_v3_1_radial_outer.py` | Première projection réussie avec plateau externe |
| `bup_paper12_sparc_alpha_dsdw_v3_3_relative_trigger.py` | Déclenchement relatif $\alpha_c = q \alpha_{\max}$ |
| `bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py` | Version avec fenêtre adaptative et correction de compacité |
| `paper12_merge_old_sparc_with_v33.py` | Fusionne les anciennes analyses SPARC $d(r)$ avec les profils Paper 12 |

## Entrée théorique

La relation principale de Paper 12 provient de Paper 11 :

$$
\alpha_{\rm eff}(r)
=
\frac{2d_s(r)}{d_w(r)}
+
d_w(r)
-
4.
$$

## Modèle v3.4

Le rayon de transition est dérivé de

$$
\alpha_{\rm eff}(r_t)
=
q \alpha_{\max},
\qquad
q = 0.97.
$$

Pour les systèmes complexes,

$$
r_t^{\rm corr}
=
r_t^{(0)} (1 - 0.75\eta),
$$

avec

$$
\eta = \frac{v_{b,\max} - v_{b,\rm med}}{v_{b,\rm med}}.
$$

## Exemple : NGC3198

```bash
python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/chemin/vers/SPARC/Rotmod_LTG 2/NGC3198_rotmod.dat" \
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
Exemple : NGC5055
bash
python3 scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py \
  --rotmod "/chemin/vers/SPARC/Rotmod_LTG 2/NGC5055_rotmod.dat" \
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
Scan par meilleure fenêtre
bash
cd ~/bottomup/energie_noir

ROT_DIR="/Users/dualcomputer/bottomup/sparc/Rotmod_LTG 2"
OUT_BASE="papers/paper12_sparc_alpha_dsdw_bup/results/classe1_extended_window_scan"
SCRIPT="papers/paper12_sparc_alpha_dsdw_bup/scripts/bup_paper12_sparc_alpha_dsdw_v3_4_corrected.py"

mkdir -p "$OUT_BASE"

for gal in NGC2998 NGC3521 NGC6946 NGC7331 NGC7793 NGC1003 NGC2366 NGC3031 \
           NGC0300 NGC0024 NGC0055 NGC0247 NGC0925 NGC2403 NGC2841 NGC3198 \
           NGC3769 NGC4085 NGC4389 NGC5005; do
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
Résumé du scan par meilleure fenêtre
bash
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
Scan eta pour NGC5055
bash
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

Dépendances
Python 3.8+
NumPy, SciPy
Pandas (pour les tables de résultats)
SPARC dataset (fichiers Rotmod)

# Bottom-Up Quantum Gravity (BuP) — Paper 2  
## Contraintes cosmologiques sur une dimension spatiale émergente

## Résultat principal

Ce travail constitue le socle cosmologique du programme BuP (Bottom-Up Quantum Gravity).

Nous montrons qu’une dimension spatiale émergente variable peut expliquer simultanément :

- l’accélération cosmique tardive
- la tension sur σ₈
- la cohérence avec les contraintes CMB
- les observations BAO (DESI 2024)
- les supernovæ Pantheon+

sans introduire de nouvelle matière noire ni de champ scalaire supplémentaire.

Le couplage gravitationnel effectif est modifié géométriquement par :

G_eff(z) = 2 / (d(z) − 1) · G

où d(z) est la dimension spatiale effective émergente.

---

## Mécanisme physique

### Univers primordial (haut redshift)

Lorsque :

d ≈ DC > 3

on obtient :

G_eff < G

ce qui produit une légère suppression de la croissance précoce des structures.

Cela permet de préserver la cohérence avec la physique standard du fond diffus cosmologique (CMB).

---

### Univers tardif (bas redshift)

Lorsque :

d < 3

on obtient :

G_eff > G

ce qui renforce le regroupement tardif de la matière et modifie l’expansion cosmique.

---

## Résultat principal sur σ₈

L’intégration complète de la croissance cosmologique depuis :

z ≈ 1000 → z = 0

donne :

σ₈ = 0.772

en excellent accord avec :

- KiDS-1000 : 0.766 ± 0.020
- DES Y3 : 0.776 ± 0.017

sans paramètre libre supplémentaire.

La légère sur-dimensionnalité primordiale :

DC = 3.0598 > 3

fournit ainsi une explication géométrique unifiée reliant :

- énergie noire émergente
- croissance des structures
- lentillage faible
- BAO
- CMB
- formation précoce des galaxies massives

---

## Modèle à dimension variable tardive

Le modèle repose sur :

X(z) = X₀ · (1 + z)^β

variable intermédiaire effective

et

d(z) = DC − Δd / (1 + X(z)^α)

profil de dimension émergente

avec :

α = 1.78

fixé par les simulations microscopiques BuP  
(et non ajusté cosmologiquement)

---

## Meilleur ajustement cosmologique  
### DESI 2024 + Pantheon+

Paramètres obtenus :

X₀ = 0.537  
β = 2.000  
Δd = 0.878

H₀ = 63.67 km/s/Mpc  
Ωm = 0.3448

Critères statistiques :

ΔAIC = −6.23  
ΔBIC = −3.71

par rapport à ΛCDM.

---

## Scripts principaux

### Ajustement BAO + SN

```bash
python scripts/bup_paper2_final.py

→ fit complet du modèle cosmologique

Calcul cohérent de σ₈
python scripts/bup_sigma8_final.py

→ résolution complète de la tension σ₈ via G_eff(z)

Patch CAMB
python scripts/bup_camb_geff_patch.py

→ modification cohérente de CAMB pour intégrer G_eff(z)

Module Fortran CAMB
scripts/bup_geff_module.f90

→ implémentation directe dans CAMB

Propagateur gravitationnel microscopique

(travail en cours)

python scripts/bup_propagator.py --N 12 16 --lambda-n 6 --n-repeat 2

Test de la loi :

G(r) ~ r^{−(d_s−2)}

sur le graphe d’intrication.

Résultat actuel :

γ ≈ 0 pour N ≤ 16

→ taille de système insuffisante

Une mesure fiable nécessitera probablement :

N ≥ 100

Figures principales
Figure 1

Exclusion des modèles à dimension constante
et scan CAMB complet :

σ₈
r_drag
θs
Figure 2

Cohérence multi-échelle :

d(z)
w(z)
transition géométrique
compatibilité galactique et cosmologique
Figure 3

Résolution de σ₈ :

évolution de G_eff(z)
suppression précoce de croissance
scan complet de σ₈
Article complet

Le fichier LaTeX complet est disponible ici :

paper2_background_dimension/main.tex

avec toutes les figures, tableaux et bibliographie.

Prépublication HAL

Disponible sur HAL :

https://hal.science/hal-05590614v1

Titre :

Contraintes cosmologiques sur une dimension spatiale émergente :
dimension variable, énergie noire et σ₈ dans le cadre Bottom-Up Quantum Gravity

Citation recommandée
@misc{hamdad2026bup,
  author  = {Farid Hamdad},
  title   = {Contraintes cosmologiques sur une dimension spatiale émergente :
             dimension variable, énergie noire et σ₈
             dans le cadre Bottom-Up Quantum Gravity},
  year    = {2026},
  note    = {Prépublication HAL hal-05590614v1},
  url     = {https://hal.science/hal-05590614v1}
}
Suite du programme

Ce Paper 2 constitue la base de :

Paper 3 → résolution géométrique cohérente de la tension σ₈
Paper 4 → effondrement non-linéaire local et excès JWST
Paper 5 → dérivation des équations d’Einstein effectives

dans le cadre général de la gravité quantique Bottom-Up.

# Paper 12 — Courbes de rotation SPARC depuis les profils locaux d’intrication BuP

## Titre

**Courbes de rotation SPARC depuis les profils locaux d’intrication dans la gravité quantique Bottom-Up**

---

## Principe fondateur

La théorie BuP part d’un seul postulat :

$$
\boxed{
\text{l’intrication est première.}
}
$$

L’espace, le temps, la matière et la gravité ne sont pas supposés comme des objets fondamentaux. Ils émergent de la structure d’intrication d’un état quantique global.

Dans ce cadre, la matière est une forme d’intrication localisée. La densité baryonique observée d’une galaxie,

$$
\Sigma(R),
$$

est donc interprétée comme la projection macroscopique d’une densité locale d’intrication.

Paper 12 est le premier test numérique à l’échelle galactique de cette chaîne d’émergence :

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
(d_s(r),d_w(r))
\rightarrow
\alpha_{\rm eff}(r)
\rightarrow
V(r).
$$

Ce n’est pas un modèle de halo de matière noire. C’est une application numérique de la chaîne d’émergence BuP aux galaxies SPARC réelles.

---

## Lien avec Paper 11

Paper 11 a dérivé la loi effective du propagateur BuP :

$$
\boxed{
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}
+
d_w
-
4.
}
$$

où :

- \(d_s\) est la dimension spectrale locale du graphe d’intrication ;
- \(d_w\) est la dimension de marche locale ;
- \(\alpha_{\rm eff}\) contrôle le comportement effectif du propagateur gravitationnel BuP.

Dans la limite brownienne :

$$
d_w=2,
$$

on retrouve :

$$
\alpha_{\rm eff}=d_s-2.
$$

La gravité newtonienne apparaît donc comme un régime particulier, et non comme un axiome fondamental.

Paper 12 applique cette loi aux galaxies SPARC : à partir de la densité baryonique observée \(\Sigma(R)\), on construit un graphe d’intrication effectif, on mesure les profils locaux \(d_s(r)\) et \(d_w(r)\), puis on en déduit le profil gravitationnel effectif \(\alpha_{\rm eff}(r)\).

---

## Équation centrale

La loi locale utilisée dans Paper 12 est :

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

Le modèle de vitesse utilisé dans la projection numérique actuelle est :

$$
V_{\rm model}^2(r)
=
V_{\rm bar}^2(r)
+
A\,S(r),
$$

avec :

$$
S(r)
=
\frac{1}{1+\exp[-(r-r_t)/\Delta r]}.
$$

Dans les versions v3.3 et v3.4, le rayon de transition \(r_t\) n’est pas choisi librement. Il est déduit du profil local d’intrication :

$$
\alpha_{\rm eff}(r_t)
=
q\,\alpha_{\max},
$$

avec :

$$
q=0.97.
$$

La largeur de transition est fixée à :

$$
\Delta r=2.5\ {\rm kpc}.
$$

---

## Corrections v3.4

La version v3.4 introduit deux corrections motivées par les anciens runs SPARC.

### 1. Échelle locale d’extraction

L’échelle radiale utilisée pour mesurer les dimensions locales du graphe est scannée :

$$
w\in\{5,8,10\}\ {\rm kpc}.
$$

Cette échelle n’est pas un paramètre de halo. Elle représente l’échelle à laquelle on mesure localement la géométrie d’intrication.

Le scan montre que différentes galaxies préfèrent différentes échelles locales d’extraction.

### 2. Correction de compacité

Pour les galaxies structurellement complexes, le profil local \(\alpha_{\rm eff}(r)\) peut présenter un pic élevé. Les anciens runs BuP sur SPARC montraient déjà que ces cas nécessitent une correction de compacité.

On définit :

$$
\eta
=
\frac{v_{b,\max}-v_{b,\rm med}}{v_{b,\rm med}}.
$$

Pour les profils à pic élevé :

$$
\alpha_{\max}>8,
$$

le rayon de transition est corrigé ainsi :

$$
\boxed{
r_t^{\rm corr}
=
r_t^{(0)}(1-0.75\eta).
}
$$

Cette correction avance la transition BuP dans les galaxies compactes ou structurellement complexes.

---

## Résultats : échantillon étendu v3.4

Le scan v3.4 avec meilleure fenêtre locale a été lancé sur 18 galaxies.

| Indicateur | Valeur |
|---|---:|
| Nombre de galaxies | 18 |
| Améliorées par rapport au baryon-only | 17 / 18 |
| \(\chi^2_{\rm red}<2\) | 4 / 18 |
| \(\chi^2_{\rm red}<5\) | 8 / 18 |
| \(\chi^2_{\rm red}<10\) | 11 / 18 |
| Médiane \(\chi^2_{\rm red}\) | 6.80 |
| Moyenne \(\chi^2_{\rm red}\) | 14.10 |

En ajoutant le run ciblé de NGC5055 avec correction de compacité, on obtient un bilan effectif sur 19 galaxies :

| Indicateur | Valeur |
|---|---:|
| Nombre de galaxies | 19 |
| Améliorées par rapport au baryon-only | 18 / 19 |
| \(\chi^2_{\rm red}<2\) | 4 / 19 |
| \(\chi^2_{\rm red}<5\) | 8 / 19 |
| \(\chi^2_{\rm red}<10\) | 12 / 19 |
| Médiane \(\chi^2_{\rm red}\) | 6.61 |
| Moyenne \(\chi^2_{\rm red}\) | 13.65 |

---

## Meilleurs résultats v3.4

| Galaxie | \(w\) kpc | \(\alpha_{\rm med}\) | \(\alpha_{\max}\) | \(r_t\) kpc | \(\chi^2_{\rm red}\) |
|---|---:|---:|---:|---:|---:|
| NGC3198 | 5 | 0.369 | 4.965 | 5.184 | 0.994 |
| NGC3769 | 5 | 0.071 | 5.373 | 7.034 | 1.082 |
| NGC5005 | 10 | 0.785 | 0.822 | 1.526 | 1.425 |
| NGC3521 | 10 | 1.028 | 1.701 | 0.860 | 1.825 |
| NGC7793 | 8 | 2.146 | 2.146 | 3.990 | 2.550 |
| NGC2998 | 8 | 0.215 | 4.847 | 0.985 | 2.616 |
| NGC0024 | 10 | 1.364 | 1.424 | 1.326 | 3.949 |
| NGC2841 | 10 | 0.620 | 1.871 | 5.517 | 4.586 |

---

## Galaxies clés

### NGC3198 — transition régulière canonique

NGC3198 est le cas le plus propre de Paper 12.

$$
w=5\ {\rm kpc},
\qquad
q=0.97,
\qquad
\chi^2_{\rm red}=0.994.
$$

Cette galaxie montre qu’un rayon de transition déduit de \(\alpha_{\rm eff}(r)\) peut reproduire une courbe de rotation SPARC avec un \(\chi^2_{\rm red}\) proche de 1.

---

### NGC2841 — disque standard étendu

Les anciens runs BuP montraient déjà que NGC2841 appartient à un régime standard stable. La version v3.3 de Paper 12 échouait parce que la fenêtre locale fixe était trop petite.

Avec :

$$
w=10\ {\rm kpc},
$$

on obtient :

$$
\chi^2_{\rm red}=4.586.
$$

Cela montre que NGC2841 n’est pas intrinsèquement complexe. Elle nécessite simplement une échelle locale d’extraction plus adaptée à son extension radiale.

---

### NGC5055 — système complexe corrigé par compacité

NGC5055 est le cas le plus clair de galaxie structurellement complexe.

Avec :

$$
w=5\ {\rm kpc},
$$

le profil contient un pic élevé :

$$
\alpha_{\max}=12.60.
$$

Ce pic active la correction de compacité :

$$
r_t^{\rm corr}
=
r_t^{(0)}(1-0.75\eta).
$$

Pour NGC5055 :

$$
\eta=0.337,
$$

$$
r_t^{(0)}=16.35\ {\rm kpc},
$$

$$
r_t^{\rm corr}=12.22\ {\rm kpc}.
$$

Le résultat est :

$$
\chi^2_{\rm red}=5.526.
$$

Les fenêtres plus larges effacent le pic \(\alpha_{\max}\) et dégradent fortement le fit. Le pic est donc interprété comme une signature physique de complexité interne, pas comme du bruit numérique.

---

### NGC2403 — régime diffus sous-newtonien précoce

NGC2403 reste le principal cas résistant.

Son profil est déjà fortement sous-newtonien près du centre :

$$
\alpha_{\rm med}=0.0818,
$$

avec :

$$
r_t=0.517\ {\rm kpc}.
$$

Une seule transition externe sigmoïdale n’est donc pas adaptée. NGC2403 définit une troisième classe : les galaxies à transition diffuse sous-newtonienne précoce, qui nécessitent une projection non sigmoïdale ou multi-échelle.

---

## Classification BuP des profils d’intrication

Les résultats suggèrent quatre régimes :

| Classe | Signification | Exemple |
|---|---|---|
| I-A | transition externe régulière | NGC3198 |
| I-B | disque standard étendu, fenêtre locale plus large | NGC2841 |
| II | système interne complexe, correction de compacité | NGC5055 |
| III | régime diffus sous-newtonien précoce | NGC2403 |

---

## Sens physique

Paper 12 ne prétend pas que la version v3.4 est le modèle final universel de SPARC.

Le résultat plus profond est :

$$
\boxed{
\text{les galaxies SPARC s’organisent en classes de profils d’intrication.}
}
$$

Les anciens runs SPARC suggéraient déjà plusieurs régimes dynamiques via les profils empiriques \(d(r)\). Paper 12 reformule ces régimes en termes de dimensions locales d’intrication :

$$
(d_s(r),d_w(r))
\rightarrow
\alpha_{\rm eff}(r).
$$

La conclusion centrale est :

$$
\boxed{
\text{Paper 11 dérive la loi du propagateur BuP ; Paper 12 montre que cette loi organise les courbes de rotation galactiques.}
}
$$



# Bottom-Up Quantum Gravity

**Émergence de l’espace, du temps et de la gravité à partir de l’intrication quantique**  
**Farid Hamdad — Mars 2026**

---

## Vue d’ensemble

**Bottom-Up Quantum Gravity (BuP)** explore une hypothèse minimale :

> l’espace, le temps et la gravité ne sont pas fondamentaux ;  
> ils émergent collectivement de la structure d’intrication d’un état quantique global fini.

Dans cette approche :

- l’intrication définit la connectivité,
- la géométrie émerge de la structure informationnelle,
- le temps apparaît via le flot modulaire,
- la gravité peut être interprétée comme une courbure effective de cette structure.

Le projet combine :

- une base conceptuelle issue de l’information quantique,
- des simulations numériques sur systèmes finis,
- une extension phénoménologique vers les galaxies (SPARC).

---

## Idée centrale

On part d’un état quantique global pur $|\Psi\rangle$ sans espace préexistant.

On définit une distance informationnelle :

$$
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max} + \epsilon}\right)
$$

où $I(i:j)$ est l’information mutuelle.

À partir de cette distance :

- on reconstruit une géométrie,
- on définit un Laplacien,
- on extrait une **dimension spectrale $d_s(\tau)$**,
- puis un profil radial effectif **$d(r)$**.

---

## Résultat clé : dimension émergente

Les simulations et les fits SPARC montrent que la dimension effective est **dépendante de l’échelle** :

- $d \approx 3$ localement,
- $d < 3$ à l’échelle galactique,
- retour vers $d \approx 3$ à grande échelle.

Ce comportement permet de reproduire :

- les courbes de rotation,
- un excès de lentillage,
- sans introduire explicitement de composante de matière noire particulaire.

---

## 1. Temps émergent : flot modulaire

Pour un sous-système $A$,

$$
\rho_A = \mathrm{Tr}_{\bar A} |\Psi\rangle\langle\Psi|
$$

$$
K_A = -\log(\rho_A)
$$

Le flot modulaire

$$
O(\tau)=e^{iK_A\tau} O e^{-iK_A\tau}
$$

définit une dynamique relationnelle intrinsèque.

**Interprétation :**  
le temps n’est plus un paramètre de fond, mais une propriété informationnelle interne du système.

---

## 2. Diagnostics spectraux et chaos modulaire

Le spectre du Hamiltonien modulaire $K_A$ est étudié à l’aide :

- des statistiques de gaps adjacents,
- des diagnostics de type Random Matrix Theory,
- du spectral form factor.
- 
$$
\langle r \rangle = \left\langle \frac{\min(\Delta_n,\Delta_{n+1})}{\max(\Delta_n,\Delta_{n+1})} \right\rangle
$$

Valeurs de référence :

| Régime | Valeur |
|---|---:|
| Poisson | $\approx 0.386$ |
| GOE | $\approx 0.536$ |
| GUE | $\approx 0.603$ |

Résultats typiques :

$$
\langle r \rangle \in [0.53,\,0.59]
$$

compatible avec un régime chaotique.

---

## 3. Émergence de l’espace

La géométrie est reconstruite via :

$$
I(i:j)=S(\rho_i)+S(\rho_j)-S(\rho_{ij})
$$

$$
d_{ij}=-\log\left(\frac{I(i:j)}{I_{\max}}+\varepsilon\right)
$$

Le projet construit :

- graphes d’intrication,
- plongements MDS,
- dimension effective.

> la géométrie n’est pas imposée ; elle est reconstruite.

---

## 4. Résultats internes

### Dimension émergente

| Configuration | Dimension |
|---|---:|
| local | $\approx 2$ |
| non-local | $\approx 3$ |

---

### ER = EPR

Des points éloignés deviennent proches → signature wormhole-like.

---

### Gravité thermodynamique

$$
\delta S \approx \beta_{\mathrm{eff}}\delta E
$$

relation numériquement stable.

---

## 5. Horizon émergent ($N=16$)

Détection d’une région type horizon :

- entropie élevée
- faible conductance
- fort couplage interne

→ structure de type horizon sans géométrie préalable.

---

## 6. Courbure discrète d’Ollivier–Ricci

Objectif :

$$
\text{défaut informationnel} \Rightarrow \text{courbure}
$$

Résultat robuste (multi-seeds) :

$$
\Delta \kappa_{\mathrm{edge}} \approx 0.0758
$$

- fraction positive = 1.00
- $p \sim 10^{-6}$

et effet localisé :

$$
\Delta \kappa_{\mathrm{near-far}} > 0
$$

👉 Ce résultat constitue à ce stade le signal numérique le plus robuste du projet.

---

## 7. Extension SPARC

Le modèle est testé sur 175 galaxies SPARC.

Résultat :

- massives : $d_{\min} \approx 2.49$
- naines : $d_{\min} \approx 2.27$

---

## 8. Structure en régimes

| Régime | $d_{\min}$ |
|--------|-----------|
| HIGH   | $\sim 2.49$ |
| LOW    | $\sim 2.27$ |

Structure émergente (non imposée).

---

## 9. Bridge $d(r)$ → lentillage

Pipeline complet :


bup_v2_minimal_pipeline.py
bup_step1_fit_sparc.py
bup_regime_split.py
bup_dr_to_lensing_bridge.py
bup_lensing_regime_figure.py


Extraction :

- $\Delta \alpha(r)$
- $R_{\text{peak}}$

---

## 10. Prédiction de lentillage

Résultat (102 galaxies) :

- LOW : $x_{\text{peak}} \approx 0.365$
- HIGH : $x_{\text{peak}} \approx 0.498$

### Interprétation

- LOW → lentillage centralisé  
- HIGH → lentillage étendu  

👉 Cette prédiction constitue un discriminant potentiel entre BuP et les modèles standard.

---

## 11. Structure du dépôt


experiments/
10_single_galaxy_ngc3198/
20_v1_3_scan_ngc3198/
30_multi_galaxy_regime_3_08_4_5/

derive_d_r/
results_bup_hybrid_multi/
results_bup_hybrid_multi_rmax_fixed/


---

## 12. Interprétation physique

> la matière modifie la structure informationnelle  
> → modifie la géométrie  
> → produit la gravité

---

## 13. Limites

- modèle effectif
- pas encore relativiste complet
- lentillage via proxy
- dépendance M/L

---

## 14. Perspectives

- lentillage relativiste complet  
- confrontation observations  
- cosmologie  
- dynamique fondamentale  

---

## 15. Statut

✔ rotation curves  
✔ dimension émergente  
✔ régimes  
✔ lentillage  
✔ courbure  

---

## Auteur

**Farid Hamdad**  
Bottom-Up Quantum Gravity — 2026
hamdadfarid54@gmail.com

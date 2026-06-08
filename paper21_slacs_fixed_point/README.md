\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{physics}
\usepackage{bm}
\usepackage{authblk}

\geometry{margin=2.5cm}

\title{
\textbf{Paper 21 --- Le Point Fixe SLACS}\\
\large Validation par Lentillage Fort du Potentiel Effectif BuP
}

\author{Farid Hamdad}
\affil{Programme Gravité Quantique Bottom-Up}
\date{2026}

\begin{document}

\maketitle

\begin{abstract}
Nous testons une prédiction du cadre Gravité Quantique Bottom-Up sur l'échantillon de lentillage fort SLACS. Le Paper 9 a introduit un potentiel gravitationnel effectif généré par une source de matière émergente via une équation de Poisson sur graphe,
\[
L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}.
\]
Nous testons ici si ce potentiel prédit les résidus de lentillage au-delà d'un modèle baryonique de référence. L'observable est
\[
C_{\rm obs}
=
\log\left(
\frac{\theta_E^{\rm obs}}
{\theta_E^{\rm baryon}}
\right).
\]
Nous trouvons une transition robuste autour de
\[
\log M_\star\simeq 11,58-11,60,
\]
où \(C_{\rm obs}\simeq0\), où le potentiel BuP surpasse localement le proxy global de masse stellaire \(\log M_\star\), et où le secteur dimensionnel satisfait \(\alpha_{\rm eff}\simeq1\) à une échelle de diffusion intermédiaire. L'utilisation des indices de Sérsic mesurés renforce le signal par rapport à un contrôle à \(n=4\) forcé. Ces résultats fournissent un test de cohérence observationnelle non trivial du potentiel effectif BuP prédit dans le Paper 9.
\end{abstract}

\section{Introduction}

Le programme Gravité Quantique Bottom-Up propose que la géométrie et la gravité émergent de la structure de l'intrication quantique. Les papiers précédents ont développé la chaîne théorique
\[
W_{ij}
\rightarrow
L_{\rm ent}
\rightarrow
d_s,d_w
\rightarrow
\alpha_{\rm eff}
\rightarrow
\text{gravité effective}.
\]

Le Paper 9 a introduit une source de matière émergente \(S_{\rm flux}\) et l'a couplée à une équation de Poisson effective sur graphe,
\[
L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}.
\]
Le présent papier teste cette prédiction sur l'échantillon de lentillage fort SLACS.

La question centrale est de savoir si \(\Phi_{\rm BuP}\) laisse une empreinte observable dans les résidus de lentillage.

\section{Prédiction du Paper 9}

Le Paper 9 prédit la chaîne
\[
S_{\rm flux}
\rightarrow
\Phi_{\rm BuP}
\rightarrow
\text{réponse gravitationnelle}.
\]

La source est
\[
S_{\rm flux}
=
T_{00}
-\frac{1}{2}T_{aa}
+\frac{1}{2}T_{\rm grad}
+
|T_{0a}|.
\]

Le potentiel effectif est obtenu par
\[
L_{\rm ent}\Phi_{\rm BuP}=S_{\rm flux}.
\]

Dans l'application au lentillage, le potentiel BuP est testé à travers la correction résiduelle
\[
C_{\rm obs}
=
\log\left(
\frac{\theta_E^{\rm obs}}
{\theta_E^{\rm baryon}}
\right).
\]

\section{Données et observables}

Nous utilisons un catalogue de travail SLACS de lentillage fort contenant les rayons de lentillage, les masses stellaires, les rayons effectifs et les informations de profil photométrique.

La prédiction baryonique de référence est notée
\[
\theta_E^{\rm baryon}.
\]

Le rayon d'Einstein observé est
\[
\theta_E^{\rm obs}.
\]

Le résidu est
\[
C_{\rm obs}
=
\log\left(
\frac{\theta_E^{\rm obs}}
{\theta_E^{\rm baryon}}
\right).
\]

Une valeur positive indique un excès par rapport à la référence baryonique ; une valeur négative indique un déficit ou une sur-normalisation selon la baseline adoptée.

\section{Construction du graphe BuP}

Pour chaque galaxie, nous construisons un graphe à partir d'un profil photométrique projeté. La densité surfacique paramétrique de base est représentée par un profil de Sérsic
\[
\Sigma(R)
\propto
\exp\left[
-b_n
\left(
\frac{R}{R_e}
\right)^{1/n}
\right].
\]

Le couplage du graphe est pris comme
\[
W_{ij}
\propto
\sqrt{\Sigma_i\Sigma_j}
\exp\left(
-\frac{d_{ij}^2}{2\xi^2}
\right).
\]

Le Laplacien d'intrication est
\[
L_{\rm ent}=D-W.
\]

Le potentiel est obtenu à partir de la structure pseudo-inverse du Laplacien du graphe et de la construction de la source du Paper 9.

\section{Secteur dimensionnel}

L'exposant gravitationnel effectif est
\[
\alpha_{\rm eff}
=
\frac{2d_s}{d_w}+d_w-4.
\]

La condition de point fixe est
\[
\alpha_{\rm eff}=1.
\]

De manière équivalente,
\[
2d_s=d_w(5-d_w),
\]
ou
\[
d_s=\frac{d_w(5-d_w)}{2}.
\]

Pour la diffusion brownienne standard,
\[
d_w=2,
\]
et le point fixe donne
\[
d_s=3.
\]

Ainsi \(\alpha_{\rm eff}=1\) correspond au point fixe newtonien ou baryonique du secteur gravitationnel effectif.

\section{Test global du résidu}

Nous comparons d'abord le potentiel BuP au résidu baryonique
\[
C_{\rm obs}.
\]

La meilleure caractéristique BuP dynamique atteint une amélioration leave-one-out d'environ \(10,49\%\) sur l'échantillon de travail complet.

Cela établit que \(\Phi_{\rm BuP}\) porte une information prédictive sur le résidu de lentillage.

\begin{table}[h]
\centering
\begin{tabular}{lcc}
\toprule
Modèle & Amélioration LOO & Interprétation \\
\midrule
Baryonique baseline & 0\% & Référence \\
Régime LOW/HIGH & \(\sim 6-8\%\) & Signal BuP discret \\
Force BuP (forme seule) & \(\sim 9\%\) & Signal du graphe morphologique \\
\(\Phi_{\rm BuP}^{\sqrt{M_\star}}\) & \(10,49\%\) & Potentiel BuP dynamique \\
\(\sqrt{M_\star}\) & \(12,25\%\) & Proxy de masse stellaire \\
\(\log M_\star\) & \(12,92\%\) & Meilleur proxy global à un paramètre \\
\bottomrule
\end{tabular}
\caption{Résumé de la hiérarchie de prédiction des résidus.}
\end{table}

\section{Contrôle par mélange dynamique}

Pour tester si le signal est causé par l'association correcte entre amplitude stellaire et structure du graphe, nous effectuons un contrôle par mélange dynamique.

L'amplitude stellaire est permutée entre les galaxies avant la construction du graphe.

Cela détruit la majeure partie du signal BuP, réduisant le meilleur potentiel mélangé à une petite contribution résiduelle.

Cela démontre que le résultat dépend de l'association amplitude-structure correcte, pas seulement de la distribution marginale des masses stellaires.

\section{Indices de Sérsic mesurés}

Une comparaison contrôlée est effectuée sur les mêmes 61 galaxies avec des indices de Sérsic mesurés.

En utilisant les \(n_i\) mesurés, la meilleure caractéristique BuP dynamique atteint une amélioration LOO de
\[
10,40\%.
\]

Pour les mêmes galaxies forcées à \(n=4\), la meilleure caractéristique BuP dynamique tombe à
\[
9,12\%.
\]

Ainsi, la morphologie photométrique mesurée renforce le signal dynamique BuP.

\begin{table}[h]
\centering
\begin{tabular}{lcc}
\toprule
Catalogue & \(N\) & Meilleure amélioration LOO BuP dynamique \\
\midrule
\(n_i\) mesurés & 61 & \(10,40\%\) \\
\(n=4\) forcé & 61 & \(9,12\%\) \\
Mesurés + secours & 70 & \(10,49\%\) \\
\bottomrule
\end{tabular}
\caption{Effet de l'utilisation des indices de Sérsic mesurés.}
\end{table}

\section{Fenêtre de transition}

Une analyse par fenêtre glissante en masse stellaire montre que le potentiel BuP surpasse localement \(\log M_\star\) autour de
\[
11,545<\log M_\star<11,645.
\]

Avec les indices de Sérsic mesurés, dans cette fenêtre,
\[
\Phi_{\rm BuP}
\]
surpasse \(\log M_\star\) pour environ
\[
81,8\%
\]
des galaxies.

Cela suggère que le potentiel BuP capture une correction dynamique localisée dans le régime de transition.

\section{Point fixe observationnel}

Le point fixe observationnel est défini par
\[
C_{\rm obs}=0.
\]

L'ajustement de \(C_{\rm obs}\) en fonction de \(\log M_\star\) donne un passage par zéro autour de
\[
\log M_\star\simeq 11,58-11,60.
\]

Cela coïncide avec la fenêtre de transition dans laquelle \(\Phi_{\rm BuP}\) surpasse le proxy global de masse stellaire.

\section{Point fixe dimensionnel}

Un balayage des fenêtres de diffusion montre que le point fixe dimensionnel
\[
\alpha_{\rm eff}=1
\]
est retrouvé à une échelle de diffusion intermédiaire.

Pour l'échantillon à Sérsic mesuré, la fenêtre optimale est
\[
t_{\min}=1,\qquad t_{\max}=21.
\]

Cela donne
\[
\langle \alpha_{\rm eff}\rangle=1,014,
\]
avec
\[
d_w\simeq4,213,
\qquad
d_s\simeq1,688.
\]

À l'intérieur de la fenêtre de masse de transition,
\[
11,545<\log M_\star<11,645,
\]
la même fenêtre de diffusion donne
\[
\langle \alpha_{\rm eff}\rangle=1,014.
\]

Ainsi le point fixe dimensionnel et le point fixe observationnel coïncident.

\section{Convergence triple}

Le résultat principal du Paper 21 est la convergence triple
\[
C_{\rm obs}=0,
\qquad
\Phi_{\rm BuP}>\log M_\star,
\qquad
\alpha_{\rm eff}\simeq1,
\]
autour de
\[
\log M_\star\simeq11,6.
\]

C'est le point fixe SLACS.

Il n'est pas introduit comme un seuil libre. Il émerge indépendamment des données de résidu de lentillage, de la fenêtre de performance du potentiel BuP et du secteur dimensionnel du graphe.

\section{Discussion}

Le résultat ne doit pas être interprété comme une validation complète de l'ensemble du programme BuP. C'est une validation d'une prédiction spécifique du Paper 9 : le potentiel effectif BuP porte une information gravitationnelle observable.

Le graphe actuel reste paramétrique. Un test plus solide utiliserait des cartes de lumière HST non-paramétriques pour construire un graphe véritablement bidimensionnel.

Néanmoins, le fait que les indices de Sérsic mesurés renforcent le résultat indique que le potentiel BuP est sensible à la morphologie photométrique réelle, pas seulement à la masse stellaire.

\section{Conclusion}

Nous avons testé le potentiel effectif BuP du Paper 9 sur l'échantillon de lentillage fort SLACS.

Les principaux résultats sont :

\begin{enumerate}
\item Le potentiel BuP prédit les résidus de lentillage avec une amélioration LOO non nulle.
\item Les contrôles par mélange dynamique suppriment le signal.
\item Les indices de Sérsic mesurés renforcent la prédiction BuP.
\item Une fenêtre de transition apparaît autour de \(\log M_\star\simeq11,6\).
\item Le résidu observé satisfait \(C_{\rm obs}=0\) près de la même masse.
\item Le secteur dimensionnel atteint \(\alpha_{\rm eff}\simeq1\) à une échelle de diffusion intermédiaire dans la même fenêtre de masse.
\end{enumerate}

Ainsi, le Paper 21 fournit une validation observationnelle de la prédiction du Paper 9 selon laquelle le potentiel effectif BuP peut encoder les corrections de lentillage gravitationnel.

\bibliographystyle{unsrt}
\bibliography{references}

\end{document}

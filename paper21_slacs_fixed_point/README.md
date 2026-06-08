\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{mathtools}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{physics}
\usepackage{bm}
\usepackage{authblk}
\usepackage{microtype}

\geometry{margin=2.5cm}

\hypersetup{
    colorlinks=true,
    linkcolor=blue!60!black,
    citecolor=blue!60!black,
    urlcolor=blue!60!black
}

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
Nous testons une prédiction du cadre Gravité Quantique Bottom-Up (BuP) sur l'échantillon de lentillage fort SLACS.

Le Paper 9 a introduit un potentiel gravitationnel effectif généré par une source de matière émergente via une équation de Poisson sur graphe :

\[
L_{\rm ent}\,\Phi_{\rm BuP} = S_{\rm flux}.
\]

Nous testons ici si ce potentiel prédit les résidus de lentillage au-delà d'un modèle baryonique de référence. L'observable est :

\[
C_{\rm obs}
= \log\left(
\frac{\theta_E^{\rm obs}}{\theta_E^{\rm baryon}}
\right).
\]

Nous trouvons une transition robuste autour de :

\[
\log M_\star \simeq 11{,}58 - 11{,}60,
\]

où :

\begin{itemize}
    \item \(C_{\rm obs} \simeq 0\),
    \item le potentiel BuP surpasse localement le proxy global de masse stellaire \(\log M_\star\),
    \item le secteur dimensionnel satisfait \(\alpha_{\rm eff} \simeq 1\) à une échelle de diffusion intermédiaire.
\end{itemize}

L'utilisation des indices de Sérsic \textbf{mesurés} renforce le signal par rapport à un contrôle à \(n = 4\) forcé.

Ces résultats fournissent un test de cohérence observationnelle non trivial du potentiel effectif BuP prédit dans le Paper 9.
\end{abstract}

\tableofcontents

\section{Introduction}

Le programme Gravité Quantique Bottom-Up (BuP) propose que la géométrie et la gravité émergent de la structure de l'intrication quantique.

Les papiers précédents ont développé la chaîne théorique :

\[
W_{ij}
\;\rightarrow\;
L_{\rm ent}
\;\rightarrow\;
d_s,\; d_w
\;\rightarrow\;
\alpha_{\rm eff}
\;\rightarrow\;
\text{gravité effective}.
\]

Le Paper 9 a introduit une source de matière émergente \(S_{\rm flux}\) et l'a couplée à une équation de Poisson effective sur graphe :

\[
L_{\rm ent}\,\Phi_{\rm BuP} = S_{\rm flux}.
\]

Le présent papier teste cette prédiction sur l'échantillon de lentillage fort SLACS.

La question centrale est de savoir si \(\Phi_{\rm BuP}\) laisse une empreinte observable dans les résidus de lentillage.

\section{Prédiction du Paper 9}

Le Paper 9 prédit la chaîne causale suivante :

\[
S_{\rm flux}
\;\rightarrow\;
\Phi_{\rm BuP}
\;\rightarrow\;
\text{réponse gravitationnelle}.
\]

\subsection{Source de matière émergente}

La source est construite à partir de quatre contributions :

\[
S_{\rm flux}
= T_{00}
- \frac{1}{2}\,T_{aa}
+ \frac{1}{2}\,T_{\rm grad}
+ |T_{0a}|.
\]

\subsection{Potentiel effectif}

Le potentiel effectif est obtenu par résolution de l'équation de Poisson discrète :

\[
\boxed{
L_{\rm ent}\,\Phi_{\rm BuP} = S_{\rm flux}
}
\]

\subsection{Résidu de lentillage}

Dans l'application au lentillage, le potentiel BuP est testé à travers la correction résiduelle :

\[
\boxed{
C_{\rm obs}
= \log\left(
\frac{\theta_E^{\rm obs}}{\theta_E^{\rm baryon}}
\right)
}
\]

\section{Données et observables}

Nous utilisons un catalogue de travail SLACS de lentillage fort contenant :

\begin{itemize}
    \item les rayons d'Einstein observés \(\theta_E^{\rm obs}\),
    \item les masses stellaires \(\log M_\star\),
    \item les rayons effectifs \(R_e\),
    \item les indices de Sérsic \(n\).
\end{itemize}

\subsection{Prédiction baryonique de référence}

La prédiction baryonique de référence est notée \(\theta_E^{\rm baryon}\).

\subsection{Résidu observationnel}

Le résidu est défini par :

\[
C_{\rm obs}
= \log\left(
\frac{\theta_E^{\rm obs}}{\theta_E^{\rm baryon}}
\right).
\]

\begin{itemize}
    \item \(C_{\rm obs} > 0\) : excès par rapport à la référence baryonique,
    \item \(C_{\rm obs} < 0\) : déficit par rapport à la référence baryonique,
    \item \(C_{\rm obs} = 0\) : la prédiction baryonique coïncide avec l'observation.
\end{itemize}

\section{Construction du graphe BuP}

Pour chaque galaxie, nous construisons un graphe à partir d'un profil photométrique projeté.

\subsection{Profil de Sérsic}

La densité surfacique paramétrique de base est représentée par un profil de Sérsic :

\[
\Sigma(R)
\propto
\exp\left[
-b_n
\left(
\frac{R}{R_e}
\right)^{1/n}
\right],
\qquad
b_n = 2n - \frac{1}{3}.
\]

\subsection{Couplage du graphe}

Le couplage entre nœuds est pris sous forme exponentielle :

\[
\boxed{
W_{ij}
\propto
\sqrt{\Sigma_i \,\Sigma_j}
\;\exp\left(
-\frac{d_{ij}^2}{2\xi^2}
\right)
}
\]

\subsection{Laplacien d'intrication}

Le Laplacien normalisé du graphe est :

\[
L_{\rm ent} = D - W,
\]

où \(D\) est la matrice diagonale des degrés :

\[
D_{ii} = \sum_j W_{ij}.
\]

\subsection{Potentiel effectif}

Le potentiel est obtenu par pseudo-inversion du Laplacien :

\[
\Phi_{\rm BuP} = L_{\rm ent}^{+} \, S_{\rm flux}.
\]

\section{Secteur dimensionnel}

L'exposant gravitationnel effectif est donné par :

\[
\boxed{
\alpha_{\rm eff}
= \frac{2d_s}{d_w} + d_w - 4
}
\]

où :
\begin{itemize}
    \item \(d_s\) est la dimension spectrale du graphe,
    \item \(d_w\) est la dimension de marche (exposant de diffusion).
\end{itemize}

\subsection{Condition de point fixe}

La condition de point fixe gravitationnel est :

\[
\alpha_{\rm eff} = 1.
\]

De manière équivalente :

\[
\frac{2d_s}{d_w} + d_w = 5,
\]

ou encore :

\[
2d_s = d_w\,(5 - d_w),
\]

\[
\boxed{
d_s = \frac{d_w\,(5 - d_w)}{2}
}
\]

\subsection{Cas particulier : diffusion brownienne}

Pour la diffusion brownienne standard, \(d_w = 2\). Le point fixe donne alors :

\[
d_s = 3.
\]

Ainsi, \(\alpha_{\rm eff} = 1\) correspond au point fixe newtonien (ou baryonique) du secteur gravitationnel effectif.

\section{Test global du résidu}

Nous comparons d'abord le potentiel BuP au résidu baryonique \(C_{\rm obs}\).

\subsection{Résultat principal}

La meilleure caractéristique BuP dynamique atteint une amélioration \textbf{leave-one-out} (LOO) d'environ :

\[
\boxed{10{,}49\%}
\]

sur l'échantillon de travail complet.

Ce résultat établit que \(\Phi_{\rm BuP}\) porte une information prédictive sur le résidu de lentillage.

\subsection{Hiérarchie des modèles}

\begin{table}[h]
\centering
\caption{Résumé de la hiérarchie de prédiction des résidus.}
\begin{tabular}{lcc}
\toprule
\textbf{Modèle} & \textbf{Amélioration LOO} & \textbf{Interprétation} \\
\midrule
Baryonique baseline & \(0\%\) & Référence \\
Régime LOW/HIGH & \(\sim 6-8\%\) & Signal BuP discret \\
Force BuP (forme seule) & \(\sim 9\%\) & Signal du graphe morphologique \\
\(\Phi_{\rm BuP}^{\sqrt{M_\star}}\) & \(10{,}49\%\) & Potentiel BuP dynamique \\
\(\sqrt{M_\star}\) & \(12{,}25\%\) & Proxy de masse stellaire \\
\(\log M_\star\) & \(12{,}92\%\) & Meilleur proxy global à un paramètre \\
\bottomrule
\end{tabular}
\end{table}

\section{Contrôle par mélange dynamique}

Pour tester si le signal est \textbf{causal}, nous effectuons un contrôle par mélange dynamique.

\subsection{Protocole}

L'amplitude stellaire \(M_\star\) est permutée \textbf{aléatoirement} entre les galaxies \emph{avant} la construction du graphe.

\subsection{Résultat}

Le mélange détruit la majeure partie du signal BuP :

\[
\text{Meilleur potentiel mélangé} \ll 10{,}49\%.
\]

\subsection{Conclusion du contrôle}

Ce contrôle démontre que :

\begin{quote}
Le résultat dépend de l'\textbf{association correcte} entre amplitude stellaire et structure du graphe, pas seulement de la distribution marginale des masses stellaires.
\end{quote}

\section{Indices de Sérsic mesurés}

Une comparaison contrôlée est effectuée sur les mêmes 61 galaxies avec des indices de Sérsic \textbf{mesurés}.

\subsection{Résultat}

\begin{table}[h]
\centering
\caption{Effet de l'utilisation des indices de Sérsic mesurés.}
\begin{tabular}{lcc}
\toprule
\textbf{Catalogue} & \(N\) & \textbf{Meilleure amélioration LOO BuP} \\
\midrule
\(n_i\) mesurés & 61 & \(10{,}40\%\) \\
\(n = 4\) forcé & 61 & \(9{,}12\%\) \\
Mesurés + secours & 70 & \(10{,}49\%\) \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Interprétation}

L'utilisation des indices de Sérsic \textbf{mesurés} renforce le signal dynamique BuP d'environ \(1{,}3\) point de pourcentage.

Ainsi, le potentiel BuP est sensible à la \textbf{morphologie photométrique réelle}, pas seulement à la masse stellaire.

\section{Fenêtre de transition}

Une analyse par \textbf{fenêtre glissante} en masse stellaire révèle un régime de transition localisé.

\subsection{Localisation de la fenêtre}

La fenêtre de transition est :

\[
\boxed{
11{,}545 < \log M_\star < 11{,}645
}
\]

\subsection{Performance dans la fenêtre}

Avec les indices de Sérsic mesurés, dans cette fenêtre :

\[
\Phi_{\rm BuP} \quad \text{surpasse} \quad \log M_\star
\]

pour environ :

\[
\boxed{81{,}8\%}
\]

des galaxies.

\subsection{Interprétation}

Cela suggère que le potentiel BuP capture une \textbf{correction dynamique localisée} dans le régime de transition.

\section{Point fixe observationnel}

Le point fixe observationnel est défini par l'annulation du résidu :

\[
C_{\rm obs} = 0.
\]

\subsection{Détermination par régression}

L'ajustement linéaire de \(C_{\rm obs}\) en fonction de \(\log M_\star\) donne un passage par zéro à :

\[
\boxed{
\log M_\star \simeq 11{,}58 - 11{,}60
}
\]

\subsection{Cohérence}

Cette valeur coïncide avec la fenêtre de transition dans laquelle \(\Phi_{\rm BuP}\) surpasse le proxy global de masse stellaire.

\section{Point fixe dimensionnel}

Un balayage des fenêtres de diffusion (\(t_{\min}, t_{\max}\)) montre que le point fixe dimensionnel est accessible.

\subsection{Fenêtre optimale}

Pour l'échantillon à Sérsic mesuré, la fenêtre de diffusion optimale est :

\[
t_{\min} = 1, \qquad t_{\max} = 21.
\]

\subsection{Valeurs obtenues}

Cette fenêtre donne :

\[
\langle \alpha_{\rm eff} \rangle = 1{,}014,
\]
\[
d_w \simeq 4{,}213,
\qquad
d_s \simeq 1{,}688.
\]

\subsection{Valeur dans la fenêtre de masse}

À l'intérieur de la fenêtre de transition \(11{,}545 < \log M_\star < 11{,}645\), la même fenêtre de diffusion donne :

\[
\langle \alpha_{\rm eff} \rangle = 1{,}014.
\]

\subsection{Conclusion}

Ainsi, le \textbf{point fixe dimensionnel} (\(\alpha_{\rm eff} = 1\)) et le \textbf{point fixe observationnel} (\(C_{\rm obs} = 0\)) \textbf{coïncident} à la même échelle de masse.

\section{Convergence triple}

Le résultat principal du Paper 21 est la \textbf{convergence triple} :

\[
\boxed{
C_{\rm obs} = 0
\qquad\text{et}\qquad
\Phi_{\rm BuP} > \log M_\star
\qquad\text{et}\qquad
\alpha_{\rm eff} \simeq 1
}
\]

au voisinage de :

\[
\boxed{
\log M_\star \simeq 11{,}6
}
\]

\subsection{Statut du point fixe}

Ce point fixe SLACS :

\begin{itemize}
    \item n'est \textbf{pas introduit} comme un seuil libre,
    \item \textbf{émerge indépendamment} :
    \begin{enumerate}
        \item des données de résidu de lentillage,
        \item de la fenêtre de performance du potentiel BuP,
        \item du secteur dimensionnel du graphe.
    \end{enumerate}
\end{itemize}

\section{Discussion}

\subsection{Portée du résultat}

Ce résultat ne doit pas être interprété comme une validation complète de l'ensemble du programme BuP.

C'est une validation d'une \textbf{prédiction spécifique du Paper 9} :

\begin{quote}
Le potentiel effectif BuP porte une information gravitationnelle observable.
\end{quote}

\subsection{Limitations actuelles}

Le graphe utilisé reste \textbf{paramétrique} (profil radial de Sérsic).

Un test plus solide utiliserait :
\begin{itemize}
    \item des cartes de lumière HST \textbf{non-paramétriques},
    \item un graphe \textbf{véritablement bidimensionnel}.
\end{itemize}

\subsection{Force du résultat}

Néanmoins, le fait que les indices de Sérsic \textbf{mesurés} renforcent le résultat indique que :

\begin{quote}
Le potentiel BuP est sensible à la morphologie photométrique \textbf{réelle}, pas seulement à la masse stellaire.
\end{quote}

\section{Conclusion}

Nous avons testé le potentiel effectif BuP du Paper 9 sur l'échantillon de lentillage fort SLACS.

\subsection{Résultats principaux}

\begin{enumerate}
    \item Le potentiel BuP prédit les résidus de lentillage avec une amélioration LOO non nulle (\(\sim 10{,}5\%\)).
    \item Les contrôles par mélange dynamique suppriment le signal (preuve de causalité).
    \item Les indices de Sérsic \textbf{mesurés} renforcent la prédiction BuP.
    \item Une fenêtre de transition apparaît autour de \(\log M_\star \simeq 11{,}6\).
    \item Le résidu observé satisfait \(C_{\rm obs} = 0\) près de la même masse.
    \item Le secteur dimensionnel atteint \(\alpha_{\rm eff} \simeq 1\) à une échelle de diffusion intermédiaire dans la même fenêtre de masse.
\end{enumerate}

\subsection{Conclusion finale}

Le Paper 21 fournit une \textbf{validation observationnelle} de la prédiction du Paper 9 :

\begin{quote}
Le potentiel effectif BuP \(\Phi_{\rm BuP}\) peut encoder les corrections de lentillage gravitationnel.
\end{quote}

La convergence triple autour de \(\log M_\star \simeq 11{,}6\) constitue un test de cohérence non trivial du cadre BuP.

\bibliographystyle{unsrt}
\bibliography{references}

\end{document}

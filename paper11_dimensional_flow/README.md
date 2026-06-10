\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{microtype}

\geometry{margin=2.5cm}

\hypersetup{
    colorlinks=true,
    linkcolor=blue!60!black,
    citecolor=blue!60!black,
    urlcolor=blue!60!black
}

\title{
\textbf{Paper 11 --- Le Flot Dimensionnel de la Gravité BuP}\\
\large Dimension Spectrale, Dimension de Marche et Exposant Gravitationnel Effectif
}

\author{Farid Hamdad}
\date{2026}

\begin{document}

\maketitle

\begin{abstract}
Ce dossier contient le matériel numérique et théorique associé au \textbf{Paper 11} du programme Gravité Quantique Bottom-Up (BuP).

Le Paper 11 introduit la loi dimensionnelle qui relie les propriétés spectrales d'un graphe d'intrication à un exposant gravitationnel effectif :

\[
\boxed{
\alpha_{\rm eff}
= \frac{2d_s}{d_w} + d_w - 4
}
\]

où :
\begin{itemize}
    \item \(d_s\) est la \textbf{dimension spectrale} du graphe d'intrication,
    \item \(d_w\) est la \textbf{dimension de marche},
    \item \(\alpha_{\rm eff}\) est l'\textbf{exposant gravitationnel effectif} qui contrôle la réponse à grande échelle émergente.
\end{itemize}

Ce papier constitue le pont entre le programme microscopique du graphe et les tests phénoménologiques réalisés plus tard dans la séquence BuP.
\end{abstract}

\tableofcontents

\section{Rôle dans le programme BuP}

Le Paper 11 fournit la colonne vertébrale dimensionnelle pour les papiers ultérieurs :

\begin{verbatim}
Paper 11 :
    W_ij → L_ent → d_s, d_w → alpha_eff

Paper 12 :
    Sigma(R) → W_ij → L_ent → alpha_eff → V(r)

Paper 14 :
    Courbes de rotation SPARC et régimes LOW/HIGH

Paper 21 :
    Point fixe de lentillage fort SLACS et alpha_eff ≈ 1
\end{verbatim}

Ainsi, le Paper 11 n'ajuste pas un catalogue de galaxies. Il établit la \textbf{loi théorique et numérique} que les papiers ultérieurs testent observationnellement.

\section{Équation centrale}

Le résultat central est :

\[
\boxed{
\alpha_{\rm eff}
= \frac{2d_s}{d_w} + d_w - 4
}
\]

Cette formule combine deux quantités de diffusion sur graphe :

\[
P(t) \sim t^{-d_s/2}
\]
et
\[
\langle r^2(t) \rangle \sim t^{2/d_w}.
\]

L'exposant gravitationnel BuP n'est donc pas introduit arbitrairement. Il est \textbf{déduit de la diffusion sur le graphe d'intrication}.

\section{Point fixe newtonien}

Le point fixe newtonien (ou baryonique) correspond à :

\[
\alpha_{\rm eff} = 1.
\]

Par conséquent :

\[
\frac{2d_s}{d_w} + d_w - 4 = 1.
\]

De manière équivalente :

\[
2d_s = d_w\,(5 - d_w),
\]

ou encore :

\[
\boxed{
d_s = \frac{d_w\,(5 - d_w)}{2}
}
\]

\subsection{Cas particulier : diffusion brownienne}

Pour la diffusion brownienne standard, \(d_w = 2\). La condition du point fixe donne alors :

\[
d_s = 3.
\]

Ainsi, le comportement newtonien tridimensionnel ordinaire apparaît comme un \textbf{point fixe particulier} du flot dimensionnel BuP.

\section{Principaux résultats numériques}

Les résultats attendus sont :

\begin{verbatim}
results/
  alpha_eff_table.csv
  finite_size_summary.csv
  alpha_predictions_vs_N.csv
  paper11_summary.json
\end{verbatim}

et les figures :

\begin{verbatim}
figures/
  fig1_pipeline_dimensional_flow.png
  fig2_ds_vs_N.png
  fig3_dw_vs_N.png
  fig4_alpha_predictions_vs_N.png
  fig5_alpha_fixed_point_curve.png
  fig6_interpretation_regimes.png
\end{verbatim}

\section{Interprétation}

Le Paper 11 montre que le comportement gravitationnel effectif est contrôlé par le couple :

\[
(d_s,\; d_w).
\]

Différents régimes correspondent à différentes réponses gravitationnelles :

\begin{table}[h]
\centering
\caption{Régimes gravitationnels du flot dimensionnel BuP}
\begin{tabular}{lcc}
\toprule
\textbf{Régime} & \textbf{Condition} & \textbf{Interprétation} \\
\midrule
Point fixe newtonien & \(\alpha_{\rm eff} = 1\) & Comportement baryonique/newtonien ordinaire \\
Sous-newtonien & \(\alpha_{\rm eff} < 1\) & Régime affaibli ou sous-couplé \\
Super-newtonien & \(\alpha_{\rm eff} > 1\) & Réponse effective renforcée \\
Régime de transition & \(\alpha_{\rm eff} \approx 1\) & Crossover entre phases du graphe \\
\bottomrule
\end{tabular}
\end{table}

\section{Reproductibilité}

Exécuter :

\begin{verbatim}
cd papers/paper11_dimensional_flow
bash scripts/run_all.sh
\end{verbatim}

Cela génère toutes les tables numériques et les figures.

\section{Statut}

Le Paper 11 doit être lu comme une \textbf{dérivation théorique et numérique} de la loi dimensionnelle BuP. Ses conséquences observationnelles sont testées plus tard dans les Papers 12, 14 et 21.

\section*{Résumé pour mémoire}

\begin{itemize}
    \item Paper 11 établit la loi \(\alpha_{\rm eff} = 2d_s/d_w + d_w - 4\)
    \item \(d_s\) et \(d_w\) sont extraits du laplacien du graphe d'intrication
    \item Le point fixe newtonien correspond à \(\alpha_{\rm eff} = 1\), soit \(d_s = 3\) lorsque \(d_w = 2\)
    \item Cette loi est testée observationnellement dans Papers 12, 14 et 21
\end{itemize}

\end{document}

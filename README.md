\documentclass[11pt,a4paper]{article}

% -------------------- Encodage / Langue --------------------
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[french]{babel}

% -------------------- Maths / Symboles --------------------
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{bm}
\usepackage{braket}

% -------------------- Mise en page / Figures --------------------
\usepackage{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{setspace}
\usepackage{booktabs}
\usepackage{authblk}
\usepackage{tcolorbox}
\usepackage{xcolor}
\usepackage{siunitx}
\usepackage{hyperref}
\usepackage{url}

\geometry{margin=2.5cm}
\setstretch{1.15}

\definecolor{myblue}{RGB}{0,82,147}
\definecolor{mygreen}{RGB}{0,128,0}
\definecolor{myred}{RGB}{180,0,0}

\hypersetup{
  colorlinks=true,
  linkcolor=myblue,
  citecolor=myblue,
  urlcolor=myblue
}

\title{\textbf{Bottom-Up Quantum Gravity :}\\
\vspace{0.3em}
\Large{Émergence de l’espace, du temps et de la gravité à partir de l’intrication quantique}}

\author{Farid Hamdad}
\date{Février 2026}

\begin{document}
\maketitle

\vspace{0.8em}
\begin{center}
\itshape
L’espace, le temps et la gravité ne sont pas fondamentaux.\\
Ils émergent collectivement de la structure d’intrication d’un état quantique global fini.
\end{center}
\vspace{1.2em}

\begin{abstract}
\begin{tcolorbox}[colback=gray!10,colframe=myblue,title=Pourquoi ce projet ?]
La physique moderne décrit avec précision la physique quantique et la gravitation classique,
mais laisse ouverte une question fondamentale : 
\textbf{pourquoi l’espace-temps existe-t-il, et pourquoi la gravité possède-t-elle une structure géométrique et thermodynamique ?}

Ce projet explore une hypothèse minimale : 
\textbf{l’espace-temps n’est pas le théâtre de la physique. Il est reconstruit à partir de l’intrication quantique.}
\end{tcolorbox}
\end{abstract}

%=============================================================================
\section{Fondement : Postulat minimal}
%=============================================================================

Il existe un état quantique global pur
\begin{equation}
\ket{\Psi} \in \mathcal{H} = \bigotimes_{i=1}^{N} \mathcal{H}_i
\end{equation}
défini sur $N$ degrés de liberté élémentaires (qubits),
\textbf{sans espace, sans temps, sans métrique préalable}.

Tout le reste — temps, espace, dimension, géométrie, gravité effective —
doit émerger exclusivement de la \textbf{structure interne de l’intrication}.

%=============================================================================
\section{Méthodologie d’émergence}
%=============================================================================

\subsection{Émergence du temps : flot modulaire}

Pour un sous-système $A$ :
\begin{equation}
\rho_A = \mathrm{Tr}_{\bar{A}} \ket{\Psi}\bra{\Psi}, \qquad
K_A = -\log \rho_A.
\end{equation}

Le \textbf{flot modulaire} :
\begin{equation}
\mathcal{O}(\tau) = e^{i K_A \tau} \mathcal{O} e^{-i K_A \tau}
\end{equation}
définit une dynamique intrinsèque relationnelle (Page–Wootters).
Le temps devient une propriété informationnelle interne.

\subsection{Chaos modulaire (nouveau résultat)}

Nous analysons le spectre de $K_A$.

\paragraph{Ratio des gaps}
\begin{equation}
\langle r \rangle = \left\langle \frac{\min(\Delta_n, \Delta_{n+1})}{\max(\Delta_n, \Delta_{n+1})} \right\rangle
\end{equation}

Références universelles :
\begin{itemize}
  \item Poisson (intégrable) $\approx 0.386$
  \item GOE $\approx 0.536$
  \item GUE $\approx 0.603$
\end{itemize}

\textbf{Résultat :} $\langle r \rangle \in [0.53, 0.59]$ \\
$\Rightarrow$ Le flot modulaire devient chaotique (régime Random Matrix Theory).

\subsection{Spectral Form Factor (SFF)}
\begin{equation}
g_2(t) = \frac{1}{d_A^2} \left| \sum_{n=1}^{d_A} e^{-i t \tilde{\kappa}_n} \right|^2
\end{equation}

Structure observée : \textbf{dip – ramp – plateau}.

\textbf{Résultat central :}
\begin{equation}
g_2^{\text{plateau}} \sim \frac{1}{d_A}
\end{equation}
Scaling universel.

\subsection{Constante modulaire topologique}

On définit :
\begin{equation}
C = d_A \times g_2^{\text{plateau}}
\end{equation}

Le scaling $1/d_A$ est universel. Mais le \textbf{préfacteur $C$ dépend de la topologie} :

\begin{table}[h]
\centering
\begin{tabular}{@{}lcc@{}}
\toprule
Topologie & $\langle r \rangle$ & $C$ \\
\midrule
Chaîne 1D   & 0.594 & 1.21 \\
Grille 3×6  & 0.575 & 1.42 \\
Graphe ER   & 0.528 & 1.63 \\
\bottomrule
\end{tabular}
\caption{Dépendance topologique de la constante modulaire ($d_A=256$).}
\end{table}

\textbf{Interprétation :} 
Le flot modulaire devient universellement chaotique.
Mais la \textbf{géométrie d’intrication} laisse une empreinte quantitative via $C$.
La structure spatiale émergente influence le générateur temporel modulaire.

%=============================================================================
\section{Émergence de l’espace}
%=============================================================================

\subsection{Reconstruction géométrique}

\textbf{Information mutuelle}
\begin{equation}
I(i:j) = S(\rho_i) + S(\rho_j) - S(\rho_{ij})
\end{equation}

\textbf{Distance informationnelle}
\begin{equation}
d_{ij} = -\log\left( \frac{I(i:j)}{I_{\max} + \epsilon} \right)
\end{equation}

\textbf{MDS} $\rightarrow$ points $x_i \in \mathbb{R}^d$

\textbf{Dimension effective} = dimension minimale stabilisant l’erreur.

%=============================================================================
\section{Résultats principaux}
%=============================================================================

\subsection{Dimension émergente}

\begin{table}[h]
\centering
\begin{tabular}{@{}lcc@{}}
\toprule
Configuration & Intrication & Dimension \\
\midrule
N=9, $\lambda\approx 0$ & Locale      & $d\approx 2$ \\
N=9, $\lambda\to 1$     & Non-locale  & $d\approx 3$ \\
N=16, $\lambda\approx 0$ & Locale      & $d\approx 2$ \\
N=16, $\lambda\to 1$     & Non-locale  & $d\approx 3$ \\
\bottomrule
\end{tabular}
\caption{Transition dimensionnelle contrôlée par l’intrication.}
\end{table}

\subsection{ER = EPR mesurable}

Des qubits \textbf{topologiquement distants} deviennent \textbf{géométriquement proches} lorsque l’intrication non-locale augmente.
Signature \textit{wormhole-like} discrète (statique).

\subsection{Gravité thermodynamique}

Test de type \textbf{Jacobson} :
\begin{equation}
\delta S \simeq \beta_{\mathrm{eff}} \, \delta E
\end{equation}
Relation stable pour $N=9$ et $N=16$.

%=============================================================================
\section{Détection d’horizon bottom-up (N = 16)}
%=============================================================================

Nous introduisons un \textbf{benchmark de détection d’horizon émergent} basé uniquement sur la structure du graphe d’intrication.

\subsection{Principe}

\begin{enumerate}
  \item Construction du graphe pondéré $W_{ij} = I(i:j)$
  \item Seuil par densité fixée ($\rho = 1/3$)
  \item Recherche d’une région de taille $k = N/2$ maximisant :
  \begin{equation}
  \text{Score}_{\mathrm{BH}} = w_S \, z(S) - w_\phi \, z(\phi) + w_{\mathrm{int}} \, z(\text{internal})
  \end{equation}
  où $S$ = entropie, $\phi$ = conductance, internal = intrication interne, $z(\cdot)$ = score normalisé.
\end{enumerate}

\subsection{Région BH-like détectée}

\begin{equation}
[0, 2, 3, 5, 6, 10, 11, 15]
\end{equation}

\begin{table}[h]
\centering
\begin{tabular}{@{}lc@{}}
\toprule
Quantité & Valeur \\
\midrule
Entropie $S$    & 7.1667 \\
Cut             & 3.2098 \\
Internal        & 5.0348 \\
Conductance $\phi$ & 0.3951 \\
Ratio $r$ (RMT) & 0.6041 \\
\bottomrule
\end{tabular}
\caption{Métriques de la région BH-like.}
\end{table}

\subsection{Comparaison via HIE fixed-size + Louvain}

\textbf{55 régions candidates} générées via :
\begin{itemize}
  \item Louvain (détection de communautés)
  \item Unions de communautés
  \item Recherche locale par swaps (Metropolis)
\end{itemize}

\subsubsection{Classement}

\begin{table}[h]
\centering
\begin{tabular}{@{}lc@{}}
\toprule
Mode & Rang BH \\
\midrule
CLASSIC        & 12 / 55 \\
HORIZON-AWARE  & 2 / 55 \\
\bottomrule
\end{tabular}
\caption{Classement de la région BH selon le mode de détection.}
\end{table}

Le mode \textbf{CLASSIC} détecte des communautés internes. \\
Le mode \textbf{HORIZON} détecte des bottlenecks entropiques.

\subsection{Percentiles vs régions aléatoires}

\begin{table}[h]
\centering
\begin{tabular}{@{}lc@{}}
\toprule
Quantité & Percentile BH \\
\midrule
$S$       & 23.7 \% \\
Cut       & 1.4 \% \\
Internal  & 98.8 \% \\
$\phi$    & 2.7 \% \\
$r$       & 62.0 \% \\
\bottomrule
\end{tabular}
\caption{Percentiles de la région BH par rapport à 3000 régions aléatoires de même taille.}
\end{table}

\textbf{Signature :}
\begin{itemize}
  \item Intrication interne \textbf{très élevée}
  \item Couplage extérieur \textbf{très faible}
  \item Bottleneck informationnel \textbf{net}
\end{itemize}

\begin{tcolorbox}[colback=mygreen!5,colframe=mygreen]
\textbf{Conclusion :} Un horizon-like peut émerger sans géométrie préalable.
\end{tcolorbox}

%=============================================================================
\section{Ce que le projet établit}
%=============================================================================

\begin{itemize}
  \item[\checkmark] Une \textbf{géométrie} peut émerger d’un état quantique fini
  \item[\checkmark] La \textbf{dimension} dépend de l’intrication
  \item[\checkmark] \textbf{ER=EPR} est mesurable
  \item[\checkmark] Une \textbf{thermodynamique d’intrication} apparaît
  \item[\checkmark] Le \textbf{flot modulaire} est chaotique
  \item[\checkmark] Le plateau \textbf{SFF} suit $1/d_A$
  \item[\checkmark] La \textbf{constante modulaire} dépend de la topologie
  \item[\checkmark] Un \textbf{horizon-like bottleneck} peut émerger bottom-up
\end{itemize}

%=============================================================================
\section{Limites}
%=============================================================================

\begin{itemize}
  \item[\texttimes] Limite continue $N \to \infty$ non démontrée
  \item[\texttimes] Pas de dynamique relativiste complète
  \item[\texttimes] Pas de dérivation des équations d’Einstein
  \item[\texttimes] Pas encore de prédiction expérimentale
\end{itemize}

%=============================================================================
\section{Organisation du dépôt}
%=============================================================================

\begin{verbatim}
bottom-up-quantum-gravity/
├── paper/
├── figures/
├── scripts/
│   └── bh_benchmark_louvain_N16.py
├── results/
│   └── examples/
│       └── N16_horizon_aware/
├── code/
├── data/
└── README.md
\end{verbatim}

%=============================================================================
\section{Perspectives}
%=============================================================================

\begin{itemize}
  \item \textbf{Finite-size scaling} $N \to 25+$
  \item États critiques / topologiques
  \item Connexion avec \textbf{SYK}
  \item Étude multi-échelle (dimension spectrale)
  \item Implémentation sur \textbf{simulateur quantique}
\end{itemize}

%=============================================================================
\section{Citation}
%=============================================================================

\begin{verbatim}
@misc{hamdad2026bottomup,
  author = {Hamdad, Farid},
  title = {Bottom-Up Quantum Gravity: Emergence of Space, 
           Time and Gravity from Quantum Entanglement},
  year = {2026},
  howpublished = {GitHub repository},
  url = {https://github.com/Farid-Hamdad/Bottom-Up-Quantum-Gravity}
}
\end{verbatim}

\begin{center}
\rule{0.5\textwidth}{0.5pt}
\textbf{contact} : \url{hamdadfarid54@gmail.com}
\end{center}

\end{document}

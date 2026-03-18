
\subsection{Prototype de lentillage optique BuP centralisé}

Afin d’aller au-delà des enveloppes phénoménologiques inverses et de rapprocher plus directement le programme Bottom-Up d’un observable de lentillage, nous avons introduit un prototype de réponse optique effective fondé sur le profil de dimension BuP calibré sur SPARC.

L’idée minimale est la suivante. Au lieu de modifier directement une densité projetée, on considère qu’un rayon lumineux se propage dans une réponse gravitationnelle effective dérivée du profil de dimension émergente. On définit d’abord un profil radial de dimension
\begin{equation}
d(r)=3-\delta_d \exp\!\left[-\frac{\ln^2(r/r_0)}{2\sigma^2}\right],
\end{equation}
puis un facteur géométrique associé
\begin{equation}
f_d(r)=\left(\frac{r}{r_{\mathrm{ref}}}\right)^{3-d(r)}.
\end{equation}
Dans une première version optique, l’accélération effective est prise sous la forme
\begin{equation}
g_{\mathrm{eff}}(r)=g_{\mathrm{bary}}(r)\,f_d(r),
\end{equation}
et la déflexion est estimée par intégration le long de la ligne de visée,
\begin{equation}
\alpha(R)\simeq \frac{2}{c^2}\int_{-\infty}^{+\infty} g_\perp(R,z)\,dz,
\qquad
g_\perp(R,z)=g_{\mathrm{eff}}(r)\frac{R}{r},
\qquad
r=\sqrt{R^2+z^2}.
\end{equation}

Les premières versions de ce schéma (v1.0--v1.1) produisent bien un excès de déflexion par rapport au cas baryonique pur, mais cet excès tend à être artificiellement porté par les grands rayons. Une étape décisive a donc consisté à introduire une \emph{localisation optique} de l’effet BuP. Dans la version centralisée retenue, l’accélération est remplacée par
\begin{equation}
g_{\mathrm{eff}}(r)=g_{\mathrm{bary}}(r)\Bigl[1+\lambda\,w_{\mathrm{tot}}(r)\,\bigl(f_d(r)-1\bigr)\Bigr],
\end{equation}
où \(w_{\mathrm{tot}}(r)\) est un poids radial combinant une localisation centrale et une coupure externe. Le centre effectif du poids est plafonné à une valeur
\begin{equation}
x_{\mathrm{center}}=\min(x_0,x_{\mathrm{cap}}),
\end{equation}
afin d’éviter que le maximum du bump SPARC, parfois situé trop loin dans le disque, ne déporte artificiellement la réponse optique vers l’extérieur.

\subsubsection{Scan de phase sur NGC3198}

Un premier scan a été effectué sur NGC3198 afin d’identifier un régime où l’effet BuP reste à la fois \emph{centralisé} et \emph{non négligeable}. Le scan portait sur les paramètres
\[
\lambda_{\mathrm{BuP}},\qquad \mathrm{width},\qquad x_{\mathrm{cap}},\qquad r_{\mathrm{cut}},
\]
avec 144 combinaisons testées.

Le meilleur compromis trouvé est
\begin{equation}
(\lambda_{\mathrm{BuP}},\mathrm{width},x_{\mathrm{cap}},r_{\mathrm{cut}})
=
(3.0,\,0.8,\,4.0,\,5.0).
\end{equation}
Dans ce régime, NGC3198 présente :
\begin{itemize}
    \item un gain robuste moyen de déflexion d’environ \(1.36\),
    \item une médiane robuste d’environ \(1.33\),
    \item un excès maximal
    \begin{equation}
    \Delta\alpha_{\max}\simeq 0.122\ \mathrm{arcsec},
    \end{equation}
    \item et un pic de l’excès situé à
    \begin{equation}
    R_{\Delta\alpha,\max}\simeq 11\ \mathrm{kpc}.
    \end{equation}
\end{itemize}

Ce résultat est important pour deux raisons. D’une part, il montre qu’un effet BuP optique \emph{centralisé} n’est pas nécessairement faible. D’autre part, il montre qu’un simple recentrage du profil sans scan (version centralisée par défaut) sous-estime fortement l’amplitude du signal.

\subsubsection{Test multi-galaxies sur spirales massives}

Le même régime
\[
(3.0,\,0.8,\,4.0,\,5.0)
\]
a ensuite été appliqué sans réajustement à plusieurs galaxies spirales massives de l’échantillon SPARC. Les cas les plus significatifs obtenus sont :

\begin{itemize}
    \item \textbf{NGC3198} :
    \(\Delta\alpha_{\max}\simeq 0.122\) arcsec,
    gain robuste moyen \(\simeq 1.36\).

    \item \textbf{NGC5055} :
    \(\Delta\alpha_{\max}\simeq 0.224\) arcsec,
    gain robuste moyen \(\simeq 1.33\).

    \item \textbf{NGC2841} :
    \(\Delta\alpha_{\max}\simeq 0.521\) arcsec,
    gain robuste moyen \(\simeq 1.39\).

    \item \textbf{NGC7331} :
    \(\Delta\alpha_{\max}\simeq 0.300\) arcsec,
    gain robuste moyen \(\simeq 1.51\).

    \item \textbf{UGC02885} :
    \(\Delta\alpha_{\max}\simeq 0.442\) arcsec,
    gain robuste moyen \(\simeq 1.49\).
\end{itemize}

Dans l’ensemble, ce régime produit donc un excès de déflexion robuste de l’ordre de \textbf{30 à 50\%} sur plusieurs spirales massives, avec une structure radiale non triviale mais encore globalement centrée sur des rayons galactiques internes à intermédiaires. Les pics de \(\Delta\alpha\) se situent typiquement entre \(\sim 8\) et \(\sim 14\) kpc pour NGC3198, NGC5055, NGC2841 et NGC7331, tandis que UGC02885 présente un maximum plus externe.

\subsubsection{Interprétation provisoire}

Ce résultat constitue, à ce stade, l’indication la plus nette en faveur d’une \emph{branche optique BuP} distincte de la simple branche de densité locale. Il suggère qu’un même profil de dimension émergente, lorsqu’il est converti en réponse optique effective avec une localisation centrale appropriée, peut induire un excès de lentillage systématique sur une sous-classe de galaxies spirales massives.

La lecture la plus prudente est la suivante :
\begin{enumerate}
    \item la branche purement locale, fondée sur une densité projetée corrigée point par point, reste insuffisante ;
    \item une branche optique construite à partir de l’accélération effective BuP est plus prometteuse ;
    \item il semble exister un \emph{régime optique BuP centralisé} quasi stable sur plusieurs spirales massives ;
    \item cette stabilité pourrait refléter une dépendance par classe galactique, en cohérence avec la structure multi-régimes déjà suggérée par l’analyse SPARC des courbes de rotation.
\end{enumerate}

\subsubsection{Limites}

Ces résultats doivent être interprétés comme une \emph{preuve de concept optique} et non comme une théorie complète du lentillage. En particulier :
\begin{itemize}
    \item la métrique optique effective n’est pas encore dérivée de premiers principes à partir de l’intrication ;
    \item la forme du poids radial \(w_{\mathrm{tot}}(r)\) reste phénoménologique ;
    \item la comparaison a été effectuée sur un petit nombre de galaxies ;
    \item l’extension aux lentilles fortes de type SLACS n’est pas encore établie.
\end{itemize}

Malgré ces limites, cette étape est importante : elle montre qu’un cadre bottom-up fondé sur la dimension émergente peut produire un excès de déflexion optique contrôlé, centralisé, et potentiellement stable sur une classe astrophysique donnée.

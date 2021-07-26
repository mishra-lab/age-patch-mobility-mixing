\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{methodsx}[2021]
% article inheritance
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
% packages
\RequirePackage[margin=3cm]{geometry}
\RequirePackage{authblk}
% references
\usepackage[
  style=numeric,
  doi=false,
  url=false,
  isbn=false,
  natbib=true,
  backend=biber,
  sorting=none,
  maxbibnames=3,
  minbibnames=3,
  giveninits=true,
]{biblatex}
% spacing
\setlength{\parskip}{1ex}
\setlength{\parindent}{0ex}
\setlength{\abovecaptionskip}{\medskipamount}
\setlength{\belowcaptionskip}{\medskipamount}
\widowpenalty10000
\clubpenalty10000
% meta-data
\newcommand{\keywords}[1]{\gdef\@keywords{#1}}
\newcommand{\submission}[1]{\gdef\@submission{#1}}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\newcommand{\grabstract}[1]{\gdef\@grabstract{#1}}
\newcommand{\area}[1]{\gdef\@area{#1}}
\newcommand{\subject}[1]{\gdef\@subject{#1}}
\newcommand{\method}[1]{\gdef\@method{#1}}
\newcommand{\resource}[1]{\gdef\@resource{#1}}
\newcommand{\funding}[1]{\gdef\@funding{#1}}
\newcommand{\acknowl}[1]{\gdef\@acknowl{#1}}
\newcommand{\contrib}[1]{\gdef\@contrib{#1}}
% helpers
\renewcommand{\paragraph}[1]{\par\textbf{#1:}}
\newcommand{\headstyle}[1]{\normalfont\textsc{\MakeLowercase{#1}}}
\newcommand{\tpart}[3][:]{%
  \hspace{0.14\linewidth}\llap{\textbf{#2#1}}\hspace{0.01\linewidth}%
  \noindent\parbox[t]{0.85\textwidth}{#3}\par\medskip}
\newcommand{\apart}[3][:]{%
  \textbf{#2#1}\par\smallskip{#3}\par\medskip}
% maketitle
\renewcommand{\maketitle}{%
  \thispagestyle{empty}
  \tpart{Title}{\@title}
  \tpart{Authors}{\@author}
  \tpart{Keywords}{\@keywords}
  \tpart{Journal}{MethodsX (preprint)}
  \tpart{Submission}{\@submission}
%  \tpart{Abstract}{\@abstract}
  \tpart{Subject Area}{\@area: \@subject}
  \tpart{Method}{\@method}
  \tpart{Resource}{\@resource}
  \tpart{Date}{\@date}
  \clearpage
  \apart{Abstract}{\@abstract}
  \apart{Graphical Abstract}{\frame{\@grabstract}}
  \clearpage\pagestyle{plain}}
\newcommand{\backmatter}{%
  \clearpage
  \subsection*{Funding}\@funding
  \subsection*{Acknowledgements}\@acknowl
  \subsection*{Contributions}\@contrib
  \clearpage
}
% appendix
\let\xappendix\appendix
\renewcommand{\appendix}{
  \xappendix
  \clearpage
  \setcounter{page}{1}
  \renewcommand*{\thepage}{A.\arabic{page}}
  \counterwithin{table}{section}
  \counterwithin{figure}{section}
  \counterwithin{equation}{section}
  \thispagestyle{plain}
  \hspace{0.15\linewidth}{\large\textsc{appendix}}
  \vskip1em
  \tpart{Title}{\@title}
  \tpart{Authors}{\@author}
  \tpart{Resource}{\@resource}
  \tpart{Date}{\@date}
  \refsection
  \clearpage
}

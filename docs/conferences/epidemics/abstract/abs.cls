\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{abs}[2021]
% article inheritance
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
% packages
\RequirePackage[margin=3cm]{geometry}
\RequirePackage{authblk}
% spacing
\setlength{\parskip}{1ex}
\setlength{\parindent}{0ex}
\def\@parboxrestore{\setlength{\parskip}{1ex}}
\frenchspacing
% meta-data
\renewcommand{\title}[1]{\gdef\@title{#1}}
\newcommand{\conference}[1]{\gdef\@conference{#1}}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\newcommand{\maxwords}[1]{\gdef\@maxwords{#1}}
% helpers
\renewcommand{\section}[1]{\par\MakeUppercase{#1.}}
\newcommand{\tstyle}[1]{\textbf{#1}}
\newcommand{\tpart}[3][:]{%
  \hspace{0.14\linewidth}\llap{\tstyle{#2#1}}\hspace{0.01\linewidth}%
  \noindent\parbox[t]{0.85\textwidth}{#3}\par\medskip}
% maketitle
\renewcommand{\maketitle}{%
  \immediate\write18{detex \@abstract.tex | wc -w > words.tex}
  \thispagestyle{empty}
  \tpart{Title}{\@title}
  \tpart{Authors}{\@author}
  \tpart{Conference}{\@conference\ (preprint)}
  \tpart{Date}{\@date}
  \tpart{Abstract}{\input{\@abstract}\par\bigskip
    Words:~\input{words}\unskip\,/\,\@maxwords}}

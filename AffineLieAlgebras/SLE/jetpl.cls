%%% JETPL -- February 2001 -- version for outer use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{jetpl}[2001/02/08 v057]
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{multicol}
\RequirePackage{tabularx}
\RequirePackage{euscript}
\RequirePackage{calc}
\RequirePackage{indentfirst}

%%%----------MAKET MACROS----------
\newcommand{\sodauthor}[1]{\global\def\@sodauthor{#1}}
\newcommand{\sodtitle}[1]{\global\def\@sodtitle{#1}}
\newcounter{issue}
\newcounter{volume}
\newcommand{\newissue}{%
  \if@maket%
    \addtocounter{issue}{1}
    \addtocontents{sod}{\protect\vfill\protect\eject}
    \addtocontents{sod}{\protect\begin{center}}
    \addtocontents{sod}{\protect\Large\protect\bfseries}
    \addtocontents{sod}{\protect Содержание\\[2mm]}
    \addtocontents{sod}{\protect\large}
    \addtocontents{sod}{\protect Выпуск\string~\theissue}
    \addtocontents{sod}{\protect\end{center}}
  \fi
  \vfill\eject
  \thispagestyle{empty}\vspace*{30mm}\hbox{}
  \cleartwopages
}
\newcommand{\newtheme}[1]{%
  \if@maket%
    \addtocontents{sod}{\protect\begin{center}}
    \addtocontents{sod}{\string{\string\bfseries\string\normalsize\ #1\string}}
    \addtocontents{sod}{\protect\end{center}}
  \fi
}
\newif\if@maket\@maketfalse
\newcounter{countart}
\setcounter{countart}{0}
\newif\if@vypusk
\newif\if@vypfp
\def\maket{\@makettrue\@vypuskfalse\@vypfpfalse}

\AtBeginDocument{%
  \setcounter{FirstPageNumber}{\value{page}}%
  \newcommand{\inputarticle}[1]{%
    \addtocounter{countart}{1}
    \def\articlename{#1}%
    \input{\theissue#1}%
    \label{\articlename:lastpagenumber}%
  }

  \if@maket%
    \addtocontents{sod}{\protect\onecolumn}
    \addtocontents{sod}{\protect\rus}
    \addtocontents{sod}{\protect\thispagestyle{plain}}
    \addtocontents{sod}{\protect\begin{center}}
    \addtocontents{sod}{\protect\Large\protect\bfseries}
    \addtocontents{sod}{\protect Содержание\\[2mm]}
    \addtocontents{sod}{\protect\large}
    \addtocontents{sod}{\protect Выпуск\string~\theissue}
    \addtocontents{sod}{\protect\end{center}}
  \fi
}

\AtEndDocument{\if@maket\@starttoc{sod}\fi}

\newcounter{mylang}
\def\rus{\setcounter{mylang}{1}\ifx\l@russian\@undefined\else
  \language\l@russian\fi\lefthyphenmin2\righthyphenmin2}
\def\lat{\setcounter{mylang}{0}\ifx\l@english\@undefined\else
  \language\l@english\fi\lefthyphenmin2\righthyphenmin3}
\newcommand{\thefigure}{\arabic{figure}}

\newcommand\@ptsize{}

\newcommand\@firstpage{}
\def\firstpage#1{\edef\@firstpage{#1}\c@page=\@firstpage}
\newcommand\@lastpage{}
\def\lastpage#1{\gdef\@lastpage{#1}}

\newif\if@restonecol
\newif\if@titlepage
\newif\if@mcol
\newif\if@wideon
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}{\setlength\paperheight{297mm}\setlength\paperwidth{210mm}}
\DeclareOption{landscape}
  {\setlength\@tempdima{\paperheight}%
  \setlength\paperheight{\paperwidth}%
  \setlength\paperwidth{\@tempdima}}
\fi
\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility\else
  \DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
  \DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
  \DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
  \DeclareOption{onecolumn}{\@twocolumnfalse\@mcolfalse\@wideontrue}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue\@mcolfalse\@wideontrue}
\DeclareOption{mcolumn}{\@mcoltrue\@twocolumnfalse}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
    \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent-\bibindent
      \listparindent\itemindent
      \parsep\z@
    }%
    \renewcommand\newblock{\par}%
  }%
}
\ExecuteOptions{a4paper,10pt,twoside,onecolumn,final}
\ProcessOptions

\input{size1\@ptsize.clo}
\AtEndOfClass{%
  \newcommand{\twocol}{\twocolumn%
    \RequirePackage{ftnright}%
    \renewcommand\footnoterule{\vspace*{-3pt}%
    \noindent\rule{.4\columnwidth}{.4pt}\vspace*{2.6pt}}%
  }%
}
\setlength\lineskip{1\p@}
\renewcommand\baselinestretch{1.05}
\setlength\parskip{0\p@\@plus.3pt}
\@lowpenalty51
\@medpenalty151
\@highpenalty301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.9}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.7}
\hoffset-25.4mm
\voffset-25.4mm
\topmargin24mm
\headsep6mm
\headheight12pt
\marginparwidth0pt
\oddsidemargin0pt
\marginparsep0pt
\footskip9mm
\footnotesep1mm
%\setlength{\skip\footins}{15\p@\@plus4\p@\@minus2\p@}
\textwidth172.3mm
\AtBeginDocument{\addtolength{\columnsep}{8pt}}
\AtBeginDocument{\addtolength{\textwidth}{10pt}}
\textheight230.3mm
\tabcolsep5pt
\oddsidemargin13.5mm
\evensidemargin13.5mm
\newlength\mycolwidth
\setlength{\mycolwidth}{\textwidth}
\addtolength{\mycolwidth}{-\columnsep}
\divide\mycolwidth by\tw@
\@addtoreset{footnote}{part}
\newskip\minskip
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength{\minskip}{1.5mm plus.15fil minus.5mm}
\setlength{\abovecaptionskip}{\minskip}
\setlength{\belowcaptionskip}{\minskip}
\setlength{\skip\footins}{.4mm plus.1fil minus.1mm}
\setlength{\abovedisplayskip}{\minskip}
\setlength{\belowdisplayskip}{\minskip}
\setlength{\abovedisplayshortskip}{\minskip}
\setlength{\belowdisplayshortskip}{\minskip}
\newcounter{FirstPageNumber}

%%%----------MULTICOLS----------
\AtEndDocument{\if@wideon\else\end{multicols}\fi}
\AtEndOfClass{%
  \if@mcol%
    \RequirePackage{multicol}\global\@wideonfalse
    \multicolsep0pt   % before/after skips
    \premulticols50pt % free space at the left of a page
    \postmulticols20pt%
    \newlength\colwidth{}\colwidth\textwidth
    \divide\colwidth\tw@
  \fi
}
\newcommand{\rulka}{}
\AtBeginDocument{%
  \renewcommand{\rulka}{\rightline{\rule[-3mm]{.1pt}{3mm}\rule{\colwidth}{.5pt}}}
  \def\wideonx{%
    \if@twocolumn%
    \else
      \if@mcol
        \if@wideon\else\global\@wideontrue%
          \par\end{multicols}\par\noindent
        \fi
      \fi
    \fi
  }
  \def\wideon{%
    \if@twocolumn%
    \else
      \if@mcol
        \if@wideon\else\global\@wideontrue%
          \par\end{multicols}%
          \leftline{\rule{\colwidth}{.1pt}\rule{.1pt}{2mm}}\par\noindent
        \fi
      \fi
    \fi
  }
  \def\wideoffx{%
    \if@twocolumn%
    \else
      \if@mcol
        \if@wideon\global\@wideonfalse%
          \par\begin{multicols}{2}[~]\noindent
        \fi
      \fi
    \fi
  }
  \def\wideoff{%
    \if@twocolumn%
    \else
      \if@mcol
        \if@wideon\global\@wideonfalse%
          \rightline{%
          \rule[-4mm]{.1pt}{2mm}\rule[-2mm]{\colwidth}{.5pt}}%
          \par\begin{multicols}{2}\noindent
        \fi
      \fi
    \fi
  }%
}

\newcommand{\eqalign}{\split}

%%%----------TITLEPAGE----------
\newcommand\@abstract{}
\newcommand{\abstract}[1]{\gdef\@abstract{\parbox{155mm}{%
  \hspace*{5.5mm}\small\ignorespaces#1}}%
}

\newcommand\@rtitle{}
\newcommand{\rtitle}[1]{\gdef\@rtitle{#1}}
\let\@rtitle=\@title
\newcommand\@rauthor{}
\newcommand\@address{}
\newcommand\@authoraddress{}
\renewcommand{\author}[1]{%
  \expandafter\def\expandafter\@authoraddress\expandafter
  {\@authoraddress%
    {\dimen0=-\prevdepth\advance\dimen0 by8mm%
      \nointerlineskip\centering%
      \vrule height\dimen0 width0pt depth3mm\relax\ignorespaces{%
      \slshape#1}\unskip\par
    }%
  }%
}
\let\@rauthor=\@authoraddress
\def\rauthor#1{\gdef\@rauthor{#1}}

\def\@PACS{}
\def\PACS#1{\gdef\@PACS{#1}}

\def\address#1{\expandafter\def\expandafter\@authoraddress\expandafter
  {\@authoraddress{\footnotesize\sl\centering\ignorespaces#1\unskip\par}}}

\newcommand{\submitted}{%
  {\ifnum\value{mylang}=0 Submitted\ \else Поступила в редакцию\ \fi}}
\newcommand{\resubmitted}{%
  {\ifnum\value{mylang}=0 Resubmitted\ \else После переработки\ \fi}}

\newcommand\@dates{}
\newcommand{\dates}[2]{\gdef\@dates{{\footnotesize\centering
  \submitted#1{\if#2*\else\par\resubmitted#2\fi}}\par
}}

\newcommand\@JETPLyear[1]{%                       Вычисляем год по тому
  \ifcase#1%
       \or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or%  1--20
    \or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or% 21--40
    \or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or\or% 41--60
    \or\or\or\or\or\or\or\or\or\or\or\or\or% 61--72
    \or2001\or2001% 73--74
    \or2002\or2002% 75--76
    \or2003\or2003% 77--78
    \or2004\or2004% 79--80
    \or2005\or2005% 81--82
    \or2006\or2006% 83--84
    \or2007\or2007% 85--86
    \else{\@latexerr{Volume \arabic{volume}~inconsistent}\@eha}%
  \fi
}

\newcommand{\rusissmonth}[1]{%
  \ifcase\number#1\or
    января\or января\or февраля\or февраля\or марта\or марта\or
    апреля\or апреля\or мая\or мая\or июня\or июня\or июля\or июля\or
    августа\or августа\or сентября\or сентября\or октября\or октября\or
    ноября\or ноября\or декабря\or декабря\fi
}
\newcommand{\engissmonth}[1]{%
  \ifcase\number#1\or
    January\or January\or February\or February\or March\or March\or
    April\or April\or May\or May\or June\or June\or July\or July\or
    August\or August\or September\or September\or October\or October\or
    November\or November\or December\or December\fi
}

\newcommand\ps@titlehead{%  Стиль titlehead первой полосы статьи
  \renewcommand\@evenhead{%
    \hbox to\textwidth{\mdseries\rmfamily\upshape\small
      \ifnum\value{mylang}=0
        Pis'ma v ZhETF%
        \if@maket%
          , vol.\,\arabic{volume}, iss.\,\arabic{issue},
          pp.\,\pageref{\articlename:firstpagenumber}\,--\,%
            \pageref{\articlename:lastpagenumber}
          \hfill\copyright\ \@JETPLyear{\value{volume}}\
          \engissmonth{\value{issue}}\ \ifodd\value{issue}10\else25\fi
        \else
          \hfill
        \fi
      \else
        Письма в ЖЭТФ%
        \if@maket%
          , том \arabic{volume}, вып.\,\arabic{issue},
          с.\,\pageref{\articlename:firstpagenumber}\,--\,%
            \pageref{\articlename:lastpagenumber}
          \hfill\copyright\ \@JETPLyear{\value{volume}}\,г.
          \ifodd\value{issue}10\else25\fi\ \rusissmonth{\value{issue}}%
        \else
          \hfill
        \fi
      \fi
    }
  }
  \let\@oddhead\@evenhead
  \renewcommand\@oddfoot{%
    \def\theList{\arabic{List}}%
    \setcounter{List}{1+(\value{page}-\value{FirstPageNumber}+3)/16}%
    \setcounter{PageInList}{(\value{page}-\value{FirstPageNumber}+3)-16*\value{List}+16}%
    \if@maket%
      \ifnum\value{List}=1%
        \rlap{\mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
        том \arabic{volume}\quad
        вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}}\hfill
        \mdseries\rmfamily\upshape\normalsize\arabic{page}\hfill
      \else
        \ifnum\value{PageInList}=1%
          \rlap{\bfseries\rmfamily\upshape\small\arabic{List}\hspace*{5mm}%
          \mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
          том \arabic{volume}\quad
          вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}}\hfill\arabic{page}\hfill
        \else
          \ifnum\value{PageInList}=3%
            \rlap{\mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
            том \arabic{volume}\quad
            вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}}\hfill\arabic{page}\hfill
            \bfseries\rmfamily\upshape\small
            \llap{\arabic{List}$\mathsurround0mm^*$}%
          \else
            \rlap{\mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
            том \arabic{volume}\quad
            вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}}\hfill
            \mdseries\rmfamily\upshape\normalsize\arabic{page}\hfill
          \fi
        \fi
      \fi
    \else
      \rlap{\mdseries\rmfamily\upshape\small Письма в ЖЭТФ}\quad\hfill
      \mdseries\rmfamily\upshape\normalsize\arabic{page}\hfill
    \fi
  }
  \renewcommand{\@evenfoot}{%
    \hfill\mdseries\rmfamily\upshape\normalsize\arabic{page}\hfill
    \llap{\mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
    том \arabic{volume}\quad
    вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}}%
  }%
}

\def\@maketitle{%
  \vbox{\vspace*{-4mm}
    \begin{center}
    {\Large\rmfamily\bfseries\upshape\@title\par\vspace*{-4mm}}%
    \vskip12pt plus.1fil minus6pt%
    \@authoraddress\par
    \vspace*{-3mm}
    \vskip6pt plus.1fil minus3pt%
    \vspace*{3mm}%
    \@dates\par
    \vskip2em plus.1fil minus1ex%
    \@abstract%
    \vskip2em plus1em minus1ex%
    \par
  \end{center}
  }%
}

\newcommand\maketitle{\par
  \setcounter{footnote}{0}%
  \begingroup
    \def\@makefnmark{$\mathsurround0mm^{\number\c@footnote)}$}
    \long\def\@makefntext##1{\parindent1em\noindent
      \hb@xt@1.8em{%
        \hss$\mathsurround0mm^{\number\c@footnote)}$}##1%
      }%
    \if@twocolumn
      \twocolumn[\@maketitle]%
    \else
      \global\@topnum\z@
      \@maketitle%\@thanks
    \fi
    \thispagestyle{titlehead}\@thanks
  \endgroup
  \if@maket%
    \label{\articlename:firstpagenumber}%
    \addtocontents{sod}{\string\def\string\art\articlename{\@sodtitle}}%
    \addtocontents{sod}{\string\def\string\aut\articlename{\@sodauthor}}%
    \addtocontents{sod}{\string\def\string\pag\articlename\string{\pageref{\articlename:firstpagenumber}\string}}%
    \addtocontents{sod}{\string\def\string\pacs\articlename{\@PACS}}%
    \addtocontents{sod}{\string\sodentry{\string\art\articlename}{\string\aut\articlename}{\string\pag\articlename}}%
  \fi
  \if@mcol\begin{multicols}{2}\fi
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@authoraddress\@empty
  \global\let\@dates\@empty
  \global\let\@title\@empty
  \global\let\@sodauthor\@empty
  \global\let\@sodtitle\@empty
  \setcounter{equation}{0}
}

%%%----------КОЛОНТИТУЛЫ----------
\newcounter{PageInList}% Номер страницы в печатном листе: 1,2,...,16
\newcounter{List}      % Номер печатного листа: 1,2,... (16 стр. в листе)
\newcounter{auxnumber}
\newbox\issuebox
\AtBeginDocument{%
  \setcounter{auxnumber}{\value{issue}}%
  \ifodd\value{issue}%
    \addtocounter{auxnumber}{1}%
    \setbox\issuebox\hbox{\arabic{issue}\,--\,\arabic{auxnumber}}%
  \else
    \addtocounter{auxnumber}{-1}%
    \setbox\issuebox\hbox{\arabic{auxnumber}\,--\,\arabic{issue}}%
  \fi
}
\def\ps@headings{%
  \def\@evenfoot{\hfill\mdseries\rmfamily\upshape\small
    Письма в ЖЭТФ\quad
    \if@maket%
      том \arabic{volume}\quad
      вып.\,\usebox\issuebox\quad\@JETPLyear{\value{volume}}%
    \fi
  }
  \def\@oddfoot{%
    \def\theList{\arabic{List}}%
    \setcounter{List}{1+(\value{page}-\value{FirstPageNumber}+3)/16}%
    \setcounter{PageInList}{(\value{page}-\value{FirstPageNumber}+3)-16*\value{List}+16}%
    \ifnum\value{List}=1%
      \mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
      \if@maket%
        том \arabic{volume}\quad
        вып.\,\usebox\issuebox\quad\@JETPLyear{\value{volume}}\hfill
      \else
        \hfill
      \fi
    \else
      \ifnum\value{PageInList}=1%
        \bfseries\rmfamily\upshape\small\arabic{List}\hspace*{5mm}%
        \mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
        \if@maket%
          том \arabic{volume}\quad
          вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}\hfill
        \else
          \hfill
        \fi
      \else
        \ifnum\value{PageInList}=3%
          \mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
          \if@maket%
            том \arabic{volume}\quad
            вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}\hfill
            \bfseries\rmfamily\upshape\small
            \arabic{List}$\mathsurround0mm^*$
          \else
            \hfill
          \fi
        \else
          \mdseries\rmfamily\upshape\small Письма в ЖЭТФ\quad
          \if@maket%
            том \arabic{volume}\quad
            вып.\,\usebox{\issuebox}\quad\@JETPLyear{\value{volume}}\hfill
          \else
            \hfill
          \fi
        \fi
      \fi
    \fi
  }
  \if\c@page=@firstpage% \@titlehead
  \else
    \def\@oddhead{\vbox{%
      \hbox to\textwidth{\hfill\normalfont\slshape\@rtitle\hfill
      \upshape\arabic{page}}%
      \smallskip\hrule width\textwidth depth.2pt height0pt%
    }}%
    \def\@evenhead{\vbox{%
      \hbox to\textwidth{\normalfont\upshape\arabic{page}\hfill
      \slshape\@rauthor\hfill}%
      \smallskip\hrule width\textwidth depth.2pt height0pt%
    }}%
  \fi
}

\pagestyle{headings}

%%%----------Math----------
\binoppenalty700
\relpenalty500
\newcommand{\qqq}[1]{#1\nobreak\discretionary{}{\hbox{\m@th$#1$}}{}}
\def\arctg{\mathop{\text{arctg}}\nolimits}
\def\Sp{\mathop{\text{Sp}}\nolimits}

\DeclareRobustCommand\lambdabar{%
  \bgroup
    \def\@tempa{%
      \hbox{%
        \raise.73\ht\z@
        \hb@xt@\z@{%
           \kern.25\wd\z@
           \vrule \@width.5\wd\z@\@height.1\p@\@depth.1\p@
           \hss
        }%
        \box\z@
      }%
    }%
    \mathchoice
      {\setbox\z@\hbox{$\displaystyle     \lambda$}\@tempa}%
      {\setbox\z@\hbox{$\textstyle        \lambda$}\@tempa}%
      {\setbox\z@\hbox{$\scriptstyle      \lambda$}\@tempa}%
      {\setbox\z@\hbox{$\scriptscriptstyle\lambda$}\@tempa}%
  \egroup
}%

%%%----------MISC----------
\newcommand{\Frac}[2]{\displaystyle\frac{#1}{#2}}

\DeclareTextSymbolDefault{\textnumero}{OT1}
\def\@makefnmark{\hbox{\@textsuperscript{\normalfont\@thefnmark)}}}

\newcommand{\centerlast}[1]{
  \addtolength{\leftskip}{#1 plus1fil}%
  \addtolength{\rightskip}{#1 plus-1fil}%
  \setlength{\parfillskip}{#1 plus2fil}%  центровка посл. строки
}

\AtBeginDocument{\renewcommand{\~}{\ensuremath{\thinspace}}}

\AtBeginDocument{
  \pretolerance-1
  \tolerance1900
  \adjdemerits10000
  \emergencystretch10mm   % Кнут, с.131
  \clubpenalty10000       % штраф за отрыв первой строки
% \widowpenalty4000       % штраф за висячую строку (перенос посл. строки абзаца)
  \widowpenalty10000
% \displaywidowpenalty4000% то же, но перед выключной формулой
  \displaywidowpenalty10000
}
\hyphenation{час-тиц час-тицы час-ти-цах рос-том воз-буж-ден-ные
  воз-буж-ден-ны-ми воз-буж-ден-ное воз-буж-де-ния воз-буж-де-нию
  воз-буж-ден-но-го воз-буж-де-ни-я-ми воз-буж-де-ний клас-тер клас-те-ра
  клас-те-ров час-ти-цы час-ти-цам час-тиц час-то час-то-та час-то-ты час-то-той
  час-то-те час-тот вы-со-ко-час-тот-ных вы-со-ко-час-тот-ной
  вы-со-ко-час-тот-но-го вы-со-ко-час-тот-ную об-суж-де-ние
}

%%%----------REFERENCES----------
\AtBeginDocument{\newcommand\refname{}}
\renewcommand{\@biblabel}[1]{{#1}.}
\newenvironment{references}%
  {\reset@font\parsep0pt\itemsep0pt\begin{thebibliography}{99}\small}%
  {\end{thebibliography}}

%%%----------FIGURE----------
%%%---------- artcle.cls remnants ----------
\newcounter{figure}
\renewcommand{\bottomfraction}{.7}
\renewcommand{\topfraction}{1}
\renewcommand{\textfraction}{0}
\renewcommand{\dbltopfraction}{.95}
\renewcommand{\floatpagefraction}{.9}
\renewcommand{\dblfloatpagefraction}{.9}

\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
  {\@float{figure}}{\end@float}
\newenvironment{figure*}
  {\@dblfloat{figure}}{\end@dblfloat}
\RequirePackage[dvips]{graphicx}
\setlength{\textfloatsep}{15pt plus2pt minus4pt}

%%%----------\ris[NO]{FILENAME}{WIDTH}{HEIGHT}---centered----------
\newif\if@captlabel\@captlabeltrue
\newcommand{\@captype}{figure}

\newcommand\ris[5][\@arabic\c@figure]{%
  \renewcommand{\@captype}{figure}%
  \ifx*#1\@captlabelfalse\else\@captlabeltrue\fi
  \def\figurename{}%
%  \renewcommand\thefigure{{\small\unskip#1.~~}}%
  \renewcommand\thefigure{}%
  \ifdim#3<.5\textwidth
    \parbox[c]{.5\textwidth-\columnsep}{%
    \centerline{\includegraphics[bb=0 0 #3 #4]{#2.pcx}}%
      \caption{#5}\vskip3mm%
    }\hspace{\columnsep}%
  \else
    \centerline{\includegraphics[bb=0 0 #3 #4]{#2.pcx}}%
    \caption{#5}\vskip3mm%
  \fi
}

\newcommand{\fig}[3][\@arabic\c@figure]{%
  \renewcommand{\@captype}{figure}%
  \ifx*#1\@captlabelfalse\else\@captlabeltrue\fi
  \def\figurename{}
  \renewcommand{\thefigure}{}%
  \ifx*#2\else\setbox0=\hbox{\includegraphics{K:/РИСУНКИ/#2.eps}}\fi
%% \ifdim\wd0<.5\textwidth%
%%   \parbox[c]{.5\textwidth-\columnsep}{\noindent
%      \centerline{\box0}%
%      \caption{#3}\vskip3mm%
%%   }\hspace{\columnsep}%
%% \else
%%   \centerline{\box0}
%%   \caption{#3}\vskip3mm%
%% \fi
  \ifdim\wd0<.5\textwidth%
    \parbox[c]{.5\textwidth-\columnsep}{\noindent
      \centerline{\box0}%
      \caption{#3}\vskip3mm%
    }\hspace{\columnsep}%
  \else
    \centerline{\box0}
    \caption{#3}\vskip3mm%
  \fi
}

\newcommand\sideris[5][\@arabic\c@figure]{%
  \renewcommand{\@captype}{figure}%
  \ifx*#1\@captlabelfalse\else\@captlabeltrue\fi
  \def\figurename{}
  \renewcommand\thefigure{}%
  \newlength{\otstup}%
  \setlength{\otstup}{\columnwidth}%
  \addtolength{\otstup}{-#3}%\addtolength{\otstup}{-4mm}%
% \setlength{\otstup}{.5\otstup}%
  \parbox[c]{#3}{\includegraphics[bb=0 0 #3 #4]{#2.pcx}}%
% \kern4mm%
  \parbox[c]{\otstup}{\caption{#5}}%
}

\newlength{\otstup}%

\newcommand\oneris[8][\@arabic\c@figure]{%
  \renewcommand{\@captype}{figure}%
  \ifx*#1\@captlabelfalse\else\@captlabeltrue\fi
  \def\figurename{}%
  \renewcommand\thefigure{}%
  \setlength{\otstup}{.5\textwidth}\addtolength{\otstup}{-\columnsep}%
  \addtolength{\otstup}{-#3}%
  \hspace*{.5\otstup}%
  \parbox[c]{#3}{\includegraphics[bb=0 0 #3 #4]{#2.pcx}}%
  \hspace*{.5\otstup}%
  \setlength{\otstup}{.5\textwidth}\addtolength{\otstup}{-\columnsep}%
  \addtolength{\otstup}{-#6}%
  \hspace*{.5\otstup}%
  \parbox[c]{#6}{\includegraphics[bb=0 0 #6 #7]{#5.pcx}}%
  \par
  \caption{#8}%
}

%%%----------TABLE----------
%%%---------- artcle.cls remnants ----------
\newcounter{table}\setcounter{table}{1}
\renewcommand\thetable{\@arabic\c@table}
\renewcommand{\tabularxcolumn}[1]{m{#1}}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}

\RequirePackage{array}
\def\arraystretch{1.2}
\RequirePackage{dcolumn}
\newcolumntype{d}[1]{D{.}{.}{#1}}

\newcommand{\tabl}[2][\@arabic\c@table\addtocounter{table}{1}]{%
  \renewcommand{\@captype}{table}%
  \renewcommand\thetable{{\sffamily\small\bfseries\unskip#1.~~}}%
  \def\tablename{%
    \ifx*#1%
    \else
      \ifnum\value{mylang}=0
        {\small Table$~$\thetable~}%
      \else
        {\small Таблица$~$\thetable~}%
      \fi
    \fi
  }%
  \addtolength{\leftskip} {.06\columnwidth plus 1fil}% центровка
  \addtolength{\rightskip}{.06\columnwidth plus-1fil}%
  \setlength{\parfillskip}{.06\columnwidth plus 2fil}%
  {\tablename{\small#2}}
}

% Средство для вертикальной центровки в таблицах (PTEMAC)
\newcommand{\cvbox}[2]{%
  \setbox0=\hbox{#2}%
  \dimen0\ht\strutbox%
  \advance\dimen0+\dp\strutbox%
  \multiply\dimen0 by#1%
  \advance\dimen0-\ht0%
  \advance\dimen0-\dp0%
  \divide\dimen0 by2
  \setbox0=\hbox{\lower\dimen0\hbox{#2}}\dp0=0pt\box0%
}

\newenvironment{table}
  {%\@centering%
%       \renewcommand{\as@makecaption}[2]{%
%       \setlength{\leftskip}{\captionmargin}%
%       \setlength{\rightskip}{\captionmargin}%
%       \addtolength{\as@captionwidth}{-2\captionmargin}%
%       \captionfont%
%       \newsavebox{\as@mycaptionbox}
%     \as@caption{\parbox{\columnwidth}{ {\hfill Табл.\thetable\par#2}}
%              \par}}%
 \@float{table}}
 {\end@float}
\newenvironment{table*}
  {\@dblfloat{table}}
  {\end@dblfloat}

%%%----------CAPTIONS----------
%%%---------- artcle.cls remnants ----------
%\newlength\abovecaptionskip
%\newlength\belowcaptionskip
%\setlength\abovecaptionskip{10\p@ plus2pt minus4pt}
%\setlength\belowcaptionskip{0\p@ plus4pt minus4pt}
\newlength{\newcaptlength}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
%  \sbox\@tempboxa{#1: #2}%
  \sbox\@tempboxa{\small#2}%
  \ifdim\wd\@tempboxa>\hsize
    \setlength{\newcaptlength}{\hsize}%
    \addtolength{\newcaptlength}{-.08\columnwidth}%
    \hbox to\hsize{\hfill\parbox{\newcaptlength}{\small#2}\hfill}%
  \else
    \global\@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
%  \ifdim \wd\@tempboxa >\hsize
%   #1: #2%\par
% \else
%  \global \@minipagefalse
%    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%  \fi
% \addtolength{\leftskip}{.06\columnwidth plus1fil}%
% \addtolength{\rightskip}{.02\columnwidth plus-1fil}%
% \setlength{\parfillskip}{0pt plus1fill}%
% {\small#2\hfill\par}
  \vskip\belowcaptionskip
}

%%%----------OLD STYLE----------
\newcount\@indentflag\global\@indentflag=1
\newcommand\@indentamount{%
  \ifcase\@indentflag 0pt\or\@centering\or0pt plus1fil\fi\relax
}
\newcommand\inc@eqnnum{\addtocounter{equation}{1}}
\newcommand\dec@eqnnum{\addtocounter{equation}{-1}}
\def\@eqnnum{\hbox{\reset@font\rm(\theequation)}}
\let\make@eqnnum=\@eqnnum
\def\eqnum#1{\dec@eqnnum\global\def\make@eqnnum{\reset@font\rm(#1)}%
  \def\@currentlabel{#1}%
}

\def\eqalign#1{\null\,\vcenter{\openup2\jot\m@th
  \ialign{\strut\hfill$\displaystyle{##}$&$\displaystyle{{}##}$\hfill
  \crcr#1\crcr}}\,%
}

\def\eqnarray{\par\vskip-\lastskip\vskip\abovedisplayskip
  \inc@eqnnum\let\@currentlabel=\theequation
  \global\@eqnswtrue\m@th
  \global\@eqcnt\z@
  \tabskip\@totalleftmargin\advance\tabskip by\@indentamount\let\\\@eqncr
  \halign to\hsize\bgroup\hskip\@centering
  $\displaystyle\tabskip\z@{##{}}$&\global\@eqcnt\@ne
  \hfill${{}##{}}$\hfill
  &\global\@eqcnt\tw@ $\displaystyle\tabskip\z@{##}$\hfil
  %\tabskip\@centering \if@eqnsw\phantom{\make@eqnnum\kern\@eqtoeqnum}\fi
  \tabskip\@centering \if@eqnsw\phantom{\make@eqnnum}\fi
  &\llap{##}\tabskip\z@\cr
}

\def\endeqnarray{%
  \@@eqncr\egroup
  \vskip\belowdisplayskip\noindent
  \dec@eqnnum\global\@indentflag=1
  \global\let\make@eqnnum\@eqnnum
  \global\@ignoretrue
}

\def\nonumber{\global\@eqnswfalse
  \def\label##1{\@revmess{error}{\string\label{##1} used
  with \string\nonumber\space before (\theequation)}}%
}

\def\endreferences{\end{thebibliography}}

%%%----------
\setcounter{secnumdepth}{3}
\newcounter{part}
\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
\newcounter{subparagraph}[paragraph]
\renewcommand\thepart{\@Roman\c@part.}
\renewcommand\thesection{\@arabic\c@section.\kern-6pt}
\renewcommand\thesubsection{\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection.}
\renewcommand\theparagraph{\thesubsubsection.\@arabic\c@paragraph.}
\renewcommand\thesubparagraph{\theparagraph.\@arabic\c@subparagraph.}
\newcommand\part{%
  \if@noskipsec\leavevmode\fi
  \par
  \addvspace{4ex}%
  \@afterindentfalse
  \secdef\@part\@spart
}

\def\@part[#1]#2{%
  \ifnum\c@secnumdepth>\m@ne
    \refstepcounter{part}%
    \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
  \else
    \addcontentsline{toc}{part}{#1}%
  \fi
  {\parindent\z@\raggedright
  \interlinepenalty\@M%
  \normalfont%
  \ifnum\c@secnumdepth>\m@ne
    \Large\bfseries\partname~\thepart
    \par\nobreak%
  \fi
  \huge\bfseries#2%
  \markboth{}{}\par}%
  \nobreak
  \vskip3ex%
  \@afterheading
}
\def\@spart#1{%
  {\parindent\z@\raggedright
    \interlinepenalty\@M%
    \normalfont
    \huge\bfseries#1\par
  }%
  \nobreak
  \vskip3ex%
  \@afterheading
}

\newcommand\section{%
  \suppressfloats[t]%
  \@startsection{section}{1}{\z@}%
  {-8mm\@plus-1ex\@minus-.2ex}%
  {2ex\@plus.3ex}%
  {\ignorespaces\normalfont\small\bfseries\centering\MakeUppercase}%
}

\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
  {-6mm\@plus-1ex\@minus-.2ex}%
  {2mm\@plus.2ex}%
  {\normalfont\bfseries\centering}%
}

\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
  {-3.25ex\@plus-1ex\@minus-.2ex}%
  {1.5ex\@plus.2ex}%
  {\normalfont\small\bfseries\MakeUppercase}%
}

\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
  {3.25ex\@plus1ex\@minus.2ex}%
  {-1em}%
  {\normalfont\bfseries}%
}

\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
  {3.25ex\@plus1ex\@minus.2ex}%
  {-1em}%
  {\normalfont\small\bfseries}%
}

\if@twocolumn
  \setlength\leftmargini{2em}
\else
  \setlength\leftmargini{2.5em}
\fi
\leftmargin\leftmargini
\setlength\leftmarginii{2.2em}
\setlength\leftmarginiii{1.87em}
\setlength\leftmarginiv{1.7em}
\if@twocolumn
  \setlength\leftmarginv{.5em}
  \setlength\leftmarginvi{.5em}
\else
  \setlength\leftmarginv{1em}
  \setlength\leftmarginvi{1em}
\fi
\setlength\labelsep{.5em}
\setlength\labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty-\@lowpenalty
\@endparpenalty-\@lowpenalty
\@itempenalty-\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
  {\list{}{\labelwidth\z@\itemindent-\leftmargin\let\makelabel\descriptionlabel}}
  {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep\normalfont\bfseries#1}

\newenvironment{quotation}
  {\list{}{\listparindent1.5em%
    \itemindent\listparindent
    \rightmargin\leftmargin
    \parsep\z@\@plus\p@}%
    \item\relax%
  }
  {\endlist}
\newenvironment{quote}
  {\list{}{\rightmargin\leftmargin}\item\relax}{\endlist}
\if@compatibility
  \newenvironment{titlepage}{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse\newpage
    \fi
    \thispagestyle{empty}%
    \setcounter{page}\z@
  }
  {\if@restonecol\twocolumn\else\newpage\fi}
\else
  \newenvironment{titlepage}
    {\if@twocolumn
      \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }
    {\if@restonecol\twocolumn\else\newpage\fi
      \if@twoside\else
        \setcounter{page}\@ne
      \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \gdef\thesection{\@Alph\c@section.}%
}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins=\skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand\theequation{\@arabic\c@equation}

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
  \section*{\contentsname
    \@mkboth{\MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
  }
\newcommand*\l@part[2]{%
  \ifnum\c@tocdepth>-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em\@plus\p@}%
    \begingroup
      \parindent\z@\rightskip\@pnumwidth
      \parfillskip-\@pnumwidth
      {\leavevmode
        \large\bfseries#1\hfil\hb@xt@\@pnumwidth{\hss#2}%
      }\par
      \nobreak
      \if@compatibility
        \global\@nobreaktrue
        \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi
}
\newcommand*\l@section[2]{%
  \ifnum\c@tocdepth>\z@
    \addpenalty\@secpenalty
    \addvspace{1em\@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent\z@\rightskip\@pnumwidth
      \parfillskip-\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip-\leftskip
      #1\nobreak\hfil\nobreak\hb@xt@\@pnumwidth{\hss#2}\par
    \endgroup
  \fi
}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
  \section*{\listfigurename
    \@mkboth{\MakeUppercase\listfigurename}{\MakeUppercase\listfigurename}%
  }%
  \@starttoc{lof}%
}
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
  \section*{\listtablename
    \@mkboth{\MakeUppercase\listtablename}{\MakeUppercase\listtablename}
  }%
  \@starttoc{lot}%
}
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
  {%\section*
    \refname\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
    \list{\@biblabel{\@arabic\c@enumiv}}%
    {\settowidth\labelwidth{\@biblabel{#1}}%
      \leftmargin\labelwidth
      \advance\leftmargin\labelsep
      \@openbib@code
      \usecounter{enumiv}%
      \let\p@enumiv\@empty
      \setlength{\itemsep}{-.9mm}
      \renewcommand\theenumiv{\@arabic\c@enumiv}%
    }%
    \clubpenalty4000
    \@clubpenalty\clubpenalty
    \widowpenalty4000%
    \sfcode`\.\@m%
    \vskip.3\minskip
    \sloppy\small
    \hbox to\mycolwidth{\hspace*{1mm}\hfill\rule{25mm}{.4pt}\hfill\hspace*{1mm}}
  }
  {\def\@noitemerr
    {\@latex@warning{Empty `thebibliography' environment}}%
  \endlist%
  }
\newcommand\newblock{\hskip.11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
  {\if@twocolumn
     \@restonecolfalse
   \else
     \@restonecoltrue
   \fi
   \columnseprule\z@
   \columnsep35\p@
   \twocolumn[\section*{\indexname}]%
   \@mkboth{\MakeUppercase\indexname}{\MakeUppercase\indexname}%
   \thispagestyle{plain}\parindent\z@
   \parskip\z@\@plus.3\p@\relax
   \let\item\@idxitem
  }
  {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent40\p@}
\newcommand\subitem{\@idxitem\hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem\hspace*{30\p@}}
\newcommand\indexspace{\par\vskip10\p@\@plus5\p@\@minus3\p@\relax}

\renewcommand\footnoterule{%
  \kern-6\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@%
}
\newcommand\@makefntext[1]{%
  \parindent1em%
  \noindent
  \hb@xt@1.8em{\hss\@makefnmark}#1%
}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}

\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi

%%%----------CALLIGRAFIC MATH----------
\newcommand\cA{\ensuremath\mathcal{A}}\newcommand\cB{\ensuremath\mathcal{B}}
\newcommand\cC{\ensuremath\mathcal{C}}\newcommand\cD{\ensuremath\mathcal{D}}
\newcommand\cE{\ensuremath\mathcal{E}}\newcommand\cF{\ensuremath\mathcal{F}}
\newcommand\cG{\ensuremath\mathcal{G}}\newcommand\cH{\ensuremath\mathcal{H}}
\newcommand\cI{\ensuremath\mathcal{I}}\newcommand\cJ{\ensuremath\mathcal{J}}
\newcommand\cK{\ensuremath\mathcal{K}}\newcommand\cL{\ensuremath\mathcal{L}}
\newcommand\cM{\ensuremath\mathcal{M}}\newcommand\cN{\ensuremath\mathcal{N}}
\newcommand\cO{\ensuremath\mathcal{O}}\newcommand\cP{\ensuremath\mathcal{P}}
\newcommand\cQ{\ensuremath\mathcal{Q}}\newcommand\cR{\ensuremath\mathcal{R}}
\newcommand\cS{\ensuremath\mathcal{S}}\newcommand\cT{\ensuremath\mathcal{T}}
\newcommand\cU{\ensuremath\mathcal{U}}\newcommand\cV{\ensuremath\mathcal{V}}
\newcommand\cW{\ensuremath\mathcal{W}}\newcommand\cX{\ensuremath\mathcal{X}}
\newcommand\cY{\ensuremath\mathcal{Y}}\newcommand\cZ{\ensuremath\mathcal{Z}}

%%%----------EULER SCRIPT MATH---euscripr.sty----------
\newcommand\sA{\ensuremath\EuScript{A}}\newcommand\sB{\ensuremath\EuScript{B}}
\newcommand\sC{\ensuremath\EuScript{C}}\newcommand\sD{\ensuremath\EuScript{D}}
\newcommand\sE{\ensuremath\EuScript{E}}\newcommand\sF{\ensuremath\EuScript{F}}
\newcommand\sG{\ensuremath\EuScript{G}}\newcommand\sH{\ensuremath\EuScript{H}}
\newcommand\sI{\ensuremath\EuScript{I}}\newcommand\sJ{\ensuremath\EuScript{J}}
\newcommand\sK{\ensuremath\EuScript{K}}\newcommand\sL{\ensuremath\EuScript{L}}
\newcommand\sM{\ensuremath\EuScript{M}}\newcommand\sN{\ensuremath\EuScript{N}}
\newcommand\sO{\ensuremath\EuScript{O}}\newcommand\sP{\ensuremath\EuScript{P}}
\newcommand\sQ{\ensuremath\EuScript{Q}}\newcommand\sR{\ensuremath\EuScript{R}}
\newcommand\sS{\ensuremath\EuScript{S}}\newcommand\sT{\ensuremath\EuScript{T}}
\newcommand\sU{\ensuremath\EuScript{U}}\newcommand\sV{\ensuremath\EuScript{V}}
\newcommand\sW{\ensuremath\EuScript{W}}\newcommand\sX{\ensuremath\EuScript{X}}
\newcommand\sY{\ensuremath\EuScript{Y}}\newcommand\sZ{\ensuremath\EuScript{Z}}

%%%----------
\def\hm#1{#1\nobreak\discretionary{}{\hbox{\m@th$#1$}}{}}
\newcommand{\aglue}{\nobreak\hskip\z@skip}
\newcommand{\myglue}{\nobreak\hskip.167em plus.08em}
\DeclareRobustCommand*{\D}{\aglue\hbox{-}\aglue}
\DeclareRobustCommand*{\d}{\aglue\hbox{-}\hskip\z@skip}
\DeclareRobustCommand*{\CH}{\myglue\hbox{--}\myglue}
\DeclareRobustCommand*{\ch}{\myglue\hbox{--}\hskip.167em}
\DeclareRobustCommand*{\T}{\myglue\hbox{--}\myglue}
\DeclareRobustCommand*{\t}{\myglue\hbox{--}\hskip.167em}
\DeclareRobustCommand*{\Д}{\aglue\hbox{-}\aglue}
\DeclareRobustCommand*{\д}{\aglue\hbox{-}\hskip\z@skip}
\DeclareRobustCommand*{\Ч}{\myglue\hbox{--}\myglue}
\DeclareRobustCommand*{\ч}{\myglue\hbox{--}\hskip.167em}
\DeclareRobustCommand*{\Т}{\myglue\hbox{--}\myglue}
\DeclareRobustCommand*{\т}{\myglue\hbox{--}\hskip.167em}
\DeclareMathOperator{\Real}{Re}
\DeclareMathOperator{\Imag}{Im}
\def\cleartwopages{\clearpage\ifodd\c@page\else\hbox{}\thispagestyle{empty}%
  \newpage\if@twocolumn\hbox{}\thispagestyle{empty}\newpage\fi\fi}
\def\thanks{\@ifnextchar[\thanksa\thanksb}
\def\thanksa[#1]#2{\footnotemark\protected@xdef\@thanks{\@thanks
  \protect\setcounter{footnote}{#1}\protect\footnotetext{#2}}}
\def\thanksb#1{\footnotemark\protected@xdef\@thanks{\@thanks
  \protect\footnotetext[\the\c@footnote]{#1}}}

%%%----------Macros to write contents-----------
\newlength\toclength
\setlength{\toclength}{\textwidth}
\addtolength{\toclength}{-12mm}

\def\sodentry#1#2#3{
  \par\raggedright
  \parbox[b]{\toclength}{{\bfseries#2}\ #1 \dotfill\ }%
  \hspace*{4mm}%
  \parbox[b]{8mm}{#3}\vspace{4mm}%
}

\newif\if@filesod\@filesodfalse

\def\thlat{\lat\thart}
\def\thrus{\rus\thart}

\def\makesod{\if@maket\newwrite\@sodfile
  \immediate\openout\@sodfile=\vypnum sod
  \def\sod{\@bsphack\begingroup
    \def\protect####1####2{\string####1####2\space}\@sanitize
    \@wrsod
  }
  \def\thart{\@bsphack\begingroup
    \def\protect####1{\string####1\space}\@sanitize
    \@wrthart
  }\fi
}

\def\@wrsod#1#2{\let\thepage\relax
  \xdef\@gtempa{
    \immediate\write\@sodfile{\string\def\string\art\artname{#1}}
    \immediate\write\@sodfile{\string\def\string\aut\artname{#2}}
    \immediate\write\@sodfile{\string\def\string\pag\artname{\thepage}}%
  }
  \endgroup\@gtempa
  \if@nobreak\ifvmode\nobreak\fi\fi\@esphack
}

\def\pacs#1{\vspace*{1.5mm}
  PACS: #1
  \vspace{1.5mm}
  \if@maket\immediate\write\@sodfile{\string\def\string\pacs\artname{#1}}
  \fi
}

\def\@wrthart#1{\let\thepage\relax
  \xdef\@gtempb{%
    \immediate\write\@sodfile{\string
    \sodentry{\string\art#1}{\string\aut#1}{\string\pag#1}}%
  }
  \endgroup\@gtempb
  \if@nobreak\ifvmode\nobreak\fi\fi\@esphack
}

\def\sodcap{%
  \immediate\write\@sodfile{\string\input{jetplrus}}
  \immediate\write\@sodfile{\string\begin{document}}
  \immediate\write\@sodfile{}
  \immediate\write\@sodfile{\string\begin{center}\string\large\string\bf}
  \immediate\write\@sodfile{Вып.\vypnum}
  \immediate\write\@sodfile{\string\vspace{2mm}\string\par}
  \immediate\write\@sodfile{СОДЕРЖАНИЕ\string\end{center}}
  \immediate\write\@sodfile{}%
}

\def\theme#1{
 \immediate\write\@sodfile{}
 \immediate\write\@sodfile{\string\begin{center}}
 \immediate\write\@sodfile{{\noexpand\bf#1}}
 \immediate\write\@sodfile{\string\end{center}}
}

\def\endvypusk{
  \immediate\write\@sodfile{\string\end{document}}
  \immediate\closeout\@sodfile
  \clearpage
  \global\@vypuskfalse
}

\def\soderzanie{
  \if@vypusk\endvypusk\else\clearpage\fi
  \def\artname{sod}
  \input\number\vypnum sod.tex%
}

\def\vypusksod#1{
  \if@vypusk\endvypusk\else\clearpage\fi
  \def\artname{sod}
  \input\number#1sod.tex%
}

%%%----------Index macros----------
\newcount\avcount
\def\sa{0}
\def\NewLtS#1{\vglue.5cm plus.3cm minus.3cm}

\def\AventS#1#2#3#4{\def\sb{#4}
  \footnotesize\par#1 {\small\bfseries#2}, #3\ifx\sa\sb\else\hskip1mm(\sb)\fi\avcount1}

\def\AventfS#1#2#3#4{\def\sb{#4}
  \footnotesize #1 {\small\bfseries#2}, #3\ifx\sa\sb\else\hskip1mm(\sb)\fi\avcount1}

\def\SAventS#1#2#3{\footnotesize\def\sb{#3}
  \ifnum\avcount=1\hskip-1.5mm;\hfill\break\else\hskip-1.5mm; \fi
  {\small\bfseries #1},~#2\ifx\sa\sb\else\hskip1mm(\sb)\fi\avcount0}

% Новые команды для текущего авторского указателя
\def\AvtVypi#1#2#3#4{\def\sb{#4}
  \hangindent10mm\noindent
  \footnotesize\par#1 {\small\bfseries #2}, #3\ifx\sa\sb\else\,(\sb)\fi}
\def\AvtVypic#1#2#3#4{\AvtVypi{#1}{#2}{#3}{#4};}
\def\AvtVypiii#1#2#3#4{\hangindent10mm\noindent
   \footnotesize\def\sb{#4}
   {\small\bfseries #2},~#3\ifx\sa\sb\else\,(\sb)\fi}
\def\AvtVypii#1#2#3#4{\hfill\break\AvtVypiii{#1}{#2}{#3}{#4}}
\def\AvtVypiic#1#2#3#4{\AvtVypii{#1}{#2}{#3}{#4};}
\def\AvtVypiiic#1#2#3#4{\AvtVypiii{#1}{#2}{#3}{#4};}

\def\indexheading#1{
% \hbox to\pagewidth
\large\bfseries\hfill#1 \hfill\break\vskip5mm}

% Theme index macros
\newcommand{\NewTL}[2]{\section*{\large\bfseries#1 #2}}
\newcommand{\NewTM}[2]{\vspace{-1.2mm}\subsection*{\bfseries#1 #2}}
\newcommand{\NewTS}[2]{\vspace{-1.2mm}\subsubsection*{\bfseries#1 #2}}
\newcommand{\NewTT}[2]{\vspace{-1.2mm}\subsubsection*{\bfseries#1 #2}}
\newcommand{\ArtentL}[4]{\par\lat{~#3 -- #4 {\bfseries#1}, #2.}\pagebreak[0]}
\newcommand{\ArtentR}[4]{\par\rus{~#3 -- #4 {\bfseries#1}, #2.}\pagebreak[0]}
\newcommand{\Artent}[4]{\par{~#3 -- #4 {\bfseries#1},~#2.}\pagebreak[0]}

% Новые команды систематического указателя
\def\IndLargeHed#1#2{\vspace{-1mm}\NewTM{#1}{#2}\vspace{1mm}}
\def\IndSmallHed#1#2{\NewTS{#1}{#2}}
\def\IndArt#1#2#3#4{\Artent{#2}{#3}{#1}{#4}}

% Author index macros
\def\enpage#1{\def\sb{#1}%
\ifx\sa\sb\else\hskip1mm(\sb)\fi}
%\newcommand{\mpar}{\hspace*{5mm}}
\newcommand{\mpar}{\par\hspace*{5mm}}
\newcommand{\AventFB}[6]{\mpar{\bf #1} -- #5; #6 {\bf#2},~#3\enpage{#4}.}
\newcommand{\AventNB}[5]{\mpar{\bf #1} -- см.\ #5 {\bf#2},~#3\enpage{#4}.}
\newcommand{\SAventFB}[5]{\mpar\rule[.1cm]{3mm}{.03cm} #4; #5
{\bf#1},~#2\enpage{#3}.}
\newcommand{\SAventNB}[4]{\mpar\rule[.1cm]{3mm}{.03cm} см. #4
{\bf#1},~#2\enpage{#3}.}
\newcommand{\NewLtB}[1]{\vspace*{-2mm}\section*{\centerline{\large\bf #1}}}
%\newcommand{\NewLtB}[1]{}

\endinput
%%
%% End of file `jetpl.cls'

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{doublecol}
              [1998/05/05 v1.3y
 Standard LaTeX document class]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
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
\DeclareOption{draft}{\setlength\overfullrule{5pt}%
  \AtBeginDocument{\PagePositionCropYes
    {7pc}{5.3pc}{5.9pc}{7pt}
    {7in}{10in}{-1in}{-1in}}}
\DeclareOption{boxdraft}{\setlength\overfullrule{5pt}%
  \AtBeginDocument{\PagePositionBoxYes
    {7pc}{5.3pc}{5.9pc}{7pt}
    {7in}{10in}{-1in}{-1in}}}
\DeclareOption{final}{\setlength\overfullrule{0pt}%
  \AtBeginDocument{\PagePositionCropYes
    {2pc}{2pc}{3.77pc}{7pt}
    {6.3in}{9.5in}{-1in}{-1in}}}
\DeclareOption{finalcro}{\setlength\overfullrule{0pt}%
  \AtBeginDocument{\PagePositionCropYes
    {7pc}{5.3pc}{5.9pc}{7pt}
    {7in}{10in}{-1in}{-1in}}}
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
\ExecuteOptions{letterpaper,10pt,twoside,twocolumn,final}
\ProcessOptions
\input{size1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
%\setcounter{topnumber}{5}
%\def\topfraction{1}
%\def\textfraction{0}
%\def\floatpagefraction{0.9}
%\setcounter{bottomnumber}{5}
%\def\bottomfraction{0.99}
%\setcounter{totalnumber}{10}
%\def\dbltopfraction{0.99}
%\def\dblfloatpagefraction{0.8}
%\setcounter{dbltopnumber}{5}
\setcounter{topnumber}{9}
\renewcommand\topfraction{1}
\setcounter{bottomnumber}{9}
\renewcommand\bottomfraction{.1}
\setcounter{totalnumber}{9}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{9}
\renewcommand\dbltopfraction{.1}
\renewcommand\dblfloatpagefraction{.1}
\if@twoside
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\slshape\leftmark}%
      \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\MakeUppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\quad
        \fi
        ##1}}}
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markright {\MakeUppercase{%
        \ifnum \c@secnumdepth >\m@ne
          \thesection\quad
        \fi
        ##1}}}}
\fi
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }
  \if@titlepage
  \newcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \null\vfil
  \vskip 60\p@
  \begin{center}%
    {\LARGE \@title \par}%
    \vskip 3em%
    {\large
     \lineskip .75em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
      \vskip 1.5em%
    {\large \@date \par}%       % Set date in \large size.
  \end{center}\par
  \@thanks
  \vfil\null
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\else
\newcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\LARGE \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    \vskip 1em%
    {\large \@date}%
  \end{center}%
  \par
  \vskip 1.5em}
\fi
\setcounter{secnumdepth}{3}

\def\LRH#1{\gdef\setLRH{\sc #1}}\def\setLRH{}%
\def\RRH#1{\gdef\setRRH{\sc #1}}\def\setRRH{}%

\LRH{author}
\RRH{short title}

\def\ps@headings{\let\@mkboth\markboth
  \def\@oddfoot{\setoddFOOT}%
  \def\@evenfoot{\setevenFOOT}%
  \def\@evenhead{\setevenRH}%
  \def\@oddhead{\setoddRH}%
  \def\sectionmark##1{}%
  \def\subsectionmark##1{}}%

\def\JOURNALNAME#1{\gdef\theJOURNALNAME{#1}}%
\def\VOL#1{\gdef\theVOL{#1}}%
\VOL{0}%
\def\PAGES#1{\gdef\thePAGES{#1}}%
\PAGES{\thepage--\pageref{LastPage}}%
\def\PUBYEAR#1{\gdef\thePUBYEAR{#1}}%
\PUBYEAR{2003}%
\def\ISSUE#1{\gdef\theISSUE{#1}}%
\ISSUE{00}%
\def\CPYRTYEAR#1{\gdef\theCPYRTYEAR{#1}}%
\CPYRTYEAR{2003}%

\def\CLline{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}
%\\[-54pt]
 \TEN\it\hskip 23pt {Progress in Computational Fluid Design
Vol. \theVOL, No. \theISSUE, \thePUBYEAR}
 & \TEN\thepage\\[-2pt]
\HD{0}{8}\vspace*{7pt}
\end{tabular*}
}

\def\bi#1{\mbox{\bf\itshape{#1}}}

\def\BottomCatch{%
%\vskip 10pt
\thispagestyle{empty}%
\begin{table}[b]%
\NINE\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lcr@{}}%
\\[-20pt]
Copyright \copyright\ 200x Inderscience Enterprises Ltd. & &%
\end{tabular*}%
\vskip -8pt%
\end{table}%
}%

%\font\titfont=heavhs at 22pt
%\font\sectfont=hlb at 9pt

\def\EGT{\footnotesize\baselineskip 10pt plus 0.23pt minus 0.21pt\relax}%
\def\NINE{\fontsize{9}{10}\selectfont}
\def\TEN{\normalsize\baselineskip12pt plus .25pt minus .24pt\relax}%
\def\TWELVE{\fontsize{12}{14}\selectfont\baselineskip14pt plus .25pt minus .24pt\relax}%
%\def\FOURTEEN{\fontsize{22}{24}\selectfont\baselineskip16pt plus .25pt minus .24pt\relax}%
\def\FOURTEEN{\fontsize{22}{24}\baselineskip25pt plus 1pt minus 1pt\relax}

\newcommand\abstractname{Abstract:\enskip}
\newcommand\keywordname{Keywords:\enskip\null}

\newbox\t@abstract
\newbox\t@receivedname
\newbox\temptbox
\newdimen\tempdime


\def\authorA#1{\gdef\@authorA{\TWELVE\bf #1}}
\def\authorB#1{\gdef\@authorB{\TWELVE\bf #1}}
\def\authorC#1{\gdef\@authorC{\TWELVE\bf #1}}
\def\authorD#1{\gdef\@authorD{\TWELVE\bf #1}}
\def\and{\lowercase{and}\ }

\def\affA#1{\gdef\@affA{\TEN\rm #1}}
\def\affB#1{\gdef\@affB{\TEN\rm #1}}
\def\affC#1{\gdef\@affC{\TEN\rm #1}}
\def\affD#1{\gdef\@affD{\TEN\rm #1}}
\def\KEY#1{\gdef\@KEY{\NINE #1}}
\def\REF#1{\gdef\@REF{\NINE #1}}
\def\BIO#1{\gdef\@BIO{\NINE #1}}

\def\HD#1#2{\vrule height #1pt depth #2pt width 0pt\relax}

\long\gdef\ARTICLEOPENING#1{\long\gdef\@maketitle{#1}}%
\ARTICLEOPENING{}%

\def\abstract{\@ifnextchar[{\@abstract}{\@abstract[]}}

\def\@abstract[#1]{%
    \hsize\textwidth%
  \global\setbox\t@abstract=\vbox\bgroup%
    \parindent10pt\leftskip 14.1pc\rightskip 0pt
    \noindent {\bf\NINE
    \abstractname}\noindent\NINE
    \ignorespaces}

\def\endabstract{\egroup}

%FFF
\def\@maketitle{%
 % \begin{Flushleft}%
    \let \footnote \thanks%
      {%
    \ifx\@docHeading\undefined\else{\TWELVE\@docHeading\par}\vskip.2\p@\fi%
       \hyphenpenalty 10000%
        \HD{40}{0}\fontsize{22}{26}\selectfont\bfseries
        \raggedright
        \leftskip 14.1pc\rightskip 0pt
        \parbox[t]{25pc}{\raggedright\@title\baselineskip23pt plus 1pt minus 1pt\relax}\HD{10}{0}%
         \par}%
    \ifx\@subtitle\undefined%
    \else%
            \TWELVE\it
            \HD{20}{0}\@subtitle%
     \fi%
%%%authorA  %GGG
    \ifx\@authorA\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{42}{0}\@authorA\HD{0}{1}\par
    \fi%
    \ifx\@affA\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{1}{0}\@affA\HD{0}{8.5}\par
    \fi%
%%%authorB  %GGG
    \ifx\@authorB\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{0}{0}\@authorB\HD{0}{1}\par
    \fi%
    \ifx\@affB\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{1}{0}\@affB\HD{0}{8.5}\par
    \fi%
%%%authorC  %GGG
     \ifx\@authorC\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{0}{0}\@authorC\HD{0}{1}\par
    \fi%
    \ifx\@affC\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{1}{0}\@affC\HD{0}{8.5}\par
    \fi%
%%%authorD  %GGG
    \ifx\@authorD\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{0}{0}\@authorD\HD{0}{1}\par
    \fi%
    \ifx\@affD\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
    \HD{1}{0}\@affD\HD{0}{8.5}\par
    \fi%
%%%
    \HD{18}{0}
    {\unvbox\t@abstract}%
        \ifx\@KEY\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
     \HD{19}{0}\NINE{\bf Keywords:}\ \@KEY \HD{0}{0}\par
    \fi%
    \ifx\@REF\undefined%
    \else%
    \leftskip 14.1pc\rightskip 0pt
     \HD{19}{0}\NINE{\bf Reference}\ \@REF \HD{0}{0}\par%
    \fi%
    \ifx\@BIO\undefined%
    \else%
    \HD{19}{0}\NINE  {\bf Biographical notes:}\ \@BIO \HD{0}{0}\par%
    \fi%
%      \ifx\@ABB\undefined%
%    \else%
%       \HD{22}{0}\@ABB \HD{0}{0}\par%
%    \fi%
%      \ifx\@NOMEN\undefined%
%    \else%
%       \HD{22}{0}\@NOMEN \HD{15}{0}\par%
%    \fi%
%   \ifx\@KEY\undefined%
%    \else%
%        \  \par%
%    \fi%
%     \end{Flushleft}%
%  \end{Flushleft}%
  \vskip 31pt plus 2pt%
  \relax%
}%


\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand \thepart {\@Roman\c@part}
\renewcommand \thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\newcommand\part{\par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname~\thepart
       \par\nobreak
     \fi
     \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
\def\firstsection#1{\vskip 12pt
\stepcounter{section}
\noindent\begin{tabular}{@{}l}
\\[-44pt]
\rule{19.5pc}{.25pt}\\
{\fontsize{9}{11}\selectfont\bfseries \thesection\enskip #1}\\[-5.5pt]
\rule{19.5pc}{.25pt}\\[-42pt]
\end{tabular}
}

\def\mcr{\\ \\[-6pt]}
\def\x{\end{document}}
%\newcommand{\pafg}[2]{\frac{\partial #1}{\partial #2}}
\def\p#1{\phantom{#1}}
\def\vvp{\vadjust{\vfill\pagebreak}}
\def\vvsp#1{\vadjust{\vspace{#1pt}\pagebreak}}
\def\vvSP#1{\HD{0}{#1}\pagebreak}
\def\tra#1{\spaceskip 0.#1em plus 0.15em minus 0.09em\relax}%
\def\vs#1{\vadjust{\vspace{#1pt}}}
\def\mathtight{
\thickmuskip2mu \thinmuskip1mu \medmuskip=1.5mu\relax}


\def\bi{\bfseries\itshape}
%\newcommand\section{\@startsection {section}{1}{\z@}%
%                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
%                                   {2.3ex \@plus.2ex}%
%                                   {\normalfont\Large\bfseries}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\large\bfseries}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
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
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

\newbox\boxH
\newbox\boxHrunin

\let\enddescription\endlist

\newenvironment{UL}{\TEN\begin{list}{}{%
 \topsep=12pt plus 3pt minus 3pt
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 \labelsep=0pt %%5pt
 \labelwidth=0pt
 \leftmargin=10pt
 \itemindent=-10pt
 \rightmargin=0pt
 \listparindent=0pt}}{\end{list}}%

\newenvironment{NLh}[1]{%
 \setbox\boxH=\hbox{{\rm #1}}%
 \begin{list}{\TEN\rm{\arabic{enumii}.}}{%
 \usecounter{enumii}
 \topsep=9pt plus 3pt minus 3pt
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth=\wd\boxH
 \leftmargin=\labelwidth
 \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0em
 \itemindent 0pt}}{\end{list}}%

%\def\KluwerBullet{\raise0.95pt\hbox{\egt\bf\char"0F}}

\def\LM{\leavevmode }
\newenvironment{BL}{%
 \begin{list}{\LM\lower0.6pt\hbox{$\bullet$}}{%
 \topsep=6pt plus 3pt minus 2pt
 \partopsep=0pt
 \itemsep=6pt %2pt
 \parsep=0pt
 \labelsep=12pt
 \labelwidth=4.5pt
 \rightmargin=0pt
 \listparindent=0pt
 \itemindent=0pt %=\leftmargin
 \leftmargin=16pt}}{\end{list}}%

\def\sectiona#1{\vskip 12pt
\noindent\begin{tabular}{@{}l}
\\[-9.5pt]
\rule{19.5pc}{.25pt}\\
{\fontsize{9}{11}\selectfont\bfseries #1}\\[-5.5pt]
\rule{19.5pc}{.25pt}\\[5pt]
\end{tabular}}

\def\NONUMBIB{%
\def\thebibliography##1{%
    \par{\sectiona{\REFERENCES\HD{6}{6}}}%
    \list{}{\settowidth\labelwidth{##1}
    \setlength{\topsep}{0pt}
    \labelsep 0pt\labelwidth 0pt
    \leftmargin=24pt
     \itemindent-24pt
    \itemsep 0pt
    \def\bibitem{\item\NINE
     \baselineskip11pt plus 2pt minus 2pt\relax}}
    \def\newblock{\hskip .11em plus .33em minus .07em}
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=1000\relax}}

\NONUMBIB % <-- default

%\def\NUMBIB{\def\@biblabel##1{\bf ##1}%
%\def\thebibliography##1{%
%  \par{\sectiona{REFERENCES\HD{0}{0}}\vs{-8}}%
%  \list
%  {\@biblabel{\arabic{enumiv}}}{\NINE
%    \baselineskip 10pt plus 0.3pt minus 0.25pt\relax
%    \settowidth\labelwidth{\@biblabel{##1}}%
%    \leftmargin\labelwidth
%    \advance\leftmargin\labelsep
%    \usecounter{enumiv}%
%    \let\p@enumiv\@empty
%    \labelsep 10pt\labelwidth 13.5pt
%    \leftmargin = 21pt
%    \itemsep  2pt
%     \itemindent0pt
%    \parskip4pt\topsep4pt
%    \def\theenumiv{\arabic{enumiv}}}%
%    \def\newblock{\hskip .11em plus.33em minus.07em}%
%    \sloppy\clubpenalty4000\widowpenalty4000
%    \sfcode`\.=1000\relax
%    \parsep0pt}}

\def\SECRR{\rightskip2pc plus 18pc\relax}

\def\firstsection#1{\vskip 12pt
\stepcounter{section}
\noindent\begin{tabular}{@{}l}
\\[-44pt]
\rule{19.5pc}{.25pt}\\
{\fontsize{9}{11}\selectfont\bfseries \thesection\enskip #1}\\[-5.5pt]
\rule{19.5pc}{.25pt}\\[-42pt]
\end{tabular}
}

\def\section#1{\vskip 12pt
\stepcounter{section}
\noindent\begin{tabular}{@{}l}
\\[-9.5pt]
\rule{19.5pc}{.25pt}\\
{\fontsize{9}{11}\selectfont\bfseries \thesection\enskip #1}\\[-5.5pt]
\rule{19.5pc}{.25pt}\\[5pt]
\end{tabular}}

\if@titlepage
%  \newenvironment{abstract}{%
%      \titlepage
%      \null\vfil
%      \@beginparpenalty\@lowpenalty
%      \begin{center}%
%        \bfseries \abstractname
%        \@endparpenalty\@M
%      \end{center}}%
%     {\par\vfil\null\endtitlepage}
%\else
%  \newenvironment{abstract}{%
%      \if@twocolumn
%        \section*{\abstractname}%
%      \else
%        \small
%        \begin{center}%
%          {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}%
%        \end{center}%
%        \quotation
%      \fi}
%      {\if@twocolumn\else\endquotation\fi}
\fi
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\if@compatibility
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\z@
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
    }
\else
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
     \if@twoside\else
        \setcounter{page}\@ne
     \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \renewcommand\thesection{\@Alph\c@section}}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand \theequation {\@arabic\c@equation}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\centerline{#1: #2}}%
  \ifdim \wd\@tempboxa >\hsize
    \centerline{#1: #2}\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
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
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename
      \@mkboth{%
          \MakeUppercase\listtablename}{\MakeUppercase\listtablename}}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\renewenvironment{thebibliography}[1]
     {\sectiona{\refname
        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\providecommand\bibsection{}
  {\renewcommand\bibsection{\sectiona{\refname
    \@mkboth{\MakeUppercase{\refname}}{\MakeUppercase{\refname}}}}}
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{REFERENCES}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
%\newcommand\abstractname{Abstract}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
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

\newdimen\markwidth
\newdimen\blth
\newdimen\marklength
\newdimen\markdistance
\newdimen\trimwidth
\newdimen\trimheight
\newdimen\auxaaa
\newdimen\auxbbb
\newdimen\auxccc
\newdimen\auxddd

\long\def\koo#1#2#3{\vbox to0pt{\hsize0pt\kern #2\hbox to0pt{\kern
#1{#3}\hss}\vss}}
\long\def\zeroCC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}}}
\long\def\zeroUC#1{\vbox to0pt{\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}}}
\long\def\zeroLL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}}}
\long\def\zeroUR#1{\vbox to0pt{\hbox to0pt{\hss #1}\vss}}
\long\def\zeroUL#1{\vbox to0pt{\hbox to0pt{#1\hss}\vss}}
\long\def\zeroCL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}\vss}}
\long\def\zeroCR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}\vss}}

\markwidth=.5pt %0.4 truept % Thickness of the mark line
\blth=0.2pt
\marklength=20 truept % Length of the mark line
\markdistance=3 truept % Distance of the real corner to the beg. of the mark
%%% \trimwidth=\textwidth+\oddsidemargin+\evensidemargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\auxaaa\marklength\advance\auxaaa\markdistance
\def\oddevenswitch{\ifodd\c@page \auxbbb=\oddsidemargin
 \else\auxbbb=\evensidemargin\fi}%

%%% First reset all page position parameters depending on the job:
% \PagePositionCropYes - should be just for testing
% \PagePositionCropNo  - should be the default behavior
% \PagePositionBoxYes  - gives the page box for testing

\def\PagePositionCropYes#1#2#3#4#5#6#7#8{%
  \oddsidemargin=#1%  % [3] inner margin / gutter
  \evensidemargin=#2% % [4] outside margin
  \topmargin=#3%      % [5] \rhsink = \topmargin + \headheight (=12pt here)
  \headheight=#4%     % [6] VERY IMPORTANT! (see above)
  \trimwidth=#5%      % [7]
  \trimheight=#6%     % [8]
  \hoffset=#7%        % = -1in to move the area in the upper left corner
  \voffset=#8%        % = -1in to move the area in the upper left corner
%%%
 \def\CROPMARKS{\oddevenswitch
 \hbox to0pt{\kern-\auxbbb
   \vbox to0bp{\kern-\topmargin\kern-\headheight
   \hbox to\trimwidth{%
     \koo{-\auxaaa}{-\markwidth}{\VRHDW{\markwidth}{0pt}{\marklength}}%
     \koo{-\markwidth}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}%
     \hfill
     \koo{\markdistance}{-\markwidth}{\VRHDW{\markwidth}{0pt}{\marklength}}%
     \koo{0pt}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}}%
   \nointerlineskip\vskip\trimheight
   \nointerlineskip
   \hbox to\trimwidth{%
     \koo{-\auxaaa}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
     \koo{-\markwidth}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}%
     \hfill
     \koo{\markdistance}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
     \koo{0pt}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}}%
   \vss}\hss}}%
}%

\def\PagePositionCropNo#1#2#3#4#5#6#7#8{%
  \oddsidemargin=#1%  % [3] inner margin / gutter
  \evensidemargin=#2% % [4] outside margin
  \topmargin=#3%      % [5] \rhsink = \topmargin + \headheight (=12pt here)
  \headheight=#4%     % [6] VERY IMPORTANT! (see above)
  \trimwidth=#5%      % [7]
  \trimheight=#6%     % [8]
  \hoffset=#7%        % = -1in to move the area in the upper left corner
  \voffset=#8%        % = -1in to move the area in the upper left corner
%%%
\def\CROPMARKS{}%
}%

\def\PagePositionBoxYes#1#2#3#4#5#6#7#8{%
  \oddsidemargin=#1%  % [3] inner margin / gutter
  \evensidemargin=#2% % [4] outside margin
  \topmargin=#3%      % [5] \rhsink = \topmargin + \headheight (=12pt here)
  \headheight=#4%     % [6] VERY IMPORTANT! (see above)
  \trimwidth=#5%      % [7]
  \trimheight=#6%     % [8]
  \hoffset=#7%        % = -1in to move the area in the upper left corner
  \voffset=#8%        % = -1in to move the area in the upper left corner
%%%
 \def\CROPMARKS{\oddevenswitch
 \hbox to0pt{\kern-\auxbbb
   \vbox to0bp{\kern-\topmargin\kern-\headheight
   \hbox to\trimwidth{%
     \auxccc=\trimwidth \advance\auxccc by 2\blth
     \auxddd=\trimheight \advance\auxddd by 2\blth
     \koo{-\blth}{0pt}{\zeroLL{\VRHDW{\blth}{0pt}{\auxccc}}}%
     \koo{-\blth}{-\blth}{\zeroUL{\VRHDW{\auxddd}{0pt}{\blth}}}%
     \koo{\trimwidth}{-\blth}{\zeroUL{\VRHDW{\auxddd}{0pt}{\blth}}}%
     \koo{-\blth}{\trimheight}{\zeroUL{\VRHDW{\blth}{0pt}{\auxccc}}}\hss}%
   \vss}\hss}}%
}%

\textheight=696pt
\textwidth=516pt

\endinput
%%
%% End of file `doublecol.cls'.

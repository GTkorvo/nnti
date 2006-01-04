#

TEX = latex

PDFTEX = pdflatex

BTEX = bibtex

BIBFILE = paper.bib

PAPER = paper

PARTS = paper.tex

PICTURES =

OTHER = acmtrans.bst acmtrans2m.cls

all: ${PAPER}.pdf

########################################################################
#
# Dependencies - 3 are required for latex/bibtex
# - a rule must exist for the .dvi, the .bbl and the .aux file
#   in order for a clean make with a new .tex file
#
########################################################################

${PAPER}.aux: 
	pdflatex ${PAPER}.tex

${PAPER}.bbl: ${PAPER}.aux
	${BTEX} ${PAPER}

all:    ${PAPER}.tex ${PAPER}.bbl
	pdflatex ${PAPER}.tex 

tarfile:
	tar chvf ${PAPER}.tar Makefile ${PARTS} ${PICTURES} ${BIBFILE} ${OTHER}
	compress ${PAPER}.tar

clean:
	/bin/rm -f ${PAPER}.pdf ${PAPER}.ps ${PAPER}.dj ${PAPER}.dvi *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z



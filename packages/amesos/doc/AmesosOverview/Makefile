#

TEX = latex

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

${PAPER}.dvi: ${PAPER}.tex ${PARTS} ${PAPER}.bbl
	${TEX} ${PAPER}
	${TEX} ${PAPER}

${PAPER}.bbl: ${PAPER}.aux
	${BTEX} ${PAPER}

${PAPER}.aux: ${PAPER}.tex ${PARTS} force
	${TEX} ${PAPER}

force:

${PAPER}.ps:    ${PAPER}.dvi
	dvips -Ppdf -o ${PAPER}.ps ${PAPER}

${PAPER}.pdf:    ${PAPER}.ps
	ps2pdf ${PAPER}.ps ${PAPER}.pdf

tarfile:
	tar chvf ${PAPER}.tar Makefile ${PARTS} ${PICTURES} ${BIBFILE} ${OTHER}
	compress ${PAPER}.tar

clean:
	/bin/rm -f ${PAPER}.pdf ${PAPER}.ps ${PAPER}.dj ${PAPER}.dvi *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z



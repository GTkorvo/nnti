##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

PSTEX  := latex
TEX    := pdflatex
BIB    := bibtex
DVI    := dvips
PDF    := ps2pdf13
SED    := sed
PSBOOK := psbook

XFIG := $(wildcard *fig)
PDFFIG := $(patsubst %.fig, %.pdf, $(XFIG))

BASEFILE := ModelEvaluatorOverview
OTHERFILES := references.bib

##---------------------------------------------------------------------------##

.PHONY: clean pdf ps realclean all denovo

all: $(BASEFILE).pdf

$(BASEFILE).pdf: $(BASEFILE).tex $(PDFFIG) $(OTHERFILES)
	$(TEX) $(BASEFILE)
	$(BIB) $(BASEFILE)
	$(TEX) $(BASEFILE)
	$(TEX) $(BASEFILE)

#set.pdf: set.fig 
#	fig2dev -L pdftex   -p set.pdf -m 0.4 set.fig set.pdf
#	fig2dev -L pdftex_t -p set.pdf -m 0.4 set.fig set.pdftex_t

clean:
	-rm -f *~
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.dvi
	-rm -f *.rel
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f *.out
	-rm -f *.toc *.lof *.lot
	-rm -f *.bak
	-rm -f TAGS

realclean: clean
	-rm -f $(BASEFILE).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##

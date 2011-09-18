
BASE_NAME = ThyraOverview
LATEX_COMMAND = pdflatex
BIBTEX_COMMAND = bibtex

SOURCEFILES= \
  $(BASE_NAME).tex \
  references.bib \
  ANA_LAL_APP_Figure.pdf \
  ThyraANAOpVecInterfaces.pdf

#
# Shortcut targets
#

default: pdf

pdf: $(BASE_NAME).pdf

clean:
	rm -f *.ps $(BASE_NAME).pdf *.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ .log *~

#
# File targets
#

$(BASE_NAME).pdf: $(SOURCEFILES)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(BIBTEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)

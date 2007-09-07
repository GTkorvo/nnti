
BASE_NAME = ASCVerticalIntegrationSummary
LATEX_COMMAND = latex
BIBTEX_COMMAND = bibtex

SOURCEFILES= \
            $(BASE_NAME).tex \
            references.bib

default: pdf

$(BASE_NAME).dvi: $(SOURCEFILES)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(BIBTEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)

ps: $(BASE_NAME).dvi
	dvips -t letter -o $(BASE_NAME).ps $(BASE_NAME)

pdf: $(BASE_NAME).dvi
	dvipdf $(BASE_NAME).dvi $(BASE_NAME).pdf

clean:
	rm -f *.ps *.pdf *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out \
           *.bbl *.blg *.lof *.lot *.tex~ .log *~

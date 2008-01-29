
BASE_NAME = ThyraCodingGuideLines
LATEX_COMMAND = latex
MAKEINDEX_COMMAND = makeindex
MAKEGLOSSARY_COMMAND = makeglossary
BIBTEX_COMMAND = bibtex

SOURCEFILES= \
  $(BASE_NAME).tex \
  guideline_defs.tex \
  CodingGuideLinesBody.tex \
  apdx_Summary.tex \
  apdx_CppCodingStandards.tex \
  apdx_MiscAmendments.tex \
  apdx_ArgumentsForConsistentStyle.tex \
  apdx_CodeReformattingGuidelines.tex \
  PassingValueObjects.tex \
  PassingReferenceObjects.tex \
  references.bib

#
# Shortcut targets
#

default: pdf

dvi: $(BASE_NAME).dvi

ps: $(BASE_NAME).ps

pdf: $(BASE_NAME).pdf

clean:
	rm -f *.ps $(BASE_NAME).pdf *.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ .log *~

#
# File targets
#

$(BASE_NAME).dvi: $(SOURCEFILES)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(MAKEINDEX_COMMAND) $(BASE_NAME)
	$(BIBTEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)

#	$(MAKEINDEX_COMMAND) -s $(BASE_NAME).ist -t $(BASE_NAME).glg -o $(BASE_NAME).gls $(BASE_NAME).glo
#	$(MAKEINDEX_COMMAND)  -t $(BASE_NAME).glg -o $(BASE_NAME).gls $(BASE_NAME).glo


$(BASE_NAME).ps: $(BASE_NAME).dvi
	dvips -t letter -o $(BASE_NAME).ps $(BASE_NAME)

$(BASE_NAME).pdf: $(BASE_NAME).ps
	ps2pdf $(BASE_NAME).ps $(BASE_NAME).pdf


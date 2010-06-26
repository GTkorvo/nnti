
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


#
# Make the sample namespace program
#
# NOTE: Define CC on the env to run these targets
#

CCFLAGS=-ansi -pedantic -Wall


ns-std:
	$(CC) $(CCFLAGS) -o NamespaceClassUsingIssues.exe NamespaceClassUsingIssues.cpp

ns-dup-class:
	$(CC) $(CCFLAGS) -DSHOW_DUPLICATE_CLASS_A -o NamespaceClassUsingIssues.exe NamespaceClassUsingIssues.cpp

ns-missing-using-decl:
	$(CC) $(CCFLAGS) -DSHOW_MISSING_USING_DECL -o NamespaceClassUsingIssues.exe NamespaceClassUsingIssues.cpp

ns-bad-using-direc:
	$(CC) $(CCFLAGS) -DSHOW_ERRONEOUS_USING_DIRECTIVE -o NamespaceClassUsingIssues.exe NamespaceClassUsingIssues.cpp



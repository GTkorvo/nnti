
LATEX=lualatex
BIBTEX=bibtex
TARGET0=PyTrilinos
TARGET1=Response

default: pytrilinos response

pytrilinos:
	$(LATEX) $(TARGET0)

response:
	$(LATEX) $(TARGET1)

bib:
	$(BIBTEX) $(TARGET0)

clean:
	@rm -f $(TARGET0)-blx.bib \
	       $(TARGET0).aux \
	       $(TARGET0).toc \
	       $(TARGET0).out \
	       $(TARGET0).auxlock \
	       $(TARGET0).bbl \
	       $(TARGET0).blg \
	       $(TARGET0).log \
	       $(TARGET0).pdf \
	       $(TARGET0).xwm \
	       $(TARGET0).snm \
	       $(TARGET0).nav
	@rm -f $(TARGET1)-blx.bib \
	       $(TARGET1).aux \
	       $(TARGET1).toc \
	       $(TARGET1).out \
	       $(TARGET1).auxlock \
	       $(TARGET1).bbl \
	       $(TARGET1).blg \
	       $(TARGET1).log \
	       $(TARGET1).pdf \
	       $(TARGET1).xwm \
	       $(TARGET1).snm \
	       $(TARGET1).nav
	@rm -f missfont.log \
	       x.log

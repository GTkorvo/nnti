name = thyra_overhead

default: $(name).pdf

$(name).dvi: $(name).tex
	latex $(name).tex

$(name).pdf: $(name).dvi
	dvipdf $(name).dvi $(name).pdf

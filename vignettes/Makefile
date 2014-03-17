all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('eiR.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
clean:
	rm -f *.md

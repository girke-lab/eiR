all:
	#"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('eiR.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('eiR.Rmd')"
	#pandoc --output eiR.html -H header.html --from markdown_github-hard_line_breaks+superscript+tex_math_dollars+pandoc_title_block --smart --table-of-contents eiR.md;
	pandoc --output eiR.html -H header.html --from markdown --smart --table-of-contents eiR.md;
clean:
	rm -f *.md

PKG_NAME=eiR

all:
	#"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('eiR.Rmd')"
	#pandoc --self-contained --output eiR.html -H header.html --from markdown --smart --table-of-contents eiR.md;

#	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('${PKG_NAME}.Rmd')"
#	pandoc --self-contained --output ${PKG_NAME}.html -H header.html --from markdown --smart --table-of-contents --template html.template ${PKG_NAME}.md;

	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitrBootstrap); knit_bootstrap('${PKG_NAME}.Rmd',show_code=TRUE)"
clean:
	rm -f *.md

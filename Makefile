PACKAGE=spict
VERSION=1.2.6
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip
R=R
SUBDIRS := $(wildcard testmore/*/.)

ifeq (testonemore,$(firstword $(MAKECMDGOALS)))
  ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARG):;@:)
endif

all:
	make doc-update
	make install
	make pdf

doc-update:
	echo "roxygen2::roxygenize('spict/', roclets=c('rd', 'collate', 'namespace'))" | R --slave
#echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | R --slave	

build-package:
	echo 'source("make.description.R")' | R --vanilla
	R CMD build --no-build-vignettes --resave-data=no $(PACKAGE)

install:
	make build-package
	R CMD INSTALL --preclean --no-multiarch $(TARBALL)
	date

unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

check:
	R CMD check $(PACKAGE)


quick-install: $(PACKAGE)/src/spict.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/spict.so: $(PACKAGE)/src/spict.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('spict.cpp','-O0 -g')" | $(R) --slave

vignette:  $(PACKAGE)/vignettes/vignette.Rmd
	R -e "rmarkdown::render('spict/vignettes/vignette.Rmd', rmarkdown::pdf_document())"

.PHONY: testmoreseq testonemore testmore $(SUBDIRS)

testmore:
	$(MAKE) -j $(NPROCS) testmoreseq

testmoreseq: $(SUBDIRS)

testonemore:
	@$(MAKE) testmore/$(ARG)/.

$(SUBDIRS):
	@cp testmore/Makefile $@
	@$(MAKE) -i -s -C $@
	@rm -f $@/Makefile


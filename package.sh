#!/bin/bash
cd spict
#echo "library(devtools); library(roxygen2); document()" | R --slave
echo "roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace'))" | R --slave
cd ..
#rm spict.pdf
#R CMD Rd2pdf spict
#cp spict.pdf spict/
R CMD build spict
#R CMD check spict
R CMD INSTALL spict

#!/bin/bash
cd spict
echo "roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace'))" | R --slave
cd ..
#rm spict.pdf
#R CMD Rd2pdf spict
#cp spict.pdf spict/
R CMD build spict
#R CMD check spict_0.1.tar.gz
R CMD INSTALL spict_0.1.tar.gz
#R CMD INSTALL spicttest_0.1.tar.gz
echo "Done!"
date
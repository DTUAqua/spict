#!/bin/bash
SECONDS=0
tmpdir=$(mktemp -d)
RED='\033[31m'

ref="master"
vignfn="spict/vignettes/spict_handbook.Rmd"

while getopts ":r:v:" opt; do
  case $opt in
    r)
      ref=$OPTARG
  ;;
    v)
      vignfn=$OPTARG
  ;;
  :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Working on temporary directory: " $tmpdir

printmsg () {
  CYAN='\033[96m'
  NC='\033[0m' # No Color
  printf "${CYAN}*************************************************\n"
  printf "$1\n"
  printf "*************************************************${NC}\n"
}


printmsg "Installing spict from github, reference:  $ref"
echo "withr::with_libpaths(new = '$tmpdir', {remotes::install_cran('ellipse', quiet = FALSE, repo = 'http://cran.rstudio.com')}, action=\"prefix\")" | R --slave
echo "withr::with_libpaths(new = '$tmpdir', {remotes::install_github('DTUAqua/spict/spict', ref = '$ref', quiet = FALSE)}, action=\"prefix\")" | R --slave


printmsg "Run vignette examples"

./runexamples.R -r $vignfn -o $tmpdir -l $tmpdir -f old.md


printmsg "Install local spict version"

echo "withr::with_libpaths(new = '$tmpdir', {remotes::install_local('spict', quiet = TRUE)}, action=\"prefix\")" | R --slave

printmsg "Run vignette examples"

./runexamples.R -r $vignfn -o $tmpdir -l $tmpdir -f new.md



printmsg "All examples are done using the two verions of R\nTime ellapsed" $SECONDS "seconds"


printmsg "Compare examples"


if [ ! -e $tmpdir/new.md ] || [ ! -e $tmpdir/old.md ]; then
  printf "${RED}The example files are not found\n"
  exit 124
fi

newhash=$(md5sum $tmpdir/new.md | cut -d " " -f1)
oldhash=$(md5sum $tmpdir/old.md | cut -d " " -f1)

if [ $newhash == $oldhash ] ; then 
  echo "The two versions of spict are producing the same results."
else 
  printf "${RED}There are differences in the two versions of spict."
  diff $tmpdir/old.md $tmpdir/new.md
  exit 1
fi

exit

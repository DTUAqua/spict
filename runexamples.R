#!/usr/bin/env Rscript
require(docopt, quietly = TRUE)
'Usage:
  runexamples [-l <LIBRARY_PATH> -r <RMD_FILE> -o <OUT_DIR> -f <OUT_FILE>]

Options:
  -l Library path
  -r Rmdfile [default: spict/vignettes/vignette.Rmd]
  -o Output directory
  -f outfilename

]' -> doc

get.slim.md <- function(rmdfile, outdir, fn) {
  tmpdir <- tempdir()
  onlyR <- knitr::purl(rmdfile, output = tempfile(), documentation = 1, quiet = TRUE)
  rmd <- knitr::spin(onlyR, knit = FALSE) 
  out <- rmarkdown::render(rmd, output_format = "md_document", output_options = list(variant = "markdown"), quiet = TRUE, output_dir = outdir, output_file = fn)
  withimg <- readLines(out)
  cat(withimg[ ! grepl("^<img", withimg)], file = out, append = FALSE, sep = "\n")
  out 
}

opts <- docopt(doc)
fn <- normalizePath(opts$r)
withr::with_libpaths(new = opts$l, {
  library(spict)
  get.slim.md(opts$r, opts$o, opts$f)
},action="prefix")



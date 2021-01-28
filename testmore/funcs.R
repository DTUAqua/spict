## General functions for testmore
if (! "cnsl" %in% ls()) cnsl <- FALSE

## Test functions
out <- function(..., sep = "", append = TRUE){
    if (cnsl) {
        cat(..., "\n", sep = sep)
    }else{
        cat(..., "\n", file = "res.out", append = append, sep = sep)
    }
}

get_nchar <- function(...){
    nchar(paste(as.character(unlist(list(...))), collapse = " "))
}

header <- function(..., append = TRUE){
    out("\n", ..., "\n", rep("=", get_nchar(...)), "\n", append = append)
}

test_this <- function(title="", expression) {
  out(title)
  tryCatch(out(capture.output(eval(expression)), sep = "\n"),
           error = function(e) out("Error:", e$message))
}

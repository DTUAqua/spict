## General functions for testmore

## Test functions
out <- function(..., sep = "", append = TRUE, cnsl = FALSE){
    if(cnsl){
        cat(..., "\n", sep = sep)
    }else{
        cat(..., "\n", file = "res.out", append = append, sep = sep)
    }
}

get_nchar <- function(...){
    nchar(paste(as.character(unlist(list(...))), collapse = " "))
}

header <- function(..., append = TRUE, cnsl = FALSE){
    out("\n", ..., "\n", rep("=", get_nchar(...)), "\n", append = append, cnsl=cnsl)
}

test_this <- function(title="", expression, cnsl = FALSE) {
  out(title, cnsl=cnsl)
  tryCatch(out(capture.output(eval(expression)), sep = "\n", cnsl=cnsl),
           error = function(e) out("Error:", e$message, cnsl=cnsl))
}

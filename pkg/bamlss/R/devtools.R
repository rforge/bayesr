## source all functions in the package once more
sbayesr <- function(R2BayesX = FALSE, dir = NULL) {
  if(is.null(dir)) {
    dir <- if(!R2BayesX) {
      "~/svn/bayesr/pkg/bamlss/R"
    } else "~/svn/bayesr/pkg/bamlss/R"
  }
  dir <- path.expand(dir)
  f <- file.path(dir, list.files(dir))
  sapply(f, source)
  invisible(NULL)
}


## open the test script
tscript <- function(file = NULL) {
  if(is.null(file))
    file <- "~/svn/bayesr/pkg/bamlss/inst/tscript.R"
  file <- path.expand(file)
  system(paste(shQuote("gedit"), shQuote(file)), wait = FALSE)
  invisible(NULL)
}

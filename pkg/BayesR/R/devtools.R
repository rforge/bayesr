## source all functions in the package once more
sbayesr <- function(R2BayesX = FALSE, dir = NULL) {
  if(is.null(dir)) {
    dir <- if(!R2BayesX) {
      "~/svn/bayesr/pkg/BayesR/R"
    } else "~/svn/bayesr/pkg/R2BayesX/R"
  }
  dir <- path.expand(dir)
  f <- file.path(dir, list.files(dir))
  sapply(f, source)
  invisible(NULL)
}


## open the test script
tscript <- function(file = NULL) {
  if(is.null(file))
    file <- "~/svn/bayesr/pkg/BayesR/inst/tscript.R"
  file <- path.expand(file)
  system(paste(shQuote("gedit"), shQuote(file)), wait = FALSE)
  invisible(NULL)
}

## source all functions in the package once more
sbayesr <- function(dir = NULL) {
  if(is.null(dir))
    dir <- "~/svn/bayesr/pkg/BayesR/R"
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

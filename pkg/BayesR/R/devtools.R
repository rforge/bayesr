## source all functions in the package once more
sbayesr <- function(dir = NULL) {
  if(is.null(dir))
    dir <- "~/sci/articles/BayesR/BayesR/R"
  dir <- path.expand(dir)
  f <- file.path(dir, list.files(dir))
  sapply(f, source)
  invisible(NULL)
}


## open the test script
tscript <- function(file = NULL) {
  if(is.null(file))
    file <- "~/sci/articles/BayesR/BayesR/inst/tscript.R"
  file <- path.expand(file)
  system(paste(shQuote("gedit"), shQuote(file)), wait = FALSE)
  invisible(NULL)
}

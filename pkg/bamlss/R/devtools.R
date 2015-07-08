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


## dynamic C code.
compile <- function(dir = NULL, tdir = NULL) {
   hold <- getwd()
   on.exit(setwd(getwd()))
   if(is.null(dir)) dir <- "~/svn/bayesr/pkg/bamlss/src"
   dir <- path.expand(dir)
   if(is.null(tdir)) {
      dir.create(tdir <- tempfile())
      on.exit(unlink(tdir), add = TRUE)
   }
   tdir <- path.expand(tdir)
   setwd(tdir)
   if(file.exists(file.path(dir, "Makevars")))
      file.copy(file.path(dir, "Makevars"), file.path(tdir, "Makevars"))
   cf <- grep(".c", dir(dir), value = TRUE, fixed = TRUE)
   for(j in cf) {
      file.copy(file.path(dir, j), file.path(tdir, j))
      system(paste("R CMD SHLIB", j))
      dyn.load(gsub(".c", ".so", j, fixed = TRUE))
   }
   setwd(hold)
}


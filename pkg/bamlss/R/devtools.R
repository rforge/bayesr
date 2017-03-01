## Source all functions in the package once more.
if(FALSE) {
sbayesr <- function(R2BayesX = FALSE, dir = NULL, svn = "svn") {
  if(is.null(dir)) {
    dir <- if(!R2BayesX) {
      paste("~/", svn, "/bayesr/pkg/bamlss/R", sep = "")
    } else paste("~/", svn, "/bayesr/pkg/bamlss/R", sep = "")
  }
  dir <- path.expand(dir)
  f <- file.path(dir, list.files(dir))
  sapply(f, source)
  invisible(NULL)
}


## Open the test script.
tscript <- function(file = NULL) {
  if(is.null(file))
    file <- "~/SVN/bayesr/pkg/bamlss/inst/tscript.R"
  file <- path.expand(file)
  system(paste(shQuote("gedit"), shQuote(file)), wait = FALSE)
  invisible(NULL)
}


## Dynamic C code.
compile <- function(dir = NULL, tdir = NULL, svn = "svn")
{
  hold <- getwd()
  on.exit(setwd(getwd()))
  if(is.null(dir)) dir <- paste("~/", svn, "/bayesr/pkg/bamlss/src", sep = "")
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
  cf <- cf[!grepl(".c~", cf, fixed = TRUE)]
  cf <- cf[!grepl("_init.c", cf, fixed = TRUE)]
  for(j in cf) {
    file.copy(file.path(dir, j), file.path(tdir, j))
    system(paste("R CMD SHLIB", j))
    dyn.load(gsub(".c", ".so", j, fixed = TRUE))
  }
  setwd(hold)
}
}


## Object size.
list.size <- function(x, j = NULL)
{
  if(is.list(x) & length(x)) {
    if(is.null(names(x)))
      names(x) <- paste("x", 1:length(x), sep = "")
    for(i in names(x)) {
      if(!is.null(j))
        cat(j, i, ".. ")
      else
        cat(i, ".. ")
      print(object.size(x[[i]]), units = "Mb")
      if(is.list(x[[i]]))
        list.size(x[[i]], j = c(j, ".. "))
    }
  }
  invisible(NULL)
}


## Derivatives helper.
compute_derivatives <- function(loglik, parameters, density = NULL, expectation = FALSE)
{
  if(!is.character(loglik))
    stop("argument loglik must be a character!")
  if(!is.character(parameters))
    stop("argument parameters must be a character!")
  if(!grepl("y", loglik))
    stop("the response y is missing in loglik!")

  if(expectation) {
    for(p in parameters) {
      v <- paste(p, ' <- ', 'rSymPy::Var("', p, '")', sep = '')
      eval(parse(text = v))
    }
    y <- rSymPy::Var("y")
    if(is.null(density))
      density <- paste("exp(", loglik, ")", sep = "")
  }

  score <- hess <- list()
  for(p in parameters) {
    dldp <- D(parse(text = loglik), p)
    score[[p]] <- dldp
    d2ld2p <- D(dldp, p)
    if(expectation) {
      int <- paste("integrate(", gsub("^", "**", gsub(" ", "",
        paste("(", paste(deparse(d2ld2p), collapse = ""), ")*(", density, ")", sep = "")), fixed = TRUE), ",y)", sep = "")
      Ep <- rSymPy::sympy(int)
      Ep <- gsub("**", "^", Ep, fixed = TRUE)
      Ep <- parse(text = Ep)
      hess[[p]] <- Ep
    } else {
      hess[[p]] <- d2ld2p
    }
  }

  return(list("score" = score, "hess" = hess))
}

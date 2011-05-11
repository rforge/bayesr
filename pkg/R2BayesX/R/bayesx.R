bayesx <-
function(formula, data, weights = NULL, subset = NULL, 
  offset = NULL, na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(...), ...)
{
  ## first check installation of bayesx
  ok <- check.install.bayesx(control$bin, verbose = FALSE)
  if(is.null(ok))
    stop("please install BayesX or verify bin path, also see function install.bayesx()")

  res <- list()
  res$formula <- formula
  weights <- deparse(substitute(weights), backtick = TRUE, width.cutoff = 500L)
  offset <- deparse(substitute(offset), backtick = TRUE, width.cutoff = 500L)
  subset <- deparse(substitute(subset), backtick = TRUE, width.cutoff = 500L)

  ## setup files for bayesx
  res$bayesx.setup <- parse.bayesx.input(formula, data,
    weights, subset, offset, na.action, contrasts, control)

  ## write prg file
  res$bayesx.prg <- write.bayesx.input(res$bayesx.setup)

  ## now estimate with BayesX
  res$bayesx.run <- run.bayesx(dir = res$bayesx.prg$file.dir, 
    prg.name = res$bayesx.prg$prg.name, verbose = res$bayesx.setup$verbose, 
    bin = control$bin)

  if(is.null(res$bayesx.setup$hlevel))
    tm <- res$bayesx.prg$model.name
  else
    tm <- s4hm(res$bayesx.prg$file.dir, control$model.name)

  ## read bayesx output files
  if(!grepl("ERROR:", res$bayesx.run$log[length(res$bayesx.run$log)]) || 
    res$bayesx.run$log[1L] == 0) {
    if(length(tm) > 1L) {
      res.h <- read.bayesx.output(res$bayesx.prg$file.dir, tm)
      for(k in 1:length(res.h)) {
        res.h[[k]]$model.fit$method = "MCMC"
        res.h[[k]]$model.fit$model.name <- tm[k] 
      }
      res.h <- res.h[length(res.h):1L]
      res.h[[1L]]$call <- match.call()
      res.h[[1L]][names(res)] <- res[names(res)]
      class(res.h) <- "bayesx"
    } else {
      res.h <- NULL
      res <- c(res, read.bayesx.model.output(res$bayesx.prg$file.dir, tm))
      res$call <- match.call()
      res$terms <- res$bayesx.setup$term.labels
      res$model.fit$family <- control$family
    }
  } else warning("an error occured during runtime of BayesX!")

  ## maybe remove output folder
  if((is.null(res$bayesx.setup$outfile)) && res$bayesx.setup$dir.rm) {
    wd <- getwd()
    try(Sys.chmod(res$bayesx.prg$file.dir, mode = "7777"), silent = TRUE)
    setwd(res$bayesx.prg$file.dir)
    files <- list.files()
    for(k in 1L:length(files)) {
      ok <- try(suppressWarnings(file.remove(files[k])), silent = TRUE)
      if(inherits(ok, "try-error"))
        try(suppressWarnings(system(paste("rm --force", files[k]))), silent = TRUE)
    }
    ok <- try(suppressWarnings(file.remove(res$bayesx.prg$file.dir)), silent = TRUE)
    if(inherits(ok, "try-error"))
      try(suppressWarnings(system(paste("rm --force", res$bayesx.prg$file.dir))), silent = TRUE)
    setwd(wd)
  }

  if(!is.null(res.h))
    return(res.h)
  else {
    eval(parse(text=paste("res <- list(", res$bayesx.prg$model.name, "= res)")))
    class(res) <- "bayesx"
    return(res)
  }
}


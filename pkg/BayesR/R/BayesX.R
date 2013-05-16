#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, chains = 1, combine = TRUE, n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, dir = NULL, ...)
{
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    on.exit(unlink(dir))
  } else dir <- path.expand(dir)

  setup <- function(x) {
    setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
      seed = seed, dir = dir, ...)
  }
  sampler <- function(x) { samplerBayesX(x, ...) }
  results <- function(x) { resultsBayesX(x, ...) }

  bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = randomize,
    setup = setup, sampler = sampler, results = results,
    cores = cores, combine = combine, ...)
}


####################################
## (2) BayesX specific functions. ##
####################################
tranformBayesX <- function(x, ...)
{
  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    for(j in seq_along(x))
      x[[j]] <- tranformBayesX(x[[j]])
  } else x <- randomize(x)
  x
}


setupBayesX <- function(x, n.iter = 1200, thin = 1, burnin = 200, seed = NULL, dir = NULL, ...)
{
  x$call <- NULL
  args <- list(...)
  name <- if(is.null(args$name)) "bayesr" else args$name
  data.name <- if(is.null(args$data.name)) "d" else args$data.name

  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  prg <- paste('% usefile', file.path(dir, paste(name, 'prg', sep = '.')))
  prg <- c(prg, paste('logopen using', file.path(dir, paste(name, 'log', sep = '.'))))

  count <- 1
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        prg <- c(prg,
          paste('dataset ', data.name, count, sep = ''),
          paste(data.name, count, '.infile using ')
        )
      }
    } else {
      prg <- c(prg,
        paste('dataset ', data.name, count, sep = ''),
        paste(data.name, count, '.infile using')
      )
    }
  }

  prg <- c(prg, paste('mcmcreg', name))

  

#dataset d1  
#d1.infile using /tmp/Rtmpybbo5u/bayesx/bayesx.estim_hlevel1_MAIN_REGRESSION.data.raw 
# 
#b.outfile = /tmp/Rtmpybbo5u/bayesx/bayesx.estim_hlevel1_MAIN_REGRESSION 
 
  eqs <- list()
  for(j in 1:n) {
    if("fake.formula" %in% names(x[[j]])) {
      prg2 <- paste("b.hregress", x[[j]]$response)
      if(all(x[[j]]$X[, 1] == 1))
        prg2 <- paste(prg2, "=", "const")
      if(length(x[[j]]$pterms))
        prg2 <- paste(prg2, "+", paste(x[[j]]$pterms, collapse = " + "))
    } else {

    }
    
  }
print(prg2)
stop()
}


samplerBayesX <- function(x, ...) { x }
resultsBayesX <- function(x, ...) { x }


## Run BayesX .prg files.
run.bayesx2 <- function(prg = NULL, verbose = TRUE, bin = NULL, ...)
{
  os.win <- .Platform$OS.type == "windows"
  bin <- if(is.null(dir)) {
    system.file(package = "BayesXsrc", "libs", .Platform$r_arch, if(os.win) "BayesX.exe" else "BayesX")
  } else path.expand(bin)
  bin <- shQuote(bin)
  if(is.null(prg)) {
    output <- file.exists("output")
    temp <- file.exists("temp")
    if(!os.win)
      bin <- paste(shQuote(file.path(R.home(),"bin","R")), "CMD", bin)
    system(bin, ...)
    if(!output && file.exists("output"))
      try(suppressWarnings(file.remove("output")), silent = TRUE)
    if(!temp && file.exists("temp"))
      try(suppressWarnings(file.remove("temp")), silent = TRUE)
    return(invisible(NULL))
  } else {
    fileo <- getwd()
    prg <- path.expand(prg)
    dir <- dirname(prg)
    prg.name <- basename(prg)
    setwd(dir)
    output <- file.exists("output")
    temp <- file.exists("temp")
    ptm <- proc.time()
    if(!verbose) {
      if(os.win) {
        log <- try(system(paste(bin, " ", dir, "/", prg.name, sep = ""), 
          intern = TRUE, show.output.on.console = FALSE, ...))
      } else {
        log <- try(system(paste(bin, " ", dir, "/", prg.name, 
          " > ", dir, "/bayesx.log", sep = ""), intern = FALSE, ...))
      }
    } else log <- try(system(paste(bin, " ", dir, "/", prg.name, sep = ""), intern = FALSE, ...))
	  if(log != 0 && !os.win)
      warning("problem processing BayesX!")
    if(!verbose && .Platform$OS.type == "unix")	
      log <- readLines(paste(dir, "/bayesx.log", sep = ""))
    now <- proc.time()
    runtime <- now - ptm
    runtime <- runtime
    if(verbose)
      cat("Total run time was:", runtime[3L], "sec\n")
    if(!output && file.exists("output"))
      try(suppressWarnings(file.remove("output")), silent = TRUE)
    if(!temp && file.exists("temp"))
      try(suppressWarnings(file.remove("temp")), silent = TRUE)
    try(setwd(fileo), silent = TRUE)
    logerror <- FALSE
    for(logline in log)
      if(grepl("ERROR", logline))
        logerror <- TRUE
    if(logerror)
      warning("an error occurred during runtime of BayesX, please check the BayesX logfile!")
    return(invisible(list(log = log, runtime = runtime)))
  }
}

#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, chains = 1, combine = TRUE, n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, dir = NULL, model.name = NULL, ...)
{
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    on.exit(unlink(dir))
  } else dir <- path.expand(dir)
  if(is.null(model.name))
    model.name <- 'bayesr'
  prg.name <- paste(model.name, 'prg', sep = '.')
  data.name <- if(!is.null(data)) {
    deparse(substitute(data), backtick = TRUE, width.cutoff = 500)
  } else "d"

  cat('%% BayesX program created by BayesR: ', as.character(Sys.time()), '\n',
    file = file.path(dir, prg.name), sep = '')
  cat('%% usefile ', file.path(dir, prg.name), '\n\n', 
    file = file.path(dir, prg.name), append = TRUE, sep = '')

  transform <- function(x) { transformBayesX(x, dir = dir, prg.name = prg.name, ...) }
  setup <- function(x) {
    setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
      seed = seed, model.name = model.name, data.name = data.name,
      dir = dir, prg.name = prg.name, ...)
  }
  sampler <- function(x) { samplerBayesX(x, ...) }
  results <- function(x) { resultsBayesX(x, ...) }

  bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transform,
    setup = setup, sampler = sampler, results = results,
    cores = cores, combine = combine, ...)
}


####################################
## (2) BayesX specific functions. ##
####################################
transformBayesX <- function(x, ...)
{
  args <- list(...)
  dir <- args$dir
  prg.name <- args$prg.name
  
  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    for(j in seq_along(x))
      x[[j]] <- transformBayesX(x[[j]], dir = dir, prg = prg.name, ...)
  } else {
    x <- randomize(x)
    if(length(x$smooth)) stop("arbitrary smooths not supported yet!")
    if(length(x$sx.smooth)) {
      for(j in seq_along(x$sx.smooth)) {
        tt <- x$sx.smooth[[j]]
        nx <- colnames(x$X)
        x$X <- cbind(x$X, x$mf[, tt$term], if(tt$by != "NA") x$mf[, tt$by] else NULL)
        colnames(x$X) <- c(nx, tt$term, if(tt$by != "NA") tt$by else NULL)
        x$sx.smooth[[j]] <- bayesx.construct(object = tt, dir = dir, prg = prg.name, data = x$mf)
      }
    }
    x$X <- cbind(x$mf[, x$response], x$X)
    colnames(x$X)[1] <- x$response
    x$X <- x$X[, !grepl('(Intercept)', colnames(x$X), fixed = TRUE)]
    x$X <- x$X[, !duplicated(colnames(x$X))]
  }

  x
}


setupBayesX <- function(x, n.iter = 1200, thin = 1, burnin = 200, seed = NULL, ...)
{
  x$call <- NULL

  args <- list(...)
  model.name <- args$model.name
  data.name <- args$data.name
  dir <- args$dir
  prg.name <- args$prg.name

  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  prg <- paste('logopen using', file.path(dir, paste(model.name, 'log', sep = '.')))

  count <- 1
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        dpath <- file.path(dir, paste(paste(data.name, count, sep = ''), "raw", sep = '.'))
        write.table(x[[j]][[i]]$X, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
        prg <- c(prg,
          paste('dataset ', data.name, count, sep = ''),
          paste(data.name, count, '.infile using ', dpath, sep = '')
        )
        count <- count + 1
      }
    } else {
      dpath <- file.path(dir, paste(paste(data.name, count, sep = ''), "raw", sep = '.'))
      write.table(x[[j]]$X, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
      prg <- c(prg,
        paste('dataset ', data.name, count, sep = ''),
        paste(data.name, count, '.infile using', dpath, sep = '')
      )
      count <- count + 1
    }
  }

  prg <- c(prg, paste('\nmcmcreg', model.name))
  prg <- c(prg, paste(model.name, '.outfile = ', file.path(dir, model.name), '\n', sep = ''))
 
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        eqn <- paste(model.name, '.hregress ', x[[j]][[i]]$response, ' = ', sep = '')
        if(length(x[[j]][[i]]$pterms))
          eqn <- paste(eqn, paste(x[[j]][[i]]$pterms, collapse = ' + '))
        if(length(x[[j]][[i]]$sx.smooth))
          eqn <- paste(eqn, paste(unlist(x[[j]][[i]]$sx.smooth), collapse = ' + '))
        prg <- c(prg, eqn)
      }
    } else {
      eqn <- paste(model.name, '.hregress ', x[[j]]$response, ' = ', sep = '')
      if(length(x[[j]]$pterms))
        eqn <- paste(eqn, paste(x[[j]]$pterms, collapse = ' + '))
      if(length(x[[j]][[i]]$sx.smooth))
        eqn <- paste(eqn, paste(unlist(x[[j]]$sx.smooth), collapse = ' + '))
      prg <- c(prg, eqn)
    }
  }

  cat(paste(prg, collapse = '\n'), file = file.path(dir, prg.name), append = TRUE)
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

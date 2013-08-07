#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, chains = 1, combine = TRUE, n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, dir = "~/tmp", model.name = NULL, ...)
{
  require("BayesXsrc")

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
  results <- function(x, ...) { resultsBayesX(x, ...) }

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
  x <- tBayesX(x, ...)
  call <- x$call; x$call <- NULL

  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  get_family <- function(x) {
    if(is.null(x$family))
      for(j in x) {
        family <- get_family(j)
    } else family <- x$family
    family
  }

  family <- get_family(x)
  family <- if(is.function(family)) family() else family
  stopifnot(!is.null(family$BayesX))
  names(x) <- rep(family$names, length.out = n)
  x$call <- call
  x$family <- family

  x
}

tBayesX <- function(x, ...)
{
  args <- list(...)
  dir <- args$dir
  prg.name <- args$prg.name
  if(inherits(x, "list") & !any(c("smooth", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx)) nx <- 1:length(x)
    if(length(unique(nx)) < length(x)) nx <- 1:length(x)
    for(j in nx)
      x[[j]] <- tBayesX(x[[j]], dir = dir, prg = prg.name, ...)
  } else {
    x <- randomize(x)
    if(length(x$smooth)) stop("arbitrary smooths not supported yet!")
    if(length(x$sx.smooth)) {
      for(j in seq_along(x$sx.smooth)) {
        x$sx.smooth[[j]] <- bayesx.construct(object = x$sx.smooth[[j]], dir = dir,
          prg = prg.name, data = x$mf)
      }
    }
  }

  x
}


BayesX.control <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "full", ...)
{
  if(is.null(seed))
    seed <- round(runif(1L) * .Machine$integer.max)
  stopifnot(burnin < n.iter)
  cvals <- list("iterations" = n.iter, "burnin" = burnin, "step" = thin, "setseed" = seed,
    "predict" = predict)
  attr(cvals, "main") <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  cvals
}

setupBayesX <- function(x, control = BayesX.control(...), ...)
{
  family <- x$family
  x$call <- x$family <- NULL

  args <- list(...)
  model.name <- args$model.name
  data.name <- args$data.name
  dir <- args$dir
  prg.name <- args$prg.name
  prg <- paste('logopen using', file.path(dir, paste(model.name, 'log\n', sep = '.')))

  BayesX_data <- function(x) {
    X <- NULL
    if(!is.null(x$X)) {
      if(ncol(x$X) > 0)
        X <- as.data.frame(x$X)
    }
    if(k2 <- length(x$sterms)) {
      X <- if(!is.null(X)) {
        cbind(X, as.data.frame(x$mf[, x$sterms]))
      } else as.data.frame(x$mf[, x$sterms])
      k1 <- ncol(X)
      colnames(X)[(k1 - k2 + 1):k1] <- x$sterms
    }
    X <- X[, !grepl("(Intercept)", colnames(X), fixed = TRUE), drop = FALSE]
    X[[x$response]] <- x$mf[[x$response]]
    X <- as.data.frame(X)
    for(j in 1:ncol(X)) {
      if(is.factor(X[, j])) {
        X[, j] <- as.integer(X[, j])
        X <- X[order(X[, j]), ]
      }
    }
    return(unique(X))
  }

  nx <- names(x)
  n <- length(x)

  count <- 1
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        dpath <- file.path(dir, paste(dname <- paste(data.name, count, sep = ''), "raw", sep = '.'))
        d <- BayesX_data(x[[j]][[i]])
        x[[j]][[i]]$hlevel <- i
        if(!is.null(d)) {
          x[[j]][[i]]$dname <- dname
          write.table(d, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
          prg <- c(prg,
            paste('dataset ', data.name, count, sep = ''),
            paste(data.name, count, '.infile using ', dpath, sep = '')
          )
          count <- count + 1
        }
      }
    } else {
      dpath <- file.path(dir, paste(dname <- paste(data.name, count, sep = ''), "raw", sep = '.'))
      d <- BayesX_data(x[[j]])
      if(!is.null(d)) {
        x[[j]]$dname <- dname
        write.table(d, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
        prg <- c(prg,
          paste('dataset ', data.name, count, sep = ''),
          paste(data.name, count, '.infile using ', dpath, sep = '')
        )
        count <- count + 1
      }
    }
  }

  prg <- c(prg, paste('\nmcmcreg', model.name))
  prg <- c(prg, paste(model.name, '.outfile = ', file.path(dir, model.name), '\n', sep = ''))

  make_eqn <- function(x, ctr = TRUE) {
    n <- length(x)
    eqn <- NULL
    for(j in n:1) {
      if(!"fake.formula" %in% names(x[[j]])) {
        eqn <- c(eqn, make_eqn(x[[j]], ctr))
      } else {
        teqn <- paste(model.name, '.hregress ', x[[j]]$response, sep = '')
        et <- x[[j]]$pterms
        fctr <- attr(x[[j]]$formula, "control")
        if(is.null(fctr)) fctr <- ""
        if(attr(terms(x[[j]]$formula), "intercept"))
          et <- c("const", et)
        if(length(x[[j]]$sx.smooth))
          et <- c(et, x[[j]]$sx.smooth)
        if(length(et))
          teqn <- paste(teqn, '=', paste(et, collapse = ' + '))
        if(ctr) {
          ok <- attr(control, "main")
          c2 <- control[!ok]
          ok <- sapply(names(c2), function(x) { !any(grepl(x, fctr)) })
          c2 <- c2[ok]
          if(length(c2)) {
            teqn <- paste(teqn, ", ", paste(names(c2), "=", unlist(c2),
              sep = "", collapse = " "), sep = "")
          }
          ctr <- FALSE
        } else teqn <- paste(teqn, ",", sep = "")
        ok <- TRUE
        if(!is.null(x[[j]]$hlevel)) {
          if(!any(grepl("hlevel", fctr))) {
            teqn <- paste(teqn, " hlevel=", x[[j]]$hlevel, sep = "")
            if(x[[j]]$hlevel > 1) {
              if(!any(grepl("family", fctr))) {
                teqn <- paste(teqn, " family=", family$BayesX$h, sep = "")
                ok <- FALSE
              }
            }
          } else ok <- grepl("1", fctr[grepl("hlevel", fctr)])
        }
        if(ok) {
          if(!any(grepl("family", fctr)))
            teqn <- paste(teqn, " family=", family$BayesX[[nx[j]]][1], sep = "")
          if(!any(grepl("equationtype", fctr)))
            teqn <- paste(teqn, " equationtype=", family$BayesX[[nx[j]]][2], sep = "")
        }
        if(length(fctr)) {
          fctr <- paste(fctr, collapse = " ")
          if(nchar(fctr) > 0)
            teqn <- paste(teqn, fctr)
        }
        if(!is.null(x[[j]]$dname)) {
          if(!any(grepl("using", fctr)))
            teqn <- paste(teqn, " using ", x[[j]]$dname, sep = "")
        }
        eqn <- c(eqn, "", teqn)
      }
    }
    eqn
  }

  prg <- c(prg, make_eqn(x))

  prg <- gsub("(random", "(hrandom", prg, fixed = TRUE)
  for(i in 1:5)
    prg <- gsub(paste("(psplinerw", i, sep = ""), "(pspline", prg, fixed = TRUE)
  prg <- c(prg, "", paste(model.name, "getsample", sep = "."))

  cat(paste(prg, collapse = '\n'), file = file.path(dir, prg.name), append = TRUE)
  x$prg <- file.path(dir, prg.name)
  x$model.name <- model.name
  x$dir <- dir

  x
}


samplerBayesX <- function(x, verbose = FALSE, ...)
{
  warn <- getOption("warn")
  options(warn = -1)
  ok <- BayesXsrc::run.bayesx(prg = x$prg, verbose = verbose, ...)
  options("warn" = warn)
  if(length(i <- grep("error", ok$log, ignore.case = TRUE))) {
    errl <- gsub("^ +", "", ok$log[i])
    errl <- gsub(" +$", "", errl)
    errl <- encodeString(errl, width = NA, justify = "left")
    errl <- paste(" *", errl)
    errm <- paste("an error occurred running the BayesX binary! The following messages are returned:\n\n",
      paste(errl, collapse = "\n", sep = ""), sep = "")
    stop(errm)
  }
  samples <- NULL
  tspecs <- grep(paste(x$model.name, "_R.r", sep = ""), dir(x$dir), value = TRUE, fixed = TRUE)
  if(length(tspecs)) {
    tspecs <- readLines(file.path(x$dir, tspecs))
    tspecs <- process_tspecs(tspecs)
    for(j in 1:nrow(tspecs)) {
      ts <- read.table(tspecs[j, "samples"], header = TRUE)
      ts$intnr <- NULL
      names(ts) <- paste(tspecs[j, "terms"], "_p[", 1:ncol(ts), "]", sep = "")
      samples <- cbind(samples, as.matrix(ts))
    }
    samples <- as.mcmc.list(list(as.mcmc(samples)))
    attr(samples, "tspecs") <- tspecs
  }
  samples
}

process_tspecs <- function(x)
{
  n <- length(x)
  samples <- eqntype <- family <- basis <- terms <- hlevel <- rep(NA, length = n)
  for(j in 1:n) {
    specs <- strsplit(x[j], ",")[[1]]
    samples[j] <- gsub("\\s", "", strsplit(grep("pathsamples", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
    eqntype[j] <- gsub("\\s", "", strsplit(grep("equationtype", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
    family[j] <- gsub("\\s", "", strsplit(grep("family", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
    basis[j] <- gsub("\\s", "", strsplit(grep("pathbasis", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
    hlevel[j] <- gsub("\\s", "", strsplit(grep("hlevel", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
    terms[j] <- gsub("\\s", "", strsplit(grep("term", specs, value = TRUE),
      "=", fixed = TRUE)[[1]][2])
  }
  if(any(i <- duplicated(terms))) {
    for(j in terms[i]) {
      tt <- grep(j, terms, value = TRUE, fixed = TRUE)
      tt <- paste(tt, 1:length(tt), sep = ":")
      terms[terms == j] <- tt
    }
  }
  rval <- cbind(terms, samples, basis, family, eqntype, hlevel)
  rval
}

resultsBayesX <- function(x, samples, ...)
{
  print(names(x))
  samples
}


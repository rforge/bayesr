#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian.BayesX, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, combine = TRUE, n.iter = 1200, thin = 1, burnin = 200, seed = NULL, ...)
{
  require("BayesXsrc")

  data.name <- if(!is.null(data)) {
    deparse(substitute(data), backtick = TRUE, width.cutoff = 500)
  } else "d"

  setup <- function(x) {
    setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
      seed = seed, data.name = data.name, cores = cores, ...)
  }

  bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transformBayesX,
    setup = setup, sampler = samplerBayesX, results = resultsBayesX,
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
  names(x) <- rep(family$names, length.out = n)
  x$call <- call
  x$family <- family

  x
}

tBayesX <- function(x, ...)
{
  if(inherits(x, "list") & !any(c("smooth", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx)) nx <- 1:length(x)
    if(length(unique(nx)) < length(x)) nx <- 1:length(x)
    for(j in nx)
      x[[j]] <- tBayesX(x[[j]], ...)
  } else {
    x <- randomize(x)
    if(length(x$smooth)) stop("arbitrary smooths not supported yet!")
    if(length(x$sx.smooth)) {
      for(j in seq_along(x$sx.smooth)) {
        x$sx.smooth[[j]] <- sx.construct(x$sx.smooth[[j]], x$mf)
      }
    }
  }

  x
}


controlBayesX <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "full", model.name = "bayesr", data.name = "d",
  prg.name = NULL, dir = NULL, verbose = FALSE, cores = NULL, ...)
{
  if(is.null(seed))
    seed <- '##seed##'
  stopifnot(burnin < n.iter)
  if(is.null(model.name))
    model.name <- 'bayesr'
  if(is.null(data.name))
    data.name <- 'd'
  if(is.null(prg.name))
    prg.name <- paste(model.name, 'prg', sep = '.')
  if(!grepl(".prg", prg.name))
    prg.name <- paste(prg.name, "prg", sep = ".")
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    attr(dir, "unlink") <- TRUE
  } else dir <- path.expand(dir)
  if(is.null(cores)) cores <- 1

  cvals <- list(
    "prg" = list(
      "iterations" = n.iter, "burnin" = burnin, "step" = thin,
      "setseed" = seed, "predict" = predict
    ),
    "setup" = list(
      "main" = rep(FALSE, 5), "model.name" = model.name, "data.name" = data.name,
      "prg.name" = prg.name, "dir" = dir, "verbose" = verbose, "cores" = cores
    )
  )

  cvals
}

setupBayesX <- function(x, control = controlBayesX(...), ...)
{
  family <- x$family
  x$call <- x$family <- NULL

  args <- list(...)
  model.name <- control$setup$model.name
  data.name <- control$setup$data.name
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  prg <- paste('logopen using', file.path(dir, paste(model.name, 'log', sep = '.')))

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
  prg <- c(prg, paste(model.name, '.outfile = ',
    if(cores < 2) file.path(dir, model.name) else '##outfile##',
    sep = ''))

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
          ok <- control$setup$main
          c2 <- control$prg[!ok]
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
                teqn <- paste(teqn, " family=", family$h, sep = "")
                ok <- FALSE
              }
            }
          } else ok <- grepl("1", fctr[grepl("hlevel", fctr)])
        }
        if(ok) {
          if(!any(grepl("family", fctr)))
            teqn <- paste(teqn, " family=", family[[nx[j]]][1], sep = "")
          if(!any(grepl("equationtype", fctr)))
            teqn <- paste(teqn, " equationtype=", family[[nx[j]]][2], sep = "")
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
  prg <- c(
    paste('%% BayesX program created by BayesR: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )

  x$prg <- prg
  x$control <- control

  x
}


samplerBayesX <- function(x, ...)
{
  if(any(grepl("##seed##", x$prg, fixed = TRUE)))
    x$prg <- gsub("##seed##", round(runif(1L) * .Machine$integer.max), x$prg, fixed = TRUE)

  dir <- x$control$setup$dir
  prg.name <- x$control$setup$prg.name
  model.name <- x$control$setup$model.name
  cores <- x$control$setup$cores

  if(cores > 1) {
    k <- 1
    while(file.exists(cdir <- file.path(dir, paste("Core", k, sep = "")))) {
      k <- k + 1
    }
    dir.create(cdir)
    x$prg <- gsub('##outfile##', prg <- file.path(cdir, model.name), x$prg, fixed = TRUE)
    prg.name <- paste(gsub(".prg", "", prg.name, fixed = TRUE), k, ".prg", sep = "")
    if(length(i <- grep('usefile', x$prg)))
      x$prg[i] <- paste('%% usefile', file.path(dir, prg.name))
    x$prg <- gsub(file.path(dir, paste(model.name, 'log', sep = '.')),
      file.path(cdir, paste(model.name, 'log', sep = '.')), x$prg, fixed = TRUE)
    x$prg <- gsub(paste(model.name, ".", sep = ""), paste(model.name, k, ".", sep = ""),
      x$prg, fixed = TRUE)
    x$prg[grep("mcmcreg", x$prg)] <- paste("mcmcreg ", model.name, k, sep = "")

    dir <- cdir
  }

  prg <- file.path(dir, prg.name)
  writeLines(x$prg, file.path(dir, prg.name))

  warn <- getOption("warn")
  options(warn = -1)
  ok <- BayesXsrc::run.bayesx(prg = prg, verbose = x$control$setup$verbose, ...)
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
  tspecs <- grep(paste(x$control$setup$model.name, "_R.r", sep = ""), dir(dir), value = TRUE, fixed = TRUE)
  if(length(tspecs)) {
    tspecs <- readLines(file.path(dir, tspecs))
    tspecs <- process_tspecs(tspecs)
    for(j in 1:nrow(tspecs$effects)) {
      ts <- read.table(tspecs$effects[j, "samples"], header = TRUE)
      ts$intnr <- NULL
      names(ts) <- paste(tspecs$effects[j, "terms"], "[", 1:ncol(ts), "]", sep = "")
      samples <- cbind(samples, as.matrix(ts))
    }
    samples <- as.mcmc.list(list(as.mcmc(samples)))
    attr(samples, "tspecs") <- tspecs
  }
  samples
}

process_tspecs <- function(x)
{
  predict <- dic <- NULL
  if(length(i <- grep("predict.res", x, ignore.case = TRUE))) {
    predict <- gsub("\\s", "", strsplit(x[i], "=", fixed = TRUE)[[1]][2])
    x <- x[-i]
  }
  if(length(i <- grep("DIC.res", x, ignore.case = TRUE))) {
    dic <- gsub("\\s", "", strsplit(x[i], "=", fixed = TRUE)[[1]][2])
    x <- x[-i]
  }
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
    terms[j] <- paste(terms[j], if(any(grepl("_", family[j]))) {
      strsplit(family[j], "_", fixed = TRUE)[[1]][2]
    } else family[j], sep = ":")
  }
  if(any(i <- duplicated(terms))) {
    for(j in terms[i]) {
      tt <- grep(j, terms, value = TRUE, fixed = TRUE)
      tt <- paste(tt, 1:length(tt), sep = ":")
      terms[terms == j] <- tt
    }
  }
  rval <- list("effects" = cbind(terms, samples, basis, family, eqntype, hlevel),
    "model" = c(predict, dic))
  rval
}

resultsBayesX <- function(x, samples, id = NULL, ...)
{
  tspecs <- attr(samples, "tspecs")
  if(!all(c("family", "formula") %in% names(x))) {
    nx <- names(x)
    nx <- nx[!(nx %in% c("call", "family"))]
    rval <- list()
    family <- x[[1]]$family
    if(is.function(family))
      family <- family()
    fn <- family$names
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- resultsBayesX(x[[nx[j]]], samples, id = fn[j])
      if(!is.null(rval[[nx[j]]]$effects)) {
        for(i in seq_along(rval[[nx[j]]]$effects)) {
          specs <- attr(rval[[nx[j]]]$effects[[i]], "specs")
          specs$label <- paste(specs$label, fn[j], sep = ":")
          attr(rval[[nx[j]]]$effects[[i]], "specs") <- specs
        }
        names(rval[[nx[j]]]$effects) <- paste(names(rval[[nx[j]]]$effects), fn[j], sep = ":")
      }
    }
    names(rval) <- fn
    class(rval) <- "bayesr"
    return(rval)
  } else {
    chains <- length(samples)
    rval <- vector(mode = "list", length = chains)
    snames <- colnames(samples[[1]])
print(snames)
print(chains)
stop()
  }
}


########################################
## (3) BayesX model term construction ##
########################################
sx.construct <- function(object, data) 
{
  UseMethod("sx.construct")
}

sx.construct.default <- function(object, data) 
{
  cl <- class(object)
  bs <- gsub(".smooth.spec", "", cl, fixed = TRUE)
  info <- if(cl == paste(bs, "smooth.spec", sep = ".")) {
    paste("with basis", sQuote(bs))
  } else {
    paste("of class", sQuote(cl))
  }
  stop(paste("BayesX does not support smooth terms ", info,
    ", it is recommended to use sx() for specifying smooth terms",
    sep = ""))
}

do.xt <- function(term, object, not = NULL, noco = FALSE)
{
  if(!is.null(object$xt)) {
    names.xt <- names(object$xt)
    if(is.null(not))
      not <- "not"
    count <- 1
    co <- ","
    if(noco)
      co <- NULL
    for(name in names.xt) {
      if(count > 1)
        co <- ","
      if(!name %in% not) {
        if(name %in% c("full", "catspecific", "center", "derivative", "nofixed") || 
          is.logical(object$xt[[name]])) {
          if(is.logical(object$xt[[name]])) {
            if(object$xt[[name]])
              term <- paste(term, co, name, sep = "")
          } else term <- paste(term, co, name, sep = "")
          count <- count + 1
        } else {
          term <- paste(term, co, name, "=", object$xt[name], sep = "")
          count <- count + 1
        }
      }
    }
  }

  return(term)
}

sx.construct.ps.smooth.spec <- sx.construct.psplinerw1.smooth.spec <-
sx.construct.psplinerw2.smooth.spec <- sx.construct.pspline.smooth.spec <-
function(object, data)
{
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 1L
  if(class(object) == "psplinerw1.smooth.spec")
    object$p.order[2L] <- 1L
  if(class(object) == "psplinerw2.smooth.spec")
    object$p.order[2L] <- 2L
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
  if(length(object$p.order) > 1L) {
    if(object$p.order[2L] > 2L) {
      warning("order of the difference penalty not supported by BayesX, set to 2!")
      object$p.order <- c(object$p.order[1L], 2L)
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  if(nrknots < 5L) {
    warning("number of inner knots smaller than 5 not supported by BayesX, set to 5!",
      call. = FALSE)
    nrknots <- 5L
  }
  termo <- object$term
  term <- paste(termo, "(psplinerw", object$p.order[2L], ",nrknots=",
    nrknots, ",degree=", object$p.order[1L], sep = "")
  object$xt[c("knots", "nrknots", "degree")] <- NULL
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)
  
  return(term)
}

sx.construct.ra.smooth.spec <- sx.construct.re.smooth.spec <-
sx.construct.random.smooth.spec <- function(object, data)
{
  term <- object$term
  if(is.null(object$ins))
    term <- paste(term, "(random", sep = "")
  else
    term <- paste(term, "(hrandom", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.rsps.smooth.spec <- function(object, data)
{
  term <- paste(object$by, "(hrandom_pspline,centermethod=meansum2", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  term <- paste(object$term, term , sep = "*")

  return(term)
}


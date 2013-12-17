#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, combine = TRUE, n.iter = 12000, thin = 10, burnin = 2000, seed = NULL, ...)
{
  require("BayesXsrc")

  data.name <- if(!is.null(data)) {
    deparse(substitute(data), backtick = TRUE, width.cutoff = 500)
  } else "d"

  setup <- function(x) {
    setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
      seed = seed, data.name = data.name, cores = cores, ...)
  }

  family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)

  rval <- bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transformBayesX,
    setup = setup, sampler = samplerBayesX, results = resultsBayesX,
    cores = cores, combine = combine, sleep = 1, ...)
  
  attr(rval, "call") <- match.call()
  
  rval
}


####################################
## (2) BayesX specific functions. ##
####################################
transformBayesX <- function(x, ...)
{
  family <- bayesr.family(attr(x, "family"))
  call <- x$call; x$call <- NULL

  if(family$cat) {
    ylevels <- attr(x, "ylevels")
    reference <- attr(x, "reference")
    grid <- attr(x, "grid")
    rn <- attr(attr(x, "model.frame"), "response.name")
    f <- as.formula(paste("~ -1 +", rn))
    rm <- as.data.frame(model.matrix(f, data = attr(x, "model.frame")))
    cn <- colnames(rm) <- rmf(colnames(rm))
    rm <- rm[, !grepl(reference, cn), drop = FALSE]
    mf <- cbind(attr(x, "model.frame"), rm)
    attr(mf, "response.name") <- rn
    x <- x[ylevels != reference]
    attr(x, "ylevels") <- ylevels
    attr(x, "reference") <- reference
    attr(x, "model.frame") <- mf
    attr(x, "grid") <- grid
  }

  names(attr(x, "model.frame")) <- rmf(names(attr(x, "model.frame")))
  attr(attr(x, "model.frame"), "response.name") <- rmf(attr(attr(x, "model.frame"), "response.name"))

  tBayesX <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- tBayesX(obj[[j]], ...)
    } else {
      obj <- randomize(obj)
      obj$response <- rmf(obj$response)
      if(!is.null(dim(obj$X)))
        colnames(obj$X) <- rmf(colnames(obj$X))
      if(length(obj$pterms))
        obj$pterms <- rmf(obj$pterms)
      if(length(obj$sterms))
        obj$sterms <- rmf(obj$sterms)

      if(length(obj$smooth)) stop("arbitrary smooths not supported yet!")
      if(length(obj$sx.smooth)) {
        for(j in seq_along(obj$sx.smooth)) {
          obj$sx.smooth[[j]]$term <- rmf(obj$sx.smooth[[j]]$term)
          tmp <- sx.construct(obj$sx.smooth[[j]],
            if(is.null(attr(obj, "model.frame"))) {
              attr(x, "model.frame")
            } else attr(obj, "model.frame"))
          attr(tmp, "specs") <- obj$sx.smooth[[j]]
          obj$sx.smooth[[j]] <- tmp
        }
      }
    }
    obj
  }

  x <- tBayesX(x, ...)

  if(inherits(x, "bayesr.input") & !any(c("formula", "fake.formula", "response") %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  names(x) <- rep(family$names, length.out = n)
  if(!is.null(family$order)) {
    mf <- attr(x, "model.frame")
    class(x) <- "list"
    x <- x[family$order]
    class(x) <- c("bayesr.input", "list")
    attr(x, "model.frame") <- mf; rm(mf)
  }
  attr(x, "call") <- call
  attr(x, "family") <- family

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
  if(!file.exists(dir)) dir.create(dir)
  if(is.null(cores)) cores <- 1

  cvals <- list(
    "prg" = list(
      "iterations" = n.iter, "burnin" = burnin, "step" = thin,
      "setseed" = seed, "predict" = predict
    ),
    "setup" = list(
      "main" = c(rep(FALSE, 4), TRUE), "model.name" = model.name, "data.name" = data.name,
      "prg.name" = prg.name, "dir" = dir, "verbose" = verbose, "cores" = cores
    )
  )

  cvals
}

setupBayesX <- function(x, control = controlBayesX(...), ...)
{
  names(attr(x, "model.frame")) <- rmf(names(attr(x, "model.frame")))

  family <- attr(x, "family")
  x$call <- x$family <- NULL

  args <- list(...)
  model.name <- control$setup$model.name
  data.name <- rmf(control$setup$data.name)
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  BayesX_data <- function(obj, h = FALSE) {
    if(!is.null(obj$cat.formula) & !h) {
      obj$response <- formula_respname(obj$cat.formula)
      obj$response.vec <- attr(x, "model.frame")[[obj$response]]
    }
    X <- NULL
    if(!is.null(obj$X)) {
      if(ncol(obj$X) > 0) {
        X <- obj$X
        if(length(i <- grep("Intercept", colnames(X), fixed = TRUE)))
          X <- obj$X[, -i, drop = FALSE]
        X <- if(ncol(X) > 0) as.data.frame(X) else NULL
      }
    }
    if(k2 <- length(obj$sterms)) {
      X <- if(!is.null(X)) {
        cbind(X, attr(x, "model.frame")[, obj$sterms, drop = FALSE])
      } else {
        attr(x, "model.frame")[, obj$sterms, drop = FALSE]
      }
      k1 <- ncol(X)
      colnames(X)[(k1 - k2 + 1):k1] <- obj$sterms
    }
    if(!is.null(X)) {
      if(ncol(X) > 0)
        X <- X[, !grepl("Intercept", colnames(X), fixed = TRUE), drop = FALSE]
      if(!is.null(obj$response)) {
        X[[obj$response]] <- obj$response.vec
      }
    } else {
      X <- list()
      if(!is.null(obj$response))
        X[[obj$response]] <- obj$response.vec
      X <- as.data.frame(X)
    }
    if(ncol(X) < 1) {
      X <- NULL
    } else {
      yf <- if(is.null(obj$response)) FALSE else is.factor(X[[obj$response]])
      X <- as.data.frame(X)
      for(j in 1:ncol(X)) {
        if(is.factor(X[, j])) {
          warn <- getOption("warn")
          options(warn = -1)
          tx <- as.integer(as.character(X[, j]))
          options("warn" = warn)
          X[, j] <- if(!any(is.na(tx))) tx else as.integer(X[, j])
          X <- X[order(X[, j]), , drop = FALSE]
        }
      }
      if(yf) {
        if(!h) {
          if(!is.null(obj$response))
            X[[obj$response]] <- as.integer(as.factor(X[[obj$response]])) - 1
        }
        X <- X[order(X[[obj$response]]), , drop = FALSE]
      }
    }
    if(h) X <- unique(X)

    return(X)
  }

  nx <- names(x)
  n <- length(x)

  count <- 1
  d <- prg <- NULL
  dpath0 <- file.path(dir, paste(dname0 <- paste(data.name, 0, sep = ''), "raw", sep = '.'))
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in 1:length(x[[j]])) {
        d2 <- NULL
        if(i < 2) {
          d2 <- BayesX_data(x[[j]][[i]])
          if(is.null(d)) {
            d <- d2
          } else {
            if(!is.null(d2))
              d <- cbind(d, d2)
          }
          d2 <- NULL
          d <- d[, unique(names(d)), drop = FALSE]
          x[[j]][[i]]$dname <- dname0
        } else d2 <- BayesX_data(x[[j]][[i]], i > 1)
        x[[j]][[i]]$hlevel <- i
        if(!is.null(d2)) {
          dpath <- file.path(dir, paste(dname <- paste(data.name, count, sep = ''),
            "raw", sep = '.'))
          x[[j]][[i]]$dname <- dname
          write.table(d2, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
          prg <- c(prg,
            paste('dataset ', data.name, count, sep = ''),
            paste(data.name, count, '.infile using ', dpath, sep = '')
          )
          count <- count + 1
        }
      }
    } else {
      d2 <- BayesX_data(x[[j]])
      if(is.null(d)) {
        d <- d2
      } else {
        if(!is.null(d2))
          d <- cbind(d, BayesX_data(x[[j]]))
      }
      d2 <- NULL
      d <- d[, unique(names(d)), drop = FALSE]
      x[[j]]$dname <- dname0
      x[[j]]$hlevel <- 1
    }
  }

  write.table(d, file = dpath0, quote = FALSE, row.names = FALSE, col.names = TRUE)
  prg <- c(
    paste('dataset', dname0),
    paste(dname0, '.infile using ', dpath0, sep = ''),
    prg
  )

  prg <- c(paste('logopen using', file.path(dir, paste(model.name, 'log', sep = '.'))), "", prg)
  prg <- c(prg, paste('\nmcmcreg', model.name))
  prg <- c(prg, paste(model.name, '.outfile = ',
    if(cores < 2) file.path(dir, model.name) else '##outfile##',
    sep = ''))

  response.name <- attr(attr(x, "model.frame"), "response.name")

  make_eqn <- function(x, ctr = TRUE, id = NULL) {
    n <- length(x)
    eqn <- NULL
    for(j in n:1) {
      ctr2 <- if(j != n) FALSE else TRUE
      if(!is.null(id)) ctr2 <- id != n
      if(!"fake.formula" %in% names(x[[j]])) {
        eqn <- c(eqn, make_eqn(x[[j]], ctr, id = j))
        ctr <- FALSE
      } else {
        teqn <- paste(model.name, '.hregress ',
          if(family$cat) {
            if(x[[j]]$hlevel > 1) {
              x[[j]]$response
            } else formula_respname(x[[j]]$cat.formula)
          } else {
            if(is.null(x[[j]]$response)) response.name else x[[j]]$response
          }, sep = '')
        et <- x[[j]]$pterms
        fctr <- attr(x[[j]]$formula, "control")
        if(is.null(fctr)) fctr <- ""
        if(length(i <- grep("Intercept", et, fixed = TRUE)))
          et[i] <- "const"
        if(length(x[[j]]$sx.smooth)) {
          et <- c(et, x[[j]]$sx.smooth)
        }
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
        eqtj <- if(family$cat & (if(is.null(id)) j else id) != 1) 3 else 2
        if(!is.null(x[[j]]$hlevel)) {
          if(!any(grepl("hlevel", fctr))) {
            teqn <- paste(teqn, " hlevel=", x[[j]]$hlevel, sep = "")
            if(x[[j]]$hlevel > 1) {
              if(!any(grepl("family", fctr))) {
                teqn <- paste(teqn, " family=", "gaussian_re", sep = "")
                teqn <- paste(teqn, " equationtype=",
                  family$bayesx[[nx[if(is.null(id)) j else id]]][eqtj], sep = "")
                ok <- FALSE
              }
            }
            if(x[[j]]$hlevel < 2 & ctr2) {
              if(any(grepl("predict", names(control$prg)))) {
                teqn <- paste(teqn, " predict=", control$prg$predict, sep = "")
              }
            }
          } else ok <- grepl("1", fctr[grepl("hlevel", fctr)])
        }
        if(ok) {
          if(!any(grepl("family", fctr)))
            teqn <- paste(teqn, " family=", family$bayesx[[nx[if(is.null(id)) j else id]]][1], sep = "")
          if(!any(grepl("equationtype", fctr)))
            teqn <- paste(teqn, " equationtype=", family$bayesx[[nx[if(is.null(id)) j else id]]][eqtj], sep = "")
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

  prg_extras <- function(x) {
    n <- length(x)
    prgex <- NULL
    for(j in n:1) {
      if(!"fake.formula" %in% names(x[[j]])) {
        prgex <- c(prgex, prg_extras(x[[j]]))
      } else {
        if(length(x[[j]]$sx.smooth)) {
          for(k in x[[j]]$sx.smooth) {
            if(!is.null(attr(k, "write"))) {
              if(!is.null(attr(k, "map.name"))) {
                if(!any(grepl(paste("map", attr(k, "map.name")), prgex)))
                  prgex <- c(prgex, attr(k, "write")(dir))
              } else prgex <- c(prgex, attr(k, "write")(dir))
            }
          }
        }
      }
    }
    if(!is.null(prgex)) prgex <- c(prgex, "")
    prgex
  }
  
  prg <- c(prg_extras(x), prg, make_eqn(x))

  prg <- gsub("(random", "(hrandom", prg, fixed = TRUE)
  for(i in 1:5)
    prg <- gsub(paste("(psplinerw", i, sep = ""), "(pspline", prg, fixed = TRUE)
  prg <- c(prg, "", paste(model.name, "getsample", sep = "."))
  prg <- c(
    paste('%% BayesX program created by BayesR: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )

  return(list("prg" = prg, "control" = control))
}


samplerBayesX <- function(x, ...)
{
  if(any(grepl("##seed##", x$prg, fixed = TRUE)))
    x$prg <- gsub("##seed##", round(runif(1L) * .Machine$integer.max), x$prg, fixed = TRUE)

  dir <- x$control$setup$dir
  if(!is.null(attr(dir, "unlink"))) {
    if(attr(dir, "unlink"))
      on.exit(unlink(dir))
  }
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
writeLines(x$prg)
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
    errm <- paste("an error occurred running the BayesX binary! The following messages are returned:\n",
      paste(errl, collapse = "\n", sep = ""), sep = "")
    warning(errm, call. = FALSE)
  }
  samples <- NULL
  mfile <- grep(paste(x$control$setup$model.name, "_R.r", sep = ""), dir(dir), value = TRUE, fixed = TRUE)
  if(length(mfile)) {
    mfile <- readLines(file.path(dir, mfile))
    mfile <- process_mfile(mfile)
    mspecs <- list("effects" = list(), "model" = list())
    for(j in 1:nrow(mfile$effects)) {
      ts <- read.table(mfile$effects[j, "samples"], header = TRUE)
      ts$intnr <- NULL
      names(ts) <- paste(mfile$effects[j, "terms"], "[", 1:ncol(ts), "]", sep = "")
      samples <- cbind(samples, as.matrix(ts))
      if(!is.na(mfile$effects[j, "varsamples"])) {
        ts <- read.table(mfile$effects[j, "varsamples"], header = TRUE)
        ts$intnr <- NULL
        names(ts) <- paste(mfile$effects[j, "terms"], ":var[", 1:ncol(ts), "]", sep = "")
        samples <- cbind(samples, as.matrix(ts))
      }
      mspecs$effects[[mfile$effects[j, "terms"]]] <- list(
        "basis" = if(!is.na(mfile$effects[j, "basis"])) {
            eval(parse(file = mfile$effects[j, "basis"]))
         } else function(x) { x },
        "family" = mfile$effects[j, "family"],
        "eqntype" = mfile$effects[j, "eqntype"],
        "filetype" = mfile$effects[j, "filetype"],
        "hlevel" =  as.integer(mfile$effects[j, "hlevel"])
      )
    }
    dic <- if(!is.null(mfile$model$dic)) read.table(mfile$model$dic, header = TRUE) else NULL
    samples <- cbind(samples, "dic" = dic$dic, "pd" = dic$pd)
    samples <- as.mcmc(samples)
    attr(samples, "model.specs") <- mspecs
  }

  samples
}

process_mfile <- function(x)
{
  predict <- dic <- NULL
  if(length(i <- grep("predict.res", x, ignore.case = TRUE))) {
    predict <- gsub(";", "", gsub("\\s", "", strsplit(x[i], "=", fixed = TRUE)[[1]][2]))
    x <- x[-i]
  }
  if(length(i <- grep("DIC.res", x, ignore.case = TRUE))) {
    dic <- gsub(";", "", gsub("\\s", "", strsplit(x[i], "=", fixed = TRUE)[[1]][2]))
    x <- x[-i]
  }
  n <- length(x)
  x <- gsub("equation=", "family=", x, fixed = TRUE)
  x <- gsub("-", "_", x, fixed = TRUE)
  samples <- varsamples <- eqntype <- family <- basis <- terms <- hlevel <- filetype <- rep(NA, length = n)

  c2c <- function(x) {
    if(any(grepl("(", x, fixed = TRUE))) {
      i <- regexpr("\\((.*)\\)", x)
      x2 <- substring(x, i + 1, i + attr(i, "match.length") - 2)
      x3 <- gsub(",", ";", x2, fixed = TRUE)
      x <- gsub(x2, x3, x, fixed = TRUE)
    }
    x
  }

  for(j in 1:n) {
    specs <- strsplit(c2c(x[j]), ",")[[1]]
    if(length(tmp <- grep("pathsamples", specs, value = TRUE))) {
      samples[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("pathvarsample", specs, value = TRUE))) {
      varsamples[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("equationtype", specs, value = TRUE))) {
      eqntype[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("filetype", specs, value = TRUE))) {
      filetype[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("family", specs, value = TRUE))) {
      family[j] <- tolower(gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2]))
    }
    if(length(tmp <- grep("pathbasis", specs, value = TRUE))) {
      basis[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("hlevel", specs, value = TRUE))) {
      hlevel[j] <- as.integer(gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2]))
    } else hlevel[j] <- 1
    if(length(tmp <- grep("term", specs, value = TRUE))) {
      tt <- strsplit(tmp, "=", fixed = TRUE)[[1]]
      tt <- paste(tt[2:length(tt)], collapse = "=")
      tt <- gsub(";", ",", tt, fixed = TRUE)
      tt <- gsub("^\\s+|\\s+$", "", tt)
      terms[j] <- gsub("\\s", "+", tt)
      
      ## FIXME: this is just a brute force fix for multinom hierarchical structures!
      ##        Searching for the correct category in the samples path.
      if(grepl("multinom", family[j])) {
        cat <- strsplit(samples[j], "_", fixed = TRUE)[[1]]
        cat <- cat[4]
        terms[j] <- paste(terms[j], cat, sep = ":")
      }

    }
  }
  for(j in 1:n) {
    if(grepl("gaussian_re", family[j], fixed = TRUE)) {
      i <- grep(eqntype[j], eqntype)
      i <- i[i != j]
      ok <- TRUE
      for(l in i) {
        if(hlevel[l] < 2 & ok) {
          family[j] <- family[l]
          ok <- FALSE
        }
      }
    }
  }

  ## FIXME: this is just a brute force fix for multinom hierarchical structures!
  ##        Searching for the correct category by matching the equationtype.
  ##        Not possible for more than 2 modeled categories!
  if(any(grepl("multinom", family))) {
    for(j in seq_along(hlevel)) {
      if(hlevel[j] > 1 & !grepl(":", terms[j], fixed = TRUE)) {
        et <- eqntype[j]
        ok <- TRUE
        for(i in seq_along(eqntype)) {
          if(i != j & et == eqntype[i]) {
            if(grepl(":", terms[i], fixed = TRUE) & ok) {
              terms[j] <- paste(terms[j], strsplit(terms[i], ":", fixed = TRUE)[[1]][2], sep = ":")
              ok <- FALSE
            }
          }
        }
      }
    }
  }

  known_paramaters <- c("mu", "sigma", "sigma2", "pi", "delta", "tau",
    "a", "b", "p", "df", "rho", "lambda", "alpha", "nu")

  ft <- sapply(strsplit(family, "_"), function(x) {
    if(length(x) > 1) {
      if(x[2] %in% known_paramaters) x[2] else x[1]
    } else x
  })
  terms <- paste(terms, ft, sep = ":")
  if(any(i <- duplicated(terms))) {
    for(j in terms[i]) {
      tt <- grep(j, terms, value = TRUE, fixed = TRUE)
      tt <- paste(tt, 1:length(tt), sep = ":")
      terms[terms == j] <- tt
    }
  }
  terms <- paste(terms, hlevel, sep = ":")

  rval <- list("effects" = cbind(terms, filetype, samples, varsamples, basis, family, eqntype, hlevel),
    "model" = list("predict" = predict, "dic" = dic))
  rval
}


resultsBayesX <- function(x, samples, ...)
{
  family <- attr(x, "family")
  grid <- attr(x, "grid")
  response.name <- attr(attr(x, "model.frame"), "response.name")
  if(is.null(grid)) grid <- 100
  if(is.function(family))
    family <- family()
  mspecs <- attr(samples, "model.specs")
  ylevels <- attr(x, "ylevels")
  reference <- attr(x, "reference")
  if(family$cat)
    ylevels <- ylevels[ylevels != reference]

  rename.p <- function(x) {
    if(family$family %in% c("beta", "gaussian", "lognormal", "binomial", "multinomial")) {
      foo <- switch(family$family,
        "gaussian" = function(x) gsub("sigma2", "sigma", x),
        "beta" = function(x) gsub("sigma2", "sigma", x),
        "lognormal" = function(x) gsub("sigma2", "sigma", x),
        "binomial" = function(x) gsub("binomial", "pi", x),
        "multinomial" = function(x) {
          if(any(grepl("):", x, fixed = TRUE))) {
            x <- strsplit(x, "):", fixed = TRUE)
            x <- sapply(x, function(x2) {
              if(length(x2) > 1) {
                x2[2] <- gsub(response.name, "", gsub("multinom:", "", x2[2]))
                x2[2] <- gsub("multinomialprobit:", "", x2[2])
                x2 <- paste(x2, collapse = "):")
              } else {
                x2 <- gsub(response.name, "", gsub("multinom:", "", x2))
                x2 <- gsub("multinomialprobit:", "", x2)
              }
              x2
            })
            x <- unlist(x)
          } else {
            x <- gsub(response.name, "", gsub("multinom:", "", x))
            x <- gsub("multinomialprobit:", "", x)
          }
          x
        }
      )
      x <- foo(x)
    }
    x
  }

  if(is.null(mspecs) & is.list(samples))
    mspecs <- attr(samples[[1]], "model.specs")
  names(mspecs$effects) <- rename.p(names(mspecs$effects))

  createBayesXresults <- function(obj, samples, id = NULL, sid = FALSE) {
    if(inherits(samples[[1]], "mcmc.list"))
      samples <- do.call("c", samples)
    chains <- length(samples)
    rval <- vector(mode = "list", length = chains)

    snames <- rename.p(colnames(samples[[1]]))
    if(is.null(snames))
      snames <- attributes(samples[[1]])$dimnames[[2]]
    for(j in 1:chains) {
      if(any(grepl("dic", snames))) {
        DIC <- mean(as.numeric(samples[[j]][, grepl("dic", snames)]))
        pd <- mean(as.numeric(samples[[j]][, grepl("pd", snames)]))
      } else {
        DIC <- pd <- NA
      }

      ## Compute model term effects.
      param.effects <- effects <- effects.hyp <- NULL
      fitted.values <- 0

      ## Parametric effects.
      if(k <- ncol(obj$X)) {
        nx <- obj$pterms
        nx <- gsub("Intercept", "const", nx, fixed = TRUE)
        pt <- paste(nx, collapse = "+")
        pt <- paste(pt, id, obj$hlevel, sep = ":")
        if(any(grepl(pt, snames, fixed = TRUE))) {
          samps <- as.matrix(samples[[j]][, grepl(pt, snames, fixed = TRUE)], ncol = k)
          nx <- gsub("const", "(Intercept)", nx, fixed = TRUE)
          colnames(samps) <- nx
          qu <- t(apply(samps, 2, quantile, probs = c(0.025, 0.5, 0.975)))
          sd <- drop(apply(samps, 2, sd))
          me <- drop(apply(samps, 2, mean))
          param.effects <- cbind(me, sd, qu)
          rownames(param.effects) <- nx
          colnames(param.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
          fitted.values <- as.vector(fitted.values + obj$X %*% param.effects[, 1])
          attr(param.effects, "samples") <- as.mcmc(samps)
        } else warning("please check the BayesX files!")
      }

      ## Smooth terms.
      if(length(i <- grep("sx(", names(mspecs$effects), fixed = TRUE))) {
        sx.smooth <- mspecs$effects[i]
        if(length(i <- grep(paste(id, obj$hlevel, sep = ":"), names(sx.smooth), fixed = TRUE))) {
          if(!is.list(effects))
            effects <- list()
          sx.smooth <- sx.smooth[i]
          for(i in names(sx.smooth)) {
            tn <- grep(i, snames, fixed = TRUE, value = TRUE)
            psamples <- as.matrix(samples[[j]][, snames %in% tn], ncol = length(tn))

            ## Possible variance parameter samples.
            if(length(vs <- grep(":var[", colnames(psamples), fixed = TRUE))) {
              vsamples <- psamples[, vs[1]]
              psamples <- psamples[, -vs, drop = FALSE]
            }

            ## Prediction matrix.
            tn0 <- strsplit(i, ":")[[1]][1]
            tn <- gsub(")", "", gsub("sx(", "", tn0, fixed = TRUE), fixed = TRUE)
            tn <- strsplit(tn, ",", fixed = TRUE)[[1]]
            tn <- tn[!grepl("by", tn)]
            tn1 <- eval(parse(text = tn0))

            basis <- sx.smooth[[i]]$basis

            if(tn1$by != "NA") {
              tn <- c(tn, tn1$by)
              get.X <- function(data) {
                diag(data[, ncol(data)]) %*% basis(data[, 1:(ncol(data) - 1)])
              }
            } else {
              get.X <- basis
            }

            X <- sx.smooth[[i]]$basis(attr(x, "model.frame")[1, tn, drop = FALSE])
            
            get.mu <- function(X, g) {
              X %*% as.numeric(g)
            }

            ## Compute final smooth term object.
            stype <- attr(X, "type")
            class(tn1) <- paste(stype, "smooth.spec", sep = ".")

            fst <- compute_term(tn1, get.X = get.X, get.mu = get.mu,
              psamples = psamples, vsamples = vsamples, asamples = NULL, FUN = NULL,
              snames = snames, effects.hyp = effects.hyp, fitted.values = fitted.values,
              data = attr(x, "model.frame")[, tn, drop = FALSE], grid = grid,
              hlevel = obj$hlevel)

            attr(fst$term, "specs")$get.mu <- get.mu
            attr(fst$term, "specs")$basis <- get.X ## sx.smooth[[i]]$basis
            if(sid)
              attr(fst$term, "specs")$label <- paste(attr(fst$term, "specs")$label, id, sep = ":")

            ## Add term to effects list.
            effects[[paste(tn0, stype, sep = ":")]] <- fst$term
            effects.hyp <- fst$effects.hyp
            fitted.values <- fst$fitted.values
            rm(fst)
          }
        }
      }

      ## Scale parameters.
      scale.m <- scale.samps.m <- NULL

      ## Compute partial residuals.
      if(!is.null(effects)) {
        if(length(obj$response)) {
          if(obj$response %in% names(attr(x, "model.frame"))) {
            if(!is.na(link <- family$links[id])) {
              link <- make.link2(link)
            } else link <- NULL
            effects <- partial.residuals(effects, attr(x, "model.frame")[[obj$response]],
              fitted.values, link)
          }
        }
      }

      ## Stuff everything together.
      rval[[j]] <- list(
        "model" = list("DIC" = DIC, "pd" = pd, "N" = nrow(attr(x, "model.frame")),
        "formula" = obj$formula), "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "scale" = scale.m, "fitted.values" = fitted.values,
        "residuals" = if(!is.factor(obj$response.vec)) obj$response.vec - fitted.values else NULL
      )

      class(rval[[j]]) <- "bayesr"
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    if(length(rval) < 2) {
      rval <- rval[[1]]
    }
    class(rval) <- "bayesr"
    return(rval)
  }

  nx <- names(x)
  nx <- nx[!(nx %in% c("call", "family"))]
  if(is(samples, "mcmc"))
    samples <- as.mcmc.list(list(samples))
  if(!all(c("formula", "fake.formula") %in% names(x))) {
    rval <- list()
    fn <- family$names
    if(is.null(names(x))) {
      nx <- paste("h", 1:length(x), sep = "")
      names(x) <- nx
    }
    if(length(unique(nx)) != length(nx)) {
      nx <- paste(rep(nx, length.out = length(nx)), 1:length(nx), sep = ":")
      names(x) <- nx
    }
    if(family$cat) {
      fn <- rmf(ylevels)
      names(x) <- nx <- fn
    }
    if(length(fn) != length(nx))
      fn <- nx
    for(j in seq_along(nx)) {
      if(!all(c("formula", "fake.formula") %in% names(x[[nx[j]]]))) {
        rval[[nx[j]]] <- list()
        for(ii in seq_along(x[[nx[j]]])) {
          nh <- paste("h", ii, sep = "")
          rval[[nx[j]]][[nh]] <- createBayesXresults(x[[nx[j]]][[ii]],
            samples, id = fn[j], sid = length(nx) > 1)
        }
        attr(rval[[nx[j]]], "hlevel") <- TRUE
      } else {
        rval[[nx[j]]] <- createBayesXresults(x[[nx[j]]], samples, id = fn[j], sid = length(nx) > 1)
      }
      class(rval[[nx[j]]]) <- "bayesr"
    }
    if(length(nx) > 1)
      names(rval) <- fn
    else
      rval <- rval[[1]]
    class(rval) <- "bayesr"
    return(rval)
  } else {
    return(createBayesXresults(x, samples, id = family$names))
  }
}


########################################
## (3) BayesX model term construction ##
########################################
sx <- function(x, z = NULL, bs = "ps", by = NA, ...)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  call <- match.call()

  available.terms <- c(
    "rw1", "rw2",
    "season",
    "ps", "psplinerw1", "psplinerw2", "pspline",
    "te", "pspline2dimrw2", "te1", "pspline2dimrw1",
    "kr", "kriging",
    "gk", "geokriging",
    "gs", "geospline",
    "mrf", "spatial",
    "bl", "baseline",
    "factor",    
    "ridge", "lasso", "nigmix",
    "re", "ra", "random",
    "cs", "catspecific",
    "offset",
    "generic",
    "rsps", "hrandom_pspline"
  )
  if(!bs %in% available.terms) stop(paste("basis type", sQuote(bs), "not supported by BayesX"))

  if(bs %in% c("rsps", "hrandom_pspline")) {
    bs <- "rsps"
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
    rcall <- paste("R2BayesX:::r(x = ", by, ", bs = ", sQuote(bs), ", by = ", x, ", ...)", sep = "")
    rval <- eval(parse(text = rcall))
  } else {
    if(length(grep("~", term <- deparse(call$x))) && bs %in% c("re", "ra", "random")) {
      x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
      rcall <- paste("R2BayesX:::r(x = ", x, ", by = ", by, ", ...)", sep = "")
      rval <- eval(parse(text = rcall))
    } else {
      k <- -1
      m <- NA
      xt <- list(...)
      if("m" %in% names(xt))
        stop("argument m is not allowed, please see function s() using this specification!")
      if("k" %in% names(xt))
        stop("argument k is not allowed, please see function s() using this specification!")
      if(!is.null(xt$xt))
        xt <- xt$xt
      warn <- getOption("warn")
      options("warn" = -1)
      if(by != "NA" && is.vector(by) && length(by) < 2L && !is.na(as.numeric(by))) {
        xt["b"] <- by
        by <- "NA"
      }
      options("warn" = warn)
      if(bs %in% c("pspline2dimrw1", "pspline2dimrw2", "te",
        "gs", "geospline", "kr", "gk", "kriging", "geokriging")) {
        if(by != "NA")
          stop(paste("by variables are not allowed for smooths of type bs = '", bs, "'!", sep = ""))
      }
      if(bs == "te1") bs <- "pspline2dimrw1"
      if(!is.null(xt$knots))
        xt$nrknots <- xt$knots
      if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
        "pspline2dimrw2", "pspline2dimrw1", "gs", "geospline")) {
        if(!is.null(xt$degree))
          m <- xt$degree
        if(!is.null(xt$order)) {
          if(is.na(m))
            m <- c(3L, xt$order)
          else
            m <- c(m[1L], xt$order)
        }
        if(length(m) < 2L && (bs %in% c("gs", "geospline")))
          m <- c(3L, 1L)
        if(length(m) < 2L && is.na(m))
          m <- c(3L, 2L)   
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 2L)
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 1L)
        m[1L] <- m[1L] - 1L
        if(!is.null(xt$nrknots))
          k <- xt$nrknots + m[1L]
        else {
          if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
            k <- 20L + m[1L]
          else
            k <- 10L + m[1L]
        }
      }
      if(bs %in% c("kr", "gk", "kriging", "geokriging")) {
        m <- c(1L, 1L)
        if(!is.null(xt$nrknots)) {
          k <- xt$nrknots
        } else k <- -1L
      }
      if(!is.null(xt$map))
        xt$map.name <- rmf(as.character(call$map))
      xt[c("degree", "order", "knots", "nrknots")] <- NULL
      if(!length(xt))
        xt <- NULL
      if(!is.null(call$z)) 
        term <- c(term, deparse(call$z))
      rval <- mgcv::s(x, z, k = k, bs = bs, m = m, xt = xt)
      rval$term <- term
      rval$dim <- length(term)
      rval$by <- by
      rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""), ")", sep = "")
    }
  }

  return(rval)
}

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

sx.construct.mrf.smooth.spec <- sx.construct.spatial.smooth.spec <- function(object, data)
{
  if(is.null(object$xt))
    stop("need to supply a map object in argument xt!")  
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  if(!is.null(object$xt$map.name))
    map.name <- object$xt$map.name
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- rmf(gsub("\\s", "", paste(map.name, sep = "", collapse = "")))

  map <- object$xt$map
  if(is.null(map)) {
    if(!is.null(object$xt$polys))
      map <- object$xt$polys
    if(!is.null(object$xt$penalty))
      map <- object$xt$penalty
  }
  if(is.null(map))
    map <- object$xt$gra
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]])) {
      if(inherits(object$xt[[1L]], "gra"))
        map <- object$xt[[1L]]
      else
        map <- object$xt
    } else map <- object$xt[[1L]]
    if(is.null(map)) {
      map <- object$xt
      if(is(map, "SpatialPolygonsDataFrame"))
        map <- SPDF2bnd(map)
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd") || !inherits(map, "gra")))
        stop("need to supply a bnd or graph file object in argument xt!")
    }
  }
  if(is(map, "nb"))
    map <- nb2gra(map)
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!inherits(map, "bnd") && !inherits(map, "gra")) {
    if(is.list(map))
      class(map) <- "bnd"
    else
      class(map) <- "gra"
  }

  write <- function(dir = NULL) {
    if(is.null(dir)) {
      dir.create(dir <- tempfile())
      on.exit(unlink(dir))
    } else dir <- path.expand(dir)
    if(!file.exists(dir)) dir.create(dir)

    counter <- NULL
    ok <- TRUE
    files <- list.files(dir)
    while(ok) {
      classm <- class(map)
      if(length(classm) > 1L)
        if("list" %in% classm)
          class(map) <- classm[classm != "list"]
      mapfile <- paste(map.name, counter, ".", class(map), sep = "")[1]
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else ok <- FALSE
    }
    mapfile <- file.path(dir, mapfile)

    prg <- paste("map", map.name)
    if(inherits(map, "bnd")) {
      if(!any(is.na(poly.names <- as.integer(names(map))))) {
        poly.names <- sort(poly.names)
        poly.names <- as.character(poly.names)
      } else poly.names <- sort(names(map))
      map <- map[poly.names]
      class(map) <- "bnd"
      write.bnd(map = map, file = mapfile, replace = TRUE)
      prg <- c(prg, paste(map.name, ".infile using ", mapfile, sep = ""))
    } else {
      if(!is.character(map)) {
        write.gra(map = map, file = mapfile, replace = TRUE)
        prg <- c(prg, paste(map.name, ".infile, graph using ", mapfile, sep = ""))
      } else {
        stopifnot(is.character(map))
        pos <- regexpr("\\.([[:alnum:]]+)$", map)
        fext <- ifelse(pos > -1L, substring(map, pos + 1L), "")
        if(fext == "gra")
          prg <- c(prg, paste(map.name, ".infile, graph using ", path.expand(map), sep = ""))
        else
          prg <- c(prg, paste(map.name, ".infile using ", path.expand(map), sep = ""))
      }
    }
    prg
  }

  term <- object$term
  term <- paste(term, "(spatial,map=", map.name, sep = "")
  term <- paste(do.xt(term, object, c("map", "polys", "penalty", "map.name")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  attr(term, "write") <- write
  attr(term, "map.name") <- map.name

  return(term)
}

make_by <- function(term, object, data)
{
  if(!missing(data) && !is.character(data)) {
    by <- data[[object$by]]
    if(is.factor(by) && nlevels(by) > 1) {
      nocenter <- paste(c("lasso", "nigmix", "ridge", "ra"), ".smooth.spec", sep = "")
      term <- paste(paste(rmf(object$by), rmf(levels(by)), sep = ""), "*", term, sep = "")
      if((k <- length(term)) > 1)
        for(j in 1:k) 
          if(!grepl("center", term[j]) && is.null(object$xt$center))
            if(!class(object) %in% nocenter)
              term[j] <- gsub(")", ",center)", term[j])
    } else term <- paste(rmf(object$by), "*", term, sep = "")
  } else term <- paste(rmf(object$by), "*", term, sep = "")

  return(term)
}

rmf <- function(x) 
{
  for(i in 1L:length(x)) {
    for(char in c("+", "-", "*", ":", "^", "/", " ", "(", ")", "]", "[",
      ",", ".", "<", ">", "?", "!", "'", "#", "~", "`", ";", "=", "&", "$", "@")) {
      x[i] <- gsub(char, "", x[i], fixed = TRUE)
    }
  }

  return(x)
}

help.map.name <- function(x)
{
  if(is.null(x))
    return("")
  x <- splitme(x)
  if(resplit(x[1L:2L]) == "s(") {
    x <- resplit(c("sfoofun", x[2L:length(x)]))
    x <- eval(parse(text = x), envir = parent.frame())
    if(is.null(x))
      x <- "map"
  } else  x <- "map"

  return(x)
}

sfoofun <- function(x, xt = NULL, ...)
{
  if(is.null(x) || is.null(xt))
    return(NULL)
  x <- rval <- deparse(substitute(xt), backtick = TRUE, width.cutoff = 500L)
  if(is.list(xt) && length(xt)>1)
    for(i in 1L:length(xt))
      if(inherits(xt[[i]], "bnd") || inherits(xt[[i]], "gra") || inherits(xt[[i]], "list")) {
        rval <- strsplit(x, ",", " ")[[1L]]
        if(length(rval) > 1L)
          rval <- rval[i]
      }
  rval <- splitme(rval)
  if(length(rval) > 5L)
    if(resplit(rval[1L:5L]) == "list(")
      rval <- rval[6L:length(rval)]
  if(rval[length(rval)] == ")")
    rval <- rval[1L:(length(rval) - 1)]
  if(any(grepl("=", rval)))
    rval <- rval[(grep("=", rval) + 2L):length(rval)]
  rval <- resplit(rval)
   
  return(rval)
}

splitme <- function(x) {
  return(strsplit(x, "")[[1L]])
}

resplit <- function(x) {
  if(!is.null(x))
    x <- paste(x, sep = "", collapse = "")
	
  return(x)
}


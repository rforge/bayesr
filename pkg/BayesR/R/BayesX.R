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

  rval <- bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transformBayesX,
    setup = setup, sampler = samplerBayesX, results = resultsBayesX,
    cores = cores, combine = combine, ...)
  
  attr(rval, "call") <- match.call()
  
  rval
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
  if(!is.null(family$order))
    x <- x[rev(family$order)]
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
        tmp <- sx.construct(x$sx.smooth[[j]], x$mf)
        attr(tmp, "specs") <- x$sx.smooth[[j]]
        x$sx.smooth[[j]] <- tmp
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
  family <- x$family
  x$call <- x$family <- NULL

  args <- list(...)
  model.name <- control$setup$model.name
  data.name <- control$setup$data.name
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  prg <- c(paste('logopen using', file.path(dir, paste(model.name, 'log', sep = '.'))), "")

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

  if(!is.null(family$all)) {
    if(n < family$k & family$all)
      stop('all parameters must have a model formula, e.g. at least an intercept only model "~ 1"!')
  }

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

  make_eqn <- function(x, ctr = TRUE, id = NULL) {
    n <- length(x)
    eqn <- NULL
    for(j in n:1) {
      if(!"fake.formula" %in% names(x[[j]])) {
        eqn <- c(eqn, make_eqn(x[[j]], ctr, id = j))
        ctr <- FALSE
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
            if(x[[j]]$hlevel < 2) {
              if(any(grepl("predict", names(control$prg)))) {
                teqn <- paste(teqn, " predict=", control$prg$predict, sep = "")
              }
            }
          } else ok <- grepl("1", fctr[grepl("hlevel", fctr)])
        }
        if(ok) {
          if(!any(grepl("family", fctr)))
            teqn <- paste(teqn, " family=", family[[nx[if(is.null(id)) j else id]]][1], sep = "")
          if(!any(grepl("equationtype", fctr)))
            teqn <- paste(teqn, " equationtype=", family[[nx[if(is.null(id)) j else id]]][2], sep = "")
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
    stop(errm, call. = FALSE)
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
      mspecs$effects[[mfile$effects[j, "terms"]]] <- list(
        "basis" = if(!is.na(mfile$effects[j, "basis"])) {
            eval(parse(file = mfile$effects[j, "basis"]))
         } else function(x) { x },
        "family" = mfile$effects[j, "family"],
        "eqntype" = mfile$effects[j, "eqntype"],
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
  samples <- eqntype <- family <- basis <- terms <- hlevel <- rep(NA, length = n)
  for(j in 1:n) {
    specs <- strsplit(x[j], ",")[[1]]
    if(length(tmp <- grep("pathsamples", specs, value = TRUE))) {
      samples[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("equationtype", specs, value = TRUE))) {
      eqntype[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("family", specs, value = TRUE))) {
      family[j] <- tolower(gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2]))
    }
    if(length(tmp <- grep("pathbasis", specs, value = TRUE))) {
      basis[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    }
    if(length(tmp <- grep("hlevel", specs, value = TRUE))) {
      hlevel[j] <- gsub("\\s", "", strsplit(tmp, "=", fixed = TRUE)[[1]][2])
    } else hlevel[j] <- 1
    if(length(tmp <- grep("term", specs, value = TRUE))) {
      tt <- strsplit(tmp, "=", fixed = TRUE)[[1]][2]
      tt <- gsub("^\\s+|\\s+$", "", tt)
      terms[j] <- gsub("\\s", "+", tt)
      terms[j] <- paste(terms[j], paste(if(any(grepl("_", family[j]))) {
        strsplit(family[j], "_", fixed = TRUE)[[1]][2]
      } else tolower(family[j]), hlevel[j], sep = ":"), sep = ":")
    }
  }
  if(any(i <- duplicated(terms))) {
    for(j in terms[i]) {
      tt <- grep(j, terms, value = TRUE, fixed = TRUE)
      tt <- paste(tt, 1:length(tt), sep = ":")
      terms[terms == j] <- tt
    }
  }
  rval <- list("effects" = cbind(terms, samples, basis, family, eqntype, hlevel),
    "model" = list("predict" = predict, "dic" = dic))
  rval
}

resultsBayesX <- function(x, samples, id = "", mspecs = NULL, ...)
{
  nx <- names(x)
  nx <- nx[!(nx %in% c("call", "family"))]
  if(length(nx) < 2)
    x <- x[[nx]]
  if(is.null(mspecs)) {
    mspecs <- attr(samples, "model.specs")
    if(is.null(mspecs) & is.list(samples))
      mspecs <- attr(samples[[1]], "model.specs")
  }
  if(is(samples, "mcmc"))
    samples <- as.mcmc.list(list(samples))
  if(!all(c("family", "formula") %in% names(x))) {
    rval <- list()
    family <- x[[1]]$family
    if(is.function(family))
      family <- family()
    fn <- family$names
    if(is.null(names(x))) {
      nx <- paste("h", 1:length(x), sep = "")
      names(x) <- nx
    }
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- resultsBayesX(x[[nx[j]]], samples, id = fn[j], mspecs = mspecs)
      if(FALSE) {
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
    if(inherits(samples[[1]], "mcmc.list"))
      samples <- do.call("c", samples)
    chains <- length(samples)
    rval <- vector(mode = "list", length = chains)
    snames <- colnames(samples[[1]])
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
      if(k <- ncol(x$X)) {
        nx <- x$pterms
        if(attr(terms(x$fake.formula), "intercept"))
          nx <- c("const", nx)
        pt <- paste(nx, collapse = "+")
        pt <- paste(pt, id, sep = ":")
        samps <- as.matrix(samples[[j]][, grepl(pt, snames, fixed = TRUE)], ncol = k)
        nx <- gsub("const", "(Intercept)", nx, fixed = TRUE)
        qu <- t(apply(samps, 2, quantile, probs = c(0.025, 0.5, 0.975)))
        sd <- drop(apply(samps, 2, sd))
        me <- drop(apply(samps, 2, mean))
        param.effects <- cbind(me, sd, qu)
        rownames(param.effects) <- nx
        colnames(param.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
        fitted.values <- as.vector(fitted.values + x$X %*% param.effects[, 1])
        attr(param.effects, "samples") <- as.mcmc(samps)
        colnames(attr(param.effects, "samples")) <- nx
      }

      ## Smooth terms.
      if(length(i <- grep("sx(", names(mspecs$effects), fixed = TRUE))) {
        sx.smooth <- mspecs$effects[i]
        if(length(i <- grep(id, names(sx.smooth)))) {
          if(!is.list(effects))
            effects <- list()
          sx.smooth <- sx.smooth[i]
          for(i in names(sx.smooth)) {
            tn <- grep(i, snames, fixed = TRUE, value = TRUE)
            psamples <- as.matrix(samples[[j]][, snames %in% tn], ncol = length(tn))

            ## Prediction matrix.
            tn0 <- strsplit(i, ":")[[1]][1]
            tn <- gsub(")", "", gsub("sx(", "", tn0, fixed = TRUE), fixed = TRUE)
            tn <- strsplit(tn, ",", fixed = TRUE)[[1]]
            X <- sx.smooth[[i]]$basis(x$mf[, tn])

            ## Possible variance parameter samples.
            vsamples <- NULL

            get.mu <- function(X, g) {
              X %*% as.numeric(g)
            }

            ## Compute samples of fitted values.
            fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })

            ## Compute final smooth term object.
            fst <- compute_term(eval(parse(text = tn0)), fsamples = fsamples, psamples = psamples,
              vsamples = vsamples, FUN = NULL, snames = snames,
              effects.hyp = effects.hyp, fitted.values = fitted.values, data = x$mf)

            attr(fst$term, "specs")$get.mu <- get.mu
            attr(fst$term, "specs")$basis <- sx.smooth[[i]]$basis

            ## Add term to effects list.
            effects[[tn0]] <- fst$term
            effects.hyp <- fst$effects.hyp

            fitted.values <- fst$fitted.values
            rm(fst)
          }
        }
      }

      ## Scale parameters.
      scale.m <- scale.samps.m <- NULL

      ## Compute partial residuals. FIXME: binomial()$linkfun()
      if(x$response %in% names(x$mf)) {
        stats <- if(id != "") {
          make.link(x$family[[paste(id, "link", sep = ".")]])
        } else {
          make.link(x$family[[grep(".link", names(x$family), fixed = TRUE)]])
        }
        for(i in seq_along(effects)) {
          e <- stats$linkfun(x$mf[[x$response]]) - (fitted.values - attr(effects[[i]], "fit"))
          if(is.null(attr(effects[[i]], "specs")$xt$center)) {
            e <- e - mean(e)
          } else {
            if(attr(effects[[i]], "specs")$xt$center)
              e <- e - mean(e)
          }
          e <- cbind(attr(effects[[i]], "x"), e)
          if(!is.null(attr(effects[[i]], "by.drop")))
            e <- e[attr(effects[[i]], "by.drop"), ]
          e <- as.data.frame(e)
          try(names(e) <- c(attr(effects[[i]], "specs")$term, "partial.resids"))
          attr(effects[[i]], "partial.resids") <- e
          attr(effects[[i]], "fit") <- NULL
          attr(effects[[i]], "x") <- NULL
          attr(effects[[i]], "by.drop") <- NULL
        }
      }

      ## Stuff everything together.
      rval[[j]] <- list("call" = x$call, "family" = x$family,
        "model" = list("DIC" = DIC, "pd" = pd, "N" = nrow(x$mf),
        "formula" = x$formula), "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "scale" = scale.m, "fitted.values" = fitted.values,
        "residuals" = x$mf[[x$response]] - fitted.values)

      class(rval[[j]]) <- "bayesr"
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    if(length(rval) < 2) {
      rval <- rval[[1]]
    }
    class(rval) <- "bayesr"
    return(rval)
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
        xt$map.name <- as.character(call$map)
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


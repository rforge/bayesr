####################################
## (1) BayesX specific functions. ##
####################################
transformBayesX <- function(x, ...)
{
  family <- bamlss.family(attr(x, "family"))
  call <- x$call; x$call <- NULL

  if(family$cat) {
    ylevels <- attr(x, "ylevels")
    reference <- attr(x, "reference")
    grid <- attr(x, "grid")
    rn <- attr(attr(x, "model.frame"), "response.name")
    if(is.null(ylevels)) {
      ylevels <- rn
      reference <- "NA"
      names(x) <- paste(names(x)[1], ylevels, sep = ":")
      for(j in seq_along(x))
        x[[j]]$cat.formula <- x[[j]]$formula
    } else {
      f <- as.formula(paste("~ -1 +", rn))
      rm <- as.data.frame(model.matrix(f, data = attr(x, "model.frame")))
      cn <- colnames(rm) <- rmf(colnames(rm))
      rm <- rm[, !grepl(reference, cn), drop = FALSE]
      mf <- cbind(attr(x, "model.frame"), rm)
      attr(mf, "response.name") <- rn
      x <- x[ylevels != reference]
      attr(x, "model.frame") <- mf
      attr(x, "grid") <- grid
    }
    attr(x, "ylevels") <- ylevels
    attr(x, "reference") <- reference
  }

  attr(attr(x, "model.frame"), "orig.names") <- names(attr(x, "model.frame"))
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

      if(length(obj$smooth)) stop("arbitrary smooths not supported yet, please use sx()!")
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

  if(inherits(x, "bamlss.input") & !any(c("formula", "fake.formula", "response") %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  names(x) <- rep(family$names, length.out = n)
  if(!is.null(family$bayesx$order)) {
    mf <- attr(x, "model.frame")
    class(x) <- "list"
    x <- x[rev(family$bayesx$order)]
    class(x) <- c("bamlss.input", "list")
    attr(x, "model.frame") <- mf; rm(mf)
  }
  attr(x, "call") <- call
  attr(x, "family") <- family

  x
}


controlBayesX <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "light", model.name = "bamlss", data.name = "d",
  prg.name = NULL, dir = NULL, verbose = FALSE, cores = NULL, ...)
{
  if(is.null(seed))
    seed <- '##seed##'
  stopifnot(burnin < n.iter)
  if(is.null(model.name))
    model.name <- 'bamlss'
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
      "main" = c(rep(FALSE, 3), rep(TRUE, 2)), "model.name" = model.name, "data.name" = data.name,
      "prg.name" = prg.name, "dir" = dir, "verbose" = verbose, "cores" = cores
    )
  )

  cvals
}

setupBayesX <- function(x, control = controlBayesX(...), ...)
{
  names(attr(x, "model.frame")) <- rmf(names(attr(x, "model.frame")))

  mf.names <- names(attr(x, "model.frame"))

  family <- attr(x, "family")
  x$call <- x$family <- NULL
  lhs <- family$bayesx$lhs

  ## Handling of weights.
  add.weights <- FALSE
  if(!is.null(family$bayesx$weights)) {
    weights0 <- attr(x, "model.frame")$weights
    rn <- attributes(attr(x, "model.frame"))$response.name
    for(j in names(family$bayesx$weights)) {
      weights1 <- family$bayesx$weights[[j]](attr(x, "model.frame")[[rn[1]]])
      if(!is.null(weights0)) weights1 <- weights1 * weights0
      attr(x, "model.frame")[[paste(j, "weights", sep = "")]] <- weights1
      rm(weights1)
    }
	  rm(weights0)
    mf.names <- names(attr(x, "model.frame"))
    add.weights <- TRUE
  }
  if(length(grep("weights", mf.names)))
    add.weights <- TRUE
  zero <- family$bayesx$zero
  if(is.null(zero)) zero <- FALSE
  if(zero) {
    attr(x, "model.frame")[["ybinom"]] <- 1 * (attr(x, "model.frame")[[rn]] > 0)
  }
  quantile <- family$bayesx$quantile

  dirichlet <- family$family == "dirichlet"
  family$bayesx[c("order", "lhs", "rm.number", "weights", "quantile", "zero")] <- NULL

  args <- list(...)
  model.name <- control$setup$model.name
  data.name <- rmf(control$setup$data.name)
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  BayesX_data <- function(obj, h = FALSE, id = NULL) {
    if(!is.null(obj$cat.formula) & !h) {
      obj$response <- formula_respname(obj$cat.formula)
      obj$response.vec <- attr(x, "model.frame")[[obj$response]]
    }
    if(zero & !is.null(attr(obj$formula, "name"))) {
      if(attr(obj$formula, "name") == "pi" & !h) {
        obj$response <- "ybinom"
        obj$response.vec <- attr(x, "model.frame")[["ybinom"]]
      }
    }
    if(h)
      obj$response.vec <- attr(x, "model.frame")[[obj$response]]
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
      for(sj in obj$sterms) {
        if(is.factor(X[, sj])) {
          f <- as.formula(paste("~ -1 +", sj))
          fX <- as.data.frame(model.matrix(f, data = X))
          X <- cbind(X, fX)
          X <- X[, unique(colnames(X))]
        }
      }
    }
    if(!is.null(X)) {
      cnX <- colnames(X)
      if(ncol(X) > 0)
        X <- X[, unique(cnX[!grepl("Intercept", cnX, fixed = TRUE)]), drop = FALSE]
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
    if(!h) {
      wi <- paste(id, "weights", sep = "")
      if(!is.null(id) & (wi %in% mf.names)) {
        X[[wi]] <- attr(x, "model.frame")[[wi]]
      } else {
        if("weights" %in% mf.names)
          X$ModelWeights <- attr(x, "model.frame")[["weights"]]
      }
      if("offset" %in% mf.names)
        X$ModelOffset <- attr(x, "model.frame")[["offset"]]
    }
    if(h) {
      X <- unique(X)
      if(!is.null(obj$response))
        X <- X[order(X[[obj$response]]), , drop = FALSE]
    }

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
          d2 <- BayesX_data(x[[j]][[i]], id = nx[j])
          if(is.null(d)) {
            d <- d2
          } else {
            if(!is.null(d2))
              d <- cbind(d, d2)
          }
          d2 <- NULL
          d <- as.data.frame(d)
          d <- d[, unique(names(d)), drop = FALSE]
          x[[j]][[i]]$dname <- dname0
        } else d2 <- BayesX_data(x[[j]][[i]], i > 1)
        if(is.null(x[[j]][[i]]$hlevel))
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
      d2 <- BayesX_data(x[[j]], id = nx[j])
      if(is.null(d)) {
        d <- d2
      } else {
        if(!is.null(d2))
          d <- cbind(d, BayesX_data(x[[j]], id = nx[j]))
      }
      d2 <- NULL
      if(!is.null(d)) {
        d <- as.data.frame(d)
        d <- d[, unique(names(d)), drop = FALSE]
      }
      x[[j]]$dname <- dname0
      x[[j]]$hlevel <- 1
    }
  }
  response.name <- attr(attr(x, "model.frame"), "response.name")
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

  make_eqn <- function(x, ctr = TRUE, id = NULL) {
    n <- length(x)
    eqn <- NULL
    for(j in n:1) {
      if(!"fake.formula" %in% names(x[[j]])) {
        eqn <- c(eqn, make_eqn(x[[j]], ctr, id = j))
        ctr <- FALSE
      } else {
        teqn <- paste(model.name, '.hregress ',
          if(family$cat | !is.null(lhs)) {
            if(x[[j]]$hlevel > 1) {
              x[[j]]$response
            } else {
              if(!is.null(lhs)) {
                if(nx[if(is.null(id)) j else id] %in% names(lhs)) {
                  lhs[nx[if(is.null(id)) j else id]]
                } else {
                  if(is.null(x[[j]]$response)) response.name[1] else x[[j]]$response
                }
              } else formula_respname(x[[j]]$cat.formula)
            }
          } else {
            if(zero & x[[j]]$hlevel < 2 & !is.null(x[[j]]$formula)) {
              if(attr(x[[j]]$formula, "name") == "pi") {
                "ybinom"
              } else {
                if(is.null(x[[j]]$response)) response.name[1] else x[[j]]$response
              }
            } else {
              if(is.null(x[[j]]$response)) response.name[1] else x[[j]]$response
            }
          }, sep = '')
        et <- x[[j]]$pterms
        fctr <- attr(x[[j]]$formula, "control")
        if(is.null(fctr)) fctr <- ""
        if(length(i <- grep("Intercept", et, fixed = TRUE)))
          et[i] <- "const"
        if(length(x[[j]]$sx.smooth)) {
          et <- c(et, unlist(x[[j]]$sx.smooth))
        }
        if(x[[j]]$hlevel < 2) {
          if("offset" %in% mf.names)
            et <- c(et, "ModelOffset(offset)")
        }
        if(length(et))
          teqn <- paste(teqn, '=', paste(et, collapse = ' + '))
        if(x[[j]]$hlevel < 2 & add.weights) {
          wi <- paste(nx[if(is.null(id)) j else id], "weights", sep = "")
          if(wi %in% mf.names) {
            teqn <- paste(teqn, "weight", wi)
          } else {
            if("weights" %in% mf.names)
              teqn <- paste(teqn, "weight ModelWeights")
          }
        }
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
            teqn <- paste(teqn, " hlevel=", if(x[[j]]$hlevel > 1) 2 else 1, sep = "")
            if(x[[j]]$hlevel > 1) {
              if(!any(grepl("family", fctr)))
                teqn <- paste(teqn, " family=gaussian_re", sep = "")
              if(!any(grepl("family", fctr))) {
                teqn <- paste(teqn, " equationtype=", family$bayesx[[nx[if(is.null(id)) j else id]]][eqtj], sep = "")
                ok <- FALSE
              }
            }
            if(x[[j]]$hlevel < 2 & (if(is.null(id)) j else id)  < 2) {
              if(any(grepl("predict", names(control$prg))) & !zero) {
                teqn <- paste(teqn, " predict=", control$prg$predict, " setseed=", control$prg$setseed, sep = "")
              }
            }
          } else ok <- grepl("1", fctr[grepl("hlevel", fctr)])
        }
        if(ok) {
          if(dirichlet) {
            teqn <- paste(teqn, " nrcat=", family$ncat, " family=dirichlet equationtype=",
              if(j < 2) "mu" else paste("alpha", j, sep = ""), sep = "")
          } else {
            if(!any(grepl("family", fctr)))
              teqn <- paste(teqn, " family=", family$bayesx[[nx[if(is.null(id)) j else id]]][1], sep = "")
            if(!any(grepl("equationtype", fctr)))
					    teqn <- paste(teqn, " equationtype=", family$bayesx[[nx[if(is.null(id)) j else id]]][eqtj], sep = "")
            if(!is.null(quantile))
              teqn <- paste(teqn, " quantile=", quantile, sep = "")
          }
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

  if(zero) {
    prg <- c(prg, "", paste(model.name, ".hregress ybinom, family=zero_adjusted setseed=",
      control$prg$setseed," predict=light using ", x[[1]]$dname, sep = "" ))
  }

  prg <- gsub("(random", "(hrandom", prg, fixed = TRUE)
  for(i in 1:5)
    prg <- gsub(paste("(psplinerw", i, sep = ""), "(pspline", prg, fixed = TRUE)
  prg <- c(prg, "", paste(model.name, "getsample", sep = "."))
  prg <- c(
    paste('%% BayesX program created by bamlss: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )

  return(list("prg" = prg, "control" = control))
}


BayesX <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  data = NULL, control = controlBayesX(...), ...)
{
  if(is.null(family$bayesx))
    stop("BayesX specifications missing in family object, cannot set up model!")

  fbx <- family$bayesx

  model.name <- control$setup$model.name
  data.name <- rmf(control$setup$data.name)
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  if(!file.exists(dir)) {
    dir.create(dir)
    on.exit(unlink(dir))
  }

  if(is.null(data)) {
    data <- try(get('bf', envir = parent.frame())$model.frame, silent = TRUE)
    if(inherits(data, "try-error"))
      stop("cannot find the model.frame for creating BayesX model term objects!")
  }
  if(is.null(data))
    stop("no data available for creating BayesX model term objects!")

  yname <- colnames(y)[1]

  single_eqn <- function(x, y, id) {
    rhs <- dfiles <- prgex <- sdata <- NULL

    if(!is.null(x$model.matrix)) {
      cn <- rmf(colnames(x$model.matrix))
      colnames(x$model.matrix) <- cn
      cn <- cn[cn != "Intercept"]
      if(length(cn)) {
        rhs <- c(rhs, cn)
        sdata <- as.data.frame(x$model.matrix[, cn, drop = FALSE])
      }
      if("Intercept" %in% colnames(x$model.matrix)) {
        sdata <- cbind("Intercept" = rep(1, nrow(y)), sdata)
        rhs <- c("const", rhs)
      }
    }

    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct)) {
        if(is.null(x$smooth.construct[[j]]$sx.construct))
          class(x$smooth.construct[[j]]) <- "userdefined.smooth.spec"
        sxc <- sx.construct(x$smooth.construct[[j]], data, id = c(id, j), dir = dir)
        if(!is.null(attr(sxc, "write")))
          prgex <- c(prgex, attr(sxc, "write")(dir))
        rhs <- c(rhs, sxc)
        sdata <- if(is.null(sdata)) {
          data[, x$smooth.construct[[j]]$term, drop = FALSE]
        } else cbind(sdata, data[, x$smooth.construct[[j]]$term, drop = FALSE])
      }
    }

    rn <- response.name(as.formula(x$formula), hierarchical = FALSE)
    if(rn %in% family$names)
      rn <- NA
    if(is.na(rn))
      rn <- yname
    eqn <- paste(rn, "=", paste(rhs, collapse = " + "))
    rval <- list("eqn" = eqn, "prgex" = prgex)

    if(!is.null(sdata)) {
      if(nrow(sdata) == nrow(y))
        sdata <- cbind(sdata, y)
      rval$dname <- paste(paste(id, collapse = "_"), data.name, sep = "_")
      write.table(sdata, file = file.path(dir, paste(rval$dname, ".raw", sep = "")),
        quote = FALSE, row.names = FALSE)
      rval$prgex <- c(
        paste("dataset", rval$dname),
        paste(rval$dname, ".infile using ",
          file.path(dir, paste(rval$dname, ".raw", sep = "")), sep = ""),
        rval$prgex
      )
    }

    rval
  }

  eqn <- list()
  prgex <- NULL
  n <- 1
  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
      eqn[[i]] <- list()
      k <- 1
      for(j in names(x[[i]])) {
        msp <- single_eqn(x[[i]][[j]], y, id = c(i, j))
        teqn <- paste(model.name, ".hregress ", msp$eqn,
          ", family=", if(k < 2) fbx[[i]][1] else "gaussian_re",
          " equationtype=", fbx[[i]][2],
          if(n == length(x) & k < 2) {
            paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " "))
          } else NULL,
          if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
        eqn[[i]][[j]] <- teqn
        prgex <- c(prgex, msp$prgex)
        k <- k + 1
      }
      eqn[[i]] <- rev(eqn[[i]])
    } else {
      msp <- single_eqn(x[[i]], y, id = i)
      teqn <- paste(model.name, ".hregress ", msp$eqn, ", family=", fbx[[i]][1],
        " equationtype=", fbx[[i]][2],
        if(n == length(x)) {
          paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " "))
        } else NULL,
        if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
      eqn[[i]] <- teqn
      prgex <- c(prgex, msp$prgex)
    }
    n <- n + 1
  }

  prg <- c(prgex, "", paste("mcmcreg", model.name), "")
  for(i in unlist(rev(eqn)))
    prg <- c(prg, i, "")

  prg <- c(prg, paste(model.name, "getsample", sep = "."))
  prg <- c(
    paste('%% BayesX program created by bamlss: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )

  if(any(grepl("##seed##", prg, fixed = TRUE)))
    prg <- gsub("##seed##", round(runif(1L) * .Machine$integer.max), prg, fixed = TRUE)

  prgf <- file.path(dir, prg.name)
  writeLines(prg, prgf)

  require("BayesXsrc")

  warn <- getOption("warn")
  options(warn = -1)
  ok <- run.bayesx(prg = prgf, verbose = control$setup$verbose)
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
  if(length(i <- grep("-nan", ok$log, ignore.case = TRUE))){
    warning("the BayesX engine returned NA samples, please check your model specification! In some cases it can be helpful to center continuous covariates!", call. = FALSE)
  }

  sfiles <- grep("_sample.raw", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)

  samples <- NULL
  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
    } else {
      if(!is.null(x[[i]]$model.matrix)) {
        sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(paste("LinearEffects", sep = ""), sfiles, fixed = TRUE)
        sf <- sfiles[sf]
        samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
        colnames(samps) <- paste(i, ".p.", colnames(x[[i]]$model.matrix), sep = "")
        samples <- cbind(samples, samps)
      }
      if(!is.null(x[[i]]$smooth.construct)) {
        for(j in seq_along(x[[i]]$smooth.construct)) {
          term <- paste("_", paste(x[[i]]$smooth.construct[[j]]$term, collapse = "_"), "_", sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & !grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]
          samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
          cn <- colnames(x[[i]]$smooth.construct[[j]]$X)
          if(is.null(cn))
            cn <- paste("b", 1:ncol(x[[i]]$smooth.construct[[j]]$X), sep = "")
          colnames(samps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", cn, sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]
          if(length(sf)) {
            vsamps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
            colnames(vsamps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", paste("tau2", 1:ncol(vsamps), sep = ""), sep = "")
            samps <- cbind(samps, vsamps)
          }
          samples <- cbind(samples, samps)
        }
      }
    }
  }

  sf <- grep("_DIC", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)
  dic <- read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE]
  samples <- cbind(samples, "DIC" = dic$dic, "pd" = dic$pd)

  as.mcmc(samples)
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
  if(length(i <- grep("-nan", ok$log, ignore.case = TRUE))){
    warning("the BayesX engine returned NA samples, please check your model specification! In some cases it can be helpful to center continuous covariates!", call. = FALSE)
  }
  samples <- NULL
  samples
}


########################################
## (2) BayesX model term construction ##
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
    "rps", "hrandom_pspline"
  )
  if(!bs %in% available.terms) stop(paste("basis type", sQuote(bs), "not supported by BayesX"))

  if(bs %in% c("rsps", "hrandom_pspline")) {
    bs <- "rsps"
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
    rcall <- paste("r(x = ", by, ", bs = ", sQuote(bs), ", by = ", x, ", ...)", sep = "")
    rval <- eval(parse(text = rcall))
  } else {
    if(length(grep("~", term <- deparse(call$x))) && bs %in% c("re", "ra", "random")) {
      x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
      rcall <- paste("r(x = ", x, ", by = ", by, ", ...)", sep = "")
      rval <- eval(parse(text = rcall))
    } else {
      k <- -1
      m <- NA
      xt <- list(...)
      if(is.null(xt$lambda))
        xt$lambda <- 100
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
      rval <- if(bs != "te") {
        s(x, z, k = k, bs = bs, m = m, xt = xt)
      } else {
        te(x, z, k = k, bs = "ps", m = m, xt = xt, mp = FALSE)
      }
      rval$by <- by
      rval$term <- term
      rval$dim <- length(term)
      rval$sx.construct <- TRUE
      if(!(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
         "re", "pspline2dimrw2", "pspline2dimrw1"))) {
        class(rval) <- c(class(rval), "no.mgcv")
      }

      rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""),
        if(by != "NA") paste(",by=", by, sep = "") else NULL, ")", sep = "")
    }
  }

  return(rval)
}

sx.construct <- function(object, data, ...)
{
  UseMethod("sx.construct")
}

sx.construct.default <- function(object, data, ...) 
{
  cl <- grep(".smooth.spec", class(object), value = TRUE, fixed = TRUE)
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

sx.construct.userdefined.smooth.spec <- function(object, data, id, dir, ...)
{
  id <- paste(rmf(id), collapse = "_")
  term <- paste(object$term, collapse = "*")
  Sn <- paste(id, "S", sep = "_")
  Xn <- paste(id, "X", sep = "_")
  if(length(object$S) > 1) {
    object$S <- do.call("+", object$S)
    object$rank <- sum(object$rank)
  }
  if(is.null(object$rank))
    object$rank <- qr(object$S[[1]])$rank
  if(is.null(object$xt$centermethod))
    object$xt$centermethod <- "meanfd"
  term <- paste(term, "(userdefined,penmatdata=", Sn, ",designmatdata=", Xn,
    ",rankK=", object$rank, sep = "")
  term <- paste(do.xt(term, object, c("center", "before")), ")", sep = "")

  write <- function(dir) {
    write.table(object$S[[1]], file = file.path(dir, paste(Sn, ".raw", sep = "")),
      quote = FALSE, row.names = FALSE)
    write.table(object$X, file = file.path(dir, paste(Xn, ".raw", sep = "")),
      quote = FALSE, row.names = FALSE)
    c(
      paste("dataset", Sn),
      paste(Sn, ".infile using ", file.path(dir, paste(Sn, ".raw", sep = "")), sep = ""),
      paste("dataset", Xn),
      paste(Xn, ".infile using ", file.path(dir, paste(Xn, ".raw", sep = "")), sep = "")
    )
  }

  attr(term, "write") <- write

  term
}

sx.construct.pspline.smooth <- sx.construct.ps.smooth.spec <- sx.construct.psplinerw1.smooth.spec <-
sx.construct.psplinerw2.smooth.spec <- sx.construct.pspline.smooth.spec <-
function(object, data, ...)
{
  if(class(object)[1] == "pspline.smooth")
    class(object) <- "psplinerw2.smooth.spec"
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
sx.construct.random.smooth.spec <- function(object, data, ...)
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

sx.construct.rps.smooth.spec <- function(object, data, ...)
{
  term <- paste(object$term, "(hrandom_pspline,centermethod=meansum2", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  term <- paste(object$by, term , sep = "*")

  return(term)
}

sx.construct.kr.smooth.spec <- sx.construct.kriging.smooth.spec <- function(object, data, ...)
{
  termo <- object$term
  if(length(termo) < 2L)
    stop("kriging method needs two terms!")
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
	nrknots <- object$bs.dim
  xt <- object$xt
  if(is.null(xt$full))
    term <- paste(termo[1L], "*", termo[2L], "(kriging,nrknots=", nrknots, sep = "")
  else {
    term <- paste(termo[1L], "*", termo[2L], "(kriging,full", sep = "")    
    object$xt$full <- NULL
  }
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

construct.shrw <- function(object, data, what)
{
  term <- object$term
  term <- paste(term, "(", what, sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.offset.smooth.spec <- function(object, data, ...)
{
  return(construct.shrw(object, data, "offset"))
}

sx.construct.mrf.smooth.spec <- sx.construct.spatial.smooth.spec <- function(object, data, ...)
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


## Create a 'bamlss.frame'.
bamlss.frame <- function(formula, data = NULL, family = gaussian.bamlss(),
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, knots = NULL, specials = NULL, reference = NULL,
  model.matrix = TRUE, smooth.construct = TRUE, ytype = c("matrix", "vector"), ...)
{
  ## Parse formula.
  if(!inherits(formula, "bamlss.formula")) {
    ## Search for additional formulas.
    formula2 <- NULL; formula3 <- list(); k <- 1
    fn <- names(fo <- formals(fun = bamlss.frame))[-1]
    fn <- fn[fn != "..."]
    for(f in fn) {
      fe <- paste("deparse(substitute(", f, "), backtick = TRUE, width.cutoff = 500)", sep = "")
      fe <- eval(parse(text = fe))
      if(any(grepl("~", fe, fixed = TRUE))) {
        fe <- as.formula(fe, env = environment(if(is.list(formula)) formula[[1]] else formula))
        formula2 <- c(formula2, as.formula(fe))
        formula3[[k]] <- fe
        eval(parse(text = paste(f, if(is.null(fo[[f]])) "NULL" else fo[[f]], sep = " = ")))
        k <- k + 1
      }
    }
    formula <- c(formula, formula2)

    ## Parse family object.
    family <- bamlss.family(family)

    ## Parse formula.
    formula <- bamlss.formula(if(length(formula) < 2) formula[[1]] else formula, family)
  } else {
    family <- bamlss.family(family)
  }
  if(!is.null(attr(formula, "orig.formula")))
    formula <- attr(formula, "orig.formula")

  ## Setup return object.
  bf <- list()
  bf$call <- match.call()

  ## Create the model frame.
  bf$model.frame <- bamlss.model.frame(formula, data, family, weights,
    subset, offset, na.action, specials, contrasts)

  ## Type of y.
  ytype <- match.arg(ytype)

  ## Process categorical responses and assign 'y'.
  cf <- bamlss.formula.cat(formula, bf$model.frame, reference)
  if(!is.null(cf)) {
    rn <- response.name(terms(bf$model.frame), hierarchical = FALSE)[1]
    orig.formula <- formula
    formula <- cf$formula
    reference <- cf$reference
    if(ytype == "matrix") {
      f <- as.formula(paste("~ -1 +", rn))
      bf$y <- bf$model.frame[rn]
      bf$y[rn] <- model.matrix(f, data = bf$model.frame)
      colnames(bf$y[[rn]]) <- rmf(gsub(rn, "", colnames(bf$y[[rn]])))
      bf$y[[rn]] <- bf$y[[rn]][, rmf(c(names(formula), reference))]
    } else {
      bf$y <- bf$model.frame[rn]
      attr(bf$y[[rn]], "reference") <- reference
    }
    family$names <- names(formula)
    family$links <- rep(family$links, length.out = length(formula))
    names(family$links) <- names(formula)
    attr(formula, "orig.formula") <- orig.formula
  } else {
    rn <- response.name(formula, hierarchical = FALSE, keep.functions = TRUE)
    rn <- rn[rn %in% names(bf$model.frame)]
    bf$y <- bf$model.frame[rn]
    for(j in rn) {
      if(is.factor(bf$y[[j]]) & (ytype == "matrix")) {
        f <- as.formula(paste("~ -1 +", j))
        bf$y[j] <- model.matrix(f, data = bf$model.frame)
      }
    }
  }
  bf$formula <- formula

  ## Add the terms object.
  bf$terms <- terms.bamlss.formula(formula, data = data, ...)

  ## Process possible score and hess functions.
  if(!is.null(score <- family$score)) {
    if(is.function(score))
      score <- list(score)
    family$score <- rep(score, length.out = length(formula))
    names(family$score) <- names(formula)
  }
  if(!is.null(hess <- family$hess)) {
    if(is.function(hess))
      hess <- list(hess)
    family$hess <- rep(hess, length.out = length(formula))
    names(family$hess) <- names(formula)
  }

  ## Add more functions to family object.
  bf$family <- complete.bamlss.family(family)

  ## Assign the 'x' master object.
  bf$x <- design.construct(bf$terms, data = bf$model.frame, knots = knots,
    model.matrix = model.matrix, smooth.construct = smooth.construct, model = NULL, ...)
  bf$knots <- knots

  ## Assign class and return.
  class(bf) <- c("bamlss.frame", "list")

  return(bf)
}


## Simple print method for 'bamlss.frame'
print.bamlss.frame <- function(x, ...)
{
  cat("'bamlss.frame' structure:", "\n")  
  nx <- c("call", "model.frame", "formula", "family", "terms", "x", "y", "knots")
  nx <- c(nx, names(x)[!(names(x) %in% nx)])
  for(i in nx) {
    if(!is.null(x[[i]])) {
      cat("  ..$", i, "\n")
      if(i == "x") {
        for(j in names(x[[i]])) {
          cat("  .. ..$", j, "\n")
          if(!all(c("formula", "fake.formula") %in% names(x[[i]][[j]]))) {
            for(k in names(x[[i]][[j]])) {
              cat("  .. .. ..$", k, "\n")
              for(d in names(x[[i]][[j]][[k]])) {
               cat("  .. .. .. ..$", d, "\n")
              }
            }
          } else {
            for(k in names(x[[i]][[j]]))
              cat("  .. .. ..$", k, "\n")
          }
        }
      }
      if(i == "y") {
        for(j in names(x[[i]])) {
          cat("  .. ..$", j, "\n")
        }
      }
    }
  }
  invisible(NULL)
}


## Compute the 'bamlss.frame' 'x' master object.
design.construct <- function(formula, data = NULL, knots = NULL,
  model.matrix = TRUE, smooth.construct = TRUE, binning = FALSE,
  before = TRUE, gam.side = TRUE, model = NULL, drop = NULL, ...)
{
  if(!model.matrix & !smooth.construct)
    return(NULL)

  if(inherits(formula, "bamlss.frame")) {
    data <- if(is.null(data)) model.frame(formula) else data
    formula <- formula(formula)
  }
  if(!inherits(formula, "bamlss.terms")) {
    if(!inherits(formula, "bamlss.formula"))
      formula <- bamlss.formula(formula, ...)
    if(inherits(formula, "bamlss.formula"))
      formula <- terms.bamlss.formula(formula, data = data, ...)
  }
  formula <- formula.bamlss.terms(formula)
  if(is.null(data))
    stop("data needs to be supplied!")
  if(!inherits(data, "data.frame"))
    data <- as.data.frame(data)
  if(!is.null(model))
    formula <- model.terms(formula, model)
  if(!binning)
    binning <- NULL

  assign.design <- function(obj, dups = NULL)
  {
    if(!is.null(dups)) {
      if(any(dups)) {
        mi <- match.index(data[, all.vars(obj$fake.formula), drop = FALSE])
        obj[names(mi)] <- mi
        data <- subset(data, !dups)
      }
    }
    obj$binning <- binning
    if(!all(c("formula", "fake.formula") %in% names(obj)))
      return(obj)
    if(model.matrix) {
      obj$model.matrix <- model.matrix(drop.terms.bamlss(obj$terms,
        sterms = FALSE, keep.response = FALSE, data = data), data = data)
    }
    if(smooth.construct) {
      tx <- drop.terms.bamlss(obj$terms,
        pterms = FALSE, keep.response = FALSE, data = data)
      sid <- unlist(attr(tx, "specials"))
      if(!length(sid))
        sid <- NULL
      if(!is.null(sid)) {
        sterms <- sterm_labels <- attr(tx, "term.labels")[sid]
        sterms <- lapply(sterms, function(x) { eval(parse(text = x)) })
        nst <- NULL
        for(j in seq_along(sterms)) {
          sl <- sterms[[j]]$label
          if(is.null(sl))
            sl <- sterm_labels[j]
          nst <- c(nst, sl)
        }
        names(sterms) <- nst
        for(tsm in sterms) {
          if(is.null(tsm$xt))
            tsm$xt <- list()
          if(is.null(tsm$xt$binning))
            tsm$xt$binning <- binning
          if(!is.null(tsm$xt$binning)) {
            if(!is.logical(tsm$xt$binning)) {
              for(tsmt in tsm$term) {
                if(!is.factor(data[[tsmt]]))
                  data[[tsmt]] <- round(data[[tsmt]], digits = tsm$xt$binning)
              }
            }
          }
        }
        no.mgcv <- NULL
        smooth <- list()
        for(tsm in sterms) {
          if(is.null(tsm$special)) {
            if(is.null(tsm$xt))
              tsm$xt <- list()
            if(is.null(tsm$xt$binning))
              tsm$xt$binning <- binning
            acons <- TRUE
            if(!is.null(tsm$xt$center))
              acons <- tsm$xt$center
            tsm$xt$center <- acons
            tsm$xt$before <- before
            if(!is.null(tsm$xt$binning)) {
              term.names <- c(tsm$term, if(tsm$by != "NA") tsm$by else NULL)
              tsm$binning <- match.index(data[, term.names, drop = FALSE])
              tsm$binning$order <- order(tsm$binning$match.index)
              tsm$binning$sorted.index <- tsm$binning$match.index[tsm$binning$order]
              smt <- smoothCon(tsm, if(before) data[tsm$binning$nodups, term.names, drop = FALSE] else data,
                knots, absorb.cons = acons)
            } else {
              smt <- smoothCon(tsm, data, knots, absorb.cons = acons)
            }
          } else {
            smt <- smooth.construct(tsm, data, knots)
            if(inherits(smt, "no.mgcv")) {
              no.mgcv <- c(no.mgcv, list(smt))
              next
            } else {
              class(smt) <- c(class(smt), "mgcv.smooth")
              smt <- list(smt)
            }
          }
          smooth <- c(smooth, smt)
        }
        if(length(smooth) > 0) {
          if(gam.side) {
            if(is.null(obj$model.matrix)) {
              Xp <- model.matrix(drop.terms.bamlss(obj$terms,
                sterms = FALSE, keep.response = FALSE, data = data), data = data)
              smooth <- try(gam.side(smooth, Xp, tol = .Machine$double.eps^.5), silent = TRUE)
            } else {
              smooth <- try(gam.side(smooth, obj$model.matrix, tol = .Machine$double.eps^.5), silent = TRUE)
            }
            if(inherits(smooth, "try-error"))
              stop("gam.side() produces an error when binning, try to set before = FALSE or set gam.side = FALSE!")
          }
          sme <- NULL
          if(smooth.construct)
            sme <- expand.t2.smooths(smooth)
          if(is.null(sme)) {
            original.smooth <- NULL
          } else {
            original.smooth <- smooth
            smooth <- sme
            rm(sme)
          }
          if(!is.null(no.mgcv))
            smooth <- c(smooth, no.mgcv)
          if(length(smooth)) {
            stl <- NULL
            for(j in seq_along(smooth))
              stl <- c(stl, smooth[[j]]$label)
            names(smooth) <- stl
          }
          obj$smooth.construct <- smooth
        }
      }
    }
    if(!is.null(drop)) {
      take <- c("model.matrix", "smooth.construct")[c(model.matrix, smooth.construct)]
      obj[!(names(obj) %in% take)] <- NULL
    }

    obj
  }

  if(!all(c("formula", "fake.formula") %in% names(formula))) {
    for(j in seq_along(formula)) {
      if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
        for(i in seq_along(formula[[j]])) {
          formula[[j]][[i]] <- assign.design(formula[[j]][[i]],
            if(i > 1) duplicated(data[, all.vars(formula[[j]][[i]]$fake.formula), drop = FALSE]) else NULL)
        }
      } else formula[[j]] <- assign.design(formula[[j]])
    }
  } else formula <- assign.design(formula)

  if(!all(c("formula", "fake.formula") %in% names(formula))) {
    for(i in seq_along(formula)) {
      if(!all(c("formula", "fake.formula") %in% names(formula[[i]]))) {
        for(j in seq_along(formula[[i]])) {
          if(!is.null(formula[[i]][[j]]$smooth.construct)) {
            for(k in seq_along(formula[[i]][[j]]$smooth.construct)) {
              if(is.null(formula[[i]][[j]]$smooth.construct[[k]]$fit.fun))
                formula[[i]][[j]]$smooth.construct[[k]]$fit.fun <- make.fit.fun(formula[[i]][[j]]$smooth.construct[[k]])
            }
          }
        }
      } else {
        if(!is.null(formula[[i]]$smooth.construct)) {
          for(j in seq_along(formula[[i]]$smooth.construct)) {
            if(is.null(formula[[i]]$smooth.construct[[j]]$fit.fun))
              formula[[i]]$smooth.construct[[j]]$fit.fun <- make.fit.fun(formula[[i]]$smooth.construct[[j]])
          }
        }
      }
    }
  } else {
    if(!is.null(formula$smooth.construct)) {
      for(j in seq_along(formula$smooth.construct)) {
        if(is.null(formula$smooth.construct[[j]]$fit.fun))
          formula$smooth.construct[[j]]$fit.fun <- make.fit.fun(formula$smooth.construct[[j]])
      }
    }
  }

  attr(formula, "specials") <- NULL
  attr(formula, ".Environment") <- NULL
  class(formula) <- "list"
  if(!is.null(drop)) {
    if(drop & (length(formula) < 2))
      formula <- formula[[1]]
  }
  return(formula)
}


## The model term fitting function.
make.fit.fun <- function(x)
{
  ff <- function(X, b, expand = TRUE) {
    if(!is.null(names(b)))
      b <- get.par(b, "b")
    f <- drop(X %*% b)
    if(!is.null(x$binning$match.index) & expand)
      f <- f[x$binning$match.index]
    return(as.numeric(f))
  }
  return(ff)
}


## Get the model.frame.
model.frame.bamlss <- model.frame.bamlss.frame <- function(formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- if(length(nargs) || is.null(formula$model.frame)) {
    fcall <- formula$call
    fcall[[1L]] <- quote(bamlss.model.frame)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$formula)
    if(is.null(env))
      env <- parent.frame()
    ft <- eval(fcall[["formula"]], env)
    if(!is.null(attr(ft, "orig.formula"))) {
      fcall["formula"] <- parse(text = paste("attr(", fcall["formula"], ", 'orig.formula')", sep = ""))
    }
    fcall["start"] <- NULL
    fcall["drop.unused.levels"] <- FALSE
    if(is.null(fcall["family"]))
      fcall["family"] <- parse(text = "gaussian.bamlss()")
    eval(fcall, env)
  } else formula$model.frame
  mf
}


## Search for parts in models, optionally extract.
model.search <- function(x, what, model = NULL, part = c("x", "formula", "terms"),
  extract = FALSE, drop = FALSE)
{
  if(!inherits(x, "bamlss.formula") & !inherits(x, "bamlss.frame"))
    stop("x must be a 'bamlss.formula' or 'bamlss.frame' object!")
  part <- match.arg(part)
  if(is.null(x[[part]]))
    return(FALSE)
  x <- model.terms(x, model = model, part = part)
  elmts <- c("formula", "fake.formula")
  nx <- names(x)
  rval <- list()
  for(i in nx) {
    if(!all(elmts %in% names(x[[i]]))) {
      rval[[i]] <- list()
      for(j in names(x[[i]])) {
        rval[[i]][[j]] <- if(is.null(x[[i]][[j]][[what]])) FALSE else TRUE
        if(extract & rval[[i]][[j]])
          rval[[i]][[j]] <- x[[i]][[j]][[what]]
      }
    } else {
      rval[[i]] <- if(is.null(x[[i]][[what]])) FALSE else TRUE
      if(extract & rval[[i]])
        rval[[i]] <- x[[i]][[what]]
    }
  }
  if(!extract) {
    rval <- unlist(rval)
  } else {
    if(drop & (length(rval) < 2))
      rval <- rval[[1]]
  }
  rval
}


## Wrapper for design construct extraction.
extract.design.construct <- function(object, data = NULL,
  knots = NULL, model = NULL, drop = TRUE, what = c("model.matrix", "smooth.construct"))
{
  if(!inherits(object, "bamlss.frame") & !inherits(object, "bamlss.formula") & !inherits(object, "bamlss.terms"))
    stop("object must be a 'bamlss.frame', 'bamlss.formula' or 'bamlss.terms' object!")
  what <- match.arg(what)
  model.matrix <- what == "model.matrix"
  smooth.construct <- what == "smooth.construct"
  if(inherits(object, "bamlss.frame")) {
    if(!is.null(data)) {
      object$model.frame <- NULL
      object <- design.construct(object, data = data, knots = knots,
        model.matrix = model.matrix, smooth.construct = smooth.construct,
        model = model, drop = drop)
    } else {
      if(!all(model.search(object, what, model, part = "x"))) {
        object <- design.construct(object, model.matrix = model.matrix,
          smooth.construct = smooth.construct, model = model, drop = TRUE)
      } else {
        object <- model.search(object, what, model, extract = TRUE, drop = drop, part = "x")
      }
    }
  } else {
    if(is.null(data))
      stop("argument data is missing!")
    object <- design.construct(object, data = data, knots = knots,
      model.matrix = model.matrix, smooth.construct = smooth.construct, model = model, drop = drop)
  }
  if(!is.null(drop)) {
    if(length(object) & drop & (length(object) < 2))
      object <- object[[1]]
  }
  if(!length(object))
    return(NULL)
  mostattributes(object) <- NULL
  attr(object, "orig.formula") <- NULL
  return(object)
}


## Model matrix extractor.
model.matrix.bamlss.frame <- model.matrix.bamlss.formula <- model.matrix.bamlss.terms <- function(object, data = NULL, model = NULL, drop = TRUE, ...)
{
  extract.design.construct(object, data = data,
    knots = NULL, model = model, drop = drop, what = "model.matrix")
}


## Extract smooth constructs.
smooth.construct <- function(object, data, knots, ...)
{
  UseMethod("smooth.construct")
}

smooth.construct.bamlss.frame <- smooth.construct.bamlss.formula <- smooth.construct.bamlss.terms <- function(object, data = NULL, knots = NULL, model = NULL, drop = TRUE, ...)
{
  extract.design.construct(object, data = data,
    knots = knots, model = model, drop = drop, what = "smooth.construct")
}


## Extract/initialize parameters.
parameters <- function(x, model = NULL, start = NULL, fill = c(0, 0.0001), list = TRUE, simple.list = FALSE)
{
  if(inherits(x, "bamlss")) {
    if(!is.null(x$parameters)) {
      if(is.null(model)) {
        if(list) return(x$parameters) else return(unlist(x$parameters))
      } else {
        if(list) return(x$parameters[model]) else return(unlist(x$parameters[model]))
      }
    }
  }
  if(inherits(x, "bamlss.frame")) {
    if(is.null(x$x)) {
      x <- design.construct(x, data = x$model.frame,
        knots = x$knots, model.matrix = TRUE, smooth.construct = TRUE, model = NULL)
    } else x <- x$x
  }
  fill <- rep(fill, length.out = 2)
  if(!is.null(start)) {
    if(is.list(start))
      start <- unlist(start)
  }
  par <- list()
  for(i in names(x)) {
    par[[i]] <- list()
    if(!all(c("formula", "fake.formula") %in% names(x[[i]]))) {
      for(j in names(x[[i]])) {
        par[[i]][[j]] <- list()
        if(!is.null(x[[i]][[j]]$model.matrix)) {
          nc <- ncol(x[[i]][[j]]$model.matrix)
          if(simple.list) {
            par[[i]][[j]]$p <- fill[1]
          } else {
            par[[i]][[j]]$p <- rep(fill[1], length = nc)
            if(is.null(cn <- colnames(x[[i]][[j]]$model.matrix)))
              cn <- paste("b", 1:nc, sep = "")
            names(par[[i]][[j]]$p) <- cn
            if(!is.null(start)) {
              if(length(ii <- grep(paste(i, j, "p", sep = "."), names(start), fixed = TRUE))) {
                spar <- start[ii]
                spn <- names(spar)
                cn2 <- paste(i, j, "p", cn, sep = ".")
                take <- which(spn %in% cn2)
                if(length(take)) {
                  par[[i]][[j]]$p[which(cn2 %in% spn)] <- spar[take]
                }
              }
            }
          }
        }
        if(!is.null(x[[i]][[j]]$smooth.construct)) {
          par[[i]][[j]]$s <- list()
          for(k in names(x[[i]][[j]]$smooth.construct)) {
            if(simple.list) {
              par[[i]][[j]]$s[[k]] <- fill[1]
            } else {
              if(!is.null(x[[i]][[j]]$smooth.construct[[k]]$rand)) {
                tpar1 <- rep(fill[1], ncol(x[[i]][[j]]$smooth.construct[[k]]$rand$Xr))
                tpar2 <- rep(fill[1], ncol(x[[i]][[j]]$smooth.construct[[k]]$Xf))
                names(tpar1) <- paste("b", 1:length(tpar1), ".re", sep = "")
                names(tpar2) <- paste("b", 1:length(tpar2), ".fx", sep = "")
                tpar <- c(tpar1, tpar2)
              } else {
                tpar <- rep(fill[1], ncol(x[[i]][[j]]$smooth.construct[[k]]$X))
                names(tpar) <- paste("b", 1:length(tpar), sep = "")
              }
              if(length(x[[i]][[j]]$smooth.construct[[k]]$S)) {
                tpar3 <- NULL
                for(kk in seq_along(x[[i]][[j]]$smooth.construct[[k]]$S)) {
                  tpar3 <- c(tpar3, fill[2])
                }
                names(tpar3) <- paste("tau2", 1:length(tpar3), sep = ".")
                tpar <- c(tpar, tpar3)
              }
              par[[i]][[j]]$s[[k]] <- tpar
              if(!is.null(start)) {
                if(length(ii <- grep(paste(i, j, "s", k, sep = "."), names(start), fixed = TRUE))) {
                  spar <- start[ii]
                  cn <- names(par[[i]][[j]]$s[[k]])
                  if(length(tau2 <- grep("tau2", names(spar)))) {
                    tau2 <- spar[tau2]
                    if(length(jj <- grep("tau2", cn, fixed = TRUE))) {
                      tau2 <- rep(tau2, length.out = length(jj))
                      par[[i]][[j]]$s[[k]][jj] <- tau2
                    }
                  }
                  if(length(b <- grep("b", names(spar)))) {
                    b <- spar[b]
                    if(length(jj <- grep("b", cn, fixed = TRUE))) {
                      b <- rep(b, length.out = length(jj))
                      par[[i]][[j]]$s[[k]][jj] <- b
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      if(!is.null(x[[i]]$model.matrix)) {
        if(simple.list) {
          par[[i]]$p <- fill[1]
        } else {
          nc <- ncol(x[[i]]$model.matrix)
          par[[i]]$p <- rep(fill[1], length = nc)
          if(is.null(cn <- colnames(x[[i]]$model.matrix)))
            cn <- paste("b", 1:nc, sep = "")
          names(par[[i]]$p) <- cn
          if(!is.null(start)) {
            if(length(ii <- grep(paste(i, "p", sep = "."), names(start), fixed = TRUE))) {
              spar <- start[ii]
              spn <- names(spar)
              cn2 <- paste(i, "p", cn, sep = ".")
              take <- which(spn %in% cn2)
              if(length(take)) {
                par[[i]]$p[which(cn2 %in% spn)] <- spar[take]
              }
            }
          }
        }
        if(!is.null(x[[i]]$smooth.construct)) {
          par[[i]]$s <- list()
          for(k in names(x[[i]]$smooth.construct)) {
            if(!is.null(x[[i]]$smooth.construct[[k]]$rand)) {
              tpar1 <- rep(fill[1], ncol(x[[i]]$smooth.construct[[k]]$rand$Xr))
              tpar2 <- rep(fill[1], ncol(x[[i]]$smooth.construct[[k]]$Xf))
              names(tpar1) <- paste("b", 1:length(tpar1), ".re", sep = "")
              names(tpar2) <- paste("b", 1:length(tpar2), ".fx", sep = "")
              tpar <- c(tpar1, tpar2)
            } else {
              tpar <- rep(fill[1], ncol(x[[i]]$smooth.construct[[k]]$X))
              names(tpar) <- paste("b", 1:length(tpar), sep = "")
            }
            if(length(x[[i]]$smooth.construct[[k]]$S)) {
              tpar3 <- NULL
              for(kk in seq_along(x[[i]]$smooth.construct[[k]]$S)) {
                tpar3 <- c(tpar3, fill[2])
              }
              names(tpar3) <- paste("tau2", 1:length(tpar3), sep = ".")
              tpar <- c(tpar, tpar3)
            }
            par[[i]]$s[[k]] <- tpar
            if(!is.null(start)) {
              if(length(ii <- grep(paste(i, "s", k, sep = "."), names(start), fixed = TRUE))) {
                spar <- start[ii]
                cn <- names(par[[i]]$s[[k]])
                if(length(tau2 <- grep("tau2", names(spar)))) {
                  tau2 <- spar[tau2]
                  if(length(jj <- grep("tau2", cn, fixed = TRUE))) {
                    tau2 <- rep(tau2, length.out = length(jj))
                    par[[i]]$s[[k]][jj] <- tau2
                  }
                }
                if(length(b <- grep("b", names(spar)))) {
                  b <- spar[b]
                  if(length(jj <- grep("b", cn, fixed = TRUE))) {
                    b <- rep(b, length.out = length(jj))
                    par[[i]]$s[[k]][jj] <- b
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(!is.null(model))
    par <- par[model]
  if(!list)
    par <- unlist(par)
  return(par)
}


## Main bamlss().
bamlss <- function(formula, family = gaussian.bamlss, data = NULL, start = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  reference = NULL, transform = NULL, optimizer = NULL, sampler = NULL, results = NULL,
  cores = NULL, sleep = NULL, combine = TRUE, model = TRUE, x = TRUE, ...)
{
  ## The environment.
  env <- get_formula_envir(formula)

  ## Search for functions in family object.
  family <- bamlss.family(family)
  if(!is.null(family$transform))
    transform <- family$transform
  if(!is.null(family$optimizer))
    optimizer <- family$optimizer
  if(!is.null(family$sampler))
    sampler <- family$sampler
  if(!is.null(family$results))
    results <- family$results

  ## Setup all processing functions.
  foo <- list("transform" = transform, "optimizer" = optimizer, "sampler" = sampler, "results" = results)
  nf <- names(foo)
  default_fun <- c("no.transform", "bfit", "GMCMC", "results.bamlss.default")
  functions <- list()
  for(j in 1:length(foo)) {
    if(is.null(foo[[j]]))
      foo[[j]] <- eval(parse(text = default_fun[j]))
    if(is.list(foo[[j]])) {
      args <- foo[[j]]
      fun <- default_fun[j]
      functions[[j]] <- function(x, ...) {
        args <- c(args, list(...))
        args$x <- x
        do.call(fun, args)
      }
    } else functions[[j]] <- foo[[j]]
    if(!is.function(functions[[j]])) {
      if(!is.logical(functions[[j]]))
        stop(paste("argument", nf[j], "is not a function!"))
    }
  }
  names(functions) <- names(foo)

  ## Create the 'bamlss.frame'.
  bf <- match.call(expand.dots = TRUE)
  bf[c("transform", "optimizer", "sampler", "results", "cores", "sleep", "combine", "model", "x")] <- NULL
  bf[[1]] <- as.name("bamlss.frame")
  bf <- eval(bf, envir = env)

  ## Transform.
  if(is.function(functions$transform)) {
    tbf <- functions$transform(bf)
    bf[names(tbf)] <- tbf
    rm(tbf)
  }

  ## Start optimizer.
  if(is.function(functions$optimizer)) {
    opt <- functions$optimizer(x = bf$x, y = bf$y, family = bf$family,
      start = start, weights = model.weights(bf$model.frame),
      offset = model.offset(bf$model.frame), ...)
    bf[names(opt)] <- opt[names(opt)]
    rm(opt)
  }

  ## Start sampling.
  if(is.function(functions$sampler)) {
    if(is.null(cores)) {
      bf$samples <- functions$sampler(x = bf$x, y = bf$y, family = bf$family,
        weights = model.weights(bf$model.frame),
        offset = model.offset(bf$model.frame),
        start = if(is.null(bf$parameters)) start else unlist(bf$parameters), ...)
    } else {
      require("parallel")
      parallel_fun <- function(j) {
        if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
        functions$sampler(x = bf$x, y = bf$y, family = bf$family,
          weights = model.weights(bf$model.frame),
          offset = model.offset(bf$model.frame),
          start = if(is.null(bf$parameters)) start else unlist(bf$parameters), ...)
      }
      bf$samples <- mclapply(1:cores, parallel_fun, mc.cores = cores)
    }
    ## Process samples.
    bf$samples <- process.chains(bf$samples, combine)
  }

  ## Compute results.
  if(is.function(functions$results))
    bf$results <- functions$results(bf, ...)

  ## Save the model frame?
  if(!model)
    bf$model.frame <- NULL

  ## Save 'x' master object?
  if(!x)
    bf$x <- NULL

  bf$call <- match.call()
  class(bf) <- c("bamlss", "bamlss.frame", "list")

  bf
}


## No transform function.
no.transform <- TRUE

## family extractor.
family.bamlss <- family.bamlss.frame <- function(object, ...)
{
  return(object$family)
}




#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
## Could be interesting: http://people.duke.edu/~neelo003/r/
##                       http://www.life.illinois.edu/dietze/Lectures2012/
xreg <- function(formula, family = gaussian.bamlss, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  reference = NULL, parse.input = bamlss.frame, transform = transformJAGS,
  setup = setupJAGS, engine = samplerJAGS, results = resultsJAGS,
  cores = NULL, sleep = NULL, combine = TRUE, model = TRUE, grid = 100, ...)
{
  ## The environment.
  ef <- environment(if(is.list(formula)) formula[[1]] else formula)

  ## Setup all processing functions.
  if(is.null(transform))
    transform <- function(x) { x }
  foo <- list("transform" = transform, "setup" = setup, "engine" = engine, "results" = results)
  nf <- names(foo)
  default_fun <- c("randomize", "setupJAGS", "samplerJAGS", "resultsJAGS")
  functions <- list()
  for(j in 1:length(foo)) {
    if(is.list(foo[[j]])) {
      args <- foo[[j]]
      fun <- default_fun[j]
      functions[[j]] <- function(x, ...) {
        args <- c(args, list(...))
        args$x <- x
        do.call(fun, args)
      }
    } else functions[[j]] <- foo[[j]]
    if(!is.function(functions[[j]])) {
      if(!is.logical(functions[[j]]))
        stop(paste("argument", nf[j], "is not a function!"))
    }
  }
  names(functions) <- names(foo)
  functions$parse.input <- if(!is.null(parse.input)) {
    stopifnot(is.function(parse.input))
    deparse(substitute(parse.input), backtick = TRUE, width.cutoff = 500)
  } else "bamlss.frame"

  ## Parse input.
  pm <- match.call(expand.dots = TRUE)
  pm$parse.input <- pm$setup <- pm$samples <- pm$results <- NULL
  pm[[1]] <- as.name(functions$parse.input)
  pm <- eval(pm, parent.frame())
  attr(pm, "environment") <- new.env(parent = ef)

  ## Transform inputs.
  pm <- functions$transform(pm)

  ## Sampling setup.
  if(is.logical(functions$setup)) {
    sc <- FALSE
  } else {
    sc <- TRUE
    ms <- functions$setup(pm)
  }

  ## Start sampling.
  if(is.null(cores)) {
    so <- functions$engine(if(sc) ms else pm)
  } else {
    require("parallel")
    parallel_fun <- function(j) {
      if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
      functions$engine(if(sc) ms else pm)
    }
    so <- mclapply(1:cores, parallel_fun, mc.cores = cores)
  }

  ## Combine samples.
  if(combine)
    so <- process.chains(so)

  ## Compute results.
  rval <- functions$results(pm, so)
  rm(so)

  ## Save the model frame?
  if(!model)
    rval$model.frame <- NULL

  rval
}


#########################
## (2) Engine stacker. ##
#########################
stacker <- function(x, optimizer = bfit0, sampler = samplerJAGS, ...)
{
  if(is.function(optimizer) | is.character(optimizer))
    optimizer <- list(optimizer)
  if(is.integer(sampler) | is.numeric(sampler)) {
    n.samples <- as.integer(sampler)
    sampler <- function(x, ...) { null.sampler(x, n.samples = n.samples) }
  }
  if(is.null(sampler))
    sampler <- null.sampler
  if(is.function(sampler) | is.character(sampler))
    sampler <- list(sampler)
  if(length(optimizer)) {
    for(j in optimizer) {
      if(is.character(j)) j <- eval(parse(text = j))
      if(!is.function(j)) stop("the optimizer must be a function!")
      x <- j(x, ...)
    }
  }
  if(length(sampler)) {
    for(j in sampler) {
      if(is.character(j)) j <- eval(parse(text = j))
      if(!is.function(j)) stop("the sampler must be a function!")
      x <- j(x, ...)
    }
  }

  x
}


#########################
## (3) BAMLSS wrapper. ##
#########################
bamlss99 <- function(formula, family = gaussian, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  optimizer = list(bfit0), sampler = list(GMCMC), results = resultsBayesG,
  engine = NULL, cores = NULL, sleep = 1, combine = TRUE,
  n.iter = 12000, thin = 10, burnin = 2000, seed = NULL, ...)
{
  ff <- try(inherits(family, "family.bamlss"), silent = TRUE)
  if(inherits(ff, "try-error")) {
    family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)
  } else {
    if(is.function(family)) {
      if(inherits(try(family(), silent = TRUE), "try-error"))
        family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)
    }
  }
  family.bamlss <- if(is.function(family)) family() else {
    if(is.character(family)) {
      if(!grepl("gF(", family, fixed = TRUE) & !grepl("gF2(", family, fixed = TRUE))
          if(!grepl("bamlss", family))
            family <- paste(family, "bamlss", sep = ".")
      family <- eval(parse(text = family[1]))
      if(is.function(family))
        family()
      else family
    } else family
  }

  if(is.null(engine)) {
    mc.cores <- cores
    transform <- if(!is.null(family.bamlss$transform)) {
      function(x) { family.bamlss$transform(x, ...) }
    } else function(x) { bamlss.setup(x, ...) }
    if(!is.null(family.bamlss$sampler)) {
      sampler <- function(x, ...) {
        family.bamlss$sampler(x, cores = mc.cores,
          n.iter = n.iter, thin = thin, burnin = burnin, seed = seed, sleep = sleep, ...)
      }
    }
    if(!is.null(family.bamlss$setup))
      setup <- function(x, ...) { family.bamlss$setup(x, ...) }
    else
      setup <- FALSE
    if(is.null(sampler))
      sampler <- function(x, ...) { null.sampler(x, ...) }
    if(!is.null(family.bamlss$engine)) {
      engine <- function(x) {
        family.bamlss$engine(x, cores = mc.cores,
          n.iter = n.iter, thin = thin, burnin = burnin, seed = seed, sleep = sleep, ...)
      }
      cores <- NULL
      xengine <- "in.family"
    } else {
      engine <- function(x) {
        stacker(x, optimizer = optimizer, sampler = sampler, cores = mc.cores,
          n.iter = n.iter, thin = thin, burnin = burnin, seed = seed, sleep = sleep, ...)
      }
      setup <- FALSE
      cores <- NULL
      xengine <- "stacker"
    }
  } else {
    xengine <- c("BayesG", "BayesX", "JAGS", "STAN")
    xengine <- xengine[pmatch(engine, xengine)]

    if(xengine == "BayesX") {
      require("BayesXsrc")

      data.name <- if(!is.null(data)) {
        deparse(substitute(data), backtick = TRUE, width.cutoff = 500)
      } else "d"

      transform <- transformBayesX
      setup <- function(x) {
        setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
          seed = seed, data.name = data.name, cores = cores, ...)
      }
      engine <- samplerBayesX
      results <- resultsBayesX
    }

    if(xengine == "BayesG") {
      transform <- transformBayesG
      setup <- FALSE
      engine <- function(x) {
        BayesG(x, n.iter = n.iter, thin = thin,
          burnin = burnin, seed = seed, ...)
      }
      results <- resultsBayesG
    }
  
    if(xengine %in% c("JAGS", "STAN")) {
      transform <- transformBUGS
      if(xengine == "JAGS") {
        require("rjags")
        setup <- setupJAGS
        engine <- function(x) {
          samplerJAGS(x, n.iter = n.iter, thin = thin,
            burnin = burnin, seed = seed, ...)
        }
      } else {
        require("rstan")
        setup <- bugs2stan
        engine <- function(x) {
          samplerSTAN(x, n.iter = n.iter, thin = thin,
            burnin = burnin, seed = seed, ...)
        }  
      }

      results <- resultsJAGS
    }
  }

  rval <- xreg(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = bamlss.frame, transform = transform,
    setup = setup, engine = engine, results = results, cores = cores,
    combine = combine, sleep = sleep, ...)
  
  rval$call <- match.call()
  rval$engine <- xengine
  
  rval
}


"[.bamlss" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  mostattributes(rval) <- attributes(x)
  rval
}


## Create the model.frame.
bamlss.model.frame <- function(formula, data, family = gaussian.bamlss(),
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  specials = NULL, contrasts.arg = NULL, drop.unused.levels = TRUE, ...)
{
  if(inherits(formula, "bamlss.frame") | inherits(formula, "bamlss")) {
    if(!is.null(formula$model.frame))
      return(formula$model.frame)
    fcall <- formula$call
    fcall[[1L]] <- quote(bamlss.model.frame)
    env <- environment(formula$formula)
    if(is.null(env))
      env <- parent.frame()
    return(eval(fcall, env))
  } else {
    family <- bamlss.family(family)
    formula <- bamlss.formula(formula, family)
    env <- environment(formula)
  }

  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  if(missing(data))
    data <- environment(formula)
  if(!is.data.frame(data))
    data <- as.data.frame(data)

  ## Make fake "Formula" object.
  fF <- make_fFormula(formula)

  ## Resulting terms object.
  mterms <- terms(formula(fF), data = data)

  ## Set up the model.frame.
  data <- list(formula = fF, data = data, subset = subset,
    na.action = na.action, drop.unused.levels = drop.unused.levels, ...)

  data <- do.call("model.frame", data)
  rownames(data) <- NULL

  ## Code from stats model.matrix()
  contr.funs <- as.character(getOption("contrasts"))
  namD <- names(data)
  for(i in namD) {
    if(is.character(data[[i]])) 
      data[[i]] <- factor(data[[i]])
  }
  isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
  isF[attr(mterms, "response")] <- FALSE
  isOF <- vapply(data, is.ordered, NA)
  for(nn in namD[isF]) {
    if(is.null(attr(data[[nn]], "contrasts"))) {
      contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
    }
  }
  if(!is.null(contrasts.arg) && is.list(contrasts.arg)) {
    if(is.null(namC <- names(contrasts.arg))) 
      stop("invalid 'contrasts' argument")
    for(nn in namC) {
      if(is.na(ni <- match(nn, namD))) 
        warning(gettextf("variable '%s' is absent, its contrast will be ignored", nn), domain = NA)
      else {
        ca <- contrasts.arg[[nn]]
        if(is.matrix(ca)) {
          contrasts(data[[ni]], ncol(ca)) <- ca
        } else {
          contrasts(data[[ni]]) <- contrasts.arg[[nn]]
        }
      }
    }
  }

  ## Process weights and offset.
  if(!is.null(weights)) {
    if(!is.list(weights))
      weights <- list(weights)
    weights <- do.call("cbind", weights)
    colnames(weights) <- names(formula)[1:ncol(weights)]
    if(!is.null(subset)) {
      weights <- if(!is.logical(subset)) {
        weights[subset, , drop = FALSE]
      } else subset(weights, subset)
    }
    if(nrow(weights) < 2)
      weights <- do.call(rbind, replicate(nrow(data), weights, simplify = FALSE))
    data[["(weights)"]] <- weights
  }
  if(!is.null(offset)) {
    if(!is.list(offset))
      offset <- list(offset)
    offset <- do.call("cbind", offset)
    colnames(offset) <- names(formula)[1:ncol(offset)]
    if(!is.null(subset)) {
      offset <- if(!is.logical(subset)) {
        offset[subset, , drop = FALSE]
      } else subset(offset, subset)
    }
    if(nrow(offset) < 2)
      offset <- do.call(rbind, replicate(nrow(data), offset, simplify = FALSE))
    data[["(offset)"]] <- offset
  }

  ## Remove inf values.
  data <- rm_infinite(data)

  ## Assign terms object.
  attr(data, "terms") <- mterms

  ## Check response.
  if(!is.null(family$valid.response)) {
    family$valid.response(model.response(data))
  }

  data
}


## Remove Inf values from data.
rm_infinite <- function(x) {
  if(is.null(dim(x))) return(x)
  if(ncol(x) > 0) {
    for(j in 1:ncol(x)) {
      if(any(class(x[, j]) %in% c("numeric", "integer"))) {
        if(any(!is.finite(x[, j]))) {
          warning("infinite values in data, removing these observations in model.frame!")
          x <- x[is.finite(x[, j]), ]
        }
      }
    }
  }
  x
}


## Parse families and get correct family object, depending on type.
bamlss.family <- function(family, type = "bamlss")
{
  family <- if(is.function(family)) family() else {
    if(is.character(family)) {
      if(!is.null(type)) {
        if(!grepl("gF(", family, fixed = TRUE) & !grepl("gF2(", family, fixed = TRUE))
          if(!grepl(type, family))
            family <- paste(family, type, sep = ".")
      }
      family <- eval(parse(text = family[1]))
      if(is.function(family))
        family()
      else family
    } else family
  }
  if(!inherits(family, "family.bamlss")) {
    if(!is.character(family)) {
      if(is.null(family$family)) stop("family is specified wrong, no family name available!")
      family <- family$family
    }
    txt <- paste(family, type, sep = if(!is.null(type)) "." else "")
    txt <- gsub("bamlss.bamlss", "bamlss", txt, fixed = TRUE)
    family <- eval(parse(text = txt[1]))
    family <- family()
  }
  if(is.null(family)) family <- list()

  family
}


complete.bamlss.family <- function(family)
{
  if(is.null(family$map2par)) {
    linkinv <- vector(mode = "list", length = length(family$names))
    for(j in family$names)
      linkinv[[j]] <- make.link2(family$links[j])$linkinv
    family$map2par <- function(eta) {
      for(j in names(eta)) {
        eta[[j]] <- linkinv[[j]](eta[[j]])
        eta[[j]][is.na(eta[[j]])] <- 0
        if(any(jj <- eta[[j]] == Inf))
          eta[[j]][jj] <- 10
        if(any(jj <- eta[[j]] == -Inf))
          eta[[j]][jj] <- -10
      }
      return(eta)
    }
  }
  if(is.null(family$mu)) {
    family$mu <- function(par) { make.link2(family$links[1])$linkinv(par[[1]]) }
  }
  if(is.null(family$loglik)) {
    if(!is.null(family$d))
      family$loglik <- function(y, par, ...) { sum(family$d(y, par, log = TRUE), na.rm = TRUE) }
  }

  return(family)
}


## Formula to list().
as.list.Formula <- function(x)
{
  if(!inherits(x, "Formula"))
    x <- as.Formula(x)
  env <- environment(x)
  lhs <- attr(x, "lhs")
  rhs <- attr(x, "rhs")
  nl <- length(lhs)
  nr <- length(rhs)
  if(nl < nr)
    lhs <- c(lhs, rep(list(NA), length = nr - nl))
  if(nr < nl)
    rhs <- c(rhs, rep(list(1), length = nl - nr))
  x <- mapply(c, lhs, rhs, SIMPLIFY = FALSE)
  formula <- list()
  for(i in seq_along(x)) {
    check <- inherits(x[[i]][[1]], "call") | inherits(x[[i]][[1]], "name")
    f <- if(check) {
      as.call(c(as.symbol("~"), x[[i]]))
    } else as.call(c(as.symbol("~"), x[[i]][[2]]))
    formula[[i]] <- eval(f, envir = env)
  }
  environment(formula) <- env
  formula
}


## Special formula parser, can deal with multi parameter models
## and hierarchical structures.
bamlss.formula <- function(formula, family = NULL)
{
  if(inherits(formula, "bamlss.formula"))
    return(formula)
  if(!is.list(formula)) {
    if(!inherits(formula, "Formula"))
      formula <- as.Formula(formula)
    if(inherits(formula, "Formula"))
      formula <- as.list.Formula(formula)
  }
  if(!is.null(family))
    family <- bamlss.family(family)
  if(!is.list(formula)) formula <- list(formula)
  if(!length(formula)) stop("formula is specified wrong!")
  
  env <- get_formula_envir(formula)

  complete_formula <- function(formula) {
    if(!is.null(family)) {
      if(length(formula) < length(family$names))
        formula <- c(formula, rep(list(), length = length(family$names) - length(formula)))
    }
    fn <- NULL
    for(j in seq_along(formula)) {
      ft <- if(!inherits(formula[[j]], "formula")) formula[[j]][[1]] else formula[[j]]
      if(!is.null(ft)) {
        yok <- attr(terms(formula(as.Formula(ft), rhs = FALSE)), "response") > 0
        fn <- c(fn, if(yok) all.vars(ft)[1] else NULL)
      }
    }
    fn[fn %in% c("1", "-1")] <- NA
    nas <- which(is.na(fn))
    if(!is.null(family)) {
      fn[nas] <- family$names[nas]
      if(is.null(family$names))
        family$names <- NA
      if(!is.na(family$names[1]))
        fn[1] <- family$names[1]
      else
        family$names <- fn
    } else fn[nas] <- paste("par", 1:length(fn[nas]), sep = ".")
    if(is.null(family)) {
      if(length(fn) < length(formula)) {
        k <- length(formula) - length(fn)
        if(k > 1)
          fn <- c(fn, paste("?par", 1:k, sep = ""))
        else
          fn <- c(fn, "?par")
      }
    }
    names(formula) <- fn
    if(!is.null(family)) {
      if(any(i <- is.na(names(formula))))
        names(formula)[i] <- family$names[i]
      for(j in family$names) {
        if(is.null(formula[[j]])) {
          formula[[j]] <- as.formula(paste(j, "1", sep = " ~ "))
          environment(formula[[j]]) <- env
        }
      }
    }
    for(j in seq_along(formula)) {
      if(!inherits(formula[[j]], "formula")) {
        if(is.null(names(formula[[j]])))
          names(formula[[j]]) <- paste("h", 1:length(formula[[j]]), sep = "")
      }
    }
    formula
  }
  formula <- formula_and(formula, env)
  formula <- formula_at(formula, env)
  formula <- complete_formula(formula_hierarchical(formula))
  formula <- formula_extend(formula, family)

  environment(formula) <- env
  class(formula) <- c("bamlss.formula", "list")

  formula
}


## Formula environment.
get_formula_envir <- function(formula)
{
  env <- environment(formula)
  if(is.null(env)) {
    get_env <- function(x) {
      if(inherits(x, "list")) {
        env <- NULL
        for(j in x)
          env <- c(env, get_env(j))
        return(env)
      } else return(environment(x))
    }
    env <- get_env(formula)
  }
  if(is.null(env)) env <- .GlobalEnv
  if(is.list(env))
    env <- env[[1]]
  return(env)
}


## Process categorical responses.
bamlss.formula.cat <- function(formula, data, reference)
{
  env <- environment(formula)

  rn <- y <- NULL
  for(j in seq_along(formula)) {
    ft <- if(!inherits(formula[[j]]$formula, "formula")) {
      formula[[j]][[1]]$formula
    } else formula[[j]]$formula
    yok <- attr(terms(formula(as.Formula(ft), rhs = FALSE)), "response") > 0
    if(yok)
      rn <- c(rn, all.vars(ft)[1])
  }
  rn2 <- rn[rn %in% names(data)]
  if(length(rn2) < 2) {
    if(is.factor(data[[rn2[1]]])) {
      if(nlevels(data[[rn2[1]]]) > 2) {
        ft <- as.formula(paste("~ -1 +", rn2[1]))
        y <- model.matrix(ft, data = data)
        colnames(y) <- rmf(gsub(rn2[1], "", colnames(y), fixed = TRUE))
        if(is.null(reference)) {
          ty <- table(data[[rn2[1]]])
          reference <- c(names(ty)[ty == max(ty)])[1]
        } else {
          ld <- levels(data[[rn2[1]]])
          reference <- ld[match(gsub(rn2[1], "", reference), ld)]
        }
        if(is.na(reference))
          stop(paste("cannot find reference category within response levels!"))
        reference <- rmf(reference)
        ylevels <- rmf(levels(data[[rn2[1]]]))
        ylevels <- ylevels[ylevels != reference]
        y <- y[, colnames(y) %in% ylevels, drop = FALSE]
        if(length(formula) < ncol(y)) {
          formula <- c(formula, rep(formula, length = ncol(y) - length(formula)))
        }
        if(!(names(formula)[[1]] %in% colnames(y))) {
          names(formula)[[1]] <- colnames(y)[1]
          ft <- if(!inherits(formula[[1]]$formula, "formula")) {
            formula[[1]][[1]]$formula
          } else formula[[1]]$formula
          env <- environment(ft)
          ft <- update(ft, as.formula(paste(colnames(y)[1], ".", sep = "~")))
          environment(ft) <- env
          if(!inherits(formula[[1]]$formula, "formula")) {
            formula[[1]][[1]]$formula <- ft
            formula[[1]][[1]]$response <- colnames(y)[1]
          } else {
            formula[[1]]$formula <- ft
            formula[[1]]$response <- colnames(y)[1]
          }
          ft <- if(!inherits(formula[[1]]$formula, "formula")) {
            formula[[1]][[1]]$fake.formula
          } else formula[[1]]$fake.formula
          ft <- update(ft, as.formula(paste(colnames(y)[1], ".", sep = "~")))
          if(!inherits(formula[[1]]$formula, "formula")) {
            formula[[1]][[1]]$fake.formula <- ft
          } else {
            formula[[1]]$fake.formula <- ft
          }
        }
        if(length(i <- !(names(formula) %in% ylevels))) {
          k <- 1
          ynot <- ylevels[!(ylevels %in% names(formula))]
          for(j in which(i)) {
            names(formula)[[j]] <- ynot[k]
            ft <- if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
              formula[[j]][[1]]$formula
            } else formula[[j]]$formula
            env <- environment(ft)
            ft <- update(ft, as.formula(paste(ynot[k], "~ .")))
            environment(ft) <- env
            if(!inherits(formula[[j]]$formula, "formula")) {
              formula[[j]][[1]]$formula <- ft
              formula[[j]][[1]]$response <- ynot[k]
            } else {
              formula[[j]]$formula <- ft
              formula[[j]]$response <- ynot[k]
            }
            attr(formula[[j]], "name") <- ynot[k]
            k <- k + 1
          }
        }
      }
    }
  }

  rval <- if(!is.null(y)) {
    class(formula) <- "bamlss.formula"
    environment(formula) <- env
    list("formula" = formula, "reference" = reference)
  } else NULL
  rval
}


## Make "Formula" object from fake.formulas.
make_fFormula <- function(formula)
{
  fF <- NULL
  for(j in seq_along(formula)) {
    if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
      for(i in seq_along(formula[[j]]))
        fF <- c(fF, formula[[j]][[i]]$fake.formula)
    } else {
      fF <- c(fF, formula[[j]]$fake.formula)
    }
  }
  fF <- do.call("as.Formula", fF)
  fF
}


all.vars.formula <- function(formula, lhs = TRUE, rhs = TRUE, specials = NULL, intercept = FALSE)
{
  env <- environment(formula)
  tf <- terms(formula, specials = unique(c("s", "te", "t2", "sx", "s2", "rs", "ti", specials)),
    keep.order = TRUE)
  sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  vars <- NULL
  if(rhs) {
    if(length(sid)) {
      vars <- tl[-sid]
      if(!length(vars))
        vars <- NULL
      for(j in tl[sid]) {
        st <- try(eval(parse(text = j), envir = env), silent = TRUE)
        if(inherits(st, "try-error"))
          st <- eval(parse(text = j), enclos = env, envir = loadNamespace("mgcv"))
        vars <- c(vars, st$term)
        if(st$by != "NA") {
          vars <- if(grepl("~", st$by)) {
            c(vars, all.vars.formula(as.formula(st$by, env = env)))
          } else c(vars, st$by)
        }
      }
    } else vars <- tl
  }
  if(lhs & (attr(tf, "response") > 0))
    vars <- c(vars, response.name(formula, keep.functions = TRUE))
  if(intercept & (attr(tf, "intercept") > 0))
    vars <- c("1", vars)
  unique(vars)
}


all.labels.formula <- function(formula, specials = NULL, full.names = FALSE)
{
  env <- environment(formula)
  tf <- terms(formula, specials = unique(c("s", "te", "t2", "sx", "s2", "rs", "ti", specials)),
    keep.order = TRUE)
  sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  labs <- NULL

  if(length(sid)) {
    labs <- tl[-sid]
    if(full.names & length(labs))
      labs <- paste("p", labs, sep = ".")
    if(!length(labs))
      labs <- NULL
    for(j in tl[sid]) {
      st <- try(eval(parse(text = j), envir = env), silent = TRUE)
      if(inherits(st, "try-error"))
        st <- eval(parse(text = j), enclos = env, envir = loadNamespace("mgcv"))
      labs <- c(labs, if(full.names) paste("s", st$label, sep = ".") else st$label)
    }
  } else labs <- if(full.names) paste("p", tl, sep = ".") else tl

  unique(labs)
}


fake.formula <- function(formula, lhs = TRUE, rhs = TRUE, specials = NULL)
{
  if(!lhs & !rhs)
    return(0 ~ 0)
  if(rhs)
    f <- paste(all.vars.formula(formula, lhs = FALSE, rhs = TRUE, specials, intercept = TRUE), collapse = "+")
  if(lhs)
    f <- paste(all.vars.formula(formula, lhs = TRUE, rhs = FALSE), "~", if(!is.null(f)) f else 0)
  else
    f <- paste("~", if(!is.null(f)) f else 0)
  f <- as.formula(f, env = environment(formula))
  f
}


## Extend formula by a fake formula with all variables
## to compute a model.frame, create smooth objects.
formula_extend <- function(formula, family)
{
  if(is.list(formula)) {
    for(j in seq_along(formula))
      formula[[j]] <- formula_extend(formula[[j]], family)
    return(formula)
  } else {
    rn <- response.name(formula)
    ff <- fake.formula(formula, lhs = !(rn %in% family$names))
    return(list("formula" = formula, "fake.formula" = ff))
  }
}


## Get response name.
response.name <- function(formula, hierarchical = TRUE, keep.functions = FALSE)
{
  rn <- NA
  if(inherits(formula, "bamlss.frame"))
    formula <- terms(model.frame(formula))
  if(!is.null(attr(formula, "terms")))
    formula <- attr(formula, "terms")
  if(inherits(formula, "formula")) {
    f <- as.Formula(formula)
    f <- formula(f, lhs = TRUE, rhs = FALSE)
    if(keep.functions) {
      cf <- as.character(formula)
      rn <- if(length(cf) < 3) character(0) else cf[2]
    } else {
      rn <- all.vars(f)
    }
  } else {
    if(inherits(formula, "list")) {
      rn <- NULL
      for(i in seq_along(formula)) {
        if(is.null(formula[[i]]$formula) & inherits(formula[[i]], "list")) {
          for(j in seq_along(formula[[i]])) {
            if(!hierarchical & (j > 1)) {
              next
            } else {
              tf <- if(is.null(formula[[i]][[j]]$formula)) {
                formula[[i]][[j]]
              } else formula[[i]][[j]]$formula
              rn <- c(rn, response.name(tf, keep.functions = keep.functions))
            }
          }
        } else {
          rn <- c(rn , response.name(formula[[i]]$formula, keep.functions = keep.functions))
        }
      }
    }
  }
  if(!length(rn))
    rn <- NA
  rn
}

## Search and process "&"
formula_and <- function(formula, env = parent.frame())
{
  if(nol <- !is.list(formula))
    formula <- list(formula)
  for(j in seq_along(formula)) {
    if(!inherits(formula[[j]], "formula")) {
      formula[[j]] <- formula_and(formula[[j]], env)
    } else {
      ft <- deparse(formula[[j]])
      if(any(grep("&", ft, fixed = TRUE))) {
        formula[[j]] <- as.list(strsplit(ft, "&", fixed = TRUE)[[1]])
        for(i in seq_along(formula[[j]])) {
          if(!any(grepl("~", formula[[j]][[i]], fixed = TRUE)))
            formula[[j]][[i]] <- paste("~", formula[[j]][[i]])
          formula[[j]][[i]] <- as.formula(formula[[j]][[i]], env = env)
        }
      }
    }
  }
  if(nol) {
    names(formula) <- response.name(formula[[1]][[1]])
    if(inherits(formula[[1]], "formula"))
      formula <- formula[[1]]
  }
  formula
}

## Search and process "@"
formula_at <- function(formula, env = parent.frame())
{
  if(nol <- !is.list(formula))
    formula <- list(formula)
  for(j in seq_along(formula)) {
    if(!inherits(formula[[j]], "formula")) {
      formula[[j]] <- formula_at(formula[[j]], env)
    } else {
      ft <- deparse(formula[[j]])
      if(any(grep("@", ft, fixed = TRUE))) {
        formula[[j]] <- strsplit(ft, "@", fixed = TRUE)[[1]]
        control <- formula[[j]][-1]
        formula[[j]] <- as.formula(formula[[j]][1], env = env)
        control <- gsub(":", "=", control)
        if(any(grepl("+", control)))
          control <- strsplit(control, "+", fixed = TRUE)[[1]]
        control <- gsub("^ +", "", control)
        control <- gsub(" +$", "", control)
        attr(formula[[j]], "control") <- gsub("using=", "using ", control, fixed = TRUE)
      }
    }
  }
  if(nol) {
    names(formula) <- response.name(formula[[1]][[1]])
    if(inherits(formula[[1]], "formula"))
      formula <- formula[[1]]
  }
  formula
}

formula_rm_at <- function(formula, env = parent.frame())
{
  ctr <- attr(formula, "control")
  attr(formula, "control") <- NULL
  if(isf <- !is.character(formula)) {
    env <- environment(formula)
    formula <- deparse(formula)
  }
  if(any(grepl("@", formula)))
    formula <- strsplit(formula, "@")[[1]][1]
  if(!isf) {
    formula <- as.formula(formula, env = env)
    formula <- deparse(formula)
  }
  if(isf)
    formula <- as.formula(formula, env = env)
  if(!is.null(ctr)) attr(formula, "control") <- ctr
  formula
}

## Hierarchical formulae.
formula_hcheck <- function(formula)
{
  if(!is.list(formula))
    return(formula)
  check <- vector(mode = "list", length = length(formula))
  for(j in seq_along(formula)) {
      for(i in seq_along(formula)) {
        if(j != i) {
          fi <- if(!is.list(formula[[i]])) list(formula[[i]]) else formula[[i]]
          rnj <- response.name(formula[[j]])
          for(jj in seq_along(fi)) {
            av <- all.vars(fi[[jj]])
            rn <- response.name(fi[[jj]])
            if(!is.na(rn))
              av <- av[av != rn]
            if(!has_dot(fi[[jj]])) {
              if(attr(terms(fi[[jj]]), "intercept") < 1) {
                av <- c(av, "-1")
              }
            }
            if(any(av %in% rnj)) {
              check[[j]] <- c(check[[j]], i)
            }
          }
        }
      }
  }
  check
}

formula_insert <- function(from, to, formula)
{
  nf <- names(formula)
  hm <- sapply(to, max)
  o <- order(hm, decreasing = TRUE)
  from <- from[o]
  to <- to[o]
  for(j in seq_along(from)) {
    for(i in seq_along(to[[j]])) {
      formula[[to[[j]][i]]] <- c(formula[[to[[j]][i]]], formula[[from[j]]])
    }
  }
  formula <- formula[take <- !(1:length(formula) %in% from)]
  names(formula) <- nf[take]
  formula
}

formula_hierarchical <- function(formula)
{
  if(!is.list(formula))
    return(formula)
  j <- formula_hcheck(formula)
  while(any(!sapply(j, is.null))) {
    i <- which(!sapply(j, is.null))
    formula <- formula_insert(i, j[i], formula)
    j <- formula_hcheck(formula)
  }
  formula
}


## Transform smooth terms to mixed model representation.
randomize <- function(x)
{
  if(!inherits(x, "bamlss.frame"))
    stop("object must be a 'bamlss.frame'!")
  if(is.null(x$x))
    stop("no 'x' object to randomize in 'bamlss.frame'!")

  vnames <- names(model.frame(x))
  x <- x$x

  rand_fun <- function(x)
  {
    if(m <- length(x$smooth.construct)) {
      for(j in 1:m) {
        if(!inherits(x$smooth.construct[[j]], "no.mgcv")) {
          if(is.null(x$smooth.construct[[j]]$rand) & is.null(x$smooth.construct[[j]]$Xf)) {
            tmp <- smooth2random(x$smooth.construct[[j]], vnames, type = 2)
            if(is.null(x$smooth.construct[[j]]$xt$nolin))
              x$smooth.construct[[j]]$Xf <- tmp$Xf
#          if(inherits(x$smooth.construct[[j]], "random.effect")) {
#            tmp$rand$Xr[tmp$rand$Xr > 0] <- 1
#            tmp$rand$Xr <- scale(tmp$rand$Xr)
#            tmp$trans.D <- rep(1, ncol(tmp$rand$Xr))
#            tmp$trans.U <- diag(1, ncol(tmp$rand$Xr))
#          }
            x$smooth.construct[[j]]$Xnore <- x$smooth.construct[[j]]$X
            x$smooth.construct[[j]]$X <- tmp$rand$Xr
            x$smooth.construct[[j]]$trans.D <- tmp$trans.D
            x$smooth.construct[[j]]$trans.U <- tmp$trans.U
            if(!is.null(x$smooth.construct[[j]]$state$parameters)) {
              g2 <- get.par(x$smooth.construct[[j]]$state$parameters, "g")
              if(!is.null(x$smooth.construct[[j]]$trans.U))
                g2 <- solve(x$smooth.construct[[j]]$trans.U) %*% g2
              g2 <- drop(g2 / x$smooth.construct[[j]]$trans.D)
              x$smooth.construct[[j]]$state$parameters <- set.par(x$smooth.construct[[j]]$state$parameters, g2, "g")
            }
          }
        }
      }
    }
    x
  }

  elmts <- c("formula", "fake.formula")
  for(j in seq_along(x)) {
    if(!all(elmts %in% names(x[[j]]))) {
      for(i in seq_along(x[[j]]))
        x[[j]][[i]] <- rand_fun(x[[j]][[i]])
    } else x[[j]] <- rand_fun(x[[j]])
  }
  return(list("x" = x))
}


## Combine sample chains.
process.chains <- function(x, combine = TRUE, drop = FALSE)
{
  if(!is.list(x))
    x <- list(x)
  model.specs <- attr(x[[1]], "model.specs")
  if(inherits(x[[1]], "mcmc.list")) {
    x <- as.mcmc.list(do.call("c", x))
  } else {
    stopifnot(inherits(x[[1]], "mcmc"))
    x <- as.mcmc.list(x)
  }
  if(combine) {
    x <- do.call("rbind", x)
    x <- as.mcmc.list(list(as.mcmc(x)))
  }
  if(drop & (length(x) < 2))
    x <- x[[1]]
  attr(x, "model.specs") <- model.specs
  return(x)
}


## Combine method for "bamlss" objects.
c.bamlss <- function(...)
{
  objects <- list(...)
  x <- NULL
  for(i in 1L:length(objects))
    x <- c(x, objects[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- c("bamlss", "cbamlss")

  return(x)
}


## Fast computation of quantiles.
quick_quantiles <- function(X, samples)
{
  rval <- .Call("quick_quantiles", X, samples)
  rval <- as.data.frame(rval)
  names(rval) <- c("2.5%", "50%", "97.5%")
  rval
}

fitted_matrix <- function(X, samples)
{
  fit <- .Call("fitted_matrix", X, as.matrix(samples))
}


## Function to compute statistics from samples of a model term.
compute_term <- function(x, get.X, fit.fun, psamples, vsamples = NULL,
  asamples = NULL, FUN = NULL, snames, s.effects.resmat, data,
  grid = 100, rug = TRUE, hlevel = 1, sx = FALSE, re.slope = FALSE, edfsamples = NULL)
{
  require("coda")

  nt <- length(x$term)

  if(x$by != "NA" | nt > 1) grid <- NA

  tterms <- NULL
  for(l in nt:1) {
    tterm <- x$term[l]
    for(char in c("(", ")", "[", "]"))
      tterm <- gsub(char, ".", tterm, fixed = TRUE)
    if(inherits(data[[tterm]], "ts"))
      data[[tterm]] <- as.numeric(data[[tterm]])
    tterms <- c(tterms, tterm)
  }

  for(char in c("(", ")", "[", "]"))
    colnames(data) <- gsub(char, ".", colnames(data), fixed = TRUE)

  ## Data for rug plotting.
  rugp <- if(length(x$term) < 2 & rug) data[[x$term]] else NULL

  if(hlevel > 1)
    data <- unique(data)

  ## Handling factor by variables.
  if(x$by != "NA") {
    if(nlevels(data[[x$by]]) > 2 & FALSE)
      x$by <- "NA"
    if(sx) {
      if(all(data[[x$by]] %% 1 == 0)) {
        if(length(unique(data[[x$by]])) < length(data[[x$by]]))
          x$by <- "NA"
      }
      if(re.slope) {
        x$by <- "NA"
      }
    }
  }

  ## New x values for which effect should
  ## be calculated, n = 100.
  if(!is.na(grid)) {
    if(length(x$term) < 2 & !is.factor(data[[tterms[1]]]) & !any(grepl("mrf", class(x))) &
      !any(grepl("re.", class(x), fixed = TRUE)) & !any(grepl("random", class(x))) & !re.slope) {
      xsmall <- TRUE
      nd <- list()
      for(j in tterms) {
        xr <- range(data[[j]], na.rm = TRUE)
        nd[[j]] <- seq(xr[1], xr[2], length = grid)
      }
      if(x$by != "NA") { ## FIXME: check by variables!
        if(!is.factor(data[[x$by]])) {
          xr <- range(data[[x$by]], na.rm = TRUE)
          nd[[x$by]] <- seq(xr[1], xr[2], length = 100)
        } else nd[[x$by]] <- rep(data[[x$by]], length.out = grid)
      }
      data0 <- data
      data <- as.data.frame(nd)
    } else xsmall <- FALSE
  } else {
    data0 <- data[, c(tterms, if(x$by != "NA") x$by else NULL), drop = FALSE]
    if(nt < 2) {
      if(x$by != "NA") {
        if(!is.factor(data[[x$by]]))
          data <- unique(data0)
      } else data <- unique(data0)
    }
    xsmall <- if((nrow(data) != nrow(data0)) & (nt < 2)) TRUE else FALSE
  }
  if(is.null(x$special)) {
    X <- get.X(data)
  } else {
    if(x$special)
      X <- get.X(data[, tterms, drop = FALSE])
    else
      X <- get.X(data)
  }

  ## Compute samples of fitted values.
  if(is.null(FUN)) {
    FUN <- function(x) {
      rval <- as.numeric(quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      names(rval) <- c("2.5%", "50%", "97.5%")
      rval
    }
  }
  if(inherits(x, "mgcv.smooth") & nrow(psamples) > 39L) {
    smf <- quick_quantiles(X, psamples)
  } else {
    if(nt < 2) {
      fsamples <- apply(psamples, 1, function(g) {
        fit.fun(X, g, expand = FALSE)
      })
      smf <- t(apply(fsamples, 1, FUN))
    } else {
      smf <- 0
      for(i in 1:nrow(psamples)) {
        smf <- smf + drop(fit.fun(X, psamples[i, ], expand = FALSE))
      }
      smf <- as.matrix(smf / nrow(psamples), ncol = 1)
      colnames(smf) <- "50%"
    }
  }

  cnames <- colnames(smf)
  smf <- as.data.frame(smf)
  for(l in 1:nt) {
    smf <- cbind(data[[tterms[l]]], smf)
  }
  names(smf) <- c(x$term, cnames)

  if(x$by != "NA") { ## FIXME: hard coded fix for plotting varying coefficient terms!
    if(!is.factor(data[[x$by]])) {
      X <- X / data[[x$by]]

      fsamples <- apply(psamples, 1, function(g) {
        fit.fun(X, g, expand = FALSE)
      })

      smf <- t(apply(fsamples, 1, FUN))

      cnames <- colnames(smf)
      smf <- as.data.frame(smf)
      for(l in 1:nt) {
        smf <- cbind(data[[tterms[l]]], smf)
      }
      names(smf) <- c(x$term, cnames)
    }
  }

  by.drop <- NULL
  if(x$by != "NA" & !is.null(x$by.level)) {
    by.drop <- (if(xsmall) data0[[x$by]] else data[[x$by]]) == x$by.level
    if(!xsmall)
      smf <- smf[by.drop, ]
  }

  ## Assign class and attributes.
  smf <- unique(smf)
  if(is.factor(data[, tterms])) {
    bbb <- 1 ## FIXME: factors!
  }
  class(smf) <- c(class(x), "data.frame")
  x[!(names(x) %in% c("term", "label", "bs.dim", "dim"))] <- NULL
  attr(smf, "specs") <- x
  class(attr(smf, "specs")) <- class(x)
  attr(smf, "x") <- if(xsmall & nt < 2) data0[, tterms] else data[, tterms]
  attr(smf, "by.drop") <- by.drop
  attr(smf, "rug") <- rugp
  colnames(psamples) <- paste(x$label, 1:ncol(psamples), sep = ".")

  ## Get samples of the variance parameter.
  edf <- FALSE
  vsamples0 <- NULL
  if(!is.null(edfsamples)) {
    vsamples0 <- vsamples
    vsamples <- edfsamples
    edf <- TRUE
  }
  if(!is.null(vsamples)) {
    if(!is.matrix(vsamples))
      vsamples <- matrix(vsamples, ncol = 1)
    ## edf <- apply(vsamples, 1, function(x) {} )
    smatfull <- NULL
    for(j in 1:ncol(vsamples)) {
      qu <- drop(quantile(vsamples[, j], probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      sd <- sd(vsamples[, j], na.rm = TRUE)
      me <- mean(vsamples[, j], na.rm = TRUE)
      smat <- matrix(c(me, sd, qu), nrow = 1)
      colnames(smat) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
      rownames(smat) <- x$label
      if(!is.null(asamples)) {
        smat <- cbind(smat, "alpha" = mean(asamples))
      }
      smatfull <- rbind(smatfull, smat)
    }
    if(!is.null(smatfull)) {
      if(any(duplicated(rownames(smatfull)))) {
        rownames(smatfull) <- paste(rownames(smatfull), 1:nrow(smatfull), sep = ":")
      }
    }
    if(!is.null(smf)) {
      attr(smf, "scale") <- smatfull
      attr(smf, "specs")$label <- gsub(")", paste(",",
        paste(formatC(me, digits = 2), collapse = ","), ")",
        sep = ""), x$label)
      colnames(vsamples) <- paste(x$label, if(is.null(edfsamples)) "tau2" else "edf", 1:nrow(smatfull), sep = ".")
      attr(smf, if(is.null(edfsamples)) "samples.scale" else "samples.edf") <- as.mcmc(vsamples)
      if(!is.null(vsamples0)) {
        if(!is.matrix(vsamples0)) vsamples0 <- matrix(vsamples0, ncol = 1)
        colnames(vsamples0) <- paste(x$label, "tau2", 1:ncol(vsamples0), sep = ".")
        attr(smf, "samples.scale") <- as.mcmc(vsamples0)
      }
      if(!is.null(asamples)) {
        asamples <- matrix(asamples, ncol = 1)
        colnames(asamples) <- paste(x$label, "alpha", sep = ".")
        attr(smf, "samples.alpha") <- as.mcmc(asamples)
      }
    }
    s.effects.resmat <- rbind(s.effects.resmat, smatfull)
    if(edf)
      attr(s.effects.resmat, "edf") <- TRUE
  }

  return(list("term" = smf, "s.effects.resmat" = s.effects.resmat))
}


## Function to add partial residuals based on weights() and score() function.
add.partial <- function(x, samples = FALSE, nsamps = 100)
{
  if(!inherits(x, "bamlss"))
    stop("x must be a 'bamlss' object!")

  nx <- names(x$terms)
  family <- x$family

  if(!is.null(family$hess) & !is.null(family$score)) {
    y <- model.response(model.frame(x))
    eta <- fitted.bamlss(x, samples = samples, nsamps = nsamps)
    if(is.null(x$model.frame))
      x$model.frame <- model.frame(x)
    for(j in seq_along(x$terms)) {
      if(!is.null(x$terms[[j]]$effects)) {
        peta <- family$map2par(eta)
        weights <- family$hess[[nx[j]]](y, peta, id = nx[j])
        score <- family$score[[nx[j]]](y, peta, id = nx[j])
        z <- eta[[nx[j]]] + 1 / weights * score
        ne <- names(x$terms[[j]]$effects)
        for(sj in seq_along(ne)) {
          f <- predict.bamlss(x, model = nx[j], term = ne[sj], nsamps = nsamps)
          term <- attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$term
          e <- z - eta[[nx[j]]] + f
          if(is.null(attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$xt$center)) {
            e <- e - mean(e)
          } else {
            if(attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$xt$center)
              e <- e - mean(e)
          }
          e <- data.frame(mf[, term], e)
          names(e) <- c(term, "partial.resids")
          attr(x$terms[[j]]$effects[[ne[sj]]], "partial.resids") <- e
        }
      }
    }
  } else {
    stop("cannot compute partial residuals, no score() and hess() function in family object!")
  }
  x
}


## A prediction method for "bamlss" objects.
## Prediction can also be based on multiple chains.
predict.bamlss <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, type = c("link", "parameter"), FUN = function(x) { mean(x, na.rm = TRUE) },
  trans = NULL, what = c("samples", "parameters"), nsamps = NULL, verbose = FALSE, drop = TRUE, ...)
{
  if(is.null(object$x))
    object$x <- smooth.construct(object)
  family <- object$family
  if(missing(newdata)) {
    newdata <- model.frame(object)
  } else {
    if(is.character(newdata)) {
      if(file.exists(newdata <- path.expand(newdata)))
        newdata <- read.table(newdata, header = TRUE, ...)
      else stop("cannot find newdata")
    }
    if(is.matrix(newdata) || is.list(newdata))
      newdata <- as.data.frame(newdata)
    ## FIXME: ??? newdata <- model.frame.bamlss.frame(object, data = newdata)
  }
  if(!is.null(attr(object, "fixed.names")))
    names(newdata) <- rmf(names(newdata))
  nn <- names(newdata)
  rn <- response.name(object, keep.functions = TRUE)
  nn <- nn[nn != rn]
  tl <- term.labels2(object, model = model, intercept = intercept)
  nx <- names(tl)
  if(!is.null(term)) {
    enames <- vector(mode = "list", length = length(nx))
    for(j in term) {
      for(i in seq_along(tl)) {
        if(!is.character(j)) {
          if(j > 0 | j < length(tl[[i]]))
            enames[[i]] <- c(enames[[i]], tl[[i]][j])
        } else {
          if(grepl("intercept", tolower(j), fixed = TRUE)) {
            if("(Intercept)" %in% tl[[i]])
              enames[[i]] <- c(enames[[i]], "(Intercept)")
          } else {
            k <- grep(j, tl[[i]], fixed = TRUE)
            if(length(k))
              enames[[i]] <- c(enames[[i]], tl[[i]][k])
          }
        }
      }
    }
    names(enames) <- nx
  } else enames <- tl
  ff <- as.formula(paste("~", paste(unlist(enames), collapse = "+")))
  vars <- all.vars.formula(ff)
  if(!all(vars[vars != "Intercept"] %in% nn))
    stop("cannot compute prediction, variables missing in newdata!")
  type <- match.arg(type)
  what <- match.arg(what)
  if(!is.null(object$samples) & what == "samples") {
    samps <- samples(object, model = model)
    if(!is.null(nsamps)) {
      i <- seq(1, nrow(samps), length = nsamps)
      samps <- samps[i, , drop = FALSE]
    }
  } else {
    if(is.null(object$parameters))
      stop("cannot find any parameters!")
    samps <- parameters(object, model = model, list = FALSE)
    cn <- names(samps)
    samps <- matrix(samps, nrow = 1)
    colnames(samps) <- cn
    samps <- as.mcmc(samps)
  }

  env <- environment(object$formula)
  enames <- lapply(enames, function(x) {
    f <- as.formula(paste("~", paste(x, collapse = "+")), env = env)
    all.labels.formula(f, full.names = TRUE)
  })

  pred <- list()
  for(i in nx) {
    pred[[i]] <- .predict.bamlss(i, object$x[[i]], samps, enames[[i]],
      intercept, FUN, trans, type, nsamps, newdata, env)
  }

  if(type != "link") {
    links <- family$links[nx]
    if(length(links) > 0) {
      for(i in nx) {
        if(links[i] != "identity") {
          linkinv <- make.link2(links[i])$linkinv
          pred[[i]] <- linkinv(pred[[i]])
        }
      }
    } else {
      warning(paste("could not compute predictions on the scale of parameter",
        ", predictions on the scale of the linear predictor are returned!", sep = ""))
    }
  }
  if(!is.null(trans)) {
    if(!is.list(trans)) {
      trans <- rep(list(trans), length = length(nx))
      names(trans) <- nx
    }
    for(i in nx) {
      if(!is.null(trans[[i]])) {
        if(!is.function(trans[[i]]))
          stop("argument trans must be a list of transformer functions!")
        pred[[i]] <- trans[[i]](pred[[i]])
      }
    }
  }
  for(i in nx) {
    pred[[i]] <- apply(pred[[i]], 1, FUN, ...)
    if(!is.null(dim(pred[[i]])))
      pred[[i]] <- t(pred[[i]])
  }

  if((length(pred) < 2) & drop)
    pred <- pred[[1]]

  return(pred)
}


.predict.bamlss <- function(id, x, samps, enames, intercept, FUN, trans, type, nsamps, data, env)
{
  if("smooth.construct" %in% names(x))
    x <- x$smooth.construct
  snames <- colnames(samps)
  enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
  has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."), snames, fixed = TRUE))
  if(intercept & has_intercept)
    enames <- c("p.(Intercept)", enames)
  enames <- unique(enames)
  ec <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][1:2], collapse = "")
  })
  enames2 <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][-c(1:2)], collapse = "")
  })
  eta <- matrix(0, nrow = nrow(data), ncol = nrow(samps))
  if(length(i <- grep("p.", ec))) {
    for(j in enames2[i]) {
      if(j != "(Intercept)") {
        f <- as.formula(paste("~", if(has_intercept) "1" else "-1", "+", j), env = env)
        X <- model.matrix(f, data = data)
        if(has_intercept)
          X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        sn <- snames[grep2(paste(id, "p", j, sep = "."), colnames(samps), fixed = TRUE)]
        eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
      } else {
        sn <- snames[grep2(paste(id, "p", j, sep = "."), colnames(samps), fixed = TRUE)]
        eta <- eta + fitted_matrix(matrix(1, nrow = nrow(data), ncol = 1), samps[, sn, drop = FALSE])
      }
    }
  }
  if(length(i <- grep("s.", ec))) {
    for(j in enames2[i]) {
      if(!inherits(x[[j]], "no.mgcv") & !inherits(x[[j]], "special")) {
        X <- PredictMat(x[[j]], data)
        sn <- snames[grep2(paste(id, "s", j, sep = "."), colnames(samps), fixed = TRUE)]
        eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
      } else {
        stop("no predictions for special terms available yet!")
      }
    }
  }
  eta
}


## Setup function for handling "special" model terms.
s2 <- function(...)
{
  rval <- s(...)
  rval$special <- TRUE
  rval$label <- gsub("s(", "s2(", rval$label, fixed = TRUE)
  rval
}


## Setup function for random scaling terms.
rsc <- function(..., by = NA)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  if(by != "NA") {
    if(!grepl("~", by, fixed = TRUE)) {
      if(by == ".") 
        stop("by=. not allowed")
      by <- paste("~", by)
    }
  }
  rval <- s(...)
  rval$by.formula <- if(by != "NA") as.formula(by) else NULL
  rval$class <- class(rval)
  rval$special <- TRUE
  class(rval) <- "rsc.smooth.spec"
  rval
}


## Smooth constructor function for random scaling terms.
smooth.construct.rsc.smooth.spec <- function(object, data, knots) {
  class(object) <- object$class
  acons <- TRUE
  if(!is.null(object$xt$center))
    acons <- object$xt$center
  rval <- smoothCon(object, data, knots, absorb.cons = acons)
  rval <- rval[[1]]
  rval$class <- class(rval)
  if(!is.null(object$by.formula)) {
    ft <- terms(object$by.formula, keep.order = TRUE)
    vars <- attr(ft, "term.labels")
    if(length(vars)) {
      rs.by <- list()
      for(j in vars) {
        rs.by[[j]] <- data[[j]]
        if(!is.factor(rs.by[[j]])) stop("random scaling by variables must be factors!")
      }
      n <- length(vars)
      g <- paste("g[", 1:n, "]", sep = "", collapse = " + ")
      fun <- paste("function(X, g) { ", "(", if(attr(ft, "intercept")) "1 + ",
        g, ") * (X %*% ", if(n > 1) paste("g[-(1:", n, ")]", sep = "") else "g[-1]", ") }", sep = "")
      rval$fit.fun <- eval(parse(text = fun))
      rval$rs.by <- rs.by
      rval$by.vars <- vars
      rval$by.formula <- object$by.formula
      rval$one <- attr(ft, "intercept")
    }
  } else {
    rval$fit.fun <- function(X, g, ...) {
      X %*% as.numeric(g)
    }
  }

  class(rval) <- "rsc.smooth"
  rval
}


## Smooth constructor function for growth curves.
smooth.construct.gc.smooth.spec <- function(object, data, knots) 
{
  object$X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  object$special <- TRUE
  if(is.null(object$xt))
    object$xt <- list("center" = FALSE)
  else
    object$xt$center <- TRUE
  object$by.done <- TRUE
  if(object$by != "NA")
    stop("by variables not implemented yet!")

  object$fit.fun <- function(X, g, ...) {
    ##f <- g[1] * exp(-g[2] * exp(-g[3] * drop(X)))
    f <- log(g[1]) - g[2] * exp(-g[3] * drop(X))
    f <- exp(f)
    if(object$xt$center)
      f <- f - mean(f)
    f
  }
  object$update <- bfit0_optim
  object$propose <- gmcmc_sm.slice
  object$prior <- function(gamma) {
    sum(dnorm(gamma, sd = 1000, log = TRUE))
  }
  object$edf <- function(x) { 3 }
  object$fixed <- TRUE
  object$np <- 3
  object$state$parameters <- rep(0, 3)
  names(object$state$parameters) <- paste("g", 1:3, sep = "")
  object$state$fitted.values <- object$fit.fun(object$X, object$state$parameters)
  object$state$edf <- 3
  object$lower <- rep(-Inf, 3)
  object$upper <- rep(Inf, 3)
  names(object$lower) <- names(object$upper) <- names(object$state$parameters)

  class(object) <- c("gc.smooth", "no.mgcv", "special")
  object
}

## Work around for the "prediction matrix" of a growth curve.
Predict.matrix.gc.smooth <- function(object, data, knots) 
{
  X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  X
}


## Rational spline constructor.
rs <- function(..., k = -1, fx = NULL, bs = "tp", m = NA, xt = NULL, link = "log", specials = NULL)
{
  smooths <- as.list(substitute(list(...)))[-1]
  if(any(grepl("~", as.character(smooths[[1]]), fixed = TRUE))) {
    stop("formulae no supported yet!")
    specials <- c(specials, "s")
    if(length(smooths) != 2) stop("there must be exactly 2 formulas!")
    sm <- list()
    for(j in seq_along(smooths)) {
      sm[[j]] <- list()
      tl <- attr(terms.formula(smooths[[j]], specials = specials), "term.labels")
      for(i in seq_along(tl)) {
        if(!grepl("s(", tl[i], fixed = TRUE)) {

        }
      }
    }
  } else {
    if(length(smooths) < 2) {
      term <- deparse(smooths[[1]], backtick = TRUE, width.cutoff = 500)
      if(!grepl("s(", term, fixed = TRUE)) {
        smooths <- s(..., k = k, fx = if(is.null(fx)) FALSE else fx, bs = bs, m = m, xt = xt)
        label <- paste("rs(", term, ")", sep = "")
      } else {
        smooths <- eval(smooths[[1]])
        label <- paste("rs(", smooths$label, ")", sep = "")
        term <- smooths$term
      }
      smooths <- rep(list(smooths), length.out = 2)
    } else {
      if(length(smooths) > 2) stop("more than two terms only possible using formula notation!")
      sm <- list(); term <- label <- NULL
      for(j in seq_along(smooths)) {
        tn <- deparse(smooths[[j]], backtick = TRUE, width.cutoff = 500)
        if(!grepl("s(", tn, fixed = TRUE)) {
          sm[[j]] <- list(
            "term" = tn,
            "param" = TRUE
          )
          term <- c(term, tn)
          label <- c(label, tn)
        } else {
          sm[[j]] <- eval(smooths[[j]])
          term <- c(term, sm[[j]]$term)
          label <- c(label, sm[[j]]$label)
        }
      }
      label <- paste("rs(", paste(label, collapse = ",", sep = ""), ")", sep = "")
      smooths <- sm
    }

    term <- unique(term)
    dim <- min(c(2, length(term)))

    k <- rep(k, length.out = 2)
    if(k[1] != -1) smooths[[1]]$bs.dim <- k[1]
    if(k[2] != -1) smooths[[2]]$bs.dim <- k[2]
    if(!is.null(fx)) {
      fx <- rep(fx, length.out = 2)
      if(!is.null(fx[1])) smooths[[1]]$fixed <- fx[1]
      if(!is.null(fx[2])) smooths[[2]]$fixed <- fx[2]
    }

    rval <- list("smooths" = smooths, "special" = TRUE,
      "term" = term, "label" = label, "dim" = dim,
      "one" = FALSE, "formula" = FALSE, "link" = link)
  }

  class(rval) <- "rs.smooth.spec"
  rval
}

smooth.construct.rs.smooth.spec <- function(object, data, knots) 
{
  lobj <- make.link2(object$link)
  link <- lobj$linkinv
  dlink <- lobj$mu.eta
  edf <- 0; fixed <- NULL
  if(!object$formula) {
    X <- interval <- list(); tau2 <- NULL
    for(j in seq_along(object$smooths)) {
      if(is.null(object$smooths[[j]]$param))
        object$smooths[[j]]$param <- FALSE
      if(object$smooths[[j]]$param) {
        X[[j]] <- data[[object$smooths[[j]]$term]]
        object$smooths[[j]]$fixed <- TRUE
        object$smooths[[j]]$xt <- list("center" = FALSE)
        edf <- edf + ncol(X[[j]])
      } else {
        stj <- object$smooths[[j]]
        acons <- TRUE
        if(!is.null(stj$xt$center))
          acons <- stj$xt$center
        stj$xt$center <- acons
        stj$xt$fixed <- stj$fixed
        fixed <- c(fixed, stj$fixed)
        stj$by.done <- TRUE
        stj <- smoothCon(stj, data, knots, absorb.cons = acons)[[1]]
        stj$class <- class(stj)
        stj$a <- if(is.null(stj$xt$a)) 1e-04 else stj$xt$a
        stj$b <- if(is.null(stj$xt$b)) 1e-04 else stj$xt$b
        if(!stj$fixed) {
          interval[[j]] <- if(is.null(stj$xt$interval)) tau2interval(stj) else stj$xt$interval
        }
        X[[j]] <- stj$X
        object$smooths[[j]] <- stj
        if(!stj$fixed) {
          sp <- if(is.null(stj$sp)) {
            if(is.null(stj$xt$lambda))
              min(c(100, mean(interval[[j]])))
            else
              1 / stj$xt$lambda
          } else stj$sp
          names(sp) <- if(j < 2) "tau2g" else "tau2w"
          tau2 <- c(tau2, sp)
          XX <- crossprod(X[[j]])
          edf <- edf + sum(diag(matrix_inv(XX + 1 / sp * stj$S[[1]]) %*% XX))
        } else {
          edf <- edf + ncol(X[[j]])
        }
        if(acons) edf <- edf - 1
        ## object$smooths[[j]][c("X", "Xf", "rand", "S")] <- NULL
      }
    }

    names(X) <- paste("g", 1:length(X), sep = "")
    X <- as.matrix(as.data.frame(X))

    object$X <- X
    object$p.save <- c("g", if(!is.null(tau2)) "tau2" else NULL)
    object$state <- list()
    object$state$g <- rep(0, ncol(X) - 1)
    object$state$edf <- edf
    object$state$tau2 <- tau2
    object$interval <- interval
    object$fixed <- all(fixed)

    object$fit.fun <- function(X, g, ...) {
      nx <- colnames(X)
      k1 <- length(grep("g1", nx, fixed = TRUE))
      w <- c(1, g[(k1 + 1):(ncol(X) - 1)])
      g <- g[1:k1]
      Z <- X[, -1 * 1:k1, drop = FALSE]
      X <- X[, 1:k1, drop = FALSE]
      f <- drop(X %*% g) / link(drop(Z %*% w))
      f <- f - mean(f, na.rm = TRUE)
      return(drop(f))
    }
  } else {
    stop("formulae no supported yet!")
  }
  object$s.colnames <- c(paste("c", 1:length(object$state$g), sep = ""),
    if(!is.null(tau2)) paste("tau2", 1:length(tau2), sep = "") else NULL)
  object$np <- c(length(object$state$g), length(tau2))

  object$prior <- function(gamma, tau2 = NULL) {
    nx <- colnames(object$X)
    k1 <- length(grep("g1", nx, fixed = TRUE))
    lp <- 0
    g <- gamma[1:k1]
    w <- c(1, gamma[(k1 + 1):(ncol(object$X) - 1)])
    lp <- if(!object$smooths[[1]]$fixed) {
      sp <- object$smooths[[1]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2g"]
        if(is.na(sp))
          sp <- tau2[1]
      }
      lp + log((1 / sp)^(object$smooths[[1]]$rank / 2)) + drop(-0.5 / sp * crossprod(g, object$smooths[[1]]$S[[1]]) %*% g) +
        log((object$smooths[[1]]$b^object$smooths[[1]]$a) / gamma(object$smooths[[1]]$a) * sp^(-object$smooths[[1]]$a - 1) * exp(-object$smooths[[1]]$b / sp))
    } else lp + sum(dnorm(g, sd = 10, log = TRUE))
    lp <- if(!object$smooths[[2]]$fixed) {
      sp <- object$smooths[[2]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2w"]
        if(is.na(sp))
          sp <- tau2[2]
      }
      lp + log((1 / sp)^(object$smooths[[2]]$rank / 2)) + drop(-0.5 / sp * crossprod(w, object$smooths[[2]]$S[[1]]) %*% w) +
        log((object$smooths[[2]]$b^object$smooths[[2]]$a) / gamma(object$smooths[[2]]$a) * sp^(-object$smooths[[2]]$a - 1) * exp(-object$smooths[[2]]$b / sp))
    } else lp + sum(dnorm(w, sd = 10, log = TRUE))
    return(lp)
  }

  object$update <- update_optim

  object$grad <- function(score, gamma, tau2 = NULL, full = TRUE) {
    nx <- colnames(object$X)
    k1 <- length(grep("g1", nx, fixed = TRUE))
    w <- c(1, gamma[(k1 + 1):(ncol(X) - 1)])
    g <- gamma[1:k1]
    Z <- object$X[, -1 * 1:k1, drop = FALSE]
    X <- object$X[, 1:k1, drop = FALSE]
    f1 <- drop(X %*% g)
    f2 <- drop(Z %*% w)
    h <- link(f2)
    Xgrad <- cbind(
      "g" = (1 / h) * X,
      "w" = as.matrix(-1 * ((f1 * Z * dlink(f2)) / (h^2)))[, -1, drop = FALSE]
    )

    grad2 <- NULL
    if(object$smooths[[1]]$fixed) {
      grad <- rep(0, length(g))
    } else {
      sp <- object$smooths[[1]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2g"]
        if(is.na(sp))
          sp <- tau2[1]
      }
      gS <- crossprod(g, object$smooths[[1]]$S[[1]])
      grad <- drop(-0.5 / sp * gS)
      if(full) {
        grad2 <- drop(-object$smooths[[1]]$rank / (2 * sp) - 1 / (2 * sp^2) * gS %*% g + (-object$smooths[[1]]$a - 1) / sp + object$smooths[[1]]$b / (sp^2))
      }
    }

    if(object$smooths[[2]]$fixed) {
      grad <- c(grad, rep(0, length(w) - 1))
    } else {
      sp <- object$smooths[[2]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2w"]
        if(is.na(sp))
          sp <- tau2[2]
      }
      wS <- crossprod(w, object$smooths[[2]]$S[[1]])
      grad <- c(grad, drop(-0.5 / sp * wS)[-1])
      if(full) {
        grad2 <- c(grad2, drop(-object$smooths[[2]]$rank / (2 * sp) - 1 / (2 * sp^2) * wS %*% w + (-object$smooths[[1]]$a - 1) / sp + object$smooths[[2]]$b / (sp^2)))
      }
    }

    gvec <- drop(crossprod(cbind(Xgrad, if(full) {
      matrix(0, nrow = nrow(Xgrad), ncol = length(tau2))
    } else NULL), score)) + c(grad, grad2)

    return(gvec)
  }

  object$propose <- function(x, family, response, eta, id, rho, ...) {
    args <- list(...)
    iter <- args$iter

    if(!is.null(args$no.mcmc)) {
      for(j in seq_along(x$state$tau2)) {
        if(!x$smooths[[j]]$fixed) {
          if(!is.null(x$smooths[[j]]$sp))
            x$state$tau2[j] <- x$smooths[[j]]$sp
          x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
            logPost = logPost5, rho = rho, lower = 0)
        }
      }
    } else {
      if(!is.null(iter)) {
        if(iter %% x$xt$step == 0) {
          if(is.null(x$state$mode)) {
            x$state <- object$update(x, family, response, eta, id, rho, ...)
          } else {
            x$state$g <- x$state$mode$g
            ## x$state$tau2 <- x$state$mode$tau2
          }
        }
      }
    }

    ## Remove fitted values.
    eta[[id]] <- eta[[id]] - x$state$fit

    for(j in seq_along(x$state$g)) {
      x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
        logPost = logPost3, rho = rho)
    }

    ## Setup return state.
    x$state$alpha <- log(1)
    x$state$fit <- x$fit.fun(x$X, x$state$g)

    for(j in seq_along(x$state$tau2)) {
      if(!x$smooths[[j]]$fixed) {
        if(!is.null(x$smooths[[j]]$sp))
          x$state$tau2[j] <- x$smooths[[j]]$sp
        x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
          logPost = logPost5, rho = rho, lower = 0)
      }
    }

    return(x$state)
  }

  object$edf <- function(x, tau2) {
    nx <- colnames(x$X)
    k1 <- length(grep("g1", nx, fixed = TRUE))

    edf1 <- edf2 <- NA
    if(!x$smooths[[1]]$fixed) {
      sp <- x$smooths[[1]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2g"]
        if(is.na(sp))
          sp <- tau2[1]
      }
      XX <- crossprod(x$X[, 1:k1, drop = FALSE])
      P <- matrix_inv(XX + 1 / sp * x$smooths[[1]]$S[[1]])
      if(!inherits(P, "try-error"))
        edf1 <- sum(diag(XX %*% P))
    } else edf1 <- k1
    if(x$smooths[[1]]$xt$center) edf1 <- edf1 - 1

    if(!x$smooths[[2]]$fixed) {
      sp <- x$smooths[[2]]$sp
      if(is.null(sp)) {
        sp <- tau2["tau2w"]
        if(is.na(sp))
          sp <- tau2[2]
      }
      XX <- crossprod(x$X[, -1 * 1:k1, drop = FALSE])
      P <- matrix_inv(XX + 1 / sp * x$smooths[[2]]$S[[1]])
      if(!inherits(P, "try-error"))
        edf2 <- sum(diag(XX %*% P))
    } else edf2 <- ncol(x$X) - k1
    if(x$smooths[[2]]$xt$center) edf2 <- edf2 - 1

    return(edf1 + edf2)
  }

  object$sample <- object$propose
  object$by <- "NA"
  object$by.done <- TRUE

  class(object) <- c("rs.smooth", "mgcv.smooth")
  object
}

Predict.matrix.rs.smooth <- function(object, data, knots)
{
  data <- as.data.frame(data)

  X <- list()
  for(j in seq_along(object$smooths)) {
    if(object$smooths[[j]]$param) {
      X[[j]] <- data[[object$smooths[[j]]$term]]
    } else {
      X[[j]] <- PredictMat(object$smooths[[j]], data)
    }
  }
  names(X) <- paste("g", 1:length(X), sep = "")
  X <- as.matrix(as.data.frame(X))

  X
}


## Penalized harmonic smooth.
smooth.construct.ha.smooth.spec <- function(object, data, knots) 
{
  x <- data[[object$term]]

  freq <- if(is.null(object$xt$frequency)) as.integer(max(x, na.rm = TRUE))
  stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), TRUE))

  if(length(object$p.order) < 2) {
    if(is.na(object$p.order))
      object$p.order <- c(2, 2)
    else
      object$p.order <- c(object$p.order, 2)
  }
  object$p.order[is.na(object$p.order)] <- 2

  order <- object$p.order[1]
  order <- min(freq, order)
  x <- x / freq
  X <- outer(2 * pi * x, 1:order)
  X <- cbind(apply(X, 2, cos), apply(X, 2, sin))
  colnames(X) <- if(order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == freq) X <- X[, -(2 * order)]
  object$X <- X

  gsin1 <- function(x) { cos(2 * pi * order * x) * 2 *pi *order }
  gsin2 <- function(x) { 4 * pi^2 * order^2 * -sin(2 * pi * order * x) }
  gcos1 <- function(x) { -sin(2 * pi * order * x) * 2 * pi * order }
  gcos2 <- function(x) { -4 * pi^2 * order^2 * cos(2 * pi * order * x) }

  if(!object$fixed) {
#    S <- outer(2 * pi * x, 1:order)
#    S <- if(object$p.order[2] < 2) {
#      cbind(apply(S, 2, gcos1), apply(S, 2, gsin1))
#    } else cbind(apply(S, 2, gcos2), apply(S, 2, gsin2))
#    object$S <- list(diag(rep(order:1, 2)))
    K <- t(diff(diag(order))) %*% diff(diag(order))
    K <- rbind(cbind(K, matrix(0, order, order)), cbind(matrix(0, order, order), K))
    object$S <- list(K)
  } else object$S <- list(diag(0, ncol(X)))

  object$frequency <- freq
  object$bs.dim <- ncol(X)
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- ncol(X)
  object$C <- matrix(nrow = 0, ncol = ncol(X))
#  object$no.rescale <- 1
#  object$side.constrain <- FALSE
  class(object) <- "harmon.smooth"
  object
}


Predict.matrix.harmon.smooth <- function(object, data, knots)
{
  x <- data[[object$term]]
  x <- x / object$frequency
  order <- object$p.order[1]
  X <- outer(2 * pi * x, 1:order)
  X <- cbind(apply(X, 2, cos), apply(X, 2, sin))
  colnames(X) <- if (order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == object$frequency) X <- X[, -(2 * order)]
  X
}


## Kriging smooth constructor.
## Evaluate a kriging
## design and penalty matrix.
krDesign1D <- function(z, knots = NULL, rho = NULL,
  phi = NULL, v = NULL, c = NULL, ...)
{
  rho <- if(is.null(rho)) {
    require("geoR")
    geoR::matern
  } else rho
  knots <- if(is.null(knots)) sort(unique(z)) else knots
  v <- if(is.null(v)) 2.5 else v
  c <- if(is.null(c)) {
    optim(1, matern, phi = 1, kappa = v, method = "L-BFGS-B", lower = 1e-10)$par
  } else c
  phi <- if(is.null(phi)) max(abs(diff(range(knots)))) / c else phi
  B <- NULL
  K <- as.matrix(dist(knots, diag = TRUE, upper = TRUE))
  for(j in seq_along(knots)) {
    h <- abs(z - knots[j])
    B <- cbind(B, rho(h, phi, v))
    K[, j] <- rho(K[, j], phi, v)
  }
  return(list("B" = B, "K" = K, "phi" = phi, "v" = v, "c" = c, "knots" = knots))
}

krDesign2D <- function(z1, z2, knots = 10, rho = NULL,
  phi = NULL, v = NULL, c = NULL, psi = NULL, delta = 1,
  isotropic = TRUE, ...)
{
  rho <- if(is.null(rho)) {
    require("geoR")
    geoR::matern
  } else rho
  if(is.null(psi)) psi <- 1
  if(is.null(delta)) delta <- 1
  if(is.null(isotropic)) isotropic <- TRUE
  if(is.null(knots)) knots <- min(c(10, nrow(unique(cbind(z1, z2)))), na.rm = TRUE)
  knots <- if(length(knots) < 2) {
    if(knots == length(z1)) {
      unique(cbind(z1, z2))
    } else {
      require("fields")
      fields::cover.design(R = unique(cbind(z1, z2)), nd = knots)
    }
  } else knots
  v <- if(is.null(v)) 2.5 else v
  c <- if(is.null(c)) {
    optim(1, rho, phi = 1, kappa = v,
      method = "L-BFGS-B", lower = 1e-10)$par
  } else c
  z <- cbind(z1, z2)
  if(class(knots) == "spatial.design")
    knots <- knots[, 1:2]
  if(!is.matrix(knots))
    knots <- matrix(knots, ncol = 2)
  nk <- nrow(knots)
  phi <- if(is.null(phi)) {
    max(abs(diff(range(knots)))) / c
  } else phi
  if(phi == 0)
    phi <- max(abs(fields::rdist(z1, z2))) / c
  K <- rho(fields::rdist(knots, knots), phi, v)
  if(isotropic) {
    B <- NULL
    for(j in 1:nk) {
      kn <- matrix(knots[j, ], nrow = 1, ncol = 2)
	    h <- fields::rdist(z, kn)
		  B <- cbind(B, rho(h, phi, v))
    }
  } else {
    B <- matrix(0, nrow(z), nk)
    R <- matrix(c(cos(psi), -1 * sin(psi),
      sin(psi), cos(psi)), 2, 2)
    D <- matrix(c(delta^(-1), 0, 0, 1), 2, 2)
    for(i in 1:nrow(z)) {
      for(j in 1:nk) {
        kn <- matrix(knots[j, ], nrow = 1, ncol = 2)
        h <- as.numeric(z[i, ] - kn)
        h <- drop(sqrt(t(h) %*% t(R) %*% D %*% R %*% h))
        B[i, j] <- rho(h, phi, v)
      }
    }
  }

  return(list("B" = B, "K" = K, "knots" = knots,
    "phi" = phi, "v" = v, "c" = c, "psi" = psi,
    "delta" = delta))
}


## Kriging smooth constructor functions.
smooth.construct.kr.smooth.spec <- function(object, data, knots)
{
  if(object$dim > 2) stop("more than 2 covariates not supported using kriging terms!")
  if(object$bs.dim < 0) object$bs.dim <- 10
  if(object$dim < 2) {
    k <- knots[[object$term]]
    x <- data[[object$term]]
    if(is.null(k))
      k <- seq(min(x), max(x), length = object$bs.dim)
    D <- krDesign1D(x, knots = k, rho = object$xt$rho,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c)
  } else {
    knots <- if(is.null(object$xt$knots)) object$bs.dim else object$xt$knots
    D <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
      knots = knots,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
      psi = object$xt$psi, delta = object$xt$delta,
      isotropic = object$xt$isotropic)
  }

  X <- D$B
  object$X <- X
  object$S <- list(D$K)
  object$rank <- qr(D$K)$rank
  object$knots <- D$knots
  object$null.space.dim <- ncol(D$K)
 
  class(object) <- "kriging.smooth"
  object
}

## Predict function for the new kriging smooth.
Predict.matrix.kriging.smooth <- function(object, data)
{
  if(object$dim < 2) {
    X <- krDesign1D(data[[object$term]], knots = object$knots, rho = object$xt$rho,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c)$B
  } else {
    X <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
      knots = object$knots,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
      psi = object$xt$psi, delta = object$xt$delta,
      isotropic = object$xt$isotropic)$B
  }
  X
}


## Smooth constructor for lag function.
## (C) Viola Obermeier.
smooth.construct.fdl.smooth.spec <- function(object, data, knots)
{
  ## Modify object so that it's fitted as a p-spline signal regression term.
  object$bs <- "ps"
  object <- smooth.construct.ps.smooth.spec(object, data, knots)

  if(!is.null(object$xt$fullrankpen) && object$xt$fullrankpen){
    ## Add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions.
    ## With same variance as difference penalty: penalty = lambda * coef' (DiffPen + RidgePen) coef.
    object$S[[1]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] <- object$S[[1]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] + 1
    object$rank <- min(object$bs.dim, object$rank + object$m[1]+2)
  }
  if(!is.null(object$xt$ridge) && object$xt$ridge){
    ## Add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions
    ## penalty = coef' (lambda_1*DiffPen + lambda_2*RidgePen) coef.
    object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
    object$S[[2]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] <- 1
    object$rank <- c(object$rank, object$m[1]+2)
  }
  if(!is.null(object$xt$constrain) && object$xt$constrain){
    ## Optionally one can constrain the last lag coefficient to be zero,
    ## not recommended as we favor a soft, data-driven shrinkage to a hard constraint!
    ## Constrain to end in zero (i.e (X%*%coefficients)[1] == 0).
    ## --> Constraint matric C = X[1,]
    object$C <- matrix(object$X[1,],nrow=1)
    object$C <- structure(object$C, always.apply=TRUE)
  }

  return(object)
}


## Plotting method for "bamlss" objects.
plot.bamlss <- function(x, model = NULL, term = NULL, which = "samples",
  ask = dev.interactive(), ...)
{
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  ## What should be plotted?
  which.match <- c("effects", "samples", "hist-resid", "qq-resid",
    "scatter-resid", "scale-resid", "max-acf", "param-samples")
  if(!is.character(which)) {
    if(any(which > 8L))
      which <- which[which <= 8L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if(which == "samples") {
    par <- if(is.null(x$parameters)) NULL else unlist(x$parameters)
    samps <- samples(x, model = model, term = term, drop = TRUE)
    snames <- colnames(samps)
    snames <- snames[!grepl(".p.edf", snames, fixed = TRUE)]
    samps <- samps[, snames, drop = FALSE]
    np <- ncol(samps)
    par(mfrow = if(np <= 4) c(np, 2) else c(4, 2))
    devAskNewPage(ask)
    tx <- as.vector(time(samps))
    for(j in 1:np) {
      al <- if(grepl("logLik", snames[j], fixed = TRUE)) {
        x$logLik
      } else {
        if(!is.null(par)) {
          if(snames[j] %in% names(par))
            par[snames[j]]
        } else NULL
      }
      lim <- range(c(al, samps[, j]), na.rm = TRUE)
      traceplot(samps[, j, drop = FALSE], main = paste("Trace of", snames[j]), ylim = lim)
      lines(lowess(tx, samps[, j]), col = "red")
      if(!is.null(al))
        abline(h = al, col = "blue")
      if(snames[j] %in% c(names(par), "logLik"))
        abline(h = mean(samps[, j], na.rm = TRUE), col = "green")
      densplot(samps[, j, drop = FALSE], main = paste("Density of", snames[j]), xlim = lim)
      if(!is.null(al))
        abline(v = al, col = "blue")
      if(snames[j] %in% c(names(par), "logLik"))
        abline(v = mean(samps[, j], na.rm = TRUE), col = "green")
    }
  }

  if(which == "effects") {
    if(is.null(x$results)) {
      plot(results.bamlss.default(x), model = model, term = term, ...)
    } else plot(x$results, model = model, term = term, ...)
  }

  return(invisible(NULL))
}


plot.bamlss.results <- function(x, model = NULL, term = NULL,
  ask = FALSE, scale = 1, spar = TRUE, ...)
{
  args <- list(...)
  cx <- class(x)

  if(is.null(args$do_par) & spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  if(!is.null(model)) {
    if(!is.character(model)) {
      if(any(model < 0 | model > length(x)))
        stop("model specified wrong!")
      model <- names(x)[model]
    } else {
      i <- grep2(model, names(x))
      if(!length(i))
        stop("model specified wrong!")
      model <- names(x)[i]
    }
    x <- x[model]
  }

  if(FALSE) {
    ## What should be plotted?
    which.match <- which <- "effects"
    if(!is.character(which)) {
      if(any(which > 8L))
        which <- which[which <= 8L]
      which <- which.match[which]
    } else which <- which.match[pmatch(tolower(which), which.match)]
    if(length(which) > length(which.match) || !any(which %in% which.match))
      stop("argument which is specified wrong!")

    args2 <- args
    args2$object <- x
    res0 <- do.call("residuals.bamlss", delete.args("residuals.bamlss", args2))
    ny <- if(is.null(dim(res0))) 1 else ncol(res0)
    if(is.null(args$do_par) & spar) {
      if(!ask) {
        par(mfrow = n2mfrow(length(which) * ny))
      } else par(ask = ask)
    }
    if(any(which %in% c("scatter-resid", "scale-resid"))) {
      fit0 <- fitted.bamlss(x, type = "parameter", samples = TRUE,
        model = if(ny < 2) 1 else NULL, nsamps = args$nsamps)
    }
    rtype <- args$type
    if(is.null(rtype)) rtype <- "quantile"
    if(rtype == "quantile2") rtype <- "quantile"
    if(rtype == "ordinary2") rtype <- "ordinary"
    for(j in 1:ny) {
      res <- if(ny > 1) res0[, j] else res0
      dropi <- !(res %in% c(Inf, -Inf)) & !is.na(res)
      res <- res[dropi]
      if(any(which %in% c("scatter-resid", "scale-resid"))) {
        fit <- if(ny < 2) {
          if(is.list(fit0)) fit0[[1]] else fit0
        } else fit0[[j]]
      }
      for(w in which) {
        args2 <- args
        if(w == "hist-resid") {
          rdens <- density(res)
          rh <- hist(res, plot = FALSE)
          args2$ylim <- c(0, max(c(rh$density, rdens$y)))
          args2$freq <- FALSE
          args2$x <- res
          args2 <- delete.args("hist.default", args2, package = "graphics")
          if(is.null(args$xlab)) args2$xlab <- paste(if(rtype == "quantile") {
              "Quantile"
            } else "Ordinary", "residuals")
          if(is.null(args$ylab)) args2$ylab <- "Density"
          if(is.null(args$main)) {
            args2$main <- "Histogramm and density"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(eval(parse(text = "graphics::hist.default")), args2))
          if(!inherits(ok, "try-error"))
            lines(rdens)
          box()
        }
        if(w == "qq-resid") {
          args2$y <- if(rtype == "quantile") (res) else (res - mean(res)) / sd(res)
          args2 <- delete.args("qqnorm.default", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$main)) {
            args2$main <- "Normal Q-Q Plot"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(qqnorm, args2))
          if(!inherits(ok, "try-error"))
  		      if(rtype == "quantile") abline(0,1) else qqline(args2$y)
        }
        if(w == "scatter-resid") {
          args2$x <- fit[dropi]
          args2$y <- res
          args2 <- delete.args("scatter.smooth", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$xlab)) args2$xlab <- "Fitted values"
          if(is.null(args$xlab)) args2$ylab <- paste(if(rtype == "quantile") {
              "Quantile"
            } else "Ordinary", "residuals")
          if(is.null(args$xlab)) {
            args2$main <- "Fitted values vs. residuals"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(scatter.smooth, args2))
          if(!inherits(ok, "try-error"))
            abline(h = 0, lty = 2)
        }
        if(w == "scale-resid") {
          args2$x <- fit[dropi]
          args2$y <- sqrt(abs((res - mean(res)) / sd(res)))
          args2 <- delete.args("scatter.smooth", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$xlab)) args2$xlab <- "Fitted values"
          if(is.null(args$ylab)) args2$ylab <- expression(sqrt(abs("Standardized residuals")))
          if(is.null(args$main)) {
            args2$main <- "Scale-location"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          try(do.call(scatter.smooth, args2))
        }
      }
    }
  } else {
    ## Get number of plots.
    get_k_n <- function(x) {
      kn <- c(0, length(x))
      ne <- pterms <- list()
      for(i in 1:kn[2]) {
        if(!any(c("s.effects", "p.effects") %in% names(x[[i]]))) {
          kn <- kn + get_k_n(x[[i]])
        } else {
          ne[[i]] <- if(!is.null(names(x[[i]]$s.effects))) names(x[[i]]$s.effects) else NA
          if(is.null(term))
            pterms[[i]] <- 1:length(ne[[i]])
          else {
            if(is.character(term)) {
              tterm <- NULL
              for(j in term)
                tterm <- c(tterm, grep(j, ne[[i]], fixed = TRUE))
              pterms[[i]] <- if(length(tterm)) tterm else NA
            } else pterms[[i]] <- term[term <= length(ne[[i]])]
          }
          if(!is.null(x[[i]]$s.effects)) {
            kn[1] <- kn[1] + length(na.omit(pterms[[i]]))
          }
        }
      }
      kn
    }

    if(any(c("s.effects", "p.effects") %in% names(x)))
      x <- list(x)

    kn <- get_k_n(x)

    if(kn[1] < 1) on.exit(warning("no terms to plot!"), add = TRUE)

    if(is.null(args$do_par) & spar) {
      if(!ask) {
        if("cbamlss" %in% cx) {
          par(mfrow = c(length(x), kn[1] / length(x)))
        } else par(mfrow = n2mfrow(kn[1]))
      } else par(ask = ask)
    }

    mmain <- if(any(c("h1", "Chain_1") %in% (nx <- names(x)))) TRUE else FALSE
    main <- args$main
    if((is.null(args$main) & mmain) | !is.null(args$mmain)) {
      main <- if(!is.null(args$main)) paste(args$main, nx, sep = "-") else nx
      args$mmain <- TRUE
    }
    if(!is.null(main)) main <- rep(main, length.out = length(x))
    for(i in seq_along(x)) {
      args[c("x", "term", "ask", "scale")] <- list(x[i], term, ask, scale)
      args$main <- if(!is.null(main)) main[i] else NULL
      if(!any(c("s.effects", "p.effects") %in% names(x[[i]]))) {
        args$do_par <- FALSE
        do.call("plot.bamlss", args)
      } else {
        args$mmain <- NULL
        do.call(".plot.bamlss.results", args)
      }
    }
  }

  invisible(NULL)
}

.plot.bamlss.results <- function(x, model = NULL, term = NULL,
  ask = FALSE, scale = 1, spar = TRUE, ...)
{
  n <- length(x)
  args <- list(...)

  ## Effect plotting.
  k <- 0; ylim <- NULL
  ylim <- args$ylim
  args$residuals <- if(is.null(args$residuals)) FALSE else args$residuals
  if(!is.null(args$ylim))
    scale <- 0
  ne <- pterms <- list()
  for(i in 1:n) {
    ne[[i]] <- if(!is.null(names(x[[i]]$s.effects))) names(x[[i]]$s.effects) else NA
    if(is.null(term))
      pterms[[i]] <- 1:length(ne[[i]])
    else {
      if(is.character(term)) {
        tterm <- NULL
        for(j in term)
          tterm <- c(tterm, grep(j, ne[[i]], fixed = TRUE))
        pterms[[i]] <- if(length(tterm)) tterm else NA
      } else pterms[[i]] <- term[term <= length(ne[[i]])]
    }
  }
  for(i in 1:n) {
    if(!is.null(x[[i]]$s.effects) & length(na.omit(pterms[[i]]))) {
      k <- k + length(na.omit(pterms[[i]]))
      if(scale > 0) {
        term <- term[1:length(x[[i]]$s.effects)]
        for(e in pterms[[i]]) {
          et <- x[[i]]$s.effects[[e]]
          de <- attr(et, "specs")$dim + 1
          ylim <- c(ylim, range(et[, de:ncol(et)], na.rm = TRUE))
          if(args$residuals) {
            if(!is.null(attr(et, "partial.resids"))) {
              res <- attr(et, "partial.resids")
              ylim <- c(ylim, range(res[, de:ncol(res)], na.rm = TRUE))
            }
          }
        }
      }
    }
  }
  if(k < 1) return(NULL)
  if(scale > 0)
    ylim <- range(ylim, na.rm = TRUE)
  args$residuals <- NULL
  for(i in 1:n) {
    if(!is.null(x[[i]]$s.effects)) {
      for(e in pterms[[i]]) {
        lim <- c("ylim", "zlim")[(attr(x[[i]]$s.effects[[e]], "specs")$dim > 1) * 1 + 1]
        setlim <- FALSE
        if(!is.null(ylim) & is.null(args[[lim]])) {
          args[[lim]]<- ylim
          setlim <- TRUE
        }
        args$x <- x[[i]]$s.effects[[e]]
        do.call("plot.bamlss.effect", args)
        if(setlim) args[[lim]] <- NULL
      }
    }
  }

  return(invisible(NULL))
}


## Generic plotting method for model terms.
plot.bamlss.effect <- function(x, ...) {
  UseMethod("plot.bamlss.effect")
}


## Default model term plotting method.
plot.bamlss.effect.default <- function(x, ...) {
  args <- list(...)

  if(attr(x, "specs")$dim > 1 & inherits(x, "rs.smooth")) {
    if(identical(x[, 1], x[, 2])) {
      cn <- colnames(x)[-2]
      xattr <- attributes(x)
      xattr$specs$dim <- 1
      x <- x[, -2, drop = FALSE]
      xattr$names <- colnames(x) <- cn
      cn <- colnames(xattr$partial.resids)[-2]
      xattr$partial.resids <- xattr$partial.resids[, -2, drop = FALSE]
      colnames(xattr$partial.resids) <- cn
      mostattributes(x) <- xattr
    }
  }

  if(length(terms <- attr(x, "specs")$term) > 2) {
    if(is.null(args$view)) args$view <- terms[1:2]
    args$view <- rep(args$view, length.out = 2)
    cn <- colnames(x)
    td <- terms[!(i <- (terms %in% args$view))]
    xattr <- attributes(x)
    xattr$specs$dim <- 2
    args$cond <- if(is.null(args$cond)) mean(x[, td], na.rm = TRUE) else args$cond
    args$cond <- x[which.min(abs(x[, td] - args$cond)), td]
    samps <- attr(x, "samples")
    specs <- attr(x, "specs")
    newdata <- x[, 1:length(attr(x, "specs")$term), drop = FALSE]
    newdata[[td]] <- args$cond

    X <- if(inherits(specs, "mgcv.smooth")) {
      PredictMat(specs, newdata)
    } else {
      if(!is.null(specs$basis)) {
        stopifnot(is.function(specs$basis))
        if(specs$by != "NA") {  ## ATTENTION: by variables with basis()!
          if(!(specs$by %in% names(newdata)))
            stop("cannot find by variable ", specs$by, " in newdata!")
          if(!is.factor(unlist(newdata[specs$by]))) {
            as.numeric(unlist(newdata[specs$by])) * specs$basis(newdata[specs$term])
          } else specs$basis(newdata[specs$term])
        } else specs$basis(newdata[specs$term])
      } else stop(paste("cannot compute design matrix for term ", specs$label, "!", sep = ""))
    }

    FUN <- function(x) {
      rval <- as.numeric(quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      names(rval) <- c("2.5%", "50%", "97.5%")
      rval
    }
    fit <- apply(samps, 1, function(g) {
      names(g) <- NULL
      specs$fit.fun(X, g, expand = FALSE)
    })
    fit <- t(apply(fit, 1, FUN))

    x <- cbind(x[, terms[i]], fit)
    xattr$names <- colnames(x)
    ##xattr$partial.resids <- xattr$partial.resids[, -td, drop = FALSE]
    mostattributes(x) <- xattr
    attr(x, "specs")$dim <- 2
  }

  args$x <- x

  lim <- c("ylim", "zlim")[(attr(x, "specs")$dim > 1) * 1 + 1]
  limNULL <- FALSE
  if(is.null(args[[lim]])) {
    limNULL <- TRUE
    args[[lim]] <- range(x[, c("2.5%", "97.5%")], na.rm = TRUE)
    if(!is.null(args$residuals)) {
      if(args$residuals & !is.null(attr(x, "partial.resids")))
        args[[lim]] <- range(c(args[[lim]], attr(x, "partial.resids")[, -1]), na.rm = TRUE)
    }
  }
  if(!is.null(args$shift))
    args[[lim]] <- args[[lim]] + args$shift
  if(attr(x, "specs")$dim < 2) {
    if(is.null(args$fill.select))
      args$fill.select <- c(0, 1, 0, 1)
    if(is.null(args$lty) & is.null(args$map))
      args$lty <- c(2, 1, 2)
    if(is.null(args$col.lines))
      args$col.lines <- c(NA, "black", NA)
    if(inherits(x, "random.effect") | inherits(x, "re.smooth.spec") |
       inherits(x, "mrf.smooth.spec") | inherits(x, "mrf.smooth") | is.factor(x[[1]])) {
      if(if(!is.null(args$density)) args$density else FALSE) {
        args$density <- NULL
        if(is.null(args$main))
          args$main <- attr(x, "specs")$label
        args$x <- density(x[, "50%"], na.rm = TRUE)
        if(!limNULL)
          args$xlim <- args$ylim
        do.call("plot", delete.args(plot.density, args, c("main", "xlim")))
      } else {
        if(!is.null(args$map)) {
          args$x <- x[, grepl("50%", colnames(x), fixed = TRUE)]
          args$id <- as.character(x[, 1])
          args$xlim <- args$ylim <- NULL
          do.call("plotmap", delete.args("plotmap", args,
            not = c("border", "lwd", "lty", names(formals("colorlegend")), "main")))
        } else {
          if(is.null(args$ylab))
            args$ylab <- attr(x, "specs")$label
          do.call("plotblock", delete.args("plotblock", args,
            c("xlim", "ylim", "pch", "main", "xlab", "ylab", "lwd", "axes", "add")))
        }
      }
    } else {
      do.call("plot2d", delete.args("plot2d", args,
        c("xlim", "ylim", "pch", "main", "xlab", "ylab", "lwd", "axes", "add")))
    }
  } else {
    if(is.null(args$c.select))
      args$c.select <- grep("50%", colnames(x), fixed = TRUE)
    if(!is.null(args$slice)) {
      do.call("sliceplot", delete.args("sliceplot", args,
        c("xlim", "ylim", "zlim", "main", "xlab", "ylab", "col", "lwd", "lty")))
    } else {
      if(inherits(x, "random.effect")) {
        do.call("bamlss_random_plot", args)
      } else {
        specs <- attr(x, "specs")
        isf <- sapply(args$x[, specs$term], is.factor)
        if(any(isf)) {
          do.call("bamlss_factor2d_plot", args)
        } else {
          do.call("plot3d", delete.args("plot3d", args,
            c("xlim", "ylim", "zlim", "pch", "main", "xlab", "ylab", "ticktype",
            "zlab", "phi", "theta", "r", "d", "scale", "range", "lrange", "pos", "image.map",
            "symmetric", "border", "lwd")))
        }
      }
    }
  }
}


bamlss_random_plot <- function(x, ...)
{
  term <- attr(x, "specs")$term
  cn <- colnames(x)
  isf <- sapply(x[, term], is.factor)
  plot(x[, "50%"] ~ x[, term[!isf]], type = "n", xlab = term[!isf], ylab = attr(x, "specs")$label)
  id <- x[, term[isf]]
  col <- rainbow_hcl(length(unique(id)))
  ii <- 1
  for(j in unique(id)) {
    d <- subset(x, x[, term[isf]] == j)
    i <- order(d[, term[!isf]])
    lines(d[i, "50%"] ~ d[i, term[!isf]], col = col[ii])
    ii <- ii + 1
  }
  return(invisible(NULL))   
}

bamlss_factor2d_plot <- function(x, ids = NULL, add = FALSE, ...)
{
  args <- list(...)
  y <- args$response
  specs <- attr(x, "specs")
  if(is.null(specs)) {
    specs <- list("term" = colnames(x)[1:2],
      label = paste("f(", colnames(x)[1], ",", colnames(x)[2], ")", sep = ""))
  }
  isf <- sapply(x[, specs$term], is.factor)
  xd <- x[, specs$term]
  fx <- unlist(x[, grepl("50", colnames(x), fixed = TRUE)])
  isf <- isf[1:length(specs$term)]
  xlab <- if(is.null(args$xlab)) colnames(xd)[!isf] else args$xlab
  ylab <- if(is.null(args$ylab)) specs$label else args$ylab
  id <- xd[, isf]
  xd <- xd[, !isf]
  if(!is.null(ids)) {
    if(!is.character(ids))
      ids <- levels(id)[as.integer(ids)]
    i <- id %in% ids
    id <- droplevels(id[i])
    xd <- xd[i]
    fx <- fx[i]
    if(!is.null(y))
      y <- y[i]
  }
  args$ylim <- args$zlim
  xlim <- if(is.null(args$xlim)) range(xd) else args$xlim
  ylim <- if(is.null(args$ylim)) range(fx) else args$ylim
  if(!add) {
    plot(1, 1, type = "n",
      xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, main = args$main)
  }
  col <- if(is.null(args$col)) rainbow_hcl(nlevels(id)) else args$col
  if(is.function(col))
     col <- col(nlevels(id))
  lwd <- if(is.null(args$lwd)) 1 else args$lwd
  lty <- if(is.null(args$lty)) 1 else args$lty
  col <- rep(col, length.out = nlevels(id))
  lwd <- rep(lwd, length.out = nlevels(id))
  lty <- rep(lty, length.out = nlevels(id))
  i <- 1
  for(j in levels(id)) {
    fid <- fx[id == j]
    tid <- xd[id == j]
    lines(fid ~ tid, col = col[i], lwd = lwd[i], lty = lty[i])
    if(!is.null(y))
      points(tid, y[id == j], col = col[i])
    i <- i + 1
  }
  rug <- if(is.null(args$rug)) TRUE else args$rug
  if(rug) {
    jitter <- if(is.null(args$jitter)) TRUE else args$jitter
    if(jitter)
      xd <- jitter(xd)
    rug(xd, col = args$rug.col)
  }
  return(invisible(NULL))   
}


## Other helping functions.
delete.args <- function(fun = NULL, args = NULL, not = NULL, package = NULL)
{
  if(is.character(fun) & !is.null(package))
    fun <- eval(parse(text = paste(package, paste(rep(":", 3), collapse = ""), fun, sep = "")))
  nf <- names(formals(fun))
  na <- names(args)
  for(elmt in na)
    if(!elmt %in% nf) {
      if(!is.null(not)) {
        if(!elmt %in% not)
          args[elmt] <- NULL
      } else args[elmt] <- NULL
    }

  return(args)
}

delete.NULLs <- function(x.list) 
{
  x.list[unlist(lapply(x.list, length) != 0)]
}


## Model summary functions.
summary.bamlss <- function(object, model = NULL, ...)
{
  call <- object$call
  family <- object$family
  object <- model.terms(object, model)
  rval <- list()
  n <- length(object)
  for(i in 1:n) {
    if(!any(c("p.effects", "s.effects.resmat") %in% names(object[[i]]))) {
      rval[[i]] <- summary.bamlss(object[[i]])
      attr(rval[[i]], "hlevel") <- TRUE
    } else {
      for(j in c("p.effects", "s.effects.resmat")) {
        if(!is.null(object[[i]][[j]]))
          attr(object[[i]][[j]], "samples") <- NULL
      }
      rval[[i]] <- with(object[[i]],
        c(list("p.effects" = p.effects,
          "s.effects.resmat" = s.effects.resmat),
          model)
      )
    }
  }
  if(n < 2)
    rval <- rval[[1]]
  else
    names(rval) <- names(object)
  attr(rval, "n") <- n
  attr(rval, "call") <- call
  attr(rval, "family") <- family
  class(rval) <- "summary.bamlss"
  rval
}

print.summary.bamlss <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  on.exit(return(invisible(x)))
  h0 <- !is.null(attr(x, "hlevel"))

  dic_out <- is.null(list(...)$dic_out)

  print_dic_pd <- function(x, ok = TRUE) {
    if(is.list(x)) {
      if(any(c("DIC", "pd") %in% names(x[[1]])))
        x <- x[[1]]
    }
    dp <- FALSE
    if(!is.null(x$IC) & !is.null(x$edf)) {
      dp <- TRUE
      if(ok) cat("\n---") else cat("---")
      cat(if(is.null(names(x$IC))) "\nlog Lik. =" else paste("\n", names(x$IC), " =", sep = ""),
        if(is.na(x$IC)) "NA" else {
            formatC(x$IC, digits = digits, flag = "-")
          }, "edf =", if(is.na(x$edf)) "NA" else {
            formatC(x$edf, digits = digits, flag = "-")
      })
    }
    if(!is.null(x$DIC) & !is.null(x$pd)) {
      dp <- TRUE
      if(ok) cat("\n---") else cat("---")
      cat("\nDIC =", if(is.na(x$DIC)) "NA" else {
            formatC(x$DIC, digits = digits, flag = "-")
          }, "pd =", if(is.na(x$pd)) "NA" else {
            formatC(x$pd, digits = digits, flag = "-")
          })
    }
    if(!is.null(x$N)) {
      dp <- TRUE
      cat(" N =", if(is.na(x$N)) "NA" else formatC(x$N, digits = digits, flag = "-"))
    }
    if(dp) cat("\n\n")
  }

  call <- attr(x, "call")
  family <- attr(x, "family")
  n <- attr(x, "n")
  nx <- NULL
  if(n < 2)
    x <- list(x)
  else
    nx <- names(x)
  cat("\n")
  if(!is.null(call)) {
    cat("Call:\n"); print(call)
    cat("\n")
  }
  if(!is.null(family)) {
    print(if(is.function(family)) family() else family)
    cat("---\n\n")
  }

  for(i in 1:n) {
    h1 <- !is.null(attr(x[[i]], "hlevel"))
    if(!is.null(nx)) {
      cat("Results for ", nx[i], ":\n", sep = "")
      if(h1) cat("---") else cat("---\n")
    }
    if(h1) {
      print.summary.bamlss(x[[i]], digits = digits, dic_out = FALSE, ...)
      if(i == n & dic_out)
        print_dic_pd(x[[i]][[1]], ok = FALSE)
    } else {
      cat("Formula:\n")
      environment(x[[i]]$formula) <- .GlobalEnv
      attr(x[[i]]$formula, "name") <- NULL
      print(x[[i]]$formula)
      if(length(x[[i]]$p.effects) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x[[i]]$p.effects, digits = digits, na.print = "NA", ...)
      }
      if(length(x[[i]]$s.effects.resmat) > 0) {
        if(is.null(attr(x[[i]]$s.effects.resmat, "edf")))
          cat("\nSmooth terms (variances):\n")
        else
          cat("\nSmooth terms (edf):\n")
        printCoefmat(x[[i]]$effects, digits = digits, na.print = "NA", ...)
      }
      if(i == n & !h0) {
        print_dic_pd(x[[1]])
      } else cat("\n")
    }
  }
}


## Simple "bamlss" print method.
print.bamlss <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  print(family(x))
  cat("*---\n")
  print(formula(x))
  if(any(c("logLik", "logPost", "IC", "edf") %in% names(x))) {
    cat("*---\n")
    sep <- ""
    if(!is.null(x$logLik)) {
      cat("logLik =", fmt(x$logLik, width = digits, digits = digits))
      sep <- " "
    }
    if(!is.null(x$logPost)) {
      cat(sep, "logPost =", fmt(x$logPost, width = digits, digits = digits))
      sep <- " "
    }
    if(!is.null(x$IC)) {
      if(!is.null(names(x$IC))) {
        cat(sep, names(x$IC), "=", fmt(x$IC, width = digits, digits = digits))
        sep <- " "
      }
    }
    if(!is.null(x$edf))
      cat(sep, "edf =", fmt(x$edf, width = 4, digits = digits))
    cat("\n")
  }
  return(invisible(NULL))
}


## More extractor functions.
DIC.bamlss <- function(object, ..., samples = TRUE, nsamps = NULL)
{
  object <- c(object, ...)
  rval <- NULL

  if(!samples) {
    for(i in 1:length(object)) {
      xs <- summary(object[[i]])
      n <- attr(xs, "n")
      if(n < 2)
        xs <- list(xs)
      rval <- rbind(rval, data.frame(
        "DIC" = if(is.null(xs[[n]]$DIC)) NA else xs[[n]]$DIC,
        "pd" = if(is.null(xs[[n]]$DIC)) NA else xs[[n]]$pd
      ))
    }
  } else {
    for(i in 1:length(object)) {
      family <- attr(object[[i]], "family")
      if(is.null(family$d)) stop("no d() function available in model family object!")
      y <- model.response2(object[[i]])
      d0 <- -2 * sum(family$d(y, family$map2par(fitted.bamlss(object[[i]], type = "link")), log = TRUE), na.rm = TRUE)
      eta <- fitted.bamlss(object[[i]], type = "link",
        samples = TRUE, nsamps = nsamps,
        FUN = function(x) { x })
      iter <- ncol(eta[[1]])
      d1 <- NULL
      for(j in 1:iter) {
        teta <- NULL
        for(ii in 1:length(eta))
          teta <- cbind(teta, eta[[ii]][, j])
        teta <- as.data.frame(teta)
        names(teta) <- names(eta)
        d1 <- c(d1, -2 * sum(family$d(y, family$map2par(teta), log = TRUE), na.rm = TRUE))
      }
      md1 <- mean(d1)
      pd <- md1 - d0
      dic <- md1 + pd
      rval <- rbind(rval, data.frame(
        "DIC" = dic,
        "pd" = pd
      ))
    }
  }

  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}


logLik.bamlss <- function(object, ..., type = 1, nsamps = NULL, FUN = NULL)
{
  object <- c(object, ...)
  rval <- if(type %in% c(1, 3)) NULL else list()

  if(type == 3) {
    for(i in 1:length(object)) {
      xs <- summary(object[[i]])
      n <- attr(xs, "n")
      if(n < 2)
        xs <- list(xs)
      rval <- rbind(rval, data.frame(
        "logLik" = if(is.null(xs[[n]]$IC)) NA else xs[[n]]$IC,
        "edf" = if(is.null(xs[[n]]$edf)) NA else xs[[n]]$edf
      ))
    }
    names(rval)[1] <- names(xs[[n]]$IC)[1]
  }
  if(type %in% c(1, 2)) {
    for(i in 1:length(object)) {
      family <- attr(object[[i]], "family")
      if(is.null(family$d)) stop("no d() function available in model family object!")
      y <- model.response2(object[[i]])
      if(type == 2) {
        eta <- fitted.bamlss(object[[i]], type = "link",
          samples = TRUE, nsamps = nsamps,
          FUN = function(x) { x })
        iter <- min(sapply(eta, ncol))
        ll <- NULL
        for(j in 1:iter) {
          teta <- NULL
          for(ii in 1:length(eta))
            teta <- cbind(teta, eta[[ii]][, j])
          teta <- as.data.frame(teta)
          names(teta) <- names(eta)
          ll <- c(ll, sum(family$d(y, family$map2par(teta), log = TRUE), na.rm = TRUE))
        }
        if(is.null(FUN)) FUN <- function(x) { x }
        rval[[i]] <- FUN(ll)
      } else {
        eta <- fitted.bamlss(object[[i]], type = "link", samples = TRUE,
          nsamps = nsamps, FUN = function(x) { mean(x, na.rm = TRUE) })
        ll <- sum(family$d(y, family$map2par(eta), log = TRUE), na.rm = TRUE)
        rval <- rbind(rval, data.frame(
          "logLik" = ll,
          "edf" = NA
        ))
      }
    }
  }

  Call <- match.call()
  if(type %in% c(1, 3)) {
    row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L])[1:nrow(rval)] else ""
  } else {
    if(length(rval) < 2) {
      rval <- rval[[1]]
    } else {
      names(rval) <- as.character(Call[-1L])[1:length(rval)]
    }
  }

  rval
}


edf <- function (object, ...) { 
  UseMethod("edf")
}

edf.bamlss <- function(object, FUN = mean, ...)
{
  family <- attr(object, "family")
  nx <- family$names
  edf0 <- 0; edf1 <- NULL
  for(j in 1:length(nx)) {
    if(!is.null(object[[nx[j]]]$p.effects))
      edf0 <- edf0 + nrow(object[[nx[j]]]$p.effects)
    if(!is.null(object[[nx[j]]]$effects)) {
      for(sj in 1:length(object[[nx[j]]]$effects)) {
        if(!is.null(attr(object[[nx[j]]]$effects[[sj]], "samples.edf"))) {
          edf1 <- cbind(edf1, attr(object[[nx[j]]]$effects[[sj]], "samples.edf"))
        } else {
          stop(paste("equivalent degrees of freedom not available for term ",
            object[[nx[j]]]$effects[[sj]]$label, "!", sep = ""))
        }
      }
    }
  }
  rval <- apply(edf1, 1, sum, na.rm = TRUE) + edf0
  FUN(rval, ...)
}


## Extract model formulas.
formula.bamlss.frame <- formula.bamlss <- function(x, model = NULL, ...)
{
  f <- model.terms(x$formula, model)
  class(f) <- "bamlss.formula"
  return(f)
}

formula.bamlss.terms <- function(x, model, ...)
{
  if(!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    x <- list(x)
    names(x) <- "formula.1"
  }
  f <- list()
  for(i in names(x)) {
    f[[i]] <- list()
    if(!inherits(x[[i]], "terms")) {
      for(j in names(x[[i]])) {
        f[[i]][[j]] <- list()
        f[[i]][[j]]$formula <- x[[i]][[j]]
        env <- environment(x[[i]][[j]])
        attributes(f[[i]][[j]]$formula) <- NULL
        environment(f[[i]][[j]]$formula) <- env
        vars <- all.vars(x[[i]][[j]])
        response <- response.name(x[[i]][[j]])
        if(all(is.na(response)))
          response <- NULL
        if(!is.null(response)) {
          response <- NULL
          vars <- vars[-1]
        }
        f[[i]][[j]]$fake.formula <- as.formula(paste(response, "~1", if(length(vars)) "+" else NULL,
          paste(vars, collapse = "+")), env = environment(x[[i]][[j]]))
        f[[i]][[j]]$terms <- x[[i]][[j]]
      }
    } else {
      f[[i]]$formula <- x[[i]]
      env <- environment(x[[i]])
      attributes(f[[i]]$formula) <- NULL
      environment(f[[i]]$formula) <- env
      vars <- all.vars(x[[i]])
      response <- response.name(x[[i]])
      if(all(is.na(response)))
        response <- NULL
      if(!is.null(response)) {
        response <- NULL
        vars <- vars[-1]
      }
      f[[i]]$fake.formula <- as.formula(paste(response, "~1", if(length(vars)) "+" else NULL,
        paste(vars, collapse = "+")), env = environment(x[[i]]))
      f[[i]]$terms <- x[[i]]
    }
  }
  class(f) <- c("bamlss.formula", "list")
  environment(f) <- environment(x)
  return(f)
}

print.bamlss.formula <- function(x, ...) {
  if(!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    print(x)
  } else {
    nx <- names(x)
    if(is.null(nx))
      nx <- as.character(1:length(x))
    for(i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if(inherits(x[[i]], "list") & "h1" %in% names(x[[i]])) {
        for(j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          attr(x[[i]][[j]], "name") <- NULL
          attr(x[[i]][[j]]$formula, ".Environment") <- NULL
          print(x[[i]][[j]]$formula, showEnv = FALSE)
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        attr(x[[i]]$formula, ".Environment") <- NULL
        if("formula" %in% names(x[[i]])) print(x[[i]]$formula, showEnv = FALSE) else print(x[[i]])
      }
      if(i < length(x))
      cat("\n")
    }
  }
  invisible(NULL)
}


## Drop terms from "bamlss.terms'.
drop.terms.bamlss <- function(f, pterms = TRUE, sterms = TRUE,
  specials = NULL, keep.response = TRUE, data = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  if(!inherits(f, "formula")) {
    if(!is.null(f$terms)) {
      f <- f$terms
    } else {
      if(!is.null(f$formula))
        f <- f$formula
    }
  }
  tx <- if(!inherits(f, "terms")) {
    terms.formula(f, specials = specials, keep.order = TRUE, data = data)
  } else f
  specials <- unique(c(names(attr(tx, "specials")), specials))
  sid <- unlist(attr(tx, "specials"))
  tl <- attr(tx, "term.labels")
  sub <- attr(tx, "response")
  if(!is.null(sid)) {
    st <- tl[sid - sub]
    pt <- tl[-1 * (sid - sub)]
  } else {
    st <- character(0)
    pt <- tl
  }
  if(!sterms & length(st)) {
    st <- paste("-", st, collapse = "")
    st <- as.formula(paste(". ~ .", st))
    tx <- terms.formula(update(tx, st), specials = specials, keep.order = TRUE, data = data)
  }
  if(!pterms & length(pt)) {
    tl <- attr(tx, "term.labels")
    sid <- unlist(attr(tx, "specials"))
    if(!is.null(sid)) {
      st <- tl[sid - 1L]
      pt <- tl[-1 * (sid - 1L)]
    } else {
      st <- character(0)
      pt <- tl
    }
    pt <- paste("-", pt, collapse = "")
    pt <- as.formula(paste(". ~ .", pt))
    tx <- terms.formula(update(tx, pt), specials = specials, keep.order = TRUE, data = data)
  }
  class(tx) <- c("formula", "terms")
  environment(tx) <- environment(f)
  if(!keep.response)
    tx <- delete.response(tx)
  tx
}

has_dot <- function(formula) {
  inherits(try(terms(formula), silent = TRUE), "try-error")
}

terms.bamlss <- terms.bamlss.frame <- terms.bamlss.formula <- function(x, specials = NULL,
  data = NULL, model = NULL, pterms = TRUE, sterms = TRUE, drop = TRUE, ...)
{
  if(inherits(x, "bamlss.frame"))
    x <- formula(x)
  if(!inherits(x, "bamlss.formula"))
    x <- bamlss.formula(x, ...)
  env <- environment(x)
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  elmts <- c("formula", "fake.formula")
  if(!any(names(x) %in% elmts) & !inherits(x, "formula")) {
    if(!is.null(model)) {
      if(is.character(model)) {
        if(all(is.na(pmatch(model[1], names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model[1]) > length(x) || is.na(model[1]) || min(model[1]) < 1) 
          stop("argument model is specified wrong!")
      }
      if(length(model) > 1)
        model <- model[1:2]
      if(length(model) < 2) {
        x <- x[model]
      } else {
        x <- x[[model[1]]]
        if(is.character(model)) {
          if(all(is.na(pmatch(model[2], names(x)))))
            stop("argument model is specified wrong!")
        } else {
          if(max(model[2]) > length(x) || is.na(model[2]) || min(model[2]) < 1) 
            stop("argument model is specified wrong!")
        }
        x <- x[model[2]]
      }
    }
  } else x <- list(x)

  rval <- list()
  if(is.null(nx <- names(x))) {
    nx <- paste("formula", 1:length(x), sep = ".")
    names(x) <- nx
  }
  for(i in seq_along(nx)) {
    if(!any(names(x[[nx[i]]]) %in% elmts) & !inherits(x[[nx[i]]], "formula")) {
      rval[[nx[i]]] <- list()
      nx2 <- names(x[[nx[i]]])
      for(j in seq_along(nx2)) {
        rval[[nx[i]]][[nx2[j]]] <- drop.terms.bamlss(x[[nx[i]]][[nx2[j]]],
          pterms = pterms, sterms = sterms, specials = specials, data = data)
      }
    } else {
      rval[[nx[i]]] <- drop.terms.bamlss(x[[nx[i]]], pterms = pterms,
        sterms = sterms, specials = specials, data = data)
    }
  }

  if(drop & (length(rval) < 2)) {
    rval <- rval[[1]]
  } else {
    class(rval) <- c("bamlss.terms", "list")
  }
  environment(rval) <- env

  rval
}


## Model terms extractor function for formulas and 'bamlss.frame'.
model.terms <- function(x, model = NULL, part = c("x", "formula", "terms"))
{
  if(!inherits(x, "bamlss.formula")) {
    if(inherits(x, "bamlss.frame")) {
      part <- match.arg(part)
      if(is.null(x[[part]]))
        stop(paste("cannot find object", part, "in 'bamlss.frame' object!"))
      x <- x[[part]]
    } else stop(paste("cannot extract parts from object of class '", class(x), "'!", sep = ""))
  }
  if(is.null(model))
    return(x)
  cx <- class(x)
  env <- environment(x)
  elmts <- c("formula", "fake.formula")
  if(!any(names(x) %in% elmts)) {
    if(is.character(model)) {
      if(all(is.na(pmatch(model[1], names(x)))))
        stop("argument model is specified wrong!")
    } else {
      if(max(model[1]) > length(x) || is.na(model[1]) || min(model[1]) < 1) 
        stop("argument model is specified wrong!")
    }
    if(length(model) > 1)
      model <- model[1:2]
    if(length(model) < 2) {
      x <- x[model]
    } else {
      x <- x[[model[1]]]
      if(is.character(model)) {
        if(all(is.na(pmatch(model[2], names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model[2]) > length(x) || is.na(model[2]) || min(model[2]) < 1) 
          stop("argument model is specified wrong!")
      }
      x <- x[model[2]]
    }
  } else x <- list(x)
  class(x) <- cx
  environment(x) <- env
  return(x)
}


## Some simple check functions for 'term' objects.
has_intercept <- function(x)
{
  if(inherits(x, "formula"))
    x <- terms(x)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(attr(x, "intercept") > 0)
}

has_response <- function(x)
{
  if(inherits(x, "formula"))
    x <- terms(x)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(attr(x, "response") > 0)
}

has_sterms <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  if(inherits(x, "formula"))
    x <- terms(x, specials = specials)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(length(unlist(attr(x, "specials"))) > 0)
}

has_pterms <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  if(inherits(x, "formula"))
    x <- terms(x, specials = specials)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  x <- drop.terms.bamlss(x, pterms = TRUE, sterms = FALSE, specials = specials, keep.response = FALSE)
  fc <- length(attr(x, "factors")) > 0
  ic <- attr(x, "intercept") > 0
  return(fc | ic)
}

get_pterms_labels <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  tl <- if(has_pterms(x, specials)) {
    x <- drop.terms.bamlss(x, pterms = TRUE, sterms = FALSE,
      keep.response = FALSE, specials = specials)
    c(attr(x, "term.labels"), if(attr(x, "intercept") > 0) "(Intercept)" else NULL)
  } else character(0)
  tl
}

get_sterms_labels <- function(x, specials = NULL)
{
  env <- environment(x)
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti"))
  if(has_sterms(x, specials)) {
    x <- drop.terms.bamlss(x, pterms = FALSE, sterms = TRUE,
      keep.response = FALSE, specials = specials)
    tl <- NULL
    for(j in attr(x, "term.labels")) {
      st <- try(eval(parse(text = j), envir = env), silent = TRUE)
      if(inherits(st, "try-error"))
        st <- eval(parse(text = j), enclos = env, envir = loadNamespace("mgcv"))
      tl <- c(tl, st$label)
    }
  } else tl <- character(0)
  tl
}


list2mat <- function(x)
{
  do.call("cbind", x)
}


## Process results with samples and bamlss.frame.
results.bamlss.default <- function(x, what = c("samples", "parameters"), grid = 100, nsamps = NULL, ...)
{
  if(!inherits(x, "bamlss.frame") & !inherits(x, "bamlss"))
    stop("x must be a 'bamlss' object!")
  if(is.null(x$samples) & is.null(x$parameters)) {
    warning("nothing to do!")
    return(NULL)
  }

  if(is.null(x$x))
    stop("cannot compute results, 'x' object is missing, see design.construct()!")

  what <- match.arg(what)
  if(!is.null(x$samples) & what == "samples") {
    samps <- samples(x)
    if(!is.null(nsamps)) {
      i <- seq(1, nrow(samps), length = nsamps)
      samps <- samps[i, , drop = FALSE]
    }
  } else {
    if(is.null(x$parameters)) {
      warning("nothing to do!")
      return(NULL)
    }
    samps <- parameters(object, model = model, list = FALSE)
    cn <- names(samps)
    samps <- matrix(samps, nrow = 1)
    colnames(samps) <- cn
    samps <- as.mcmc(samps)
  }

  family <- x$family
  snames <- colnames(samps)
  mf <- model.frame(x)

  make_results <- function(obj, id = NULL)
  {
    DIC <- pd <- NA
    if(any(grepl("deviance", snames))) {
      DIC <- as.numeric(samps[, grepl("deviance", snames)])
      pd <- var(DIC, na.rm = TRUE) / 2
      DIC <- mean(DIC, na.rm = TRUE)
    }
    if(any(grepl("logLik", snames))) {
      DIC <- -2 * as.numeric(samps[, grepl("logLik", snames)])
      pd <- var(DIC, na.rm = TRUE) / 2
      DIC <- mean(DIC, na.rm = TRUE)
    }
    IC <- c("DIC" = DIC, "pd" = pd)

    ## Compute model term effects.
    p.effects <- s.effects <- s.effects.resmat <- NULL

    ## Parametric effects.
    if(has_pterms(obj$terms)) {
      tl <- get_pterms_labels(obj$terms)
      sn <- paste(id, "p", tl, sep = ".")
      i <- grep2(sn, snames, fixed = TRUE)
      if(length(i)) {
        psamples <- as.matrix(samps[, snames[i], drop = FALSE])
        nas <- apply(psamples, 1, function(x) { any(is.na(x)) } )
        psamples <- psamples[!nas, , drop = FALSE]
        qu <- t(apply(psamples, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
        sd <- drop(apply(psamples, 2, sd, na.rm = TRUE))
        me <- drop(apply(psamples, 2, mean, na.rm = TRUE))
        p.effects <- cbind(me, sd, qu)
        rownames(p.effects) <- gsub(paste(id, "p.", sep = "."), "", snames[i], fixed = TRUE)
        colnames(p.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
      } else stop(paste("cannot find samples for terms: ", paste(tl, sep = ", "), "!", sep = ""))
    }

    ## Smooth effects.
    if(has_sterms(obj$terms)) {
      tl <- get_sterms_labels(obj$terms)
      sn <- paste(id, "s", tl, sep = ".")
      i <- grep2(sn, snames, fixed = TRUE)
      if(length(i)) {
        for(j in tl) {
          sn <- paste(id, "s", j, sep = ".")
          psamples <- as.matrix(samps[, snames[grep2(sn, snames, fixed = TRUE)], drop = FALSE])
          nas <- apply(psamples, 1, function(x) { any(is.na(x)) } )
          psamples <- psamples[!nas, , drop = FALSE]
       
          ## FIXME: retransform!
          if(!is.null(obj$smooth.construct[[j]]$Xf) & FALSE) {
            stop("no randomized terms supported yet!")
            kx <- ncol(obj$smooth.construct[[j]]$Xf)
            if(kx) {
              pn <- paste(paste(id, ":h1:linear.",
                paste(paste(obj$smooth.construct[[j]]$term, collapse = "."), "Xf", sep = "."), sep = ""),
                1:kx, sep = ".")
              xsamps <- matrix(samples[[j]][, snames %in% pn], ncol = kx)
              psamples <- cbind("ra" = psamples, "fx" = xsamps)
              re_trans <- function(g) {
                g <- obj$smooth.construct[[j]]$trans.D * g
                if(!is.null(obj$smooth.construct[[j]]$trans.U))
                  g <- obj$smooth.construct[[j]]$trans.U %*% g
                g
              }
              psamples <- t(apply(psamples, 1, re_trans))
            }
          }

          ## Possible variance/edf parameter samples.
          vsamples <- NULL
          tau2 <- paste(id, "s", j, "tau2", sep = ".")
          if(length(tau2 <- grep(tau2, snames, fixed = TRUE))) {
            vsamples <- samps[, tau2, drop = FALSE]
            vsamples <- vsamples[!nas, , drop = FALSE]
          } else {
            if(obj$smooth.construct[[j]]$fixed)
              vsamples <- rep(.Machine$double.eps, nrow(samps))
          }
          edfsamples <- NULL
          edf <- paste(id, "s", j, "edf", sep = ".")
          if(length(edf <- grep(edf, snames, fixed = TRUE))) {
            edfsamples <- samps[, edf, drop = FALSE]
            edfsamples <- edfsamples[!nas, , drop = FALSE]
          }

          ## Acceptance probalities.
          asamples <- NULL
          alpha <- paste(id, "s", j, "alpha", sep = ".")
          if(length(alpha <- grep(alpha, snames, fixed = TRUE))) {
            asamples <- as.numeric(samps[, alpha])
            asamples <- asamples[!nas]
          }

          ## Prediction matrix.
          get.X <- function(x) { ## FIXME: time(x)
            X <- PredictMat(obj$smooth.construct[[j]], x)
            X
          }

          ## Compute effect.
          if(!is.list(s.effects))
            s.effects <- list()
          if(length(s.effects)) {
            if(obj$smooth.construct[[j]]$label %in% names(s.effects)) {
              ct <- gsub(".smooth.spec", "", class(obj$smooth.construct[[j]]))[1]
              if(ct == "random.effect") ct <- "re"
              obj$smooth.construct[[j]]$label <- paste(obj$smooth.construct[[j]]$label, ct, sep = ".")
            }
          }
          if(is.null(obj$smooth.construct[[j]]$fit.fun)) {
            obj$smooth.construct[[j]]$fit.fun <- function(X, b, ...) {
              drop(X %*% b)
            }
          }

          b <- grep2(paste(id, "s", j, "b", sep = "."), colnames(psamples), fixed = TRUE)

          fst <- compute_term(obj$smooth.construct[[j]], get.X = get.X,
            fit.fun = obj$smooth.construct[[j]]$fit.fun,
            psamples = psamples[, b, drop = FALSE], vsamples = vsamples, asamples = asamples,
            FUN = NULL, snames = snames, s.effects.resmat = s.effects.resmat,
            data = mf[, obj$smooth.construct[[j]]$term, drop = FALSE],
            grid = grid, edfsamples = edfsamples)

          ## Add term to effects list.
          s.effects[[obj$smooth.construct[[j]]$label]] <- fst$term
          s.effects.resmat <- fst$s.effects.resmat
          remove(fst)
        }
      } else stop(paste("cannot find samples for terms: ", paste(tl, sep = ", "), "!", sep = ""))
    }

    rval <- list(
      "model" = list("formula" = obj$formula,
        "DIC" = DIC, "pd" = pd, "N" = nrow(mf)),
      "p.effects" = p.effects, "s.effects" = s.effects,
      "s.effects.resmat" = s.effects.resmat
    )

    class(rval) <- "bamlss.results"
    return(rval)
  }

  rval <- list()
  nx <- names(x$terms)
  for(j in nx) {
    rval[[j]] <- make_results(x$x[[j]], id = j)
    if(!is.null(rval[[j]]$s.effects)) {
      for(i in seq_along(rval[[j]]$s.effects)) {
        specs <- attr(rval[[j]]$s.effects[[i]], "specs")
        specs$label <- paste(specs$label, j, sep = ".")
        attr(rval[[j]]$s.effects[[i]], "specs") <- specs
      }
    }
  }

  class(rval) <- "bamlss.results"
  return(rval)
}


## Fitted values/terms extraction
fitted.bamlss <- function(object, model = NULL, term = NULL,
  type = c("link", "parameter"), samples = FALSE, FUN = mean,
  nsamps = NULL, ...)
{
  type <- match.arg(type)
  family <- object$family

  if(type != "parameter" & !samples)
    object <- model.terms(object, model)

  h1check <- any(grepl("h1", names(object$terms)))

  elmts <- c("formula", "fake.formula", "model", "p.effects",
    "effects", "fitted.values", "residuals")

  one <- FALSE
  if(any(elmts %in% names(object))) {
    if(!samples)
      object <- list(object)
    one <- TRUE
  }

  rval <- vector(mode = "list", length = length(object))
  names(rval) <- names(object)
  nrval <- if(is.null(names(rval))) 1:length(object) else names(rval)

  if(!samples) {
    for(j in nrval) {
      if(!any(elmts %in% names(object[[j]]))) {
        rval[[j]] <- fitted.bamlss(object[[j]], term = term, ...)
      } else {
        if(is.null(term))
          rval[[j]] <- object[[j]]$fitted.values
        else {
          if(!is.null(object[[j]]$effects)) {
            fe <- list()
            ne <- names(object[[j]]$effects)
            for(i in seq_along(term)) {
              if(length(e <- grep(term[i], ne, fixed = TRUE)))
                fe[[ne[e[1]]]] <- object[[j]]$effects[[e[1]]]
            }
            if(length(fe))
              rval[[j]] <- fe
          }
        }
      }
    }
    if(type != "link")
      stop("to compute fitted values on the parameter scale use option samples = TRUE!")
  } else {
    if(is.null(model))
      model <- nrval
    if(!is.character(model))
      model <- nrval[model]
    ind <- if(one) 1 else seq_along(model)
    for(j in ind) {
      if(!one & !all(c("model", "fitted.values") %in% names(object[[j]]))) {
        rval[[j]] <- 0
        for(jj in seq_along(object[[j]])) {
          if(is.null(object[[j]][[jj]]$p.effects) & is.null(object[[j]][[jj]]$effects)) next
          rval[[j]] <- rval[[j]] + predict.bamlss(object, model = c(j, jj), term = term,
            FUN = function(x) { x }, nsamps = nsamps, ...)
        }
      } else {
        if(is.null(object[[j]]$p.effects) & is.null(object[[j]]$effects)) next
        rval[[j]] <- predict.bamlss(object, model = if(one) NULL else j, term = term,
          FUN = function(x) { x }, nsamps = nsamps, ...)
      }
      if(type != "link" & !h1check)
        rval[[j]] <- apply(rval[[j]], 2, make.link2(family$links[if(one) 1 else nrval[j]])$linkinv)
      if(!h1check) {
        if(!is.null(dim(rval[[j]])))
          rval[[j]] <- t(apply(rval[[j]], 1, FUN))
        if(!is.null(dim(rval[[j]]))) {
          if(nrow(rval[[j]]) == 1) {
            rval[[j]] <- drop(rval[[j]])
          } else {
            if(ncol(rval[[j]]) == 1)
              rval[[j]] <- drop(rval[[j]])
          }
        }
      }
    }
    if(h1check) {
      rval <- delete.NULLs(rval)
      rval <- do.call("+", rval)
      if(type != "link")
        rval <- apply(rval, 2, make.link2(family$links[1])$linkinv)
      rval <- t(apply(rval, 1, FUN))
      if(!is.null(dim(rval))) {
        if(nrow(rval) == 1) {
          rval <- drop(rval)
        } else {
          if(ncol(rval[[j]]) == 1)
            rval <- drop(rval)
        }
      }
    }
  }

  rval <- delete.NULLs(rval)
  if(length(rval) < 2)
    rval <- rval[[1]]

  rval
}

## Functions for model samples
grep2 <- function(pattern, x, ...) {
  i <- NULL
  for(p in pattern)
    i <- c(i, grep(p, x, ...))
  sort(unique(i))
}

samples <- function(x, model = NULL, term = NULL, combine = TRUE, drop = TRUE, ...)
{
  if(!inherits(x, "bamlss") & !inherits(x, "bamlss.frame"))
    stop("x is not a 'bamlss' object!")
  if(is.null(x$samples))
    stop("no samples to extract!")
  tx <- terms(x)
  x <- x$samples

  x <- process.chains(x, combine, drop = FALSE)

  if(is.null(model) & is.null(term)) {
    if(length(x) < 2)
      x <- x[[1]]
    return(x)
  }

  snames <- colnames(x[[1]])
  nx <- names(tx)

  if(!is.null(model)) {
    model <- model[1]
    i <- if(is.character(model)) {
      pmatch(model, nx)
    } else {
      if(length(model) > length(tx)) NA else model
    }
    if(is.na(i))
      stop("cannot find model!")
    j <- grep(nx[i], snames, fixed = TRUE, value = TRUE)
    for(k in seq_along(x)) {
      x[[k]] <- x[[k]][, j]
    }
    tx <- tx[nx[i]]
    snames <- colnames(x[[1]])
  }

  if(!is.null(term)) {
    term <- term[1]
    if(!is.character(term)) {
      if(term < 1)
        term <- "(Intercept)"
    }
    rval <- vector(mode = "list", length = length(x))
    nx <- names(tx)
    for(i in seq_along(tx)) {
      tl <- all.labels.formula(tx[[i]], full.names = TRUE)
      if(attr(tx[[i]], "intercept") > 0)
        tl <- c(tl, "(Intercept)")
      if(is.character(term)) {
        j <- grep(term, tl, fixed = TRUE)
        if(!length(j))
          j <- NA
      } else {
        j <- if(length(term) > length(tl)) NA else term
      }
      if(is.na(j))
        next
      specials <- unlist(attr(tx[[i]], "specials"))
      jj <- grep(tl[j], snames, fixed = TRUE, value = TRUE)
      for(ii in jj) {
        for(k in seq_along(x))
          rval[[k]] <- cbind(rval[[k]], x[[k]][, ii, drop = FALSE])
      }
    }
    for(k in seq_along(x))
      rval[[k]] <- as.mcmc(rval[[k]], start = start(x[[k]]), end = end(x[[k]]))
    x <- as.mcmc.list(rval)
  }

  if(drop & (length(x) < 2))
    x <- x[[1]]
  
  return(x)
}


## Credible intervals of coefficients.
confint.bamlss <- function(object, parm, level = 0.95, model = NULL, ...)
{
  args <- list(...)
  if(!is.null(args$term))
    parm <- args$term
  if(missing(parm))
    parm <- term.labels(object, ne = TRUE, id = FALSE)
  samps <- samples(object, model = model, term = parm)
  np <- colnames(samps)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  apply(samps, 2, quantile, probs = probs)
}


## Extract model coefficients.
coef.bamlss <- function(object, model = NULL, term = NULL, FUN = mean, ...)
{
  object <- model.terms(object)
  if(is.null(term))
    term <- term.labels(object, ne = TRUE, id = FALSE)
  samps <- samples(object, model = model, term = term)
  apply(samps, 2, function(x) { FUN(na.omit(x), ...) })
}


## Get all terms names used.
term.labels <- function(x, model = NULL, pterms = TRUE, sterms = TRUE,
  intercept = TRUE, list = TRUE, ...)
{
  if(inherits(x, "bamlss") | inherits(x, "bamlss.frame")) {
    x <- terms(x)
  } else {
    if(!inherits(x, "bamlss.terms")) {
      if(inherits(x, "terms")) {
        x <- list("p" = x)
      } else stop("x must be a 'terms' or 'bamlss.terms' object!")
    }
  }

  nx <- names(x)

  if(is.null(model)) {
    model <- nx
  } else {
    if(!is.character(model)) {
      if(max(model) > length(nx) | min(model) < 1)
        stop("model is specified wrong")
    } else {
      for(j in seq_along(model)) {
        mm <- pmatch(model[j], nx)
        if(is.na(mm))
          stop("model is specified wrong")
        model[j] <- nx[mm]
      }
    }
  }

  x <- x[model]
  nx <- names(x)
  rval <- vector(mode = "list", length = length(nx))
  for(j in seq_along(x)) {
    rval[[j]] <- list()
    txj <- drop.terms.bamlss(x[[j]], pterms = pterms, sterms = sterms, keep.response = FALSE)
    tl <- attr(txj, "term.labels")
    specials <- unlist(attr(txj, "specials"))
    if(length(specials)) {
      sub <- if(attr(txj, "response") > 0) 1 else 0
      rval[[j]]$p <- tl[-1 * c(specials - sub)]
      rval[[j]]$s <- tl[specials - sub]
    } else {
      rval[[j]]$p <- tl
    }
    if(intercept & (attr(x[[j]], "intercept") > 0)) {
      rval[[j]]$p <- if(!length(tl)) "(Intercept)" else c("(Intercept)", rval[[j]]$p)
    }
  }
  names(rval) <- nx
  if(!list)
    rval <- unlist(rval)

  rval
}

term.labels2 <- function(x, model = NULL, pterms = TRUE, sterms = TRUE,
  intercept = TRUE, list = TRUE, ...)
{
  if(inherits(x, "bamlss") | inherits(x, "bamlss.frame")) {
    x <- terms(x)
  } else {
    if(!inherits(x, "bamlss.terms")) {
      if(inherits(x, "terms")) {
        x <- list("p" = x)
      } else stop("x must be a 'terms' or 'bamlss.terms' object!")
    }
  }

  nx <- names(x)

  if(is.null(model)) {
    model <- nx
  } else {
    if(!is.character(model)) {
      if(max(model) > length(nx) | min(model) < 1)
        stop("model is specified wrong")
    } else {
      for(j in seq_along(model)) {
        mm <- pmatch(model[j], nx)
        if(is.na(mm))
          stop("model is specified wrong")
        model[j] <- nx[mm]
      }
    }
  }

  x <- x[model]
  nx <- names(x)
  rval <- vector(mode = "list", length = length(nx))
  for(j in seq_along(x)) {
    txj <- drop.terms.bamlss(x[[j]], pterms = TRUE, sterms = TRUE, keep.response = FALSE)
    rval[[j]] <- attr(txj, "term.labels")
    if(intercept & (attr(txj, "intercept") > 0))
      rval[[j]] <- c(rval[[j]], "(Intercept)")
  }
  names(rval) <- nx
  rval
}


## Scores for model comparison.
score <- function(x, limits = NULL, FUN = function(x) { mean(x, na.rm = TRUE) },
  type = c("mean", "samples"), kfitted = TRUE, nsamps = NULL, ...)
{
  stopifnot(inherits(x, "bamlss"))
  family <- attr(x, "family")
  stopifnot(!is.null(family$d))
  type <- match.arg(type)
  y <- model.response2(x)
  n <- if(is.null(dim(y))) length(y) else nrow(y)
  maxy <- max(y, na.rm = TRUE)

  if(is.null(family$nscore)) {
    nscore <- function(eta) {
      integrand <- function(x) {
        int <- family$d(x, family$map2par(eta))^2
        int[int == Inf | int == -Inf] <- 0
        int
      }
      rval <- if(is.null(limits)) {
          try(integrate(integrand, lower = -Inf, upper = Inf), silent = TRUE)
        } else try(integrate(integrand, lower = limits[1], upper = limits[2]), silent = TRUE)
      if(inherits(rval, "try-error")) {
        rval <- try(integrate(integrand, lower = min(y, na.rm = TRUE),
          upper = max(y, na.rm = TRUE)))
      }
      rval <- if(inherits(rval, "try-error")) NA else rval$value
      rval
    }
  } else {
    nscore <- function(eta) {
	    integrand <- function(x) {
        family$d(x, family$map2par(eta))^2
	    }
	    rval <- sum(integrand(seq(0, maxy)))
	    rval
    }

	  nscore2 <- function(y, eta) {
	    integrand <- function(x) {
         -sum(((x == y) * 1 - family$d(x, family$map2par(eta)))^2)
	    }
	    rval <- (integrand(seq(0, maxy)))
	    rval
    }
  }

  scorefun <- function(eta) {
    norm <- rep(0, n)
    for(i in 1:n) {
      ni <- try(nscore(eta[i, , drop = FALSE]), silent = TRUE)
      if(inherits(ni, "try-error")) ni <- NA
      norm[i] <- ni
    }
    pp <- family$d(y, family$map2par(eta))
    pp[pp == Inf | pp == -Inf] <- 0
    loglik <- log(pp)
    if(is.null(family$nscore)) {
      quadratic <- 2 * pp - norm
    } else {
      quadratic <- rep(0, n)
      for(i in 1:n) {
        ni <- try(nscore2(y[i], eta[i, , drop = FALSE]), silent = TRUE)
        if(inherits(ni, "try-error")) ni <- NA
        quadratic[i] <- ni
      }
    }
    spherical <- pp / sqrt(norm)

    return(data.frame(
      "log" = FUN(loglik),
      "quadratic" = FUN(quadratic),
      "spherical" = FUN(spherical)
    ))
  }

  if(type == "mean") {
    eta <- if(kfitted) {
      kfitted(x, nsamps = nsamps,
        FUN = function(x) { mean(x, na.rm = TRUE) }, ...)
    } else fitted(x, samples = if(!is.null(h_response(x))) TRUE else FALSE)
    if(!inherits(eta, "list")) {
      eta <- list(eta)
      names(eta) <- family$names[1]
    }
    eta <- as.data.frame(eta)
    res <- unlist(scorefun(eta))
  } else {
    nx <- names(x)
    eta <- if(kfitted) {
      kfitted(x, FUN = function(x) { x }, nsamps = nsamps, ...)
    } else fitted(x, samples = TRUE, FUN = function(x) { x }, nsamps = nsamps)
    if(!inherits(eta, "list")) {
      eta <- list(eta)
      names(eta) <- family$names[1]
    }
    for(j in nx) {
      colnames(eta[[j]]) <- paste("i",
        formatC(1:ncol(eta[[j]]), width = nchar(ncol(eta[[1]])), flag = "0"),
        sep = ".")
    }
    nc <- ncol(eta[[1]])
    eta <- as.data.frame(eta)
    res <- list()
    for(i in 1:nc) {
      eta2 <- eta[, grep(ni <- paste(".i.",
        formatC(i, width = nchar(nc), flag = "0"), sep = ""),
        names(eta)), drop = FALSE]
      names(eta2) <- gsub(ni, "", names(eta2))
      res[[i]] <- scorefun(eta2)
    }
    res <- do.call("rbind", res)
  }

  res
}


## Compute fitted values with dropping data.
kfitted <- function(x, k = 5, weighted = FALSE, random = FALSE,
  engine = NULL, verbose = TRUE, FUN = mean, nsamps = NULL, ...)
{
  if(!inherits(x, "bamlss")) stop('argument x is not a "bamlss" object!')
  if(is.null(engine))
    engine <- attr(x, "engine")
  if(is.null(engine)) stop("please choose an engine!")
  mf <- model.frame(x)
  i <- rep(1:k, length.out = nrow(mf))
  if(random)
    i <- sample(i)
  k <- sort(unique(i))
  f <- formula(x)
  family <- family(x)
  ny <- length(unique(attr(mf, "response.name")))
  rval <- NULL
  jj <- 1
  for(j in k) {
    if(verbose) cat("subset:", jj, "\n")
    drop <- mf[i == j, ]
    if(!weighted) {
      take <- mf[i != j, ]
      bcv <- bamlss(f, data = take, family = family,
        engine = engine, verbose = verbose, ...)
    } else {
      w <- 1 * (i != j)
      bcv <- bamlss(f, data = mf, family = family,
        engine = engine, verbose = verbose, weights = w, ...)
    }
    if(!is.null(attr(mf, "orig.names")))
      names(drop) <- rmf(names(drop))
    fit <- fitted.bamlss(bcv, newdata = drop, samples = TRUE, FUN = FUN, nsamps = nsamps)
    if(!inherits(fit, "list")) {
      fit <- list(fit)
      names(fit) <- family$names
    }
    if(is.null(rval)) {
      rval <- list()
      for(ii in names(fit)) {
        rval[[ii]] <- matrix(NA, nrow = nrow(mf),
          ncol = if(is.null(dim(fit[[ii]]))) 1 else ncol(fit[[ii]]))
      }
    }
    for(ii in names(fit)) {
      rval[[ii]][i == j, ] <- fit[[ii]]
    }
    jj <- jj + 1
  }

  for(ii in names(fit)) {
    rval[[ii]] <- if(ncol(rval[[ii]]) > 1) {
      as.data.frame(rval[[ii]])
    } else drop(rval[[ii]])
  }

  if(length(rval) < 2)
    rval <- rval[[1]]

  rval
}


## Modified p() and d() functions.
create.dp <- function(family)
{
  if(is.null(names(family$links)))
    names(family$links) <- family$names
  links <- list()
  for(j in names(family$links))
    links[[j]] <- make.link2(family$links[j])$linkinv
  d <- function(y, eta, ...) {
    for(j in names(eta))
      eta[[j]] <- links[[j]](eta[[j]])
    family$d(y, eta, ...)
  }
  p <- function(y, eta, ...) {
    for(j in names(eta))
      eta[[j]] <- links[[j]](eta[[j]])
    family$p(y, eta, ...)
  }
  return(list("d" = d, "p" = p))
}


## Extract model residuals.
residuals.bamlss <- function(object, type = c("quantile", "ordinary", "quantile2", "ordinary2"),
  samples = FALSE, FUN = mean, nsamps = NULL)
{
  type <- match.arg(type)
  y <- model.response2(object)
  if(!is.null(dim(y))) {
    if(!is.null(hrn <- h_response(object))) {
      if(!samples) stop("need to set argument samples = TRUE using hierarchical models!")
      y <- y[, !(colnames(y) %in% hrn)]
    }
  }
  family <- attr(object, "family")
  if(is.null(family$type)) family$type <- 1
  if(type == "ordinary") {
    if(is.factor(y)) y <- as.integer(y) - 1
  } else stopifnot(!is.null(family$p))

  if(type %in% c("quantile2", "ordinary2")) {
    samples <- TRUE
    FUN2 <- function(x) { x }
  } else FUN2 <- FUN

  res <- fitted(object, type = "link", samples = samples, FUN = FUN2, nsamps = nsamps)
  if(!is.list(res)) {
    res <- list(res)
    names(res) <- family$names
  }

  if(samples & type %in% c("quantile2", "ordinary2")) {
    nc <- ncol(res[[1]])
    res2 <- matrix(NA, nrow(res[[1]]), nc)
    if(!is.null(dim(y)))
      res2 <- rep(list(res2), length = ncol(y))
    for(j in seq_along(res)) {
      colnames(res[[j]]) <- paste("i",
        formatC(1:ncol(res[[j]]), width = nchar(nc), flag = "0"),
        sep = ".")
    }
    res <- as.data.frame(res)
    for(i in 1:nc) {
      eta2 <- res[, grep(ni <- paste(".i.",
        formatC(i, width = nchar(nc), flag = "0"), sep = ""),
        names(res)), drop = FALSE]
      names(eta2) <- gsub(ni, "", names(eta2))
      tres <- if(type == "quantile2") {
        if(family$type == 3) {
          le <- family$p(y - 1, family$map2par(eta2))
          ri <- family$p(y, family$map2par(eta2))
          qnorm(runif(length(y), min = le, max = ri))
        } else qnorm(family$p(y, family$map2par(eta2)))
      } else family$mu(family$map2par(eta2))
      if(is.null(dim(y))) {
        res2[, i] <- unlist(tres)
      } else {
        for(j in 1:ncol(y))
          res2[[j]][, i] <- tres[, j]
      }
    }
    if(is.null(dim(y))) {
      res <- apply(res2, 1, FUN)
    } else {
      res <- list()
      for(j in 1:ncol(y))
        res[[j]] <- apply(res2[[j]], 1, FUN)
      names(res) <- names(y)
      res <- as.data.frame(res)
    }
  }

  if(type %in% c("quantile", "ordinary")) {
    eta <- res
    res <- if(type == "quantile") {
      if(family$type == 3) {
        le <- family$p(y - 1, family$map2par(eta))
        ri <- family$p(y, family$map2par(eta))
        qnorm(runif(length(y), min = le, max = ri))
      } else qnorm(family$p(y, family$map2par(eta)))
    } else family$mu(family$map2par(eta))
  }

  if(type %in% c("ordinary", "ordinary2")) {
    if(is.null(dim(y))) {
      res <- y - res
    } else {
      res <- as.data.frame(res)
      colnames(res) <- colnames(y)
      for(j in 1:ncol(y))
        res[[j]] <- y[, j] - res[[j]]
    }
  } else {
    if(!is.null(dim(res))) {
      res <- as.data.frame(res)
      colnames(res) <- colnames(y)
    }
  }

  res
}


## Extract the model response.
model.response2 <- function(data, hierarchical = FALSE, ...)
{
  if(!inherits(data, "data.frame")) {
    f <- if(inherits(data, "bamlss")) formula(data) else NULL
    data <- model.frame(data)
    if(!is.null(f)) {
      if("h1" %in% names(f)) {
        rn <- all.vars(f$h1)[1]
        attr(data, "response.name") <- rn
      } else {
        rn <- NULL
        for(j in seq_along(f)) {
          if(is.list(f[[j]])) {
            if("h1" %in% names(f[[j]]))
              rn <- c(rn, all.vars(f[[j]]$h1)[1])
          }
        }
        rn <- rn[rn %in% names(data)]
        if(length(rn))
          attr(data, "response.name") <- rn
      }
    }
  }
  rn <- attr(data, "response.name")
  y <- if(is.null(rn)) {
    model.response(data, ...)
  } else data[, unique(rn), ...]
  y
}

## find hierarchical responses
h_response <- function(x)
{
  rval <- NULL
  if(!all(c("model", "fitted.values") %in% names(x))) {
    for(j in seq_along(x))
      rval <- c(rval, h_response(x[[j]]))
  } else {
    if(!is.null(x$model$hlevel)) {
      if(x$model$hlevel > 1)
        rval <- response.name(x$model$formula)
    }
  }
  rval
}


## Create the inverse of a matrix.
matrix_inv <- function(x)
{
  if(length(x) < 2)
    return(1 / x)
  rn <- rownames(x)
  cn <- colnames(x)
  p <- try(chol(x), silent = TRUE)
  p <- if(inherits(p, "try-error")) {
    try(solve(x), silent = TRUE)
  } else {
    try(chol2inv(p), silent = TRUE)
  }
  if(inherits(p, "try-error")) {
    diag(x) <- jitter(diag(x), amount = 1e-5)
    p <- try(solve(x), silent = TRUE)
  }
  rownames(p) <- rn
  colnames(p) <- cn
  return(p)
}


## Compute matching index for duplicates in data.
match.index <- function(x)
{
  if(!is.vector(x)) {
    if(!inherits(x, "matrix") & !inherits(x, "data.frame"))
      stop("x must be a matrix or a data.frame!")
    x <- if(inherits(x, "matrix")) {
      apply(x, 1, paste, sep = "\r")
    } else do.call("paste", c(x, sep = "\r"))
  }
  nodups <- which(!duplicated(x))
  ind <- match(x, x[nodups])
  return(list("match.index" = ind, "nodups" = nodups))
}


XinY <-
    function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by,
             notin = FALSE, incomparables = NULL,
             ...)
{
    fix.by <- function(by, df)
    {
        ## fix up 'by' to be a valid set of cols by number: 0 is row.names
        if(is.null(by)) by <- numeric(0L)
        by <- as.vector(by)
        nc <- ncol(df)
        if(is.character(by))
            by <- match(by, c("row.names", names(df))) - 1L
        else if(is.numeric(by)) {
            if(any(by < 0L) || any(by > nc))
                stop("'by' must match numbers of columns")
        } else if(is.logical(by)) {
            if(length(by) != nc) stop("'by' must match number of columns")
            by <- seq_along(by)[by]
        } else stop("'by' must specify column(s) as numbers, names or logical")
        if(any(is.na(by))) stop("'by' must specify valid column(s)")
        unique(by)
    }

    nx <- nrow(x <- as.data.frame(x)); ny <- nrow(y <- as.data.frame(y))
    by.x <- fix.by(by.x, x)
    by.y <- fix.by(by.y, y)
    if((l.b <- length(by.x)) != length(by.y))
        stop("'by.x' and 'by.y' specify different numbers of columns")
    if(l.b == 0L) {
        ## was: stop("no columns to match on")
        ## returns x
        return(x)
    }
    else {
        if(any(by.x == 0L)) {
            x <- cbind(Row.names = I(row.names(x)), x)
            by.x <- by.x + 1L
        }
        if(any(by.y == 0L)) {
            y <- cbind(Row.names = I(row.names(y)), y)
            by.y <- by.y + 1L
        }
        ## create keys from 'by' columns:
        if(l.b == 1L) {                  # (be faster)
            bx <- x[, by.x]; if(is.factor(bx)) bx <- as.character(bx)
            by <- y[, by.y]; if(is.factor(by)) by <- as.character(by)
        } else {
            ## Do these together for consistency in as.character.
            ## Use same set of names.
            bx <- x[, by.x, drop=FALSE]; by <- y[, by.y, drop=FALSE]
            names(bx) <- names(by) <- paste("V", seq_len(ncol(bx)), sep="")
            bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
            bx <- bz[seq_len(nx)]
            by <- bz[nx + seq_len(ny)]
        }
        comm <- match(bx, by, 0L)
        if (notin) {
            res <- x[comm == 0,]
        } else {
            res <- x[comm > 0,]
        }
    }
    ## avoid a copy
    ## row.names(res) <- NULL
    attr(res, "row.names") <- .set_row_names(nrow(res))
    res
}


XnotinY <-
    function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by,
             notin = TRUE, incomparables = NULL,
             ...)
{
    XinY(x,y,by,by.x,by.y,notin,incomparables)
}


#############################
## (13) Utility functions. ##
#############################
TODOs <- NA
class(TODOs) <- "TODOs"

print.TODOs <- function(x, ...)
{
  todos <- .TODOs(...)
  print(todos, row.names = FALSE)
  invisible(todos)
}

.TODOs <- function(file = NULL)
{
  require("bamlss")
  if(is.null(file))
    file <- "~/svn/bayesr/pkg/bamlss/R/families.R"
  file <- path.expand(file)
  env <- new.env()
  source(file, local = env)
  fun <- grep(".bamlss", ls(env), fixed = TRUE, value = TRUE)
  fun <- fun[!grepl("print", fun)]
  tab <- NULL
  for(i in seq_along(fun)) {
    fe <- try(eval(parse(text = paste(fun[i], "()", sep = "")), envir = env), silent = TRUE)
    if(inherits(fe, "family.bamlss")) {
      if(!is.null(fe$family)) {
        dgp <- try(get(paste("dgp", fe$family, sep = "_")), silent = TRUE)
        dgp <- if(!inherits(dgp, "try-error")) "yes" else "no"
      } else dgp <- "no"
      tab <- rbind(tab, cbind(
        "family" = if(!is.null(fe$family)) fe$family else "na",
        "type" = if(!is.null(fe$type)) fe$type else "na",
        "loglik" = if(!is.null(fe$loglik)) "yes" else "no",
        "scorefun" = if(!is.null(fe$score)) "yes" else "no",
        "weightfun" = if(!is.null(fe$weights)) "yes" else "no",
        "d" = if(!is.null(fe$d)) "yes" else "no",
        "p" = if(!is.null(fe$p)) "yes" else "no",
        "mu" = if(!is.null(fe$mu)) "yes" else "no",
        "dgp" = dgp,
        "BayesX" = if(!is.null(fe$bayesx)) "yes" else "no",
        "JAGS" = if(!is.null(fe$jagstan)) "yes" else "no",
        "IWLS" = if(!is.null(fe$score) & !is.null(fe$weights) & !is.null(fe$d)) "yes" else "no"
      ))
    }
  }
  as.data.frame(tab)
}

#.First.lib <- function(lib, pkg)
#{
#  library.dynam("bamlss", pkg, lib)
#}


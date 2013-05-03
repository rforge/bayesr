################################################
## (1) BayesR main model fitting constructor. ##
################################################
bayesr <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  parse.input = parse.input.bayesr, transform = randomize, setup = setupJAGS,
  sampler = samplerJAGS, results = resultsJAGS,
  cores = NULL, combine = TRUE, ...)
{
  ## Setup all processing functions.
  if(is.null(transform))
    transform <- function(x) { x }
  foo <- list("transform" = transform, "setup" = setup, "sampler" = sampler, "results" = results)
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
    if(!is.function(functions[[j]]))
      stop(paste("argument", nf[j], "is not a function!"))
  }
  names(functions) <- names(foo)
  functions$parse.input <- if(!is.null(parse.input)) {
    stopifnot(is.function(parse.input))
    deparse(substitute(parse.input), backtick = TRUE, width.cutoff = 500)
  } else "parse.input.bayesr"

  ## Parse input.
  pm <- match.call(expand.dots = FALSE)
  pm$parse.input <- pm$setup <- pm$samples <- pm$results <- NULL
  pm[[1]] <- as.name(functions$parse.input)
  pm <- eval(pm, parent.frame())
  pm$call <- match.call()

  ## Transform inputs.
  pm <- functions$transform(pm)

  ## Sampling setup.
  ms <- functions$setup(pm)

  ## Start sampling.
  if(is.null(cores)) {
    so <- functions$sampler(ms)
  } else {
    require("parallel")
    parallel_fun <- function(j) {
      functions$sampler(ms)
    }
    so <- mclapply(1:cores, parallel_fun, mc.cores = cores)
    if(length(so) > 1) {
      so2 <- so[[1]]
      for(j in 2:length(so)) {
        so2 <- c(so2, so[[j]])
      }
      so <- so2
      rm(so2)
    } else so <- so[[1]]
  }

  ## Combine samples.
  if(combine)
    so <- combine_chains(so)

  ## Compute results.
  rval <- functions$results(pm, so)
  attr(rval, "functions") <- functions
  rm(so)

  ## Assign more model information.
  attr(rval, "model.frame") <- pm$mf
  attr(rval, "call") <- pm$call
  attr(rval, "functions") <- functions

  rval
}


##########################################################
## (2) Parsing all input using package mgcv structures. ##
##########################################################
parse.input.bayesr <- function(formula, data, family = gaussian(),
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, knots = NULL, specials = NULL, ...)
{
  formula <- parse.formula.bayesr(formula)
  if(missing(data))
    data <- environment(formula)
  if(is.matrix(data))
    data <- as.data.frame(data)
  specials <- unique(c("s", "te", "t2", "sx", "s2", "rs", specials))
  mt <- terms(formula, specials = specials, keep.order = TRUE)
  tl <- attr(mt, "term.labels")
  sm <- grepl("(", tl, fixed = TRUE)
  pterms <- tl[!sm]
  pf <- paste("~ ", if(attr(mt, "intercept")) 1 else -1,
    if(length(pterms)) paste(" +", paste(pterms, collapse = " + ")), sep = "")
  sterms <- NULL
  if(length(sm <- tl[sm])) {
    for(j in sm) {
      smj <- eval(parse(text = j))
      sterms <- c(sterms, smj$term, smj$by)
      if(!is.null(smj$by.formula)) {
        sterms <- c(sterms, attr(terms(smj$by.formula, keep.order = TRUE), "term.labels"))
      }
    }
  }
  sterms <- unique(sterms[sterms != "NA"])
  response <- as.character(formula)[2]
  fake.formula <- as.formula(paste(response, "~ 1 +", paste(c(pterms, sterms), collapse = " + ")))
  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  mf <- list(formula = fake.formula, data = data, weights = weights, subset = subset,
    offset = offset, na.action = na.action, drop.unused.levels = TRUE)
  mf <- do.call("model.frame", mf)
  for(j in 1:ncol(mf)) {
    if(class(mf[, j]) %in% c("numeric", "integer"))
      if(any(!is.finite(mf[, j]))) {
        warning("infinite values in data, removing these observations in model frame!")
        mf <- mf[is.finite(mf[, j]), ]
      }
  }
  X <- model.matrix(as.formula(pf), data = mf, contrasts.arg = contrasts, ...)
  smooth <- list()
  if(length(sm)) {
    for(j in sm) {
      tsm <- eval(parse(text = j))
      if(is.null(tsm$special)) {
        acons <- TRUE
        if(!is.null(tsm$xt$center))
          acons <- tsm$xt$center
        smt <- smoothCon(tsm, mf, knots, absorb.cons = acons)
      } else {
        smt <- smooth.construct(tsm, mf, knots)
        class(smt) <- c(class(smt), "mgcv.smooth")
        smt <- list(smt)
      }
      smooth <- c(smooth, smt)
    }
    smooth <- mgcv:::gam.side(smooth, X, tol = .Machine$double.eps^.5)
    sme <- mgcv:::expand.t2.smooths(smooth)
    if(is.null(sme)) {
      original.smooth <- NULL
    } else {
      original.smooth <- smooth
      smooth <- sme
      rm(sme)
    }
  }

  rval <- list("formula" = formula, "fake.formula" = fake.formula, "response" = response,
    "X" = X, "smooth" = smooth, "pterms" = pterms, "sterms" = sterms, "mf" = mf, "family" = family)

  rval
}

parse.formula.bayesr <- function(formula)
{
  env <- environment(formula)
  if(is.null(env)) env <- .GlobalEnv
  if(!is.list(formula)) {
    tf <- tempfile()
    capture.output(print(formula), file = tf)
    ft <- paste(readLines(tf), collapse = "")
    on.exit(unlink(tf))
    if(any(grep("|", ft, fixed = TRUE))) {
      formula <- as.list(strsplit(ft, "|", fixed = TRUE)[[1]])
      for(j in seq_along(formula))
        formula[[j]] <- as.formula(formula[[j]], env = env)
    }
  }
  if(is.list(formula)) {
    nf <- names(formula)
    if(is.null(nf)) nf <- rep("", length(formula))
    nf2 <- rep(NA, length(formula))
    for(j in seq_along(formula)) {
      tf <- terms(formula[[j]])
      if(attr(tf, "response") < 1) {
        if(is.null(nf))
          stop("formulae need named responses when supplied as unnamed list!")
        else {
          nf2[j] <- nf[j]
          ft <- paste(nf[j], paste(as.character(formula[[j]]), collapse = " "))
          formula[[j]] <- as.formula(ft)
        }
      } else {
        if(is.null(nf) | nf[j] == "") {
          nf2[j] <- as.character(formula[[j]])[2]
        }
      }
      environment(formula[[j]]) <- env
    }
    names(formula) <- nf2
    check <- as.list(rep(NA, length = length(nf2) - 1))
    for(j in 2:length(nf2)) {
      for(i in seq_along(nf2)) {
        if(j != i) {
          if(p <- length(grep(nf2[j], as.character(formula[[i]]), fixed = TRUE))) {
            if(p > 1) stop(paste("variable", nf2[j], "specified in too often!"))
              formula[[i]] <- c(formula[[i]], formula[[j]])
              formula[[j]] <- NA
          }
        }
      }
    }
    if(any(is.na(unlist(formula)))) {
      formula <- lapply(formula, function(x) { if(length(x) < 2 & is.na(x[1])) NULL else x })
   	  formula <- formula[unlist(lapply(formula, length) != 0)]
    }
  }

  formula
}


###########################
## (3) Utility functions ##
###########################
## Transform smooth terms to mixed model representation.
randomize <- function(x)
{
  if(m <- length(x$smooth)) {
    for(j in 1:m) {
      if(!inherits(x$smooth[[j]], "no.mgcv")) {
        tmp <- mgcv:::smooth2random(x$smooth[[j]], names(x$mf), type = 2)
        if(is.null(x$smooth[[j]]$xt$nolin))
          x$smooth[[j]]$Xf <- tmp$Xf
#        if(inherits(x$smooth[[j]], "random.effect")) {
#          tmp$rand$Xr[tmp$rand$Xr > 0] <- 1
#          tmp$rand$Xr <- scale(tmp$rand$Xr)
#          tmp$trans.D <- rep(1, ncol(tmp$rand$Xr))
#          tmp$trans.U <- diag(1, ncol(tmp$rand$Xr))
#        }
        x$smooth[[j]]$rand <- tmp$rand
        x$smooth[[j]]$trans.D <- tmp$trans.D
        x$smooth[[j]]$trans.U <- tmp$trans.U
      }
    }
  }

  x
}


## Combine sample chains.
combine_chains <- function(x)
{
  rval <- NULL
  for(j in 1:length(x)) {
    rval <- rbind(rval, x[[j]])
  }
  rval <- as.mcmc.list(list(as.mcmc(rval)))
  rval
}


## Combine method for "bayesr" objects.
c.bayesr <- function(...)
{
  objects <- list(...)
  x <- NULL
  for(i in 1L:length(objects))
    x <- c(x, objects[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- c("bayesr", "bayesx")

  return(x)
}


## Model extractor function.
get.model <- function(x, model)
{
  cx <- class(x)
  elmts <- c("formula", "fitted.values", "residuals", "response", "effects", "smooth.hyp", 
    "fixed.effects", "variance", "model.fit", "call")
  if(!any(names(x) %in% elmts)) {
    if(!is.null(model)) {
      if(is.character(model)) {
        if(all(is.na(model <- pmatch(model, names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model) > length(x) || is.na(model) || min(model) < 1) 
          stop("argument model is specified wrong!")
      }
      x <- x[model]
    }
  } else x <- list(x)
  class(x) <- cx

  return(x)
}


## Model frame extractor function.
model.frame.bayesr <- function(formula, parse.input = NULL, ...)
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  if(length(nargs) || is.null(attr(formula, "model.frame"))) {
    fcall <- attr(formula, "call")
    if(is.null(fcall))
      fcall <- formula$call
    if(is.null(fcall)) {
      pin <- if(!is.null(parse.input)) {
        stopifnot(is.function(parse.input))
        deparse(substitute(parse.input), backtick = TRUE, width.cutoff = 500)
      } else "parse.input.bayesr"
      dots$formula <- formula
      rval <- do.call(pin, dots)
      rval <- rval$mf
    } else {
      fcall[[1L]] <- as.name(attr(formula, "functions")$parse.input)
      fcall[names(nargs)] <- nargs
      env <- environment(terms(formula))
      if(is.null(env)) 
        env <- parent.frame()
      rval <- eval(fcall, env)
      rval <- rval$mf
    }
  } else {
    rval <- attr(formula, "model.frame")
  }
  rval
}


## Function to compute statistics from samples of a model term.
compute_term <- function(x, fsamples, psamples, vsamples = NULL,
  FUN = NULL, snames, smooth.hyp, fitted.values, data)
{
  require("coda")
  if(is.null(FUN)) {
    FUN <- function(x) {
      rval <- as.numeric(quantile(x, probs = c(0.025, 0.5, 0.975)))
      names(rval) <- c("2.5%", "50%", "97.5%")
      rval
    }
  }
  smf <- t(apply(fsamples, 1, FUN))
  cnames <- colnames(smf)
  nt <- length(x$term)
  for(l in nt:1) {
    smf <- cbind(data[[x$term[l]]], smf)
  }
  smf <- data.frame(smf)
  names(smf) <- c(x$term, cnames)

  ## Compute new linear predictor.
  fit <- rep(0, nrow(smf))
  if(any(im <- grepl("50%", tolower(colnames(smf)), fixed = TRUE))) {
    im <- c(1:ncol(smf))[im]
    fit <- smf[, im[1]]
  }
  by.drop <- NULL
  if(x$by != "NA") {
    by.drop <- data[[x$by]] == x$by.level
    fit[!by.drop] <- 0
    smf <- smf[by.drop, ]
  }
  fitted.values <- if(!is.null(fitted.values)) fitted.values + fit else fit

  ## Assign class and attributes.
  smf <- as.data.frame(unique(as.matrix(smf)))
  class(smf) <- c(class(x), "data.frame")
  x["X"] <- NULL
  attr(smf, "specs") <- x
  attr(smf, "specs")[c("X", "Xf", "rand", "trans.D", "trans.U")] <- NULL
  class(attr(smf, "specs")) <- class(x)
  attr(smf, "fit") <- fit
  attr(smf, "x") <- data[, x$term]
  attr(smf, "by.drop") <- by.drop
  colnames(psamples) <- paste(x$label, 1:ncol(psamples), sep = ".")

  ## Get samples of the variance parameter.
  if(!is.null(vsamples)) {
    if(!is.matrix(vsamples))
      vsamples <- matrix(vsamples, ncol = 1)
    smatfull <- NULL
    for(j in 1:ncol(vsamples)) {
      qu <- drop(quantile(vsamples[, j], probs = c(0.025, 0.5, 0.975)))
      sd <- sd(vsamples[, j])
      me <- mean(vsamples[, j])
      ra <- range(vsamples[, j])
      smat <- matrix(c(me, sd, qu, ra), nrow = 1)
      colnames(smat) <- c("Mean", "Sd", "2.5%", "50%", "97.5%", "Min", "Max")
      rownames(smat) <- x$label
      smatfull <- rbind(smatfull, smat)
    }
    if(!is.null(smf)) {
      attr(smf, "scale") <- smatfull
      attr(smf, "specs")$label <- gsub(")", paste(",",
        paste(formatC(me, digits = 2), collapse = ","), ")",
        sep = ""), x$label)
      colnames(vsamples) <- paste(x$label, "tau", 1:nrow(smatfull), sep = ".")
      attr(smf, "samples.scale") <- as.mcmc(vsamples)
    }
    smooth.hyp <- rbind(smooth.hyp, smatfull)
  }

  ## Assign samples.
  attr(smf, "samples") <- as.mcmc(psamples)

  return(list("term" = smf, "smooth.hyp" = smooth.hyp, "fitted.values" = fitted.values))
}


#####################
## (4) Prediction. ##
#####################
## A prediction method for "bayesr" objects.
## Prediction can also be based on multiple chains.
predict.bayesr <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = mean, ...)
{
  object <- get.model(object, model)
  k <- length(object)
  enames <- list()
  for(j in 1:k) {
    enames[[j]] <- names(object[[j]]$effects)
  }
  if(length(diff) & !all(diff(sapply(enames, length)) == 0))
    stop("the number of terms in the models is not identical, cannot compute prediction!")
  enames <- data.frame(enames)
  if(!all(apply(enames, 1, function(x) length(unique(x))) == 1))
    stop("different terms in the supplied models, cannot compute prediction!")
  enames <- names(object[[1L]]$effects)
  term <- if(!is.null(enames)) {
    if(is.null(term)) enames else {
      if(is.character(term)) {
        enames[grepl(gsub("[[:space:]]", "", term), enames, fixed = TRUE)]
      } else enames[term]
    }
  } else NULL
  term <- term[!is.na(term)]
  if(!length(term)) term <- NULL
  rval <- NULL

  if(missing(newdata))
    newdata <- model.frame(object[[1L]])
  if(is.character(newdata)) {
    if(file.exists(newdata <- path.expand(newdata)))
      newdata <- read.table(newdata, header = TRUE, ...)
  }
  if(is.matrix(newdata) || is.list(newdata))
    newdata <- as.data.frame(newdata)
  nn <- names(newdata)
  m.samples <- m.designs <- m.specials <- list()
  for(j in 1:k) {
    if(!is.null(term)) {
      for(i in term) {
        specs <- attr(object[[j]]$effects[[i]], "specs")
        if(!all(specs$term %in% nn))
          stop(paste("cannot find variables", specs$term, "in newdata!"))
        if(is.null(specs$is.factor)) specs$is.factor <- FALSE
        tmp <- attr(object[[j]]$effects[[i]], "samples")
        if(is.null(dim(tmp))) {
          tmp <- matrix(tmp, ncol = 1)
        }
        if(!is.null(specs$special)) {
          m.specials[[i]] <- list("X" = PredictMat(specs, newdata),
            "get.mu" = specs$get.mu, "samples" = tmp)
        } else {
          m.samples[[i]] <- rbind(m.samples[[i]], tmp)
          if(inherits(object[[j]]$effects[[i]], "linear.bayesx")) {
            if(specs$is.factor & !is.character(newdata)) {
              hi <- if(!is.null(object[[j]]$fixed.effects)) {
                any(grepl("(Intercept)", rownames(object[[j]]$fixed.effects), fixed = TRUE))
              } else FALSE
              nl <- nlevels(newdata[[i]]) - if(hi) 1 else 0
              if(nl != ncol(m.samples[[i]]))
                stop(paste("levels of factor variable", i, "in newdata not identical to model levels!"))
              f <- as.formula(paste("~", if(hi) "1" else "-1", "+", i))
              if(j < 2) {
                m.designs[[i]] <- model.matrix(f, data = newdata)
                if(hi) m.designs[[i]] <- m.designs[[i]][, -1]
              }
            } else {
              if(j < 2)
                m.designs[[i]] <- newdata[[i]]
            }
          } else {
            if(j < 2) {
              m.designs[[i]] <- if(inherits(specs, "mgcv.smooth")) {
                PredictMat(specs, newdata)
              } else newdata[[specs$label]]
            }
          }
          attr(m.samples[[i]], "is.factor") <- specs$is.factor
        }
      }
    }

    if(intercept & j < 2) {
      sami <- attr(object[[j]]$fixed.effects, "samples")
      if(!is.null(sami)) {
        sami <- sami[, grep("(Intercept)", colnames(sami), fixed = TRUE)]
        if(length(sami)) {
          m.samples$Intercept <- sami
          i <- if(is.null(term)) 1 else {
            length(m.designs) + 1
          }
          m.designs[[i]] <- matrix(1, nrow = nrow(newdata), ncol = 1)
        }
      }
    }
  }
  if(length(m.samples) || length(m.specials)) {
    warn <- getOption("warn")
    options("warn" = -1)
    if(length(m.samples)) {
      m.samples <- as.data.frame(m.samples)
      m.designs <- as.data.frame(m.designs)
      options("warn" = warn)
      get.mu <- function(X, b) {
        as.matrix(X) %*% as.numeric(b)
      }
      rval <- apply(m.samples, 1, function(x) { get.mu(m.designs, x) })
    } else rval <- 0
    if(length(m.specials)) {
      for(i in m.specials) {
        rval <- rval + apply(i$samples, 1, function(x) { i$get.mu(i$X, x) })
      }
    }
    rval <- apply(rval, 1, FUN)
    if(!is.null(dim(rval))) {
      if(nrow(rval) != nrow(newdata))
        rval <- as.data.frame(t(rval))
    }
  } else stop("no model terms selected for prediction!")

  rval
}


####################################
## (5) Creating new smooth terms. ##
####################################
## Setup function for handling "special" model terms.
s2 <- function(...)
{
  rval <- s(...)
  rval$special <- TRUE
  rval$label <- gsub("s(", "s2(", rval$label, fixed = TRUE)
  rval
}


## Setup function for random scaling terms.
rs <- function(..., by = NA)
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
  class(rval) <- "rs.smooth.spec"
  rval
}


## Smooth constructor function for random scaling terms.
smooth.construct.rs.smooth.spec <- function(object, data, knots) {
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
      rval$get.mu <- eval(parse(text = fun))
      rval$rs.by <- rs.by
      rval$by.vars <- vars
      rval$by.formula <- object$by.formula
      rval$one <- attr(ft, "intercept")
    }
  } else {
    rval$get.mu <- function(X, g) {
      X %*% as.numeric(g)
    }
  }

  class(rval) <- "rs.smooth"
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
    object$xt$center <- FALSE
  object$by.done <- TRUE
  if(object$by != "NA") {
    by <- data[[object$by]]
    if(!is.factor(by))
      by <- as.factor(data[[object$by]])
    object$by.levels <- levels(by)
    object$fid <- as.integer(by)
    object$byname <- object$by
    object$by <- "NA"
    object$get.mu <- function(X, g) {
      (g[4] + g[1]) * exp(-(g[5] + g[2]) * exp(-(g[6] + g[3]) * X))
    }
  } else {
    object$get.mu <- function(X, g) {
      g[1] * exp(-g[2] * exp(-g[3] * X))
    }
  }
  class(object) <- c("gc.smooth", "no.mgcv")
  object
}


## Work around for the "prediction matrix" of a growth curve.
Predict.matrix.gc.smooth <- function(object, data, knots) 
{
  X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  X
}


###################
## (6) Plotting. ##
###################
## Plotting method for "bayesr" objects.
plot.bayesr <- function(x, model = NULL, term = NULL, which = 1, ask = FALSE, scale = 1, ...)
{
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  ## What should be plotted?
  which.match <- c("effects", "samples", "hist-resid", "qq-resid",
    "scatter-resid", "scale-resid", "scale-samples", "max-acf")
  if(!is.character(which)) {
    if(any(which > 8L))
      which <- which[which <= 8L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  x <- get.model(x, model)
  n <- length(x)
  args <- list(...)

  ## Effect plotting.
  if(which %in% c("effects", "samples")) {
    k <- 0; ylim <- NULL
    ylim <- args$ylim
    args$residuals <- if(is.null(args$residuals)) FALSE else args$residuals
    if(!is.null(args$ylim) || which == "samples")
      scale <- 0
    ne <- names(x[[1]]$effects)
    if(is.null(term))
      term <- 1:length(ne)
    else {
      if(is.character(term)) {
        tterm <- NULL
        for(j in term)
          tterm <- c(tterm, grep(j, ne, fixed = TRUE))
        term <- na.omit(tterm)
      } else term <- term[term <= length(ne)]
    }
    for(i in 1:n) {
      if(!is.null(x[[i]]$effects)) {
        k <- k + length(term)
        if(scale > 0) {
          for(e in term) {
            et <- x[[i]]$effects[[e]]
            de <- attr(et, "specs")$dim + 1
            ylim <- c(ylim, range(et[, de:ncol(et)]))
            if(args$residuals) {
              res <- attr(et, "partial.resids")
              ylim <- c(ylim, range(res[, de:ncol(res)]))
            }
          }
        }
      }
    }
    if(k < 0) stop("no terms to plot in model object!")
    if(scale > 0)
      ylim <- range(ylim)
    args$ylim <- ylim
    args$which <- which
    if(which == "effects") {
      if(!ask) par(mfrow = n2mfrow(k)) else par(ask = ask)
    } else args$residuals <- NULL
    for(i in 1:n) {
      if(!is.null(x[[i]]$effects)) {
        for(e in term) {
          args$x <- x[[i]]$effects[[e]]
          do.call("plot.bayesr.effect", args)
        }
      }
    }
  }
  if(which == "scale-samples") {
    if(!ask) par(mfrow = n2mfrow(n)) else par(ask = ask)
    for(i in 1:n) {
      args$x <- attr(x[[i]]$variance, "samples")
      do.call("plot", args)
    }
  }

  return(invisible(NULL))
}


## Generic plotting method for model terms.
plot.bayesr.effect <- function(x, which = "effects", ...) {
  if(which == "effects") {
    UseMethod("plot.bayesr.effect")
  } else {
    require("coda")
    args <- list(...)
    args$x <- attr(x, "samples")
    if(!is.null(attr(x, "samples.scale")))
      args$x <- as.mcmc(cbind(as.matrix(args$x), as.matrix(attr(x, "samples.scale"))))
    par(ask = TRUE)
    do.call("plot", args)
  }
}


## Default model term plotting method.
plot.bayesr.effect.default <- function(x, ...) {
  args <- list(...)
  args$x <- x
  if(attr(x, "specs")$dim < 2) {
    if(is.null(args$fill.select))
      args$fill.select <- c(0, 1, 0, 1)
    if(is.null(args$lty))
      args$lty <- c(2, 1, 2)
    if(is.null(args$col.lines))
      args$col.lines <- c(NA, "black", NA)
    do.call("plot2d", args)
  } else {
    do.call("plot3d", args)
  }
}


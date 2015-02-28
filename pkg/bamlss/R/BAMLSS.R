################################################
## (1) BAMLSS main model fitting constructor. ##
################################################
## Could be interesting: http://people.duke.edu/~neelo003/r/
##                       http://www.life.illinois.edu/dietze/Lectures2012/
xreg <- function(formula, family = gaussian.bamlss, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  reference = NULL, parse.input = parse.input.bamlss, transform = transformJAGS,
  setup = setupJAGS, engine = samplerJAGS, results = resultsJAGS,
  cores = NULL, sleep = NULL, combine = TRUE, model = TRUE, grid = 100, ...)
{
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
  } else "parse.input.bamlss"

  ## Parse input.
  pm <- match.call(expand.dots = TRUE)
  pm$parse.input <- pm$setup <- pm$samples <- pm$results <- NULL
  pm[[1]] <- as.name(functions$parse.input)
  pm <- eval(pm, parent.frame())
  formula <- attr(pm, "formula0")

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
    so <- combine_chains(so)

  ## Compute results.
  rval <- functions$results(pm, so)
  rm(so)

  ## Assign more model information.
  attr(rval, "functions") <- functions
  if(model)
    attr(rval, "model.frame") <- attr(pm, "model.frame")
  attr(rval, "family") <- attr(pm, "family")
  attr(rval, "formula") <- formula
  attr(rval, "call") <- match.call()

  rval
}


#########################
## (2) Engine stacker. ##
#########################
stacker <- function(x, optimizer = bfit0, sampler = samplerJAGS,
  cores = NULL, sleep = NULL, ...)
{
  if(is.function(optimizer) | is.character(optimizer))
    optimizer <- list(optimizer)
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
      if(is.null(cores)) {
        x <- j(x, ...)
      } else {
        require("parallel")
        parallel_fun <- function(i) {
          if(i > 1 & !is.null(sleep)) Sys.sleep(sleep)
          j(x, ...)
        }
        x <- mclapply(1:cores, parallel_fun, mc.cores = cores)
      }
    }
  }

  x
}


#########################
## (3) BAMLSS wrapper. ##
#########################
## Using the stacker.
bamlss0 <- function(formula, family = gaussian2, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  optimizer = list(bfit0), sampler = list(GMCMC), cores = NULL, combine = TRUE,
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

  transform <- function(x) { optimizer.setup(x, ...) }
  engine <- function(x) {
    stacker(x, optimizer = optimizer, sampler = sampler, mc.cores = cores,
      n.iter = n.iter, thin = thin, burnin = burnin, seed = seed, ...)
  }
  results <- resultsBayesG

  rval <- xreg(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bamlss, transform = transform,
    setup = FALSE, engine = engine, results = results, cores = NULL,
    combine = combine, sleep = 1, ...)
  
  attr(rval, "engine") <- "stacker"
  attr(rval, "call") <- match.call()
  
  rval
}


## First version, no stacking.
bamlss <- function(formula, family = gaussian2, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  engine = c("BayesG", "BayesX", "JAGS", "STAN"), cores = NULL, combine = TRUE,
  n.iter = 12000, thin = 10, burnin = 2000, seed = NULL, ...)
{
  xengine <- match.arg(engine)
  ff <- try(inherits(family, "family.bamlss"), silent = TRUE)
  if(inherits(ff, "try-error")) {
    family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)
  } else {
    if(is.function(family)) {
      if(inherits(try(family(), silent = TRUE), "try-error"))
        family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)
    }
  }

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

  rval <- xreg(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bamlss, transform = transform,
    setup = setup, engine = engine, results = results, cores = cores,
    combine = combine, sleep = 1, ...)
  
  attr(rval, "engine") <- xengine
  attr(rval, "call") <- match.call()
  
  rval
}


##########################################################
## (4) Parsing all input using package mgcv structures. ##
##########################################################
parse.input.bamlss <- function(formula, data = NULL, family = gaussian2.bamlss,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, knots = NULL, specials = NULL, reference = NULL,
  grid = 100, binning = FALSE, ...)
{
  ## Search for additional formulas
  formula2 <- NULL; formula3 <- list(); k <- 1
  fn <- names(fo <- formals(fun = parse.input.bamlss))[-1]
  fn <- fn[fn != "..."]
  for(f in fn) {
    fe <- eval(parse(text = f))
    if(is.character(fe)) {
      if(grepl("~", fe, fixed = TRUE))
        fe <- as.formula(fe)
    }
    if(inherits(fe, "formula")) {
      formula2 <- c(formula2, fe)
      formula3[[k]] <- fe
      eval(parse(text = paste(f, if(is.null(fo[[f]])) "NULL" else fo[[f]], sep = " = ")))
      k <- k + 1
    }
  }

  ## Parse family object.
  family <- bamlss.family(family)

  ## Parse formula
  formula <- bamlss.formula(c(formula, formula2), specials, family)
  formula0 <- attr(formula, "formula0")

  ## Create the model frame.
  mf <- bamlss.model.frame(formula, data, family, weights,
    subset, offset, na.action, specials)
  response.name <- attr(mf, "response.name")

  ## For categorical responses, extend formula object.
  ylevels <- NULL
  if((length(response.name) < 2) | (cat <- if(is.null(family$cat)) FALSE else family$cat)) {
    if(is.factor(mf[[response.name[1]]])) {
      if(cat & nlevels(mf[[response.name[1]]]) > 1) {
        if(is.null(reference)) {
          ty <- table(mf[[response.name[1]]])
          reference <- c(names(ty)[ty == max(ty)])[1]
        }
        reference <- rmf(reference)
        ylevels <- rmf(levels(mf[[response.name[1]]]))
        ylevels <- ylevels[ylevels != reference]
        if(length(formula) != (n <- length(ylevels)))
          formula <- rep(formula, length.out = n)
        names(formula) <- nf <- paste(response.name[1], ylevels, sep = "")
        for(j in seq_along(formula)) {
          uf <- eval(parse(text = paste(nf[j], " ~ .", sep = "")))
          if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
            formula[[j]][[1]]$cat.formula <- update(formula[[j]][[1]]$formula, uf)
          } else formula[[j]]$cat.formula <- update(formula[[j]]$formula, uf)
        }
      }
    } else mf[[response.name[1]]] <- as.numeric(mf[[response.name[1]]]) ## FIXME: matrices?
  } else {
    for(y in response.name)
      mf[[y]] <- as.numeric(mf[[y]])
  }

  ## Assign all design matrices and the hierarchical level, if any.
  rval <- bamlss.design(formula, mf, contrasts, knots, binning, ...)
  rval <- bamlss.hlevel(rval)

  ## Dirichlet specials.
  if(family$family == "dirichlet") {
    family$ncat <- length(rval)
    family$names <- paste(family$names, 1:family$ncat, sep = "")
    names(rval) <- family$names
  }

  attr(rval, "family") <- family
  attr(rval, "reference") <- reference
  attr(rval, "ylevels") <- ylevels
  attr(rval, "grid") <- grid
  attr(rval, "model.frame") <- mf
  attr(rval, "formula0") <- formula0

  class(rval) <- c("bamlss.input", "list")

  rval
}

"[.bamlss.input" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  xattr <- attributes(x)
  mostattributes(rval) <- attributes(x)
  rval
}

"[.bamlss" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  mostattributes(rval) <- attributes(x)
  rval
}


## Assign all designs matrices.
bamlss.design <- function(x, data, contrasts = NULL, knots = NULL, binning = FALSE,
  before = TRUE, ...)
{
  if(!binning)
    binning <- NULL
  assign.design <- function(obj, mf) {
    if(!all(c("formula", "fake.formula", "response") %in% names(obj)))
      return(obj)
    pf <- paste("~ ", if(obj$intercept) 1 else -1,
      if(length(obj$pterms)) paste(" +", paste(obj$pterms, collapse = " + ")), sep = "")
    obj$X <- model.matrix(as.formula(pf),
      data = mf, contrasts.arg = contrasts, ...)
    obj$pterms <- colnames(obj$X)
    obj$binning <- binning
    rn <- obj$response
    obj$response.vec <- if(!is.null(rn)) mf[[rn]] else NULL
    no.mgcv <- NULL
    if(length(obj$smooth)) {
      smooth <- list()
      for(j in obj$smooth) {
        tsm <- eval(parse(text = j))
        if(is.null(tsm$xt))
          tsm$xt <- list()
        if(is.null(tsm$xt$xbin))
          tsm$xt$xbin <- binning
        if(!is.null(tsm$xt$xbin)) {
          if(!is.logical(tsm$xt$xbin)) {
            for(tsmt in tsm$term) {
              if(!is.factor(mf[[tsmt]]))
                mf[[tsmt]] <- round(mf[[tsmt]], digits = tsm$xt$xbin)
            }
          }
        }
      }
      for(j in obj$smooth) {
        if((rt <- grepl('):s(', j, fixed = TRUE)) | grepl(')/s(', j, fixed = TRUE)) {
          j <- strsplit(j, if(rt) ":" else "/", fixed = TRUE)[[1]]
          j <- paste('rs(', j[1], ',', j[2], ',link="', if(rt) "inverse" else "log", '")', sep = "")
          tsm <- eval(parse(text = j))
          tsm$label <- paste(tsm$smooths[[1]]$label, tsm$smooths[[2]]$label, sep = if(rt) ":" else "/")
        } else tsm <- eval(parse(text = j))
        if(is.null(tsm$special)) {
          if(is.null(tsm$xt))
            tsm$xt <- list()
          if(is.null(tsm$xt$xbin))
            tsm$xt$xbin <- binning
          acons <- TRUE
          if(!is.null(tsm$xt$center))
            acons <- tsm$xt$center
          tsm$xt$center <- acons
          tsm$before <- before
          if(!is.null(tsm$xt$xbin)) {
            term.names <- c(tsm$term, if(tsm$by != "NA") tsm$by else NULL)
            ind <- lapply(mf[, term.names, drop = FALSE], function(x) {
              if(!is.character(x) & !is.factor(x)) {
                if(!is.logical(tsm$xt$xbin))
                  rval <- format(x, digits = tsm$xt$xbin, nsmall = tsm$xt$xbin)
                else rval <- sprintf("%.48f", x)
              } else rval <- x
              rval
            })
            ind <- as.vector(apply(do.call("cbind", ind), 1, paste, collapse = ",", sep = ""))
            uind <- unique(ind)
            tsm$xbin.take <- !duplicated(ind)
            tsm$xbin.ind <- rep(NA, nrow(mf))
            xbin.uind <- seq_along(uind)
            for(ii in xbin.uind)
              tsm$xbin.ind[ind == uind[ii]] <- ii
            tsm$xbin.order <- order(tsm$xbin.ind)
            tsm$xbin.k <- length(xbin.uind)
            tsm$xbin.sind <- tsm$xbin.ind[tsm$xbin.order]
            smt <- smoothCon(tsm, if(before) mf[tsm$xbin.take, term.names, drop = FALSE] else mf, knots, absorb.cons = acons)
          } else {
            smt <- smoothCon(tsm, mf, knots, absorb.cons = acons)
          }
        } else {
          smt <- smooth.construct(tsm, mf, knots)
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
        smooth <- try(gam.side(smooth, obj$X, tol = .Machine$double.eps^.5), silent = TRUE)
        if(inherits(smooth, "try-error"))
          stop("gam.side() produces an error when binning, try to set before = FALSE!?")
        sme <- mgcv:::expand.t2.smooths(smooth)
        if(is.null(sme)) {
          original.smooth <- NULL
        } else {
          original.smooth <- smooth
          smooth <- sme
          rm(sme)
        }
      }
      if(!is.null(no.mgcv))
        smooth <- c(smooth, no.mgcv)
      obj$smooth <- smooth
    }
    if(length(obj$sx.smooth)) {
      sx.smooth <- list()
      for(j in obj$sx.smooth) {
        sx.smooth <- c(sx.smooth, list(eval(parse(text = j))))
      }
      obj$sx.smooth <- sx.smooth
    }

    obj
  }

  if(!all(c("formula", "fake.formula", "response") %in% names(x))) {
    for(j in seq_along(x)) {
      if(!all(c("formula", "fake.formula", "response") %in% names(x[[j]]))) {
        for(i in seq_along(x[[j]])) {
          x[[j]][[i]] <- assign.design(x[[j]][[i]],
            if(i > 1) unique(data[all.vars(x[[j]][[i]]$fake.formula)]) else data)
        }
      } else x[[j]] <- assign.design(x[[j]], data)
    }
  } else x <- assign.design(x, data)

  x
}


## Create the model.frame.
bamlss.model.frame <- function(formula, data, family, weights = NULL,
  subset = NULL, offset = NULL, na.action = na.omit, specials = NULL)
{
  family <- bamlss.family(family)
  formula <- bamlss.formula(formula, specials, family)

  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  if(missing(data))
    data <- environment(formula)
  if(!is.data.frame(data))
    data <- as.data.frame(data)

  ## Make fake "Formula" object.
  fF <- make_fFormula(formula)

  if(!is.null(offset)) {
    if(length(offset) != nrow(data))
      offset <- rep(offset, nrow(data))
  }

  ## Set up the model.frame.
  mf <- list(formula = fF, data = data, weights = weights,
    subset = subset, offset = offset, na.action = na.action,
    drop.unused.levels = TRUE)
  mf <- do.call("model.frame", mf)
  rownames(mf) <- NULL

  ## Remove inf values
  mf <- rm_infinite(mf)

  ## assign response names
  tf <- terms(formula(fF, rhs = 0))
  rn <- as.character(attr(tf, "variables"))[2]
  rn <- strsplit(rn, " | ", fixed = TRUE)[[1]]
  attr(mf, "response.name") <- unique(rn)

  ## Check response.
  if(!is.null(family$valid.response)) {
    family$valid.response(mf[, rn[1]])
  }

  mf
}


## Remove Inf values from data.
rm_infinite <- function(x) {
  if(is.null(dim(x))) return(x)
  if(ncol(x) > 0) {
    for(j in 1:ncol(x)) {
      if(class(x[, j]) %in% c("numeric", "integer")) {
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
      if(is.null(family$family)) stop("family is specidied wrong, no family name available!")
      family <- family$family
    }
    txt <- paste(family, type, sep = if(!is.null(type)) "." else "")
    txt <- gsub("bamlss.bamlss", "bamlss", txt, fixed = TRUE)
    family <- eval(parse(text = txt[1]))
    family <- family()
  }
  if(is.null(family$cat))
    family$cat <- FALSE
  if(is.null(family$mu)) {
    family$mu <- function(x) { make.link2(family$links[1])$linkinv(x[[1]]) }
  }
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
  if(is.null(family$loglik)) {
    if(!is.null(family$d))
      family$loglik <- function(y, eta) { sum(family$d(y, eta, log = TRUE), na.rm = TRUE) }
  }
  if(is.null(family)) family <- list()
  if(is.null(family$score) & FALSE) {
    nf <- family$names
    score <- list()
    for(j in family$names) {
      score[[j]] <- function(y, eta, ...) {
        call <- paste('num_deriv(y, eta, family = family, id = "',  j, '", d = 1)', sep = '')
        eval(parse(text = call))
      }
    }
    family$score <- score
  }
  if(is.null(family$weights) & FALSE) {
    nf <- family$names
    weights <- list()
    for(j in family$names) {
      weights[[j]] <- function(y, eta, ...) {
        call <- paste('num_deriv(y, eta, family = family, id = "',  j, '", d = 2)', sep = '')
        -1 * eval(parse(text = call))
      }
    }
    family$weights <- weights
  }

  family
}


## Special formula parser, can deal with multi parameter models
## and hierarchical structures.
bamlss.formula <- function(formula, specials = NULL, family = gaussian.bamlss())
{
  if(inherits(formula, "bamlss.formula"))
    return(formula)

  specials <- unique(c("s", "te", "t2", "sx", "s2", "rs", specials))

  env <- environment(formula)
  if(is.null(env)) env <- .GlobalEnv
  if(!is.list(formula)) formula <- list(formula)
  if(!length(formula)) stop("formula is specified wrong!")

  complete_formula <- function(formula) {
    if(length(formula) < length(family$names))
      formula <- c(formula, rep(list(), length = length(family$names) - length(formula)))
    names(formula) <- family$names
    if(any(i <- is.na(names(formula))))
      names(formula)[i] <- family$names[i]
    for(j in family$names) {
      if(is.null(formula[[j]])) {
        formula[[j]] <- as.formula(paste(j, "1", sep = " ~ "))
        environment(formula[[j]]) <- env
      }
      attr(formula[[j]], "name") <- j
    }
    formula
  }
  formula <- formula_and(formula, env)
  formula <- formula_at(formula, env)
  formula <- formula0 <- complete_formula(formula_hierarchical(formula))
  formula <- formula_extend(formula, specials, family)

  environment(formula) <- environment(formula0) <- env
  class(formula) <- class(formula0) <- c("bamlss.formula", "list")
  for(j in seq_along(formula0)) {
    if(!inherits(formula0[[j]], "formula")) {
      if(is.null(names(formula0[[j]])))
        names(formula0[[j]]) <- paste("h", 1:length(formula0[[j]]), sep = "")
    }
  }
  attr(formula0, "specials") <- specials
  attr(formula, "formula0") <- formula0
  attr(formula, "raw.formula") <- TRUE

  formula
}


## Make "Formula" object from fake.formulas.
make_fFormula <- function(formula)
{
  fF <- NULL
  for(j in seq_along(formula)) {
    if(!all(c("formula", "intercept", "fake.formula") %in% names(formula[[j]]))) {
      for(i in seq_along(formula[[j]]))
        fF <- c(fF, formula[[j]][[i]]$fake.formula)
    } else {
      fF <- c(fF, formula[[j]]$fake.formula)
    }
  }
  fF <- do.call("as.Formula", fF)
  fF
}


## Extend formula by a fake formula with all variables
## to compute a model.frame, create smooth objects.
formula_extend <- function(formula, specials = NULL, family)
{
  if(is.list(formula)) {
    for(j in seq_along(formula))
      formula[[j]] <- formula_extend(formula[[j]], specials, family)
    return(formula)
  } else {
    specials <- unique(c("s", "te", "t2", "sx", "s2", "rs", specials))
    mt <- terms(formula, specials = specials, keep.order = TRUE)

    get.term.labels <- function(formula) {
      tl <- attr(mt, "term.labels")
      tl2 <- NULL
      for(j in as.character(formula)) {
        if(FALSE & grepl("s(", j , fixed = TRUE) & grepl("/", j, fixed = TRUE)) {
          j2 <- strsplit(j, "/", fixed = TRUE)[[1]]
          for(jj in j2)
            tl <- tl[tl != jj]
          tl <- tl[tl != gsub("/", ":", j)]
          tl2 <- c(tl2, j)
        }
      }
      c(tl, tl2)
    }

    tl <- get.term.labels(formula)

    sm <- rep(NA, length = length(tl))
    for(j in seq_along(tl)) {
      iss <- FALSE
      for(i in specials) {
        if(grepl(paste(i, "(", sep = ""), tl[j], fixed = TRUE))
          iss <- TRUE
      }
      sm[j] <- iss
    }
    pterms <- if(length(tl[!sm])) tl[!sm] else NULL
    intercept <- attr(mt, "intercept") > 0
    sterms <- NULL
    if(length(sm <- tl[sm])) {
      for(j in sm) {
        if(FALSE & ((rt <- grepl(":", j, fixed = TRUE)) | grepl("/", j, fixed = TRUE))) {
          for(jj in strsplit(j, if(rt) ":" else "/", fixed = TRUE)[[1]]) {
            smj <- eval(parse(text = jj))
            sterms <- c(sterms, smj$term, smj$by)
            if(!is.null(smj$by.formula)) {
              sterms <- c(sterms, attr(terms(smj$by.formula, keep.order = TRUE), "term.labels"))
            }
          }
        } else {
          smj <- eval(parse(text = j))
          sterms <- c(sterms, smj$term, smj$by)
          if(!is.null(smj$by.formula)) {
            sterms <- c(sterms, attr(terms(smj$by.formula, keep.order = TRUE), "term.labels"))
          }
        }
      }
    }
    sterms <- unique(sterms[sterms != "NA"])
    response <- formula_respname(formula)
    if(is.na(response) | response %in% family$names) response <- NULL
    fake.formula <- as.formula(paste(response, "~ 1", if(length(c(pterms, sterms))) " + " else NULL,
      paste(c(pterms, sterms), collapse = " + ")), env = environment(formula))

    smooth <- sx.smooth <- NULL
    if(length(sm)) {
      for(j in sm) {
        if("sx(" == paste(strsplit(j, "")[[1]][1:3], collapse = ""))
          sx.smooth <- c(sx.smooth, j)
        else
          smooth <- c(smooth, j)
      }
    }

    return(list("formula" = formula, "intercept" = intercept, "fake.formula" = fake.formula,
      "response" = response, "pterms" = pterms, "sterms" = sterms, "smooth" = smooth,
      "sx.smooth" = sx.smooth))
  }
}

## Get response name.
formula_respname <- function(formula)
{
  rn <- if(!inherits(formula, "list")) {
    f <- as.Formula(formula)
    f <- formula(f, lhs = TRUE, rhs = FALSE)
    if(length(f) < 3) NA else deparse(f[[2]])
  } else NA
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
    names(formula) <- formula_respname(formula[[1]][[1]])
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
    names(formula) <- formula_respname(formula[[1]][[1]])
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
  nf <- sapply(formula, formula_respname)
  snf <- seq_along(nf)
  check <- vector(mode = "list", length = length(formula))
  for(j in snf) {
      for(i in snf) {
        if(j != i) {
          fi <- if(!is.list(formula[[i]])) list(formula[[i]]) else formula[[i]]
          for(jj in seq_along(fi)) {
            av <- all.vars(fi[[jj]])
            rn <- formula_respname(fi[[jj]])
            if(!is.na(rn))
              av <- av[av != rn]
            if(attr(terms(fi[[jj]]), "intercept") < 1) {
              av <- c(av, "-1")
            }
            if(any(av %in% nf[j])) {
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

bamlss.hlevel <- function(x)
{
  do <- TRUE
  while(do) {
    if(!all(c("formula", "fake.formula", "response") %in% names(x))) {
      for(j in seq_along(x)) {
        if(!all(c("formula", "fake.formula", "response") %in% names(x[[j]]))) {
          for(i in seq_along(x[[j]])) {
            rn <- formula_respname(x[[j]][[i]]$formula)
            ind <- seq_along(x[[j]])
            ind <- ind[ind != i]
            base <- TRUE
            for(iii in ind) {
              if(rn %in% all.vars((x[[j]][[iii]]$formula)))
                base <- FALSE
            }
            if(base) {
              x[[j]][[i]]$hlevel <- 1
            } else {
              for(iii in ind) {
                if(rn %in% all.vars((x[[j]][[iii]]$formula)) & !is.null(x[[j]][[iii]]$hlevel))
                  x[[j]][[i]]$hlevel <- x[[j]][[iii]]$hlevel + 1
              }
            }
          }
        } else x[[j]]$hlevel <- 1
      }
    } else x$hlevel <- 1

    do <- drop(unlist(sapply(x, function(x) {
      rval <- NULL
      if(!all(c("formula", "fake.formula", "response") %in% names(x))) {
        for(jj in seq_along(x)) {
          if(!all(c("formula", "fake.formula", "response") %in% names(x[[jj]]))) {
            for(ii in seq_along(x[[jj]])) {
              rval <- c(rval, is.null(x[[jj]][[ii]]$hlevel))
            }
          } else rval <- c(rval, is.null(x[[jj]]$hlevel))
        }
      } else rval <- is.null(x$hlevel)
      rval
    })))

    do <- all(do == TRUE)
  }

  x
}


###########################
## (5) Utility functions ##
###########################
## Transform smooth terms to mixed model representation.
randomize <- function(x)
{
  if(inherits(x, "bamlss.input") & !any(c("smooth", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx)) nx <- 1:length(x)
    if(length(unique(nx)) < length(x)) nx <- 1:length(x)
    for(j in nx)
      x[[j]] <- randomize(x[[j]])
  } else {
    if(m <- length(x$smooth)) {
      for(j in 1:m) {
        if(!inherits(x$smooth[[j]], "no.mgcv")) {
          tmp <- mgcv:::smooth2random(x$smooth[[j]], names(attr(x, "model.frame")), type = 2)
          if(is.null(x$smooth[[j]]$xt$nolin))
            x$smooth[[j]]$Xf <- tmp$Xf
#          if(inherits(x$smooth[[j]], "random.effect")) {
#            tmp$rand$Xr[tmp$rand$Xr > 0] <- 1
#            tmp$rand$Xr <- scale(tmp$rand$Xr)
#            tmp$trans.D <- rep(1, ncol(tmp$rand$Xr))
#            tmp$trans.U <- diag(1, ncol(tmp$rand$Xr))
#          }
          x$smooth[[j]]$rand <- tmp$rand
          x$smooth[[j]]$trans.D <- tmp$trans.D
          x$smooth[[j]]$trans.U <- tmp$trans.U
        }
      }
    }
  }

  x
}


## Assign weights.
assign.weights <- function(x, weights = NULL)
{
  if(!is.null(attr(x, "model.frame"))) {
    if(length(i <- grep("(weights)", names(attr(x, "model.frame")))))
      weights <- attr(x, "model.frame")[, i]
  }

  if(!is.null(weights)) {
    if(inherits(x, "bamlss.input") & !any(c("smooth", "response") %in% names(x))) {
      nx <- names(x)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(x)
      if(length(unique(nx)) < length(x)) nx <- 1:length(x)
      for(j in nx)
        x[[j]] <- assign.weights(x[[j]], weights = weights)
    } else {
      if(!is.null(x$X))
        x$X <- x$X * weights
      if(m <- length(x$smooth)) {
        for(j in 1:m) {
          if(!inherits(x$smooth[[j]], "no.mgcv")) {
            x$smooth[[j]]$X <- x$smooth[[j]]$X * weights
          } 
        }
      }
    }
  }

  x
}


## Combine sample chains.
combine_chains <- function(x)
{
  if(!is.list(x))
    return(x)
  model.specs <- attr(x[[1]], "model.specs")
  if(inherits(x[[1]], "mcmc.list")) {
    x <- as.mcmc.list(do.call("c", x))
  } else {
    stopifnot(inherits(x[[1]], "mcmc"))
    x <- as.mcmc.list(x)
  }
  rval <- NULL
  for(j in 1:length(x)) {
    rval <- rbind(rval, x[[j]])
  }
  rval <- as.mcmc.list(list(as.mcmc(rval)))
  attr(rval, "model.specs") <- model.specs
  rval
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


## Function to compute statistics from samples of a model term.
compute_term <- function(x, get.X, get.mu, psamples, vsamples = NULL,
  asamples = NULL, FUN = NULL, snames, effects.hyp, fitted.values, data,
  grid = 100, rug = TRUE, hlevel = 1, sx = FALSE, re.slope = FALSE)
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
    if(nt < 2)
      data <- unique(data0)
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
  fsamples <- apply(psamples, 1, function(g) {
    get.mu(X, g, expand = FALSE)
  })

  if(is.null(FUN)) {
    FUN <- function(x) {
      rval <- as.numeric(quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      names(rval) <- c("2.5%", "50%", "97.5%")
      rval
    }
  }
  smf <- t(apply(fsamples, 1, FUN))

  cnames <- colnames(smf)
  smf <- as.data.frame(smf)
  for(l in 1:nt) {
    smf <- cbind(data[[tterms[l]]], smf)
  }
  names(smf) <- c(x$term, cnames)

  ## Compute new linear predictor.
  fit <- rep(0, nrow(smf))
  if(any(im <- grepl("50%", tolower(colnames(smf)), fixed = TRUE))) {
    im <- c(1:ncol(smf))[im]
    fit <- smf[, im[1]]
    if(xsmall & !all(is.na(fit)))
      fit <- approx(data[[tterms]], fit, xout = data0[[tterms]])$y
  }

  by.drop <- NULL
  if(x$by != "NA" & !is.null(x$by.level)) {
    by.drop <- (if(xsmall) data0[[x$by]] else data[[x$by]]) == x$by.level
    fit[!by.drop] <- 0
    if(!xsmall)
      smf <- smf[by.drop, ]
  }
  if(is.null(x$xt$center)) {
    mean_fit <- mean(fit, na.rm = TRUE)
    fit <- fit - mean_fit
    if(any(grepl("50%", colnames(smf))))
      smf[, c("2.5%", "50%", "97.5%")] <- smf[, c("2.5%", "50%", "97.5%")] - mean_fit
  } else {
    if(x$xt$center)
      fit <- fit - mean(fit, na.rm = TRUE)
  }

  if(!is.null(fitted.values)) {
    if((length(fit) == length(fitted.values)) | (length(fitted.values) < 2 & fitted.values[1] == 0))
      fitted.values <- fitted.values + fit
  } else fitted.values <- fit

  ## Assign class and attributes.
  smf <- unique(smf)
  if(is.factor(data[, tterms])) {
    bbb <- 1 ## FIXME: factors!
  }
  class(smf) <- c(class(x), "data.frame")
  x["X"] <- NULL
  attr(smf, "specs") <- x
  attr(smf, "specs")[c("X", "Xf", "rand", "trans.D", "trans.U")] <- NULL
  class(attr(smf, "specs")) <- class(x)
  attr(smf, "fit") <- fit
  attr(smf, "x") <- if(xsmall & nt < 2) data0[, tterms] else data[, tterms]
  attr(smf, "by.drop") <- by.drop
  attr(smf, "rug") <- rugp
  colnames(psamples) <- paste(x$label, 1:ncol(psamples), sep = ".")

  ## Get samples of the variance parameter.
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
      colnames(vsamples) <- paste(x$label, "tau", 1:nrow(smatfull), sep = ".")
      attr(smf, "samples.scale") <- as.mcmc(vsamples)
      if(!is.null(asamples)) {
        asamples <- matrix(asamples, ncol = 1)
        colnames(asamples) <- paste(x$label, "alpha", sep = ".")
        attr(smf, "samples.alpha") <- as.mcmc(asamples)
      }
    }
    effects.hyp <- rbind(effects.hyp, smatfull)
  }

  ## Assign samples.
  attr(smf, "samples") <- as.mcmc(psamples)

  return(list("term" = smf, "effects.hyp" = effects.hyp, "fitted.values" = fitted.values))
}


## Function to compute partial residuals.
partial.residuals <- function(effects, response, fitted.values, family)
{
  if(!is.null(response)) {
    if(length(mulink <- family$links[grep("mu", names(family$links))]) < 1)
      mulink <- family$links[1]
    linkfun <- make.link2(mulink[1])$linkfun
    for(i in seq_along(effects)) {
      if(is.factor(response)) response <- as.integer(response) - 1
      e <- linkfun(response) - fitted.values + attr(effects[[i]], "fit")
      if(is.null(attr(effects[[i]], "specs")$xt$center)) {
        e <- e - mean(e, na.rm = TRUE)
      } else {
        if(attr(effects[[i]], "specs")$xt$center)
          e <- e - mean(e, na.rm = TRUE)
      }
      e <- if(is.factor(attr(effects[[i]], "x"))) {
        warn <- getOption("warn")
        options(warn = -1)
        tx <- as.integer(as.character(attr(effects[[i]], "x")))
        options("warn" = warn)
        cbind(if(!any(is.na(tx))) tx else as.integer(attr(effects[[i]], "x")), e)
      } else {
        cbind(attr(effects[[i]], "x"), e)
      }
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
  
  effects
}


## Function to add partial residuals based on weights() and score() function.
add.partial <- function(x, samples = FALSE, nsamps = 100) {
  stopifnot(!is.null(names(x)))
  nx <- names(x)
  family <- attr(x, "family")
  if(is.null(family)) stop("cannot compute partial residuals, no iwls score and weights functions supplied with family object!")
  if(!is.null(family$weights) & !is.null(family$score)) {
    y <- model.response2(x)
    eta <- fitted.bamlss(x, samples = samples, nsamps = nsamps)
    mf <- model.frame(x)
    for(j in seq_along(x)) {
      if(!is.null(x[[j]]$effects)) {
        peta <- family$map2par(eta)
        weights <- family$weights[[nx[j]]](y, peta)
        score <- family$score[[nx[j]]](y, peta)
        z <- eta[[nx[j]]] + 1 / weights * score
        ne <- names(x[[j]]$effects)
        for(sj in seq_along(ne)) {
          f <- predict(x, model = nx[j], term = ne[sj], nsamps = nsamps)
          term <- attr(x[[j]]$effects[[ne[sj]]], "specs")$term
          e <- z - eta[[nx[j]]] + f
          if(is.null(attr(x[[j]]$effects[[ne[sj]]], "specs")$xt$center)) {
            e <- e - mean(e)
          } else {
            if(attr(x[[j]]$effects[[ne[sj]]], "specs")$xt$center)
              e <- e - mean(e)
          }
          e <- data.frame(mf[, term], e)
          names(e) <- c(term, "partial.resids")
          attr(x[[j]]$effects[[ne[sj]]], "partial.resids") <- e
        }
      }
    }
  }
  x
}


#####################
## (6) Prediction. ##
#####################
## A prediction method for "bamlss" objects.
## Prediction can also be based on multiple chains.
predict.bamlss <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = function(x) { mean(x, na.rm = TRUE) }, trans = NULL, MARGIN = 1,
  type = c("link", "parameter"), nsamps = NULL, verbose = FALSE, ...)
{
  family <- attr(object, "family")
  if(missing(newdata))
    newdata <- model.frame(object)
  if(is.character(newdata)) {
    if(file.exists(newdata <- path.expand(newdata)))
      newdata <- read.table(newdata, header = TRUE, ...)
  }
  if(is.matrix(newdata) || is.list(newdata))
    newdata <- as.data.frame(newdata)  
  if(!is.null(attr(object, "fixed.names")))
    names(newdata) <- rmf(names(newdata))
  nn <- names(newdata)

  object <- get.model(object, model)
  if(any(c("effects", "param.effects") %in% names(object)))
    object <- list(object)
  k <- length(object)
  enames <- list()
  for(j in 1:k)
    enames[[j]] <- all.terms(object[[j]], j, ne = TRUE, id = FALSE, intercept = FALSE)
  if(!all(diff(sapply(enames, length)) == 0))
    stop("the number of terms in the models is not identical, cannot compute prediction!")
  enames <- data.frame(enames)
  if(!all(apply(enames, 1, function(x) length(unique(x))) == 1))
    stop("different terms in the supplied models, cannot compute prediction!")
  enames <- all.terms(object[[1L]], ne = TRUE, id = FALSE, intercept = FALSE)
  term <- if(!is.null(enames)) {
    if(is.null(term)) enames else {
      if(is.character(term)) {
        unlist(sapply(term, function(x) {
          enames[grepl(gsub("[[:space:]]", "", x), enames, fixed = TRUE)]
        }))
      } else enames[term]
    }
  } else NULL
  term <- term[!is.na(term)]
  if(!length(term)) term <- NULL
  if(!is.null(term) & verbose)
    cat("terms used for prediction:", paste(term, collapse = ", "), "\n")
  rval <- NULL
  m.samples <- m.designs <- m.specials <- list()
  for(j in 1:k) {
    hi <- if(!is.null(object[[j]]$param.effects)) {
      any(grepl("(Intercept)", rownames(object[[j]]$param.effects), fixed = TRUE))
    } else FALSE
    if(!is.null(term)) {
      for(i in term) {
        specs <- attr(object[[j]]$effects[[i]], "specs")
        if(!is.null(specs)) {
          if(!all(specs$term %in% nn))
            stop(paste("cannot find variables", specs$term, "in newdata!"))
          if(is.null(specs$is.factor)) specs$is.factor <- FALSE
          tmp <- attr(object[[j]]$effects[[i]], "samples")
          if(is.null(dim(tmp))) {
            tmp <- matrix(tmp, ncol = 1)
          } else tmp <- as.matrix(tmp)
          if(!is.null(specs$special)) {
            m.specials[[i]] <- list("X" = PredictMat(specs, newdata), ## FIXME: also allow basis()?
              "get.mu" = specs$get.mu, "samples" = tmp)
          } else {
            m.samples[[i]] <- rbind(m.samples[[i]], tmp)
            if(inherits(object[[j]]$effects[[i]], "linear.bamlss")) {
              if(specs$is.factor & !is.character(newdata)) {
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
              }
            }
            attr(m.samples[[i]], "is.factor") <- specs$is.factor
          }
        } else {
          if(any(grepl(i, rownames(object[[j]]$param.effects), fixed = TRUE))) {
            ij <- which(colnames(attr(object[[j]]$param.effects, "samples")) %in% i)
            m.samples[[i]] <- cbind(m.samples[[i]],
              matrix(attr(object[[j]]$param.effects, "samples")[, ij, drop = FALSE],
              ncol = length(ij)))
            if(is.null(newdata[[i]])) {
              enames2 <- all.terms(object[[1L]], ne = FALSE, id = FALSE,
                intercept = FALSE, what = "parametric")
              tte <- NULL
              for(en in enames2) {
                if(grepl(rmf(en), rmf(i), fixed = TRUE))
                  tte <- en
              }
              if(is.null(tte)) stop(paste("cannot find term", i, "in newdata!"))
              if(is.null(newdata[[tte]])) {
                f <- as.formula(paste("~", if(hi) "-1" else "1", "+", tte))
                tmm <- model.matrix(f, data = newdata)
                m.designs[[i]] <- tmm[, grep(rmf(i), rmf(colnames(tmm)), fixed = TRUE)]
              } else {
                if(is.factor(newdata[[tte]])) {
                  nl <- nlevels(newdata[[tte]]) - if(hi) 1 else 0
                  f <- as.formula(paste("~", if(hi) "1" else "-1", "+", tte))
                  tmm <- model.matrix(f, data = newdata)
                  m.designs[[i]] <- tmm[, i, drop = FALSE]
                } else m.designs[[i]] <- newdata[[tte]]
              }
              if(length(m.designs)) {
                if(is.null(dim(m.designs[[i]]))) m.designs[[i]] <- matrix(m.designs[[i]], ncol = 1)
              }
            } else {
              if(is.factor(newdata[[i]])) {
                nl <- nlevels(newdata[[i]]) - if(hi) 1 else 0
                if(nl != ncol(m.samples[[i]]))
                  stop(paste("levels of factor variable", i, "in newdata not identical to model levels!"))
                f <- as.formula(paste("~", if(hi) "1" else "-1", "+", i))
                m.designs[[i]] <- model.matrix(f, data = newdata)
                if(hi) m.designs[[i]] <- m.designs[[i]][, -1]
              } else m.designs[[i]] <- newdata[[i]]
            }
          }
        }
        if(!is.null(nsamps) & length(m.samples)) {
          if(!is.null(dim(m.samples[[i]])))
            m.samples[[i]] <- m.samples[[i]][seq.int(1:ncol(m.samples[[i]]), length = nsamps), , drop = FALSE]
        }
        if(!is.null(m.designs[[i]])) {
          if(is.null(dim(m.designs[[i]])))
            m.designs[[i]] <- matrix(m.designs[[i]], ncol = 1)
        }
      }
    }

    if(intercept & j < 2) {
      sami <- attr(object[[j]]$param.effects, "samples")
      if(!is.null(sami)) {
        sami <- sami[, grep("(Intercept)", colnames(sami), fixed = TRUE)]
        if(!is.null(nsamps)) sami <- sami[seq.int(1:length(sami), length = nsamps)]
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
      sn <- sapply(m.samples, function(x) {
        if(is.null(dim(x))) length(x) else nrow(x)
      })
      if(!all(sn == max(sn))) {
        sn <- max(sn)
        m.samples <- lapply(m.samples, function(x) {
          if(is.null(dim(x))) c(x, rep(NA, sn - length(x))) else rbind(x, matrix(NA, sn - nrow(x), ncol(x)))
        })
        warning("different numbers of samples of model terms, please see also argument accept.only!")
      }
      m.samples <- as.data.frame(m.samples)
      m.designs <- as.data.frame(m.designs)
      options("warn" = warn)
      get.mu <- function(X, b, ...) {
        as.matrix(X) %*% as.numeric(b)
      }
      rval <- apply(m.samples, 1, function(x) { get.mu(m.designs, x, expand = FALSE) })
    } else rval <- 0
    if(length(m.specials)) {
      for(i in m.specials) {
        rval <- rval + apply(i$samples, 1, function(x) { i$get.mu(i$X, x, expand = FALSE) })
      }
    }
    type <- match.arg(type)
    if(type != "link") {
      link <- family$links[if(!is.null(model)) grep(model, names(family$links)) else 1]
      if(length(link) > 0) {
        linkinv <- make.link2(link)$linkinv
        rval <- t(apply(rval, 1, linkinv))
      } else {
        warning(paste("could not compute predictions on the scale of parameter",
          model, ", predictions on the scale of the linear predictor are returned!", sep = ""))
      }
    }
    if(!is.null(trans)) {
      rval <- apply(rval, MARGIN, trans)
      if(MARGIN < 2)
        rval <- t(rval)
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
## (8) Creating new smooth terms. ##
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
      rval$get.mu <- eval(parse(text = fun))
      rval$rs.by <- rs.by
      rval$by.vars <- vars
      rval$by.formula <- object$by.formula
      rval$one <- attr(ft, "intercept")
    }
  } else {
    rval$get.mu <- function(X, g, ...) {
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
  if(object$by != "NA") {
    by <- data[[object$by]]
    if(!is.factor(by))
      by <- as.factor(data[[object$by]])
    object$by.levels <- levels(by)
    object$fid <- as.integer(by)
    object$byname <- object$by
    object$by <- "NA"
    object$get.mu <- function(X, g, ...) {
      (g[4] + g[1]) * exp(-(g[5] + g[2]) * exp(-(g[6] + g[3]) * X))
    }
  } else {
    object$get.mu <- function(X, g, ...) {
      f <- g[1] / (1 + exp(g[2]) * (exp(g[3]) / (1 + exp(g[3])))^(drop(X)))
      if(object$xt$center)
        f <- f - mean(f)
      f
    }
    object$update <- update_optim2
    object$propose <- propose_slice
    object$prior <- function(gamma, tau2 = NULL) {
      sum(dnorm(gamma, sd = 1000, log = TRUE))
    }
    object$grad <- FALSE
    object$edf <- function(...) { 3 }
    object$fixed <- TRUE
    object$np <- 3
    object$p.save <- "g"
    object$state <- list("g" = rep(0, 3))
    object$s.colnames = c("c1", "c2", "c3")
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

    object$get.mu <- function(X, g, ...) {
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
    x$state$fit <- x$get.mu(x$X, x$state$g)

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


###################
## (9) Plotting. ##
###################
## Plotting method for "bamlss" objects.
plot.bamlss <- function(x, model = NULL, term = NULL, which = 1,
  ask = FALSE, scale = 1, spar = TRUE, ...)
{
  family <- attr(x, "family")
  args <- list(...)
  cx <- class(x)

  if(is.null(args$do_par) & spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  x <- get.model(x, model, drop = FALSE)

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

  if(all(which %in% c("hist-resid", "qq-resid", "scatter-resid", "scale-resid"))) {
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
        if(!any(c("effects", "param.effects") %in% names(x[[i]]))) {
          kn <- kn + get_k_n(x[[i]])
        } else {
          ne[[i]] <- if(!is.null(names(x[[i]]$effects))) names(x[[i]]$effects) else NA
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
          if(!is.null(x[[i]]$effects)) {
            kn[1] <- kn[1] + length(na.omit(pterms[[i]]))
          }
        }
      }
      kn
    }

    if(any(c("effects", "param.effects") %in% names(x)))
      x <- list(x)

    kn <- get_k_n(x)

    if(which == "effects" & kn[1] < 1) on.exit(warning("no terms to plot in model object!"), add = TRUE)

    if(is.null(args$do_par) & spar) {
      if(!ask) {
        if("cbamlss" %in% cx) {
          par(mfrow = c(length(x), kn[1] / length(x)))
        } else par(mfrow = n2mfrow(kn[if(which == "effects") 1 else 2]))
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
      args[c("x", "term", "which", "ask", "scale")] <- list(x[[i]], term, which, ask, scale)
      args$main <- if(!is.null(main)) main[i] else NULL
      if(!any(c("effects", "param.effects") %in% names(x[[i]]))) {
        args$do_par <- FALSE
        do.call("plot.bamlss", args)
      } else {
        args$mmain <- NULL
        do.call(".plot.bamlss", args)
      }
    }
  }

  invisible(NULL)
}

.plot.bamlss <- function(x, model = NULL, term = NULL, which = 1,
  ask = FALSE, scale = 1, spar = TRUE, ...)
{
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
    ne <- pterms <- list()
    for(i in 1:n) {
      ne[[i]] <- if(!is.null(names(x[[i]]$effects))) names(x[[i]]$effects) else NA
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
      if(!is.null(x[[i]]$effects) & length(na.omit(pterms[[i]]))) {
        k <- k + length(na.omit(pterms[[i]]))
        if(scale > 0) {
          term <- term[1:length(x[[i]]$effects)]
          for(e in pterms[[i]]) {
            et <- x[[i]]$effects[[e]]
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
    args$which <- which
    if(which != "effects")
      args$residuals <- NULL
    for(i in 1:n) {
      if(!is.null(x[[i]]$effects)) {
        for(e in pterms[[i]]) {
          lim <- c("ylim", "zlim")[(attr(x[[i]]$effects[[e]], "specs")$dim > 1) * 1 + 1]
          setlim <- FALSE
          if(!is.null(ylim) & is.null(args[[lim]])) {
            args[[lim]]<- ylim
            setlim <- TRUE
          }
          args$x <- x[[i]]$effects[[e]]
          do.call("plot.bamlss.effect", args)
          if(setlim) args[[lim]] <- NULL
        }
      }
    }
  }
  if(which == "param-samples") {
    for(i in 1:n) {
      args$main <- NULL
      args$x <- attr(x[[i]]$param.effects, "samples")
      do.call("plot", args)
    }
  }

  return(invisible(NULL))
}


## Generic plotting method for model terms.
plot.bamlss.effect <- function(x, which = "effects", ...) {
  if(which == "effects") {
    UseMethod("plot.bamlss.effect")
  } else {
    require("coda")
    args <- list(...)
    args$main <- NULL
    args$x <- attr(x, "samples")
    if(!is.null(attr(x, "samples.scale")))
      args$x <- as.mcmc(cbind(as.matrix(args$x), as.matrix(attr(x, "samples.scale"))))
    if(!is.null(attr(x, "samples.alpha")))
      args$x <- as.mcmc(cbind(as.matrix(args$x), as.matrix(attr(x, "samples.alpha"))))
    acf <- if(is.null(args$acf)) FALSE else args$acf
    args$acf <- NULL
    if(is.null(args$ask)) args$ask <- TRUE
    if(acf) {
      do.call("autocorr.plot", args)
    } else {
      do.call("plot", args)
    }
  }
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
      specs$get.mu(X, g, expand = FALSE)
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
       inherits(x, "mrf.smooth.spec") | inherits(x, "mrf.smooth")) {
      if(if(!is.null(args$density)) args$density else FALSE) {
        args$density <- NULL
        if(is.null(args$main))
          args$main <- attr(x, "specs")$label
        args$x <- density(x[, "50%"])
        if(!limNULL)
          args$xlim <- args$ylim
        do.call("plot", delete.args(stats:::plot.density, args, c("main", "xlim")))
      } else {
        if(!is.null(args$map)) {
          args$x <- x[, grepl("50%", colnames(x), fixed = TRUE)]
          args$id <- as.character(x[, 1])
          args$xlim <- args$ylim <- NULL
          do.call("plotmap", delete.args("plotmap", args,
            not = c("border", "lwd", "lty", names(formals("colorlegend")), "main")))
        } else do.call("plotblock", args)
      }
    } else {
      do.call("plot2d", delete.args("plot2d", args,
        c("xlim", "ylim", "pch", "main", "xlab", "ylab", "lwd", "axes")))
    }
  } else {
    if(is.null(args$c.select))
      args$c.select <- grep("50%", colnames(x), fixed = TRUE)
    if(!is.null(args$slice)) {
      do.call("sliceplot", delete.args("sliceplot", args,
        c("xlim", "ylim", "zlim", "main", "xlab", "ylab", "col", "lwd", "lty")))
    } else {
      do.call("plot3d", delete.args("plot3d", args,
        c("xlim", "ylim", "zlim", "pch", "main", "xlab", "ylab", "ticktype",
        "zlab", "phi", "theta", "r", "d", "scale", "range", "lrange", "pos", "image.map")))
    }
  }
}


###################################
## (10) Other helping functions. ##
###################################
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



###################################
## (11) Model summary functions. ##
###################################
summary.bamlss <- function(object, model = NULL, ...)
{
  call <- attr(object, "call")
  family <- attr(object, "family")
  object <- get.model(object, model)
  rval <- list()
  n <- length(object)
  for(i in 1:n) {
    if(!any(c("param.effects", "effects.hyp") %in% names(object[[i]]))) {
      rval[[i]] <- summary.bamlss(object[[i]])
      attr(rval[[i]], "hlevel") <- TRUE
    } else {
      for(j in c("param.effects", "effects.hyp")) {
        if(!is.null(object[[i]][[j]]))
          attr(object[[i]][[j]], "samples") <- NULL
      }
      rval[[i]] <- with(object[[i]],
        c(list("param.effects" = param.effects,
          "effects.hyp" = effects.hyp),
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
      cat("\nlog Lik. =", if(is.na(x$IC)) "NA" else {
            formatC(x$IC, digits = digits, flag = "-")
          }, "edf =", if(is.na(x$ed)) "NA" else {
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
      if(length(x[[i]]$param.effects) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x[[i]]$param.effects, digits = digits, na.print = "NA", ...)
      }
      if(length(x[[i]]$effects.hyp) > 0) {
        cat("\nSmooth effects variances:\n")
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
  on.exit(return(invisible(x)))
  xs <- summary(x)

  h0 <- !is.null(attr(x, "hlevel"))

  print_dic_pd <- function(x, ok = TRUE) {
    if(is.list(x)) {
      if(any(c("DIC", "pd") %in% names(x[[1]])))
        x <- x[[1]]
    }
    dp <- FALSE
    if(!is.null(x$DIC) & !is.null(x$pd)) {
      dp <- TRUE
      cat("---\n")
      cat("DIC =", if(is.na(x$DIC)) "NA" else {
          formatC(x$DIC, digits = digits, flag = "-")
        }, "pd =", if(is.na(x$pd)) "NA" else {
          formatC(x$pd, digits = digits, flag = "-")
        })
    }
    if(!is.null(x$N)) {
      if(!dp) cat("---\n")
      dp <- TRUE
      cat(" N =", if(is.na(x$N)) "NA" else formatC(x$N, digits = digits, flag = "-"))
    }
    if(ok & dp) cat("\n\n")
  }

  family <- attr(x, "family")
  n <- attr(xs, "n")
  nx <- NULL
  if(n < 2)
    xs <- list(xs)
  else
    nx <- names(xs)

  pdic <- TRUE
  
  cat("\n")
  print(if(is.function(family)) family() else family)
  cat("---\n")
  for(i in 1:n) {
    if(!is.null(nx)) {
      cat("Formula ", nx[i], ":\n", sep = "")
    } else cat("Formula:\n")
    if(!is.null(attr(xs[[i]], "hlevel"))) {
      nh <- names(xs[[i]])
      for(j in seq_along(xs[[i]])) {
        attr(xs[[i]][[j]]$formula, "name") <- NULL
        cat(nh[j], ": ", sep = ""); print(xs[[i]][[j]]$formula)
      }
      if(i < n) cat("---\n")
      if(i == n & pdic) {
        print_dic_pd(xs[[i]][[1]])
        pdic <- FALSE
      }
    } else {
      attr(xs[[i]]$formula, "name") <- NULL
      print(xs[[i]]$formula)
      if(i == n & pdic) {
        print_dic_pd(xs[[1]])
        pdic <- FALSE
      }
    }

    if(i < n & h0) cat("---\n")

    if(i == n & h0 & pdic) {
      print_dic_pd(xs[[1]])
    }
  }
}


####################################
## (12) More extractor functions. ##
####################################
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


#edf <- function(x, samples = TRUE, nsamps = NULL)
#{
#  family <- attr(x, "family")
#  if(is.null(family$weights))
#    stop("cannot compute edf without the weights() function in the model family object!")
#  nx <- family$names
#  edf0 <- 0; edf1 <- NULL
#  y <- model.response2(x)
#  eta <- fitted(x, type = "link", samples = samples, nsamps = nsamps)
#  mf <- model.frame(x)
#  for(j in 1:length(nx)) {
#    if(!is.null(x[[nx[j]]]$param.effects))
#      edf0 <- edf0 + ncol(x[[nx[j]]]$param.effects)
#    if(!is.null(x[[nx[j]]]$effects)) {
#      weights <- family$weights[[nx[j]]](y, eta)
#      for(sj in 1:length(x[[nx[j]]]$effects)) {
#        specs <- attr(x[[nx[j]]]$effects[[sj]], "specs")
#        X <- if(inherits(specs, "mgcv.smooth")) {
#          PredictMat(specs, mf)
#        } else {
#          if(!is.null(specs$basis)) {
#            stopifnot(is.function(specs$basis))
#            specs$basis(mf[, specs$term])
#          } else stop(paste("cannot compute design matrix for term ", specs$label, "!", sep = ""))
#        }
#      }
#    }
#  }
#}


## Extract model formulas.
formula.bamlss <- function(x, model = NULL, ...)
{
  if(all(c("model", "effects") %in% names(x))) {
    f <- x$model$formula
  } else {
    f <- attr(x, "formula")
    if(inherits(f, "list")) {
      if(length(f) < 2) {
        f <- f[[1]]
      } else {
        if(!is.null(model)) {
          for(j in model) {
            f <- f[[j]]
          }
        }
      }
    }
  }
  f
}

print.bamlss.formula <- function(x, ...) {
  if(!inherits(x, "list")) {
    print(x)
  } else {
    nx <- names(x)
    if(is.null(nx))
      nx <- as.character(1:length(x))
    for(i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if(inherits(x[[i]], "list") & is.null(attr(x, "raw.formula"))) {
        for(j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          attr(x[[i]][[j]], "name") <- NULL
          print(x[[i]][[j]])
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        if(is.null(attr(x, "raw.formula"))) print(x[[i]]) else print(x[[i]]$formula)
      }
      if(i < length(x))
      cat("\n")
    }
  }
  invisible(NULL)
}


## Extract formula terms.
terms.bamlss <- function(x, ...) {
  terms.bamlss.formula(formula(x, ...))
}

terms.bamlss.formula <- function(x)
{
  if(!inherits(x, "list")) {
    tx <- terms(x)
  } else {
    tx <- list()
    for(i in seq_along(x)) {
      if(inherits(x[[i]], "list")) {
        tx[[i]] <- terms.bamlss.formula(x[[i]])
        names(tx[[i]]) <- paste("h", 1:length(x[[i]]), sep = "")
      } else tx[[i]] <- terms(x[[i]])
    }
    names(tx) <- names(x)
    if(length(tx) < 2)
      tx <- tx[[1]]
  }
  tx
}


## Model extractor function.
get.model <- function(x, model = NULL, drop = TRUE)
{
  if(!drop)
    mf <- model.frame(x)
  if(length(model) > 1) {
    for(j in model)
      x <- get.model(x, j, drop = drop)
  } else {
    family <- attr(x, "family")
    cx <- class(x)
    elmts <- c("formula", "fake.formula", "model", "param.effects",
      "effects", "fitted.values", "residuals")
    if(!any(names(x) %in% elmts)) {
      if(!is.null(model)) {
        if(is.character(model)) {
          if(all(is.na(model <- pmatch(model, names(x)))))
            stop("argument model is specified wrong!")
        } else {
          if(max(model) > length(x) || is.na(model) || min(model) < 1) 
            stop("argument model is specified wrong!")
        }
        x <- if(drop) x[[model]] else {
          class(x) <- "list"
          x[model]
        }
      }
    } else x <- list(x)
    class(x) <- cx
    attr(x, "family") <- family
  }
  if(!drop)
    attr(x, "model.frame") <- mf

  return(x)
}


## Fitted values/terms extraction
fitted.bamlss <- function(object, model = NULL, term = NULL,
  type = c("link", "parameter"), samples = FALSE, FUN = mean,
  nsamps = NULL, ...)
{
  type <- match.arg(type)
  family <- attr(object, "family")

  if(type != "parameter" & !samples)
    object <- get.model(object, model, drop = FALSE)

  h1check <- any(grepl("h1", names(object)))

  elmts <- c("formula", "fake.formula", "model", "param.effects",
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
    for(j in seq_along(object)) {
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
              if(length(e <- grep(term[i], ne)))
                fe[[ne[e]]] <- object[[j]]$effects[[e]]
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
          if(is.null(object[[j]][[jj]]$param.effects) & is.null(object[[j]][[jj]]$effects)) next
          rval[[j]] <- rval[[j]] + predict.bamlss(object, model = c(j, jj), term = term,
            FUN = function(x) { x }, nsamps = nsamps, ...)
        }
      } else {
        if(is.null(object[[j]]$param.effects) & is.null(object[[j]]$effects)) next
        rval[[j]] <- predict.bamlss(object, model = if(one) NULL else j, term = term,
          FUN = function(x) { x }, nsamps = nsamps, ...)
      }
      if(type != "link" & !h1check)
        rval[[j]] <- apply(rval[[j]], 2, make.link2(family$links[if(one) 1 else nrval[j]])$linkinv)
      if(!h1check) {
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
samples <- function(x, ...) {
  UseMethod("samples")
}

grep2 <- function(pattern, x, ...) {
  i <- NULL
  for(p in pattern)
    i <- c(i, grep(p, x, ...))
  sort(i)
}

samples.bamlss <- function(x, model = NULL, term = NULL, ...)
{
  x <- get.model(x, model)
  if(!is.null(nx <- names(x))) {
    if(all(c("effects", "param.effects") %in% nx))
      x <- list(x)
  }
  if(is.null(nx))
    nx <- paste("p", 1:length(x), sep = "")
  args <- list(...)
  id <- args$id
  samps <- list(); k <- 1
  for(j in seq_along(x)) {
    if(!all(c("effects", "param.effects") %in% names(x[[j]]))) {
      samps[[j]] <- samples.bamlss(x[[j]], term = term, id = nx[j])
    } else {
      if(!is.null(x[[j]]$param.effects)) {
        if(nrow(x[[j]]$param.effects) > 0) {
          if(length(i <- grep2(term, rownames(x[[j]]$param.effects), fixed = TRUE))) {
            samps[[k]] <- attr(x[[j]]$param.effects, "samples")[, i, drop = FALSE]
            if(length(x) > 1) {
              colnames(samps[[k]]) <- paste(colnames(samps[[k]]),
                if(!is.null(id)) paste(id, ".", sep = ""), ":", nx[j], sep = "")
            }
            k <- k + 1
          }
        }
      }
      if(!is.null(x[[j]]$effects)) {
        if(length(i <- grep2(term, names(x[[j]]$effects), fixed = TRUE))) {
          for(e in i) {
            samps[[k]] <- attr(x[[j]]$effects[[e]], "samples")
            if(length(x) > 1) {
              colnames(samps[[k]]) <- paste(colnames(samps[[k]]),
                if(!is.null(id)) paste(id, ".", sep = ""), ":", nx[j], sep = "")
            }
            k <- k + 1
          }
        }
      }
    }
  }
  if(k < 2) stop("no samples could be extracted, please check term and model selection!")
  n <- max(sapply(samps, length))
  samps <- lapply(samps, function(sm) {
    if(nrow(sm) < n) {
      NAs <- matrix(NA, nrow = n - nrow(sm), ncol = ncol(sm))
      colnames(NAs) <- colnames(sm)
      sm <- rbind(sm, NAs)
    }
    sm
  })
  samps <- if(length(samps) < 1) NULL else as.mcmc(do.call("cbind", samps))
  samps <- samps[, unique(colnames(samps)), drop = FALSE]
  samps
}


## Credible intervals of coefficients.
confint.bamlss <- function(object, parm, level = 0.95, model = NULL, ...)
{
  args <- list(...)
  if(!is.null(args$term))
    parm <- args$term
  if(missing(parm))
    parm <- all.terms(object, ne = TRUE, id = FALSE)
  samps <- samples(object, model = model, term = parm)
  np <- colnames(samps)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  apply(samps, 2, quantile, probs = probs)
}


## Extract model coefficients.
coef.bamlss <- function(object, model = NULL, term = NULL, FUN = mean, ...)
{
  object <- get.model(object)
  if(is.null(term))
    term <- all.terms(object, ne = TRUE, id = FALSE)
  samps <- samples(object, model = model, term = term)
  apply(samps, 2, function(x) { FUN(na.omit(x), ...) })
}


## Get all terms names used.
all.terms <- function(x, model = NULL, ne = TRUE, what = c("parametric", "smooth"), ...)
{
  what <- match.arg(what, several.ok = TRUE)
  args <- list(...)
  nx <- names(x)
  if(!ne) {
    tx <- terms.bamlss(x, model)
    if(!inherits(tx, "list"))
      tx <- list(tx)
    tl <- NULL
    for(j in tx) {
      if(inherits(j, "list")) {
        for(i in j) {
          if(attr(i, "intercept"))
            tl <- c(tl, "(Intercept)")
          tl <- c(tl, attr(i, "term.labels"))
        }
      } else {
        if(attr(j, "intercept"))
          tl <- c(tl, "(Intercept)")
        tl <- c(tl, attr(j, "term.labels"))
      }
    }
  } else {
    x <- get.model(x, model)
    if(all(c("effects", "param.effects") %in% names(x))) {
      x <- list(x)
      if(!is.null(model)) {
        if(!is.character(model))
          nx <- nx[as.integer(model)]
        else nx <- grep(model, nx, value = TRUE)
      }
    }
    tl <- NULL
    for(j in seq_along(x)) {
      if(!"effects" %in% names(x[[j]]))
        tl <- c(tl, all.terms(x[[j]], ne = TRUE, ...))
      else {
        if(!is.null(x[[j]]$param.effects)) {
          tl <- c(tl, paste(rownames(x[[j]]$param.effects),
            if(is.null(args$id)) paste(":", nx[j], sep = ""), sep = ""))
        }
        tl <- c(tl, names(x[[j]]$effects))
      }
    }
  }
  if(!is.null(args$intercept))
    tl <- tl[!grepl("Intercept", tl)]
  if(!("smooth" %in% what)) {
    specials <- attr(formula(x), "specials")
    for(sp in specials) {
      if(any(dte <- grepl(paste(sp, "(", sep = ""), tl, fixed = TRUE)))
        tl <- tl[!dte]
    }
  }
  if(!("parametric" %in% what)) {
    specials <- attr(formula(x), "specials")
    for(sp in specials) {
      if(any(dte <- grepl(paste(sp, "(", sep = ""), tl, fixed = TRUE)))
        tl <- tl[dte]
    }
  }

  tl
}


## Get the model.frame.
model.frame.bamlss <- function(formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- if(length(nargs) || is.null(attr(formula, "model.frame"))) {
    fcall <- attr(formula, "call")
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(bamlss.model.frame)
    fcall[names(nargs)] <- nargs
    env <- environment(attr(formula, "formula"))
    if(is.null(env)) 
      env <- parent.frame()
    eval(fcall, env)
  } else attr(formula, "model.frame")
  if(!is.null(attr(mf, "orig.names")))
    names(mf) <- attr(mf, "orig.names")
  mf
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
model.response2 <- function(data, ...)
{
  if(!inherits(data, "data.frame"))
    data <- model.frame(data)
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
        rval <- formula_respname(x$model$formula)
    }
  }
  rval
}


## Create the inverse of a matrix.
matrix_inv <- function(x)
{
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
  return(p)
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


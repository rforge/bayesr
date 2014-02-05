################################################
## (0) BayesR main model fitting constructor. ##
################################################
## Could be interesting: http://people.duke.edu/~neelo003/r/
##                       http://www.life.illinois.edu/dietze/Lectures2012/
xreg <- function(formula, family = gaussian.BayesR, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  reference = NULL, parse.input = parse.input.bayesr, transform = transformJAGS,
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
  } else "parse.input.bayesr"

  ## Parse input.
  pm <- match.call(expand.dots = FALSE)
  pm$parse.input <- pm$setup <- pm$samples <- pm$results <- NULL
  pm[[1]] <- as.name(functions$parse.input)
  pm <- eval(pm, parent.frame())

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
  attr(rval, "formula") <- attr(pm, "formula0")
  attr(rval, "call") <- match.call()

  rval
}


#########################
## (1) BayesR wrapper. ##
#########################
bayesr <- function(formula, family = gaussian, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  engine = c("IWLS", "BayesX", "JAGS", "STAN"), cores = NULL, combine = TRUE,
  n.iter = 12000, thin = 10, burnin = 2000, seed = NULL, ...)
{
  xengine <- match.arg(engine)
  family <- deparse(substitute(family), backtick = TRUE, width.cutoff = 500)

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

  if(xengine == "IWLS") {
    transform <- transformIWLS
    setup <- FALSE
    engine <- function(x) {
      samplerIWLS(x, n.iter = n.iter, thin = thin,
        burnin = burnin, seed = seed, ...)
    }
    results <- resultsIWLS
  }
  
  if(xengine %in% c("JAGS", "STAN")) {
    transform <- transformJAGS
    if(xengine == "JAGS") {
      require("rjags")
      setup <- setupJAGS
      engine <- function(x) {
        samplerJAGS(x, n.iter = n.iter, thin = thin,
          burnin = burnin, seed = seed, ...)
      }
    } else {
      require("rstan")
      setup <- jags2stan
      engine <- function(x) {
        samplerSTAN(x, n.iter = n.iter, thin = thin,
          burnin = burnin, seed = seed, ...)
      }  
    }

    results <- resultsJAGS
  }

  rval <- xreg(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transform,
    setup = setup, engine = engine, results = results, cores = cores,
    combine = combine, sleep = 1, ...)
  
  attr(rval, "call") <- match.call()
  
  rval
}


##########################################################
## (2) Parsing all input using package mgcv structures. ##
##########################################################
parse.input.bayesr <- function(formula, data = NULL, family = gaussian.BayesR,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, knots = NULL, specials = NULL, reference = NULL,
  grid = 100, ...)
{
  ## Search for additional formulas
  formula2 <- NULL
  fn <- names(fo <- formals(fun = parse.input.bayesr))[-1]
  fn <- fn[fn != "..."]
  for(f in fn) {
    fe <- eval(parse(text = f))
    if(is.character(fe)) {
      if(grepl("~", fe, fixed = TRUE))
        fe <- as.formula(fe)
    }
    if(inherits(fe, "formula")) {
      formula2 <- c(formula2, fe)
      eval(parse(text = paste(f, if(is.null(fo[[f]])) "NULL" else fo[[f]], sep = " = ")))
    }
  }

  ## Parse family object.
  family <- bayesr.family(family)

  ## Parse formula
  formula <- bayesr.formula(c(formula, formula2), specials, family)
  formula0 <- attr(formula, "formula0")

  ## Create the model frame.
  mf <- bayesr.model.frame(formula, data, family, weights,
    subset, offset, na.action, specials)
  response.name <- attr(mf, "response.name")

  ## For categorical responses, extend formula object.
  ylevels <- NULL
  if(length(response.name) < 2) {
    if(is.factor(mf[[response.name]])) {
      cat <- if(is.null(family$cat)) FALSE else family$cat
      if(cat & nlevels(mf[[response.name]]) > 1) {
        if(is.null(reference)) {
          ty <- table(mf[[response.name]])
          reference <- c(names(ty)[ty == max(ty)])[1]
        }
        reference <- rmf(reference)
        ylevels <- rmf(levels(mf[[response.name]]))
        ylevels <- ylevels[ylevels != reference]
        if(length(formula) != (n <- length(ylevels)))
          formula <- rep(formula, length.out = n)
        names(formula) <- nf <- paste(response.name, ylevels, sep = "")
        for(j in seq_along(formula)) {
          uf <- eval(parse(text = paste(nf[j], " ~ .", sep = "")))
          if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
            formula[[j]][[1]]$cat.formula <- update(formula[[j]][[1]]$formula, uf)
          } else formula[[j]]$cat.formula <- update(formula[[j]]$formula, uf)
        }
      }
    } else mf[[response.name]] <- as.numeric(mf[[response.name]])
  } else {
    for(y in response.name)
      mf[[y]] <- as.numeric(mf[[y]])
  }

  ## Assign all design matrices and the hierarchical level, if any.
  rval <- bayesr.design(formula, mf, contrasts, knots, ...)
  rval <- bayesr.hlevel(rval)

  attr(rval, "family") <- family
  attr(rval, "reference") <- reference
  attr(rval, "ylevels") <- ylevels
  attr(rval, "grid") <- grid
  attr(rval, "model.frame") <- mf
  attr(rval, "formula0") <- formula0

  class(rval) <- c("bayesr.input", "list")

  rval
}

"[.bayesr.input" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  xattr <- attributes(x)
  mostattributes(rval) <- attributes(x)
  rval
}

"[.bayesr" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  mostattributes(rval) <- attributes(x)
  rval
}


## Assign all designs matrices.
bayesr.design <- function(x, data, contrasts = NULL, knots = NULL, ...)
{
  assign.design <- function(obj, mf) {
    if(!all(c("formula", "fake.formula", "response") %in% names(obj)))
      return(obj)
    pf <- paste("~ ", if(obj$intercept) 1 else -1,
      if(length(obj$pterms)) paste(" +", paste(obj$pterms, collapse = " + ")), sep = "")
    obj$X <- model.matrix(as.formula(pf),
      data = mf, contrasts.arg = contrasts, ...)
    obj$pterms <- colnames(obj$X)
    rn <- obj$response
    obj$response.vec <- if(!is.null(rn)) mf[[rn]] else NULL
    if(length(obj$smooth)) {
      smooth <- list()
      for(j in obj$smooth) {
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
      smooth <- gam.side(smooth, obj$X, tol = .Machine$double.eps^.5)
      sme <- mgcv:::expand.t2.smooths(smooth)
      if(is.null(sme)) {
        original.smooth <- NULL
      } else {
        original.smooth <- smooth
        smooth <- sme
        rm(sme)
      }
      obj$smooth <- smooth
    }
    if(length(obj$sx.smooth)) {
      sx.smooth <- list()
      for(j in obj$sx.smooth)
        sx.smooth <- c(sx.smooth, list(eval(parse(text = j))))
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
bayesr.model.frame <- function(formula, data, family, weights = NULL,
  subset = NULL, offset = NULL, na.action = na.omit, specials = NULL)
{
  family <- bayesr.family(family)
  formula <- bayesr.formula(formula, specials, family)

  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  if(missing(data))
    data <- environment(formula)
  if(is.matrix(data))
    data <- as.data.frame(data)

  ## Make fake "Formula" object.
  fF <- make_fFormula(formula)

  ## Set up the model.frame.
  mf <- list(formula = fF, data = data, weights = weights,
    subset = subset, offset = offset, na.action = na.action,
    drop.unused.levels = TRUE)
  mf <- do.call("model.frame", mf)

  ## Remove inf values
  mf <- rm_infinite(mf)

  ## assign response names
  attr(mf, "response.name") <- rn <- all.vars(formula(fF, rhs = 0))

  ## Check response.
  if(!is.null(family$valid.response)) {
    family$valid.response(mf[, rn])
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
          warning("infinite values in data, removing these observations in model frame!")
          x <- x[is.finite(x[, j]), ]
        }
      }
    }
  }
  x
}


## Parse families and get correct family object, depending on type.
bayesr.family <- function(family, type = "BayesR")
{
  family <- if(is.function(family)) family() else {
    if(is.character(family)) {
      if(!is.null(type)) {
        if(!grepl(type, family))
          family <- paste(family, type, sep = ".")
      }
      family <- eval(parse(text = family))
      if(is.function(family))
        family()
      else family
    } else family
  }
  if(!inherits(family, "family.BayesR")) {
    if(!is.character(family)) {
      if(is.null(family$family)) stop("family is specidied wrong, no family name available!")
      family <- family$family
    }
    family <- eval(parse(text = paste(family, type, sep = if(!is.null(type)) "." else "")))
    family <- family()
  }
  if(is.null(family$cat))
    family$cat <- FALSE
  if(is.null(family$mu))
    family$mu <- function(x) { x }
  family
}


## Special formula parser, can deal with multi parameter models
## and hierarchical structures.
bayesr.formula <- function(formula, specials, family)
{
  if(inherits(formula, "bayesr.formula"))
    return(formula)

  env <- environment(formula)
  if(is.null(env)) env <- .GlobalEnv
  if(!is.list(formula)) formula <- list(formula)
  if(!length(formula)) stop("formula is specified wrong!")
  names(formula) <- family$names[1:length(formula)]

  complete_formula <- function(formula) {
    if(any(i <- is.na(names(formula))))
      names(formula)[i] <- family$names[i]
    for(j in family$names) {
      if(is.null(formula[[j]])) {
        formula[[j]] <- as.formula(paste(j, "1", sep = " ~ "))
        environment(formula[[j]]) <- env
      }
    }
    formula
  }

  formula <- formula_and(formula, env)
  formula <- formula_at(formula, env)
  formula <- formula0 <- complete_formula(formula_hierarchical(formula))
  formula <- formula_extend(formula, specials, family)

  environment(formula) <- environment(formula0) <- env
  class(formula) <- class(formula0) <- c("bayesr.formula", "list")
  for(j in seq_along(formula0)) {
    if(!inherits(formula0[[j]], "formula")) {
      if(is.null(names(formula0[[j]])))
        names(formula0[[j]]) <- paste("h", 1:length(formula0[[j]]), sep = "")
    }
  }
  attr(formula, "formula0") <- formula0

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
    tl <- attr(mt, "term.labels")
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
        smj <- eval(parse(text = j))
        sterms <- c(sterms, smj$term, smj$by)
        if(!is.null(smj$by.formula)) {
          sterms <- c(sterms, attr(terms(smj$by.formula, keep.order = TRUE), "term.labels"))
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
  if(length(formula) < 3) NA else deparse(formula[[2]])
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

bayesr.hlevel <- function(x)
{
  if(!all(c("formula", "fake.formula", "response") %in% names(x))) {
    for(j in seq_along(x)) {
      if(!all(c("formula", "fake.formula", "response") %in% names(x[[j]]))) {
        for(i in seq_along(x[[j]])) {
          x[[j]][[i]]$hlevel <- i
        }
      } else x[[j]]$hlevel <- 1
    }
  } else x$hlevel <- 1
  x
}


###########################
## (3) Utility functions ##
###########################
## Transform smooth terms to mixed model representation.
randomize <- function(x)
{
  if(inherits(x, "bayesr.input") & !any(c("smooth", "response") %in% names(x))) {
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


## Combine method for "bayesr" objects.
c.bayesr <- function(...)
{
  objects <- list(...)
  x <- NULL
  for(i in 1L:length(objects))
    x <- c(x, objects[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- c("bayesr", "cbayesr")

  return(x)
}


## Function to compute statistics from samples of a model term.
compute_term <- function(x, get.X, get.mu, psamples, vsamples = NULL,
  asamples = NULL, FUN = NULL, snames, effects.hyp, fitted.values, data,
  grid = 100, rug = TRUE, hlevel = 1)
{
  require("coda")

  ## Data for rug plotting.
  rugp <- if(length(x$term) < 2 & rug) data[[x$term]] else NULL

  if(hlevel > 1)
    data <- unique(data)

  ## New x values for which effect should
  ## be calculated, n = 100.
  if(!is.na(grid)) {
    if(length(x$term) < 2 & !is.factor(data[[x$term[1]]]) & !any(grepl("mrf", class(x))) &
      !any(grepl("re.", class(x), fixed = TRUE)) & !any(grepl("random", class(x)))) {
      xsmall <- TRUE
      nd <- list()
      for(j in x$term) {
        xr <- range(data[[j]], na.rm = TRUE)
        nd[[j]] <- seq(xr[1], xr[2], length = grid)
      }
      if(x$by != "NA") { ## FIXME: check by variables!
        if(!is.factor(data[[x$by]])) {
          xr <- range(data[[x$by]], na.rm = TRUE)
          nd[[x$by]] <- seq(xr[1], xr[2], length = 100)
        } else nd[[x$by]] <- data[[x$by]]
      }
      data0 <- data
      data <- as.data.frame(nd)
    } else xsmall <- FALSE
  } else {
    data0 <- data[, c(x$term, if(x$by != "NA") x$by else NULL), drop = FALSE]
    data <- unique(data0)
    xsmall <- if(nrow(data) != nrow(data0)) TRUE else FALSE
  }

  X <- get.X(data)

  ## Compute samples of fitted values.
  fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })

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
  nt <- length(x$term)
  for(l in nt:1) {
    smf <- cbind(data[[x$term[l]]], smf)
  }
  names(smf) <- c(x$term, cnames)

  ## Compute new linear predictor.
  fit <- rep(0, nrow(smf))
  if(any(im <- grepl("50%", tolower(colnames(smf)), fixed = TRUE))) {
    im <- c(1:ncol(smf))[im]
    fit <- smf[, im[1]]
    if(xsmall)
      fit <- approx(data[[x$term]], fit, xout = data0[[x$term]])$y
  }
  by.drop <- NULL
  if(x$by != "NA") {
    by.drop <- (if(xsmall) data0[[x$by]] else data[[x$by]]) == x$by.level
    fit[!by.drop] <- 0
    if(!xsmall)
      smf <- smf[by.drop, ]
  }

  fitted.values <- if(!is.null(fitted.values)) fitted.values + fit else fit

  ## Assign class and attributes.
  smf <- unique(smf)
  if(is.factor(data[, x$term])) {
    bbb <- 1 ## FIXME: factors!
  }
  class(smf) <- c(class(x), "data.frame")
  x["X"] <- NULL
  attr(smf, "specs") <- x
  attr(smf, "specs")[c("X", "Xf", "rand", "trans.D", "trans.U")] <- NULL
  class(attr(smf, "specs")) <- class(x)
  attr(smf, "fit") <- fit
  attr(smf, "x") <- if(xsmall) data0[, x$term] else data[, x$term]
  attr(smf, "by.drop") <- by.drop
  attr(smf, "rug") <- rugp
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
      smat <- matrix(c(me, sd, qu), nrow = 1)
      colnames(smat) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
      rownames(smat) <- x$label
      if(!is.null(asamples)) {
        smat <- cbind(smat, "alpha" = mean(asamples))
      }
      smatfull <- rbind(smatfull, smat)
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
        e <- e - mean(e)
      } else {
        if(attr(effects[[i]], "specs")$xt$center)
          e <- e - mean(e)
      }
      e <- if(is.factor(attr(effects[[i]], "x"))) {
        warn <- getOption("warn")
        options(warn = -1)
        tx <- as.integer(as.character(attr(effects[[i]], "x")))
        options("warn" = warn)
        cbind(if(!any(is.na(tx))) tx else as.integer(attr(effects[[i]], "x")), e)
      } else cbind(attr(effects[[i]], "x"), e)
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


#####################
## (4) Prediction. ##
#####################
## A prediction method for "bayesr" objects.
## Prediction can also be based on multiple chains.
predict.bayesr <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = mean, trans = NULL, type = c("link", "parameter"),
  nsamps = NULL, ...)
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
          }
          if(!is.null(specs$special)) {
            m.specials[[i]] <- list("X" = PredictMat(specs, newdata), ## FIXME: also allow basis()?
              "get.mu" = specs$get.mu, "samples" = tmp)
          } else {
            m.samples[[i]] <- rbind(m.samples[[i]], tmp)
            if(inherits(object[[j]]$effects[[i]], "linear.bayesx")) {
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
                    specs$basis(newdata[specs$term])
                  } else stop(paste("cannot compute design matrix for term ", specs$label, "!", sep = ""))
                }
              }
            }
            attr(m.samples[[i]], "is.factor") <- specs$is.factor
          }
        } else {
          if(any(grepl(i, rownames(object[[j]]$param.effects), fixed = TRUE))) {
            ij <- grep(i, colnames(attr(object[[j]]$param.effects, "samples")), fixed = TRUE)
            m.samples[[i]] <- cbind(m.samples[[i]],
              matrix(attr(object[[j]]$param.effects, "samples")[, ij, drop = FALSE],
              ncol = length(ij)))
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
        if(!is.null(nsamps)) {
          m.samples[[i]] <- m.samples[[i]][seq.int(1:ncol(m.samples[[i]]), length = nsamps), , drop = FALSE]
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
    type <- match.arg(type)
    if(type != "link") {
      link <- family$links[grep(model, names(family$links))]
      if(length(link) > 0) {
        linkinv <- make.link2(link)$linkinv
        rval <- t(apply(rval, 1, linkinv))
      } else {
        warning(paste("could not compute predictions on the scale of parameter",
          model, ", predictions on the scale of the linear predictor are returned!", sep = ""))
      }
    }
    if(!is.null(trans))
      rval <- t(apply(rval, 1, trans))
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
    rval$get.mu <- function(X, g) {
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


## Rational spline constructor.
rs <- function(...)
{
  rval <- s(...)
  rval$class <- class(rval)
  rval$special <- TRUE
  class(rval) <- "rs.smooth.spec"
  rval
}

smooth.construct.rs.smooth.spec <- function(object, data, knots) 
{
  class(object) <- object$class
  acons <- FALSE
  if(!is.null(object$xt$center))
    acons <- object$xt$center
  object$xt$center <- acons
  object$fixed <- TRUE
  if(!is.null(object$xt$fixed))
    object$fixed <- object$xt$fixed
  object$xt$fixed <- object$fixed
  object <- smoothCon(object, data, knots, absorb.cons = acons)[[1]]
  object$by.done <- TRUE
  object$get.mu <- function(X, g) {
    k <- ncol(X)
    w <- c(1, g[(k + 1):(2 * k - 1)])
    g <- g[1:k]
    drop(X %*% g) / exp(drop(X %*% w))
  }
  object$class <- class(object)
  class(object) <- if(object$fixed) c("rs.smooth", "no.mgcv") else "rs.smooth"
  object
}

Predict.matrix.rs.smooth <- function(object, data, knots)
{
  class(object) <- object$class
  Predict.matrix(object, data) 
}


###################
## (6) Plotting. ##
###################
## Plotting method for "bayesr" objects.
plot.bayesr <- function(x, model = NULL, term = NULL, which = 1,
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
    res0 <- do.call("residuals.bayesr", delete.args("residuals.bayesr", args2))
    ny <- if(is.null(dim(res0))) 1 else ncol(res0)
    if(is.null(args$do_par) & spar) {
      if(!ask) {
        par(mfrow = n2mfrow(length(which) * ny))
      } else par(ask = ask)
    }
    if(any(which %in% c("scatter-resid", "scale-resid"))) {
      fit0 <- fitted.bayesr(x, type = "parameter", samples = TRUE,
        model = if(ny < 2) 1 else NULL, nsamps = args$nsamps)
    }
    rtype <- args$type
    if(is.null(rtype)) rtype <- "quantile"
    for(j in 1:ny) {
      res <- if(ny > 1) res0[, j] else res0
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
          if(is.null(args$xlab)) args2$xlab <- "Residuals"
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
          args2$y <- if(rtype == "quantile") qnorm(res) else (res - mean(res)) / sd(res)
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
          args2$x <- fit
          args2$y <- res
          args2 <- delete.args("scatter.smooth", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$xlab)) args2$xlab <- "Fitted values"
          if(is.null(args$xlab)) args2$ylab <- "Residuals"
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
          args2$x <- fit
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
        if("cbayesr" %in% cx) {
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
        do.call("plot.bayesr", args)
      } else do.call(".plot.bayesr", args)
    }
  }

  invisible(NULL)
}

.plot.bayesr <- function(x, model = NULL, term = NULL, which = 1,
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
          do.call("plot.bayesr.effect", args)
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
plot.bayesr.effect <- function(x, which = "effects", ...) {
  if(which == "effects") {
    UseMethod("plot.bayesr.effect")
  } else {
    require("coda")
    args <- list(...)
    args$main <- NULL
    args$x <- attr(x, "samples")
    if(!is.null(attr(x, "samples.scale")))
      args$x <- as.mcmc(cbind(as.matrix(args$x), as.matrix(attr(x, "samples.scale"))))
    if(!is.null(attr(x, "samples.alpha")))
      args$x <- as.mcmc(cbind(as.matrix(args$x), as.matrix(attr(x, "samples.alpha"))))
    par(ask = TRUE)
    do.call("plot", args)
  }
}


## Default model term plotting method.
plot.bayesr.effect.default <- function(x, ...) {
  args <- list(...)
  args$x <- x

  lim <- c("ylim", "zlim")[(attr(x, "specs")$dim > 1) * 1 + 1]
  if(is.null(args[[lim]])) {
    args[[lim]] <- range(x[, c("2.5%", "97.5%")], na.rm = TRUE)
    if(!is.null(args$residuals)) {
      if(args$residuals & !is.null(attr(x, "partial.resids")))
        args[[lim]] <- range(c(args[[lim]], attr(x, "partial.resids")[, -1]), na.rm = TRUE)
    }
  }
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
        do.call("plot", delete.args(stats:::plot.density, args, "main"))
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
        c("xlim", "ylim", "pch", "main", "xlab", "ylab")))
    }
  } else {
    do.call("plot3d", delete.args("plot3d", args,
      c("xlim", "ylim", "zlim", "pch", "main", "xlab", "ylab", "zlab")))
  }
}


##################################
## (7) Other helping functions. ##
##################################
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



##################################
## (8) Model summary functions. ##
##################################
summary.bayesr <- function(object, model = NULL, ...)
{
  call <- attr(object, "call")
  family <- attr(object, "family")
  object <- get.model(object, model)
  rval <- list()
  n <- length(object)
  for(i in 1:n) {
    if(!any(c("param.effects", "effects.hyp", "scale") %in% names(object[[i]]))) {
      rval[[i]] <- summary.bayesr(object[[i]])
      attr(rval[[i]], "hlevel") <- TRUE
    } else {
      for(j in c("param.effects", "effects.hyp", "scale")) {
        if(!is.null(object[[i]][[j]]))
          attr(object[[i]][[j]], "samples") <- NULL
      }
      rval[[i]] <- with(object[[i]],
        c(list("param.effects" = param.effects,
          "effects.hyp" = effects.hyp, "scale" = scale),
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
  class(rval) <- "summary.bayesr"
  rval
}

print.summary.bayesr <- function(x, digits = max(3, getOption("digits") - 3), ...)
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
      print.summary.bayesr(x[[i]], digits = digits, dic_out = FALSE, ...)
      if(i == n & dic_out)
        print_dic_pd(x[[i]][[1]], ok = FALSE)
    } else {
      cat("Formula:\n")
      print(x[[i]]$formula)
      if(length(x[[i]]$param.effects) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x[[i]]$param.effects, digits = digits, na.print = "NA", ...)
      }
      if(length(x[[i]]$effects.hyp) > 0) {
        cat("\nSmooth effects variances:\n")
        printCoefmat(x[[i]]$effects, digits = digits, na.print = "NA", ...)
      }
      if(!is.function(x[[i]]$scale)) {
        if(!is.null(x[[i]]$scale)) {
          cat("\nScale estimate:\n")
          printCoefmat(x[[i]]$scale, digits = digits, na.print = "NA", ...)
        }
      }
      if(i == n & !h0) {
        print_dic_pd(x[[1]])
      } else cat("\n")
    }
  }
}


## Simple "bayesr" print method.
print.bayesr <- function(x, digits = max(3, getOption("digits") - 3), ...) 
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
        cat(nh[j], ": ", sep = ""); print(xs[[i]][[j]]$formula)
      }
      if(i < n) cat("---\n")
      if(i == n & pdic) {
        print_dic_pd(xs[[i]][[1]])
        pdic <- FALSE
      }
    } else {
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


###################################
## (9) More extractor functions. ##
###################################
DIC <- function(object, ...)
{
  UseMethod("DIC")
}

DIC.bayesr <- function(object, ...)
{
  object <- c(object, ...)
  rval <- NULL
  for(i in 1:length(object)) {
    xs <- summary(object[[i]])
    n <- attr(xs, "n")
    if(n < 2)
      xs <- list(xs)
    rval <- rbind(rval, data.frame(
      "DIC" = xs[[n]]$DIC,
      "pd" = xs[[n]]$pd
    ))
  }
  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}


## Extract model formulas.
formula.bayesr <- function(x, model = NULL, ...)
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

print.bayesr.formula <- function(x, ...) {
  if(!inherits(x, "list")) {
    print(x)
  } else {
    nx <- names(x)
    if(is.null(nx))
      nx <- as.character(1:length(x))
    for(i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if(inherits(x[[i]], "list")) {
        for(j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          print(x[[i]][[j]])
        }
      } else print(x[[i]])
      if(i < length(x))
      cat("\n")
    }
  }
  invisible(NULL)
}


## Extract formula terms.
terms.bayesr <- function(x, ...) {
  terms.bayesr.formula(formula(x, ...))
}

terms.bayesr.formula <- function(x)
{
  if(!inherits(x, "list")) {
    tx <- terms(x)
  } else {
    tx <- list()
    for(i in seq_along(x)) {
      if(inherits(x[[i]], "list")) {
        tx[[i]] <- terms.bayesr.formula(x[[i]])
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
fitted.bayesr <- function(object, model = NULL, term = NULL,
  type = c("link", "parameter"), samples = FALSE, FUN = mean,
  nsamps = NULL, ...)
{
  type <- match.arg(type)
  family <- attr(object, "family")

  if(type != "parameter" & !samples)
    object <- get.model(object, model, drop = FALSE)

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
        rval[[j]] <- fitted.bayesr(object[[j]], term = term, ...)
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
      rval[[j]] <- predict.bayesr(object, model = if(one) NULL else j, term = term,
        FUN = function(x) { x }, nsamps = nsamps, ...)
      if(type != "link")
        rval[[j]] <- apply(rval[[j]], 2, make.link2(family$links[if(one) 1 else nrval[j]])$linkinv)
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

samples.bayesr <- function(x, model = NULL, term = NULL, ...)
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
      samps[[j]] <- samples.bayesr(x[[j]], term = term, id = nx[j])
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
  samps <- if(length(samps) < 1) NULL else as.mcmc(do.call("cbind", samps))
  samps <- samps[, unique(colnames(samps)), drop = FALSE]
  samps
}


## Credible intervals of coefficients.
confint.bayesr <- function(object, parm, level = 0.95, model = NULL, ...)
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
coef.bayesr <- function(object, model = NULL, term = NULL, FUN = mean, ...)
{
  object <- get.model(object)
  if(is.null(term))
    term <- all.terms(object, ne = TRUE, id = FALSE)
  samps <- samples(object, model = model, term = term)
  apply(samps, 2, FUN, ...)
}


## Get all terms names used.
all.terms <- function(x, model = NULL, ne = TRUE, ...)
{
  args <- list(...)
  nx <- names(x)
  if(!ne) {
    tx <- terms.bayesr(x, model)
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

  tl
}


## Get the model.frame.
model.frame.bayesr <- function(formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  if(length(nargs) || is.null(attr(formula, "model.frame"))) {
    fcall <- attr(formula, "call")
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(bayesr.model.frame)
    fcall[names(nargs)] <- nargs
    env <- environment(attr(formula, "formula"))
    if(is.null(env)) 
      env <- parent.frame()
    eval(fcall, env)
  }
  else attr(formula, "model.frame")
}


## Scores for model comparison.
score <- function(x, limits = NULL, FUN = function(x) { mean(x, na.rm = TRUE) },
  type = c("mean", "samples"), nsamps = NULL, ...)
{
  stopifnot(inherits(x, "bayesr"))
  family <- attr(x, "family")
  stopifnot(!is.null(family$d))
  type <- match.arg(type)
  y <- model.response2(x)
  n <- length(y)
  maxy <- max(y, na.rm = TRUE)

  if(is.null(family$score.norm)) {
    score.norm <- function(eta) {
      integrand <- function(x) {
        family$d(x, eta)^2
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
    score.norm <- function(eta) {
	    integrand <- function(x) {
        family$d(x, eta)^2
	    }
	    rval <- sum(integrand(seq(0, maxy)))
	    rval
    }

	  score.norm2 <- function(y, eta) {
	    integrand <- function(x) {
         -sum(((x == y) * 1 - family$d(x, eta))^2)
	    }
	    rval <- (integrand(seq(0, maxy)))
	    rval
    }
  }

  scorefun <- function(eta) {
    norm <- rep(0, n)
    for(i in 1:n) {
      norm[i] <- score.norm(eta[i, , drop = FALSE])
    }

    pp <- family$d(y, eta)
    loglik <- log(pp)
    if(is.null(family$score.norm)) {
      quadratic <- 2 * pp - norm
    } else {
      quadratic <- rep(0, n)
      for(i in 1:n)
        quadratic[i] <- score.norm2(y[i], eta[i, , drop = FALSE])
    }
    spherical <- pp / sqrt(norm)

    return(data.frame(
      "logLik" = FUN(loglik),
      "quadratic" = FUN(quadratic),
      "spherical" = FUN(spherical)
    ))
  }

  if(type == "mean") {
    eta <- list()
    for(j in family$names)
      eta[[j]] <- fitted(x, model = j)
    eta <- as.data.frame(eta)
    res <- unlist(scorefun(eta))
  } else {
    nx <- names(x)
    eta <- list()
    for(j in nx) {
      eta[[j]] <- predict.bayesr(x, model = j, FUN = function(x) { x }, nsamps = nsamps)
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


## Extract model residuals.
residuals.bayesr <- function(object, type = c("quantile", "ordinary"),
  FUN = mean, nsamps = NULL)
{
  type <- match.arg(type)
  y <- model.response2(object)
  family <- attr(object, "family")
  if(!is.null(family$type)) family$type <- 1
  if(type == "ordinary") {
    if(is.factor(y)) y <- as.integer(y) - 1
  } else stopifnot(!is.null(family$p))
  res <- fitted(object, type = "link", samples = TRUE, FUN = function(x) { x }, nsamps = nsamps)
  if(!is.list(res)) {
    res <- list(res)
    names(res) <- family$names
  }
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
    tres <- if(type == "quantile") {
      if(family$type == 3) {
        le <- family$p(y - 1, eta2)
        ri <- family$p(y, eta2)
        qnorm(runif(length(y), min = le, max = ri))
      } else family$p(y, eta2)
    } else family$mu(eta2)
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
  if(type == "ordinary") {
    if(is.null(dim(y))) {
      res <- y - res
    } else {
      for(j in 1:ncol(y))
        res[[j]] <- y[, j] - res[[j]]
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
  } else data[, rn]
  y
}


#############################
## (10) Utility functions. ##
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
  require("BayesR")
  if(is.null(file))
    file <- "~/svn/bayesr/pkg/BayesR/R/families.R"
  file <- path.expand(file)
  env <- new.env()
  source(file, local = env)
  fun <- grep(".BayesR", ls(env), fixed = TRUE, value = TRUE)
  fun <- fun[!grepl("print", fun)]
  tab <- NULL
  for(i in seq_along(fun)) {
    fe <- try(eval(parse(text = paste(fun[i], "()", sep = "")), envir = env), silent = TRUE)
    if(inherits(fe, "family.BayesR")) {
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
        "IWLS" = if(!is.null(fe$score) & !is.null(fe$weights) & !is.null(fe$loglik)) "yes" else "no"
      ))
    }
  }
  as.data.frame(tab)
}

#.First.lib <- function(lib, pkg)
#{
#  library.dynam("BayesR", pkg, lib)
#}


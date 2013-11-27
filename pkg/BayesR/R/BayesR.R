################################################
## (1) BayesR main model fitting constructor. ##
################################################
## Could be interesting: http://people.duke.edu/~neelo003/r/
##                       http://www.life.illinois.edu/dietze/Lectures2012/
bayesr <- function(formula, family = gaussian.BayesR, data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  reference = NULL, parse.input = parse.input.bayesr, transform = transformJAGS,
  setup = setupJAGS, sampler = samplerJAGS, results = resultsJAGS,
  cores = NULL, sleep = NULL, combine = TRUE, model = TRUE, grid = 100, ...)
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
      if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
      functions$sampler(ms)
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


##########################################################
## (2) Parsing all input using package mgcv structures. ##
##########################################################
parse.input.bayesr <- function(formula, data, family = gaussian.BayesR,
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
    if(inherits(fe, "formula")) {
      formula2 <- c(formula2, fe)
      eval(parse(text = paste(f, if(is.null(fo[[f]])) "NULL" else fo[[f]], sep = " = ")))
    }
  }

  ## Parse family object.
  family <- bayesr.family(family)

  ## Parse formula
  formula <- bayesr.formula(c(formula, formula2), specials, family)

  ## Create the model frame.
  mf <- bayesr.model.frame(formula, data, family, weights,
    subset, offset, na.action, specials)
  response.name <- attr(mf, "response.name")

  ## For categorical responses, extend formula object.
  ylevels <- NULL
  if(is.factor(mf[[response.name]])) {
    cat <- if(is.null(family$cat)) FALSE else family$cat
    if(cat & nlevels(mf[[response.name]]) > 2) {
      if(is.null(reference)) {
        ty <- table(mf[[response.name]])
        reference <- c(names(ty)[ty == max(ty)])[1]
      }
      reference <- rmf(reference)
      ylevels <- rmf(levels(mf[[response.name]]))
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

  ## Assign all design matrices and the hierarchical level, if any.
  rval <- bayesr.design(formula, mf, contrasts, knots, ...)
  rval <- bayesr.hlevel(rval)

  attr(rval, "family") <- family
  attr(rval, "reference") <- reference
  attr(rval, "ylevels") <- ylevels
  attr(rval, "grid") <- grid
  attr(rval, "model.frame") <- mf

  class(rval) <- c("bayesr.input", "list")

  rval
}

"[.bayesr.input" <- function(x, ...) {
  rval <- NextMethod("[")
  xattr <- attributes(x)
  mostattributes(rval) <- attributes(x)
  rval
}

"[.bayesr" <- function(x, ...) {
  rval <- NextMethod("[")
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

  ## assign main response name
  ok <- all(c("formula", "response", "fake.formula") %in% names(formula[[1]]))
  rn <- attr(mf, "response.name") <- if(!ok) formula[[1]][[1]]$response else formula[[1]]$response

  ## Check response.
  if(!is.null(family$valid.response))
    family$valid.response(mf[[rn]])

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
  formula <- complete_formula(formula_hierarchical(formula))
  formula <- formula_extend(formula, specials, family)

  environment(formula) <- env
  class(formula) <- "bayesr.formula"

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
  o <- order(from, decreasing = TRUE)
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
    hm <- sapply(j[i], max)
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
  class(x) <- c("bayesr", "bayesx")

  return(x)
}


## Function to compute statistics from samples of a model term.
compute_term <- function(x, get.X, get.mu, psamples, vsamples = NULL,
  FUN = NULL, snames, effects.hyp, fitted.values, data,
  grid = 100, rug = TRUE)
{
  require("coda")

  ## Data for rug plotting.
  rugp <- if(length(x$term) < 2 & rug) data[[x$term]] else NULL

  ## Compute new data set for which effects should
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
    effects.hyp <- rbind(effects.hyp, smatfull)
  }

  ## Assign samples.
  attr(smf, "samples") <- as.mcmc(psamples)

  return(list("term" = smf, "effects.hyp" = effects.hyp, "fitted.values" = fitted.values))
}


#####################
## (4) Prediction. ##
#####################
## A prediction method for "bayesr" objects.
## Prediction can also be based on multiple chains.
predict.bayesr <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = mean, trans = NULL, ...)
{
  object <- get.model(object, model)
  if(any(c("effects", "param.effects") %in% names(object)))
    object <- list(object)
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
    newdata <- model.frame(object)
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
          m.specials[[i]] <- list("X" = PredictMat(specs, newdata), ## FIXME: also allow basis()?
            "get.mu" = specs$get.mu, "samples" = tmp)
        } else {
          m.samples[[i]] <- rbind(m.samples[[i]], tmp)
          if(inherits(object[[j]]$effects[[i]], "linear.bayesx")) {
            if(specs$is.factor & !is.character(newdata)) {
              hi <- if(!is.null(object[[j]]$param.effects)) {
                any(grepl("(Intercept)", rownames(object[[j]]$param.effects), fixed = TRUE))
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
      }
    }

    if(intercept & j < 2) {
      sami <- attr(object[[j]]$param.effects, "samples")
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
    if(!is.null(trans))
      rval <- apply(rval, 1, trans)
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
  acons <- TRUE
  if(!is.null(object$xt$center))
    acons <- object$xt$center
  object <- smoothCon(object, data, knots, absorb.cons = acons)[[1]]
  object$by.done <- TRUE
  object$get.mu <- function(X, g) {
    k <- ncol(X)
    w <- c(1, g[(k + 1):(2 * k - 1)])
    g <- g[1:k]
    R <- diag(1 / drop(X %*% w)) %*% X %*% diag(w)
    R %*% g
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
  args <- list(...)

  if(is.null(args$do_par) & spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  x <- get.model(x, model)

  ## What should be plotted?
  which.match <- c("effects", "samples", "hist-resid", "qq-resid",
    "scatter-resid", "scale-resid", "scale-samples", "max-acf", "param-samples")
  if(!is.character(which)) {
    if(any(which > 8L))
      which <- which[which <= 8L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

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
    if(!ask) par(mfrow = n2mfrow(kn[if(which == "effects") 1 else 2])) else par(ask = ask)
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
  if(which == "scale-samples") {
    for(i in 1:n) {
      args$x <- attr(x[[i]]$scale, "samples")
      do.call("plot", args)
    }
  }
  if(which == "param-samples") {
    for(i in 1:n) {
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
delete.args <- function(fun = NULL, args = NULL, not = NULL)
{
  nf <- names(formals(fun))
  na <- names(args)
  for(elmt in na)
    if(!elmt %in% nf) {
      if(!is.null(not)) {
        if(!elmt %in% not)
          args[elmt] <- NULL
      } else args[elmt] <- NULL
    }
  args
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


formula.bayesr <- function(x, ...) { attr(x, "formula") }


## Model extractor function.
get.model <- function(x, model = NULL)
{
  if(length(model) > 1) {
    for(j in model)
      x <- get.model(x, j)
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
        x <- x[[model]]
      }
    } else x <- list(x)
    class(x) <- cx
    attr(x, "family") <- family
  }

  return(x)
}


#############################
## (10) Utility functions. ##
#############################
scale2 <- function(x, lower = 0, upper = 1)
{
  x <- if(length(unique(x)) > 1) {
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (upper - lower) + lower
  } else x
  x
}


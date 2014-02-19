#########################################
## (1) Special JAGS transform function ##
#########################################
transformJAGS <- function(x)
{
  family <- attr(x, "family")
  cat <- if(is.null(family$cat)) FALSE else family$cat

  if(cat) {
    reference <- attr(x, "reference")

    xr <- list(
      "formula" = as.formula(paste(reference, "~ 1", sep = "")),
      "fake.formula" = as.formula(paste(reference, "~ 1", sep = "")),
      "cat.formula" = as.formula(paste(reference, "~ 1", sep = "")),
      "intercept" = TRUE,
      "X" = matrix(1, nrow = nrow(attr(x, "model.frame")), ncol = 1)
    )

    nx <- names(x)
    x[[length(x) + 1]] <- xr
    names(x) <- c(nx, reference)
    attr(x, "ylevels") <- c(attr(x, "ylevels"), reference)

    tJAGS <- function(obj) {
      if(inherits(obj, "bayesr.input") & !any(c("smooth", "response") %in% names(obj))) {
        no <- names(obj)
        no <- no[no != "call"]
        if(is.null(no)) no <- 1:length(obj)
        if(length(unique(no)) < length(obj)) no <- 1:length(obj)
        for(j in no)
          obj[[j]] <- tJAGS(obj[[j]])
      }
      obj
    }

    x <- tJAGS(x)
  }

  x <- randomize(x)
  
  x
}


######################################################################
## (2) General JAGS setup functions, model code, initials and data. ##
######################################################################
## Setup the model structure needed for the sampler.
## These functions create the model code, initials and data
## to run a JAGS sampler.
## Examples: http://sourceforge.net/projects/mcmc-jags/files/
##           http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/
## Default linear predictor and model setup functions.
JAGSeta <- function(x, id = NULL, zero = FALSE, ...) {
  setup <- list()
  if(is.null(zero)) zero <- FALSE
  setup$inits <- list()

  ## Parametric part.
  if(k <- ncol(x$X)) {
    id2 <- if(is.null(id)) 1 else id
    setup$data[[paste("X", id2, sep = "")]] <- x$X
    for(j in 1:k) {
      setup$param <- c(setup$param, paste(if(k < 2) paste("beta", id2, sep = "") else paste("beta", id2, "[", j, "]",
        sep = ""), "*X", id2, "[i, ", j, "]", sep = ""))
    }
    setup$param <- paste(setup$param, collapse = " + ")
    setup$param <- paste("    param", id2, "[i] <- ", setup$param, sep = "")
    setup$loops <- k
    setup$priors.coef <- if(k > 1) {
      paste("    beta", id2, if(zero) "[j] <- 0.0" else "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  beta", id2, if(zero) " <- 0.0" else " ~ dnorm(0, 1.0E-6)", sep = "")
    if(!zero)
      setup$inits[[paste("beta", id2, sep = "")]] <- runif(k)
    setup$psave <- c(setup$psave, paste("beta", id2, sep = ""))
    setup$eta <- paste("param", id2, "[i]", sep = "")
  }

  ## Smooth part.
  if(m <- length(x$smooth)) {
    for(i in 1:m) {
      setup <- if(!is.null(x$smooth[[i]]$special)) {
        buildJAGS.smooth.special(x$smooth[[i]], setup, paste(i, id, sep = ""), zero)
      } else {
        buildJAGS.smooth(x$smooth[[i]], setup, paste(i, id, sep = ""), zero)
      }
    }
  }

  ## Final touch ups.
  if(!is.null(x$response.vec))
    setup$data$response <- x$response.vec

  setup
}

## Get link functions.
JAGSlinks <- function(x)
{
  switch(x,
    "identity" = "eta",
    "log" = "exp(eta)",
    "exp" = "log(eta)",
    "inverse" = "1 / (eta^2)",
    "logit" = "1 / (1 + exp(-(eta)))",
    "probit" = "phi(eta)",
    "cloglog" = "log(-log(1 - eta))",
    "pow" = "pow(eta, -2)"
  )
}

## Construct the final model code.
JAGSmodel <- function(x, family, cat = FALSE, is.stan = FALSE, ...) {
  if(is.function(family))
    family <- family()
  k <- if(all(c("inits", "data", "psave") %in% names(x))) {
    x <- list(x)
    1
  } else length(x)
  if(k > length(family$names) & !cat) {
    stop(paste("more parameters specified than existing in family ",
      family$family, ".BayesR()!", sep = ""), call. = FALSE)
  }

  model <- "model {"
  for(j in 1:k) {
    model <- c(model, x[[j]]$start)
  }
  pn <- family$names
  if(!is.null(family$jagstan$reparam))
    pn[repi <- match(names(family$jagstan$reparam), pn)] <- paste("rp", 1:length(family$jagstan$reparam), sep = "")
  if(is.null(pn)) pn <- paste("theta", 1:k, sep = "")
  if(length(pn) < 2 & length(pn) != k)
    pn <- paste(pn, 1:k, sep = "")

  pn[1:k] <- paste(pn[1:k], "[i]", sep = "")
  on <- if(cat) family$names else NULL
  links <- family[[grep("links", names(family), fixed = TRUE, value = TRUE)]]
  links <- rep(sapply(links, JAGSlinks), length.out = k)
  model <- c(model,  "  for(i in 1:n) {",
    paste("    response[i] ~ ", family$jagstan$dist, "(",
      paste(if(is.null(on)) pn else paste(on, ## if(cat) "n" else NULL,
      "[i, 1:", k, "]", sep = ""),
      collapse = ", "), ")", sep = ""))
#  if(cat) {
#    npn <- if(is.null(on)) pn else on
#    model <- c(model, paste("    ", npn, "n[i, ", 1:k, "] <- ",
#      npn, "[i, ", 1:k, "] / sum(", npn, "[i, 1:", k, "])", sep = ""))
#  }

  if(!is.null(family$jagstan$reparam)) {
    reparam <- NULL
    for(j in seq_along(family$jagstan$reparam))
      reparam <- c(reparam, paste("    rp", j, "[i] <- ", family$jagstan$reparam[j], sep = ""))
    for(j in family$names)
      reparam <- gsub(j, paste(j, "[i]", sep = ""), reparam)
    model <- c(model, reparam)
    pn[repi] <- paste(family$names[repi], "[i]", sep = "")
  }
  if(!is.null(family$jagstan$addparam)) {
    for(j in family$jagstan$addparam)
      model <- c(model, paste("   ", j))
  }
  if(!is.null(family$jagstan$addvalues)) {
    for(j in names(family$jagstan$addvalues))
      model <- gsub(j, family$jagstan$addvalues[[j]], model)
  }

  for(j in 1:k) {
    model <- c(model, paste("    ", if(is.null(on)) pn[j] else paste(on, "[i, ", j, "]", sep = ""),
      " <- ", gsub("eta", x[[j]]$eta, links[[j]]), sep = ""))
  }
  for(j in 1:k)
    model <- c(model, x[[j]]$adds)
  for(j in 1:k)
    model <- c(model, x[[j]]$param, x[[j]]$smooth)
  model <- c(model, "  }")

  for(i in 1:k)
    model <- c(model, x[[i]]$close1)

  for(i in 1:k) {
    lp <- list()
    if(!is.null(x[[i]]$loops)) {
      for(j in 1:length(x[[i]]$loops))
        lp[[paste(x[[i]]$loops[j])]] <- c(lp[[paste(x[[i]]$loops[j])]], x[[i]]$priors.coef[j])
    }
    if(length(lp)) {
      for(j in names(lp)) {
        if(j != 1)
          tmp <- c(paste("  for(j in 1:", j, ") {", sep = ""), lp[[j]], "  }")
        else
          tmp <- lp[[j]]
        model <- c(model, tmp)
      }
    }
    model <- c(model, x[[i]]$priors.scale, x[[i]]$close2, x[[i]]$close3)
  }
  model <- c(model, "}")

  model
}

## Create final model setup.
setupJAGS <- function(x)
{
  is.stan <- if(is.null(attr(x, "is.stan"))) FALSE else TRUE
  family <- attr(x, "family")
  reference <- attr(x, "reference")
  ylevels <- attr(x, "ylevels")
  if(is.function(family))
    family <- family()
  ncat <- NULL
  if(!all(c("formula", "fake.formula", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx))
      nx <- 1:length(x)
    rval <- list()
    fn <- family$names
    cat <- if(!is.null(family$cat)) family$cat else FALSE
    if(cat)
      fn <- gsub(attr(attr(x, "model.frame"), "response.name"), "", names(x))
    if(length(fn) < length(x))
      fn <- paste(fn, 1:length(nx), sep = "")
    for(i in seq_along(nx)) {
      rval[[nx[i]]] <- family$jagstan$eta(x[[i]], fn[i],
        zero = if(!is.null(reference)) ylevels[i] == reference else NULL)
    }
  } else {
    rval <- family$jagstan$eta(x)
  }
  
  ## Create model code.
  model <- family$jagstan$model(rval, family, cat = !is.null(reference), is.stan)

  ## Collect data.
  if(all(c("inits", "data", "psave") %in% names(rval)))
    rval <- list(rval)
  data <- inits <- psave <- NULL
  for(j in seq_along(rval)) {
    data <- c(data, rval[[j]]$data)
    inits <- c(inits, rval[[j]]$inits)
    psave <- c(psave, rval[[j]]$psave)
  }
  data <- data[unique(names(data))]
  inits <- inits[unique(names(inits))]
  psave <- unique(psave)
  data$n <- nrow(attr(x, "model.frame"))
  psave <- c(psave, attr(model, "psave"))

  if(is.factor(data$response)) {
    nl <- nlevels(data$response)
    data$response <- as.integer(data$response)
    if(nl < 3)
      data$response <- data$response - 1
  }

  rval <- list("model" = model, "data" = data,
    "inits" = inits, "psave" = psave)
  
  return(rval)
}


## Build the JAGS model code for a smooth term. 
buildJAGS.smooth <- function(smooth, setup, i, zero) {
  fall <- NULL
  kr <- if(is.null(smooth$rand$Xr)) 0 else ncol(smooth$rand$Xr)
  kx <- if(is.null(smooth$Xf)) 0 else ncol(smooth$Xf)
  hcauchy <- if(is.null(smooth$xt$hcauchy)) FALSE else smooth$xt$hcauchy
  if(kx > 0) {
    fall <- c(fall, paste("b", i, if(kx > 1) paste("[", 1:kx, "]", sep = ""),
      "*Xf", i, "[i, ", 1:kx, "]", sep = ""))
    setup$data[[paste("Xf", i, sep = "")]] <- smooth$Xf
    tmp <- if(kx > 1) {
        paste("    b", i, if(zero) "[j] <- 0.0" else "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  b", i, if(zero) " <- 0.0" else " ~ dnorm(0, 1.0E-6)", sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kx)
    if(!zero)
      setup$inits[[paste("b", i, sep = "")]] <- runif(kx, 0.1, 0.2)
    setup$psave <- c(setup$psave, paste("b", i, sep = ""))
  }
  if(kr > 0) {
    fall <- c(fall, paste("g", i, if(kr > 1) paste("[", 1:kr, "]", sep = ""), "*Xr",
      i, "[i, ", 1:kr, "]", sep = ""))
    setup$data[[paste("Xr", i, sep = "")]] <- smooth$rand$Xr
    taug <- paste("taug", if(is.null(smooth$id)) i else smooth$id, sep = "")
    tmp <- if(kr > 1) {
      paste("    g", i, if(zero) "[j] <- 0.0" else paste("[j] ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    } else paste("g", i, if(zero) " <- 0.0" else paste(" ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kr)
    if(is.null(setup$priors.scale) || !any(grepl(taug, setup$priors.scale))) {
      if(!hcauchy) {
        setup$priors.scale <- c(setup$priors.scale, paste("  ", taug,
          if(zero) " <- 0.0" else " ~ dgamma(1.0E-6, 1.0E-6)", sep = ""))
      } else {
        setup$priors.scale <- c(setup$priors.scale, paste("  ", taug,
          if(zero) " <- 0.0" else " <- abs(", taug, 0, ")", sep = ""),
          paste("  ", taug, 0, "~ dt(0, 10, 1)", sep = ""))
      }
      if(!zero)
        	setup$inits[[taug]] <- runif(1, 0.1, 0.2)
      setup$psave <- c(setup$psave, taug)
    }
    if(!zero)
      setup$inits[[paste("g", i, sep = "")]] <- runif(kr, 0.1, 0.2)
    setup$psave <- c(setup$psave, paste("g", i, sep = ""))
  }

  setup$smooth <- c(setup$smooth, paste("    sm", i, "[i] <- ",
    paste(fall, collapse = " + ", sep = ""), sep = ""))
  setup$eta <- paste(setup$eta, paste("sm", i, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


## For special terms, e.g. growth curves, this function
## builds the model code.
buildJAGS.smooth.special <- function(smooth, setup, i, zero)
{
  UseMethod("buildJAGS.smooth.special")
}


## Default special model term builder.
buildJAGS.smooth.special.default <- function(smooth, setup, i, zero)
{
  buildJAGS.smooth(smooth, setup, i, zero)
}


## JAGS random scaling model term constructor.
buildJAGS.smooth.special.rsc.smooth <- function(smooth, setup, i, zero)
{
  smooth$special <- FALSE
  setup <- buildJAGS.smooth(smooth, setup, i, zero)

  if(!is.null(smooth$by.formula)) {
    st <- setup$smooth[si <- grep(paste("sm", i, "[i] <-", sep = ""), setup$smooth, fixed = TRUE)]
    st <- strsplit(st, " <- ", fixed = TRUE)[[1]]
    n <- length(smooth$rs.by)
    rs <- paste("rs", i, 1:n, "[rsid", i, 1:n, "[i]]", sep = "", collapse = " + ")
    st[2] <- paste("(", if(smooth$one) "1 + ", rs, ")*(", st[2], ")", sep = "")
    setup$smooth[si] <- paste(st[1], "<-", st[2])
    for(j in 1:n) {
      tmp <- paste("    rs", i, j, "[j] ~ dnorm(0, taugrs", i, j, ")", sep = "")
      setup$priors.coef <- c(setup$priors.coef, tmp)
      setup$loops <- c(setup$loops, nlevels(smooth$rs.by[[j]]))
      setup$priors.scale <- c(setup$priors.scale,
        paste("  ", "taugrs", i, j, " ~ dgamma(1.0E-6, 1.0E-6)", sep = ""))
      setup$psave <- c(setup$psave, paste("rs", i, j, sep = ""),
        paste("taugrs", i, j, sep = ""))
      setup$data[[paste("rsid", i, j, sep = "")]] <- as.integer(smooth$rs.by[[j]])
      setup$inits[[paste("rs", i, j, sep = "")]] <- rnorm(nlevels(smooth$rs.by[[j]]))
    }
  }

  setup
}


## Special code builder for growth curve terms.
buildJAGS.smooth.special.gc.smooth <- function(smooth, setup, i, zero)
{
  center <- if(is.null(smooth$xt$center)) TRUE else smooth$xt$center
  pn <- paste("g", i, sep = "")

  setup$data[[paste("X", pn, sep = "")]] <- as.numeric(smooth$X)
  setup$inits[[paste(pn, "0", sep = "")]] <- runif(3, 0.1, 0.2)
  setup$psave <- c(setup$psave, pn)

  setup$close2 <- c(setup$close2,
    "  for(j in 1:3) {",
    paste("    ", pn, "[j] <- exp(", pn, "0[j])", sep = ""),
    paste("    ", pn, "0[j] ~ dnorm(0, 1.0E-6)", sep = "")
  )

  ## logistic
  ## fall <- paste("g", i, "[1] / (1 + g", i, "[2]*exp(g", i, "[3]*Xgc", i, "[i]))", sep = "")

  ## Gompertz growth function.
  if(is.null(smooth$by.levels)) {
    fall <- paste(pn, "[1]*exp(-", pn, "[2]*exp(-", pn, "[3]*X", pn, "[i]))", sep = "")
  } else {
    setup$data[[paste(pn, "id", sep = "")]] <- as.integer(smooth$fid)
    fall <- paste("(", pn, "r[", pn , "id[i], 1] + ", pn, "[1])*exp(-(", pn, "r[", pn,
      "id[i], 2] + ", pn, "[2])*exp(-(", pn, "r[", pn, "id[i], 3] + ", pn, "[3])*X", pn, "[i]))", sep = "")
    setup$close2 <- c(setup$close2,
    paste("    for(k in 1:", length(smooth$by.levels), ") {", sep = ""),
      paste("      ", pn, "r[k, j] <- exp(", pn, "r0[k, j])", sep = ""),
      paste("      ", pn, "r0[k, j] ~ dnorm(0, taug", i, "[j])", sep = ""),
      "    }",
      paste("    taug", i, "[j] ~ dgamma(1.0E-3, 1.0E-3)", sep = "")
    )
    setup$psave <- c(setup$psave, paste(pn, "r", sep = ""), paste("taug", i, sep = ""))
    setup$inits[[paste(pn, "r0", sep = "")]] <- matrix(runif(length(smooth$by.levels) * 3, 0.1, 0.2), ncol = 3)
    setup$inits[[paste("taug", i, sep = "")]] <- runif(3, 0.1, 0.2)
  }

  setup$close2 <- c(setup$close2, "  }")

  prefix <- if(center) "0" else NULL
  if(!center) {
    setup$smooth <- c(setup$smooth, paste("    sm", i, "[i] <- ",
      paste(fall, collapse = " + ", sep = ""), sep = ""))
  } else {
    setup$close1 <- c(setup$close1, paste("  sm", i,  1, " <- sm", i, 0, " - mean(sm", i, 0, ")", sep = ""))
    setup$close1 <- c(setup$close1,
      paste("  for(i in 1:n) {", sep = ""),
      paste("    sm", i, 0, "[i] <- ",
        paste(fall, collapse = " + ", sep = ""), sep = ""), "  }")
  }
  setup$eta <- paste(setup$eta, paste("sm", i, if(center) 1 else NULL, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


## Special code builder for rational splines.
buildJAGS.smooth.special.rs.smooth <- function(smooth, setup, i, zero)
{
  fall <- fall0 <- NULL
  kr <- if(is.null(smooth$rand$Xr)) 0 else ncol(smooth$rand$Xr)
  kx <- if(is.null(smooth$Xf)) 0 else ncol(smooth$Xf)
  kx <- if(kr < 1 & kx < 1) ncol(smooth$X) else kx

  if(kx > 0) {
    fall0 <- c(fall0, paste("Xf", i, "[i, ", 1:kx, "]", sep = ""))
    fall <- c(fall, paste("b", i, if(kx > 1) paste("[", 1:kx, "]", sep = ""),
      "*Xf", i, "[i, ", 1:kx, "]", sep = ""))
    setup$data[[paste("Xf", i, sep = "")]] <- if(smooth$fixed) smooth$X else smooth$Xf
    tmp <- if(kx > 1) {
        paste("    b", i, if(zero) "[j] <- 0.0" else "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  b", i, if(zero) " <- 0.0" else " ~ dnorm(0, 1.0E-6)", sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kx)
    if(!zero)
      setup$inits[[paste("b", i, sep = "")]] <- runif(kx, 0.1, 0.2)
    setup$psave <- c(setup$psave, paste("b", i, sep = ""))
  }
  if(kr > 0) {
    fall0 <- c(fall0, paste("Xr", i, "[i, ", 1:kr, "]", sep = ""))
    fall <- c(fall, paste("g", i, if(kr > 1) paste("[", 1:kr, "]", sep = ""), "*Xr",
      i, "[i, ", 1:kr, "]", sep = ""))
    setup$data[[paste("Xr", i, sep = "")]] <- smooth$rand$Xr
    taug <- paste("taug", if(is.null(smooth$id)) i else smooth$id, sep = "")
    tmp <- if(kr > 1) {
      paste("    g", i, if(zero) "[j] <- 0.0" else paste("[j] ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    } else paste("g", i, if(zero) " <- 0.0" else paste(" ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kr)
    if(is.null(setup$priors.scale) || !any(grepl(taug, setup$priors.scale))) {
      setup$priors.scale <- c(setup$priors.scale, paste("  ", taug,
        if(zero) " <- 0.0" else " ~ dgamma(1.0E-6, 1.0E-6)", sep = ""))
      if(!zero)
        	setup$inits[[taug]] <- runif(1, 0.1, 0.2)
      setup$psave <- c(setup$psave, taug)
    }
  }

#  tmp <- if((kw <- length(fall)) > 1) {
#    paste("    w", i, if(zero) "[j] <- 0.0" else "[j] ~ dgamma(1.0E-6, 1.0E-6)", sep = "")
#  } else paste("  w", i, if(zero) " <- 0.0" else " ~ dgamma(1.0E-6, 1.0E-6)", sep = "")

  tmp <- if((kw <- length(fall)) > 1) {
    paste("    w", i, if(zero) "[j] <- 0.0" else "[j] ~ dunif(1.0E-6, 1)", sep = "")
  } else paste("  w", i, if(zero) " <- 0.0" else " ~ dunif(1.0E-6, 1)", sep = "")

  ## setup$adds <- paste("  w2", i, " <- 1 / sum(w", i, ")", sep = "")
  setup$priors.coef <- c(setup$priors.coef, tmp)
  setup$loops <- c(setup$loops, kw)
  if(!zero)
    setup$inits[[paste("w", i, sep = "")]] <- runif(kw, 0.1, 0.2)
  setup$psave <- c(setup$psave, paste("w", i, sep = ""))

  fall0 <- paste(fall0, c(1, paste("w", i, "[", 1:(length(fall0) - 1), "]", sep = "")), sep = "*")
  ## fall <- paste(fall, c(1, paste("w", i, "[", 1:(length(fall) - 1), "]", sep = "")), sep = "*")

  center <- if(is.null(smooth$xt$center)) TRUE else smooth$xt$center

  if(!center) {
    fall0 <- paste("    sm0", i, "[i] <- 1 / exp(", paste(fall0, collapse = " + "), ")", sep = "")
    fall <- paste("    sm", i, "[i] <- sm0", i, "[i] * (", paste(fall, collapse = " + "), ")", sep = "")
    setup$smooth <- c(setup$smooth, fall, fall0)
  } else {
    fall0 <- paste("    sm0", i, "[i] <- 1 / exp(", paste(fall0, collapse = " + "), ")", sep = "")
    fall <- paste("    sm", i, 0, "[i] <- sm0", i, "[i] * (", paste(fall, collapse = " + "), ")", sep = "")
    setup$start <- c(setup$start,
      paste("  sm", i, " <- sm", i, 0, " - mean(sm", i, 0, ")", sep = ""))
    setup$close <- c(setup$close,
      paste("  for(i in 1:n) {", sep = ""), fall0, fall, "  }")
  }

  setup$eta <- paste(setup$eta, paste("sm", i, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


########################################
## (3) Interface to the JAGS sampler. ##
########################################
samplerJAGS <- function(x, tdir = NULL,
  n.chains = 1, n.adapt = 100,
  n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = TRUE, set.inits = FALSE, ...)
{
  require("rjags")

  ## Temporary directory handling.
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  } else tdir <- path.expand(tdir)
  if(!file.exists(tdir))
    dir.create(tdir)

  ## Write the model code.
  writeLines(paste(x$model, collapse = "\n"), mfile <- file.path(tdir, "model.txt"))

  ## Set the seed of the random number generator.
  if(is.null(seed))
    seed <- floor(runif(n.chains) * .Machine$integer.max)
  inits <- rep(list(x$inits), n.chains)
  for(j in seq_along(inits)) {
    inits[[j]][[".RNG.name"]] <- "base::Super-Duper"
    inits[[j]][[".RNG.seed"]] <- seed[j]
  }

  ## Sampling.
  load.module("dic"); load.module("glm")
  
  if(verbose) writeLines(x$model)
  
  if(set.inits) {
    jmodel <- jags.model(mfile, data = x$data, inits = inits,
      n.chains = n.chains, n.adapt = n.adapt, ...)
  } else {
    jmodel <- jags.model(mfile, data = x$data,
      n.chains = n.chains, n.adapt = n.adapt, ...)
  }
  jsamples <- coda.samples(jmodel, variable.names = c(x$psave, "deviance"),
    n.iter = n.iter, thin = thin, ...)

  ## Remove burnin.
  if(is.null(burnin))
    burnin <- floor(n.iter * 0.2)
  jsamples <- window(jsamples, start = burnin)

  jsamples
}


###################################################################
## (4) Functions to compute summary statistics, plots, etc. from ##
##     samples returned the JAGS sampler.                        ##
###################################################################
## Function to extract all results obtained by running the JAGS
## sampler. The function uses BayesR structures to represent fitted
## model terms etc.
resultsJAGS <- function(x, samples)
{
  family <- attr(x, "family")
  grid <- attr(x, "grid")
  if(is.null(grid)) grid <- 100
  if(is.function(family))
    family <- family()

  createJAGSresults <- function(obj, samples, id = NULL)
  {
    if(inherits(samples[[1]], "mcmc.list"))
      samples <- do.call("c", samples)
    chains <- length(samples)
    rval <- vector(mode = "list", length = chains)
    snames <- colnames(samples[[1]])
    for(j in 1:chains) {
      if(any(grepl("deviance", snames))) {
        DIC <- as.numeric(samples[[j]][, grepl("deviance", snames)])
        pd <- var(DIC) / 2
        DIC <- mean(DIC)
      } else {
        DIC <- pd <- NA
      }

      ## Compute model term effects.
      param.effects <- effects <- effects.hyp <- NULL
      fitted.values <- 0

      ## Parametric effects.
      if(k <- ncol(obj$X)) {
        samps <- as.matrix(samples[[j]][, grepl(paste("beta", id, sep = ""), snames)], ncol = k)
        nx <- colnames(obj$X)
        qu <- t(apply(samps, 2, quantile, probs = c(0.025, 0.5, 0.975)))
        sd <- drop(apply(samps, 2, sd))
        me <- drop(apply(samps, 2, mean))
        param.effects <- cbind(me, sd, qu)
        rownames(param.effects) <- nx
        colnames(param.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
        fitted.values <- as.vector(fitted.values + obj$X %*% param.effects[, 1])
        attr(param.effects, "samples") <- as.mcmc(samps)
        colnames(attr(param.effects, "samples")) <- nx
      }

      ## Smooth terms.
      if(length(obj$smooth)) {
        if(!is.list(effects))
          effects <- list()
        for(i in 1:length(obj$smooth)) {
          if(!is.null(obj$smooth[[i]]$special)) {
            fst <- resultsJAGS.special(obj$smooth[[i]], samples[[j]], attr(x, "model.frame"), i, id = id)
            if(is.null(attr(fst, "by"))) {
              effects[[obj$smooth[[i]]$label]] <- fst$term
              effects.hyp <- rbind(effects.hyp, fst$effects.hyp)
              fitted.values <- fitted.values + fst$fitted.values
            } else {
              tjl <- names(fst)
              for(l in tjl) {
                effects[[l]] <- fst[[l]]$term
                effects.hyp <- rbind(effects.hyp, fst[[l]]$effects.hyp)
                fitted.values <- fitted.values + fst[[l]]$fitted.values
              }
            }
            rm(fst)
          } else {
            ## Get coefficient samples of smooth term.
            xsamples <- rsamples <- NULL
            kr <- if(is.null(obj$smooth[[i]]$rand$Xr)) 0 else ncol(obj$smooth[[i]]$rand$Xr)
            kx <- if(is.null(obj$smooth[[i]]$Xf)) 0 else ncol(obj$smooth[[i]]$Xf)
            kw <- 0
            if(kx) {
              pn <- grep(paste("b", i, id, sep = ""), snames, value = TRUE, fixed = TRUE)
              pn <- pn[!grepl(paste("tau", id, sep = ""), pn)]
              xsamples <- matrix(samples[[j]][, snames %in% pn], ncol = kx)
            }
            if(kr) {
              pn <- grep(paste("g", i, id, sep = ""), snames, value = TRUE, fixed = TRUE)
              pn <- pn[!grepl(paste("taug", i, id, sep = ""), pn)]
              rsamples <- as.matrix(samples[[j]][, snames %in% pn], ncol = kr)
            }
            psamples <- cbind("ra" = rsamples, "fx" = xsamples)
  
            ## Retransform parameter samples.
            if(kr) {
              re_trans <- function(g) {
                g <- obj$smooth[[i]]$trans.D * g
                if(!is.null(obj$smooth[[i]]$trans.U))
                  g <- obj$smooth[[i]]$trans.U %*% g
                g
              }
              psamples <- t(apply(psamples, 1, re_trans))
            }

            ## Possible variance parameter samples.
            vsamples <- NULL
            taug <- paste("taug", if(is.null(obj$smooth[[i]]$id)) i else obj$smooth[[i]]$id, id, sep = "")
            if(taug %in% snames) {
              vsamples <- as.numeric(samples[[j]][, snames %in% taug])
            }

            get.mu <- function(X, g) {
              X %*% as.numeric(g)
            }

            ## Prediction matrix.
            get.X <- function(x) {
              X <- PredictMat(obj$smooth[[i]], x)
              X
            }

            ## Compute final smooth term object.
            tn <- c(obj$smooth[[i]]$term, if(obj$smooth[[i]]$by != "NA") obj$smooth[[i]]$by else NULL)

            if(length(effects)) {
              if(obj$smooth[[i]]$label %in% names(effects)) {
                ct <- gsub(".smooth.spec", "", class(obj$smooth[[i]]))[1]
                if(ct == "random.effect") ct <- "re"
                obj$smooth[[i]]$label <- paste(obj$smooth[[i]]$label, ct, sep = ":")
              }
            }

            fst <- compute_term(obj$smooth[[i]], get.X = get.X, get.mu = get.mu,
              psamples = psamples, vsamples = vsamples, FUN = NULL, snames = snames,
              effects.hyp = effects.hyp, fitted.values = fitted.values,
              data = attr(x, "model.frame")[, tn, drop = FALSE], grid = grid)

            attr(fst$term, "specs")$get.mu <- get.mu 

            ## Add term to effects list.
            effects[[obj$smooth[[i]]$label]] <- fst$term
            effects.hyp <- fst$effects.hyp

            fitted.values <- fst$fitted.values
            rm(fst)
          }
        }
      }

      ## Scale parameters.
      scale.m <- scale.samps.m <- NULL
      sn <- if(is.null(id)) {
        family <- if(is.function(family)) family() else family
        family$names
      } else id
      for(snj in sn) {
        if(snj %in% snames) {
          samps <- 1 / as.numeric(samples[[j]][, snames %in% snj])
          qu <- drop(quantile(samps, probs = c(0.025, 0.5, 0.975)))
          sd <- sd(drop(samps))
          me <- mean(samps)
          scale <- matrix(c(me, sd, qu), nrow = 1)
          colnames(scale) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
          rownames(scale) <- snj
          samps <- matrix(samps, ncol = 1)
          colnames(samps) <- snj
          scale.samps.m <- cbind(scale.samps.m, samps)
          scale.m <- rbind(scale.m, scale)
          attr(scale.m, "samples") <- as.mcmc(scale.samps.m)
        }
      }

      ## Compute partial residuals.
      if(!is.null(effects)) {
        if(length(obj$response)) {
          if(obj$response %in% names(attr(x, "model.frame"))) {
            effects <- partial.residuals(effects, attr(x, "model.frame")[[obj$response]],
              fitted.values, family)
          }
        }
      }

      ## Stuff everything together.
      rval[[j]] <- list(
        "model" = list("DIC" = DIC, "pd" = pd,
          "N" = nrow(attr(x, "model.frame")), "formula" = obj$formula),
        "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "scale" = scale.m, "fitted.values" = fitted.values
      )
      
#      ## Clean.
#      rval[[j]] <- delete.NULLs(rval[[j]])

      class(rval[[j]]) <- "bayesr"
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    if(length(rval) < 2) {
      rval <- rval[[1]]
    }
    class(rval) <- "bayesr"
    return(rval)
  }

  if(inherits(x, "bayesr.input") & !all(c("formula", "fake.formula", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    rval <- list()
    fn <- family$names
    cat <- if(!is.null(family$cat)) family$cat else FALSE
    if(cat)
      fn <- gsub(attr(attr(x, "model.frame"), "response.name"), "", names(x))
    if(length(fn) != length(nx))
      fn <- paste(fn, 1:length(nx), sep = "")
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- createJAGSresults(x[[nx[j]]], samples, id = fn[j])
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
    if(cat) {
      reference <- attr(x, "reference")
      rval <- rval[!grepl(reference, fn)]
    }
    attr(rval, "family") <- family
    class(rval) <- "bayesr"
    return(rval)
  } else {
    return(createJAGSresults(x, samples))
  }
}


## Result extractor function for special terms.
resultsJAGS.special <- function(x, samples, data, i, ...) 
{
  UseMethod("resultsJAGS.special")
}


## Default special term results method.
resultsJAGS.special.default <- function(x, samples, data, i, ...)
{
  snames <- colnames(samples)

  ## Get coefficient samples of smooth term.
  xsamples <- rsamples <- NULL
  kr <- if(is.null(x$rand$Xr)) 0 else ncol(x$rand$Xr)
  kx <- if(x$fixed) {
    ncol(x$X)
  } else {
    if(is.null(x$Xf)) 0 else ncol(x$Xf)
  }
  if(kx) {
    pn <- grep(paste("b", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]
    xsamples <- matrix(samples[, snames %in% pn], ncol = kx)
  }
  if(kr) {
    pn <- grep(paste("g", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]
    rsamples <- as.matrix(samples[, snames %in% pn], ncol = kr)
  }
  psamples <- cbind("ra" = rsamples, "fx" = xsamples)
  
  ## Retransform parameter samples.
  if(kr & !x$fixed) {
    re_trans <- function(g) {
      g <- x$trans.D * g
      if(!is.null(x$trans.U))
        g <- x$trans.U %*% g
      g
    }
    psamples <- t(apply(psamples, 1, re_trans))
  }

  ## Prediction matrix.
  get.X <- function(data) {
    X <- PredictMat(x, data)
    X
  }

  ## Possible variance parameter samples.
  vsamples <- NULL
  taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
  if(taug %in% snames) {
    vsamples <- as.numeric(samples[, snames %in% taug])
  }

  if(is.null(x$get.mu)) {
    get.mu <- function(X, g) {
      X %*% as.numeric(g)
    }
  } else {
    get.mu <- x$get.mu
  }

  ## Get weight parameters for rational splines.
  if(inherits(x, "rs.smooth")) {
    pn <- grep(paste("w", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]
    wsamples <- matrix(samples[, snames %in% pn], ncol = (kx + kr))
    psamples <- cbind(psamples, "w" = wsamples)
  }

  ## Compute final smooth term object.
  fst <- compute_term(x, get.X = get.X, get.mu = get.mu,
    psamples = psamples, vsamples = vsamples, FUN = NULL, snames = snames,
    effects.hyp = NULL, fitted.values = NULL,
    data = data, grid = 100)

  attr(fst$term, "specs")$get.mu <- get.mu 

  rval <- list("term" = fst$term, "effects.hyp" = fst$effects.hyp, "fitted.values" = fst$fitted.values)
}


## Random scaling results function.
resultsJAGS.special.rsc.smooth <- function(x, samples, data, i, ...)
{
  snames <- colnames(samples)

  ## Prediction matrix
  class(x) <- x$class
  X <- PredictMat(x, data)

  ## Get coefficient samples of smooth term.
  xsamples <- rsamples <- NULL
  kr <- if(is.null(x$rand$Xr)) 0 else ncol(x$rand$Xr)
  kx <- if(is.null(x$Xf)) 0 else ncol(x$Xf)
  kw <- 0
  if(kx) {
    pn <- grep(paste("b", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]
    xsamples <- matrix(samples[, snames %in% pn], ncol = kx)
  }
  if(kr) {
    pn <- grep(paste("g", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]
    rsamples <- as.matrix(samples[, snames %in% pn], ncol = kr)
  }
  psamples <- cbind("ra" = rsamples, "fx" = xsamples)
  
  ## Retransform parameter samples.
  if(kr) {
    re_trans <- function(g) {
      g <- x$trans.D * g
      if(!is.null(x$trans.U))
        g <- x$trans.U %*% g
      g
    }
    psamples <- t(apply(psamples, 1, re_trans))
  }

  vsamples <- NULL
  taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
  if(taug %in% snames) {
    vsamples <- as.numeric(samples[, snames %in% taug])
  }

  get.mu <- x$get.mu

  if(!is.null(x$by.formula)) {

    rval <- list()
    for(j in seq_along(x$rs.by)) {
      x$by <- x$by.vars[j]
      by.levels <- levels(x$rs.by[[j]])
      for(jj in 1:length(by.levels)) {
        pn <- paste("rs", i, j, "[", jj, "]", sep = "")
        rsamples <- as.matrix(samples[, snames %in% pn], ncol = length(pn))
        rsamples <- cbind(rsamples, psamples)
        fsamples <- apply(rsamples, 1, function(g) { get.mu(X, g) })
        x2 <- x
        x2$label <- paste(x$label, ":", x$by.vars[j], by.levels[jj], sep = "")
        x2$by.level <- by.levels[jj]
        fst <- compute_term(x2, fsamples = fsamples, psamples = rsamples,
          vsamples = vsamples, FUN = NULL, snames = snames,
          effects.hyp = NULL, fitted.values = NULL, data = data)

        rval[[x2$label]] <- fst
      }
      attr(rval, "by") <- x$by.vars[1]
    }

  } else {
    fst <- compute_term(x2, fsamples = fsamples, psamples = psamples,
      vsamples = NULL, FUN = NULL, snames = snames,
      effects.hyp = NULL, fitted.values = NULL, data = data)

    rval <- list("term" = fst$term, "effects.hyp" = fst$effects.hyp, "fitted.values" = fst$fitted.values)
  }

  rval
}


## Special results function for growth curves.
resultsJAGS.special.gc.smooth <- function(x, samples, data, i, ...)
{
  snames <- colnames(samples)

  ## Prediction matrix.
  X <- PredictMat(x, data)
  get.mu <- x$get.mu
  vsamples <- NULL

  args <- list(...)
  id <- args$id

  if(!is.null(x$by.levels)) {
    ## Possible variance parameter samples.
    taug <- paste("taug", if(is.null(x$id)) i else x$id, id, sep = "")
    if(any(grepl(taug, snames, fixed = TRUE))) {
      vsamples <- as.matrix(samples[, grepl(taug, snames, fixed = TRUE)])
    }
    rval <- list()
    x$by <- x$byname
    x$byname <- NULL
    ## by <- as.factor(data[[x$by]])
    for(j in seq_along(x$by.levels)) {
      pn <- c(paste("g", i, id, "[", 1:3, "]", sep = ""), paste("g", i, id, "r[", j, ",", 1:3, "]", sep = ""))
      psamples <- as.matrix(samples[, snames %in% pn], ncol = length(pn))
      x2 <- x
      x2$label <- paste(x$label, ":", x$by, x$by.levels[j], sep = "")
      x2$by.levels <- NULL
      x2$by.level <- x$by.levels[j]

      fst <- compute_term(x2, get.X = function(x) { as.matrix(x) }, get.mu = get.mu,
        psamples = psamples, vsamples = vsamples, FUN = NULL, snames = snames,
        effects.hyp = NULL, fitted.values = NULL,
        data = data, grid = NA)

      rval[[x2$label]] <- fst
    }
    
    attr(rval, "by") <- x$by
  } else {
    pn <- grep(paste("g", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]

    psamples <- as.matrix(samples[, snames %in% pn], ncol = length(pn))

    ## Possible variance parameter samples.
    taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
    if(taug %in% snames) {
      vsamples <- as.matrix(samples[, snames %in% taug])
    }

    fst <- compute_term(x, get.X = function(x) { as.matrix(x) }, get.mu = get.mu,
      psamples = psamples, vsamples = vsamples, FUN = NULL, snames = snames,
      effects.hyp = NULL, fitted.values = NULL,
      data = data[, x$term, drop = FALSE], grid = 100)

    rval <- list("term" = fst$term, "effects.hyp" = fst$effects.hyp, "fitted.values" = fst$fitted.values)
  }
  rval
}


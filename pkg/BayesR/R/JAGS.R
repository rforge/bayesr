######################################################################
## (1) General JAGS setup functions, model code, initials and data. ##
######################################################################
## Setup the model structure needed for the sampler.
## The function creates the model code, initials and data
## to run a JAGS sampler.
setupJAGS <- function(x)
{
  family <- x$family
  if(!is.character(family))
    family <- as.character(family)[1]
  for(char in c("(", ")"))
    family <- gsub(char, "", family, fixed = TRUE)

  ## Start actual setup.
  setup <- list()
  setup$inits <- list("sigma" = runif(1))
  setup$psave <- "sigma"
 
  ## Parametric part.
  if(k <- ncol(x$X)) {
    setup$data$X <- x$X
    if(all(x$X[, 1] == 1))
      attr(setup, "intercept") <- TRUE
    for(j in 1:k) {
      setup$param <- c(setup$param, paste(if(k < 2) "beta" else paste("beta[", j, "]", sep = ""),
        " * X[i, ", j, "]", sep = ""))
    }
    setup$param <- paste(setup$param, collapse = " + ")
    setup$param <- paste("    param[i] <-", setup$param)
    setup$loops <- k
    setup$priors.coef <- if(k > 1) "    beta[j] ~ dnorm(0, sigma)" else "  beta ~ dnorm(0, sigma)"
    setup$inits$beta <- runif(k)
    setup$psave <- c(setup$psave, "beta")
    setup$eta <- "param[i]"
  }
  attr(setup, "intercept") <- attr(terms(x$fake.formula), "intercept")
  attr(setup, "k") <- k
  
  ## Smooth part.
  if(m <- length(x$smooth)) {
    for(i in 1:m) {
      setup <- if(!is.null(x$smooth[[i]]$special)) {
        buildJAGS.smooth.special(x$smooth[[i]], setup, i)
      } else {
        buildJAGS.smooth(x$smooth[[i]], setup, i)
      }
    }
  }

  ## Setup complete model.
  setup$eta <- paste("    eta[i] <-", setup$eta)
  if(!is.na(pmatch(tolower(family), "gaussian"))) {
    response <- "    response[i] ~ dnorm(eta[i], sigma)"
    setup$priors.scale <- c(setup$priors.scale, "  sigma ~ dgamma(1.0E-6, 1.0E-6)")
  }
  model <- c(
    "model {",
    setup$start,
    "  for(i in 1:n) {",
    setup$param, setup$smooth,
    if(is.null(setup$adds)) c(setup$eta, response),
    "  }"
  )
  if(!is.null(setup$adds)) {
    model <- c(model,
      "  for(i in 1:n) {",
      setup$adds,
      setup$eta,
      response,
      "  }"
    )
  }

  lp <- list()
  for(j in 1:length(setup$loops))
    lp[[paste(setup$loops[j])]] <- c(lp[[paste(setup$loops[j])]], setup$priors.coef[j])
  for(j in names(lp)) {
    if(j != 1)
      tmp <- c(paste("  for(j in 1:", j, ") {", sep = ""), lp[[j]], "  }")
    else
      tmp <- lp[[j]]
    model <- c(model, tmp)
  }
  model <- c(model, setup$priors.scale, setup$close, setup$close2, "}")

  ## writeLines(model)

  ## Final touch ups.
  setup$data$response <- x$mf[[x$response]]
  setup$data$n <- length(setup$data$response)

  rval <- list("model" = model, "data" = setup$data,
    "inits" = setup$inits, "psave" = setup$psave)

  rval
}


## Build the JAGS model code for a smooth term. 
buildJAGS.smooth <- function(smooth, setup, i) {
  fall <- NULL
  kr <- if(is.null(smooth$rand$Xr)) 0 else ncol(smooth$rand$Xr)
  kx <- if(is.null(smooth$Xf)) 0 else ncol(smooth$Xf)
  if(kx > 0) {
    fall <- c(fall, paste("b", i, if(kx > 1) paste("[", 1:kx, "]", sep = ""),
      "*Xf", i, "[i, ", 1:kx, "]", sep = ""))
    setup$data[[paste("Xf", i, sep = "")]] <- smooth$Xf
    tmp <- if(kx > 1) {
        paste("    b", i, "[j] ~ dnorm(0, sigma)", sep = "")
    } else paste("  b", i, " ~ dnorm(0, sigma)", sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kx)
    setup$inits[[paste("b", i, sep = "")]] <- runif(kx)
    setup$psave <- c(setup$psave, paste("b", i, sep = ""))
  }
  if(kr > 0) {
    fall <- c(fall, paste("g", i, if(kr > 1) paste("[", 1:kr, "]", sep = ""), "*Xr",
      i, "[i, ", 1:kr, "]", sep = ""))
    setup$data[[paste("Xr", i, sep = "")]] <- smooth$rand$Xr
    taug <- paste("taug", if(is.null(smooth$id)) i else smooth$id, sep = "")
    tmp <- if(kr > 1) {
      paste("    g", i, "[j] ~ dnorm(0, ", taug, ")", sep = "")
    } else paste("g", i, " ~ dnorm(0, ", taug, ")", sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kr)
    if(is.null(setup$priors.scale) || !any(grepl(taug, setup$priors.scale))) {
      setup$priors.scale <- c(setup$priors.scale, paste("  ", taug, " ~ dgamma(1.0E-6, 1.0E-6)", sep = ""))
      setup$inits[[taug]] <- runif(1, 0.00001, 0.0001)
      setup$psave <- c(setup$psave, taug)
    }
    setup$inits[[paste("g", i, sep = "")]] <- runif(kr)
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
buildJAGS.smooth.special <- function(smooth, setup, i)
{
  UseMethod("buildJAGS.smooth.special")
}


## Default special model term builder.
buildJAGS.smooth.special.default <- function(smooth, setup, i)
{
  buildJAGS.smooth(smooth, setup, i)
}


## JAGS random scaling model term constructor.
buildJAGS.smooth.special.rs.smooth <- function(smooth, setup, i)
{
  smooth$special <- FALSE
  setup <- buildJAGS.smooth(smooth, setup, i)

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
buildJAGS.smooth.special.gc.smooth <- function(smooth, setup, i)
{
  pn <- paste("g", i, sep = "")

  setup$data[[paste("X", pn, sep = "")]] <- as.numeric(smooth$X)
  setup$inits[[paste(pn, "0", sep = "")]] <- runif(3)
  setup$psave <- c(setup$psave, pn)

  setup$close <- c(setup$close,
    "  for(j in 1:3) {",
    paste("    ", pn, "[j] <- exp(", pn, "0[j])", sep = ""),
    paste("    ", pn, "0[j] ~ dnorm(0, sigma)", sep = "")
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
    setup$close <- c(setup$close,
    paste("    for(k in 1:", length(smooth$by.levels), ") {", sep = ""),
      paste("      ", pn, "r[k, j] <- exp(", pn, "r0[k, j])", sep = ""),
      paste("      ", pn, "r0[k, j] ~ dnorm(0, taug", i, "[j])", sep = ""),
      "    }",
      paste("    taug", i, "[j] ~ dgamma(1.0E-3, 1.0E-3)", sep = "")
    )
    setup$psave <- c(setup$psave, paste(pn, "r", sep = ""), paste("taug", i, sep = ""))
    setup$inits[[paste(pn, "r0", sep = "")]] <- matrix(runif(length(smooth$by.levels) * 3), ncol = 3)
    setup$inits[[paste("taug", i, sep = "")]] <- runif(3, 0.0001, 0.001)
  }

  setup$close <- c(setup$close, "  }")

  setup$smooth <- c(setup$smooth, paste("    sm", i, "[i] <- ",
    paste(fall, collapse = " + ", sep = ""), sep = ""))
  setup$eta <- paste(setup$eta, paste("sm", i, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


########################################
## (2) Interface to the JAGS sampler. ##
########################################
samplerJAGS <- function(x, tdir = NULL,
  n.chains = 1, n.adapt = 100,
  n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, ...)
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
  
  jmodel <- jags.model(mfile, data = x$data, inits = inits,
    n.chains = n.chains, n.adapt = n.adapt, ...)
  jsamples <- coda.samples(jmodel, variable.names = c(x$psave, "deviance"),
    n.iter = n.iter, thin = thin, ...)

  ## Remove burnin.
  if(is.null(burnin))
    burnin <- floor(n.iter * 0.2)
  jsamples <- window(jsamples, start = burnin)

  jsamples
}


###################################################################
## (3) Functions to compute summary statistics, plots, etc. from ##
##     samples returned the JAGS sampler.                        ##
###################################################################
## Function to extract all results obtained by running the JAGS
## sampler. The function uses BayesR structures to represent fitted
## model terms etc.
resultsJAGS <- function(x, samples)
{
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
    fixed.effects <- effects <- smooth.hyp <- NULL
    fitted.values <- 0

    ## Fixed effects.
    if(k <- ncol(x$X)) {
      samps <- as.matrix(samples[[j]][, grepl("beta", snames)], ncol = k)
      nx <- colnames(x$X)
      for(i in seq_along(nx)) {
        if(nx[i] != "(Intercept)") {
          if(is.null(effects))
            effects <- list()
          Z <- matrix(x$X[, i])
          if(grepl("(Intercept)", nx[i])) {
            Z <- cbind(x$X[, 1], Z)
            f <- apply(matrix(cbind(samps[, 1], samps[, i]), ncol = 1), 1, function(x) Z %*% x)
          } else {
            f <- apply(matrix(samps[, i], ncol = 1), 1, function(x) Z %*% x)
          }
          f <- t(apply(f, 1, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)))
          fm <- apply(f, 1, mean)
          fitted.values <- fitted.values + fm
          f <- f - mean(fm)
          fm <- fm - mean(fm)
          effects[[nx[i]]] <- as.data.frame(unique(cbind(Z, fm, f)))
          colnames(effects[[nx[i]]]) <- c(nx[i], "Mean", "2.5%", "10%", "50%", "90%", "97.5%")
          class(effects[[nx[i]]]) <- "data.frame"
          attr(effects[[nx[i]]], "sample") <- samps[, i]
          attr(effects[[nx[i]]], "specs") <- list(dim = 1, label = nx[i])
          attr(effects[[nx[i]]], "fit") <- fm
          attr(effects[[nx[i]]], "x") <- x$X[, i]
        }
      }
      qu <- t(apply(samps, 2, quantile, probs = c(0.025, 0.5, 0.975)))
      sd <- drop(apply(samps, 2, sd))
      me <- drop(apply(samps, 2, mean))
      fixed.effects <- cbind(me, sd, qu)
      rownames(fixed.effects) <- nx
      colnames(fixed.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
      if(grepl("(Intercept)", nx[i])) {
        fitted.values <- fitted.values + fixed.effects[1, 1]
      }
      attr(fixed.effects, "samples") <- samps
      colnames(attr(fixed.effects, "sample")) <- nx
    }
  
    ## Smooth terms.
    if(length(x$smooth)) {
      if(!is.list(effects))
        effects <- list()
      for(i in 1:length(x$smooth)) {
        if(!is.null(x$smooth[[i]]$special)) {
          fst <- resultsJAGS.special(x$smooth[[i]], samples[[j]], x$mf, i)
          if(is.null(attr(fst, "by"))) {
            effects[[x$smooth[[i]]$label]] <- fst$term
            smooth.hyp <- rbind(smooth.hyp, fst$smooth.hyp)
            fitted.values <- fitted.values + fst$fitted.values
          } else {
            tjl <- names(fst)
            for(l in tjl) {
              effects[[l]] <- fst[[l]]$term
              smooth.hyp <- rbind(smooth.hyp, fst[[l]]$smooth.hyp)
              fitted.values <- fitted.values + fst[[l]]$fitted.values
            }
          }
          rm(fst)
        } else {
          ## Get coefficient samples of smooth term.
          xsamples <- rsamples <- NULL
          kr <- if(is.null(x$smooth[[i]]$rand$Xr)) 0 else ncol(x$smooth[[i]]$rand$Xr)
          kx <- if(is.null(x$smooth[[i]]$Xf)) 0 else ncol(x$smooth[[i]]$Xf)
          kw <- 0
          if(kx) {
            pn <- grep(paste("b", i, sep = ""), snames, value = TRUE, fixed = TRUE)
            pn <- pn[!grepl("tau", pn)]
            xsamples <- matrix(samples[[j]][, snames %in% pn], ncol = kx)
          }
          if(kr) {
            pn <- grep(paste("g", i, sep = ""), snames, value = TRUE, fixed = TRUE)
            pn <- pn[!grepl("tau", pn)]
            rsamples <- as.matrix(samples[[j]][, snames %in% pn], ncol = kr)
          }
          psamples <- cbind("ra" = rsamples, "fx" = xsamples)
  
          ## Retransform parameter samples.
          if(kr) {
            re_trans <- function(g) {
              g <- x$smooth[[i]]$trans.D * g
              if(!is.null(x$smooth[[i]]$trans.U))
                g <- x$smooth[[i]]$trans.U %*% g
              g
            }
            psamples <- t(apply(psamples, 1, re_trans))
          }

          ## Prediction matrix.
          X <- PredictMat(x$smooth[[i]], x$mf)

          ## Possible variance parameter samples.
          vsamples <- NULL
          taug <- paste("taug", if(is.null(x$smooth[[i]]$id)) i else x$smooth[[i]]$id, sep = "")
          if(taug %in% snames) {
            vsamples <- as.numeric(samples[[j]][, snames %in% taug])
          }

          get.mu <- function(X, g) {
            X %*% as.numeric(g)
          }

          ## Compute samples of fitted values.
          fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })

          ## Compute final smooth term object.
          fst <- compute_term(x$smooth[[i]], fsamples = fsamples, psamples = psamples,
            vsamples = vsamples, FUN = NULL, snames = snames,
            smooth.hyp = smooth.hyp, fitted.values = fitted.values, data = x$mf)

          attr(fst$term, "specs")$get.mu <- get.mu 

          ## Add term to effects list.
          effects[[x$smooth[[i]]$label]] <- fst$term
          smooth.hyp <- fst$smooth.hyp

          fitted.values <- fst$fitted.values
          rm(fst)
        }
      }
    }

    ## Scale parameter.
    samps <- 1 / as.numeric(samples[[j]][, snames %in% "sigma"])
    qu <- drop(quantile(samps, probs = c(0.025, 0.5, 0.975)))
    sd <- sd(drop(samps))
    me <- mean(samps)
    variance <- matrix(c(me, sd, qu), nrow = 1)
    colnames(variance) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
    rownames(variance) <- "Sigma2"
    samps <- matrix(samps, ncol = 1)
    colnames(samps) <- "scale"
    attr(variance, "samples") <- as.mcmc(samps)

    ## Compute partial residuals.
    for(i in seq_along(effects)) {
      e <- x$mf[[x$response]] - (fitted.values - attr(effects[[i]], "fit"))
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

    ## Stuff everything together.
    rval[[j]] <- list("fixed.effects" = fixed.effects, "effects" = effects,
      "smooth.hyp" = smooth.hyp, "model.fit" = list("method" = "MCMC", "DIC" = DIC, "pd" = pd,
      "N" = nrow(x$mf), "formula" = x$formula), "variance" = variance, "call" = x$call,
      "fitted.values" = fitted.values, "residuals" = x$mf[[x$response]] - fitted.values)

    class(rval[[j]]) <- c("bayesr", "bayesx")
  }
  names(rval) <- paste("Chain", 1:chains, sep = "_")
  if(length(rval) < 2) {
    rval <- rval[[1]]
  }
  class(rval) <- c("bayesr", "bayesx")

  rval
}


## Result extractor function for special terms.
resultsJAGS.special <- function(x, samples, data, i, ...) 
{
  UseMethod("resultsJAGS.special")
}


## Default special term results method.
resultsJAGS.special.default <- function(x, samples, data, i, ...)
{
  warning("Please fixme!")
  snames <- colnames(samples)

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

  ## Prediction matrix.
  X <- PredictMat(x, data)

  ## Possible variance parameter samples.
  vsamples <- NULL
  taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
  if(taug %in% snames) {
    vsamples <- as.numeric(samples[, snames %in% taug])
  }

  get.mu <- function(X, g) {
    X %*% as.numeric(g)
  }

  ## Compute samples of fitted values.
  fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })

  ## Compute final smooth term object.
  fst <- compute_term(x, fsamples = fsamples, psamples = psamples,
    vsamples = vsamples, FUN = NULL, snames = snames,
    smooth.hyp = NULL, fitted.values = NULL, data = data)

  attr(fst$term, "specs")$get.mu <- get.mu 

  rval <- list("term" = fst$term, "smooth.hyp" = fst$smooth.hyp, "fitted.values" = fst$fitted.values)
}


## Random scaling results function.
resultsJAGS.special.rs.smooth <- function(x, samples, data, i, ...)
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
          smooth.hyp = NULL, fitted.values = NULL, data = data)

        rval[[x2$label]] <- fst
      }
      attr(rval, "by") <- x$by.vars[1]
    }

  } else {
    fst <- compute_term(x2, fsamples = fsamples, psamples = psamples,
      vsamples = NULL, FUN = NULL, snames = snames,
      smooth.hyp = NULL, fitted.values = NULL, data = data)

    rval <- list("term" = fst$term, "smooth.hyp" = fst$smooth.hyp, "fitted.values" = fst$fitted.values)
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

  if(!is.null(x$by.levels)) {
    ## Possible variance parameter samples.
    taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
    if(any(grepl(taug, snames, fixed = TRUE))) {
      vsamples <- as.matrix(samples[, grepl(taug, snames, fixed = TRUE)])
    }
    rval <- list()
    x$by <- x$byname
    x$byname <- NULL
    ## by <- as.factor(data[[x$by]])
    for(j in seq_along(x$by.levels)) {
      pn <- c(paste("g", i, "[", 1:3, "]", sep = ""), paste("g", i, "r[", j, ",", 1:3, "]", sep = ""))
      psamples <- as.matrix(samples[, snames %in% pn], ncol = length(pn))
      fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })
      x2 <- x
      x2$label <- paste(x$label, ":", x$by, x$by.levels[j], sep = "")
      x2$by.levels <- NULL
      x2$by.level <- x$by.levels[j]
      fst <- compute_term(x2, fsamples = fsamples, psamples = psamples,
        vsamples = NULL, FUN = NULL, snames = snames,
        smooth.hyp = NULL, fitted.values = NULL, data = data)

      rval[[x2$label]] <- fst
    }
    
    attr(rval, "by") <- x$by
  } else {
    pn <- grep(paste("g", i, sep = ""), snames, value = TRUE, fixed = TRUE)
    pn <- pn[!grepl("tau", pn)]

    psamples <- as.matrix(samples[, snames %in% pn], ncol = length(pn))

    ## Compute samples of fitted values.
    fsamples <- apply(psamples, 1, function(g) { get.mu(X, g) })

    ## Possible variance parameter samples.
    taug <- paste("taug", if(is.null(x$id)) i else x$id, sep = "")
    if(taug %in% snames) {
      vsamples <- as.matrix(samples[, snames %in% taug])
    }

    fst <- compute_term(x, fsamples = fsamples, psamples = psamples,
      vsamples = vsamples, FUN = NULL, snames = snames,
      smooth.hyp = NULL, fitted.values = NULL, data = data)

    rval <- list("term" = fst$term, "smooth.hyp" = fst$smooth.hyp, "fitted.values" = fst$fitted.values)
  }
  rval
}


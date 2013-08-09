######################################################################
## (1) General JAGS setup functions, model code, initials and data. ##
######################################################################
## Setup the model structure needed for the sampler.
## The function creates the model code, initials and data
## to run a JAGS sampler.
restructure_x <- function(x) {
  if(!all(c("family", "formula", "response") %in% names(x))) {
    call <- x$call
    x <- x[names(x) != "call"]
    x2 <- list(); k <- 1
    for(i in seq_along(x)) {
      if(!all(c("family", "formula", "response") %in% names(x[[i]]))) {
        x2[[k]] <- x[[i]]
        x[[i]] <- NULL
      }
    }
    if(length(x2))
      x <- c(x, x2)
    if(any(duplicated(nx <- names(x))))
      names(x) <- paste(nx, 1:length(nx), sep = "")
    x$call <- call
  }
  x
}

## Transform design and penalty matrices to iid random effects structure.
transformJAGS <- function(x) {
  x <- randomize(x)
  x <- restructure_x(x)
  x
}

## Default linear predictor and model setup functions.
JAGSeta <- function(x, id = NULL, ...) {
  setup <- list()
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
      paste("    beta", id2, "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  beta", id2, " ~ dnorm(0, 1.0E-6)", sep = "")
    setup$inits[[paste("beta", id2, sep = "")]] <- runif(k)
    setup$psave <- c(setup$psave, paste("beta", id2, sep = ""))
    setup$eta <- paste("param", id2, "[i]", sep = "")
  }

  ## Smooth part.
  if(m <- length(x$smooth)) {
    for(i in 1:m) {
      setup <- if(!is.null(x$smooth[[i]]$special)) {
        buildJAGS.smooth.special(x$smooth[[i]], setup, paste(i, id, sep = ""))
      } else {
        buildJAGS.smooth(x$smooth[[i]], setup, paste(i, id, sep = ""))
      }
    }
  }

  ## Final touch ups.
  if(!is.null(x$mf[[x$response]]))
    setup$data$response <- x$mf[[x$response]]

  setup
}

## Get link functions.
JAGSlinks <- function(x)
{
  switch(x,
    "identity" = "eta",
    "log" = "exp(eta)",
    "exp" = "log(eta)",
    "inverse" = "1 / (eta)",
    "logit" = "1 / (1 + exp(-(eta)))",
    "probit" = "phi(eta)"
  )
}

## Construct the final model code.
JAGSmodel <- function(x, family, ...) {
  if(is.function(family))
    family <- family()
  k <- if(all(c("inits", "data", "psave") %in% names(x))) {
    x <- list(x)
    1
  } else length(x)
  if(k > family$k) {
    stop(paste("more parameters specified than existing in family ",
      family$family, ".BayesR()!", sep = ""), call. = FALSE)
  }

  model <- "model {"
  for(j in 1:k) {
    model <- c(model, x[[j]]$start)
  }
  pn <- family$names
  if(is.infinite(family$k))
    family$k <- k
  if(is.null(pn)) pn <- paste("theta", 1:family$k, sep = "")
  if(length(pn) < 2 & length(pn) != k)
    pn <- paste(pn, 1:k, sep = "")
  pn[1:k] <- paste(pn[1:k], "[i]", sep = "")
  on <- family$oname
  links <- family[grep("link", names(family), fixed = TRUE, value = TRUE)]
  links <- rep(sapply(links, JAGSlinks), length.out = k)
  model <- c(model,  "  for(i in 1:n) {",
    paste("    response[i] ~ ", family$dist, "(",
      paste(if(is.null(on)) pn else paste(on, "[i, ]", sep = ""), collapse = ", "), ")", sep = ""))
  for(j in 1:k) {
    model <- c(model, paste("    ", if(is.null(on)) pn[j] else paste(on, "[i, ", j, "]", sep = ""),
      " <- ", gsub("eta", x[[j]]$eta, links[[j]]), sep = ""))
  }
  for(j in 1:k)
    model <- c(model, x[[j]]$adds)
  for(j in 1:k)
    model <- c(model, x[[j]]$param, x[[j]]$smooth)
  model <- c(model, "  }")

  for(i in 1:k) {
    lp <- list()
    for(j in 1:length(x[[i]]$loops))
      lp[[paste(x[[i]]$loops[j])]] <- c(lp[[paste(x[[i]]$loops[j])]], x[[i]]$priors.coef[j])
    for(j in names(lp)) {
      if(j != 1)
        tmp <- c(paste("  for(j in 1:", j, ") {", sep = ""), lp[[j]], "  }")
      else
        tmp <- lp[[j]]
      model <- c(model, tmp)
    }
    model <- c(model, x[[i]]$priors.scale, x[[i]]$close, x[[i]]$close2)
  }
  if(k < family$k) {
    model <- c(model, paste(" ", family$default.prior[(k + 1):family$k]))
  }
  model <- c(model, "}")
  if(k < family$k) {
    attr(model, "psave") <- family$names[(k + 1):family$k]
  }

  model
}

## Create final model setup.
setupJAGS <- function(x)
{
  x <- restructure_x(x)
  ncat <- NULL
  if(!all(c("family", "formula", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx))
      nx <- 1:length(x)
    rval <- list()
    family <- if(is.function(x[[1]]$family)) x[[1]]$family() else x[[1]]$family
    fn <- family$names
    for(i in seq_along(nx)) rval[[nx[i]]] <- family$eta(x[[i]], fn[i])
  } else {
    family <- if(is.function(x$family)) x$family() else x$family
    rval <- family$eta(x)
  }
  
  ## Create model code.
  model <- family$model(rval, family)

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
  data$n <- nrow(if(all(c("family", "formula", "response") %in% names(x))) x$mf else x[[1]]$mf)
  psave <- c(psave, attr(model, "psave"))
#  if(!is.null(fhamily$oname))
#    inits[[family$oname]] <- matrix(0, data$n, length(rval))

  rval <- list("model" = model, "data" = data,
    "inits" = inits, "psave" = psave)
  
  return(rval)
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
        paste("    b", i, "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  b", i, " ~ dnorm(0, 1.0E-6)", sep = "")
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
  seed = NULL, verbose = TRUE, ...)
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
resultsJAGS <- function(x, samples, id = NULL)
{
  if(!all(c("family", "formula") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    rval <- list()
    family <- x[[1]]$family
    if(is.function(family))
      family <- family()
    fn <- family$names
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- resultsJAGS(x[[nx[j]]], samples, id = fn[j])
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
    class(rval) <- "bayesr"
    return(rval)
  } else {
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
      if(k <- ncol(x$X)) {
        samps <- as.matrix(samples[[j]][, grepl(paste("beta", id, sep = ""), snames)], ncol = k)
        nx <- colnames(x$X)
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
      if(length(x$smooth)) {
        if(!is.list(effects))
          effects <- list()
        for(i in 1:length(x$smooth)) {
          if(!is.null(x$smooth[[i]]$special)) {
            fst <- resultsJAGS.special(x$smooth[[i]], samples[[j]], x$mf, i)
            if(is.null(attr(fst, "by"))) {
              effects[[x$smooth[[i]]$label]] <- fst$term
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
            kr <- if(is.null(x$smooth[[i]]$rand$Xr)) 0 else ncol(x$smooth[[i]]$rand$Xr)
            kx <- if(is.null(x$smooth[[i]]$Xf)) 0 else ncol(x$smooth[[i]]$Xf)
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
            taug <- paste("taug", if(is.null(x$smooth[[i]]$id)) i else x$smooth[[i]]$id, id, sep = "")
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
              effects.hyp = effects.hyp, fitted.values = fitted.values, data = x$mf)

            attr(fst$term, "specs")$get.mu <- get.mu 

            ## Add term to effects list.
            effects[[x$smooth[[i]]$label]] <- fst$term
            effects.hyp <- fst$effects.hyp

            fitted.values <- fst$fitted.values
            rm(fst)
          }
        }
      }

      ## Scale parameters.
      scale.m <- scale.samps.m <- NULL
      sn <- if(is.null(id)) {
        family <- if(is.function(x$family)) x$family() else x$family
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
      if(x$response %in% names(x$mf)) {
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
      }

      ## Stuff everything together.
      rval[[j]] <- list("call" = x$call, "family" = x$family,
        "model" = list("DIC" = DIC, "pd" = pd, "N" = nrow(x$mf),
        "formula" = x$formula), "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "scale" = scale.m, "fitted.values" = fitted.values,
        "residuals" = x$mf[[x$response]] - fitted.values)
      
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
    effects.hyp = NULL, fitted.values = NULL, data = data)

  attr(fst$term, "specs")$get.mu <- get.mu 

  rval <- list("term" = fst$term, "effects.hyp" = fst$effects.hyp, "fitted.values" = fst$fitted.values)
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
        effects.hyp = NULL, fitted.values = NULL, data = data)

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
      effects.hyp = NULL, fitted.values = NULL, data = data)

    rval <- list("term" = fst$term, "effects.hyp" = fst$effects.hyp, "fitted.values" = fst$fitted.values)
  }
  rval
}


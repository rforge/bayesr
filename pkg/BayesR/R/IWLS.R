## Setup for IWLS sampler, handling
## sampling functions.
transformIWLS <- function(x, ...)
{
  call <- x$call; x$call <- NULL

  tIWLS <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- tIWLS(obj[[j]], ...)
    } else {
      if(!is.null(dim(obj$X))) {
        if(is.null(obj$smooth)) obj$smooth <- list()
        obj$smooth[["parametric"]] <- list(
          "X" = obj$X,
          "S" = list(diag(0, ncol(obj$X))),
          "rank" = ncol(obj$X),
          "term" = "linear",
          "label" = "linear",
          "bs.dim" = ncol(obj$X),
          "fixed" = TRUE,
          "by" = "NA",
          "is.linear" = TRUE
        )
        obj$sterms <- c(obj$strems, "parametric")
        obj$X <- NULL
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]] <- smooth.IWLS(obj$smooth[[j]])
        }
      }
    }
    obj
  }

  x <- tIWLS(x, ...)

  attr(x, "call") <- call
  attr(x, "response.vec") <- attr(x, "model.frame")[, attr(attr(x, "model.frame"), "response.name")]

  x
}


## Function to setup IWLS smooths.
smooth.IWLS <- function(x, ...) {
  UseMethod("smooth.IWLS")
}

smooth.IWLS.default <- function(x, ...)
{
  if(is.null(x$get.mu)) {
    x$get.mu <- if(is.null(x$get.mu)) {
      function(X, g) {
        X %*% as.numeric(g)
      }
    } else x$get.mu
  }

  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b

  if(is.null(x$state)) {
    x$p.save <- c("g", "tau2")
    x$state <- list()
    x$state$g <- runif(ncol(x$X), 1.001, 1.002)
    if(!x$fixed)
      x$state$tau2 <- if(is.null(x$sp)) runif(1, 0.99, 1) else x$sp
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("c", 1:length(x$state$g), sep = ""),
        if(!x$fixed) "tau2" else NULL)
    } else x$s.colnames
    x$np <- length(x$s.colnames)
  }

  if(is.null(x$propose)) {
    if(is.null(x$rand)) {
      x$propose <- function(x, family, response, eta, id, ...) {
        ## Compute weights.
        weights <- family$weights[[id]](response, eta)

        ## Score.
        score <- family$score[[id]](response, eta)

        ## Compute working observations.
        z <- eta[[id]] + 1 / weights * score

#        ## Save old predictor.
#        eta0 <- eta[[id]]

        ## Compute old log likelihood and old log coefficients prior.
        pibeta <- family$loglik(response, eta)
#        p1 <- if(x$fixed) {
#          0
#        } else drop(-0.5 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)

        ## Compute partial predictor.
        eta[[id]] <- eta[[id]] - x$state$fit

#        ## Compute mean and precision.
#        XW <- t(x$X * weights)
#        P <- if(x$fixed) {
#          chol2inv(chol(P0 <- XW %*% x$X))
#        } else chol2inv(chol(P0 <- XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
#        M <- P %*% (XW %*% (z - eta[[id]]))

        ## Save old coefficients
        g0 <- x$state$g

#        ## Sample new parameters.
#        x$state$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

        ## Start Nadja.
        coeff <- g0
        pr <- t(x$X) %*% (x$X * weights) + if(!x$fixed) 1 / x$state$tau2 * x$S[[1]] else 0
        cholpr <- chol(pr, symmetric = TRUE)
        zufall <- rnorm(ncol(x$X), mean = 0, sd = 1)
        coeffprop <- solve(cholpr, zufall) ## N(0, P^-1)
        vec <- solve(t(cholpr), t(x$X * weights) %*% (z - eta[[id]]))
        mu <- solve(cholpr, vec)
        coeffprop <- coeffprop + mu
        if(!x$fixed) {
          hilfsvec <- solve(t(cholpr), t(rep(1, length(score)) %*% x$X))
          v <- solve(cholpr, hilfsvec)
          w <- drop(rep(1, length(score)) %*% x$X %*% v)
          x$state$g <- drop(coeffprop - t(1/w * t(v)) %*% (rep(1, length(score)) %*% x$X %*% coeffprop))
        }
        x$state$fit <- x$X %*% x$state$g
        ## End Nadja.

#        ## Compute log priors.
#        p2 <- if(x$fixed) {
#          0
#        } else drop(-0.5 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)
#        g1 <- x$state$g - M
#        qbetaprop <- 0.5 * sum(log((diag(chol(P0, symmetric = TRUE))^2))) -0.5 * crossprod(g1, P0) %*% g1
#        qbetaprop2 <- dmvnorm(x$state$g, mean = M, sigma = P, log = TRUE)

#        ## Compute fitted values.        
#        x$state$fit <- drop(x$X %*% x$state$g)
#        x$state$fit <- x$state$fit - mean(x$state$fit)

        ## Set up new predictor.
        eta[[id]] <- eta[[id]] + x$state$fit

        ## Compute new log likelihood.
        pibetaprop <- family$loglik(response, eta)

        ## Compute new weights
        weights <- family$weights[[id]](response, eta)

        ## New score.
        score <- family$score[[id]](response, eta)

        ## New working observations.
        z <- eta[[id]] + 1 / weights * score

#        ## Compute partial predictor.
#        eta[[id]] <- eta[[id]] - x$state$fit

        ## Start Nadja.      
        prprop <- t(x$X) %*% (x$X * weights) + if(!x$fixed) 1 / x$state$tau2 * x$S[[1]] else 0
        cholprprop <- chol(prprop, symmetric = TRUE)
        vecprop <- solve(t(cholprprop), t(x$X * weights) %*% (z - (eta[[id]] - x$state$fit)))
        muprop <- solve(cholprprop, vecprop)
        qbeta <- 0.5 * sum(log((diag(cholprprop)^2))) - 0.5 * t(coeff - muprop) %*% prprop %*% (coeff - muprop)
        qbetaprop <- 0.5 * sum(log((diag(cholpr)^2))) - 0.5 * t(x$state$g - mu) %*% pr %*% (x$state$g - mu)
        p1 <- if(!x$fixed) - 1 / (2 * x$state$tau2) * t(coeff) %*% x$S[[1]] %*% coeff else 0
        p2 <- if(!x$fixed) - 1 / (2 * x$state$tau2) * t(x$state$g) %*% x$S[[1]] %*% x$state$g else 0
        ## End Nadja.

#        ## Compute mean and precision.
#        XW <- t(x$X * weights)
#        P2 <- if(x$fixed) {
#          chol2inv(chol(P0 <- XW %*% x$X))
#        } else chol2inv(chol(P0 <- XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
#        M2 <- P2 %*% (XW %*% (z - eta[[id]]))

#        ## Get the log prior.
#        g2 <- g0 - M2
#        qbeta <- 0.5 * sum(log((diag(chol(P0, symmetric = TRUE))^2))) -0.5 * crossprod(g2, P0) %*% g2
#        qbeta2 <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

        ## Sample variance parameter.
        if(!x$fixed & is.null(x$sp)) {
          a <- x$rank / 2 + x$a
          b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
          x$state$tau2 <- 1 / rgamma(1, a, b)
        }

        ## Compute acceptance probablity.
        x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

        return(x$state)
      }
    }
  }

  x
}


## Sampler based on IWLS proposals.
samplerIWLS <- function(x, n.iter = 1200, thin = 1, burnin = 200,
  verbose = TRUE, step = 100, tdir = NULL, ...)
{
  require("mvtnorm")

  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  response <- attr(x, "response.vec")

  ## Actual number of samples to save.
  if(n.iter < burnin) stop("argument burnin exceeds n.iter!")
  if(thin > (n.iter - burnin)) stop("argument thin is set too large!")
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))
  n.save <- length(iterthin)
  
  ## Add accptance rate and fitted values vectors.
  smIWLS <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- smIWLS(obj[[j]], ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]]$s.alpha <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.samples <- matrix(0, nrow = n.save, ncol = obj$smooth[[j]]$np)
          obj$smooth[[j]]$state$fit <- rep(0, nrow(obj$smooth[[j]]$X))
        }
      }
    }
    obj
  }

  x <- smIWLS(x, ...)

  ## Number of parameters
  np <- length(nx)

  ## Set up predictors.
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np)
    eta[[j]] <- rep(0.1, length(response))

  deviance <- rep(0, length(iterthin))

  ## Start sampling
  for(i in 1:n.iter) {
    if(save <- i %in% iterthin)
      js <- which(iterthin == i)
    
    ## Cycle through all parameters
    for(j in 1:np) {
      ## And all terms.
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        ## Get proposed states.
        p.state <- x[[nx[j]]]$smooth[[sj]]$propose(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j])

        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= min(c(p.state$alpha, 1))

print(accepted)
plot(p.state$fit ~ dat$x)
Sys.sleep(2)
accepted <- TRUE

        if(accepted) {
          eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit
          x[[nx[j]]]$smooth[[sj]]$state <- p.state 
        }

        ## Save the samples and acceptance.
        if(save) {
          x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- accepted
          x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$state[x[[nx[j]]]$smooth[[sj]]$p.save])
          deviance[js] <- -2 * family$loglik(response, eta)
        }
      }
    }
    if(verbose) {
      if(i %% step == 0)
        cat("iteration:", i, "\n")
    }
  }

  ## Return all samples as mcmc matrix.
  ## (1) Write out all samples to tdir.
tdir <- "~/tmp"
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  }
  tdir <- path.expand(tdir)

  samplesIWLS <- function(obj, id = NULL, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- paste("p", 1:length(obj), sep = "")
      if(length(unique(nx)) < length(obj)) nx <- paste("p", 1:length(obj), sep = "")
      for(j in nx)
        samplesIWLS(obj[[j]], id = j, ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h",
            obj$hlevel, ":", obj$smooth[[j]]$label, ".raw", sep = ""))
          colnames(obj$smooth[[j]]$s.samples) <- obj$smooth[[j]]$s.colnames
          write.table(obj$smooth[[j]]$s.samples, file = fn, row.names = FALSE, quote = FALSE)
        }
      }
    }
    NULL
  }

  samplesIWLS(x)

  ## (2) Remove old object.
  rm(x)

  ## (3) Read back in samples
  sf <- grep(".raw", dir(tdir), value = TRUE)
  samples <- NULL
  for(j in sf) {
    st <- as.matrix(read.table(file.path(tdir, j), header = TRUE))
    colnames(st) <- paste(gsub(".raw", "", j), gsub("c", "", colnames(st)), sep = ".")
    samples <- cbind(samples, st)
  }
  samples <- cbind(samples, "deviance" = deviance)

  return(as.mcmc(samples))
}


resultsIWLS <- function(x, samples)
{
  family <- attr(x, "family")
  grid <- attr(x, "grid")
  if(is.null(grid)) grid <- 100
  if(is.function(family))
    family <- family()

  createIWLSresults <- function(obj, samples, id = NULL)
  {
    if(inherits(samples[[1]], "mcmc.list"))
      samples <- do.call("c", samples)
    else samples <- as.mcmc.list(list(samples))
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

      ## Smooth terms.
      if(length(obj$smooth)) {
        for(i in 1:length(obj$smooth)) {
          ## Get coefficient samples of smooth term.
          psamples <- NULL
          k <- ncol(obj$smooth[[i]]$X)
          pn <- grep(paste(id, "h1", obj$smooth[[i]]$label, sep = ":"), snames, value = TRUE, fixed = TRUE) ## FIXME: hlevels!
          pn <- pn[!grepl("tau2", pn)]
          psamples <- matrix(samples[[j]][, snames %in% pn], ncol = k)

          ## Possible variance parameter samples.
          vsamples <- NULL
          tau2 <- paste(id, "h1", paste(obj$smooth[[i]]$label, "tau2", sep = "."), sep = ":")
          if(length(tau2 <- grep(tau2, snames, fixed = TRUE))) {
            vsamples <- as.numeric(samples[[j]][, tau2])
          }

          ## Prediction matrix.
          get.X <- function(x) {
            X <- PredictMat(obj$smooth[[i]], x)
            X
          }

          ## Compute final smooth term object.
          tn <- c(obj$smooth[[i]]$term, if(obj$smooth[[i]]$by != "NA") obj$smooth[[i]]$by else NULL)

          if(is.null(obj$smooth[[i]]$is.linear)) {
            if(!is.list(effects))
              effects <- list()
            if(length(effects)) {
              if(obj$smooth[[i]]$label %in% names(effects)) {
                ct <- gsub(".smooth.spec", "", class(obj$smooth[[i]]))[1]
                if(ct == "random.effect") ct <- "re"
                obj$smooth[[i]]$label <- paste(obj$smooth[[i]]$label, ct, sep = ":")
              }
            }

            fst <- compute_term(obj$smooth[[i]], get.X = get.X, get.mu = obj$smooth[[i]]$get.mu,
              psamples = psamples, vsamples = vsamples, FUN = NULL, snames = snames,
              effects.hyp = effects.hyp, fitted.values = fitted.values,
              data = attr(x, "model.frame")[, tn, drop = FALSE], grid = grid)

            attr(fst$term, "specs")$get.mu <- obj$smooth[[i]]$get.mu

            ## Add term to effects list.
            effects[[obj$smooth[[i]]$label]] <- fst$term
            effects.hyp <- fst$effects.hyp

            fitted.values <- fst$fitted.values
            rm(fst)
          } else {
            nx <- colnames(obj$smooth[[i]]$X)
            qu <- t(apply(psamples, 2, quantile, probs = c(0.025, 0.5, 0.975)))
            sd <- drop(apply(psamples, 2, sd))
            me <- drop(apply(psamples, 2, mean))
            param.effects <- cbind(me, sd, qu)
            rownames(param.effects) <- nx
            colnames(param.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
            fitted.values <- as.vector(fitted.values + obj$smooth[[i]]$X %*% param.effects[, 1])
            attr(param.effects, "samples") <- as.mcmc(psamples)
            colnames(attr(param.effects, "samples")) <- nx
          }
        }
      }

      ## Compute partial residuals.
      if(!is.null(effects)) {
        if(length(obj$response)) {
          if(obj$response %in% names(attr(x, "model.frame"))) {
            effects <- partial.residuals(effects, attr(x, "model.frame")[[obj$response]],
              fitted.values, NULL)
          }
        }
      }

      ## Stuff everything together.
      rval[[j]] <- list(
        "model" = list("DIC" = DIC, "pd" = pd, "N" = nrow(attr(x, "model.frame")), "formula" = obj$formula),
        "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "fitted.values" = fitted.values,
        "residuals" = if(is.factor(obj$response.vec)) {
            as.integer(obj$response.vec) - 1
          } else {
            obj$response.vec
          } - fitted.values
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
      rval[[nx[j]]] <- createIWLSresults(x[[nx[j]]], samples, id = fn[j])
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
    return(createIWLSresults(x, samples))
  }
}


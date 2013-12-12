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
          "term" = "parametric",
          "bs.dim" = ncol(obj$X),
          "fixed" = TRUE
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
  attr(x, "model.frame") <- NULL

  x
}


## Function to setup IWLS smooths.
smooth.IWLS <- function(x, ...) {
  UseMethod("smooth.IWLS")
}

smooth.IWLS.default <- function(x, ...)
{
  if(is.null(x$get.mu)) {
    x$get.mu <- if(is.null(x$xt$get.mu)) {
      function(X, g) {
        X %*% as.numeric(g)
      }
    } else x$xt$get.mu
  }

  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b

  if(is.null(x$prop)) {
    x$p.save <- if(x$fixed) "g" else c("g", "tau2")
    x$prop <- list()
    x$prop$g <- runif(ncol(x$X), 0.001, 0.002)
    if(!x$fixed)
      x$prop$tau2 <- runif(1, 0.99, 1)
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("g[", 1:length(x$prop$g), "]", sep = ""), if(!x$fixed) "tau2[1]" else NULL)
    } else x$s.colnames
    if(x$fixed) {
      x$prop$log.priors <- c(
        0,
        dmvnorm(x$prop$g, mean = rep(0.001, ncol(x$X)), sigma = solve(crossprod(x$X)), log = TRUE)
      )
    } else {
      x$prop$log.priors <- c(
        drop(-1 / x$prop$tau2 * crossprod(x$prop$g, x$S[[1]]) %*% x$prop$g),
        dmvnorm(x$prop$g, mean = rep(0.001, ncol(x$X)), sigma = solve(crossprod(x$X)), log = TRUE)
      )
    }
    x$np <- ncol(x$X) + if(x$fixed) 0 else 1
  }

  if(is.null(x$update)) {
    if(is.null(x$rand)) {
      ## Is a list, if length sample > 1 means additional accept check in MCMC step.
      if(x$fixed) {
        x$update <- function(x, z, eta, weights, sample = TRUE, ...) {
          XW <- crossprod(x$X, diag(weights))
          P <- chol2inv(chol(XW %*% x$X))
          M <- P %*% (XW %*% (z - eta))
          if(sample) {
            x$prop$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))
            x$prop$fs <- drop(x$X %*% x$prop$g)
            x$prop$log.priors <- c(
              0,
              dmvnorm(x$prop$g, mean = M, sigma = P, log = TRUE)
            )
            return(x$prop)
          } else return(function(x) { dmvnrom(x, mean = M, sigma = P, log = TRUE) })
        }
      } else {
        x$update <- function(x, z, eta, weights, sample = TRUE, ...) {
          XW <- crossprod(x$X, diag(weights))
          P <- chol2inv(chol(XW %*% x$X + 1 / x$prop$tau2 * x$S[[1]]))
          M <- P %*% (XW %*% (z - eta))
          if(sample) {
            x$prop$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))
            a <- x$rank / 2 + x$a
            b <- 0.5 * crossprod(x$prop$g, x$S[[1]]) %*% x$prop$g + x$b
            x$prop$tau2 <- 1 / rgamma(1, a, b)
            x$prop$fs <- drop(x$X %*% x$prop$g)
            x$prop$log.priors <- c(
              drop(-1 / x$prop$tau2 * crossprod(x$prop$g, x$S[[1]]) %*% x$prop$g),
              dmvnorm(x$prop$g, mean = M, sigma = P, log = TRUE)
            )
            return(x$prop)
          } else return(function(x) { dmvnrom(x, mean = M, sigma = P, log = TRUE) })
        }
      }
    }
  }

  x
}


## Sampler based on IWLS proposals.
samplerIWLS <- function(x, n.iter = 1000, thin = 2, burnin = 200,
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
      if(!is.null(dim(obj$X))) {
        obj$s.alpha <- rep(0, n.save)
        obj$s.samples <- matrix(0, nrow = n.save, ncol = ncol(obj$X))
        obj$fp <- rep(0, nrow(obj$X))
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]]$s.alpha <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.samples <- matrix(0, nrow = n.save, ncol = obj$smooth[[j]]$np)
          obj$smooth[[j]]$samples$fs <- rep(0, nrow(obj$smooth[[j]]$X))
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

  ## IWLS family specs.
  iwls <- family$iwls

  ## Start sampling
  for(i in 1:n.iter) {
    if(save <- i %in% iterthin)
      js <- which(iterthin == i)
      
    for(j in 1:np) {
      ## Sample smooth effects
      if(length(x[[nx[j]]]$smooth)) {
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## PART (1)
          ## Compute weights.
          weights <- iwls[[nx[j]]]$weights(response, eta)

          ## Score.
          score <- iwls[[nx[j]]]$score(response, eta)

          ## Compute working observations.
          z <- eta[[nx[j]]] + 1 / weights * score

          ## Save old predictor.
          eta0 <- eta[[nx[j]]]

          ## Compute old log likelihood.
          p0 <- iwls$loglik(response, eta)

          ## Compute partial predictor.
          eta[[nx[j]]] <- eta0 - x[[nx[j]]]$smooth[[sj]]$samples$fs

          ## Propose new parameters.
          sms <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]], z, eta[[nx[j]]], weights)

          ## Set up new predictor.
          eta[[nx[j]]] <- eta[[nx[j]]] + sms$fs

          ## Compute new log likelihood.
          p1 <- iwls$loglik(response, eta)

          ## PART (2)
          ## Compute new weights
          weights <- iwls[[nx[j]]]$weights(response, eta)

          ## New score.
          score <- iwls[[nx[j]]]$score(response, eta)

          ## New working observations.
          z <- eta[[nx[j]]] + 1 / weights * score

          ## Compute partial predictor.
          eta[[nx[j]]] <- eta0 - sms$fs

          ## Obtain function to evaluate log proposal probablity
          lppfun <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]], z, eta[[nx[j]]], weights, sample = FALSE)

          ## Compute log proposal probablity
          lpp <- lppfun(x[[nx[j]]]$smooth[[sj]]$samples)

          ## Compute acceptance probablity.
          alpha <- (p1 + sms$log.priors[1] + x[[nx[j]]]$smooth[[sj]]$samples$log.priors[2]) -
            (p0 + x[[nx[j]]]$smooth[[sj]]$samples$log.priors[1] + sms$log.priors[2])

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha
print(accepted)
          if(accepted) {
            x[[nx[j]]]$smooth[[sj]]$samples <- sms
          } else eta[[nx[j]]] <- eta0

          ## Save the samples and acceptance.
          if(save) {
            x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- accepted
            x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$samples[x[[nx[j]]]$smooth[[sj]]$p.save])
          }
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
      if(!is.null(dim(obj$X))) {
        fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h", obj$hlevel, ":", "param", ".raw", sep = ""))
        colnames(obj$s.samples) <- paste("b[", 1:ncol(obj$X), "]", sep = "")
        write.table(obj$s.samples, file = fn, row.names = FALSE, quote = FALSE)
      }
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
    colnames(st) <- paste(gsub(".raw", "", j), colnames(st), sep = ":")
    samples <- cbind(samples, st)
  }

  return(as.mcmc(samples))
}


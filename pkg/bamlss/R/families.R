######################
## BAMLSS Families. ##
######################
## Print method.
print.family.bamlss <- function(x, ...)
{
  cat("Family:", x$family, "\n")
  links <- paste(names(x$links), x$links, sep = " = ")
  links <- paste(links, collapse = ", ")
  cat(if(length(links) > 1) "Link functions:" else "Link function:", links, sep = " ")
  cat("\n")
}

## Extract method.
family.bamlss <- function(object, ...)
{
  attr(object, "family")
}

## Second make.link function.
make.link2 <- function(link)
{
  link0 <- link
  mu.eta2 <- function(x) {
    if(link0 == "identity") {
      x$mu.eta2 <- function(eta) rep.int(0, length(eta))
      return(x)
    }
    if(link0 == "log") {
      x$mu.eta2 <- function(eta) exp(eta)
      return(x)
    }
    if(link0 == "logit") {
      x$mu.eta2 <- function(eta) {
        eta <- exp(eta)
        return(-eta * (eta - 1) / (eta + 1)^3)
      }
      return(x)
    }
    if(link0 == "probit") {
      x$mu.eta2 <- function(eta) {
        -eta * dnorm(eta, mean = 0, sd = 1)
      }
      return(x)
    }
    if(link0 == "inverse") {
      x$mu.eta2 <- function(eta) {
        2 / (eta^3)
      }
      return(x)
    }
    if(link0 == "1/mu^2") {
      x$mu.eta2 <- function(eta) {
        0.75 / eta^(2.5)
      }
      return(x)
    }
    if(link0 == "sqrt") {
      x$mu.eta2 <- function(eta) { rep(2, length = length(eta)) }
      return(x)
    }
    stop(paste('higher derivatives of link "', link, '" not available!', sep = ''), call. = FALSE)
  }

  if(link %in% c("logit", "probit", "cauchit", "cloglog", "identity",
                 "log", "sqrt", "1/mu^2", "inverse")) {
    rval <- make.link(link)
  } else {
    rval <- switch(link,
                   "rhogit" = list(
                     "linkfun" = function(mu) { mu / sqrt(1 - mu^2) },
                     "linkinv" = function(eta) { eta / sqrt(1 + eta^2) }, 
                     "mu.eta" = function(eta) { 1 / (1 + eta^2)^1.5 }
                   ),
                   "cloglog2" = list(
                     "linkfun" = function(mu) { log(-log(mu)) },
                     "linkinv" = function(eta) {
                       pmax(pmin(1 - expm1(-exp(eta)), .Machine$double.eps), .Machine$double.eps)
                     },
                     "mu.eta" = function(eta) {
                       eta <- pmin(eta, 700)
                       pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
                     }
                   )
    )
  }
  
  rval <- mu.eta2(rval)
  rval$name <- link
  rval
}

parse.links <- function(links, default.links, ...)
{
  dots <- list(...)
  nl <- names(default.links)
  if(length(dots))
    links <- as.character(dots)
  if(is.null(names(links)))
    names(links) <- rep(nl, length.out = length(links))
  links <- as.list(links)
  for(j in nl) {
    if(is.null(links[[j]]))
      links[[j]] <- default.links[j]
  }
  links <- links[nl]
  links <- as.character(links)
  names(links) <- nl
  links
}


## http://stats.stackexchange.com/questions/41536/how-can-i-model-a-proportion-with-bugs-jags-stan
beta.bamlss <- function(links = c(mu = "logit", sigma2 = "logit"), ...)
{
  rval <- list(
    "family" = "beta",
    "names" = c("mu", "sigma2"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit"), ...),
    "valid.response" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in (0, 1)!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("beta", "mean"),
      "sigma2" = c("beta", "sigma2")
    ),
    "bugs" = list(
      "dist" = "dbeta",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(
        mu = "mu * (1 / sigma2)",
        sigma2 = "(1 - mu) * (1 / sigma2)"
      )
    ),
    "score" = list(
      "mu" = function(y, eta, ...) {
        a <- eta$mu
        b <- eta$sigma2
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(a * hilfs2 * log(y) - a * hilfs2 * log(1 - y) + ((1 - b) / b) * a * (1 - a) * (-digamma(hilfs) + digamma(hilfs2)))
      },
      "sigma2" = function(y, eta, ...) {
        a <- eta$mu
        b <- eta$sigma2
        hilfs <- a*(1-b)/b
        hilfs2 <- (1-a)*(1-b)/b
        drop(-(1 - b) / (b) * ( -a * digamma(hilfs) - (1 - a) * digamma(hilfs2) + digamma((1 - b) / (b)) + a * log(y) + (1 - a) * log(1 - y)))
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        a <- eta$mu
        b <- eta$sigma2
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * a^2 * (1 - a)^2 * (trigamma(hilfs) + trigamma(hilfs2)))
      },
      "sigma2" = function(y, eta, ...) {
        a <- eta$mu
        b <- eta$sigma2
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * (a^2 * trigamma(hilfs) + (1 - a)^2 * trigamma(hilfs2) - trigamma((1 - b) / (b))))
      }
    ),
    "mu" = function(eta, ...) {
      eta$mu
    },
	  "d" = function(y, eta, log = FALSE) {
       mu <- eta$mu
       sigma2 <- eta$sigma2
		   a <- mu * (1 - sigma2) / (sigma2)
		   b <- a * (1 - mu) / mu
		   dbeta(y, shape1 = a, shape2 = b, log = log)
	  },
	  "p" = function(y, eta, ...) {
       mu <- eta$mu
       sigma2 <- eta$sigma2
		   a <- mu * (1 - sigma2) / (sigma2)
		   b <- a * (1 - mu) / mu
		   pbeta(y, shape1 = a, shape2 = b, ...)
	  },
    "type" = 1
  )
  class(rval) <- "family.bamlss"
  rval
}


betazoi.bamlss <- function(links = c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...)
{
  rval <- list(
    "family" = "betazoi",
    "names" = c("mu", "sigma2", "nu", "tau"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0 & x <= 1)) stop("response values not in [0, 1]!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf", "mu"),
      "sigma2" = c("betainf", "sigma2"),
      "nu" = c("betainf", "nu"),
      "tau" = c("betainf", "tau"),
      "order" = 1:4,
      "weights" = list(
        "mu" = function(x) { 1 * ((x != 1) & (x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 1) & (x != 0)) }
      )
    ),
    "mu" = function(eta, ...) {
       eta$mu * (1 - (eta$nu + eta$tau) / (1 + eta$nu + eta$tau)) + eta$tau / (1 + eta$nu + eta$tau)
    },
	  "d" = function(y, eta, log = FALSE) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  a <- mu * (1 - sigma) / (sigma)
		  b <- a * (1 - mu) / mu
		  d <- ifelse(y == 0, eta$nu / (1 + eta$nu + eta$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + eta$nu + eta$tau))
		  ifelse (y==1, eta$tau / (1 + eta$nu + eta$tau), d)
      if(log) d <- log(d)
      d
	  },
	  "p" = function(y, eta, ...) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  a <- mu * (1 - sigma) / (sigma)
		  b <- a * (1 - mu) / mu
		  hilfs <- eta$nu / (1 + eta$nu + eta$tau)
		  hilfs2 <- eta$tau / (1 + eta$nu + eta$tau)
		  cdf <- ifelse(y == 0, hilfs, hilfs + (1 - (hilfs + hilfs2)) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
		  ifelse(y == 1, 1, cdf)
	  }
  )
  class(rval) <- "family.bamlss"
  rval
}


betazi.bamlss <- function(links = c(mu = "logit", sigma2 = "logit", nu = "log"), ...)
{
  rval <- list(
    "family" = "betazi",
    "names" = c("mu", "sigma2", "nu"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", nu = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0 & x < 1)) stop("response values not in [0, 1)!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf0", "mu"),
      "sigma2" = c("betainf0", "sigma2"),
      "nu" = c("betainf0", "nu"),
      "order" = 1:3,
      "weights" = list(
        "mu" = function(x) { 1 * ((x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 0)) }
      )
    ),
    "mu" = function(eta, ...) {
      eta$mu * (1 - (eta$nu) / (1 + eta$nu))
    },
	  "d" = function(y, eta, log = FALSE) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  a <- mu * (1 - sigma) / (sigma)
		  b <- a * (1 - mu) / mu
		  d <- ifelse(y == 0, eta$nu / (1 + eta$nu), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + eta$nu))
      if(log) d <- log(d)
      d
	  },
	  "p" = function(y, eta, ...) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  a <- mu * (1 - sigma) / (sigma)
		  b <- a * (1 - mu) / mu
		  hilfs <- eta$nu / (1 + eta$nu)
		  ifelse(y == 0, hilfs, hilfs + (1 - hilfs) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
	  }
  )
  class(rval) <- "family.bamlss"
  rval
}


betaoi.bamlss <- function(links = c(mu = "logit", sigma2 = "logit", tau = "log"), ...)
{
  rval <- list(
    "family" = "betazi",
    "names" = c("mu", "sigma2", "tau"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0 & x <= 1)) stop("response values not in (0, 1]!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf1", "mu"),
      "sigma2" = c("betainf1", "sigma2"),
      "tau" = c("betainf1", "tau"),
      "order" = 1:3,
      "weights" = list(
        "mu" = function(x) { 1 * ((x != 1)) },
        "sigma2" = function(x) { 1 * ((x != 1)) }
      )
    ),
    "mu" = function(eta, ...) {
      eta$mu * (1 - eta$tau / (1 + eta$tau)) +
        eta$tau / (1 + eta$tau)
    },
    "d" = function(y, eta, log = FALSE) {
      mu <- eta$mu
      sigma <- eta$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 1, eta$tau / (1 + eta$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + eta$tau))
      if(log) d <- log(d)
      d
	  },
    "p" = function(y, eta, ...) {
      mu <- eta$mu
      sigma <- eta$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      hilfs <- eta$tau / (1 + eta$tau)
      ifelse(y == 1, hilfs, hilfs + (1 - hilfs) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
	  }
  )
  class(rval) <- "family.bamlss"
  rval
}


binomial.bamlss <- function(link = "logit", ...)
{
  rval <- list(
    "family" = "binomial",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "logit"), ...),
    "valid.response" = function(x) {
      if(!is.factor(x)) {
        if(length(unique(x)) > 2)
          stop("response has more than 2 levels!", call. = FALSE)
      } else {
        if(nlevels(x) > 2)
          stop("more than 2 levels in factor response!", call. = FALSE)
      }
      TRUE
    },
    "bayesx" = list(
      "pi" = c(paste("binomial", link, sep = "_"), "mu")
    ),
    "bugs" = list(
      "dist" = "dbern",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(eta, ...) {
      eta$pi
    },
	  "d" = function(y, eta, log = FALSE) {
      if(is.factor(y)) y <- as.integer(y) - 1
		  dbinom(y, size = 1, prob = eta$pi, log = log)
	  },
	  "p" = function(y, eta, ...) {
      y <- as.integer(y) - 1
		  pbinom(y, size = 1, prob = eta$pi, ...)
	  },
    "score" = list(
      "pi" = function(y, eta, ...) {
        if(is.factor(y)) y <- as.integer(y) - 1
        (eta$pi - y) / ((eta$pi - 1) * eta$pi)
      }
    ),
    "weights" = list(
      "pi" = function(y, eta, ...) {
        if(is.factor(y)) y <- as.integer(y) - 1
        -1 * ((eta$pi^2 + y - 2 * eta$pi * y) / ((-1 + eta$pi^2) * eta$pi^2))
      }
    ),
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}

cloglog.bamlss <- function(link = "cloglog", ...)
{
  rval <- list(
    "family" = "cloglog",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "cloglog"), ...),
    "valid.response" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      if(nlevels(x) > 2) stop("more than 2 levels in factor response!", call. = FALSE)
      TRUE
    },
    "bayesx" = list(
      "pi" = c(paste("cloglog", link, sep = "_"), "mean")
    ),
    "mu" = function(eta, ...) {
      eta$pi
    },
	  "d" = function(y, eta, log = FALSE) {
		  dbinom(y, size = 1, prob = eta$pi, log = log)
	  },
	  "p" = function(y, eta, ...) {
		  pbinom(y, size = 1, prob = eta$pi, ...)
	  }
  )

  class(rval) <- "family.bamlss"
  rval
}

gaussian.bamlss <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal2", "mu"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal2", "sigma"),
        "logit" = c("normal_sigma_logit", "scale")
      )
    ),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(sigma = "1 / sqrt(sigma)")
    ),
    "score" = list(
      "mu" = function(y, eta, ...) { drop((y - eta$mu) / (eta$sigma^2)) },
      "sigma" = function(y, eta, ...) { drop(-1 + (y - eta$mu)^2 / (eta$sigma^2)) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { drop(1 / (eta$sigma^2)) },
      "sigma" = function(y, eta, ...) { rep(2, length(y)) }
    ),
    "gradient" = list(
      "mu" = function(g, y, eta, x, ...) {
        eta$mu <- eta$mu + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop((y - eta$mu) / (exp(eta$sigma)^2)))
      },
      "sigma" = function(g, y, eta, x, ...) {
        eta$sigma <- eta$sigma + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop(-1 + (y - eta$mu)^2 / (exp(eta$sigma)^2)))
      }
    ),
    "hessian" = list(
      "mu" = function(g, y, eta, x, ...) {
        if(!is.null(x$xbin.ind)) {
          return(crossprod(x$X[x$xbin.ind, , drop = FALSE], x$X[x$xbin.ind, , drop = FALSE] * drop(1 / exp(eta$sigma)^2)))
        } else {
          return(crossprod(x$X, x$X * drop(1 / exp(eta$sigma)^2)))
        }
      },
      "sigma" = function(g, y, eta, x, ...) {
        if(is.null(x$XX)) return(crossprod(x$X)) else return(x$XX)
      }
    ),
    "loglik" = function(y, eta, ...) {
      sum(dnorm(y, eta$mu, eta$sigma, log = TRUE))
    },
    "mu" = function(eta, ...) {
      eta$mu
    },
    "d" = function(y, eta, log = FALSE) {
      dnorm(y, mean = eta$mu, sd = eta$sigma, log = log)
    },
    "p" = function(y, eta, ...) {
      pnorm(y, mean = eta$mu, sd = eta$sigma, ...)
    },
    "type" = 1
  )
  
  class(rval) <- "family.bamlss"
  rval
}


gaussian2.bamlss <- function(links = c(mu = "identity", sigma2 = "log"), ...)
{
  rval <- list(
    "family" = "gaussian2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "bayesx" = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal", "mu"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma2" = switch(links["sigma2"],
        "log" = c("normal", "sigma2"),
        "logit" = c("normal_sigma2_logit", "scale")
      )
    ),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(sigma2 = "1 / sigma2")
    ),
    "score" = list(
      "mu" = function(y, eta, ...) { drop((y - eta$mu) / eta$sigma2) },
      "sigma2" = function(y, eta, ...) { drop(-0.5 + (y - eta$mu)^2 / (2 * eta$sigma2)) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { drop(1 / eta$sigma2) },
      "sigma2" = function(y, eta, ...) { rep(0.5, length(y)) }
    ),
    "loglik" = function(y, eta, ...) {
      sum(dnorm(y, eta$mu, sqrt(eta$sigma2), log = TRUE))
    },
    "mu" = function(eta, ...) {
      eta$mu 
    },
    "d" = function(y, eta, log = FALSE) {
      dnorm(y, mean = eta$mu, sd = sqrt(eta$sigma2), log = log)
    },
    "p" = function(y, eta, ...) {
      pnorm(y, mean = eta$mu, sd = sqrt(eta$sigma2), ...)
    },
    "type" = 1
  )
 
  class(rval) <- "family.bamlss"
  rval
}


truncgaussian2.bamlss <- function(links = c(mu = "identity", sigma2 = "log"), ...)
{
  rval <- list(
    "family" = "truncgaussian2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal", "mu"),
      "sigma2" = c("truncnormal", "sigma2")
    ),
    "mu" = function(eta, ...) {
      mu <-  eta$mu
	    sigma <- sqrt(eta$sigma2)
	    arg <- - mu / sigma
	    mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, eta, log = FALSE) {
	    sigma <- sqrt(eta$sigma2)
	    arg <- - eta$mu / sigma
	    d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
	    sigma <- sqrt(eta$sigma2)
	    arg <- - eta$mu / sigma
	    2 * (pnorm(y / sigma + arg) - pnorm(arg))
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}

truncgaussian.bamlss <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  rval <- list(
    "family" = "truncgaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal2", "mu"),
      "sigma" = c("truncnormal2", "sigma")
    ),
    "mu" = function(eta, ...) {
      mu <-  eta$mu
	    sigma <- eta$sigma
	    arg <- - mu / sigma
	    mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, eta, log = FALSE) {
      mu <-  eta$mu
	    sigma <- eta$sigma
	    arg <- - mu / sigma
	    d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
      mu <-  eta$mu
	    sigma <- eta$sigma
	    arg <- - mu / sigma
	    2 * (pnorm(y / sigma + arg) - pnorm(arg))
    },
    "score" = list(
      "mu" = function(y, eta, ...) {
        rval <- with(eta, (y - mu) / (sigma^2) - (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      },
      "sigma" = function(y, eta, ...) {
        rval <- with(eta, -1 + (y - mu)^2 / (sigma^2) + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        rval <- with(eta, 1 / (sigma^2) - (mu / sigma^2) * (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))
          - ((1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))^2)
        return(drop(rval))
      },
      "sigma" = function(y, eta, ...) {
        rval <- with(eta, 2 - (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)) * 
          (1 + (mu / sigma)^2 + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))))
        return(drop(rval))
      }
    ),
    "loglik" = function(y, eta, ...) {
      rval <- with(eta, sum(-0.5 * log(2 * pi) - log(sigma) - (y - mu)^2 / (2*sigma^2) - log(pnorm(mu / sigma))))
      return(rval)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


trunc.bamlss <- function(links = c(mu = "identity", sigma = "log"),
  direction = "left", point = 0, ...)
{
  tgrad <- function(y, eta, what = "mu") {
    eta$sigma[eta$sigma < 1] <- 1
    eta$mu[eta$mu < -10] <- -10
    eta$mu[eta$mu > 10] <- 10
    resid <- y - eta$mu
    if(direction == "left") {
      trunc <- eta$mu - point
      sgn <- 1
    } else {
      trunc <- point - eta$mu
      sgn <- -1
    }
    mills <- dnorm(trunc / eta$sigma) / pnorm(trunc / eta$sigma)
    g <- if(what == "mu") {
      resid / eta$sigma^2 - sgn / eta$sigma * mills
    } else {
      resid^2 / eta$sigma^3 - 1 / eta$sigma + trunc / eta$sigma^2 * mills
    }
    return(g)
  }

  rval <- list(
    "family" = "truncreg",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "d" = function(y, eta, log = FALSE) {
      eta$sigma[eta$sigma < 1] <- 1
      eta$mu[eta$mu < -10] <- -10
      eta$mu[eta$mu > 10] <- 10
      resid <- y - eta$mu
      if(direction == "left") {
        trunc <- eta$mu - point
      } else {
        trunc <- point - eta$mu
      }
      ll <- log(dnorm(resid / eta$sigma)) - log(eta$sigma) - log(pnorm(trunc / eta$sigma))
      if(!log) ll <- exp(ll)
      ll
    },
    "p" = function(y, eta) {
      require("msm")
      ptnorm(y, mean = eta$mu, sd = eta$sigma,
        lower = if(direction == "left") point else -Inf,
        upper = if(direction == "right") point else Inf)
    },
    "score" = list(
      "mu" = function(y, eta, ...) {
        as.numeric(tgrad(y, eta, what = "mu"))
      },
      "sigma" = function(y, eta, ...) {
        as.numeric(tgrad(y, eta, what = "sigma"))
      }
    )
  )
  
  class(rval) <- "family.bamlss"
  rval
}


cnorm.bamlss <- function(...)
{
  f <- cens.bamlss(left = 0)
  f$transform <- function(x, ...) {
    x <- bamlss.setup(x, ...)
    y <- attr(x, "response.vec")
    check <- as.integer(y <= 0)
    attr(y, "check") <- check
    attr(x, "response.vec") <- y
    x
  }
  f$engine <- function(x, ...) {
    sampler <- function(x, ...) { GMCMC(x, propose = "iwls", ...) }
    stacker(x, optimizer = bfit_cnorm, sampler = sampler, ...)
  }
  f$score <- list(
    "mu" = function(y, eta, ...) {
      .Call("cnorm_score_mu",
        as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, eta, ...) {
      .Call("cnorm_score_sigma",
        as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$weights <- list(
    "mu" = function(y, eta, ...) {
      .Call("cnorm_weights_mu",
        as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, eta, ...) {
      .Call("cnorm_weights_sigma",
        as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$loglik <- function(y, eta, ...) {
    .Call("cnorm_loglik",
      as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma),
      as.integer(attr(y, "check")))
  }
  f$d <- function(y, eta, log = FALSE) {
    ifelse(y <= 0, pnorm(-eta$mu / eta$sigma, log.p = log),
      dnorm((y - eta$mu) / eta$sigma, log = log) / eta$sigma^(1 - log) - log(eta$sigma) * log)
  }
  f$p <- function(y, eta, log = FALSE) {
    ifelse(y < 0, 0, pnorm((y - eta$mu) / eta$sigma, log = log))
  }
  f$q <- function(y, eta, ...) {
    rval <- qnorm(y) * eta$sigma + eta$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$r <- function(n, y, eta) {
    rval <- rnorm(n) * eta$sigma + eta$mu
    pmax(pmin(rval, Inf), 0)
  }
  f
}

pcnorm.bamlss <- function(alpha = NULL, ...)
{
  f <- cens.bamlss(left = 0)
  if(is.null(alpha)) {
    f$names <- c(f$names, "alpha")
    f$links <- c(f$links, "alpha" = "log")
  }
  f$transform <- function(x, ...) {
    x <- bamlss.setup(x, ...)
    y <- attr(x, "response.vec")
    check <- as.integer(y <= 0)
    attr(y, "check") <- check
    attr(x, "response.vec") <- y
    if("alpha" %in% names(x)) {
      if(length(x$alpha$smooth)) {
        for(j in seq_along(x$alpha$smooth))
          x$alpha$smooth[[j]]$propose <- gmcmc_sm.slice
      }
    }
    x
  }
  f$engine <- function(x, ...) {
    optimizer <- if(is.null(alpha)) {
      function(x, ...) { opt0(x, gradient = FALSE, ...) }
    } else {
      function(x, ...) { bfit0(x, ...) }
    }
    sampler <- function(x, ...) { GMCMC(x, propose = "iwls", ...) }
    stacker(x, optimizer = optimizer, sampler = null.sampler, ...)
  }
  f$score <- list(
    "mu" = function(y, eta, ...) {
      if(!is.null(alpha))
        eta$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_score_mu",
        as.numeric(y^(1 / eta$alpha)), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, eta, ...) {
      if(!is.null(alpha))
        eta$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_score_sigma",
        as.numeric(y^(1 / eta$alpha)), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$weights <- list(
    "mu" = function(y, eta, ...) {
      if(!is.null(alpha))
        eta$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_weights_mu",
        as.numeric(y^(1 / eta$alpha)), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, eta, ...) {
      if(!is.null(alpha))
        eta$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_weights_sigma",
        as.numeric(y^(1 / eta$alpha)), as.numeric(eta$mu), as.numeric(eta$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$loglik <- function(y, eta, ...) {
    if(!is.null(alpha))
      eta$alpha <- rep(alpha, length = length(y))
    .Call("cnorm_power_loglik",
      as.numeric(y), as.numeric(eta$mu), as.numeric(eta$sigma), as.numeric(eta$alpha),
      as.integer(attr(y, "check")))
  }
  f$d <- function(y, eta, log = FALSE) {
    if(!is.null(alpha))
      eta$alpha <- rep(alpha, length = length(y))
    ifelse(y <= 0, pnorm(-eta$mu / eta$sigma, log.p = log),
      dnorm((y^(1 / eta$alpha) - eta$mu) / eta$sigma, log = log) / eta$sigma^(1 - log) - log(eta$sigma) * log
      -1 * log(eta$alpha) - (1 / eta$alpha - 1) * log(y))
  }
  f$p <- function(y, eta, log = FALSE) {
    if(!is.null(alpha))
      eta$alpha <- rep(alpha, length = length(y))
    ifelse(y < 0, 0, pnorm((y^(1 / eta$alpha) - eta$mu) / eta$sigma, log = log))
  }
  f$q <- function(y, eta, ...) {
    if(!is.null(alpha))
      eta$alpha <- rep(alpha, length = length(y))
    rval <- qnorm(y^(1 / eta$alpha)) * eta$sigma + eta$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$r <- function(n, y, eta) {
    rval <- rnorm(n) * eta$sigma + eta$mu
    pmax(pmin(rval, Inf), 0)
  }
  f
}


cens.bamlss <- function(links = c(mu = "identity", sigma = "log", df = "log"),
  left = 0, right = Inf, dist = "gaussian", ...)
{
  dist <- match.arg(dist, c("student", "gaussian", "logistic"))

  ddist <- switch(dist,
    "student"  = function(x, location, scale, df, log = TRUE) 
      dt((x - location)/scale, df = df, log = log)/scale^(1-log) - 
      log*log(scale),
    "gaussian" = function(x, location, scale, df, log = TRUE) 
      dnorm((x - location)/scale, log = log)/scale^(1-log) - 
      log*log(scale),
    "logistic" = function(x, location, scale, df, log = TRUE) 
      dlogis((x - location)/scale, log = log)/scale^(1-log) - 
      log*log(scale)
  )
  pdist <- switch(dist,
    "student"  = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) pt((x - location)/scale, df = df, lower.tail = lower.tail,
      log.p = log.p),
    "gaussian" = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) pnorm((x - location)/scale, lower.tail = lower.tail, 
      log.p = log.p),
    "logistic" = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) plogis((x - location)/scale, lower.tail = lower.tail, 
      log.p = log.p)
  )

  gradfun <- function(y, eta, type = "gradient", name = "mu") {
    ## functions used to evaluate gradient and hessian
    mills <- function(y, lower.tail = TRUE) {
      with(eta, sigma * ddist(y, mu, sigma, df, log = FALSE)/
      pdist(y, mu, sigma, df, log.p = FALSE, lower.tail = lower.tail))
    }

    ## ddensity/dmu
    d1 <- with(eta, switch(dist,
      "student"  = function(x)  
        (x - mu)/sigma^2 * (df + 1) / (df + (x - mu)^2/sigma^2),
      "gaussian" = function(x) 
        (x - mu)/sigma^2,
      "logistic" = function(x)  
        (1 - 2 * pdist(-x, - mu, sigma, log.p = FALSE))/sigma)
    )
    
    ## ddensity/dsigma
    d2 <- function(x) with(eta, d1(x) * (x-mu))

    ## d^2density/dmu^2
    d3 <- with(eta, switch(dist,
      "student"  = function(x)  
        (df + 1)*((x - mu)^2 - df*sigma^2) / (df*sigma^2 + (x - mu)^2)^2,
      "gaussian" = function(x) 
        - 1/sigma^2,
      "logistic" = function(x)  
        - 2/sigma * ddist(x, mu, sigma, log = FALSE))
    )    
    
    ## d^2density/dsigma^2
    d5 <- with(eta, switch(dist,
      "student"  = function(x)  
        - (x - mu)^2 * (df + 1) / (df*sigma^2 + (x - mu)^2)^2*2*df*sigma^2,
      "gaussian" = function(x) 
        2 * d3(x) * (x-mu)^2,
      "logistic" = function(x)  
        - d2(x) - 2*(x-mu)^2/sigma*ddist(x,mu,sigma, log = FALSE)
    ))
      
    ## d^2density/dmudsigma
    d4 <- with(eta, switch(dist,
      "student"  = function(x)  
          d5(x) / (x - mu),
      "gaussian" = function(x) 
        2 * d3(x) * (x-mu),
      "logistic" = function(x)  
        - d1(x) + (x-mu)*d3(x)
    ))

    ## compute gradient
    if(type == "gradient") {
      if(name == "mu") {
        rval <- with(eta, ifelse(y <= left, 
          - mills(left)/sigma,
          ifelse(y >= right, 
            mills(right, lower.tail = FALSE)/sigma,
            d1(y)
          )))
      } else {
        rval <- with(eta, ifelse(y <= left, 
          - mills(left) * (left - mu)/sigma,
          ifelse(y >= right, 
            mills(right, lower.tail = FALSE) * (right - mu)/sigma,
            d2(y) - 1
          )))
      }

    ## compute hessian
    } else {
      if(name == "mu") {
        rval <- with(eta, ifelse(y <= left, 
          -d1(left)/sigma * mills(left) - mills(left)^2/sigma^2,
          ifelse(y >= right, 
            d1(right)/sigma * mills(right, lower.tail = FALSE) - 
              mills(right, lower.tail = FALSE)^2/sigma^2,
            d3(y)
          )))
      } else {
        rval <- with(eta, ifelse(y <= left, 
          ((left-mu)/sigma - (left-mu)*d2(left))*mills(left) - 
            (left - mu)^2/sigma^2 * mills(left)^2,
          ifelse(y >= right, 
            (-(right-mu)/sigma + (right-mu)*d2(right))*
              mills(right, lower.tail = FALSE) 
              - (right - mu)^2/sigma^2 * mills(right, lower.tail = FALSE)^2,
            d5(y)
          )))
      }
    }
    return(rval)
  }

  score <- list(
    "mu" =  function(y, eta, ...) {
      gradmu <- gradfun(y, eta, type = "gradient", name = "mu")
      return(drop(gradmu))
    },
    "sigma" =  function(y, eta, ...) {
      gradsigma <- gradfun(y, eta, type = "gradient", name = "sigma")
      return(drop(gradsigma))
    }
  )

  weights <- list(
    "mu" =  function(y, eta, ...) {
      wmu <- -1 * gradfun(y, eta, type = "weights", name = "mu")
      return(drop(wmu))
    },
    "sigma" =  function(y, eta, ...) {
      wsigma <- -1 * gradfun(y, eta, type = "weights", name = "sigma")
      return(drop(wsigma))
    }
  )

  names <- switch(dist,
    "student" = c("mu", "sigma", "df"),
    "gaussian" = c("mu", "sigma"),
    "logistic" = c("mu", "sigma")
  )
  
  i <- 1:length(names)

  rval <- list(
    "family" = "cens",
    "names" = names,
    "links" = parse.links(links[i], c(mu = "identity", sigma = "log", df = "log")[i], ...),
    "d" = function(y, eta, log = FALSE, ...) {
      ll <- with(eta, ifelse(y <= left,
        pdist(left, mu, sigma, df, lower.tail = TRUE, log = TRUE),
        ifelse(y >= right,
          pdist(right, mu, sigma, df, lower.tail = FALSE, log = TRUE),
          ddist(y, mu, sigma, df, log = TRUE))))
      if(!log) ll <- exp(ll)
      return(ll)
    },
    "p" = function(y, eta, log = FALSE, ...) {
      with(eta, pdist(y, mu, sigma, df, lower.tail = TRUE, log = log))
    },
    "score" = score,
    "weights" = weights,
    "type" = 1
  )
 
  class(rval) <- "family.bamlss"
  rval
}


cens0.bamlss <- function(links = c(mu = "identity", sigma = "log", df = "log"),
  left = 0, right = Inf, dist = "gaussian", ...)
{
  dist <- match.arg(dist, c("student", "gaussian", "logistic"))

  ddist <- switch(dist,
    "student"  = function(x, location, scale, df, log = TRUE) 
      dt((x - location)/scale, df = df, log = log)/scale^(1-log) - 
      log*log(scale),
    "gaussian" = function(x, location, scale, df, log = TRUE) 
      dnorm((x - location)/scale, log = log)/scale^(1-log) - 
      log*log(scale),
    "logistic" = function(x, location, scale, df, log = TRUE) 
      dlogis((x - location)/scale, log = log)/scale^(1-log) - 
      log*log(scale)
  )
  pdist <- switch(dist,
    "student"  = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) pt((x - location)/scale, df = df, lower.tail = lower.tail,
      log.p = log.p),
    "gaussian" = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) pnorm((x - location)/scale, lower.tail = lower.tail, 
      log.p = log.p),
    "logistic" = function(x, location, scale, df, lower.tail = TRUE, 
      log.p = TRUE) plogis((x - location)/scale, lower.tail = lower.tail, 
      log.p = log.p)
  )

  dddist <- switch(dist,
    "student"  = function(x, location, scale, df) 
      - ddist(x, location, scale, df, log = FALSE) * 
      (x - location)/scale^2 * (df + 1) / (df + (x - location)^2/scale^2),
    "gaussian" = function(x, location, scale, df) 
      - (x - location) * ddist(x, location, scale, log = FALSE)/scale^2,
    "logistic" = function(x, location, scale, df) 
      ddist(x, location, scale, df, log = FALSE)/scale * 
      (- 1 + 2 * pdist(-x, - location, scale, log.p = FALSE))
  )

  if(dist == "student") {
    score <- NULL
  } else {
    score <- list(
      "mu" =  function(y, eta) {
         gradmu <- with(eta, ifelse(y <= left, 
          - ddist(left, mu, sigma, df, log = FALSE) /
            pdist(left, mu, sigma, df, log.p = FALSE),
          ifelse(y >= right, 
          ddist(right, mu, sigma, df, log = FALSE) /
            pdist(right, mu, sigma, df, lower.tail = FALSE, log.p = FALSE),
          - dddist(y, mu, sigma, df)/ddist(y, mu, sigma, df, log = FALSE))))
        return(drop(gradmu))
      },
      "sigma" =  function(y, eta) {
        gradsigma <- with(eta, ifelse(y <= left, 
          - ddist(left, mu, sigma, df, log = FALSE) /
            pdist(left, mu, sigma, df, log.p = FALSE) * (left - mu),
          ifelse(y >= right, 
          ddist(right, mu, sigma, df, log = FALSE)/
            pdist(right, mu, sigma, df, lower.tail = FALSE, log.p = FALSE)*
            (right - mu),
          - dddist(y, mu, sigma, df) * (y - mu)/
            ddist(y, mu, sigma, df, log = FALSE) - 1)))
         return(drop(gradsigma))
      }
    )
  }

  names <- switch(dist,
    "student" = c("mu", "sigma", "df"),
    "gaussian" = c("mu", "sigma"),
    "logistic" = c("mu", "sigma")
  )
  
  i <- 1:length(names)

  rval <- list(
    "family" = "cens",
    "names" = names,
    "links" = parse.links(links[i], c(mu = "identity", sigma = "log", df = "log")[i], ...),
    "d" = function(y, eta, log = FALSE, ...) {
      ll <- with(eta, ifelse(y <= left,
        pdist(left, mu, sigma, df, lower.tail = TRUE, log = TRUE),
        ifelse(y >= right,
          pdist(right, mu, sigma, df, lower.tail = FALSE, log = TRUE),
          ddist(y, mu, sigma, df, log = TRUE))))
      if(!log) ll <- exp(ll)
      return(ll)
    },
    "score" = score,
    "type" = 1
  )
 
  class(rval) <- "family.bamlss"
  rval
}


## Truncated distributions.
dtrunc <- function(x, spec, a = 1, b = Inf, ...) {
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...) / (G(b, ...) - G(a, ...))
  return(tt)
}

ptrunc <- function(x, spec, a = -Inf, b = Inf, ...)
{
  tt <- x
  aa <- rep(a, length(x))
  bb <- rep(b, length(x))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
  tt <- tt - G(aa, ...)
  tt <- tt / (G(bb, ...) - G(aa, ...))
  return(tt)
}

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p * (G(b, ...) - G(a, ...)), ...)
  return(tt)
}

trunc2.bamlss <- function(links = c(mu = "identity", sigma = "log"),
  name = "norm", a = -Inf, b = Inf, ...)
{
  rval <- list(
    "family" = "trunc2",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "d" = function(y, eta, log = FALSE) {
      if(name == "gamma") {
		    a2 <- eta$sigma
		    s2 <- eta$mu / eta$sigma
        d <- dtrunc(y, name, a = a, b = b, shape = s2, scale = s2, log = log)
      } else {
        d <- dtrunc(y, name, a = a, b = b, mean = eta$mu, sd = eta$sigma, log = log)
      }
      return(d)
    },
    "p" = function(y, eta, ...) {
      if(name == "gamma") {
		    a2 <- eta$sigma
		    s2 <- eta$mu / eta$sigma
        p <- ptrunc(y, name, a = a, b = b, shape = a2, scale = s2, ...)
      } else {
        p <- ptrunc(y, name, a = a, b = b, mean = eta$mu, sd = eta$sigma, ...)
      }
      return(p)
    },
    "q" = function(y, eta, ...) {
      q <- qtrunc(y, name, a = a, b = b, mean = eta$mu, sd = eta$sigma, ...)
      return(q)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


t.bamlss <- function(links = c(mu = "identity", sigma2 = "log", df = "log"), ...)
{
  rval <- list(
    "family" = "t",
    "names" = c("mu", "sigma2", "df"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log", df = "log"), ...),
    "bayesx" = list(
      "mu" = c("t", "mu"),
      "sigma2" =  c("t", "sigma2"),
	    "df" = c("t", "df")
    ),
    "bugs" = list(
      "dist" = "dt",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(eta, ...) {
      rval <- eta$mu
      rval[eta$df <= 1] <- 0
      rval
    },
    "d" = function(y, eta, log = FALSE) {
      arg <- (y - eta$mu) / sqrt(eta$sigma2)
      dt(arg, df = eta$df, log = log)
    },
    "p" = function(y, eta, ...) {
      arg <- (y - eta$mu) / sqrt(eta$sigma2)
      pt(arg, df = eta$df, ...)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


invgaussian.bamlss <- function(links = c(mu = "log", sigma2 = "log"), ...)
{
  rval <- list(
    "family" = "invgaussian",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "log", sigma2 = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu"  = c("invgaussian", "mu"),
      "sigma2" = c("invgaussian", "sigma2")
    ),
    "score" = list(
      "mu" = function(y, eta, ...) {
        mu <- eta$mu
        (y - mu) / (mu^2 * eta$sigma2)
      },
      "sigma2" = function(y, eta, ...) {
        mu <- eta$mu
        -0.5 + (y - mu)^2 / (2 * y * mu^2 * eta$sigma2)
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { 1 / (eta$mu * eta$sigma2) },
      "sigma2" = function(y, eta) { rep(0.5, length(y)) }
    ),
    "mu" = function(eta, ...) {
      eta$mu
    },
	  "d" = function (y, eta, log = FALSE) {
		  mu <- eta$mu
		  sigma <- sqrt(eta$sigma2)
		  d <- exp( -0.5 * log(2 * pi) - log(sigma) - (3 / 2) * log(y) - ((y - mu)^2) / (2 * sigma^2 * (mu^2) * y))
      if(log) d <- log(d)
      d
	  },
	  "p" = function (y, eta, ...) {
		  mu <- eta$mu
		  lambda <- 1 / sqrt(eta$sigma2)
		  lq <- sqrt(lambda / y)
		  qm <- y / mu
		  pnorm(lq * (qm - 1)) + exp(2 * lambda / mu) * pnorm(-lq * (qm + 1), ...)
	  }
  )

  class(rval) <- "family.bamlss"
  rval
}

weibull.bamlss <- function(links = c(lambda = "log", alpha = "log"), ...)
{
  rval <- list(
    "family" = "weibull",
    "names" = c("lambda", "alpha"),
    "links" = parse.links(links, c(lambda = "log", alpha = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "lambda"  = c("weibull", "lambda"),
      "alpha" = c("weibull", "alpha")
    ),
    "bugs" = list(
      "dist" = "dweib",
      "eta" = BUGSeta,
      "model" = BUGSmodel,  ## FIXME: order of parameters?
      "order" = 2:1
    ),
    "mu" = function(eta, ...) {
      alpha <-  eta$alpha
      lambda <- eta$lambda
      alpha * gamma(1 + 1 / lambda)
    },
    "d" = function (y, eta, log = FALSE) {
      alpha <-  eta$alpha
      lambda <- eta$lambda
      dweibull(y, scale = lambda, shape = alpha, log = log)
    },
    "p" = function (y, eta, ...) {
      alpha <- eta$alpha
      lambda <- eta$lambda
      rweibull(y, scale = lambda, shape = alpha, ...)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


pareto.bamlss <- function(links = c(b = "log", p = "log"), ...)
{
  rval <- list(
    "family" = "pareto",
    "names" = c("b", "p"),
    "links" = parse.links(links, c(b = "log", p = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "b"  = c("pareto", "b"),
      "p" = c("pareto", "p")
    ),
    "bugs" = list(
      "dist" = "dpar",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(eta, ...) {
      p <- eta$p
      b <- eta$b
      p * gamma(1 + 1 / b)
    },
    "d" = function (y, eta, log = FALSE) {
      p <- eta$p
      b <- eta$b
      d <- p * b^p * (y + b)^(-p - 1)
      if(log) d <- log(d)
      d
    },
    "p" = function (y, eta, ...) {
      p <- eta$p
      b <- eta$b
      1 - ((b/(b + y))^(p))
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


gamma.bamlss <- function(links = c(mu = "log", sigma = "log"), ...)
{
  rval <- list(
    "family" = "gamma",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "log", sigma = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("gamma", "mu"),
      "sigma" = c("gamma", "sigma")
    ),
    "bugs" = list(
      "dist" = "dgamma",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "score" = list(
      "mu" = function(y, eta, ...) {
        sigma <- eta$sigma
        sigma * (-1 + y / eta$mu)
      },
      "sigma" = function(y, eta, ...) {
        mu <- eta$mu
        sigma <- eta$sigma
        sigma * (log(sigma) + 1 - log(mu) - digamma(sigma) + log(y) - y / mu)
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { eta$sigma },
      "sigma" = function(y, eta, ...) {
        sigma <- eta$sigma
        sigma^2 * trigamma(sigma) - sigma
      }
    ),
    "loglik" = function(y, eta, ...) {
		  a <- eta$sigma
		  s <- eta$mu / eta$sigma 
		  sum(dgamma(y, shape = a, scale = s, log = TRUE), na.rm = TRUE)
    },
    "mu" = function(eta, ...) {
      eta$mu
    },
	  "d" = function(y, eta, log = FALSE) {
		  a <- eta$sigma
		  s <- eta$mu / eta$sigma
		  dgamma(y, shape = a, scale = s, log = log)
	  },
	  "p" = function(y, eta, lower.tail = TRUE, log.p = FALSE) {
		  a <- eta$sigma
		  s <- eta$mu / eta$sigma
		  pgamma(y, shape = a, scale = s, lower.tail = lower.tail, log.p = log.p)
	  },
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}


gengamma.bamlss <- function(links = c(mu = "log", sigma = "log", tau = "log"), ...)
{
  rval <- list(
    "family" = "gengamma",
    "names" = c("mu", "sigma", "tau"),
    "links" = parse.links(links, c(mu = "log", sigma = "log", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("gengamma", "mu"),
      "sigma" = c("gengamma", "sigma"),
	    "tau" = c("gengamma", "tau")
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


lognormal.bamlss <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  rval <- list(
    "family" = "lognormal",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("lognormal2", "mu"),
      "sigma" = c("lognormal2", "sigma")
    ),
    "score" = list(
      "mu" = function(y, eta, ...) { (log(y) - eta$mu) / (eta$sigma^2) },
      "sigma" = function(y, eta, ...) { -1 + (log(y) - eta$mu)^2 / (eta$sigma^2) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { 1 / (eta$sigma^2) },
      "sigma" = function(y, eta, ...) { rep(2, length(y)) }
    ),
	  "mu" = function(eta, ...) {
      exp(eta$mu + 0.5 * (eta$sigma)^2)
    },
    "d" = function(y, eta, log = FALSE) {
      dlnorm(y, meanlog = eta$mu, sdlog = eta$sigma, log = log)
    },
    "p" = function(y, eta, ...) {
      plnorm(y, meanlog = eta$mu, sdlog = eta$sigma, ...)
    },
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}


lognormal2.bamlss <- function(links = c(mu = "identity", sigma2 = "log"), ...)
{
  rval <- list(
    "family" = "lognormal2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("lognormal", "mu"),
      "sigma2" = c("lognormal", "sigma2")
    ),
	  "score" = list(
      "mu" = function(y, eta, ...) { (log(y) - eta$mu) / (eta$sigma2) },
      "sigma2" = function(y, eta, ...) { -0.5 + (log(y) - eta$mu)^2 / (2 * eta$sigma2) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { 1 / (eta$sigma2) },
      "sigma2" = function(y, eta, ...) { rep(0.5, length(y)) }
    ),
	  "mu" = function(eta, ...) {
      exp(eta$mu + 0.5 * (eta$sigma2))
    },
    "d" = function(y, eta, log = FALSE) {
      dlnorm(y, meanlog = eta$mu, sdlog = sqrt(eta$sigma2), log = log)
    },
    "p" = function(y, eta, ...) {
      plnorm(y, meanlog = eta$mu, sdlog = sqrt(eta$sigma2), ...)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


dagum.bamlss <- function(links = c(a = "log", b = "log", p = "log"), ...)
{
  rval <- list(
    "family" = "dagum",
    "names" = c("a", "b", "p"),
    "links" = parse.links(links, c(a = "log", b = "log", p = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "a" = c("dagum", "a"),
      "b" = c("dagum", "b"),
      "p" = c("dagum", "p")
    ),
	  "mu" = function(eta, ...) {
	    a <- eta$a
      b <- eta$b
      p <- eta$p
      -(b/a) * (gamma(- 1 / a) * gamma(p + 1 / a)) / (gamma(p))
    },
	  "d" = function(y, eta, log = FALSE) {
		  a <- eta$a
		  b <- eta$b
		  p <- eta$p
		  ap <- a * p
		  d <-ap * y^(ap -1) / (b^ap * (1 + (y / b)^a)^(p + 1))
      if(log) d <- log(d)
      d
	  },
	  "p" = function(y, eta, ...) {
		  a <- eta$a
      b <- eta$b
      p <- eta$p
		  (1 + (y / b)^(-a))^(-p)
	  }
  )

  class(rval) <- "family.bamlss"
  rval
}


BCCG2.bamlss <- function(links = c(mu = "log", sigma = "log", nu = "identity"), ...)
{
  rval <- list(
    "family" = "BCCG",
    "names" = c("mu", "sigma", "nu"),
    "links" = parse.links(links, c(mu = "log", sigma = "log", nu = "identity"), ...),
	  "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("BCCG", "mu"),
      "sigma" =  c("BCCG", "sigma"),
	    "nu" = c("BCCG", "nu"),
	    "order" = c("nu", "sigma", "mu")
    ),
	  "d" = function(y, eta, log = FALSE) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  nu <- eta$nu
		  z <- ifelse(nu == 0, log(y/mu)/sigma, (((y/mu)^nu - 1)/(nu * sigma)))
		  d <- (1 / (sqrt(2 * pi) * sigma)) * (y^(nu - 1) / mu^nu) * exp(-z^2 / 2)
      if(log) d <- log(d)
      d
	  },
	  "p" = function(y, eta, ...) {
		  mu <- eta$mu
		  sigma <- eta$sigma
		  nu <- eta$nu
		  z <- ifelse(nu == 0, log(y/mu)/sigma, (((y/mu)^nu - 1)/(nu * sigma)))
		  FYy1 <- pnorm(z, ...)
      FYy2 <- ifelse(nu > 0, pnorm(-1/(sigma * abs(nu))), 0)
		  FYy3 <- pnorm(1/(sigma * abs(nu)), ...)
      (FYy1 - FYy2)/FYy3
	  },
	  "type" = 1
  )
  
  class(rval) <- "family.bamlss"
  rval
}


mvn.bamlss <- function(links = c(mu1 = "identity", mu2 = "identity",
  sigma1 = "log", sigma2 = "log", rho = "rhogit"), ...)
{
  rval <- list(
    "family" = "mvn",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity",
      sigma1 = "log", sigma2 = "log", rho = "rhogit"), ...),
    "bayesx" = list(
      "mu1" = c("bivnormal_mu", "mean"),
      "mu2" = c("bivnormal_mu", "mu"),
      "sigma1" = c("bivnormal_sigma", "scale1"),
      "sigma2" = c("bivnormal_sigma", "scale2"),
      "rho" = c("bivnormal_rho", "rho"),
      "order" = 5:1,
      "rm.number" = TRUE
    ),
	  "mu" = function(eta, ...) {
      cbind(eta$mu1, eta$mu2)
    },
    "d" = function(y, eta, log = FALSE) {
      cbind(
        dnorm(y[, 1], mean = eta$mu1, sd = eta$sigma1, log = log),
        dnorm(y[, 2], mean = eta$mu2, sd = eta$sigma2, log = log)
      )
    },
    "p" = function(y, eta, ...) {
      cbind(
        pnorm(y[, 1], mean = eta$mu1, sd = eta$sigma1, ...),
        pnorm(y[, 2], mean = eta$mu2, sd = eta$sigma2, ...)
      )
    },
    "type" = 2
  )

  class(rval) <- "family.bamlss"
  rval
}


bivprobit.bamlss <- function(links = c(mu1 = "identity", mu2 = "identity", rho = "rhogit"), ...)
{
  rval <- list(
    "family" = "bivprobit",
    "names" = c("mu1", "mu2", "rho"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity", rho = "rhogit"), ...),
    "bayesx" = list(
      "mu1" = c("bivprobit", "mu"),
      "mu2" = c("bivprobit", "mu"),
      "rho" = c("bivprobit", "rho"),
      "order" = 3:1,
      "rm.number" = TRUE
    ),
	  "mu" = function(eta, ...) {
      c(eta$mu1, eta$mu2)
    },
    "type" = 2
  )

  class(rval) <- "family.bamlss"
  rval
}


bivlogit.bamlss <- function(links = c(p1 = "logit", p2 = "logit", psi = "log"), ...)
{
  rval <- list(
    "family" = "bivlogit",
    "names" = c("p1", "p2", "psi"),
    "links" = parse.links(links, c(p1 = "logit", p2 = "logit", psi = "log"), ...),
    "bayesx" = list(
      "mu1" = c("bivlogit", "mu"),
      "mu2" = c("bivlogit", "mu"),
      "psi" = c("bivlogit", "oddsratio"),
      "order" = 3:1,
      "rm.number" = TRUE
    ),
	  "mu" = function(eta, ...) {
      c(eta$mu1, eta$mu2)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


mvt.bamlss <- function(links = c(mu1 = "identity", mu2 = "identity",
  sigma1 = "log", sigma2 = "log", rho = "fisherz", df = "log"), ...)
{
  rval <- list(
    "family" = "mvt",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho", "df"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity",
      sigma1 = "log", sigma2 = "log", rho = "fisherz", df = "log"), ...),
    "bayesx" = list(
      "mu1" = c("bivt", "mu"),
      "mu2" = c("bivt", "mu"),
      "sigma1" = c("bivt", "sigma"),
      "sigma2" = c("bivt", "sigma"),
      "rho" = c("bivt", "rho"),
      "df" = c("bivt", "df"),
      "order" = 6:1,
      "rm.number" = TRUE
    ),
	  "mu" = function(eta, ...) {
      c(eta$mu1, eta$mu2)
    },
    "type" = 2
  )
  class(rval) <- "family.bamlss"
  rval
}


dirichlet.bamlss <- function(link = "logit", ...)
{
  rval <- list(
    "family" = "dirichlet",
    "names" = "alpha",
    "links" = parse.links(link, c(pi = "logit"), ...),
    "bayesx" = list(
      "alpha" = c("dirichlet", "mu", "alpha")
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


multinomial.bamlss <- function(link = "logit", ...)
{
  rval <- list(
    "family" = "multinomial",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "probit"), ...),
    "cat" = TRUE,
    "bugs" = list(
      "dist" = "dcat",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "bayesx" = list(
      "pi" = c(paste("multinom", link, sep = "_"), "mu", "meanservant")
    )
  )

  rval$d <- switch(link,
    "logit" = function(y, eta, log = FALSE) {
      eta <- as.matrix(as.data.frame(eta))
      eta <- cbind(eta, exp(0))
      eta <- eta / rowSums(eta)
      d <- dcat(y, eta, log = log)
      return(d)
    },
    "probit" = function(...) stop("multinomial probit not supported yet!")
  )

  class(rval) <- "family.bamlss"
  rval
}


## Count Data distributions
poisson.bamlss <- function(links = c(lambda = "log"), ...)
{
  rval <- list(
    "family" = "poisson",
    "names" = c("lambda"),
    "links" = parse.links(links, c(lambda = "log"), ...),
    "bayesx" = list(
      "lambda" = c("poisson", "lambda")
    ),
    "bugs" = list(
      "dist" = "dpois",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
	  "mu" = function(eta, ...) {
       eta$lambda
    },
	  "d" = function(y, eta, log = FALSE) {
      dpois(y, lambda = eta$lambda, log = log)
    },
    "p" = function(y, eta, ...) {
      ppois(y, lambda = eta$lambda, ...)
    },
	  "score" = list(
      "lambda" = function(y, eta, ...) {
        y / eta$lambda - 1
      }
    ),
	  "weights" = list(
      "lambda" = function(y, eta, ...) {
        1 / eta$lambda
      }
    ),
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


zip.bamlss <- function(links = c(lambda = "log", pi = "logit"), ...)
{
  rval <- list(
    "family" = "zip",
    "names" = c("lambda", "pi"),
    "links" = parse.links(links, c(lambda = "log", pi = "logit"), ...),
    "bayesx" = list(
      "lambda" = c("zip", "lambda"),
      "pi" = switch(links["pi"],
        "logit" = c("zip", "pi"),
        "cloglog2" = c("zip", "pi")
      ) 
    ),
	  "mu" = function(eta, ...) {
      eta$lambda * (1 - eta$pi)
    },
    "d" = function(y, eta, log = FALSE) {
      d <- ifelse(y == 0, eta$pi + (1 - eta$pi) * dpois(y, lambda = eta$lambda), 
				(1 - eta$pi) * dpois(y, lambda = eta$lambda))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
      ifelse(y < 0, 0, eta$pi + (1 - eta$pi) * ppois(y, lambda = eta$lambda))
    },
    "nscore" = TRUE,
    "type" = 3
  )
  if(rval$bayesx[[2]][[1]] == "zip_pi_cloglog")
    rval$bayesx[[1]][[1]] <- "zip_lambda_cloglog"

  class(rval) <- "family.bamlss"
  rval
}


hurdleP.bamlss <- function(links = c(lambda = "log", pi = "logit"), ...)
{ 
  rval <- list(
    "family" = "hurdle",
    "names" = c("lambda", "pi"),
    "links" = parse.links(links, c(lambda = "log", pi = "logit"), ...),
    "bayesx" = list(
      "lambda" = c("hurdle", "lambda"),
      "pi" = c("hurdle", "pi"),
	    "weights" = list(
        "lambda" = function(x) { 1 * (x != 0)}
      )
    ),
	  "mu" = function(eta, ...) {
      (1 - eta$pi) * eta$lambda / (1 - exp(-eta$lambda))
    },
    "d" = function(y, eta, log = FALSE) {
      d <- ifelse(y == 0, eta$pi, 
				(1 - eta$pi) * dpois(y, lambda = eta$lambda) / (1 - exp(-eta$lambda)))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
		  cdf1 <- ppois(y, lambda = eta$lambda)
		  cdf2 <- ppois(0, lambda = eta$lambda)
		  cdf3 <- eta$pi + ((1 - eta$pi) * (cdf1 - cdf2)/(1 - cdf2))
		  cdf <- ifelse((y == 0), eta$pi, cdf3)
      cdf
    },
    "nscore" = TRUE,
    "type" = 3
  )
 
  class(rval) <- "family.bamlss"
  rval
}

negbin.bamlss <- function(links = c(mu = "log", delta = "log"), ...)
{
  rval <- list(
    "family" = "negbin",
    "names" = c("mu", "delta"),
    "links" = parse.links(links, c(mu = "log", delta = "log"), ...),
    "bayesx" = list(
      "mu" = c("negbin", "mu"),
      "delta" = c("negbin", "delta")
    ),
	  "mu" = function(eta, ...) {
      eta$mu
    },
    "d" = function(y, eta, log = FALSE) {
      dnbinom(y, mu = eta$mu, size = eta$delta, log = log)
    },
    "p" = function(y, eta, ...) {
      pnbinom(y, mu = eta$mu, size = eta$delta)
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


zinb.bamlss <- function(links = c(mu = "log", pi = "logit", delta = "log"), ...)
{
  rval <- list(
    "family" = "zinb",
    "names" = c("mu", "pi", "delta"),
    "links" = parse.links(links, c(mu = "log", pi = "logit", delta = "log"), ...),
    "bayesx" = list(
      "mu" = c("zinb", "mu"),
      "pi" = c("zinb", "pi"),
      "delta" = c("zinb", "delta")
    ),
	  "mu" = function(eta, ...) {
      eta$mu * (1 - eta$pi)
    },
    "d" = function(y, eta, log = FALSE) {
      d <- ifelse(y == 0, eta$pi + (1 - eta$pi) * dnbinom(y, mu = eta$mu, size = eta$delta), 
				(1 - eta$pi) * dnbinom(y, mu = eta$mu, size = eta$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
      ifelse(y<0, 0, eta$pi + (1 - eta$pi) * pnbinom(y, size = eta$delta, mu = eta$mu))
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}

hurdleNB.bamlss <- function(links = c(mu = "log", pi = "logit", delta = "log"), ...)
{
  rval <- list(
    "family" = "hurdleNB",
    "names" = c("mu", "delta", "pi"),
    "links" = parse.links(links, c(mu = "log", delta = "log", pi = "logit"), ...),
    "bayesx" = list(
      "mu" = c("hurdle", "mu"),
      "delta" = c("hurdle", "delta"),
	    "pi" = c("hurdle", "pi"),
	     "weights" = list(
         "mu" = function(x) { 1 * (x != 0)},
		     "delta" = function(x) { 1 * (x != 0)}
       )
    ),
	  "mu" = function(eta, ...) {
      (1 - eta$pi) * eta$mu / (1 - (eta$delta) / (eta$delta + eta$mu)^eta$delta)
    },
    "d" = function(y, eta, log = FALSE) {
      d <- ifelse(y == 0, eta$pi + (1 - eta$pi) * dnbinom(y, mu = eta$mu, size = eta$delta), 
				(1 - eta$pi) * dnbinom(y, mu = eta$mu, size = eta$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, eta, ...) {
		  cdf1 <- pnbinom(y, size = eta$delta, mu = eta$mu)
		  cdf2 <- pnbinom(0, size = eta$delta, mu = eta$mu)
		  cdf3 <- eta$pi + ((1 - eta$pi) * (cdf1 - cdf2)/(1 - cdf2))
		  cdf <- ifelse((y == 0), eta$pi, cdf3)
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


## http://stats.stackexchange.com/questions/17672/quantile-regression-in-jags
quant.bamlss <- function(links = c(mu = "identity"), prob = 0.5, ...)
{
  rval <- list(
    "family" = "quant",
    "names" = "mu",
    "links" = parse.links(links, c(mu = "identity"), ...),
    "bayesx" = list(
      "mu" = c("quantreg", "mean"),
      "quantile" = prob
    )
  )
  class(rval) <- "family.bamlss"
  rval
}


quant2.bamlss <- function(links = c(mu = "identity", sigma = "log"), prob = 0.5, ...)
{
  rval <- list(
    "family" = "quant2",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(
        "mu" = "(1 - 2 * prop) / (prop * (1 - prop)) * w[i] + mu",
        "sigma" = "(prop * (1 - prop) * (1 / sigma)) / (2 * w[i])"
      ),
      "addparam" = list("w[i] ~ dexp(1 / sigma[i])"),
      "addvalues" = list("prop" = prob)
    )
  )
  class(rval) <- "family.bamlss"
  rval
}


## Zero adjusted families.
zero.bamlss <- function(pi = "logit", g = invgaussian)
{
  gg <- try(inherits(g, "family.bamlss"), silent = TRUE)
  if(inherits(gg, "try-error")) {
    g <- deparse(substitute(g), backtick = TRUE, width.cutoff = 500)
  } else {
    if(is.function(g)) {
      if(inherits(try(g(), silent = TRUE), "try-error"))
        g <- deparse(substitute(g), backtick = TRUE, width.cutoff = 500)
    }
  }
  g <- bamlss.family(g)
  g[c("mu", "map2par", "loglik")] <- NULL
  g0 <- g
  np <- g$names
  g$links <- c(g$links, "pi" = pi)
  g$family <- paste("zero-adjusted | ", g$family, ", ", "binomial", sep = "")
  g$names <- c(g$names, "pi")
  g$valid.response <- function(x) {
    if(any(x < 0)) stop("response includes values smaller than zero!")
    TRUE
  }
  dfun <- g0$d
  pfun <- g0$p
  if(is.function(dfun)) {
    g$d <- function(y, eta, log = TRUE) {
      d <- dfun(y, eta, log = FALSE) * eta$pi * (y > 0) + (1 - eta$pi) * (y == 0)
      if(log)
        d <- log(d)
      d
    }
  }
  if(is.function(pfun)) {
    g$p <- function(y, eta, log = FALSE) {
      pfun(y, eta, log = log) * eta$pi + (1 - eta$pi)
    }
  }
  g$bayesx <- c(g$bayesx, list("pi" = c(paste("binomial", pi, sep = "_"), "meanservant")))
  g$bayesx$weights <- list()
  for(j in np) {
    g$bayesx$weights[[j]] <- function(x) { 1 * (x > 0) }
    if(grepl("mean", g$bayesx[[j]][2]))
      g$bayesx[[j]][2] <- "meanservant"
  }
  g$bayesx$zero <- TRUE
  class(g) <- "family.bamlss"
  g
}


## General bamlss family creator.
gF <- function(x, ...) {
  x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  F <- get(paste(x, "bamlss", sep = "."), mode = "function")
  F(...)
}

gF2 <- function(x, ...) {
  x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  F <- get(paste(x, "bamlss", sep = "."), mode = "function")
  x <- F(...)
  linkinv <- vector(mode = "list", length = length(x$names))
  for(j in x$names)
    linkinv[[j]] <- make.link2(x$links[j])$linkinv
  x$map2par <- function(eta) {
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
  if(is.null(x$loglik))
    x$loglik <- function(y, eta, ...) { sum(x$d(y, eta, log = TRUE), na.rm = TRUE) }
  x
}


## Function to transform gamlss.family objects.
tF <- function(x, ...)
{
  require("gamlss")

  if(is.function(x)) x <- x()
  if(!inherits(x, "gamlss.family")) stop('only "gamlss.family" objects can be transformed!')

  args <- list(...)
  k <- x$nopar
  nx <- c("mu", "sigma", "nu", "tau")[1:k]
  nf <- names(x)
  de <- c("m", "d", "v", "t")[1:k]
  score <- weights <- list()

  args <- if(!is.null(names(args))) {
    paste(', ', paste(names(args), "=", unlist(args), sep = '', collapse = ', '), sep = '')
  } else ''

  mu.linkfun <- make.link.gamlss(x$mu.link)$linkfun
  score$mu  <- function(y, eta) {
    call <- paste('x$dldm(y, ', paste('eta$', nx, sep = '', collapse = ', '), args, ')', sep = "")
    eval(parse(text = call))
  }
  weights$mu <- function(y, eta) {
    fo <- names(formals(x$d2ldm2))
    call <- paste('x$d2ldm2(', paste('eta$', fo, sep = '', collapse = ', '), ')', sep = "")
    weights <- eval(parse(text = call))
    dlink <- 1 / x$mu.dr(mu.linkfun(eta$mu))
    weights <- -(weights / (dlink * dlink)) * dlink
    weights
  }
  if(k > 1) {
    sigma.linkfun <- make.link.gamlss(x$sigma.link)$linkfun
    score$sigma  <- function(y, eta) {
      call <- paste('x$dldd(y, ', paste('eta$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    weights$sigma <- function(y, eta) {
      fo <- names(formals(x$d2ldd2))
      call <- paste('x$d2ldd2(', paste('eta$', fo, sep = '', collapse = ', '), ')', sep = "")
      weights <- eval(parse(text = call))
      dlink <- 1 / x$sigma.dr(sigma.linkfun(eta$sigma))
      weights <- -(weights / (dlink * dlink)) * dlink
      weights
    }
  }
  if(k > 2) {
    nu.linkfun <- make.link.gamlss(x$nu.link)$linkfun
    score$nu  <- function(y, eta) {
      call <- paste('x$dldv(y, ', paste('eta$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    weights$nu <- function(y, eta) {
      fo <- names(formals(x$d2ldv2))
      call <- paste('x$d2ldv2(', paste('eta$', fo, sep = '', collapse = ', '), ')', sep = "")
      weights <- eval(parse(text = call))
      dlink <- 1 / x$nu.dr(nu.linkfun(eta$sigma))
      weights <- -(weights / (dlink * dlink)) * dlink
      weights
    }
  }
  if(k > 3) {
    tau.linkfun <- make.link.gamlss(x$tau.link)$linkfun
    score$tau  <- function(y, eta) {
      call <- paste('x$dldt(y, ', paste('eta$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    weights$tau <- function(y, eta) {
      fo <- names(formals(x$d2ldt2))
      call <- paste('x$d2ldt2(', paste('eta$', fo, sep = '', collapse = ', '), ')', sep = "")
      weights <- eval(parse(text = call))
      dlink <- 1 / x$tau.dr(tau.linkfun(eta$sigma))
      weights <- -(weights / (dlink * dlink)) * dlink
      weights
    }
  }

  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- get(paste("p", x$family[1], sep = ""))

  rval <- list(
    "family" = x$family[1],
    "names" = nx,
    "links" = unlist(x[paste(nx, "link", sep = ".")]),
    "score" = score,
    "weights" = weights,
    "d" = function(y, eta, log = FALSE, ...) {
      call <- paste('dfun(y, ', paste('eta$', nx, sep = '', collapse = ', '), ', ...)', sep = "")
      d <- try(eval(parse(text = call)), silent = TRUE)
      if(inherits(d, "try-error")) {
        d <- rep(0, length(y))
      } else {
        if(log) d <- log(d)
      }
      d
    },
    "p" = function(y, eta, ...) {
      call <- paste('pfun(y, ', paste('eta$', nx, sep = '', collapse = ', '), ', ...)', sep = "")
      eval(parse(text = call))
    }
  )
  names(rval$links) <- nx

  class(rval) <- "family.bamlss"
  rval
}


### Ordered logit.
### http://staff.washington.edu/lorenc2/bayesian/ologit.R

## Categorical distribution.
dcat <- function(x, p, log=FALSE)
{
  if(is.vector(x) & !is.matrix(p))
    p <- matrix(p, length(x), length(p), byrow = TRUE)
  if(is.matrix(x) & !is.matrix(p))
    p <- matrix(p, nrow(x), length(p), byrow = TRUE)
  if(is.vector(x) & {length(x) == 1}) {
    temp <- rep(0, ncol(p))
    temp[x] <- 1
    x <- t(temp)
  } else if(is.vector(x) & (length(x) > 1)) x <- as.indicator.matrix(x)
  if(!identical(nrow(x), nrow(p))) stop("The number of rows of x and p differ.")
  if(!identical(ncol(x), ncol(p))) {
    x.temp <- matrix(0, nrow(p), ncol(p))
    x.temp[,as.numeric(colnames(x))] <- x
    x <- x.temp
  }
  dens <- x*p
  if(log == TRUE) dens <- x*log(p)
  dens <- as.vector(rowSums(dens))
  return(dens)
}

qcat <- function(pr, p, lower.tail = TRUE, log.pr = FALSE)
{
  if(!is.vector(pr)) pr <- as.vector(pr)
  if(!is.vector(p)) p <- as.vector(p)
  if(log.pr == FALSE) {
    if(any(pr < 0) | any(pr > 1))
      stop("pr must be in [0,1].")
  } else if(any(!is.finite(pr)) | any(pr > 0)) stop("pr, as a log, must be in (-Inf,0].")
  if(sum(p) != 1) stop("sum(p) must be 1.")
  if(lower.tail == FALSE) pr <- 1 - pr
  breaks <- c(0, cumsum(p))
  if(log.pr == TRUE) breaks <- log(breaks)
  breaks <- matrix(breaks, length(pr), length(breaks), byrow = TRUE)
  x <- rowSums(pr > breaks)
  return(x)
}

rcat <- function(n, p)
{
  if(is.vector(p)) {
    x <- as.vector(which(rmultinom(n, size = 1, prob = p) == 1, arr.ind = TRUE)[, "row"])
  } else {
    n <- nrow(p)
    x <- apply(p, 1, function(x) {
      as.vector(which(rmultinom(1, size = 1, prob = x) == 1, arr.ind = TRUE)[, "row"])
    })
  }
  return(x)
}

as.indicator.matrix <- function(x)
{
  n <- length(x)
  x <- as.factor(x)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x) - 1)] <- 1
  dimnames(X) <- list(names(x), levels(x))
  X
}


#####################################
## New family specifictaion setup. ##
#####################################
##############################################
## (1) Score, hessian and fisher functions. ##
##############################################
##########################
## Normal distribution. ##
##########################
snorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  dxm <- x - mean
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      score <- cbind(score, dxm / sd2)
    if(w == "sigma")
      score <- cbind(score, (dxm^2 - sd2) / sd^3)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

hnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  n <- length(x)
  hess <- list()
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- rep(-1 / sd2, length.out = n)
    if(w == "sigma")
      hess[[w]] <- (sd2 - 3 * (x - mean)^2) / sd2^2
    if(w == "mu.sigma")
      hess[[w]] <- -2 * (x - mean) / sd^3
    if(w == "sigma.mu")
      hess[[w]] <- 2 * (mean - x) / sd^3
  }
  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

fnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  n <- length(x)
  fish <- list()
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      fish[[w]] <- rep(1 / sd2, length.out = n)
    if(w == "sigma")
      fish[[w]] <- rep(2 / sd2, length.out = n)
    if(w %in% c("mu.sigma", "sigma.mu"))
      fish[[w]] <- 0
  }
  fish <- do.call("cbind", fish)
  colnames(fish) <- gsub("mu", "dmu", colnames(fish))
  colnames(fish) <- gsub("sigma", "dsigma", colnames(fish))
  colnames(fish)[colnames(fish) == "dmu"] <- "d2mu"
  colnames(fish)[colnames(fish) == "dsigma"] <- "d2sigma"
  fish
}


###################################
## (2) Family creator functions. ##
###################################
get.dist <- function(distribution = "norm")
{
  ddist <- get(paste("d", distribution, sep = ""))
  pdist <- try(get(paste("p", distribution, sep = "")), silent = TRUE)
  if(inherits(pdist, "try-error"))
    pdist <- NULL
  qdist <- try(get(paste("q", distribution, sep = "")), silent = TRUE)
  if(inherits(qdist, "try-error"))
    qdist <- NULL
  rdist <- try(get(paste("r", distribution, sep = "")), silent = TRUE)
  if(inherits(pdist, "try-error"))
    rdist <- NULL
  sdist <- try(get(paste("s", distribution, sep = "")), silent = TRUE)
  if(inherits(sdist, "try-error"))
    sdist <- NULL
  hdist <- try(get(paste("h", distribution, sep = "")), silent = TRUE)
  if(inherits(hdist, "try-error"))
    hdist <- NULL
  fdist <- try(get(paste("f", distribution, sep = "")), silent = TRUE)
  if(inherits(fdist, "try-error"))
    fdist <- NULL
  dist <- list("d" = ddist, "p" = pdist, "q" = qdist, "r" = rdist,
    "s" = sdist, "h" = hdist, "f" = fdist)
  return(dist)
}


gaussian5.bamlss <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  links <- parse.links(links, c(mu = "identity", sigma = "log"), ...)
  lfun <- list()
  for(j in names(links))
    lfun[[j]] <- make.link2(links[j])

  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = links,
    "score" = list(
      "mu" = function(y, eta, ...) {
        mu <- lfun$mu$linkfun(eta$mu)
        drop(snorm(y, mean = eta$mu, sd = eta$sigma, which = 1) * lfun$mu$mu.eta(mu))
      },
      "sigma" = function(y, eta, ...) {
        sigma <- lfun$sigma$linkfun(eta$sigma)
        drop(snorm(y, mean = eta$mu, sd = eta$sigma, which = 2) * lfun$sigma$mu.eta(sigma))
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        mu <- lfun$mu$linkfun(eta$mu)
        w <- -1 * drop(snorm(y, mean = eta$mu, sd = eta$sigma, which = 1) * lfun$mu$mu.eta2(mu) +
          hnorm(y, mean = eta$mu, sd = eta$sigma, which = 1) * lfun$mu$mu.eta(mu)^2)
        w
      },
      "sigma" = function(y, eta, ...) {
        sigma <- lfun$sigma$linkfun(eta$sigma)
        w <- -1 * drop(snorm(y, mean = eta$mu, sd = eta$sigma, which = 2) * lfun$sigma$mu.eta2(sigma) +
          hnorm(y, mean = eta$mu, sd = eta$sigma, which = 2) * lfun$sigma$mu.eta(sigma)^2)
        w
      }
    ),
    "mu" = function(eta, ...) {
      eta$mu
    },
    "d" = function(y, eta, log = FALSE) {
      dnorm(y, mean = eta$mu, sd = eta$sigma, log = log)
    },
    "p" = function(y, eta, ...) {
      pnorm(y, mean = eta$mu, sd = eta$sigma, ...)
    },
    "type" = 1
  )
  
  class(rval) <- "family.bamlss"
  rval
}


###############
##  Survival ##
###############
cox.bamlss <- function(links = c(lambda = "identity", mu = "identity"), ...)
{
  require("survival")
  rval <- list(
    "family" = "cox",
    "names" = c("lambda", "mu"),
    "links" = parse.links(links, c(lambda = "log", mu = "identity"), ...),
    "transform" = surv.transform,
    "loglik" = function(y, eta, ...) {
      n <- attr(y, "subdivisions")
      eeta <- exp(eta_Surv_timegrid)
      int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
      ll <- (eta$lambda + eta$mu) * y[, "status"] - exp(eta$mu) * int
      sum(ll)
    },
    "score" = list(
      "mu" = function(y, eta, ...) {
        n <- attr(y, "subdivisions")
        eeta <- exp(eta_Surv_timegrid)
        int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
        y[, "status"] - exp(eta$mu) * int
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        n <- attr(y, "subdivisions")
        eeta <- exp(eta_Surv_timegrid)
        int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
        exp(eta$mu) * int
      }
    ),
    "gradient" = list(
      "lambda" = function(g, y, eta, x, ...) {
        n <- attr(y, "subdivisions")
        X <- x$get.mu_timegrid(NULL)
        eeta <- eta_Surv_timegrid + x$get.mu_timegrid(g)
        eeta <- exp(eeta)
        dummy <- vector("list", ncol(x$X))
        for(i in 1:ncol(x$X)) {
          dummy[[i]] <- matrix(X[, i], nrow = nrow(eeta), ncol = ncol(eeta), byrow = TRUE)
          dummy[[i]] <- dummy[[i]] * eeta
          dummy[[i]] <- attr(y, "width") * (0.5 * (dummy[[i]][, 1] + dummy[[i]][, n]) + apply(dummy[[i]][, 2:(n - 1)], 1, sum))
        }
        dummy <- sapply(dummy, cbind)
        dummy <- dummy * exp(eta$mu)
        int <- apply(dummy, 2, sum)
        xgrad <- drop(t(y[, "status"]) %*% x$X - int)
        return(xgrad)
      }
    ),
    "hessian" = list(
      "lambda" = function(g, y, eta, x, ...) {
        n <- attr(y, "subdivisions")
        X <- x$get.mu_timegrid(NULL)
        eeta <- eta_Surv_timegrid + x$get.mu_timegrid(g)
        eeta <- exp(eeta)
        nobs <- nrow(y)
        dummy <- vector("list", nobs)
        width <- attr(y, "width")
        xhess <- matrix(0, ncol = ncol(x$X), nrow = ncol(x$X))
        for(i in 1:nobs) {
          forward <- n * (i - 1)
          dummy[[i]] <- matrix(0, ncol = ncol(X), nrow = ncol(X))
          for(j in 1:n) {
            MAT <- X[j + forward,] %o% X[j + forward,] * eeta[i, j]
            if(j == 1 || j == n){
              dummy[[i]] <- dummy[[i]] + 0.5 * MAT
            } else {
              dummy[[i]] <- dummy[[i]] + MAT
            }
          }
          dummy[[i]] <- dummy[[i]] * width[i]
          xhess <- xhess + exp(eta$mu[i]) * dummy[[i]]
        }
        return(xhess)
      }
    )
  )
  class(rval) <- "family.bamlss"
  rval
}


Surv2 <- function(..., obs = NULL, subdivisions = 100)
{
  require("survival")
  rval <- cbind(as.matrix(Surv(...)), "obs" = obs)
  class(rval) <- c("Surv")
  rval
}


rSurvTime2 <- function (lambda, x, cens_fct, upper = 1000, ..., file = NULL,
  subdivisions = 1000) 
{
  if(!is.matrix(x)) 
    x <- cbind(x)
  time <- rep(NA, nrow(x))
  Lambda <- function(lambda, x, time) {
    integrate(lambda, 0, time, x = x, subdivisions = subdivisions)$value
  }
  InvLambda <- function(Lambda, lambda, x) {
    negLogU <- -log(runif(1, 0, 1))
    rootfct <- function(time) {
      negLogU - Lambda(lambda, x, time)
    }
    return(uniroot(rootfct, interval = c(0, upper))$root)
  }
  for(i in 1:nrow(x)) {
    time[i] = InvLambda(Lambda, lambda, x[i, ])
  }
  time_event = cens_fct(time, ...)
  data = data.frame(time = time_event[, 1], event = time_event[, 2], x = x)
  names(data) <- gsub("x.", "x", names(data), fixed = TRUE)
  if(!is.null(file)) {
    save(data, file = file)
    invisible(data)
  } else {
    return(data)
  }
}


integrate2 <- function(f, a, b, ...)
{
  nodes <- c(
    -0.999713726773441,
    -0.998491950639596,
    -0.996295134733125,
    -0.993124937037443,
    -0.988984395242991,
    -0.983877540706057,
    -0.977809358486919,
    -0.970785775763707,
    -0.962813654255816,
    -0.953900782925493,
    -0.944055870136256,
    -0.93328853504308,
    -0.921609298145334,
    -0.90902957098253,
    -0.895561644970728,
    -0.881218679385019,
    -0.866014688497165,
    -0.849964527879591,
    -0.833083879888401,
    -0.815389238339177,
    -0.796897892390315,
    -0.777627909649496,
    -0.757598118519707,
    -0.73682808980202,
    -0.715338117573057,
    -0.693149199355802,
    -0.670283015603141,
    -0.64676190851413,
    -0.622608860203708,
    -0.597847470247179,
    -0.572501932621382,
    -0.546597012065094,
    -0.520158019881763,
    -0.493210789208192,
    -0.465781649773359,
    -0.437897402172032,
    -0.409585291678302,
    -0.38087298162463,
    -0.351788526372422,
    -0.32236034390053,
    -0.292617188038472,
    -0.262588120371503,
    -0.232302481844974,
    -0.201789864095736,
    -0.171080080538604,
    -0.140203137236114,
    -0.109189203580061,
    -0.0780685828134366,
    -0.0468716824215914,
    -0.0156289844215428,
    0.0156289844215435,
    0.0468716824215913,
    0.0780685828134369,
    0.109189203580061,
    0.140203137236114,
    0.171080080538603,
    0.201789864095736,
    0.232302481844974,
    0.262588120371504,
    0.292617188038472,
    0.322360343900529,
    0.351788526372422,
    0.38087298162463,
    0.409585291678302,
    0.437897402172032,
    0.465781649773358,
    0.493210789208191,
    0.520158019881763,
    0.546597012065094,
    0.572501932621381,
    0.597847470247179,
    0.622608860203708,
    0.64676190851413,
    0.670283015603142,
    0.693149199355802,
    0.715338117573057,
    0.736828089802021,
    0.757598118519707,
    0.777627909649496,
    0.796897892390315,
    0.815389238339176,
    0.833083879888401,
    0.849964527879592,
    0.866014688497164,
    0.881218679385019,
    0.895561644970727,
    0.90902957098253,
    0.921609298145334,
    0.93328853504308,
    0.944055870136256,
    0.953900782925492,
    0.962813654255816,
    0.970785775763707,
    0.977809358486919,
    0.983877540706057,
    0.988984395242992,
    0.993124937037444,
    0.996295134733125,
    0.998491950639596,
    0.999713726773442
  )

  weights <- c(
    0.000734634490505196,
    0.00170939265351828,
    0.00268392537155396,
    0.00365596120132663,
    0.00462445006342169,
    0.00558842800386553,
    0.00654694845084563,
    0.00749907325546386,
    0.00844387146966807,
    0.00938041965369529,
    0.0103078025748694,
    0.0112251140231855,
    0.0121314576629791,
    0.0130259478929708,
    0.0139077107037196,
    0.0147758845274415,
    0.015629621077546,
    0.0164680861761455,
    0.0172904605683233,
    0.0180959407221274,
    0.0188837396133761,
    0.019653087494435,
    0.0204032326462095,
    0.0211334421125275,
    0.0218430024162464,
    0.0225312202563371,
    0.0231974231852537,
    0.0238409602659681,
    0.0244612027079576,
    0.0250575444815794,
    0.0256294029102085,
    0.0261762192395456,
    0.0266974591835708,
    0.0271926134465774,
    0.0276611982207921,
    0.0281027556591007,
    0.0285168543223954,
    0.0289030896011252,
    0.0292610841106386,
    0.0295904880599137,
    0.0298909795933315,
    0.0301622651051686,
    0.0304040795264544,
    0.0306161865839809,
    0.030798379031153,
    0.0309504788504913,
    0.0310723374275658,
    0.03116383569621,
    0.0312248842548494,
    0.0312554234538637,
    0.0312554234538625,
    0.0312248842548495,
    0.0311638356962099,
    0.0310723374275674,
    0.030950478850491,
    0.0307983790311533,
    0.0306161865839802,
    0.0304040795264556,
    0.0301622651051682,
    0.0298909795933326,
    0.0295904880599136,
    0.0292610841106384,
    0.0289030896011253,
    0.0285168543223953,
    0.0281027556591007,
    0.0276611982207923,
    0.0271926134465769,
    0.0266974591835716,
    0.0261762192395464,
    0.025629402910208,
    0.0250575444815787,
    0.024461202707957,
    0.0238409602659674,
    0.0231974231852531,
    0.0225312202563353,
    0.0218430024162472,
    0.0211334421125276,
    0.0204032326462092,
    0.0196530874944354,
    0.0188837396133751,
    0.0180959407221291,
    0.0172904605683237,
    0.0164680861761457,
    0.0156296210775464,
    0.0147758845274413,
    0.0139077107037179,
    0.0130259478929714,
    0.0121314576629803,
    0.011225114023186,
    0.0103078025748689,
    0.00938041965369477,
    0.00844387146966944,
    0.00749907325546371,
    0.00654694845084634,
    0.00558842800386576,
    0.00462445006342231,
    0.00365596120132682,
    0.00268392537155354,
    0.00170939265351811,
    0.000734634490505862
  )

  fn <- function(x) {
    (b - a) / 2 * f((b - a) / 2 * x + (a + b) / 2, ...)
  }

  sum(fn(nodes) * weights)
}


integrate3 <- function(f, a, b, m = 2, n = 40, nx = 1000, ret.fun = FALSE, plot = FALSE) {
  require("splines")
  xl <- a
  xu <- b
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(n - 1)
  k <- seq(xl - dx * (m + 1), xu + dx * (m + 1), length = n + 2 * m + 2)
  x <- seq(a, b, length = nx)
  X <- spline.des(k, x, m + 2, rep(0, nx))$design
  X1 <- spline.des(k, x, m + 2, rep(1, nx))$design
  X <- X[, -1]
  X1 <- X1[, -1]
  y <- f(x)
  bf <- lm.fit(X1, y)
  y2 <- drop(X %*% coef(bf))
  f2 <- splinefun(x, y2)
  if(plot) {
    x2 <- abs(f(x))
    i <- which(x2 <= quantile(x2, prob = 0.01))
    par(mfrow = c(2, 1))
    plot(x, f2(x), type = "l", ylab = "int(f(x))")
    abline(v = x[i], lty = 2)
    plot(x, f(x), type = "l", col = "red", lwd = 2, ylab = "f(x)")
    lines(x, drop(X1 %*% coef(bf)), col = "black")
    abline(h = 0, lty = 2)
    abline(v = x[i], lty = 2)
  }
  if(ret.fun) {
    return(f2)
  } else {
    return(f2(b) - f2(a))
  }
}

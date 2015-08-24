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
beta.bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit")

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
      "mu" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(a * h2 * log(y) - a * h2 * log(1 - y) + ((1 - b) / b) * a * (1 - a) * (-digamma(h1) + digamma(h2)))
      },
      "sigma2" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a*(1-b)/b
        h2 <- (1-a)*(1-b)/b
        drop(-(1 - b) / (b) * ( -a * digamma(h1) - (1 - a) * digamma(h2) + digamma((1 - b) / (b)) + a * log(y) + (1 - a) * log(1 - y)))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * a^2 * (1 - a)^2 * (trigamma(h1) + trigamma(h2)))
      },
      "sigma2" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * (a^2 * trigamma(h1) + (1 - a)^2 * trigamma(h2) - trigamma((1 - b) / (b))))
      }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
       mu <- par$mu
       sigma2 <- par$sigma2
       a <- mu * (1 - sigma2) / (sigma2)
       b <- a * (1 - mu) / mu
       dbeta(y, shape1 = a, shape2 = b, log = log)
    },
    "p" = function(y, par, ...) {
       mu <- par$mu
       sigma2 <- par$sigma2
       a <- mu * (1 - sigma2) / (sigma2)
       b <- a * (1 - mu) / mu
       pbeta(y, shape1 = a, shape2 = b, ...)
    },
    "type" = 1
  )
  class(rval) <- "family.bamlss"
  rval
}


betazoi.bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log")

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
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 1) & (x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 1) & (x != 0)) }
      )
    ),
    "mu" = function(par, ...) {
       par$mu * (1 - (par$nu + par$tau) / (1 + par$nu + par$tau)) + par$tau / (1 + par$nu + par$tau)
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 0, par$nu / (1 + par$nu + par$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$nu + par$tau))
      ifelse (y==1, par$tau / (1 + par$nu + par$tau), d)
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$nu / (1 + par$nu + par$tau)
      h2 <- par$tau / (1 + par$nu + par$tau)
      cdf <- ifelse(y == 0, h1, h1 + (1 - (h1 + h2)) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
      ifelse(y == 1, 1, cdf)
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


betazi.bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", nu = "log")

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
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 0)) }
      )
    ),
    "mu" = function(par, ...) {
      par$mu * (1 - (par$nu) / (1 + par$nu))
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 0, par$nu / (1 + par$nu), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$nu))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$nu / (1 + par$nu)
      ifelse(y == 0, h1, h1 + (1 - h1) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


betaoi.bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", tau = "log")

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
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 1)) },
        "sigma2" = function(x) { 1 * ((x != 1)) }
      )
    ),
    "mu" = function(par, ...) {
      par$mu * (1 - par$tau / (1 + par$tau)) +
        par$tau / (1 + par$tau)
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 1, par$tau / (1 + par$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$tau))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$tau / (1 + par$tau)
      ifelse(y == 1, h1, h1 + (1 - h1) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


binomial.bamlss <- function(...)
{
  link <- "logit"

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
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      if(is.factor(y)) y <- as.integer(y) - 1
      dbinom(y, size = 1, prob = par$pi, log = log)
    },
    "p" = function(y, par, ...) {
      y <- as.integer(y) - 1
      pbinom(y, size = 1, prob = par$pi, ...)
    },
    "score" = list(
      "pi" = function(y, par, ...) {
        if(is.factor(y)) y <- as.integer(y) - 1
        (par$pi - y) / ((par$pi - 1) * par$pi)
      }
    ),
    "hess" = list(
      "pi" = function(y, par, ...) {
        if(is.factor(y)) y <- as.integer(y) - 1
        -1 * ((par$pi^2 + y - 2 * par$pi * y) / ((-1 + par$pi^2) * par$pi^2))
      }
    ),
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}

cloglog.bamlss <- function(...)
{
  link <- "cloglog"

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
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      dbinom(y, size = 1, prob = par$pi, log = log)
    },
    "p" = function(y, par, ...) {
      pbinom(y, size = 1, prob = par$pi, ...)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

gaussian.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

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
      "mu" = function(y, par, ...) { drop((y - par$mu) / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { drop(-1 + (y - par$mu)^2 / (par$sigma^2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
    "gradient" = list(
      "mu" = function(g, y, par, x, ...) {
        par$mu <- par$mu + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop((y - par$mu) / (exp(par$sigma)^2)))
      },
      "sigma" = function(g, y, par, x, ...) {
        par$sigma <- par$sigma + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop(-1 + (y - par$mu)^2 / (exp(par$sigma)^2)))
      }
    ),
    "hessian" = list(
      "mu" = function(g, y, par, x, ...) {
        if(!is.null(x$xbin.ind)) {
          return(crossprod(x$X[x$xbin.ind, , drop = FALSE], x$X[x$xbin.ind, , drop = FALSE] * drop(1 / exp(par$sigma)^2)))
        } else {
          return(crossprod(x$X, x$X * drop(1 / exp(par$sigma)^2)))
        }
      },
      "sigma" = function(g, y, par, x, ...) {
        if(is.null(x$XX)) return(crossprod(x$X)) else return(x$XX)
      }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "type" = 1
  )
  
  class(rval) <- "family.bamlss"
  rval
}


gaussian2.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

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
      "mu" = function(y, par, ...) { drop((y - par$mu) / par$sigma2) },
      "sigma2" = function(y, par, ...) { drop(-0.5 + (y - par$mu)^2 / (2 * par$sigma2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, sqrt(par$sigma2), log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu 
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = sqrt(par$sigma2), log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = sqrt(par$sigma2), ...)
    },
    "type" = 1
  )
 
  class(rval) <- "family.bamlss"
  rval
}


truncgaussian2.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

  rval <- list(
    "family" = "truncgaussian2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal", "mu"),
      "sigma2" = c("truncnormal", "sigma2")
    ),
    "mu" = function(par, ...) {
      mu <-  par$mu
	    sigma <- sqrt(par$sigma2)
	    arg <- - mu / sigma
	    mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, par, log = FALSE) {
	    sigma <- sqrt(par$sigma2)
	    arg <- - par$mu / sigma
	    d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
	    sigma <- sqrt(par$sigma2)
	    arg <- - par$mu / sigma
	    2 * (pnorm(y / sigma + arg) - pnorm(arg))
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}

truncgaussian.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "truncgaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal2", "mu"),
      "sigma" = c("truncnormal2", "sigma")
    ),
    "mu" = function(par, ...) {
      mu <-  par$mu
	    sigma <- par$sigma
	    arg <- - mu / sigma
	    mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, par, log = FALSE) {
      mu <-  par$mu
	    sigma <- par$sigma
	    arg <- - mu / sigma
	    d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <-  par$mu
	    sigma <- par$sigma
	    arg <- - mu / sigma
	    2 * (pnorm(y / sigma + arg) - pnorm(arg))
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        rval <- with(par, (y - mu) / (sigma^2) - (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      },
      "sigma" = function(y, par, ...) {
        rval <- with(par, -1 + (y - mu)^2 / (sigma^2) + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        rval <- with(par, 1 / (sigma^2) - (mu / sigma^2) * (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))
          - ((1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))^2)
        return(drop(rval))
      },
      "sigma" = function(y, par, ...) {
        rval <- with(par, 2 - (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)) * 
          (1 + (mu / sigma)^2 + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))))
        return(drop(rval))
      }
    ),
    "loglik" = function(y, par, ...) {
      rval <- with(par, sum(-0.5 * log(2 * pi) - log(sigma) - (y - mu)^2 / (2*sigma^2) - log(pnorm(mu / sigma))))
      return(rval)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


trunc.bamlss <- function(direction = "left", point = 0, ...)
{
  links <- c(mu = "identity", sigma = "log")

  tgrad <- function(y, par, what = "mu") {
    par$sigma[par$sigma < 1] <- 1
    par$mu[par$mu < -10] <- -10
    par$mu[par$mu > 10] <- 10
    resid <- y - par$mu
    if(direction == "left") {
      trunc <- par$mu - point
      sgn <- 1
    } else {
      trunc <- point - par$mu
      sgn <- -1
    }
    mills <- dnorm(trunc / par$sigma) / pnorm(trunc / par$sigma)
    g <- if(what == "mu") {
      resid / par$sigma^2 - sgn / par$sigma * mills
    } else {
      resid^2 / par$sigma^3 - 1 / par$sigma + trunc / par$sigma^2 * mills
    }
    return(g)
  }

  rval <- list(
    "family" = "truncreg",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "d" = function(y, par, log = FALSE) {
      par$sigma[par$sigma < 1] <- 1
      par$mu[par$mu < -10] <- -10
      par$mu[par$mu > 10] <- 10
      resid <- y - par$mu
      if(direction == "left") {
        trunc <- par$mu - point
      } else {
        trunc <- point - par$mu
      }
      ll <- log(dnorm(resid / par$sigma)) - log(par$sigma) - log(pnorm(trunc / par$sigma))
      if(!log) ll <- exp(ll)
      ll
    },
    "p" = function(y, par, ...) {
      require("msm")
      ptnorm(y, mean = par$mu, sd = par$sigma,
        lower = if(direction == "left") point else -Inf,
        upper = if(direction == "right") point else Inf)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        as.numeric(tgrad(y, par, what = "mu"))
      },
      "sigma" = function(y, par, ...) {
        as.numeric(tgrad(y, par, what = "sigma"))
      }
    )
  )
  
  class(rval) <- "family.bamlss"
  rval
}


cnorm.bamlss <- function(...)
{
  f <- list(
    "family" = "cnorm",
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "identity", "sigma" = "log")
  )
  f$transform <- function(x, ...) {
    x <- bamlss.setup(x, ...)
    y <- attr(x, "response.vec")
    check <- as.integer(y <= 0)
    attr(y, "check") <- check
    attr(x, "response.vec") <- y
    x
  }
  f$engine <- function(x, ...) {
    sampler <- if(is.null(list(...)$no.mcmc)) {
      function(x, ...) { GMCMC(x, propose = "iwls", ...) }
    } else null.sampler
    optimizer <- if(is.null(list(...)$no.bfit)) bfit0 else opt0
    stacker(x, optimizer = optimizer, sampler = sampler, ...)
  }
  f$score <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_score_mu",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_score_sigma",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$hess <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_hess_mu",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_hess_sigma",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$loglik <- function(y, par, ...) {
    .Call("cnorm_loglik",
      as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
      as.integer(attr(y, "check")))
  }
  f$d <- function(y, par, log = FALSE) {
    ifelse(y <= 0, pnorm(-par$mu / par$sigma, log.p = log),
      dnorm((y - par$mu) / par$sigma, log = log) / par$sigma^(1 - log) - log(par$sigma) * log)
  }
  f$p <- function(y, par, log = FALSE) {
    ifelse(y < 0, 0, pnorm((y - par$mu) / par$sigma, log = log))
  }
  f$q <- function(y, par, ...) {
    rval <- qnorm(y) * par$sigma + par$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$r <- function(n, y, par, ...) {
    rval <- rnorm(n) * par$sigma + par$mu
    pmax(pmin(rval, Inf), 0)
  }
  class(f) <- "family.bamlss"
  f
}


pcnorm.bamlss <- function(alpha = NULL, start = 1.5, ...)
{
  f <- list(
    "family" = "pcnorm",
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "identity", "sigma" = "log")
  )
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
        if(!is.null(x$alpha$smooth$parametric)) {
          x$alpha$smooth$parametric$state$parameters["g1"] <- log(start)
          x$alpha$smooth$parametric$state$fitted.values <- x$alpha$smooth$parametric$get.mu(x$alpha$smooth$parametric$X, x$alpha$smooth$parametric$state$parameters)
        }
      }
    }
    x
  }
  f$engine <- function(x, ...) {
    optimizer <- if(is.null(list(...)$no.opt)) opt0 else bfit_cnorm
    sampler <- if(is.null(list(...)$no.mcmc)) {
      function(x, ...) { GMCMC(x, propose = "iwls", ...) }
    } else null.sampler
    stacker(x, optimizer = optimizer, sampler = sampler, ...)
  }
  f$score <- list(
    "mu" = function(y, par, ...) {
      if(!is.null(alpha))
        par$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_score_mu",
        as.numeric(y^(1 / par$alpha)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, par, ...) {
      if(!is.null(alpha))
        par$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_score_sigma",
        as.numeric(y^(1 / par$alpha)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    },
    "alpha" = function(y, par, ...) {
      .Call("cnorm_power_score_alpha",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.numeric(par$alpha), as.integer(attr(y, "check")))
    }
  )
  f$hess <- list(
    "mu" = function(y, par, ...) {
      if(!is.null(alpha))
        par$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_hess_mu",
        as.numeric(y^(1 / par$alpha)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    },
    "sigma" = function(y, par, ...) {
      if(!is.null(alpha))
        par$alpha <- rep(alpha, length = length(y))
      .Call("cnorm_hess_sigma",
        as.numeric(y^(1 / par$alpha)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")))
    }
  )
  f$loglik <- function(y, par, ...) {
    if(!is.null(alpha))
      par$alpha <- rep(alpha, length = length(y))
    .Call("cnorm_power_loglik",
      as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma), as.numeric(par$alpha),
      as.integer(attr(y, "check")))
  }
  f$d <- function(y, par, log = FALSE) {
    if(!is.null(alpha))
      par$alpha <- rep(alpha, length = length(y))
    dy <- ifelse(y <= 0, pnorm(0, par$mu, par$sigma, log.p = TRUE),
      dnorm(y^(1 / par$alpha), par$mu, par$sigma, log = TRUE) -
      log(par$alpha) + (1 / par$alpha - 1) * log(y))
    if(!log)
      dy <- exp(dy)
    dy
  }
  f$p <- function(y, par, log = FALSE) {
    if(!is.null(alpha))
      par$alpha <- rep(alpha, length = length(y))
    ifelse(y <= 0, 0, pnorm(y^(1 / par$alpha), par$mu, par$sigma, log = log))
  }
  f$q <- function(y, par, ...) {
    if(!is.null(alpha))
      par$alpha <- rep(alpha, length = length(y))
    rval <- qnorm(y^(1 / par$alpha), par$mu, par$sigma)
    pmax(pmin(rval, Inf), 0)
  }
  f$r <- function(n, y, par, ...) {
    rval <- rnorm(n, par$mu, par$sigma)
    pmax(pmin(rval, Inf), 0)
  }
  class(f) <- "family.bamlss"
  f
}



cens.bamlss <- function(links = c(mu = "identity", sigma = "log", df = "log"),
  left = 0, right = Inf, dist = "gaussian", ...)
{
  dist <- match.arg(dist, c("student", "gaussian", "logistic"))

  ddist <- switch(dist,
    "student"  = function(x, location, scale, df, log = TRUE) 
      dt((x - location)/scale, df = df, log = log)/scalexp(1-log) - 
      log*log(scale),
    "gaussian" = function(x, location, scale, df, log = TRUE) 
      dnorm((x - location)/scale, log = log)/scalexp(1-log) - 
      log*log(scale),
    "logistic" = function(x, location, scale, df, log = TRUE) 
      dlogis((x - location)/scale, log = log)/scalexp(1-log) - 
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

  gradfun <- function(y, par, type = "gradient", name = "mu") {
    ## functions used to evaluate gradient and hessian
    mills <- function(y, lower.tail = TRUE) {
      with(par, sigma * ddist(y, mu, sigma, df, log = FALSE)/
      pdist(y, mu, sigma, df, log.p = FALSE, lower.tail = lower.tail))
    }

    ## ddensity/dmu
    d1 <- with(par, switch(dist,
      "student"  = function(x)  
        (x - mu)/sigma^2 * (df + 1) / (df + (x - mu)^2/sigma^2),
      "gaussian" = function(x) 
        (x - mu)/sigma^2,
      "logistic" = function(x)  
        (1 - 2 * pdist(-x, - mu, sigma, log.p = FALSE))/sigma)
    )
    
    ## ddensity/dsigma
    d2 <- function(x) with(par, d1(x) * (x-mu))

    ## d^2density/dmu^2
    d3 <- with(par, switch(dist,
      "student"  = function(x)  
        (df + 1)*((x - mu)^2 - df*sigma^2) / (df*sigma^2 + (x - mu)^2)^2,
      "gaussian" = function(x) 
        - 1/sigma^2,
      "logistic" = function(x)  
        - 2/sigma * ddist(x, mu, sigma, log = FALSE))
    )    
    
    ## d^2density/dsigma^2
    d5 <- with(par, switch(dist,
      "student"  = function(x)  
        - (x - mu)^2 * (df + 1) / (df*sigma^2 + (x - mu)^2)^2*2*df*sigma^2,
      "gaussian" = function(x) 
        2 * d3(x) * (x-mu)^2,
      "logistic" = function(x)  
        - d2(x) - 2*(x-mu)^2/sigma*ddist(x,mu,sigma, log = FALSE)
    ))
      
    ## d^2density/dmudsigma
    d4 <- with(par, switch(dist,
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
        rval <- with(par, ifelse(y <= left, 
          - mills(left)/sigma,
          ifelse(y >= right, 
            mills(right, lower.tail = FALSE)/sigma,
            d1(y)
          )))
      } else {
        rval <- with(par, ifelse(y <= left, 
          - mills(left) * (left - mu)/sigma,
          ifelse(y >= right, 
            mills(right, lower.tail = FALSE) * (right - mu)/sigma,
            d2(y) - 1
          )))
      }

    ## compute hessian
    } else {
      if(name == "mu") {
        rval <- with(par, ifelse(y <= left, 
          -d1(left)/sigma * mills(left) - mills(left)^2/sigma^2,
          ifelse(y >= right, 
            d1(right)/sigma * mills(right, lower.tail = FALSE) - 
              mills(right, lower.tail = FALSE)^2/sigma^2,
            d3(y)
          )))
      } else {
        rval <- with(par, ifelse(y <= left, 
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
    "mu" =  function(y, par, ...) {
      gradmu <- gradfun(y, par, type = "gradient", name = "mu")
      return(drop(gradmu))
    },
    "sigma" =  function(y, par, ...) {
      gradsigma <- gradfun(y, par, type = "gradient", name = "sigma")
      return(drop(gradsigma))
    }
  )

  hess <- list(
    "mu" =  function(y, par, ...) {
      wmu <- -1 * gradfun(y, par, type = "hess", name = "mu")
      return(drop(wmu))
    },
    "sigma" =  function(y, par, ...) {
      wsigma <- -1 * gradfun(y, par, type = "hess", name = "sigma")
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
    "d" = function(y, par, log = FALSE, ...) {
      ll <- with(par, ifelse(y <= left,
        pdist(left, mu, sigma, df, lower.tail = TRUE, log = TRUE),
        ifelse(y >= right,
          pdist(right, mu, sigma, df, lower.tail = FALSE, log = TRUE),
          ddist(y, mu, sigma, df, log = TRUE))))
      if(!log) ll <- exp(ll)
      return(ll)
    },
    "p" = function(y, par, log = FALSE, ...) {
      with(par, pdist(y, mu, sigma, df, lower.tail = TRUE, log = log))
    },
    "score" = score,
    "hess" = hess,
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
      dt((x - location)/scale, df = df, log = log)/scalexp(1-log) - 
      log*log(scale),
    "gaussian" = function(x, location, scale, df, log = TRUE) 
      dnorm((x - location)/scale, log = log)/scalexp(1-log) - 
      log*log(scale),
    "logistic" = function(x, location, scale, df, log = TRUE) 
      dlogis((x - location)/scale, log = log)/scalexp(1-log) - 
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
      "mu" =  function(y, par, ...) {
         gradmu <- with(par, ifelse(y <= left, 
          - ddist(left, mu, sigma, df, log = FALSE) /
            pdist(left, mu, sigma, df, log.p = FALSE),
          ifelse(y >= right, 
          ddist(right, mu, sigma, df, log = FALSE) /
            pdist(right, mu, sigma, df, lower.tail = FALSE, log.p = FALSE),
          - dddist(y, mu, sigma, df)/ddist(y, mu, sigma, df, log = FALSE))))
        return(drop(gradmu))
      },
      "sigma" =  function(y, par, ...) {
        gradsigma <- with(par, ifelse(y <= left, 
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
    "d" = function(y, par, log = FALSE, ...) {
      ll <- with(par, ifelse(y <= left,
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
    "d" = function(y, par, log = FALSE) {
      if(name == "gamma") {
		    a2 <- par$sigma
		    s2 <- par$mu / par$sigma
        d <- dtrunc(y, name, a = a, b = b, shape = s2, scale = s2, log = log)
      } else {
        d <- dtrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, log = log)
      }
      return(d)
    },
    "p" = function(y, par, ...) {
      if(name == "gamma") {
		    a2 <- par$sigma
		    s2 <- par$mu / par$sigma
        p <- ptrunc(y, name, a = a, b = b, shape = a2, scale = s2, ...)
      } else {
        p <- ptrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, ...)
      }
      return(p)
    },
    "q" = function(y, par, ...) {
      q <- qtrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, ...)
      return(q)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


t.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log", df = "log")

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
    "mu" = function(par, ...) {
      rval <- par$mu
      rval[par$df <= 1] <- 0
      rval
    },
    "d" = function(y, par, log = FALSE) {
      arg <- (y - par$mu) / sqrt(par$sigma2)
      dt(arg, df = par$df, log = log)
    },
    "p" = function(y, par, ...) {
      arg <- (y - par$mu) / sqrt(par$sigma2)
      pt(arg, df = par$df, ...)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


invgaussian.bamlss <- function(...)
{
  links <- c(mu = "log", sigma2 = "log")

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
      "mu" = function(y, par, ...) {
        mu <- par$mu
        (y - mu) / (mu^2 * par$sigma2)
      },
      "sigma2" = function(y, par, ...) {
        mu <- par$mu
        -0.5 + (y - mu)^2 / (2 * y * mu^2 * par$sigma2)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$mu * par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- sqrt(par$sigma2)
      d <- exp( -0.5 * log(2 * pi) - log(sigma) - (3 / 2) * log(y) - ((y - mu)^2) / (2 * sigma^2 * (mu^2) * y))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      lambda <- 1 / sqrt(par$sigma2)
      lq <- sqrt(lambda / y)
      qm <- y / mu
      pnorm(lq * (qm - 1)) + exp(2 * lambda / mu) * pnorm(-lq * (qm + 1), ...)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

weibull.bamlss <- function(...)
{
  links <- c(lambda = "log", alpha = "log")

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
    "mu" = function(par, ...) {
      alpha <-  par$alpha
      lambda <- par$lambda
      alpha * gamma(1 + 1 / lambda)
    },
    "d" = function(y, par, log = FALSE) {
      alpha <-  par$alpha
      lambda <- par$lambda
      dweibull(y, scale = lambda, shape = alpha, log = log)
    },
    "p" = function(y, par, ...) {
      alpha <- par$alpha
      lambda <- par$lambda
      rweibull(y, scale = lambda, shape = alpha, ...)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


pareto.bamlss <- function(...)
{
  links <- c(b = "log", p = "log")

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
    "mu" = function(par, ...) {
      p <- par$p
      b <- par$b
      p * gamma(1 + 1 / b)
    },
    "d" = function(y, par, log = FALSE) {
      p <- par$p
      b <- par$b
      d <- p * b^p * (y + b)^(-p - 1)
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      p <- par$p
      b <- par$b
      1 - ((b/(b + y))^(p))
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


gamma.bamlss <- function(...)
{
  links <- c(mu = "log", sigma = "log")

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
      "mu" = function(y, par, ...) {
        sigma <- par$sigma
        sigma * (-1 + y / par$mu)
      },
      "sigma" = function(y, par, ...) {
        mu <- par$mu
        sigma <- par$sigma
        sigma * (log(sigma) + 1 - log(mu) - digamma(sigma) + log(y) - y / mu)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { par$sigma },
      "sigma" = function(y, par, ...) {
        sigma <- par$sigma
        sigma^2 * trigamma(sigma) - sigma
      }
    ),
    "loglik" = function(y, par, ...) {
		  a <- par$sigma
		  s <- par$mu / par$sigma 
		  sum(dgamma(y, shape = a, scale = s, log = TRUE), na.rm = TRUE)
    },
    "mu" = function(par, ...) {
      par$mu
    },
	  "d" = function(y, par, log = FALSE) {
		  a <- par$sigma
		  s <- par$mu / par$sigma
		  dgamma(y, shape = a, scale = s, log = log)
	  },
	  "p" = function(y, par, lower.tail = TRUE, log.p = FALSE) {
		  a <- par$sigma
		  s <- par$mu / par$sigma
		  pgamma(y, shape = a, scale = s, lower.tail = lower.tail, log.p = log.p)
	  },
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}


gengamma.bamlss <- function(...)
{
  links <- c(mu = "log", sigma = "log", tau = "log")

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


lognormal.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

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
      "mu" = function(y, par, ...) { (log(y) - par$mu) / (par$sigma^2) },
      "sigma" = function(y, par, ...) { -1 + (log(y) - par$mu)^2 / (par$sigma^2) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$sigma^2) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
	  "mu" = function(par, ...) {
      exp(par$mu + 0.5 * (par$sigma)^2)
    },
    "d" = function(y, par, log = FALSE) {
      dlnorm(y, meanlog = par$mu, sdlog = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      plnorm(y, meanlog = par$mu, sdlog = par$sigma, ...)
    },
    "type" = 1
  )

  class(rval) <- "family.bamlss"
  rval
}


lognormal2.bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

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
      "mu" = function(y, par, ...) { (log(y) - par$mu) / (par$sigma2) },
      "sigma2" = function(y, par, ...) { -0.5 + (log(y) - par$mu)^2 / (2 * par$sigma2) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(par, ...) {
      exp(par$mu + 0.5 * (par$sigma2))
    },
    "d" = function(y, par, log = FALSE) {
      dlnorm(y, meanlog = par$mu, sdlog = sqrt(par$sigma2), log = log)
    },
    "p" = function(y, par, ...) {
      plnorm(y, meanlog = par$mu, sdlog = sqrt(par$sigma2), ...)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


dagum.bamlss <- function(...)
{
  links <- c(a = "log", b = "log", p = "log")

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
    "mu" = function(par, ...) {
      a <- par$a
      b <- par$b
      p <- par$p
      -(b/a) * (gamma(- 1 / a) * gamma(p + 1 / a)) / (gamma(p))
    },
    "d" = function(y, par, log = FALSE) {
      a <- par$a
      b <- par$b
      p <- par$p
      ap <- a * p
      d <- ap * y^(ap - 1) / (b^ap * (1 + (y / b)^a)^(p + 1))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      a <- par$a
      b <- par$b
      p <- par$p
      (1 + (y / b)^(-a))^(-p)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


BCCG2.bamlss <- function(...)
{
  links <- c(mu = "log", sigma = "log", nu = "identity")

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
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      nu <- par$nu
      z <- ifelse(nu == 0, log(y/mu)/sigma, (((y/mu)^nu - 1)/(nu * sigma)))
      d <- (1 / (sqrt(2 * pi) * sigma)) * (y^(nu - 1) / mu^nu) * exp(-z^2 / 2)
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      nu <- par$nu
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


mvn.bamlss <- function(...)
{
  links <- c(mu1 = "identity", mu2 = "identity",
    sigma1 = "log", sigma2 = "log", rho = "rhogit")

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
    "mu" = function(par, ...) {
      cbind(par$mu1, par$mu2)
    },
    "d" = function(y, par, log = FALSE) {
      cbind(
        dnorm(y[, 1], mean = par$mu1, sd = par$sigma1, log = log),
        dnorm(y[, 2], mean = par$mu2, sd = par$sigma2, log = log)
      )
    },
    "p" = function(y, par, ...) {
      cbind(
        pnorm(y[, 1], mean = par$mu1, sd = par$sigma1, ...),
        pnorm(y[, 2], mean = par$mu2, sd = par$sigma2, ...)
      )
    },
    "type" = 2
  )

  class(rval) <- "family.bamlss"
  rval
}


bivprobit.bamlss <- function(...)
{
  links <- c(mu1 = "identity", mu2 = "identity", rho = "rhogit")

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
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    },
    "type" = 2
  )

  class(rval) <- "family.bamlss"
  rval
}


bivlogit.bamlss <- function(...)
{
  links <- c(p1 = "logit", p2 = "logit", psi = "log")

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
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


mvt.bamlss <- function(...)
{
  links <- c(mu1 = "identity", mu2 = "identity",
    sigma1 = "log", sigma2 = "log", rho = "fisherz", df = "log")
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
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    },
    "type" = 2
  )
  class(rval) <- "family.bamlss"
  rval
}


dirichlet.bamlss <- function(...)
{
  link <- "logit"

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


multinomial.bamlss <- multinom.bamlss <- function(...)
{
  link <- "log"

  rval <- list(
    "family" = "multinomial",
    "names" = "pi",
    "links" = link,
    "bugs" = list(
      "dist" = "dcat",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "bayesx" = list(
      "pi" = c(paste("multinom", link, sep = "_"), "mu", "meanservant")
    ),
    "score" = function(y, par, id, ...) {
      pi <- par[[id]] / (1 + rowSums(do.call("cbind", par)))
      if(is.factor(y))
        1 * (y == id)
      else
        return(y[, id] - pi)
    },
    "hess" = function(y, par, id, ...) {
      pi <- par[[id]] / (1 + rowSums(do.call("cbind", par)))
      return(pi * (1 - pi))
    },
    "d" = function(y, par, log = FALSE) {
      if(is.factor(y))
        y <- model.matrix(~ y - 1)
      par <- cbind(do.call("cbind", par), 1)
      d1 <- rowSums(y * log(par))
      d2 <- log(rowSums(par))
      d <- d1 - d2
      if(!log)
        d <- exp(d)
      return(d)
    },
    "loglik" = function(y, par, ...) {
      if(is.factor(y))
        y <- model.matrix(~ y - 1)
      par <- cbind(do.call("cbind", par), 1)
      d1 <- rowSums(y * log(par))
      d2 <- log(rowSums(par))
      d <- d1 - d2
      return(sum(d, na.rm = TRUE))
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


## Count Data distributions
poisson.bamlss <- function(...)
{
  links <- c(lambda = "log")

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
    "mu" = function(par, ...) {
       par$lambda
    },
    "d" = function(y, par, log = FALSE) {
      dpois(y, lambda = par$lambda, log = log)
    },
    "p" = function(y, par, ...) {
      ppois(y, lambda = par$lambda, ...)
    },
    "score" = list(
      "lambda" = function(y, par, ...) {
        y / par$lambda - 1
      }
    ),
    "hess" = list(
      "lambda" = function(y, par, ...) {
        1 / par$lambda
      }
    ),
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


zip.bamlss <- function(...)
{
  links <- c(lambda = "log", pi = "logit")

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
	  "mu" = function(par, ...) {
      par$lambda * (1 - par$pi)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dpois(y, lambda = par$lambda), 
				(1 - par$pi) * dpois(y, lambda = par$lambda))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      ifelse(y < 0, 0, par$pi + (1 - par$pi) * ppois(y, lambda = par$lambda))
    },
    "nscore" = TRUE,
    "type" = 3
  )
  if(rval$bayesx[[2]][[1]] == "zip_pi_cloglog")
    rval$bayesx[[1]][[1]] <- "zip_lambda_cloglog"

  class(rval) <- "family.bamlss"
  rval
}


hurdleP.bamlss <- function(...)
{
  links <- c(lambda = "log", pi = "logit")

  rval <- list(
    "family" = "hurdle",
    "names" = c("lambda", "pi"),
    "links" = parse.links(links, c(lambda = "log", pi = "logit"), ...),
    "bayesx" = list(
      "lambda" = c("hurdle", "lambda"),
      "pi" = c("hurdle", "pi"),
      "hess" = list(
        "lambda" = function(x) { 1 * (x != 0)}
      )
    ),
    "mu" = function(par, ...) {
      (1 - par$pi) * par$lambda / (1 - exp(-par$lambda))
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi, 
        (1 - par$pi) * dpois(y, lambda = par$lambda) / (1 - exp(-par$lambda)))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      cdf1 <- ppois(y, lambda = par$lambda)
      cdf2 <- ppois(0, lambda = par$lambda)
      cdf3 <- par$pi + ((1 - par$pi) * (cdf1 - cdf2)/(1 - cdf2))
      cdf <- ifelse((y == 0), par$pi, cdf3)
      cdf
    },
    "nscore" = TRUE,
    "type" = 3
  )
 
  class(rval) <- "family.bamlss"
  rval
}

negbin.bamlss <- function(...)
{
  links <- c(mu = "log", delta = "log")

  rval <- list(
    "family" = "negbin",
    "names" = c("mu", "delta"),
    "links" = parse.links(links, c(mu = "log", delta = "log"), ...),
    "bayesx" = list(
      "mu" = c("negbin", "mu"),
      "delta" = c("negbin", "delta")
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnbinom(y, mu = par$mu, size = par$delta, log = log)
    },
    "p" = function(y, par, ...) {
      pnbinom(y, mu = par$mu, size = par$delta)
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


zinb.bamlss <- function(...)
{
  links <- c(mu = "log", pi = "logit", delta = "log")

  rval <- list(
    "family" = "zinb",
    "names" = c("mu", "pi", "delta"),
    "links" = parse.links(links, c(mu = "log", pi = "logit", delta = "log"), ...),
    "bayesx" = list(
      "mu" = c("zinb", "mu"),
      "pi" = c("zinb", "pi"),
      "delta" = c("zinb", "delta")
    ),
	  "mu" = function(par, ...) {
      par$mu * (1 - par$pi)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta), 
				(1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      ifelse(y<0, 0, par$pi + (1 - par$pi) * pnbinom(y, size = par$delta, mu = par$mu))
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}

hurdleNB.bamlss <- function(...)
{
  links <- c(mu = "log", pi = "logit", delta = "log")

  rval <- list(
    "family" = "hurdleNB",
    "names" = c("mu", "delta", "pi"),
    "links" = parse.links(links, c(mu = "log", delta = "log", pi = "logit"), ...),
    "bayesx" = list(
      "mu" = c("hurdle", "mu"),
      "delta" = c("hurdle", "delta"),
      "pi" = c("hurdle", "pi"),
      "hess" = list(
         "mu" = function(x) { 1 * (x != 0)},
         "delta" = function(x) { 1 * (x != 0)}
      )
    ),
    "mu" = function(par, ...) {
      (1 - par$pi) * par$mu / (1 - (par$delta) / (par$delta + par$mu)^par$delta)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta), 
        (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      cdf1 <- pnbinom(y, size = par$delta, mu = par$mu)
      cdf2 <- pnbinom(0, size = par$delta, mu = par$mu)
      cdf3 <- par$pi + ((1 - par$pi) * (cdf1 - cdf2)/(1 - cdf2))
      cdf <- ifelse((y == 0), par$pi, cdf3)
    },
    "nscore" = TRUE,
    "type" = 3
  )

  class(rval) <- "family.bamlss"
  rval
}


## http://stats.stackexchange.com/questions/17672/quantile-regression-in-jags
quant.bamlss <- function(prob = 0.5, ...)
{
  links <- c(mu = "identity")
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


quant2.bamlss <- function(prob = 0.5, ...)
{
  links <- c(mu = "identity", sigma = "log")
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
zero.bamlss <- function(g = invgaussian)
{
  pi <- "logit"
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
    g$d <- function(y, par, log = TRUE) {
      d <- dfun(y, par, log = FALSE) * par$pi * (y > 0) + (1 - par$pi) * (y == 0)
      if(log)
        d <- log(d)
      d
    }
  }
  if(is.function(pfun)) {
    g$p <- function(y, par, log = FALSE) {
      pfun(y, par, log = log) * par$pi + (1 - par$pi)
    }
  }
  g$bayesx <- c(g$bayesx, list("pi" = c(paste("binomial", pi, sep = "_"), "meanservant")))
  g$bayesx$hess <- list()
  for(j in np) {
    g$bayesx$hess[[j]] <- function(x) { 1 * (x > 0) }
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
    x$loglik <- function(y, par, ...) { sum(x$d(y, par, log = TRUE), na.rm = TRUE) }
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
  score <- hess <- list()

  args <- if(!is.null(names(args))) {
    paste(', ', paste(names(args), "=", unlist(args), sep = '', collapse = ', '), sep = '')
  } else ''

  mu.linkfun <- make.link.gamlss(x$mu.link)$linkfun
  score$mu  <- function(y, par, ...) {
    call <- paste('x$dldm(y, ', paste('par$', nx, sep = '', collapse = ', '), args, ')', sep = "")
    eval(parse(text = call))
  }
  hess$mu <- function(y, par, ...) {
    fo <- names(formals(x$d2ldm2))
    call <- paste('x$d2ldm2(', paste('par$', fo, sep = '', collapse = ', '), ')', sep = "")
    hess <- eval(parse(text = call))
    dlink <- 1 / x$mu.dr(mu.linkfun(par$mu))
    hess <- -(hess / (dlink * dlink)) * dlink
    hess
  }
  if(k > 1) {
    sigma.linkfun <- make.link.gamlss(x$sigma.link)$linkfun
    score$sigma  <- function(y, par, ...) {
      call <- paste('x$dldd(y, ', paste('par$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    hess$sigma <- function(y, par, ...) {
      fo <- names(formals(x$d2ldd2))
      call <- paste('x$d2ldd2(', paste('par$', fo, sep = '', collapse = ', '), ')', sep = "")
      hess <- eval(parse(text = call))
      dlink <- 1 / x$sigma.dr(sigma.linkfun(par$sigma))
      hess <- -(hess / (dlink * dlink)) * dlink
      hess
    }
  }
  if(k > 2) {
    nu.linkfun <- make.link.gamlss(x$nu.link)$linkfun
    score$nu  <- function(y, par, ...) {
      call <- paste('x$dldv(y, ', paste('par$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    hess$nu <- function(y, par, ...) {
      fo <- names(formals(x$d2ldv2))
      call <- paste('x$d2ldv2(', paste('par$', fo, sep = '', collapse = ', '), ')', sep = "")
      hess <- eval(parse(text = call))
      dlink <- 1 / x$nu.dr(nu.linkfun(par$sigma))
      hess <- -(hess / (dlink * dlink)) * dlink
      hess
    }
  }
  if(k > 3) {
    tau.linkfun <- make.link.gamlss(x$tau.link)$linkfun
    score$tau  <- function(y, par, ...) {
      call <- paste('x$dldt(y, ', paste('par$', nx, sep = '', collapse = ', '), ')', sep = "")
      eval(parse(text = call))
    }
    hess$tau <- function(y, par, ...) {
      fo <- names(formals(x$d2ldt2))
      call <- paste('x$d2ldt2(', paste('par$', fo, sep = '', collapse = ', '), ')', sep = "")
      hess <- eval(parse(text = call))
      dlink <- 1 / x$tau.dr(tau.linkfun(par$sigma))
      hess <- -(hess / (dlink * dlink)) * dlink
      hess
    }
  }

  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- get(paste("p", x$family[1], sep = ""))

  rval <- list(
    "family" = x$family[1],
    "names" = nx,
    "links" = unlist(x[paste(nx, "link", sep = ".")]),
    "score" = score,
    "hess" = hess,
    "d" = function(y, par, log = FALSE, ...) {
      call <- paste('dfun(y, ', paste('par$', nx, sep = '', collapse = ', '), ', ...)', sep = "")
      d <- try(eval(parse(text = call)), silent = TRUE)
      if(inherits(d, "try-error")) {
        d <- rep(0, length(y))
      } else {
        if(log) d <- log(d)
      }
      d
    },
    "p" = function(y, par, ...) {
      call <- paste('pfun(y, ', paste('par$', nx, sep = '', collapse = ', '), ', ...)', sep = "")
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
  dens <- x * p
  if(log) dens <- x * log(p)
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
## New family specification setup. ##
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
      "mu" = function(y, par, ...) {
        mu <- lfun$mu$linkfun(par$mu)
        drop(snorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta(mu))
      },
      "sigma" = function(y, par, ...) {
        sigma <- lfun$sigma$linkfun(par$sigma)
        drop(snorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta(sigma))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        mu <- lfun$mu$linkfun(par$mu)
        w <- -1 * drop(snorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta2(mu) +
          hnorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta(mu)^2)
        w
      },
      "sigma" = function(y, par, ...) {
        sigma <- lfun$sigma$linkfun(par$sigma)
        w <- -1 * drop(snorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta2(sigma) +
          hnorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta(sigma)^2)
        w
      }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "type" = 1
  )
  
  class(rval) <- "family.bamlss"
  rval
}


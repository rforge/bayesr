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
      "mu" = c("beta_mu", "mean"),
      "sigma2" = c("beta_sigma2", "scale")
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
      "mu" = c("betainf_mu", "location"),
      "sigma2" = c("betainf_sigma2", "scale"),
      "nu" = c("betainf_nu", "shape"),
      "tau" = c("betainf_tau", "mean"),
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
      "mu" = c("beta_mu", "location"),
      "sigma2" = c("beta_sigma2", "scale"),
      "nu" = c("betainf0_nu", "mean"),
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
      "mu" = c("beta_mu", "location"),
      "sigma2" = c("beta_sigma2", "scale"),
      "tau" = c("betainf1_tau", "mean"),
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
      "pi" = c(paste("binomial", link, sep = "_"), "mean")
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
        "identity" = c("normal2_mu", "mean"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal2_sigma", "scale"),
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
        "identity" = c("normal_mu", "mean"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma2" = switch(links["sigma2"],
        "log" = c("normal_sigma2", "scale"),
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
      "mu" =  c("truncnormal_mu", "mean"),
      "sigma2" = c("truncnormal_sigma2", "scale")
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
      "mu" =  c("truncnormal2_mu", "mean"),
      "sigma" = c("truncnormal2_sigma", "scale")
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
      "mu" = function(y, eta) {
        as.numeric(tgrad(y, eta, what = "mu"))
      },
      "sigma" = function(y, eta) {
        as.numeric(tgrad(y, eta, what = "sigma"))
      }
    )
  )
  
  class(rval) <- "family.bamlss"
  rval
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
    "mu" =  function(y, eta) {
      gradmu <- gradfun(y, eta, type = "gradient", name = "mu")
      return(drop(gradmu))
    },
    "sigma" =  function(y, eta) {
      gradsigma <- gradfun(y, eta, type = "gradient", name = "sigma")
      return(drop(gradsigma))
    }
  )

  weights <- list(
    "mu" =  function(y, eta) {
      wmu <- -1 * gradfun(y, eta, type = "weights", name = "mu")
      return(drop(wmu))
    },
    "sigma" =  function(y, eta) {
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
      "mu" = c("t_mu", "mean"),
      "sigma2" =  c("t_sigma2", "scale"),
	    "df" = c("t_df", "df")
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
      "mu"  = c("invgaussian_mu", "mean"),
      "sigma2" = c("invgaussian_sigma2", "scale")
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
      "lambda"  = c("weibull_lambda", "mean"),
      "alpha" = c("weibull_alpha", "shape")
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
      "b"  = c("pareto_b", "mean"),
      "p" = c("pareto_p", "shape")
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
      "mu" = c("gamma_mu", "mean"),
      "sigma" = c("gamma_sigma", "shape")
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
      "mu" = c("gengamma_mu", "mean"),
      "sigma" = c("gengamma_sigma", "shape1"),
	    "tau" = c("gengamma_tau", "shape2")
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
      "mu" = c("lognormal2_mu", "mean"),
      "sigma" = c("lognormal2_sigma", "scale")
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
      "mu" = c("lognormal_mu", "mean"),
      "sigma2" = c("lognormal_sigma2", "scale")
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
      "a" = c("dagum_a", "mean"),
      "b" = c("dagum_b", "scale"),
      "p" = c("dagum_p", "shape2")
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
      "mu" = c("BCCG_mu", "mean"),
      "sigma" =  c("BCCG_sigma", "scale"),
	    "nu" = c("BCCG_nu", "nu"),
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
      "mu1" = c("bivprobit_mu", "mean"),
      "mu2" = c("bivprobit_mu", "mu"),
      "rho" = c("bivprobit_rho", "rho"),
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
      "mu1" = c("bivlogit_mu", "mean"),
      "mu2" = c("bivlogit_mu", "mu"),
      "psi" = c("bivlogit_or", "oddsratio"),
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
      "mu1" = c("bivt_mu", "mean"),
      "mu2" = c("bivt_mu", "mu"),
      "sigma1" = c("bivt_sigma", "scale1"),
      "sigma2" = c("bivt_sigma", "scale2"),
      "rho" = c("bivt_rho", "rho"),
      "df" = c("bivt_df", "df"),
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
      "alpha" = c("dirichlet", "mean", "alpha")
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
      "pi" = c(paste("multinom", link, sep = "_"), "mean", "meanservant")
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
      "lambda" = c("poisson", "mean")
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
      "lambda" = c("zip_lambda", "mean"),
      "pi" = switch(links["pi"],
        "logit" = c("zip_pi", "pi"),
        "cloglog2" = c("zip_pi_cloglog", "pi")
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
      "lambda" = c("hurdle_lambda", "mean"),
      "pi" = c("hurdle_pi", "pi"),
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
      "mu" = c("negbin_mu", "mean"),
      "delta" = c("negbin_delta", "delta")
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
      "mu" = c("zinb_mu", "mean"),
      "pi" = c("zinb_pi", "pi"),
      "delta" = c("zinb_delta", "delta")
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
      "mu" = c("hurdle_mu", "mean"),
      "delta" = c("hurdle_delta", "delta"),
	    "pi" = c("hurdle_pi", "pi"),
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


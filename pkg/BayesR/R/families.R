######################
## BayesR Families. ##
######################
## Print method.
print.family.BayesR <- function(x, ...)
{
  cat("Family:", x$family, "\n")
  links <- paste(names(x$links), x$links, sep = " = ")
  links <- paste(links, collapse = ", ")
  cat(if(length(links) > 1) "Link functions:" else "Link function:", links, sep = " ")
  cat("\n")
}

## Extract method.
family.bayesr <- function(object, ...)
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
beta.BayesR <- function(links = c(mu = "logit", sigma = "log"), ...)
{
  links <- parse.links(links, c(mu = "logit", sigma = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "beta",
    "names" = c("mu", "sigma"),
    "links" =  links,
    "valid.response" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in (0, 1)!", call. = FALSE)
      ok
    },
    bayesx = list(
      "mu" = c("beta_mu", "mean"),
      "sigma" = c("beta_sigma2", "scale")
    ),
    jagstan = list(
      "dist" = "dbeta",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(
        mu = "mu * (1 / sigma)",
        sigma = "(1 - mu) * (1 / sigma)"
      )
    ),
    "loglik" = function(y, eta, ...) {
      a <- linkinv$mu(eta$mu)
      b <- linkinv$sigma(eta$sigma)
      hilfs <- a * (1 - b) / b
      hilfs2 <- (1 - a) * (1 - b) / b
      sum((hilfs - 1) * log(y) + (hilfs2 - 1) * log(1 - y) - lgamma(hilfs) - lgamma(hilfs2) + lgamma((1 - b) / b))
    },
    "score" = list(
      "mu" = function(y, eta, ...) {
        a <- linkinv$mu(eta$mu)
        b <- linkinv$sigma(eta$sigma)
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(a * hilfs2 * log(y) - a * hilfs2 * log(1 - y) +
               ((1 - b) / b) * a * (1 - a) * (-digamma(hilfs) + digamma(hilfs2)))
      },
      "sigma" = function(y, eta, ...) {
        a <- linkinv$mu(eta$mu)
        b <- linkinv$sigma(eta$sigma)
        hilfs <- a*(1-b)/b
        hilfs2 <- (1-a)*(1-b)/b
        drop(-(1 - b) / (b) * ( -a * digamma(hilfs) - (1 - a) * digamma(hilfs2) +
                                  digamma((1 - b) / (b)) + a * log(y) + (1 - a) * log(1 - y)))
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        a <- linkinv$mu(eta$mu)
        b <- linkinv$sigma(eta$sigma)
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * a^2 * (1 - a)^2 * (trigamma(hilfs) + trigamma(hilfs2)))
      },
      "sigma" = function(y, eta, ...) {
        a <- linkinv$mu(eta$mu)
        b <- linkinv$sigma(eta$sigma)
        hilfs <- a * (1 - b) / b
        hilfs2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * (a^2 * trigamma(hilfs) + (1 - a)^2 * trigamma(hilfs2) - trigamma((1 - b) / (b))))
      }
    ),
    "mu" = function(eta, ...) {
      linkinv$mu(eta$mu)
    },
	  "d" = function(y, eta) {
       mu <- linkinv$mu(eta$mu)
       sigma <- linkinv$sigma(eta$sigma)
		   a <- mu * (1 - sigma) / (sigma)
		   b <- a * (1 - mu) / mu
		   dbeta(y, shape1 = a, shape2 = b, ncp = 0)
	  },
	  "p" = function(y, eta) {
       mu <- linkinv$mu(eta$mu)
       sigma <- linkinv$sigma(eta$sigma)
		   a <- mu * (1 - sigma) / (sigma)
		   b <- a * (1 - mu) / mu
		   pbeta(y, shape1 = a, shape2 = b, ncp = 0)
	  }
  )
  class(rval) <- "family.BayesR"
  rval
}


betazoi.BayesR <- function(links = c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...)
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
    bayesx = list(
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
      linkinv$mu(eta$mu) * (1 - (linkinv$nu(eta$nu) + linkinv$tau(eta$tau)) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau))) +
        linkinv$tau(eta$tau) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau))
    },
	"d" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		d <- ifelse(y == 0, linkinv$nu(eta$nu) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau)), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau)))
		ifelse (y==1, linkinv$tau(eta$tau) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau)), d)
	},
	"p" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		hilfs <- linkinv$nu(eta$nu) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau))
		hilfs2 <- linkinv$tau(eta$tau) / (1 + linkinv$nu(eta$nu) + linkinv$tau(eta$tau))
		cdf <- ifelse(y == 0, hilfs, hilfs + (1 - (hilfs + hilfs2)) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
		ifelse(y == 1, 1, cdf)
	}
  )
  class(rval) <- "family.BayesR"
  rval
}


betazi.BayesR <- function(links = c(mu = "logit", sigma2 = "logit", nu = "log"), ...)
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
    bayesx = list(
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
      linkinv$mu(eta$mu) * (1 - (linkinv$nu(eta$nu)) / (1 + linkinv$nu(eta$nu)))
    },
	"d" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		ifelse(y == 0, linkinv$nu(eta$nu) / (1 + linkinv$nu(eta$nu)), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + linkinv$nu(eta$nu)))
	},
	"p" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		hilfs <- linkinv$nu(eta$nu) / (1 + linkinv$nu(eta$nu))
		ifelse(y == 0, hilfs, hilfs + (1 - hilfs) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
	}
  )
  class(rval) <- "family.BayesR"
  rval
}


betaoi.BayesR <- function(links = c(mu = "logit", sigma2 = "logit", tau = "log"), ...)
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
    bayesx = list(
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
      linkinv$mu(eta$mu) * (1 - linkinv$tau(eta$tau) / (1 + linkinv$tau(eta$tau))) +
        linkinv$tau(eta$tau) / (1 + linkinv$tau(eta$tau))
    },
	"d" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		ifelse(y == 1, linkinv$tau(eta$tau) / (1 + linkinv$tau(eta$tau)), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + linkinv$tau(eta$tau)))
	},
	"p" = function(y, eta) {
		mu <- linkinv$mu(eta$mu)
		sigma <- linkinv$sigma(eta$sigma)
		a <- mu * (1 - sigma) / (sigma)
		b <- a * (1 - mu) / mu
		hilfs <- linkinv$tau(eta$tau) / (1 + linkinv$tau(eta$tau))
		ifelse(y == 1, hilfs, hilfs + (1 - hilfs) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
	}
  )
  class(rval) <- "family.BayesR"
  rval
}


binomial.BayesR <- function(link = "logit", ...)
{
  rval <- list(
    "family" = "binomial",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "logit"), ...),
    "valid.response" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      if(nlevels(x) > 2) stop("more than 2 levels in factor response!", call. = FALSE)
      TRUE
    },
    bayesx = list(
      "pi" = c(paste("binomial", link, sep = "_"), "mean")
    ),
    jagstan = list(
      "dist" = "dbern",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    ),
    "mu" = function(eta, ...) {
      linkinv$pi(eta$pi)
    },
	"d" = function(y, eta) {
		dbinom(y, size = 1, prob = linkinv$pi(eta$pi))
	},
	"p" = function(y, eta) {
		pbinom(y, size = 1, prob = linkinv$pi(eta$pi))
	}
  )
  class(rval) <- "family.BayesR"
  rval
}




gaussian.BayesR <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  links <- parse.links(links, c(mu = "identity", sigma = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = links,
    bayesx = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal2_mu", "mean"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal2_sigma", "scale"),
        "logit" = c("normal_sigma_logit", "scale")
      )
    ),
    jagstan = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(sigma = "1 / sqrt(sigma)")
    ),
    "loglik" = function(y, eta, ...) {
      sum(dnorm(y, eta$mu, linkinv$sigma(eta$sigma), log = TRUE))
    },
    "score" = list(
      "mu" = function(y, eta, ...) { drop((y - eta$mu) / (linkinv$sigma(eta$sigma)^2)) },
      "sigma" = function(y, eta, ...) { drop(-0.5 + (y - eta$mu)^2 / (linkinv$sigma(eta$sigma)^2)) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { drop(1 / (linkinv$sigma(eta$sigma)^2)) },
      "sigma" = function(y, eta, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(eta, ...) {
      eta$mu
    },
    "d" = function(y, eta) {
      dnorm(y, mean = eta$mu, sd = linkinv$sigma(eta$sigma))
    },
    "p" = function(y, eta) {
      pnorm(y, mean = eta$mu, sd = linkinv$sigma(eta$sigma))
    }
  )
  
  class(rval) <- "family.BayesR"
  rval
}


gaussian2.BayesR <- function(links = c(mu = "identity", sigma2 = "log"), ...)
{
  links <- parse.links(links, c(mu = "identity", sigma2 = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "gaussian2",
    "names" = c("mu", "sigma2"),
    "links" = links,
    bayesx = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal_mu", "mean"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma2" = switch(links["sigma2"],
        "log" = c("normal_sigma2", "scale"),
        "logit" = c("normal_sigma2_logit", "scale")
      )
    ),
    jagstan = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(sigma2 = "1 / sigma2")
    ),
    "loglik" = function(y, eta, ...) {
      sum(dnorm(y, eta$mu, sqrt(linkinv$sigma2(eta$sigma2)), log = TRUE))
    },
    "score" = list(
      "mu" = function(y, eta, ...) { drop((y - eta$mu) / linkinv$sigma2(eta$sigma2)) },
      "sigma2" = function(y, eta, ...) { drop(-0.5 + (y - eta$mu)^2 / (2 * linkinv$sigma2(eta$sigma2))) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { drop(1 / linkinv$sigma2(eta$sigma2)) },
      "sigma2" = function(y, eta, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(eta, ...) {
      eta$mu 
    },
    "d" = function(y, eta) {
      dnorm(y, mean = eta$mu, sd = sqrt(linkinv$sigma2(eta$sigma2)))
    },
    "p" = function(y, eta) {
      pnorm(y, mean = eta$mu, sd = sqrt(linkinv$sigma2(eta$sigma2)))
    }
  )
  
  class(rval) <- "family.BayesR"
  rval
}


invgaussian.BayesR <- function(links = c(mu = "log", sigma2 = "log"), ...)
{
  links <- parse.links(links, c(mu = "log", sigma2 = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "invgaussian",
    "names" = c("mu", "sigma2"),
    "links" = links,
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    bayesx = list(
      "mu"  = c("invgaussian_mu", "mean"),
      "sigma2" = c("invgaussian_sigma2", "scale")
    ),
    "mu" = function(eta, ...) {
      linkinv$mu(eta$mu) 
    },
	"d" = function (y, eta) 
	{
		mu <- linkinv$mu(eta$mu)
		sigma <- sqrt(linkinv$sigma(eta$sigma))
		exp( -0.5 * log(2 * pi) - log(sigma) - (3 / 2) * log(x) - ((x - mu)^2) / (2 * sigma^2 * (mu^2) * x))
	},
	"p" = function (y, eta) 
	{
		mu <- linkinv$mu(eta$mu)
		lambda <- 1 / sqrt(linkinv$sigma(eta$sigma))
		lq <- sqrt(lambda / y)
		qm <- y / mu
		pnorm(lq * (qm - 1)) + exp(2 * lambda / mu) * pnorm(-lq * (qm + 1))
	}
  )
  class(rval) <- "family.BayesR"
  rval
}

gamma.BayesR <- function(links = c(mu = "log", sigma = "log"), ...)
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
    bayesx = list(
      "mu" = c("gamma_mu", "mean"),
      "sigma" = c("gamma_sigma", "shape")
    ),
    jagstan = list(
      "dist" = "dgamma",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    ),
    "mu" = function(eta, ...) {
      linkinv$mu(eta$mu) 
    },
	"d" = function (y, eta, log = FALSE) {
		a <- linkinv$sigma(eta$sigma) 
		s <- linkinv$mu(eta$mu) / linkinv$sigma(eta$sigma) 
		dgamma(y, shape = a, scale = s, log = log)
	},
	"d" = function (y, eta, lower.tail = TRUE, log.p = FALSE) {
		a <- linkinv$sigma(eta$sigma) 
		s <- linkinv$mu(eta$mu) / linkinv$sigma(eta$sigma) 
		pgamma(y, shape = a, scale = s, lower.tail = lower.tail, log.p = log.p)
	}
  )
  class(rval) <- "family.BayesR"
  rval
}


lognormal.BayesR <- function(links = c(mu = "log", sigma = "log"), ...)
{
  rval <- list(
    "family" = "lognormal",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "log", sigma = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    bayesx = list(
      "mu" = c("lognormal2_mu", "mean"),
      "sigma" = c("lognormal2_sigma", "scale")
    ),
    "integrand" = function(y, eta) {
      dlnorm(y, meanlog = eta$mu, sdlog = (linkinv$sigma(eta$sigma)))^2
    },
	"mu" = function(eta, ...) {
      exp(eta$mu + 0.5 * (linkinv$sigma(eta))^2)
    },
    "d" = function(y, eta) {
      dlnorm(y, meanlog = eta$mu, sdlog = linkinv$sigma(eta$sigma))
    },
    "p" = function(y, eta) {
      plnorm(y, meanlog = eta$mu, sdlog = linkinv$sigma(eta$sigma))
    }
  )
  class(rval) <- "family.BayesR"
  rval
}


lognormal2.BayesR <- function(links = c(mu = "log", sigma2 = "log"), ...)
{
  rval <- list(
    "family" = "lognormal2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "log", sigma2 = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    bayesx = list(
      "mu" = c("lognormal_mu", "mean"),
      "sigma2" = c("lognormal_sigma2", "scale")
    ),
	"mu" = function(eta, ...) {
      exp(linkinv$mu(eta) + 0.5 * (linkinv$sigma(eta)))
    },
    "d" = function(y, eta) {
      dlnorm(y, meanlog = eta$mu, sdlog = sqrt(linkinv$sigma(eta$sigma)))
    },
    "p" = function(y, eta) {
      plnorm(y, meanlog = eta$mu, sdlog = sqrt(linkinv$sigma(eta$sigma)))
    }
  )
  class(rval) <- "family.BayesR"
  rval
}


dagum.BayesR <- function(links = c(a = "log", b = "log", p = "log"), ...)
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
    bayesx = list(
      "a" = c("dagum_a", "mean"),
      "b" = c("dagum_b", "scale"),
      "p" = c("dagum_p", "shape2")
    ),
	"mu" = function(eta, ...) {
      - (linkinv$b(eta$b)/linkinv$a(eta$a)) * (gamma(- 1 / linkinv$a(eta$a)) * gamma(linkinv$p(eta$p) + 1 / linkinv$a(eta$a))) / (gamma(linkinv$p(eta$p)))
    },
	"d" = function(y, eta) {
		ap <- linkinv$a(eta$a) * linkinv$p(eta$p)
		ap * y^(ap -1) / (b^ap * (1 + (y / linkinv$b(eta$b))^linkinv$a(eta$a))^(linkinv$p(eta$p) + 1))
	},
	"p" = function(y, eta) {
		( 1 + (y / linkinv$b(eta$b))^(-linkinv$a(eta$a)))^(-linkinv$p(eta$p))
	}
  )
  class(rval) <- "family.BayesR"
  rval
}



mvn.BayesR <- function(links = c(mu1 = "identity", mu2 = "identity",
                                 sigma1 = "log", sigma2 = "log", rho = "rhogit"), ...)
{  
  links <- parse.links(links, c(mu1 = "identity", mu2 = "identity",
                                sigma1 = "log", sigma2 = "log", rho = "rhogit"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "mvn",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho"),
    "links" = links,
    bayesx = list(
      "mu1" = c("bivnormal_mu", "mean"),
      "mu2" = c("bivnormal_mu", "mu"),
      "sigma1" = c("bivnormal_sigma", "scale1"),
      "sigma2" = c("bivnormal_sigma", "scale2"),
      "rho" = c("bivnormal_rho", "rho"),
      "order" = 5:1,
      "rm.number" = TRUE
    ),
	"mu" = function(eta, ...) {
      c(eta$mu1, eta$mu2)
    }
  )
  class(rval) <- "family.BayesR"
  rval
}


mvt.BayesR <- function(links = c(mu1 = "identity", mu2 = "identity",
                                 sigma1 = "log", sigma2 = "log", rho = "fisherz", df = "log"), ...)
{
  rval <- list(
    "family" = "mvt",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho", "df"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity",
                                   sigma1 = "log", sigma2 = "log", rho = "fisherz", df = "log"), ...),
    bayesx = list(
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
      c(linkinv$mu1(eta), linkinv$mu2(eta))
    }
  )
  class(rval) <- "family.BayesR"
  rval
}




multinomial.BayesR <- function(link = "probit", ...)
{
  rval <- list(
    "family" = "multinomial",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "probit"), ...),
    "cat" = TRUE,
    "valid.response" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      TRUE
    },
    jagstan = list(
      "dist" = "dcat",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    ),
    bayesx = list(
      "pi" = c(paste("multinom", link, sep = "_"), "mean", "meanservant")
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


## Count Data distributions

poisson.BayesR <- function(links = c(lambda = "log"), ...)
{
  links <- parse.links(links, c(lambda = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "poisson",
    "names" = c("lambda"),
    "links" = links,
    bayesx = list(
      "lambda" = c("poisson", "mean")
    ),
	"mu" = function(eta, ...) {
      linkinv$lambda(eta)
    }
  )
  class(rval) <- "family.BayesR"
  rval
}

zip.BayesR <- function(links = c(lambda = "log", pi = "logit"), ...)
{
  links <- parse.links(links, c(lambda = "log", pi = "logit"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "zip",
    "names" = c("lambda", "pi"),
    "links" = links,
    bayesx = list(
      "lambda" = c("zip_lambda", "mean"),
      "pi" = switch(links["pi"],
                    "logit" = c("zip_pi", "pi"),
                    "cloglog2" = c("zip_pi_cloglog", "pi")
      ) 
    ),
	"mu" = function(eta, ...) {
      linkinv$lambda(eta) * (1 - linkinv$pi(eta))
    }
  )
  if(rval$bayesx[[2]][[1]] == "zip_pi_cloglog")
    rval$bayesx[[1]][[1]] <- "zip_lambda_cloglog"
  class(rval) <- "family.BayesR"
  rval
}


negbin.BayesR <- function(links = c(mu = "log", delta = "log"), ...)
{
  links <- parse.links(links, c(mu = "log", delta = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "negbin",
    "names" = c("mu", "delta"),
    "links" = links,
    bayesx = list(
      "mu" = c("negbin_mu", "mean"),
      "delta" = c("negbin_delta", "delta")
    ),
	"mu" = function(eta, ...) {
      linkinv$mu(eta)
    }
  )
  class(rval) <- "family.BayesR"
  rval
}


zinb.BayesR <- function(links = c(mu = "log", "pi" = "logit", delta = "log"), ...)
{
  links <- parse.links(links, c(mu = "log", "pi" = "logit", delta = "log"), ...)
  linkinv <- list()
  for(j in names(links))
    linkinv[[j]] <- make.link2(links[[j]])$linkinv
  
  rval <- list(
    "family" = "zinb",
    "names" = c("mu", "pi", "delta"),
    "links" = links,
    bayesx = list(
      "mu" = c("zinb_mu", "mean"),
      "pi" = c("zinb_pi", "pi"),
      "delta" = c("zinb_delta", "delta")
    ),
	"mu" = function(eta, ...) {
      linkinv$mu(eta) * (1 - linkinv$pi(eta))
    } 
  )
  class(rval) <- "family.BayesR"
  rval
}

## http://stats.stackexchange.com/questions/17672/quantile-regression-in-jags
quant.BayesR <- function(links = c(mu = "identity", sigma = "log"), prob = 0.5, ...)
{
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    bayesx = list(
      "mu" = c("normal_mu", "mean"),
      "sigma" = c("normal_sigma2", "scale")
    ),
    jagstan = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(
        mu = "(1 - 2 * prop) / (prop * (1 - prop)) * w[i] + mu",
        sigma = "(prop * (1 - prop) * (1 / sigma)) / (2 * w[i])"
      ),
      "addparam" = list("w[i] ~ dexp(1 / sigma[i])"),
      "addvalues" = list("prop" = prob)
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


### Ordered logit.
### http://staff.washington.edu/lorenc2/bayesian/ologit.R


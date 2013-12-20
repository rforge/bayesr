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
  } else if(link == "cloglog2") {
  
	
	cloglog2 <- function() {
        linkfun <- function(mu) log(-log(mu))
        linkinv <- function(eta) pmax(pmin(1-expm1(-exp(eta)), 
             .Machine$double.eps), .Machine$double.eps)
        mu.eta <- function(eta) {
            eta <- pmin(eta, 700)
            pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
        }
	}
	rval <- cloglog2()
  } else {
    stop("unknown link function!")
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
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in [0, 1]!", call. = FALSE)
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
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


betazoi.BayesR <- function(links = c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...)
{
  rval <- list(
    "family" = "betainflated",
    "names" = c("mu", "sigma2", "nu", "tau"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in [0, 1]!", call. = FALSE)
      ok
    },
    "order" = 4:1,
    bayesx = list(
      "mu" = c("betainf_mu", "location"),
      "sigma2" = c("betainf_sigma2", "scale"),
      "nu" = c("betainf_nu", "shape"),
      "tau" = c("betainf_tau", "mean")
    )
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
    )
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
        "identity" = c("normal_mu", "mean"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal_sigma2", "scale"),
        "logit" = c("normal_sigma2_logit", "scale")
      )
    ),
    jagstan = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(sigma = "1 / sigma")
    ),
    "loglik" = function(y, eta, ...) {
      sum(dnorm(y, eta$mu, sqrt(linkinv$sigma(eta$sigma)), log = TRUE))
    },
    "score" = list(
      "mu" = function(y, eta, ...) { drop((y - eta$mu) / linkinv$sigma(eta$sigma)) },
      "sigma" = function(y, eta, ...) { drop(-0.5 + (y - eta$mu)^2 / (2 * linkinv$sigma(eta$sigma))) }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) { drop(1 / linkinv$sigma(eta$sigma)) },
      "sigma" = function(y, eta, ...) { rep(0.5, length(y)) }
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


mvn.BayesR <- function(links = c(mu1 = "identity", mu2 = "identity",
  sigma1 = "log", sigma2 = "log", rho = "identity"), ...)
{
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity",
       sigma1 = "log", sigma2 = "log", rho = "identity"), ...),
    bayesx = list(
      "mu1" = c("bivnormal_mu", "mean"),
      "mu2" = c("bivnormal_mu", "mu"),
      "sigma1" = c("bivnormal_sigma", "scale"),
      "sigma2" = c("bivnormal_sigma", "scale"),
      "rho" = c("bivnormal_rho", "rho"),
      "order" = 5:1,
      "lhs" = c("sigma1" = "sigma1", "sigma2" = "sigma2", "rho" = "rho")
    )
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
    bayesx = list(
      "mu" = c("gamma_mu", "mean"),
      "sigma" = c("gamma_sigma", "shape")
    ),
    jagstan = list(
      "dist" = "dgamma",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


lognormal.BayesR <- function(links = c(mu = "log", sigma = "log"), ...)
{
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "log", sigma = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    bayesx = list(
      "mu" = c("lognormal_mu", "mean"),
      "sigma" = c("lognormal_sigma2", "scale")
    ))
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
	  
    )
  )
  if(rval$bayesx[[2]][[1]] == "zip_pi_cloglog")
		rval$bayesx[[1]][[1]] <- "zip_lambda_cloglog"
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


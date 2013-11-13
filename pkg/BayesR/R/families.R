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


make.link2 <- function(link)
{
  if(link %in% c("logit", "probit", "cauchit", "cloglog", "identity",
    "log", "sqrt", "1/mu^2", "inverse")) {
    rval <- make.link(link)
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
beta.BayesR <- function(links = c(mu = "logit", phi = "log"), ...)
{
  rval <- list(
    "family" = "beta",
    "names" = c("mu", "phi"),
    "links" =  parse.links(links, c(mu = "logit", phi = "log"), ...),
    "valid.response" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in [0, 1]!", call. = FALSE)
      ok
    },
    jags = list(
      "dist" = "dbeta",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(
        alpha = "mu * phi",
        beta = "(1 - mu) * phi"
      )
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
    jags = list(
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
  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    bayesx = list(
      "mu" = c("normal_mu", "mean"),
      "sigma" = c("normal_sigma2", "scale")
    ),
    jags = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"

  rval
}


multinomial.BayesR <- function(link = "logit", ...)
{
  rval <- list(
    "family" = "multinomial",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "logit"), ...),
    "cat" = TRUE,
    "valid.response" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      TRUE
    },
    jags = list(
      "dist" = "dcat",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
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
    jags = list(
      "dist" = "dnorm",
      "eta" = JAGSeta,
      "model" = JAGSmodel,
      "reparam" = c(
        m1 = "(1 - 2 * prop) / (prop * (1 - prop)) * w[i] + mu",
        s1 = "(prop * (1 - prop) * sigma) / (2 * w[i])"
      ),
      "addparam" = list("w[i] ~ dexp(sigma[i])"),
      "addvalues" = list("prop" = prob)
    )
  )
  class(rval) <- "family.BayesR"

  rval
}


### Ordered logit.
### http://staff.washington.edu/lorenc2/bayesian/ologit.R

#################
### (b) BayesX ##
#################
#binomial.BayesX <- function(link = "logit")
#{
#  rval <- list(
#    "family" = "binomial",
#    "k" = 1,
#    "mu.link" = link,
#    "names" = "binomial",
#    "binomial" = c(paste("binomial", link, sep = "_"), "mean"),
#    "all" = TRUE,
#    "h" = "gaussian_re",
#    "factor" = TRUE
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}

#gaussian.BayesX <- function(...)
#{
#  rval <- list(
#    "family" = "gaussian",
#    "k" = 1,
#    "mu.link" = "identity",
#    "names" = "gaussian",
#    "gaussian" = c("gaussian", "mean"),
#    "all" = TRUE,
#    "h" = "gaussian_re"
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}

#normal.BayesX <- function(...)
#{
#  rval <- list(
#    "family" = "normal",
#    "k" = 2,
#    "mu.link" = "identity",
#    "sigma2.link" = "log",
#    "names" = c("mu", "sigma2"),
#    "mu" = c("normal_mu", "mean"),
#    "sigma2" = c("normal_sigma2", "scale"),
#    "all" = TRUE,
#    "h" = "gaussian_re"
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}

#lognormal.BayesX <- function(...)
#{
#  rval <- list(
#    "family" = "lognormal",
#    "k" = 2,
#    "mu.link" = "log",
#    "sigma2.link" = "log",
#    "names" = c("mu", "sigma2"),
#    "mu" = c("lognormal_mu", "mean"),
#    "sigma2" = c("lognormal_sigma2", "scale"),
#    "all" = TRUE,
#    "h" = "gaussian_re",
#    "valid.response" = function(x) {
#      if(is.factor(x)) return(FALSE)
#      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
#      ok
#    }
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}

#beta.BayesX <- function(...)
#{
#  rval <- list(
#    "family" = "beta",
#    "k" = 2,
#    "mu.link" = "logit",
#    "sigma2.link" = "logit",
#    "names" = c("mu", "sigma2"),
#    "mu" = c("beta_mu", "mean"),
#    "sigma2" = c("beta_sigma2", "scale"),
#    "all" = TRUE,
#    "h" = "gaussian_re",
#    "valid.response" = function(x) {
#      if(is.factor(x)) return(FALSE)
#      if(ok <- !all(x > 0 & x < 1)) stop("response values not in [0, 1]!", call. = FALSE)
#      ok
#    }
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}

#betainflated.BayesX <- function(...)
#{
#  rval <- list(
#    "family" = "betainflated",
#    "k" = 4,
#    "mu.link" = "logit",
#    "sigma2.link" = "logit",
#    "nu.link" = "log",
#    "tau.link" = "log",
#    "names" = c("mu", "sigma2", "nu", "tau"),
#    "mu" = c("betainf_mu", "location"),
#    "sigma2" = c("betainf_sigma2", "scale"),
#    "nu" = c("betainf_nu", "shape"),
#    "tau" = c("betainf_tau", "mean"),
#    "all" = TRUE,
#    "h" = "gaussian_re",
#    "order" = 4:1,
#    "valid.response" = function(x) {
#      if(is.factor(x)) return(FALSE)
#      if(ok <- !all(x > 0 & x < 1)) stop("response values not in [0, 1]!", call. = FALSE)
#      ok
#    }
#  )
#  class(rval) <- "family.BayesR"
#  rval
#}


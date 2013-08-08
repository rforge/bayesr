###################
## (1) Families. ##
###################
## Print method.
print.family.BayesR <- function(x)
{
  cat("Family:", x$family, "\n")
  links <- NULL
  for(j in grep("link", names(x), value = TRUE))
    links <- c(links, paste(gsub(".link", "", j, fixed = TRUE), " = \"", x[[j]], "\"", sep = ""))
  links <- paste(links, collapse = ", ")
  cat("Links:", links, sep = " ")
  cat("\n")
}


##############
## (a) JAGS ##
##############
gaussian.JAGS <- function(mu.link = "identity", sigma.link = "log")
{
  rval <- list(
    "family" = "gaussian",
    "k" = 2,
    "mu.link" = mu.link,
    "sigma.link" = sigma.link,
    "names" = c("mu", "sigma"),
    "dist" = "dnorm",
    "default.prior" = c("mu ~ dnorm(0, 1.0E-6)", "sigma ~ dgamma(1.0E-6, 1.0E-6)"),
    "trans" <- list("mu" = function(x) { x }, "sigma" = function(x) { 1 / exp(x) }),
    "eta" = JAGSeta,
    "model" = JAGSmodel
  )
  class(rval) <- "family.BayesR"
  rval
}

beta.JAGS <- function(shape1.link = "log", shape2.link = "log")
{
  rval <- list(
    "family" = "beta",
    "k" = 2,
    "shape1.link" = shape1.link,
    "shape2.link" = shape2.link,
    "names" = c("shape1", "shape2"),
    "valideta" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("values not in [0, 1] using beta.BayesR!", call. = FALSE)
      ok
    },
    "dist" = "dbeta",
    "default.prior" = c("shape1 ~ dgamma(1.0E-6, 1.0E-6)", "shape2 ~ dgamma(1.0E-6, 1.0E-6)"),
    "eta" = JAGSeta,
    "model" = JAGSmodel
  )
  class(rval) <- "family.BayesR"
  rval
}

binomial.JAGS <- function(link = "logit")
{
  rval <- list(
    "family" = "binomial",
    "k" = 1,
    "link" = link,
    "names" = "pi",
    "valideta" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      if(nlevels(x) > 2) stop("more than 2 levels in factor response!", call. = FALSE)
      TRUE
    },
    "dist" = "dbern",
    "default.prior" = "pi ~ dnorm(0, 1.0E-6)",
    "eta" = JAGSeta,
    "model" = JAGSmodel
  )
  class(rval) <- "family.BayesR"
  rval
}

multinomial.JAGS <- function(link = "logit")
{
  rval <- list(
    "family" = "multinomial",
    "k" = Inf,
    "link" = link,
    "names" = "pi",
    "valideta" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      TRUE
    },
    "cat" = TRUE,
    "dist" = "dcat",
    "default.prior" = "pi ~ dnorm(0, 1.0E-6)",
    "eta" = JAGSeta,
    "model" = JAGSmodel
  )
  class(rval) <- "family.BayesR"
  rval
}

## Ordered logit.
## http://staff.washington.edu/lorenc2/bayesian/ologit.R


################
## (b) BayesX ##
################
gaussian.BayesX <- function(...)
{
  rval <- list(
    "family" = "gaussian",
    "k" = 1,
    "mu.link" = "identity",
    "names" = "mu",
    "mu" = c("gaussian", "mean"),
    "all" = TRUE,
    "h" = "gaussian_re"
  )
  class(rval) <- "family.BayesR"
  rval
}

lognormal.BayesX <- function(...)
{
  rval <- list(
    "family" = "lognormal",
    "k" = 2,
    "mu.link" = "log",
    "sigma2.link" = "log",
    "names" = c("mu", "sigma2"),
    "mu" = c("lognormal_mu", "mean"),
    "sigma2" = c("lognormal_sigma2", "scale"),
    "all" = TRUE,
    "h" = "gaussian_re"
  )
  class(rval) <- "family.BayesR"
  rval
}

beta.BayesX <- function(...)
{
  rval <- list(
    "family" = "beta",
    "k" = 2,
    "mu.link" = "logit",
    "sigma2.link" = "logit",
    "names" = c("mu", "sigma2"),
    "mu" = c("beta_mu", "mean"),
    "sigma2" = c("beta_sigma2", "scale"),
    "all" = TRUE,
    "h" = "gaussian_re"
  )
  class(rval) <- "family.BayesR"
  rval
}

betainflated.BayesX <- function(...)
{
  rval <- list(
    "family" = "betainflated",
    "k" = 4,
    "mu.link" = "logit",
    "sigma2.link" = "logit",
    "nu.link" = "log",
    "tau.link" = "log",
    "names" = c("mu", "sigma2", "nu", "tau"),
    "mu" = c("betainf_mu", "location"),
    "sigma2" = c("betainf_sigma2", "scale"),
    "nu" = c("betainf_nu", "shape"),
    "tau" = c("betainf_tau", "mean"),
    "all" = TRUE,
    "h" = "gaussian_re",
    "order" = 1:4
  )
  class(rval) <- "family.BayesR"
  rval
}


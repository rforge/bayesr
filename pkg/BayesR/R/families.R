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


##  Continuous responses.
gaussian.BayesR <- function(mu.link = "identity", sigma.link = "log")
{
  rval <- list(
    "family" = "gaussian",
    "k" = 2,
    "mu.link" = mu.link,
    "sigma.link" = sigma.link,
    "names" = c("mu", "sigma"),
    "JAGS" = list(
      "dist" = "dnorm",
      "default.prior" = c("mu ~ dnorm(0, 1.0E-6)", "sigma ~ dgamma(1.0E-6, 1.0E-6)"),
      "trans" <- list("mu" = function(x) { x }, "sigma" = function(x) { 1 / exp(x) }),
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"
  rval
}

beta.BayesR <- function(shape1.link = "log", shape2.link = "log")
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
    "JAGS" = list(
      "dist" = "dbeta",
      "default.prior" = c("shape1 ~ dgamma(1.0E-6, 1.0E-6)", "shape2 ~ dgamma(1.0E-6, 1.0E-6)"),
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


## Categorical responses.
binomial.BayesR <- function(link = "logit")
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
    "JAGS" = list(
      "dist" = "dbern",
      "default.prior" = "pi ~ dnorm(0, 1.0E-6)",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"
  rval
}

multinomial.BayesR <- function(link = "logit")
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
    "JAGS" = list(
      "dist" = "dmulti",
      "default.prior" = "pi ~ dnorm(0, 1.0E-6)",
      "eta" = JAGSeta,
      "model" = JAGSmodel
    )
  )
  class(rval) <- "family.BayesR"
  rval
}


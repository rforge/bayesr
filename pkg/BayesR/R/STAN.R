############################
## (1) STAN installation. ##
############################
install.stan <- function() {
  require("Rcpp")
  require("inline")
  options("repos" = c(getOption("repos"), rstan = "http://wiki.rstan-repo.googlecode.com/git/"))
  install.packages("rstan", type = "source")
}


################################
## (2) STAN helper functions. ##
################################
list2stan <- function(x)
{
  if(is.null(nx <- names(x))) stop("the data list must be a named list!")

  dtype <- function(type) {
    switch(type,
      'integer' = 'int',
      'numeric' = 'vector',
      'logical' = NA
    )
  }

  rval <- NULL
  for(j in seq_along(x)) {
    rval <- if(!is.matrix(x[[j]])) {
      tx <- dtype(class(x[[j]]))
      if(tx == 'vector') {
        c(rval, paste(tx, '[', length(x[[j]]), '] ', nx[j], ';', sep = ''))
      } else c(rval, paste(tx, ' ', nx[j], '[', length(x[[j]]), ']', ';', sep = ''))
    } else c(rval, paste('matrix[', nrow(x[[j]]), ',', ncol(x[[j]]), '] ', nx[j], ';', sep = ''))
  }

  rval
}


## Sets up data and model code for fitting with STAN.
jags2stan <- function(x)
{
  x <- setupJAGS(x)

  jdata <- list2stan(x$data)
  jpar <- list2stan(x$inits)
  
  data <- c(
    'data {',
    paste('  ', jdata, sep = ''),
    '}'
  )

  parameters <- c(
    'parameters {',
    paste('  ', jpar, sep = ''),
    '}'
  )

  dist2stan <- function(x) {
    x <- gsub('dbeta(', 'beta(', x, fixed = TRUE)
    x <- gsub('dbern(', 'bernoulli(', x, fixed = TRUE)
    x <- gsub('dnorm(', 'normal(', x, fixed = TRUE)
    x <- gsub('dgamma(', 'gamma(', x, fixed = TRUE)
    x
  }

  model <- dist2stan(c(data, parameters, x$model))

  x$model <- model

  x
}


########################################
## (3) Interface to the STAN sampler. ##
########################################
samplerSTAN <- function(x, tdir = NULL,
  n.chains = 1, n.adapt = 100,
  n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = TRUE, set.inits = FALSE, ...)
{
  require("rstan")

  ## Temporary directory handling.
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  } else tdir <- path.expand(tdir)
  if(!file.exists(tdir))
    dir.create(tdir)

  ## Set the seed of the random number generator.
  if(is.null(seed))
    seed <- floor(runif(n.chains) * .Machine$integer.max)

  ## Write the model code.
  writeLines(paste(x$model, collapse = "\n"), mfile <- file.path(tdir, "model.txt"))

  if(verbose) writeLines(x$model)

  smodel <- stan(mfile, fit = NA, data = x$data, chains = n.chains, iter = n.iter,
    thin = thin, warmup = burnin, seed = seed, verbose = verbose, ...)

  smodel
}

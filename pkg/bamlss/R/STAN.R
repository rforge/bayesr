################################
## (2) STAN helper functions. ##
################################
## Sets up data and model code for fitting with STAN.
bugs2stan <- function(x)
{
  STAN_data <- function(x)
  {
    if(is.null(nx <- names(x))) stop("the data list must be a named list!")

    dtype <- function(type) {
      switch(type,
        'integer' = 'int',
        'numeric' = 'real',
        'logical' = NA
      )
    }

    rval <- NULL
    for(j in seq_along(x)) {
      rval <- if(!is.matrix(x[[j]])) {
        tx <- dtype(class(x[[j]]))
        if(tx == 'vector') {
          c(rval, paste(tx, '[', length(x[[j]]), '] ', nx[j], ';', sep = ''))
        } else {
          if(length(x[[j]]) > 1)
            c(rval, paste(tx, ' ', nx[j], '[', length(x[[j]]), ']', ';', sep = ''))
          else
            c(rval, paste(tx, ' ', nx[j], ';', sep = ''))
        }
      } else c(rval, paste('matrix[', nrow(x[[j]]), ',', ncol(x[[j]]), '] ', nx[j], ';', sep = ''))
    }

    rval
  }

  STAN_model_data <- function(x, n) {
    d <- grep("<-", x, fixed = TRUE, value = TRUE)
    d <- sapply(strsplit(d, "<-", fixed = TRUE), function(x) { x[1] })
    if(length(d)) {
      d <- gsub("[i]", "", d, fixed = TRUE)
      d <- gsub("\\s", "", d)
      d <- paste("  real ", d, "[", n, "]", sep = "")
    } else d <- NULL
    d
  }

  attr(x, "is.stan") <- TRUE
  x <- setupJAGS(x)
  
  data <- c(
    'data {',
    paste('  ', STAN_data(x$data), sep = ''),
    '}'
  )

  parameters <- c(
    'parameters {',
    paste('  ', STAN_data(x$inits), sep = ''),
    '}'
  )

  ## See http://www.jrnold.me/blog/jagsopenbugs-to-stan-distributions.html
  dist2stan <- function(x) {
    x <- gsub('dbeta(', 'beta(', x, fixed = TRUE)
    x <- gsub('dcat(', 'categorical(', x, fixed = TRUE)
    x <- gsub('dbern(', 'bernoulli(', x, fixed = TRUE)
    if(length(i <- grep('dnorm(', x, fixed = TRUE))) {
      for(j in i) {
        xs <- strsplit(x[j], 'dnorm(', fixed = TRUE)[[1]]
        ri <- strsplit(xs[2], ',', fixed = TRUE)[[1]]
        ri[2] <- gsub(';', '', gsub(')', '', ri[2], fixed = TRUE), fixed = TRUE)
        ri[2] <- gsub(' ', '', ri[2])
        ri[2] <- paste('pow(', ri[2], ', -0.5)', sep = '')
        ri <- paste('normal(', ri[1], ', ', ri[2], ');', sep = '')
        xs <- paste(xs[1], ri, sep = '')
        x[j] <- xs
      }
    }
    x <- gsub('dgamma(', 'gamma(', x, fixed = TRUE)
    x <- gsub('dpois(', 'poisson(', x, fixed = TRUE)
    x <- gsub('dunif(', 'uniform(', x, fixed = TRUE)
    x
  }

  i <- grep("model {", x$model, fixed = TRUE)
  for(j in (i + 1):length(x$model)) {
    if(grepl("}", x$model[j], fixed = TRUE))
      break
  }

  i2 <- i + 2; j2 <- j - 1
  x$model[i2:j2] <- rev(x$model[i2:j2])
  x$model <- c(x$model[i], x$model[(j + 1):(length(x$model) - 1)],
    x$model[(i + 1):j], x$model[length(x$model)])

  x$model <- c(x$model[1:i], STAN_model_data(x$model, x$data$n), x$model[(i + 1):length(x$model)])
  i <- grepl("{", x$model, fixed = TRUE) |  grepl("}", x$model, fixed = TRUE)
  x$model[!i] <- paste(x$model[!i], ";", sep = "") 

  model <- dist2stan(c(data, parameters, x$model))
  model <- gsub('i in 1:n', paste('i in 1:', x$data$n, sep = ''), model, fixed = TRUE)

  x$model <- model

  x
}


########################################
## (3) Interface to the STAN sampler. ##
########################################
STAN <- function(x, tdir = NULL,
  n.chains = 1, n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = FALSE, show.model = TRUE, ...)
{
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
  writeLines(paste(x$model, collapse = "\n"), mfile <- file.path(tdir, "STANmodel.txt"))

  if(show.model) writeLines(paste(x$model, collapse = "\n"))

  if(verbose) writeLines(x$model)

  smodel <- stan(mfile, fit = NA, data = x$data, chains = n.chains, iter = n.iter,
    thin = thin, warmup = burnin, seed = seed, verbose = verbose, ...)

  samples <- slot(smodel, "sim")$samples
  for(j in seq_along(samples))
    samples[[j]] <- as.mcmc(do.call("cbind", samples[[j]]))
  samples <- as.mcmc.list(samples)
  samples <- window(samples, start = ceiling(burnin / thin))

  samples
}

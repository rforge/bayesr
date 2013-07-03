#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, chains = 1, combine = TRUE, n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, dir = NULL, model.name = NULL, ...)
{
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    on.exit(unlink(dir))
  } else dir <- path.expand(dir)
  if(is.null(model.name))
    model.name <- 'bayesr'
  prg.name <- paste(model.name, 'prg', sep = '.')
  data.name <- if(!is.null(data)) {
    deparse(substitute(data), backtick = TRUE, width.cutoff = 500)
  } else "d"

  cat('%% BayesX program created by BayesR: ', as.character(Sys.time()), '\n',
    file = file.path(dir, prg.name), sep = '')
  cat('%% usefile ', file.path(dir, prg.name), '\n\n', 
    file = file.path(dir, prg.name), append = TRUE, sep = '')

  transform <- function(x) { transformBayesX(x, dir = dir, prg.name = prg.name, ...) }
  setup <- function(x) {
    setupBayesX(x, n.iter = n.iter, thin = thin, burnin = burnin,
      seed = seed, model.name = model.name, data.name = data.name,
      dir = dir, prg.name = prg.name, ...)
  }
  sampler <- function(x) { samplerBayesX(x, ...) }
  results <- function(x, ...) { resultsBayesX(x, ...) }

  bayesr(formula, family = family, data = data, knots = knots,
    weights = weights, subset = subset, offset = offset, na.action = na.action,
    contrasts = contrasts, parse.input = parse.input.bayesr, transform = transform,
    setup = setup, sampler = sampler, results = results,
    cores = cores, combine = combine, ...)
}


####################################
## (2) BayesX specific functions. ##
####################################
transformBayesX <- function(x, ...)
{
  args <- list(...)
  dir <- args$dir
  prg.name <- args$prg.name
  if(inherits(x, "list") & !any(c("smooth", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    if(is.null(nx)) nx <- 1:length(x)
    for(j in nx)
      x[[j]] <- transformBayesX(x[[j]], dir = dir, prg = prg.name, ...)
  } else {
    x <- randomize(x)
    if(length(x$smooth)) stop("arbitrary smooths not supported yet!")
    if(length(x$sx.smooth)) {
      for(j in seq_along(x$sx.smooth)) {
        x$sx.smooth[[j]] <- bayesx.construct(object = x$sx.smooth[[j]], dir = dir,
          prg = prg.name, data = x$mf)
      }
    }
  }

  x
}


BayesX.control <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "full", ...)
{
  if(is.null(seed))
    seed <- round(runif(1L) * .Machine$integer.max)
  stopifnot(burnin < n.iter)
  cvals <- list("iterations" = n.iter, "burnin" = burnin, "step" = thin, "setseed" = seed,
    "predict" = predict)
  cvals
}

setupBayesX <- function(x, control = BayesX.control(...), ...)
{
  x$call <- NULL

  args <- list(...)
  model.name <- args$model.name
  data.name <- args$data.name
  dir <- args$dir
dir <- path.expand("~/tmp")
  prg.name <- args$prg.name

  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  prg <- paste('logopen using', file.path(dir, paste(model.name, 'log\n', sep = '.')))

  BayesX_data <- function(x) {
    X <- NULL
    if(!is.null(x$X)) {
      if(ncol(x$X) > 0)
        X <- as.data.frame(x$X)
    }
    if(k2 <- length(x$sterms)) {
      X <- if(!is.null(X)) {
        cbind(X, as.data.frame(x$mf[, x$sterms]))
      } else as.data.frame(x$mf[, x$sterms])
      k1 <- ncol(X)
      colnames(X)[(k1 - k2 + 1):k1] <- x$sterms
    }
    X <- X[, !grepl("(Intercept)", colnames(X), fixed = TRUE), drop = FALSE]
    X[[x$response]] <- x$mf[[x$response]]
    X <- as.data.frame(X)
    for(j in 1:ncol(X)) {
      if(is.factor(X[, j])) {
        X[, j] <- as.integer(X[, j])
        X <- X[order(X[, j]), ]
      }
    }
    return(unique(X))
  }

  get_family <- function(x) {
    if(is.null(x$family))
      for(j in x) {
        family <- get_family(j)
    } else family <- x$family
    family
  }

  family <- get_family(x)
  family <- if(is.function(family)) family() else family
  stopifnot(!is.null(family$BayesX))
  names(x) <- rep(family$names, length.out = n)
  nx <- names(x)

  count <- 1
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        dpath <- file.path(dir, paste(dname <- paste(data.name, count, sep = ''), "raw", sep = '.'))
        d <- BayesX_data(x[[j]][[i]])
        x[[j]][[i]]$hlevel <- i
        if(!is.null(d)) {
          x[[j]][[i]]$dname <- dname
          write.table(d, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
          prg <- c(prg,
            paste('dataset ', data.name, count, sep = ''),
            paste(data.name, count, '.infile using ', dpath, sep = '')
          )
          count <- count + 1
        }
      }
    } else {
      dpath <- file.path(dir, paste(dname <- paste(data.name, count, sep = ''), "raw", sep = '.'))
      d <- BayesX_data(x[[j]])
      if(!is.null(d)) {
        x[[j]]$dname <- dname
        write.table(d, file = dpath, quote = FALSE, row.names = FALSE, col.names = TRUE)
        prg <- c(prg,
          paste('dataset ', data.name, count, sep = ''),
          paste(data.name, count, '.infile using', dpath, sep = '')
        )
        count <- count + 1
      }
    }
  }

  prg <- c(prg, paste('\nmcmcreg', model.name))
  prg <- c(prg, paste(model.name, '.outfile = ', file.path(dir, model.name), '\n', sep = ''))

  count <- 1
  for(j in n:1) {
    if(!"fake.formula" %in% names(x[[j]])) {
      for(i in length(x[[j]]):1) {
        eqn <- paste(model.name, '.hregress ', x[[j]][[i]]$response, sep = '')
        et <- x[[j]][[i]]$pterms
        if(attr(terms(x[[j]][[i]]$formula), "intercept"))
          et <- c("const", et)
        if(length(x[[j]][[i]]$sx.smooth))
          et <- c(et, x[[j]][[i]]$sx.smooth)
        if(length(et))
          eqn <- paste(eqn, '=', paste(et, collapse = ' + '))
        if(count < 2) {
          eqn <- paste(eqn, ", ", paste(names(control), "=", unlist(control),
            sep = "", collapse = " "), sep = "")
        } else eqn <- paste(eqn, ", ", sep = "")
        if(!is.null(x[[j]][[i]]$hlevel)) {
          eqn <- paste(eqn, " hlevel=", x[[j]][[i]]$hlevel, sep = "")
          if(x[[j]][[i]]$hlevel > 1)
            eqn <- paste(eqn, " family=", family$BayesX$h, sep = "")
        }
        if(i < 2)
          eqn <- paste(eqn, " family=", family$BayesX[[nx[j]]][1], sep = "")
        eqn <- paste(eqn, " equationtype=", family$BayesX[[nx[j]]][2], sep = "")
        if(!is.null(x[[j]][[i]]$dname))
          eqn <- paste(eqn, " using ", x[[j]][[i]]$dname, sep = "")
        prg <- c(prg, eqn)
        count <- count + 1
      }
    } else {
      eqn <- paste(model.name, '.hregress ', x[[j]]$response, sep = '')
      et <- x[[j]]$pterms
      if(attr(terms(x[[j]]$formula), "intercept"))
        et <- c("const", et)
      if(length(x[[j]]$sx.smooth))
        et <- c(et, x[[j]]$sx.smooth)
      if(length(et))
        eqn <- paste(eqn, '=', paste(et, collapse = ' + '))
      if(!is.null(x[[j]]$dname))
        eqn <- paste(eqn, " using ", x[[j]]$dname, sep = "")
      prg <- c(prg, eqn)
      count <- count + 1
    }
  }

  prg <- gsub("(random", "(hrandom", prg, fixed = TRUE)
  for(i in 1:3)
    prg <- gsub(paste("(psplinerw", i, sep = ""), "(pspline", prg, fixed = TRUE)

writeLines(prg)

  cat(paste(prg, collapse = '\n'), file = file.path(dir, prg.name), append = TRUE)
}


samplerBayesX <- function(x, ...) { x }
resultsBayesX <- function(x, ...) { x }



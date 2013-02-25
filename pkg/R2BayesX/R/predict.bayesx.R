predict.bayesx <- function(object, newdata = NULL, ...) {
  if(!is.null(newdata))
    stop("predicting for newdata is currently not implemented!")
  return(fitted.bayesx(object, ...))
}

.predict.bayesx <- function(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = mean, ...)
{
  object <- get.model(object, model)
  k <- length(object)
  types <- NULL
  enames <- list()
  for(j in 1:k) {
    types <- c(types, object[[j]]$model.fit$method)
    enames[[j]] <- names(object[[j]]$effects)
  }
  if(length(diff) & !all(diff(sapply(enames, length)) == 0))
    stop("the number of terms in the models is not identical, cannot compute prediction!")
  enames <- data.frame(enames)
  if(!all(apply(enames, 1, function(x) length(unique(x))) == 1))
    stop("different terms in the supplied models, cannot compute prediction!")
  enames <- names(object[[1L]]$effects)
  term <- if(!is.null(enames)) {
    if(is.null(term)) enames else {
      if(is.character(term)) {
        enames[pmatch(gsub("[[:space:]]", "", term), enames)]
      } else enames[term]
    }
  } else NULL
  term <- term[!is.na(term)]
  if(!length(term)) term <- NULL
  types <- tolower(types)
  if(is.null(types) | all(is.na(types)))
    stop("unknown estimation method used, cannot predict!")
  if(length(ut <- unique(types)) > 1)
    stop("cannot combine predictions of models using different estimation methods!")
  rval <- NULL
  if(grepl("mcmc", ut)) {
    if(missing(newdata))
      newdata <- model.frame(object[[1L]])
    if(is.character(newdata)) {
      if(file.exists(newdata <- path.expand(newdata)))
        newdata <- read.table(newdata, header = TRUE, ...)
    }
    if(is.matrix(newdata) || is.list(newdata))
      newdata <- as.data.frame(newdata)
    nn <- names(newdata)
    m.samples <- m.designs <- list()
    for(j in 1:k) {
      if(!is.null(term)) {
        for(i in term) {
          specs <- attr(object[[j]]$effects[[i]], "specs")
          if(inherits(object[[j]]$effects[[i]], "geo.bayesx"))
            specs$term <- specs$term[1]
          if(!all(specs$term %in% nn))
            stop(paste("cannot find variables", specs$term, "in newdata!"))
          if(is.null(specs$is.factor)) specs$is.factor <- FALSE
          tmp <- samples(object[[j]], term = i)
          tmp <- tmp[[1]]$Coef
          if(is.null(dim(tmp))) {
            tmp <- matrix(tmp, ncol = 1)
          }
          m.samples[[i]] <- rbind(m.samples[[i]], tmp)
          if(inherits(object[[j]]$effects[[i]], "linear.bayesx")) {
            if(specs$is.factor & !is.character(newdata)) {
              hi <- if(!is.null(object[[j]]$fixed.effects)) {
                any(grepl("(Intercept)", rownames(object[[j]]$fixed.effects), fixed = TRUE))
              } else FALSE
              nl <- nlevels(newdata[[i]]) - if(hi) 1 else 0
              if(nl != ncol(m.samples[[i]]))
                stop(paste("levels of factor variable", i, "in newdata not identical to model levels!"))
              f <- as.formula(paste("~", if(hi) "1" else "-1", "+", i))
              if(j < 2) {
                m.designs[[i]] <- model.matrix(f, data = newdata)
                if(hi) m.designs[[i]] <- m.designs[[i]][, -1]
              }
            } else {
              if(j < 2)
                m.designs[[i]] <- newdata[[i]]
            }
          } else {
            if(j < 2)
              m.designs[[i]] <- Predict.matrix.bayesx(object[[j]]$effects[[i]], newdata)
          }
          attr(m.samples[[i]], "is.factor") <- specs$is.factor
        }
      }
      if(intercept) {
        sami <- attr(object[[j]]$fixed.effects, "sample")
        if(!is.null(sami)) {
          sami <- sami[, grep("(Intercept)", colnames(sami), fixed = TRUE)]
          if(length(sami)) {
            m.samples$Intercept <- sami
            i <- if(is.null(term)) 1 else {
              length(m.designs) + 1
            }
            m.designs[[i]] <- matrix(1, nrow = nrow(newdata), ncol = 1)
          }
        }
      }
    }
    if(length(m.samples)) {
      warn <- getOption("warn")
      options("warn" = -1)
      m.samples <- as.data.frame(m.samples)
      m.designs <- as.data.frame(m.designs)
      options("warn" = warn)
      get.mu <- function(X, b) {
        as.matrix(X) %*% as.numeric(b)
      }
      rval <- apply(m.samples, 1, function(x) { get.mu(m.designs, x) })
      rval <- apply(rval, 1, FUN)
      if(!is.null(dim(rval))) {
        if(nrow(rval) != nrow(newdata))
          rval <- as.data.frame(t(rval))
      }
    } else stop("no model terms selected for prediction!")
  } else {
    if(!missing(newdata))
      stop(paste("out of sample prediction for", toupper(ut), "models currently not available!"))
    rval <- fitted.bayesx(object, term = term)
  }

  rval
}


Predict.matrix.bayesx <- function(object, data) 
{
  UseMethod("Predict.matrix.bayesx")
}


Predict.matrix.bayesx.sm.bayesx <- function(object, data)
{
  specs <- attr(object, "specs")
  if(is.null(specs$type)) stop("cannot compute predict matrix, unknown term type!")
  if(is.null(specs$call)) stop("cannot compute predict matrix, unknown term call!")
  if(specs$dim < 2) {
    if(specs$type == "ps") {
      ## P-splines
      require("splines")
      term <- eval(parse(text = specs$call))
      p.order <- term$p.order[1L] + 1L
      nrknots <- term$bs.dim - p.order
      xr <- range(object[, term$term], na.rm = TRUE)
      dxr <- diff(xr)
      xr <- xr + dxr * c(-0.01, 0.01)
      step <- dxr / nrknots
      knots <- seq(xr[1] - p.order * step, xr[2] + p.order * step, by = step)
      X <- splineDesign(knots, data[[term$term]], ord = p.order + 1L, outer.ok = TRUE)
    } else stop("predict type not available!")
  } else {
    ## tensor splines
    require("splines")
    term <- eval(parse(text = specs$call))
    p.order <- term$p.order[1L] + 1L
    Xm <- list()
    for(j in term$term) {
      nrknots <- term$bs.dim - p.order
      xr <- range(object[, j], na.rm = TRUE)
      dxr <- diff(xr)
      xr <- xr + dxr * c(-0.01, 0.01)
      step <- dxr / nrknots
      knots <- seq(xr[1] - p.order * step, xr[2] + p.order * step, by = step)
      Xm[[j]] <- splineDesign(knots, data[[j]], ord = p.order + 1L, outer.ok = TRUE)
    }
    X <- mgcv:::tensor.prod.model.matrix(Xm)
  }

  return(X)
}


Predict.matrix.bayesx.mrf.bayesx <- Predict.matrix.bayesx.random.bayesx <- function(object, data)
{
  specs <- attr(object, "specs")
  if(is.null(specs$call)) stop("cannot compute predict matrix, unknown term call!")
  term <- eval(parse(text = specs$call))
  x <- unique(object[[term$term]])
  ndx <- data[[term$term]]
  if(is.factor(ndx))
    ndx <- f2int(ndx)
  if(!all(ndx %in% x)) stop(paste("newdata variable", term$term, "has unknown levels!"))
  nl <- length(x)
	X <- diag(nl)[factor(ndx, levels = x), ]

  return(X)
}


Predict.matrix.bayesx.geo.bayesx <- function(object, data)
{
  specs <- attr(object, "specs")
  if(specs$type == "gk")
    stop("predict matrix for geokriging term not supported yet!")
  if(is.null(specs$call)) stop("cannot compute predict matrix, unknown term call!")
  require("splines")
  term <- eval(parse(text = specs$call))
  p.order <- term$p.order[1L] + 1L
  t.data <- f2int(data[[term$term]])
  xy <- NULL
  if(!all(t.data %in% object[, term$term]))
    stop(paste("not all values in newdata variable", term$term, "in model term!"))
  for(i in unique(t.data)) {
    x <- object[object[, 1L] == i, "x"]
    y <- object[object[, 1L] == i, "x"]
    n <- sum(t.data == i)
    x <- rep(x, length.out = n)
    y <- rep(y, length.out = n)
    xy <- rbind(xy, cbind("x" = x, "y" = y))
  }
  Xm <- list()
  for(j in c("x", "y")) {
    xr <- range(xy[, j], na.rm = TRUE) + c(-0.001, 0.001)
    nrknots <- term$bs.dim - p.order
    step <- diff(xr) / nrknots
    knots <- seq(xr[1] - p.order * step, xr[2] + p.order * step, by = step)
    Xm[[j]] <- splineDesign(knots, xy[, j], ord = p.order + 1L, outer.ok = TRUE)
  }
  X <- matrix(0, nrow(data), 0)
  for(j in 1:ncol(Xm[[1L]]))
    X <- cbind(X, Xm[[1L]][, j] * Xm[[2L]])

  return(X)
}


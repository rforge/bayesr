GRstat <- function(object, term = NULL) {
  if((n <- length(object)) < 2L)
    stop("at least two models ar needed for calculation!")
  os <- samples(object, model = NULL, term)
  rval <- vector(mode = "list", length = length(os[[1L]]))
  for(i in 1:length(os[[1L]])) {
    chains <- list()
    for(j in 1:n) {
      chains[[j]] <- os[[j]][[i]]
    }
    if(is.list(chains[[1L]])) {
      for(j in 1:n)
        chains[[j]] <- as.data.frame(chains[[j]])
    }
    for(k in 1:ncol(chains[[1L]])) {
      grs <- NULL
      for(j in 1:n)
        grs <- cbind(grs, chains[[j]][, k])
      rval[[i]] <- c(rval[[i]], GR_calc(grs))
    }
    names(rval[[i]]) <- colnames(chains[[1L]])
  }
  names(rval) <- names(os[[1L]])

  rval
}


GR_calc <- function(theta) {
  if(!is.matrix(theta))
    stop("theta must be a matrix")

  nObs <- nrow(theta)
  nCols <- ncol(theta)
  n1 <- floor(nObs * 0.5)
  n2 <- nObs - n1

  if(nObs < 100)
    stop("There must be at least 100 observations from each chain")

  if(nCols < 2)
    stop("There must be at least two chains")

  theta <- theta[-(1:n1), ]

  vars <- apply(theta, 2, var)
  means  <- apply(theta, 2, mean)
  mBar <- mean(means)

  B <- n2*sum((means - mBar)^2) / (nCols - 1)
  W <- sum(vars) / nCols
  sigmaSq <- ((n2 - 1) * W + B) / (n2)
  vHat <- sigmaSq + B / (n2 * nCols)
  df <- n2
  R <- sqrt(vHat / W * (df / (df - 2)))

  return(R)
}


GRstats <- function(object, term = NULL, combine = TRUE, ...) {
  if(!inherits(object, "bayesx.hpc"))
    stop("cannot compute Gelman Rubin statistics of this object, object does not contain of parallel chains!")
  require("coda")
  if((n <- length(object)) < 2L)
    stop("at least two models ar needed for calculation!")
  os <- samples(object, model = NULL, term)
  nos <- names(os[[1L]])
  rval <- vector(mode = "list", length = length(os[[1L]]))
  if(combine)
    chains.combine <- vector(mode = "list", length = n)
  for(i in 1:length(os[[1L]])) {
    chains <- list()
    for(j in 1:n) {
      chains[[j]] <- os[[j]][[i]]
    }
    if(is.list(chains[[1L]])) {
      for(j in 1:n)
        chains[[j]] <- as.data.frame(chains[[j]])
    }
    if(combine) {
      for(j in 1:n) {
        if(is.null(chains.combine[[j]]))
          chains.combine[[j]] <- data.frame("id" = 1:nrow(chains[[j]]))
        if(is.null(dim(chains[[j]]))) {
          chains[[j]] <- data.frame(chains[[j]])
          names(chains[[j]]) <- NULL
        }
        colnames(chains[[j]]) <- paste(nos[i], colnames(chains[[j]]), sep = ".")
        chains.combine[[j]] <- cbind(chains.combine[[j]], chains[[j]])
      }
    } else {
      for(j in 1:n) {
        chains[[j]] <- mcmc(chains[[j]])
      }
      rval[[i]] <- gelman.diag(mcmc.list(chains), ...)
    }
  }
  if(combine) {
    for(j in 1:n) {
      chains.combine[[j]]$id <- NULL
      chains.combine[[j]] <- mcmc(chains.combine[[j]])
    }
    rval <- gelman.diag(mcmc.list(chains.combine), ...)
  } else {
    names(rval) <- names(os[[1L]])
  }

  rval
}


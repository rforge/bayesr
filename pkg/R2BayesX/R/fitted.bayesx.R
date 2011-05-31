fitted.bayesx <-
function(object, model = NULL, term = NULL, ...)
{    
  object <- get.model(object, model)
  k <- length(object)
  mn <- rep("model", length.out = k)
  if(is.null(term)) {
    if(length(object) > 1L) {
      rval <- vector("list", length = k)
      for(i in 1L:k) {
        rval[[i]] <- object[[i]]$fitted.values
        if(!is.null(object[[i]]$bayesx.setup$model.name))
          mn[i] <- object[[i]]$bayesx.setup$model.name
      }
      mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
      names(rval) <- mn
    } else rval <- object[[1L]]$fitted.values
  } else {
    if(length(object) > 1L) {
      rval <- list()
      for(i in 1L:k) {
        rval[[i]] <- object[[i]]$effects[term]
        if(!is.null(object[[i]]$bayesx.setup$model.name))
          mn[i] <- object[[i]]$bayesx.setup$model.name
      }
      mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
      names(rval) <- mn
      class(rval) <- "fit.bayesx"
    } else {
      if(length(term) > 1L)
        rval <- object[[1L]]$effects[term]
      else
        rval <- object[[1L]]$effects[[term]]
    }
  }
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}


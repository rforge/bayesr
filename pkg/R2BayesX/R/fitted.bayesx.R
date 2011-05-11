fitted.bayesx <-
function(object, model = NULL, term = NULL, ...)
{    
  object <- get.model(object, model)
  if(is.null(term)) {
    if(length(object) > 1L) {
      rval <- vector("list", length = length(object))
      for(i in 1L:length(object))
        rval[[i]] <- object[[i]]$fitted.values
      names(rval) <- names(object)
    } else rval <- object[[1L]]$fitted.values
  } else {
    if(length(object) > 1L) {
      rval <- list()
      for(i in 1L:length(object))
        rval[[i]] <- object[[i]]$effects[term]
      names(rval) <- names(object)
    } else rval <- object[[1L]]$effects[term]
    class(rval) <- "fit.bayesx"
  }
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}


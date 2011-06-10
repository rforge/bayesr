fitted.bayesx <-
function(object, model = NULL, term = NULL, fit.attributes = FALSE, ...)
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
        if(!fit.attributes) {
          attr(rval[[i]], "sample") <- NULL
          attr(rval[[i]], "variance.sample") <- NULL
          attr(rval[[i]], "specs") <- NULL
          attr(rval[[i]], "variance") <- NULL
          attr(rval[[i]], "partial.resids") <- NULL
          class(rval[[i]]) <- "matrix"
        }
        if(!is.null(object[[i]]$bayesx.setup$model.name))
          mn[i] <- object[[i]]$bayesx.setup$model.name
      }
      mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
      names(rval) <- mn
    } else {
      if(length(term) > 1L) {
        rval <- object[[1L]]$effects[term]
        if(!fit.attributes) {
          for(i in 1L:length(rval)) {
            attr(rval[[i]], "sample") <- NULL
            attr(rval[[i]], "variance.sample") <- NULL
            attr(rval[[i]], "specs") <- NULL
            attr(rval[[i]], "variance") <- NULL
            attr(rval[[i]], "partial.resids") <- NULL
            class(rval[[i]]) <- "matrix"
          }
        }
      } else {
        rval <- object[[1L]]$effects[[term]]
        if(!fit.attributes) {
          attr(rval, "sample") <- attr(rval, "variance.sample") <- attr(rval, "specs") <- NULL
          attr(rval, "variance") <- NULL
          attr(rval, "specs") <- NULL
          attr(rval, "variance") <- NULL
          attr(rval, "partial.resids") <- NULL
          class(rval) <- "matrix"
        }
      }
    }
  }
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA
  else {
    if(is.list(rval) && any(grepl("bayesx", class(rval[[1L]]))))
      class(rval) <- "fit.bayesx"
  }
  if(any(is.na(rval)))
    warning("fitted values are missing in object!")

  return(rval)
}


"[.fit.bayesx" <- function(x, term)
{
  if(is.list(x))
    return(x[[term]])
  else return(x)
}


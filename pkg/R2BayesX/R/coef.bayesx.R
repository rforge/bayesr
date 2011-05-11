coef.bayesx <- function(object, model = NULL, ...)
{
  object <- get.model(object, model)
  if(length(object) > 1L) {
    rval <- vector("list", length = length(object))
    for(i in 1L:length(object))
      rval[[i]] <- object[[i]]$fixed.effects
    names(rval) <- names(object)
  } else rval <- object[[1L]]$fixed.effects

  return(rval)
}

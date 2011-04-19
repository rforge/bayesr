fitted.bayesx <- function(object, model = NULL, term = NULL, ...)
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


residuals.bayesx <- function(object, model = NULL, ...)
{
  object <- get.model(object, model)
  if(length(object) > 1L) {
    rval <- vector("list", length = length(object))
    for(i in 1L:length(object))
      rval[[i]] <- object[[i]]$residuals
    names(rval) <- names(object)
  } else rval <- object[[1L]]$residuals
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}


plot.fit.bayesx <- function(x, ...)
{
  nx <- length(x)
  if(nx > 1L) {
    pval <- list()
    for(i in 1L:nx)
      pval[[i]] <- list(effects = x[[i]])
    class(pval) <- "bayesx"
  } else pval <- x[[1L]]
  plot(pval, ...)

  return(invisible(NULL))
}


terms.bayesx <- function(x, model = NULL, ...)
{
  x <- get.model(x, model)
  if(length(x) > 1L) {
    rval <- vector("list", length = length(x))
    for(i in 1L:length(x))
      rval[[i]] <- x[[i]]$bayesx.setup$terms
    names(rval) <- names(x)
  } else rval <- x[[1L]]$bayesx.setup$terms
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}

terms.bayesx <-
function(x, model = NULL, ...)
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


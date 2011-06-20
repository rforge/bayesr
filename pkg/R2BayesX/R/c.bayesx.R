c.bayesx <-
function(...)
{
  b <- list(...)
  x <- nx <- NULL
  for(i in 1L:length(b)) {
    x <- c(x, b[i])
    if(!is.null(b[[i]]$model.fit$model.name))
      nx <- c(nx, b[[i]]$model.fit$model.name)
    else
      nx <- c(nx, paste(i))
  }
  names(x) <- nx
  class(x) <- "bayesx"

  return(x)
}


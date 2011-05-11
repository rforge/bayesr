c.bayesx <- function(...)
{
  b <- list(...)
  x <- NULL
  for(i in 1L:length(b))
    x <- c(x, b[[i]])
  class(x) <- "bayesx"

  return(x)
}

c.bayesx <-
function(...)
{
  b <- list(...)
print(length(b))
  x <- NULL
  for(i in 1L:length(b))
    x <- c(x, b[i])
print(length(x))
  class(x) <- "bayesx"

  return(x)
}


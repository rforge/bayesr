f2int <- function(x, type = 1L) 
{
  if(!is.factor(x))
    x <- as.factor(x)
  levels(x) <- nl <- 1:nlevels(x)
  x <- factor(x, levels = nl, labels = nl)
  x <- as.integer(x)
  if(type != 2L)
    x <- x - 1L
  if(min(x, na.rm = TRUE) < 0)
    x <- x + abs(min(x, na.rm = TRUE))

  return(x)
}

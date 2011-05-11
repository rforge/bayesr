df2m <-
function(x)
{
  x$intnr <- NULL
  x <- as.matrix(x)
  rownames(x) <- 1L:nrow(x)

  return(x)
}


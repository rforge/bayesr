df2m <-
function(x)
{
  x$intnr <- NULL
  cn <- names(x)
  x <- as.matrix(x)
  rownames(x) <- 1L:nrow(x)
  colnames(x) <- cn

  return(x)
}


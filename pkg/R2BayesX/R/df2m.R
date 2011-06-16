df2m <-
function(x)
{
  xattr <- attributes(x)
  nxa <- names(xattr)
  x$intnr <- NULL
  cn <- names(x)
  x <- as.matrix(x)
  rownames(x) <- 1L:nrow(x)
  colnames(x) <- cn
  for(k in 1L:length(nxa)) 
    if(all(nxa[k] != c("dim", "dimnames", "class")))
      attr(x, nxa[k]) <- xattr[[k]]

  return(x)
}


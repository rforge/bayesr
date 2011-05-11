plot.fit.bayesx <-
function(x, ...)
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


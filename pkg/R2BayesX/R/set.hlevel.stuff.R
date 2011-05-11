set.hlevel.stuff <-
function(x, outfile)
{
  for(k in 1L:length(x)) {
    if(is.null(x[[k]]$hlevel))
      x[[k]]$hlevel <- 1L
    x[[k]]$hlevel <- x[[k]]$hlevel + 1L
    x[[k]]$first <- FALSE
    x[[k]]$outfile <- outfile
    if(!is.null(x[[k]]$h.random))
      x[[k]]$h.random <- set.hlevel.stuff(x[[k]]$h.random, outfile)
  }

  return(x)
}


is.rt <-
function(x)
{
  n <- length(x)
  isrt <- rep(FALSE, n)
  for(i in 1L:length(x)) {
    split <- splitme(x[i])
    if(length(split) > 2L) {
      if(resplit(split[1L:2L]) == c("r("))
        isrt[i] <- TRUE
    }
  }

  return(isrt)
}


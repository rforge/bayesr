is.f <-
function(x)
{
  if(length(x)) {
    n <- length(x)
    issm <- rep(FALSE, n)
    for(i in 1L:length(x)) {
      split <- splitme(x[i])
      if(length(split) > 2L) {
        if(resplit(split[1L:2L]) %in% "f(")
          issm[i] <- TRUE
      }
      if(length(split) > 4L) {
        if(resplit(split[1L:4L]) %in% "fbx(")
          issm[i] <- TRUE
      }
    }
  } else issm <- FALSE

  return(issm)
}


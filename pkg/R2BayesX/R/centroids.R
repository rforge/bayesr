centroids <-
function(map) 
{
  n <- length(map)
  cp <- matrix(0, n, 2L)
  for(i in 1L:n)
    cp[i,] <- centroidpos(map[[i]])

  return(cp)
}


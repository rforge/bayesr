f2int <-
function(x) 
{
  if(is.factor(x))
    levels(x) <- 1:nlevels(x)
  x <- as.integer(as.numeric(as.character(x)))

  return(x)
}


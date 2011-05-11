rrd <-
function(x) 
{
  x <- splitme(x)
  go <- TRUE
  rval <- NULL
  for(i in 1L:length(x)) {
    if(x[i] == ",")
      go <- FALSE
    if(go)
      rval <- c(rval, x[i])
  }
  rval <- resplit(c(rval, ")"))
  rval <- gsub("))", ")", rval)
  
  return(rval)
}


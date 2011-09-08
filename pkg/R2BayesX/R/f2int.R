f2int <-
function(x) 
{
  warn <- getOption("warn")
  options(warn = -1L)
  if(any(is.na(as.numeric(as.character(x)))))
    if(is.factor(x))
      levels(x) <- 1:nlevels(x)
  x <- as.integer(as.numeric(as.character(x))) - 1
  options(warn = warn)

  return(x)
}


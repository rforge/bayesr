print.bayesx <-
function(x, which = NULL, ...)
{
  n <- length(x)
  if(!is.null(which)) {
    start <- stop <- which
    n <- 1L
  } else {
    start <- 1L
    stop <- n
  }
  ncheck <- n > 1L
  for(i in start:stop) {
    if(ncheck)
      cat("###", i, "\n")
    .print_bayesx(x[[i]])
  }
  if(ncheck) {
    cat("###\n")
    cat("Object contains of", n, "models\n")
  }

  return(invisible(NULL))	
}


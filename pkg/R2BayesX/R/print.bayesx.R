print.bayesx <- function(x, which = NULL, ...)
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


.print_bayesx <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("---\n")
  if(!is.null(x$model.fit)) {
    mfn <- names(x$model.fit)
    step <- 5L
    for(i in 1L:length(mfn)) {
      if(i < step) {
        if(!is.null(x$model.fit[[i]]) && x$model.fit[[i]] != "") {
          if(mfn[i] != "step.final.model")
            cat(mfn[i], "=", x$model.fit[[i]]," ")
          else {
            cat("\n---\n")
            cat("Stepwise final model:\n")
            cat("-\n")
            cat(x$model.fit[[i]])
          }
        }
      }
      if(i == step) {
        if(i != length(mfn))
          cat("\n")
        step <- step + step
      }
    }
  cat("\n")
  }

  return(invisible(NULL))
}

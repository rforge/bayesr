getscript <-
function(object, file = NULL)
{
  if(!inherits(object, "bayesx"))
    stop("object must be of class \'bayesx\'")
  on <- deparse(substitute(object), backtick = TRUE, width.cutoff = 500L)
  mn <- names(object)
  object <- get.model(object, NULL)
  n <- length(object)
  script <- file("script", "w")
  cat("## The following R code may be applied on object ", on, ".\n", sep = "", file = script)
  cat("## Note that for map effect plots the corresponding\n", file = script)
  cat("## boundary object needs to be provided to function\n", file = script)
  cat("## plot.bayesx()!\n", file = script)
  cat("## Also see the help files for all available plots!\n\n", file = script)
  if(n > 1L)
    cat("## Note: object contains of", n, "models:\n\n", file = script)
  cat("## Summary statistics are provided with:\n", file = script)
  cat("summary(", on, ")", sep = "", file = script) 
  for(i in 1L:length(object)) {
    if(!is.null(object[[i]]$model.fit$method) && (object[[i]]$model.fit$method == "MCMC" ||
      object[[i]]$model.fit$method == "HMCMC")) {
      if(n > 1L)
        cat("\n\n## MCMC model diagnostic plots of model", mn[i], "\n", file = script)
      else
        cat("\n\n## MCMC model diagnostic plots:\n", file = script)
      cat("plot(", on, ", model = ", i, ", which = \"max-samples\")\n", sep = "", file = script)
      cat("plot(", on, ", model = ", i, ", which = \"max-samples\", acf = TRUE)", sep = "", file = script)
    }
    if(!is.null(object[[i]]$effects)) {
      tn <- names(object[[i]]$effects)
      if(n > 1L) {
        cat("\n\n## Plots of estimated effects of model: ", mn[i], sep = "", file = script) 
      } else {
        cat("\n\n## Plots of estimated effects", sep = "", file = script)
      } 
      for(k in 1L:length(tn)) {
        cat("\n## Plots for term:", tn[k], "\n", file = script)
        cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], "\")\n", sep = "", 
          file = script)
        if(!is.null(attr(object[[i]]$effects[[k]], "sample"))) {
          cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
            "\", which = \"coef-samples\")\n", sep = "", file = script)
          cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
            "\", which = \"coef-samples\", acf = TRUE)\n", sep = "", file = script)
        }
        if(!is.null(attr(object[[i]]$effects[[k]], "variance.sample"))) {
          cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
            "\", which = \"var-samples\")\n", sep = "", file = script)
          cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
            "\", which = \"var-samples\", acf = TRUE)\n", sep = "", file = script)
        }
      }
    }
  }
  close(script)
  script <- rval <- readLines("script")
  unlink("script")
  if(!is.null(file)) {
    file <- path.expand(file)
    cat("", file = file)    
    writeLines(script, con = file)
  }
  class(rval) <- "bayesx.script"

  return(invisible(rval))
}


print.bayesx.script <- function(x, ...)
{
  writeLines(x)
  return(invisible(NULL))
}


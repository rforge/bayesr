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
  cat("## Also see the help files of function plot.bayesx()\n## for all available plots!\n\n", 
    file = script)
  if(n > 1L)
    cat("## Note: object contains of", n, "models:\n\n", file = script)
  cat("## Summary statistics are provided with:\n", file = script)
  cat("summary(", on, ")\n", sep = "", file = script) 
  for(i in 1L:length(object)) {
    if(!is.null(object[[i]]$model.fit$method) && (object[[i]]$model.fit$method == "MCMC" ||
      object[[i]]$model.fit$method == "HMCMC")) {
      cat("\n", file = script)
      if(n > 1L)
        cat("## MCMC model diagnostic plots of model", mn[i], "\n", file = script)
      else
        cat("## MCMC model diagnostic plots:\n", file = script)
      cat("plot(", on, ", model = ", i, ", which = \"all-acf\")\n", sep = "", file = script)
    }
    cat("\n", file = script)
    if(!is.null(object[[i]]$effects)) {
      tn <- names(object[[i]]$effects)
      if(n > 1L) {
        cat("## Plots of estimated effects of model: ", mn[i], "\n", sep = "", file = script) 
      } else {
        cat("## Plots of estimated effects\n", sep = "", file = script)
      } 
      ccheck <- vcheck <- FALSE
      for(k in 1L:length(tn)) {
        cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], "\")\n", sep = "", 
          file = script)
        if(!is.null(attr(object[[i]]$effects[[k]], "sample")))
          ccheck <- TRUE
        if(!is.null(attr(object[[i]]$effects[[k]], "variance.sample")))
          vcheck <- TRUE
      }
      if(ccheck) {
        cat("\n## Coefficient sampling paths\n", file = script)
        for(k in 1L:length(tn)) {
          if(!is.null(attr(object[[i]]$effects[[k]], "sample"))) {
            cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
              "\", which = \"coef-samples\")\n", sep = "", file = script)
            cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
              "\", which = \"coef-samples\", acf = TRUE)\n", sep = "", file = script)
          }
        }
      }
      if(vcheck) {
        cat("\n## Variance sampling paths\n", file = script)
        for(k in 1L:length(tn)) {
          if(!is.null(attr(object[[i]]$effects[[k]], "variance.sample"))) {
            cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
              "\", which = \"var-samples\")\n", sep = "", file = script)
            cat("plot(", on, ", model = ", i, ", term = ", "\"", tn[k], 
              "\", which = \"var-samples\", acf = TRUE)\n", sep = "", file = script)
          }
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


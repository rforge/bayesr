coeffun <- function(x, args, diagnostics)
{
  if(is.null(x))
    return(invisible(NULL))			
  args$var <- FALSE
  if(diagnostics == 2) {
    args$var <- TRUE
    args$x <- attr(x, "variance.sample")
  } else args$x <- attr(x, "sample")
  if(!is.null(args$x)) {
    args$selected <- attr(x, "specs")$label
    do.call("plotsamples", args)
  }

  return(invisible(NULL))
}

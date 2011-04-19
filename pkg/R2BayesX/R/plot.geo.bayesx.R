plot.geo.bayesx <- function(x, diagnostics = FALSE, ...)
{
  if(!is.null(x)) {
    args <- list(...)
    if(diagnostics == FALSE) {
      p3d <- TRUE
      if(is.null(args$map)) {
        plottype <- "plot3d"
        args$x <- x[, 2L:ncol(x)]
      } else {
        p3d <- FALSE
        plottype <- "plotmap"
        args$x <- x[,c(1L, 4L:ncol(x))]
        class(args$x) <- "mrf.bayesx"
      }
      an <- names(attributes(x))
      for(att in an)
        if(!att %in% c("dim","dimnames")) {
          if(p3d)
            attr(args$x, att) <- attr(x, att)
          else
            attr(args$x, att) <- attr(x, att)
        }
      do.call(plottype, args)
    } else coeffun(x, args, diagnostics)
  }

  return(invisible(NULL))
}

plot.bayesx <-
function(x, model = NULL, term = NULL, which = 1L, ask = FALSE, ...)
{
  op <- par(no.readonly = TRUE)
  which.match <- c("effect", "coef-samples", "var-samples", "intcpt-samples", 
    "hist-resid", "qq-resid", "scatter-resid", "scale-resid", "scale-samples")
  if(!is.character(which)) {
    if(any(which > 9))
      which <- which[which <= 9]
    which <- which.match[which]
  } else which <- which.match[pmatch(which, which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")
  x <- get.model(x, model)
  nx <- length(x)
  if((!"effect" %in% which) && (!"coef-samples" %in% which) 
    && (!"var-samples" %in% which) && (!"intcpt-samples" %in% which)) {
    model.names <- names(x)
    for(i in 1L:nx)
      which.plots(x[[i]], which, ask, model.names[i], nx, ...)
    if(nx > 1L || length(which) > 1L)		
      par(op)
  } else {
    if(is.null(term) && !ask) {
      nt <- 0L
      for(i in 1L:nx)
        nt <- nt + length(x[[i]]$effects)
      if(nt > 1L)
        setmfrow(nt)
    } else {
      nt <- neffects(x, term)
      if(!ask)
        if(nt > 1L) 
          setmfrow(nt)
    }
    args <- list(...)
    for(i in 1L:nx) {
      if("intcpt-samples" %in% which) {
        if(!is.null(attr(x[[i]]$fixed.effects, "sample"))) {
          par(oma = c(1, 1, 2, 1))
          par(mfrow = c(1, 1))
          args$x <- attr(x[[i]]$fixed.effects, "sample")[,1L]
          args$selected <- "(Intercept)"
          args$var <- FALSE
          do.call("plotsamples", args)	
        }
      } else {
        if(is.null(term))
          ts <- 1:length(x[[i]]$effects)
        else
          ts <- term
        ne <- names(x[[i]]$effects)
        if(is.null(ne) || !is.character(ts))
          ne <- 1L:length(x[[i]]$effects)
        for(j in 1L:length(ts)) {
          if(is.character(ts[j])) {
            tmp <- splitme(ts[j])
            tmp <- resplit(tmp[tmp != " "])
            take <- pmatch(tmp, ne)
          } else take <- match(ts[j], ne)
          if(length(take) > 0L && length(x[[i]]$effects) > 0L && !is.na(take)) {
            args$x <- x[[i]]$effects[[take]]
            args$diagnostics <- FALSE
            if("coef-samples" %in% which) {
              par(oma = c(1, 1, 2, 1))
              args$diagnostics <- 1L
            }
            if("var-samples" %in% which) {
              par(oma = c(1, 1, 2, 1))
              args$diagnostics <- 2L
            }
            if(!is.null(args$x)) {
              args$ask <- ask
              do.call("plot", args)	
              if(j == 1L)
                if(ask)
                  par(ask = TRUE)
            }
          }
        }
      }	
    }		
    if(nt > 1L)		
      par(op)
  }

  return(invisible(NULL))
}


plot.bayesx <- function(x, model = NULL, term = NULL, which = 1L, ask = FALSE, ...)
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


which.plots <- function(x, which, ask, model.name, nx, ...)
{
  args <- list(...)
  nw <- length(which)
  if((!ask || nx > 1L) && nw > 1L)
    par(mfrow = n2mfrow(nw))
  else
    ask <- TRUE  
  if(nx > 1L)
    ask <- TRUE
  if(nw > 1L && ask && nx > 1L) {
    par(oma = c(1, 1, 5, 1))
    if(is.null(args$main))
      main <- paste("Diagnostic plots for model", model.name)
    else
      main <- args$main
  }
  residuals <- x$residuals
  if(is.matrix(residuals))
    residuals <- residuals[,1L]
  residuals <- as.numeric(residuals)
  fitvalues <- x$fitted.values
  if(is.matrix(fitvalues))
    fitvalues <- fitvalues[,1L]
  fitvalues <- as.numeric(fitvalues)
  k <- 0L
  ok <- FALSE
  for(ww in which) {
    k <- k + 1L
    if(ww == "hist-resid" && length(residuals) > 0L) {
      args2 <- args
      dens <- density(residuals)
      hst <- hist(residuals, plot = FALSE)
      args2$ylim <- c(0, max(c(hst$density, dens$y)))
      args2$xlab <- "Residuals"
      args2$ylab <- "Density"
      args2$main <- "Histogramm and density"
      args2$freq <- FALSE
      args2$x <- residuals
      args2 <- delete.args(graphics::hist.default, args2)
      args2$cex <- args$cex
      do.call(graphics::hist.default, args2)
      graphics::lines(dens)
      box()
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "qq-resid" && length(residuals) > 0L) {
      args2 <- args
      r2 <- (residuals - mean(residuals))/sd(residuals)
      args2$y <- r2
      args2$ylab <- "Standardized residuals"
      args2$xlab <- "Theoretical quantiles"
      args2$main <- "Normal Q-Q Plot"
      args2$ylim <- NULL
      args2$xlim <- NULL
      args2 <- delete.args(stats::qqnorm.default, args2)
      args2$cex <- args$cex
      do.call(stats::qqnorm, args2)
      stats::qqline(r2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "scatter-resid" && length(residuals) > 0L  && length(fitvalues) > 0L) {
      args2 <- args
      args2$y <- residuals
      args2$x <- fitvalues
      args2 <- delete.args(stats::scatter.smooth, args2)
      args2$xlab <- "Fitted values"
      args2$ylab <- "Residuals"
      args2$main <- "Fitted values vs. residuals"
      args2$cex <- args$cex
      do.call(stats::scatter.smooth, args2)
      abline(h = 0, lty = 2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
    if(ww == "scale-resid" && length(residuals) > 0L  && length(fitvalues) > 0L) {
      args2 <- args
      args2$y <- sqrt(abs((residuals - mean(residuals))/sd(residuals)  ))
      args2$x <- fitvalues
      args2 <- delete.args(stats::scatter.smooth, args2)
      args2$xlab <- "Fitted values"
      args2$ylab <- expression(sqrt(abs("Standardized residuals")))
      args2$main <- "Scale-location" 
      args2$cex <- args$cex
      do.call(stats::scatter.smooth, args2)
      if(k == 1L && nx < 2L)
        par(ask = ask) 
      ok <- TRUE
    }
##!!! FIXME    if(ww == "scale-samples") {

  }
  if(nw > 1L && ask && nx > 1L && ok)
    mtext(main, side = 3L, line = 2L, outer = TRUE, font = 2, cex = 1)
  if(k == 1L && nx > 1L)
    par(ask = ask) 

  return(invisible(NULL))
}


neffects <- function(x, term)
{
  n <- 0L
  for(i in 1L:length(x)) {
    ne <- names(x[[i]]$effects)
    if(is.null(ne) || !is.character(term))
      ne <- 1L:length(x[[i]]$effects)
    for(j in 1L:length(term)) {
      if(is.character(term[j])) {
        tmp <- splitme(term[j])
        tmp <- resplit(tmp[tmp != " "])
        take <- pmatch(tmp, ne)
      } else take <- match(term[j], ne)
      if(length(take) > 0L && length(x[[i]]$effects) > 0L && !is.na(take))
        n <- n + 1L
    }
  }

  return(n)
}

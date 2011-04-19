smoothone <- function(x, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, main = NULL, 
                      colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, cex = 1, ...)
	{
	if(col[1] != FALSE)
		color <- col
	else
		color <- rgb(0.01,0.01,0.01,0.2)
	args <- as.list(substitute(list(...)))[-1]
	by <- attr(x,"specs")$by
	if(is.null(by))
		by <- "NA"
	name <- attr(x,"specs")$term
	if(is.null(name))
		name <- colnames(x)[1]
	f <- x[,2]
	f01 <- x[,3]
	f011 <- x[,4]
	f021 <- x[,6]
	f02 <- x[,7]
	e <- x[,10]
	ceffect <- attr(x,"ceffect")
        x <- x[,1]
	icheck <- FALSE
	if(by[1]!="NA")
		{
		if(any(by == 0))
			{
			icheck <- TRUE
			ind <- 1:length(f)
			ind <- ind[by!=0]
			be <- e[ind]
			bx <- x[ind]
			f <- f[ind]
			f01 <- f01[ind]
			f02 <- f02[ind]
			f011 <- f011[ind]
			f021 <- f021[ind]
			e <- e[ind]
			x <- x[ind]
			}
		#else
		#	{
			#f <- f/by
			#f01 <- f01/by
			#f02 <- f02/by
			#f011 <- f011/by
			#f021 <- f021/by
		#	}	
		if(length(name)>1)	
			byname <- name[length(name)]
		else
			byname <- by
		name <- name[1]
		}
	if(const)
		{
		if(!is.null(ceffect))
			{
			f <- f + ceffect
			f01 <- f01 + ceffect
			f02 <- f02 + ceffect
			f011 <- f011 + ceffect
			f021 <- f021 + ceffect
			e <- e + ceffect
			}
		}
				
	fmax <- max(f)
	fmin <- min(f)
	ymax01 <- max(f01)
	ymax02 <- max(f02)
	ymin01 <- min(f01)
	ymin02 <- min(f02)
	emax <- max(e)
	emin <- min(e)
	
	if(length(name) > 1)
		name <- name[1]
	if(is.null(xlab))
		xlab <- name
	if(is.null(ylab))
		{
		ylab <- paste("Effect of", xlab)
		if(by[1]!="NA")
			ylab <- paste(ylab,":",byname,sep="")
		}

	if(resid != FALSE || resid != F)
		{
		maxy <- max(c(ymax01,ymax02,emax,fmax))
		miny <- min(c(ymin01,ymin02,emin,fmin))

		if(is.null(args$ylim))
			args$ylim <- substitute(c(miny,maxy))
	
		args$y <- substitute(e)
		args$x <- substitute(x)
		args$xlab <- xlab
		args$ylab <- ylab
		args$main <- main
		args$cex <- cex
		args$type <- "n"

		do.call("plot",args)
		if(colored != FALSE || colored != F)
			{
			if(length(color) == 1)
				color <- c(color,color)
			xo <- x[order(x)]
			xc <- c(xo,xo[length(xo):1])
			f02o <- f02[order(x)]
			yc <- c(f01[order(x)],f02o[length(f02o):1])
			f02o2 <- f021[order(x)]
			yc2 <- c(f011[order(x)],f02o2[length(f02o2):1])
			polygon(xc, yc, col=color[1], lty=0)
			polygon(xc, yc2, col=color[2], lty=0)
			}
		lines(f[order(x)]~x[order(x)], lwd=lwdc)
		if(lwdconf!=0)
			{
			lines(f01[order(x)]~x[order(x)],lty=3,lwd=lwdconf)
			lines(f02[order(x)]~x[order(x)],lty=3,lwd=lwdconf)
			lines(f011[order(x)]~x[order(x)],lty=2,lwd=lwdconf)
			lines(f021[order(x)]~x[order(x)],lty=2,lwd=lwdconf)
			}
		# if(icheck)
		#	points(x=bx,y=be,cex=cex)
		# else
			points(x=x,y=e,cex=cex)
		}
	else
		{
		maxy <- max(c(ymax01,ymax02,fmax))
		miny <- min(c(ymin01,ymin02,fmin))
		#maxy <- maxy + 0.15*abs(maxy)
		#miny <- miny - 0.15*abs(miny)

		if(is.null(args$ylim))
			args$ylim <- substitute(c(miny,maxy))

		args$y <- substitute(f)
		args$x <- substitute(x)
		args$xlab <- xlab
		args$ylab <- ylab
		args$main <- main
		args$type <- "n"

		do.call("plot",args)
		if(colored != FALSE || colored != F)
			{
			if(length(color) == 1)
				color <- c(color,color)
			xo <- x[order(x)]
			xc <- c(xo,xo[length(xo):1])
			f02o <- f02[order(x)]
			yc <- c(f01[order(x)],f02o[length(f02o):1])
			f02o2 <- f021[order(x)]
			yc2 <- c(f011[order(x)],f02o2[length(f02o2):1])
			polygon(xc, yc, col=color[1], lty=0)
			polygon(xc, yc2, col=color[2], lty=0)
			}
		lines(f[order(x)]~x[order(x)], lwd=lwdc)
		if(lwdconf!=0)
			{
			lines(f01[order(x)]~x[order(x)],lty=3,lwd=lwdconf)
			lines(f02[order(x)]~x[order(x)],lty=3,lwd=lwdconf)
			lines(f011[order(x)]~x[order(x)],lty=2,lwd=lwdconf)
			lines(f021[order(x)]~x[order(x)],lty=2,lwd=lwdconf)
			}
		if(jit)
			rug(jitter(x))
		}
	}

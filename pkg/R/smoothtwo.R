smoothtwo <- function(X,grid,xlab,ylab,zlab,main,theta,phi,image,const,col,colored,resid,cex,border,nrc,pal,...)
	{
	x <- X[,1]
	z <- X[,2]
	e <- X[,11]
	specs <- attr(X,"smooth.specs")
	names <- specs$term
	by <- specs$by
	if(is.null(by))
		by <- "NA"
	else
		{
		if(length(specs$term)>2)
			by <- specs$term[length(specs$term)]
		}
	call <- specs$call
	draws <- attr(X,"coef.draws.utr")
	secheck <- FALSE
	args <- as.list(substitute(list(...)))[-1]
	ipcheck <- TRUE
	if(!is.null(args$ip) || is.null(draws))
		{
		args$ip <- NULL
		ipcheck <- TRUE
		}
	if(!is.null(args$se))
		{
		secheck <- TRUE
		args$se <- NULL
		}
	ceffect <- attr(X,"ceffect")
	xn <- seq(min(x),max(x),length=grid)
	zn <- seq(min(z),max(z),length=grid)
	if(ipcheck)
		{
		fitted <- list()
		fitted$mean <- akima::interp(x,z,X[,3],xo=xn,yo=zn,duplicat="strip")$z
		if(secheck)
			{
			fitted$q1 <- akima::interp(x,z,X[,4],xo=xn,yo=zn)$z
			fitted$q2 <- akima::interp(x,z,X[,8],xo=xn,yo=zn)$z
			}
		names <- colnames(X)[1:2]
		}
	else
		{
		XN <- rep(xn,grid)
		ZN <- rep(zn,rep(grid,grid))
		dat <- eval(parse(text=paste("data.frame(",specs$term[1],"=XN,",specs$term[2],"=ZN)")))
		if(!is.null(specs$xt$geo))
			specs$xt$geo <- NULL
		X <- Predict.matrix(attr(specs,"smooth.construct"),dat)
		if(is.null(const))
			const <- FALSE
		fitted <- qhelp(draws,X)
		}
	fit01 <- fit02 <- NA
	if(const)
		{
		fit <- fitted$mean + ceffect
		if(secheck)
			{
			fit01 <- fitted$q1 + ceffect
			fit02 <- fitted$q2 + ceffect
			}
		e <- e + ceffect
		}
	else
		{
		fit <- fitted$mean 
		if(secheck)
			{
			fit01 <- fitted$q1
			fit02 <- fitted$q2 
			}
		}

	if(resid != FALSE || resid != F)
		zlimit <- range(c(fit,fit01,fit02,e),na.rm=TRUE)
	else
		zlimit <- range(c(fit,fit01,fit02),na.rm=TRUE)
	if(is.null(xlab))
		xlab <- names[1]
	if(is.null(ylab))
		ylab <- names[2]
	if(is.null(zlab))
		{
		zlab <- "fitted"
		if(by!="NA")
			zlab <- paste(zlab,":",by,sep="")
		}
	if(secheck)	
		myfit <- matrix(fit02, grid, grid)
	else
		myfit <- matrix(fit, grid, grid)
	if(!image)
		{
		if(colored)
			{
			if(col == FALSE || col == F)
				{
				if(is.null(pal) || is.na(pal))
					pal <- "heat"
				col <- ret.pal.colors(nrc,pal,myfit)
				}
			}
		else 
			col <- "white"
		args$y <- substitute(zn)
		args$x <- substitute(xn)
		args$z <- substitute(myfit)
		args$xlab <- xlab
		args$ylab <- ylab
		args$zlab <- zlab
		args$main <- main
		args$col <- col
		if(is.null(args$zlim))
			args$zlim <- zlimit
		args$theta <- theta
		args$phi <- phi
		args$border <- border
		if(secheck)
			{
			args$border <- "green"
			args$col <- NA
			pmat <- do.call("persp",args)

			par(new = TRUE)
			args$border <- "black"
			myfit <- matrix(fit, grid, grid)
			args$z <- substitute(myfit)
			pmat <- do.call("persp",args)

			myfit <- matrix(fit01, grid, grid)
			args$z <- substitute(myfit)			
			args$border <- "red"
			par(new = TRUE)
			pmat <- do.call("persp",args)
			}
		else
			pmat <- do.call("persp",args)
		if(resid != FALSE || resid != F) 
			{
			t3d <- trans3d(x,z,e,pmat)
			points(x=t3d$x,y=t3d$y,cex=cex)
			}
		}
	if(image)
		{
		if(colored)
			{
			if(col != F || col != FALSE)
				colfill <- col
			else
				{
				colfill <- heat.colors(100)
				colfill <- colfill[100:1]
				}
			}
		if(is.null(args$xlim))
			xlim <- range(xn)
		else
			xlim <- args$xlim
		if(is.null(args$ylim))
			ylim <- range(zn)
		else
			ylim <- args$ylim
		image.plot(xn,zn,myfit,xlab=xlab,ylab=ylab,main=main,col=colfill,xlim=xlim,ylim=ylim)
		if(!is.null(list(...)$contour))
			if(list(...)$contour)
				contour(xn,zn,myfit, add = TRUE, lty = "solid")
		}
	}

plot.mrf.gibbs <- function(x, resid = FALSE, map = NULL, names = FALSE, values = NULL, colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0,
		           range = NULL, pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, dgts = 4, 
		           cex = 1, lpos = c(0.2,0), const = FALSE, diagnostics = FALSE, acf = FALSE, p3d = FALSE, theta = 0, phi = 25, pcat = NULL, border = NULL, ...)
	{
	if(diagnostics == FALSE)
		{
		if(inherits(x,"random.gibbs"))
			{
			values <- x[,2:8]
			region <- x[,1]
			}
		else
			{
			values <- attr(x,"smooth.coef")
			region <- as.integer(levels(as.factor(x[,1])))
			}
		if(const)
			{
			if(inherits(x,"random.gibbs"))
				{
				values <- values + attr(x,"random.ceffect")
				x[,2:7] <- x[,2:7] + attr(x,"random.ceffect")
				}
			else
				{
				values <- values + attr(x,"smooth.ceffect")
				x[,2:7] <- x[,2:7] + attr(x,"smooth.ceffect")
				}
			}
		if(is.null(map))
			{
			if(is.null(range))
				{
				up <- 0.3
				dow <- 0.3	
				}
			else
				{
				if(length(range) < 2)
					range <- c(range,range)
				dow <- range[1]
				up <- range[2]	
				}
			nl <- nrow(values)
			coef1 <- matrix(rep(values[,2],3),nl,3)
			coef2 <- matrix(rep(values[,6],3),nl,3)
			coef11 <- matrix(rep(values[,3],3),nl,3)
			coef21 <- matrix(rep(values[,5],3),nl,3)
			coef <- matrix(rep(values[,1],3),nl,3)
			if(inherits(x,"random.gibbs"))
				e <- attr(x,"random.partial.resid")
			else
				{
				res <- x[,8]
				flev <- as.integer(levels(as.factor(x[,1])))
				e <- vector("list",length=nl)
				for(j in 1:nl)
					{
					e[[j]] <- res[x[,1]==flev[j]]
					if(const)
						e[[j]] <- e[[j]] + attr(x,"smooth.ceffect")
					}
				}

			fnames <- levels(as.factor(x[,1]))
			xtmp <- xtmp2 <- matrix(0,nl,3)
			limchecke <- NULL
			limcheck <- c(coef1[,1],coef2[,1],coef11[,1],coef21[,1])
			for(j in 1:nl)
				{
				limchecke <- c(limchecke,e[[j]])
				e[[j]] <- cbind(rep(j,length(e[[j]])),e[[j]])
                 	 	xtmp[j,] <- c(j-dow,j,j+up)
                  	xtmp2[j,] <- c(j+up,j,j-dow)
                  	}
            	if(resid != F || resid != FALSE)
				y <- c(min(limchecke),max(limchecke))
			else
				y <- c(min(limcheck),max(limcheck))
			xp <- c(0.5,nl+0.5)	

			if(colored != FALSE)
				{
				if(col != FALSE)
					color <- col
				else
					color <- rgb(0.01,0.01,0.01,0.2)
				}
			else
				color <- "white"
			args <- as.list(substitute(list(...)))[-1]
			if(length(color) == 1)
				color <- c(color,color)
			
			if(is.null(args$ylim))
				args$ylim <- substitute(y)
			xnam <- colnames(x)[1]
			if(is.null(xlab))
				xlab <- xnam
			if(is.null(ylab))
				{
				if(inherits(x,"random.gibbs"))
					ylab <- paste("Estimated coefficients of",attr(x,"random.term"))
				else
					ylab <- paste("Estimated coefficients of",xnam)
				}
			args$y <- substitute(y)
			args$x <- substitute(xp)
			args$xlab <- xlab
			args$ylab <- ylab
			args$main <- main
			args$cex <- cex
			args$type <- "n"
			args$axes <- FALSE
			do.call("plot",args)
			box()
			axis(2)
			axis(1,at=1:nl,labels=fnames)
			
			for(k in 1:nl)
				{
				f1 <- cbind(c(xtmp[k,],xtmp2[k,]),c(coef2[k,],coef1[k,]))
				f2 <- cbind(c(xtmp[k,],xtmp2[k,]),c(coef21[k,],coef11[k,]))
				polygon(f1,col=color[1],lty=0)
				polygon(f2,col=color[2],lty=0)
				if(resid!=F || resid != FALSE)
					points(e[[k]],cex=cex)
				lines(coef[k,]~xtmp[k,],lwd=lwdc)
				if(lwdconf!=0)
					{
					lines(coef1[k,]~xtmp[k,],lwd=lwdconf,lty=3)
					lines(coef2[k,]~xtmp[k,],lwd=lwdconf,lty=3)
					lines(coef11[k,]~xtmp[k,],lwd=lwdconf,lty=2)
					lines(coef21[k,]~xtmp[k,],lwd=lwdconf,lty=2)
					}
				}
			}
		else
			{
			drawmap(map,names,x,colored,range,pal,legend,scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos,p3d,theta,phi,pcat,border,resid=resid,...)		
			if(inherits(x,"random.gibbs"))
				cat("Estimated random effects:\n")
			else
				cat("Estimated regional effects:\n")
			rn <- rep("",nrow(values))
			rownames(values) <- region
			printCoefmat(values)		
			}
		}
	else
		{
		if(!inherits(x,"random.gibbs"))
			{
			if(diagnostics == 2)
				coefplot(attr(x,"smooth.variance.draws"),xlab,ylab,main,attr(x,"smooth.specs")$label,acf,TRUE,...)
			else
				coefplot(attr(x,"smooth.coef.draws"),xlab,ylab,main,attr(x,"smooth.specs")$label,acf,FALSE,...)
			}
		}
	return(invisible(NULL))
	}

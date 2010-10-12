plot.linear.gibbs <- function(x, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, 
                              main = NULL, colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, range = NULL,
                              cex = 1, ask = FALSE, ...) 
	{
	if(is.null(attr(x,"factorcheck")))
		attr(x,"factorcheck") <- "nonfactor"
	if(attr(x,"factorcheck") == "nonfactor")
		{
		for(j in 2:ncol(x))
			x[,j] <- x[,j] - mean(x[,j])
		nams <- colnames(x)
		attr(x,"specs") <- list(term=nams[1],by="NA")
		if(diagnostics == FALSE)
			smoothone(x,resid,jit,const,diagnostics,acf,xlab,ylab,main, 
		    	    	    colored,col,lwdc,lwdconf,cex,...)
		else
			coefplot(attr(x,"coef.draws.utr"),xlab,ylab,main,colnames(x)[1],acf,FALSE,...)
		}
	else
		{
		if(diagnostics == FALSE)
			{
			colfill <- col
			xnam <- attr(x,"factorname")
			x <- x[[1]]
			if(is.null(range))
				{
				dow <- 0.3
				up <- 0.3
				}
			else
				{
				dow <- range[1]
				up <- range[2]	
				}
			lx <- length(x)
			xt <- ft <- et <- f01t <- f02t <- f011t <- f021t <- vector("list",lx)
			fnames <- limcheck <- limchecke <- NULL
			for(j in 1:lx)
				{
				X <- x[[j]]
				many <- sum(X[,1])
				xtmp <- c(j-dow,j,j+up)
				xtmp2 <- c(j+up,j,j-dow)
				check <- X[,1]!=0
				f <- X[,2][check]
				f <- rep(f[1],3)
				f01 <- X[,2][check]
				f01 <- rep(f01[1],3)
				f02 <- X[,7][check]
				f02 <- rep(f02[1],3)
				f011 <- X[,4][check]
				f011 <- rep(f011[1],3)
				f021 <- X[,6][check]
				f021 <- rep(f021[1],3)
				e <- X[,10][check]
				limcheck <- c(limcheck,f02[1],f021[1],f01[1],f011[1])
				limchecke <- c(limchecke,e)
				xt[[j]] <- xtmp	
				ft[[j]] <- cbind(xtmp,f)
				et[[j]] <- cbind(e,rep(j,many))
				f01t[[j]] <- cbind(xtmp2,f01)
				f02t[[j]] <- cbind(xtmp,f02)
				f011t[[j]] <- cbind(xtmp2,f011)
				f021t[[j]] <- cbind(xtmp,f021)
				tmpn <- colnames(X)[1]
				fnames <- c(fnames,tmpn)
				}
			if(resid!=F || resid != FALSE)
				y <- c(min(limchecke),max(limchecke))
			else
				y <- c(min(limcheck),max(limcheck))
			xt <- c(0.5,lx+0.5)
			
			if(colored != FALSE)
				{
				if(colfill != FALSE)
					color <- colfill
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
			if(is.null(xlab))
				xlab <- xnam
			if(is.null(ylab))
				ylab <- paste("Effects of factor",xnam)
			args$y <- substitute(y)
			args$x <- substitute(xt)
			args$xlab <- xlab
			args$ylab <- ylab
			args$main <- main
			args$cex <- cex
			args$type <- "n"
			args$axes <- FALSE

			do.call("plot",args)
			box()
			axis(2)
			axis(1,at=1:lx,labels=fnames)
			
			for(k in 1:lx)
				{
				f1 <- rbind(f02t[[k]],f01t[[k]])
				f2 <- rbind(f021t[[k]],f011t[[k]])
				polygon(f1,col=color[1],lty=0)
				polygon(f2,col=color[2],lty=0)
				if(resid!=F || resid != FALSE)
					points(et[[k]][,1]~et[[k]][,2],cex=cex)
				lines(ft[[k]],lwd=lwdc)
				if(lwdconf!=0)
					{
					lines(f01t[[k]][,2]~ft[[k]][,1],lwd=lwdconf,lty=3)
					lines(f02t[[k]][,2]~ft[[k]][,1],lwd=lwdconf,lty=3)
					lines(f011t[[k]][,2]~ft[[k]][,1],lwd=lwdconf,lty=2)
					lines(f021t[[k]][,2]~ft[[k]][,1],lwd=lwdconf,lty=2)
					}
				}
			}
		else
			{
			x <- x[[1]]
			lx <- length(x)
			nx  <- length(attr(x[[1]],"coef.draws"))
			dmat <- matrix(0,lx,nx)
			for(j in 1:lx)
				dmat[j,] <- attr(x[[1]],"coef.draws.utr")
			coefplot(dmat,xlab,ylab,main,strsplit(colnames(X)[1],":")[[1]][1],acf,FALSE,...)
			}
		}
	}

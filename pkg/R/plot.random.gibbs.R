plot.random.gibbs <- function(x, resid = FALSE, const = FALSE, diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, zlab = NULL,
                              main = NULL, colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, map = NULL, names = FALSE, 
                              range = NULL, pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, dgts = 2, 
                              lpos = c(0.2,0), cex = 1, ask = FALSE, p3d = FALSE, theta = 0, phi = 25, pcat = NULL, border = NULL, ...)
	{	
	if(diagnostics == FALSE)
		{
		if(!is.null(map))
			{
			values <- x
			plot.mrf.gibbs(x,resid,map,names,values,colored,col,lwdc,lwdconf,
					range,pal,legend,scale,nrc,xlab,ylab,zlab,main,dgts,
					cex,lpos,const,diagnostics,acf,p3d,theta,phi,pcat,border,...)					
			}
		else
			{
			cefo <- 0
			if(const)
				cefo <- attr(x,"random.ceffect")
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
			nl <- nrow(x)
			coef1 <- matrix(rep(x[,4]+cefo,3),nl,3)
			coef2 <- matrix(rep(x[,8]+cefo,3),nl,3)
			coef11 <- matrix(rep(x[,5]+cefo,3),nl,3)
			coef21 <- matrix(rep(x[,7]+cefo,3),nl,3)
			coef <- matrix(rep(x[,2]+cefo,3),nl,3)
			e <- attr(x,"random.partial.resid")
			fnames <- levels(as.factor(x[,1]))
			yescheck <- FALSE
			if(!is.null(attr(x,"byplots")))
				{
				fnames <- rownames(x)
				yescheck <- TRUE
				}
			xtmp <- xtmp2 <- matrix(0,nl,3)
			limchecke <- NULL
			limcheck <- c(coef1[,1],coef2[,1],coef11[,1],coef21[,1])
			for(j in 1:nl)
				{
				limchecke <- c(limchecke,e[[j]]+cefo,coef1+cefo,coef2+cefo)
				e[[j]] <- cbind(rep(j,length(e[[j]])),e[[j]]+cefo)
                 	 	xtmp[j,] <- c(j-dow,j,j+up)
                  	xtmp2[j,] <- c(j+up,j,j-dow)
                  	}
            	if(resid!=F || resid != FALSE)
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
				if(yescheck)
					ylab <- paste("Estimated coefficients of",attr(x,"random.term"))
				else
					ylab <- paste("Random effects of",xnam)
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
		}
	else
		{
		if(diagnostics == 2)
			coefplot(attr(x,"random.variance.draws"),xlab,ylab,main,attr(x,"random.term"),acf,TRUE,...)
		else
			coefplot(attr(x,"random.coefs.draws"),xlab,ylab,main,attr(x,"random.term"),acf,FALSE,...)
		}
	}

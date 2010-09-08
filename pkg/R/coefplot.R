coefplot <- function(draws, xlab, ylab, main, selected, acf, var, ...)
	{
	yes <- FALSE
	if(!is.matrix(draws))
		{
		yes <- TRUE
		draws <- matrix(draws,1,length(draws))
		}
	nr <- nrow(draws)
	args <- as.list(substitute(list(...)))[-1]
	setmfrow(nr)

	for(k in 1:nr)
		{
		par(mar=c(0,0,0,0))
		if(!acf)
			{
			plot(draws[k,], type="l", main="",xlab="", lwd=1,...)
			box()
			}
		else
			{
			tims <- ts(data=draws[k,])
			acf(tims,lag.max=args$lag.max,main="",xlab="",verbose=FALSE,...)
			}
		if(!is.null(xlab))
			mtext(xlab, side=1, line=3, outer=TRUE, font=1)
		if(is.null(ylab))
			mtext(ylab, side=2, line=3, outer=TRUE, font=1)
		if(!is.null(main))
			mtext(main, side=3, line=2, outer=TRUE, font=2, cex=1)
		if(is.null(main))
			{
			if(var)
				mtext(paste("Sampling path of the variance of term",selected), side=3, line=2, outer=TRUE, font=2, cex=1)
			else
				{
				if(yes)
					{
					if(!acf)
						mtext(paste("Sampling path of coefficient of term",selected), side=3, line=2, outer=TRUE, font=2, cex=1)
					else
						mtext(paste("Autocorrelation function of sampled coefficient of term",selected), side=3, line=2, outer=TRUE, font=2, cex=1)
					}
				else
					{
					if(!acf)
						mtext(paste("Sampling paths of coefficients of term",selected), side=3, line=2, outer=TRUE, font=2, cex=1)
					else
						mtext(paste("Autocorrelation functions of sampled coefficients of term",selected), side=3, line=2, outer=TRUE, font=2, cex=1)
					}
				}
			}
		}
	}

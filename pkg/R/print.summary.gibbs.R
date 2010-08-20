print.summary.gibbs <- function(x,...)
	{
	cat("Call:\n")
	print(x$call)

	tcheck <- rownames(x$coeflin)
	tcheck2 <- as.character(x$coeflin)
	if(x$K > 0)
		{
		tcheck <- c(tcheck,rownames(x$smoothhyp))
		tcheck2 <- c(tcheck2,as.character(x$smoothhyp))
		}
	if(x$R > 0)
		{
		if(!is.null(x$ra$coeflin))
			{
			tcheck <- c(tcheck,rownames(x$ra$coeflin))
			tcheck2 <- c(tcheck2,as.character(x$ra$coeflin))
			}
		if(!is.null(x$ra$smoothhyp))
			{
			tcheck <- c(tcheck,rownames(x$ra$smoothhyp))
			tcheck2 <- c(tcheck2,as.character(x$ra$smoothhyp))
			}
		}
	if(x$M > 0)
		{
		if(!is.null(x$coeflinmu))
			{
			tcheck <- c(tcheck,rownames(x$coeflinmu))
			tcheck2 <- c(tcheck2,as.character(x$coeflinmu))
			}
		if(!is.null(x$smoothmuhyp))
			{
			tcheck <- c(tcheck,rownames(x$smoothmuhyp))
			tcheck2 <- c(tcheck2,as.character(x$smoothmuhyp))
			}
		}
	ltcheck <- length(tcheck)
	maxsize <- rep(0,ltcheck)
	for(j in 1:ltcheck)
		maxsize[j] <- length(strsplit(tcheck[j],"","")[[1]])
	maxsize <- max(maxsize)

	ltcheck2 <- length(tcheck2)
	maxsize2 <- rep(0,ltcheck2)
	for(j in 1:ltcheck2)
		maxsize2[j] <- length(strsplit(tcheck2[j],"","")[[1]])
	maxsize2 <- max(maxsize2)

	if(maxsize2 > 7)
		maxsize <- maxsize + maxsize2*7 + 7	
	else
		maxsize <- maxsize + maxsize2*7 + 7
	liner <- paste(rep("-",maxsize),collapse="")

	fc <- TRUE
	if(nrow(x$coeflin) < 2)
		{
		if(all(x$coeflin[1,]==0))
			fc <- FALSE
		}
	else
		{
		if(all(x$coeflin[1,]==0))
			{
			m <- ncol(x$coeflin)
			nc <- colnames(x$coeflin)
			nr <- rownames(x$coeflin)[2:nrow(x$coeflin)]
			x$coeflin <- matrix(x$coeflin[2:nrow(x$coeflin),],ncol=m)
			colnames(x$coeflin) <- nc
			rownames(x$coeflin) <- nr
			}
		}
	if(fc || (x$K > 0))
		{
		cat(liner,"\n")
		cat("Fixed effects estimation results:\n")
		cat("---\n")
		}
	if(fc)
		{
		cat("Parametric Coefficients:\n")
		printCoefmat(x$coeflin)
		}
	
	if(x$K > 0)
		{
		if(fc)
			cat("-\n")
		cat("Smooth terms variances:\n")
		ls <- ncol(x$smoothhyp)
		terms <- colnames(x$smoothhyp)	
		rn <- rownames(x$smoothhyp)
		#for(j in 1:ls)
			#{
			#cat(terms[j],"\n")
			#newmat <- matrix(x$smoothhyp[j,],nrow=1,ncol=ncol(x$smoothhyp))
			#rownames(newmat) <- ""
			#colnames(newmat) <- rn
			#printCoefmat(newmat)
			#}
		printCoefmat(x$smoothhyp)
		}
	cat(liner,"\n")

	if(x$R > 0)
		{
		cat("Random effects estimation results:\n")
		cat("---\n")

		if(!is.null(x$ra$coeflin))
			{
			cat("Parametric Coefficients:\n")
			printCoefmat(x$ra$coeflin)
			cat("---\n")
			}
		if(!is.null(x$ra$smoothhyp))
			{
			cat("Smooth terms variances:\n")
			printCoefmat(x$ra$smoothhyp)
			cat("---\n")
			}
		cat("Random variance estimation results:\n")
		cat("---\n")
		printCoefmat(x$ra$s2r)		
		cat(liner,"\n")
		}
	if(x$M > 0)
		{
		cat("Multilevel effects estimation results:\n")
		if(!is.null(x$coeflinmu))
			{
			cat("---\n")
			cat("Parametric Coefficients:\n")
			printCoefmat(x$coeflinmu)
			}
		if(!is.null(x$smoothmuhyp))
			{
			cat("---\n")
			cat("Smooth terms variances:\n")
			printCoefmat(x$smoothmuhyp)
			}
		cat(liner,"\n")
		}
	cat("Global variance estimation results:\n")
	printCoefmat(x$S2)
	cat(liner,"\n")
	cat("DIC = ",x$DIC$DIC,"  pd = ",x$DIC$pd,"  n = ",x$N,"  samptime = ",x$samptime,"sec\n", sep="")
	cat("Iterations =",x$iter," burnin =",x$burn," thinning =",x$thin,"\n")
	}

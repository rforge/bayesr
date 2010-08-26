extract <- function(cdraws,lambda,response,eta,Z,C,K,check,new,mcdraws,S2,sm.specs,resp.ind,data)
	{
	smooth.out <- vector("list", K)
	smooth.mat <- matrix(0,K,7)
	matrownames <- rep("",K)
	for(k in 1:K)
		{
		if(check)
			{
			fout <- qhelp2(cdraws[[k]],Z[[k]]$tildeZ,response,eta,resp.ind[[k]])
			redraws <- redraw(Z[[k]]$RQ,cdraws[[k]])
			}
		X <- NULL
		if(class(sm.specs[[k]])!="mrf.smooth.spec")
			{
			if(!is.null(sm.specs[[k]]$xt$geo))
				{
				tmp <- get("geospline.centroid.data",envir=.GlobalEnv)
				X <- cbind(tmp[[1]],tmp[[2]])
				tmp <- cbind(X,tmp[[3]])
				tmp <- unique(tmp)
				colnames(tmp) <- names(tmp)
				assign("geospline.centroid.data",tmp,envir=.GlobalEnv)
				}
			else
				{
				for(j in 1:sm.specs[[k]]$dim)
					{
					if(is.null(sm.specs[[k]]$data))
						X <- cbind(X,eval(parse(text=sm.specs[[k]]$term[j]),envir=data))
					else
						X <- cbind(X,eval(parse(text=sm.specs[[k]]$term[j]),envir=sm.specs[[k]]$data))
					}
				}
			}
		else
			{
			if(is.null(sm.specs[[k]]$specs$data))
				X <- eval(parse(text=sm.specs[[k]]$term[1]),envir=data)
			else
				X <- eval(parse(text=sm.specs[[k]]$term[1]),envir=sm.specs[[k]]$specs$data)
			}
		fout <- cbind(X,fout)
		if(sm.specs[[k]]$by!="NA" || (!is.null(sm.specs[[k]]$loc) && (length(sm.specs[[k]]$terms) == 3)))
			{
			tmp <- sm.specs[[k]]$term
			sm.specs[[k]]$term <- c(sm.specs[[k]]$term,sm.specs[[k]]$by)
			sm.specs[[k]]$by <- eval(parse(text=sm.specs[[k]]$by),envir=data)
			if(is.factor(sm.specs[[k]]$by))
				sm.specs[[k]]$by <- 1*(sm.specs[[k]]$by==levels(sm.specs[[k]]$by)[1])
			sm.specs[[k]]$by <- sm.specs[[k]]$by[order(fout[,1])]
			}
		else
			{
			if(class(sm.specs[[k]])!="mrf.smooth.spec")
				tmp <- sm.specs[[k]]$term
			else
				tmp <- sm.specs[[k]]$term[1]
			}
		if(!is.null(sm.specs[[k]]$data$byokcheck))
			fout <- fout[sm.specs[[k]]$data$byokcheck,]
		fout <- fout[order(fout[,1]),]
		colnames(fout) <- c(tmp,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")

		gamma <- sfout(redraws,type=2)
		gcn <- paste(sm.specs[[k]]$label,":",1:nrow(gamma),sep="")
		rownames(gamma) <- gcn
		matrownames[k] <- sm.specs[[k]]$label
		attr(fout,"smooth.coef") <- gamma
		attr(fout,"s") <- Z[[k]]$s
		attr(fout,"rq") <- Z[[k]]$RQ
		attr(fout,"smooth.hyper") <- sfout(lambda[k,])
		attr(fout,"smooth.hyper.draws") <- lambda[k,]
		attr(fout,"smooth.variance") <- sfout(S2/lambda[k,])
		attr(fout,"smooth.variance.draws") <- S2/lambda[k,]
		attr(fout,"smooth.ceffect") <- mean(C[k,])
		attr(fout,"smooth.ceffect.draws") <- C[k,]
		attr(fout,"smooth.coef.draws") <- cdraws[[k]]
		attr(fout,"smooth.coef.mean") <- mcdraws[[k]]
		attr(fout,"smooth.coef.draws.utr") <- redraws
		attr(fout,"smooth.specs") <- sm.specs[[k]]
		attr(fout,"term.type") <- "smooth"
		
		if(class(sm.specs[[k]])!="mrf.smooth.spec")
			class(fout) <- "sm.gibbs"
		else
			class(fout) <- "mrf.gibbs"
		smooth.out[[k]] <- fout
		smooth.mat[k,] <- attr(fout,"smooth.hyper")
		}
	rownames(smooth.mat) <- matrownames
	colnames(smooth.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")

	return(list(smooth.out,smooth.mat))
	}

extract <- function(cdraws,lambda,response,eta,Z,C,K,check,new,mcdraws,S2,sm.specs,resp.ind,data)
	{
	smooth.out <- vector("list", K)
	smooth.mat <- matrix(0,K,7)
	matrownames <- rep("",K)
	for(k in 1:K)
		{
		if(check)
			{
			fitted <- qhelp(cdraws[[k]],Z[[k]]$tildeZ)
			redraws <- redraw(Z[[k]]$RQ,cdraws[[k]])
			pqu2p5 <- fitted$q2
			pqu97p5 <- fitted$q1
			pqu10 <- fitted$q21
			pqu90 <- fitted$q11
			pmean <- fitted$mean
			pmed <- fitted$median
			pmean <- pmean[new[[k]]]
			pmed <- pmed[new[[k]]]
			pqu2p5 <- pqu2p5[new[[k]]]
			pqu97p5 <- pqu97p5[new[[k]]]
			pqu10 <- pqu10[new[[k]]]
			pqu90 <- pqu90[new[[k]]]
			partial.resid <- as.vector(response[resp.ind[[k]]] - (eta[resp.ind[[k]]] - pmean))
			}
		pcat80 <- (pqu10 < 0 & pqu90 < 0)*(-1) + (pqu10 <= 0 & pqu90 >= 0)*0 + (pqu10 > 0 & pqu90 > 0)*1
		pcat95 <- (pqu2p5 < 0 & pqu97p5 < 0)*(-1) + (pqu2p5 <= 0 & pqu97p5 >= 0)*0 + (pqu2p5 > 0 & pqu97p5 > 0)*1
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
		fout <- cbind(X,pmean,pqu2p5,pqu10,pmed,pqu90,pqu97p5,partial.resid,pcat95,pcat80)
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

		gmean <- apply(redraws, 1, mean)
		gpqu2p5 <- apply(redraws, 1, quantile, probs=0.025)
		gpqu10 <- apply(redraws, 1, quantile, probs=0.1)
		gpmed <- apply(redraws, 1, quantile, probs=0.5)
		gpqu90 <- apply(redraws, 1, quantile, probs=0.9)
		gpqu97p5 <- apply(redraws, 1, quantile, probs=0.975)
		gamma <- cbind(gmean,gpqu2p5,gpqu10,gpmed,gpqu90,gpqu97p5) 
		colnames(gamma) <- c("pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
		gcn <- paste(sm.specs[[k]]$label,":",1:nrow(gamma),sep="")
		rownames(gamma) <- gcn
		attr(fout,"smooth.coef") <- gamma
		attr(fout,"s") <- Z[[k]]$s
		attr(fout,"rq") <- Z[[k]]$RQ

		lamed <- median(lambda[k,])
		lamean <- mean(lambda[k,])
		lasd <- sd(lambda[k,])
		lapqu2p5 <- quantile(lambda[k,], probs = 0.025)
		lapqu10 <- quantile(lambda[k,], probs = 0.1)
		lapqu90 <- quantile(lambda[k,], probs = 0.9)
		lapqu97p5 <- quantile(lambda[k,], probs = 0.975)
		la <- cbind(lamean,lasd,lapqu2p5,lapqu10,lamed,lapqu90,lapqu97p5)
		colnames(la) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
		rownames(la) <- NULL
		attr(fout,"smooth.hyper") <- la
		attr(fout,"smooth.hyper.draws") <- lambda[k,]

		taumed <- median(S2/lambda[k,])
		taumean <- mean(S2/lambda[k,])
		tausd <- sd(S2/lambda[k,])
		taupqu2p5 <- c(quantile(S2/lambda[k,], probs = 0.025))
		taupqu10 <- c(quantile(S2/lambda[k,], probs = 0.1))
		taupqu90 <- c(quantile(S2/lambda[k,], probs = 0.9))
		taupqu97p5 <- c(quantile(S2/lambda[k,], probs = 0.975))
		tau2 <- cbind(taumean,tausd,taupqu2p5,taupqu10,taumed,taupqu90,taupqu97p5)
		colnames(tau2) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
		matrownames[k] <- sm.specs[[k]]$label
		attr(fout,"smooth.variance") <- tau2
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
		smooth.mat[k,] <- la
		}
	rownames(smooth.mat) <- matrownames
	colnames(smooth.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")

	return(list(smooth.out,smooth.mat))
	}

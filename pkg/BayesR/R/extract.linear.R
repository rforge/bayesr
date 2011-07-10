extract.linear <- function(L,K,draws,mdraws,X,LIN,response,eta,cL,tnames,plus,intcheck,itcpt)
	{
	pmed <- psd <- pmean <- pqu2p5 <- pqu10 <- pqu90 <- pqu97p5 <- rep(0,(L+plus))
	fits <- NULL
	redraws <- redraw(X$RQ,draws)
	if(itcpt > 0 || L > 0)
		{
		for(k in 1:K)
			{
 			pmed[k] <- median(redraws[k,])
			psd[k] <- sd(redraws[k,])
			pmean[k] <- mean(redraws[k,])
 			pqu2p5[k] <- quantile(redraws[k,], probs = (0.05/2))
 			pqu97p5[k] <- quantile(redraws[k,], probs = (1-0.05/2))
 			pqu10[k] <- quantile(redraws[k,], probs = 0.1)
 			pqu90[k] <- quantile(redraws[k,], probs = 0.9)
			}
		}	
	beta <- cbind(pmean,psd,pqu2p5,pqu10,pmed,pqu90,pqu97p5)
	bcolnam <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	if(itcpt < 1)
		beta[1,] <- rep(0,7)
	if(L > 0)
		{
		fits <- vector("list", L)
		tmpceffect <- apply(cL, 1, mean)
		for(k in 1:L)
			{
			fpqu2p5 <- LIN[[k]]*pqu2p5[k+plus]
			fpqu10 <- LIN[[k]]*pqu10[k+plus]
			fpqu90 <- LIN[[k]]*pqu90[k+plus]
			fpqu97p5 <- LIN[[k]]*pqu97p5[k+plus]
			fpmean <- LIN[[k]]*pmean[k+plus]
			fpmed <- LIN[[k]]*pmed[k+plus]
			e <- response - (eta - fpmed)
			pcat80 <- (fpqu10 < 0 & fpqu90 < 0)*(-1) + (fpqu10 <= 0 & fpqu90 >= 0)*0 + (fpqu10 > 0 & fpqu90 > 0)*1
			pcat95 <- (fpqu2p5 < 0 & fpqu97p5 < 0)*(-1) + (fpqu2p5 <= 0 & fpqu97p5 >= 0)*0 + (fpqu2p5 > 0 & fpqu97p5 > 0)*1
			fout <- cbind(LIN[[k]],fpmean,fpqu2p5,fpqu10,fpmed,fpqu90,fpqu97p5,e,pcat95,pcat80)
			fout <- fout[order(fout[,1]),]
			colnames(fout) <- c(tnames[k],"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
			attr(fout,"linear.coef") <- beta[k+plus,]
			attr(fout,"linear.ceffect") <- tmpceffect[k]
			attr(fout,"linear.coef.draws") <- draws[k+plus,]
			attr(fout,"linear.coef.draws.utr") <- redraws[k+plus,]
			attr(fout,"linear.coef.mean") <- mdraws[k+plus,]
			class(fout) <- "linear.gibbs"
			fits[[k]] <- fout
			}
		if(intcheck)
			brownam <- c("(Intercept)",tnames)
		else
			brownam <- tnames
		}
	else
		{
		beta <- as.matrix(beta,1,7)
		brownam <- c("(Intercept)",tnames)
		}
	attr(beta,"linear.coef.draws") <- draws
	attr(beta,"linear.coef.draws.utr") <- redraws
	attr(beta,"linear.coef.mean") <- mdraws
	attr(beta,"dimnames") <- list(brownam,bcolnam)
	class(beta) <- "gibbsfit"

	if(L > 0)
		return(list(beta=beta,fits=fits))
	else
		return(list(beta=beta))
	}

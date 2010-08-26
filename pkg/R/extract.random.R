extract.random <- function(R,RR,RRl,RRs,RRi,fit,RA,lambda,Z,drawsc,mdrawsc,dPR,response,eta,drawb,drawbm,S2R,ins,
                           drawTZRl,drawTZRlM,desRl,design,metaR,newR,icheck,iZ,cRMS,cRSS,cRLL,S2o,ZRA,which,data,stetaR)
	{
	refout <- iNS <- vector("list", R)
	which <- which + 1	

	for(k in 1:R)
		{
		varmat <- NULL
		# pmean <- apply(drawb[[k]], 1, mean)
		# psd <- apply(drawb[[k]], 1, sd)
		# pqu97p5 <- apply(drawb[[k]], 1, quantile, probs = 0.975)
		# pqu2p5 <- apply(drawb[[k]], 1, quantile, probs = 0.025)
		# pqu90 <- apply(drawb[[k]], 1, quantile, probs = 0.9)
		# pqu10 <- apply(drawb[[k]], 1, quantile, probs = 0.1)
		# pmed <- apply(drawb[[k]], 1, median)

		fitted <- qhelp(drawb[[k]],ZRA[[k]]$RQ)
		redraws <- redraw(ZRA[[k]]$RQ,drawb[[k]])
		# fitted <- apply(redraws,1,quantile,probs=c(0.025,0.975,0.1,0.9,0.5))
		pqu2p5 <- fitted$q2
		pqu97p5 <- fitted$q1
		pqu10 <- fitted$q21
		pqu90 <- fitted$q11
		pmean <- fitted$mean
		pmed <- fitted$median
		psd <- apply(drawb[[k]],1,sd)

		partial.resid <- vector("list",ncol(RA[[k]]$Z))
		for(j in 1:ncol(RA[[k]]$Z))
			{
			tmp <- pmean[j]*RA[[k]]$Z[,j]
			if(icheck)
				{
				tmp <- tmp[iZ]
				ztmp <- RA[[k]]$Z[iZ,]
				}
			tmp <- response - eta + tmp
			if(icheck)
				partial.resid[[j]] <- tmp[ztmp[,j]!=0]
			else
				partial.resid[[j]] <- tmp[RA[[k]]$Z[,j]!=0]
			}

		var <- sfout(S2R[k,])
		rownames(var) <- RA[[k]]$facname
		if(is.null(RA[[k]]$byfacid))
			faclev <- as.integer(levels(as.factor(RA[[k]]$facorig)))
		else
			faclev <- RA[[k]]$byfacid

		pcat80 <- (pqu10 < 0 & pqu90 < 0)*(-1) + (pqu10 <= 0 & pqu90 >= 0)*0 + (pqu10 > 0 & pqu90 > 0)*1
		pcat95 <- (pqu2p5 < 0 & pqu97p5 < 0)*(-1) + (pqu2p5 <= 0 & pqu97p5 >= 0)*0 + (pqu2p5 > 0 & pqu97p5 > 0)*1

		effects <- cbind(faclev,pmean,psd,pqu2p5,pqu10,pmed,pqu90,pqu97p5,pcat95,pcat80)
		colnames(effects) <- c(RA[[k]]$facname,"pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5","pcat95","pcat80")

		if(is.null(RA[[k]]$byfac))
			rownames(effects) <- paste(RA[[k]]$facname,":",faclev,sep="")
		else
			{
			rownames(effects) <- RA[[k]]$byfac
			nlby <- length(RA[[k]]$bywhat)
			byplots <- vector("list",length=nlby)
			specs <- attr(RA[[k]]$bywhat,"specs")
			lxby <- length(specs)
			names <- vector("list",length=lxby)

			for(j in 1:nlby)
				{
				check <- attr(RA[[k]]$bywhat[[j]],"used")
				byplots[[j]] <- vector("list",length=lxby)
				for(jj in 1:lxby)
					{
					if(RA[[k]]$rintc)
						ok <- RA[[k]]$bywhat[[j]][jj+1,1]:RA[[k]]$bywhat[[j]][jj+1,2]
					else
						ok <- RA[[k]]$bywhat[[j]][jj,1]:RA[[k]]$bywhat[[j]][jj,2]
					Ztmp <- RA[[k]]$Z[,ok]
					# Ztmp <- ZRA[[k]]$tildeZ[,ok]
					dtmp <- drawb[[k]][ok,]	
					redtmp <- redraws[ok,]
					if(is.vector(Ztmp))
						{
						xtmp <- specs[[jj]]$x[check]
						tmp <- effects[ok,]
						pqu2p5 <- tmp[4]*xtmp
						pqu97p5 <- tmp[8]*xtmp
						pqu10 <- tmp[5]*xtmp
						pqu90 <- tmp[7]*xtmp
						pmean <- tmp[2]*xtmp
						pmed <- tmp[6]*xtmp
						by.partial <- response[check] - eta[check] + pmean
						pcat80 <- (pqu10 < 0 & pqu90 < 0)*(-1) + (pqu10 <= 0 & pqu90 >= 0)*0 + (pqu10 > 0 & pqu90 > 0)*1
						pcat95 <- (pqu2p5 < 0 & pqu97p5 < 0)*(-1) + (pqu2p5 <= 0 & pqu97p5 >= 0)*0 + (pqu2p5 > 0 & pqu97p5 > 0)*1
						fitted <- cbind(pmean,pqu2p5,pqu10,pmed,pqu90,pqu97p5,by.partial,pcat95,pcat80)
						byplots[[j]][[jj]] <- cbind(xtmp,fitted)[order(xtmp),]
						}
					else
						{
						fitted <- qhelp2(drawb[[k]],ZRA[[k]]$tildeZ,response,eta)
						if(class(specs[[jj]]) == "mrf.smooth")
							xtmp <- specs[[jj]]$x[check]
						else
							xtmp <- specs[[jj]]$x[check,]
						byplots[[j]][[jj]] <- cbind(xtmp,fitted[check,])[order(xtmp),]
						}
					if(nlby > 1)
						tmpnam <- paste(specs[[jj]]$term,".",j,sep="")
					else
						tmpnam <- specs[[jj]]$term
					colnames(byplots[[j]][[jj]]) <- c(tmpnam,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
					gamma <- matrix(effects[ok,c(2,4,5,6,7,8)],ncol=6)
					colnames(gamma) <- c("pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
					rownames(gamma) <- RA[[k]]$byfac[ok]
					if(specs[[jj]]$bytype == "linear")
						{
						attr(byplots[[j]][[jj]],"term.type") <- "linear"
						attr(byplots[[j]][[jj]],"linear.ceffect") <- mean(cRMS[k,])
						attr(byplots[[j]][[jj]],"linear.coef.draws") <- dtmp
						attr(byplots[[j]][[jj]],"linear.coef.draws.utr") <- redraws[ok,]
						attr(byplots[[j]][[jj]],"linear.coef.mean") <- drawbm[[k]][ok,]
						}
					else
						{
						attr(byplots[[j]][[jj]],"smooth.specs") <- specs[[jj]]$specs
						attr(byplots[[j]][[jj]],"smooth.specs")$term <- tmpnam
						attr(byplots[[j]][[jj]],"smooth.specs")$label <- rownames(RA[[k]]$bywhat[[j]])[jj]
						attr(byplots[[j]][[jj]],"term.type") <- "smooth"
						attr(byplots[[j]][[jj]],"smooth.coef") <- gamma
						attr(byplots[[j]][[jj]],"smooth.coef.draws") <- dtmp
						attr(byplots[[j]][[jj]],"smooth.variance.draws") <- S2R[k,]
						attr(byplots[[j]][[jj]],"smooth.ceffect") <- mean(cRMS[k,])
						attr(byplots[[j]][[jj]],"smooth.coef.mean") <- drawbm[[k]][ok,]
						attr(byplots[[j]][[jj]],"smooth.coef.draws.utr") <- redraws[ok,]
						}
					if(class(specs[[jj]])!="mrf.smooth.spec")
						class(byplots[[j]][[jj]]) <- "sm.gibbs"
					else
						class(byplots[[j]][[jj]]) <- "mrf.gibbs"
					}
				}
			attr(effects,"byplots") <- byplots
			}
		attr(effects,"random.variance") <- var
		attr(effects,"random.variance.draws") <- S2R[k,]
		attr(effects,"random.coefs.draws") <- drawb[[k]]
		attr(effects,"random.coefs.mean") <- drawbm[[k]]
		attr(effects,"random.partial.resid") <- partial.resid
		attr(effects,"random.ceffect") <- mean(cRMS[k,])
		attr(effects,"random.term") <- RA[[k]]$call
		attr(effects,"s") <- 1/diag(crossprod(RA[[k]]$Z))
		varmat <- rbind(varmat,var)
		class(effects) <- "random.gibbs"

		if(RA[[k]]$hspec)
			{
			fitted <- qhelp(stetaR[[k]],ZRA[[k]]$RQ)
			pqu2p5 <- fitted$q2
			pqu97p5 <- fitted$q1
			pqu10 <- fitted$q21
			pqu90 <- fitted$q11
			pmean <- fitted$mean
			pmed <- fitted$median
			psd <- apply(stetaR[[k]],1,sd)
			
			partial.resid <- vector("list",ncol(RA[[k]]$Z))
			for(j in 1:ncol(RA[[k]]$Z))
				{
				tmp <- pmean[j]*RA[[k]]$Z[,j]
				if(icheck)
					{
					tmp <- tmp[iZ]
					ztmp <- RA[[k]]$Z[iZ,]
					}
				tmp <- response - eta + tmp
				if(icheck)
					partial.resid[[j]] <- tmp[ztmp[,j]!=0]
				else
					partial.resid[[j]] <- tmp[RA[[k]]$Z[,j]!=0]
				}
			faclev <- as.integer(levels(as.factor(RA[[k]]$facorig)))

			pcat80 <- (pqu10 < 0 & pqu90 < 0)*(-1) + (pqu10 <= 0 & pqu90 >= 0)*0 + (pqu10 > 0 & pqu90 > 0)*1
			pcat95 <- (pqu2p5 < 0 & pqu97p5 < 0)*(-1) + (pqu2p5 <= 0 & pqu97p5 >= 0)*0 + (pqu2p5 > 0 & pqu97p5 > 0)*1

			effects.i <- cbind(faclev,pmean,psd,pqu2p5,pqu10,pmed,pqu90,pqu97p5,pcat95,pcat80)
			colnames(effects.i) <- c(RA[[k]]$facname,"pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5","pcat95","pcat80")

			attr(effects.i,"random.variance") <- var
			attr(effects.i,"random.variance.draws") <- S2R[k,]
			attr(effects.i,"random.coefs.draws") <- stetaR[[k]]
			attr(effects.i,"random.coefs.mean") <- drawbm[[k]]
			attr(effects.i,"random.partial.resid") <- partial.resid
			attr(effects.i,"random.ceffect") <- mean(cRMS[k,])
			attr(effects.i,"random.term") <- RA[[k]]$call
			attr(effects.i,"s") <- 1/diag(crossprod(RA[[k]]$Z))
			class(effects.i) <- "random.gibbs"
			attr(effects,"effects.nototal") <- effects.i

			betas <- betastmp <- smooth.mat.tmp <- NULL
			ra.fout <- list()
			mapind <- 1:RA[[k]]$length
			if(RRs[k] > 0)
				{
				smooth.mat <- matrix(0,RRs[k],7)
				matrownames <- rep("",RRs[k])
				tmpcheck <- 1:ncol(RA[[k]]$indspec)

				for(j in 1:RRs[k])
					{
					jj <- tmpcheck[RA[[k]]$indspec[2,]==j]
					if(is.null(RA[[k]][[jj]]$rwcheck))
						{
						fout <- qhelp2(drawsc[[k]][[j]],Z[[k]][[j]]$tildeZ,response,eta,newR[[k]][[j]])
						fout <- fout[RA[[k]]$fac,]
						xtmp <- NULL
						for(i in 1:RA[[k]][[jj]]$specs$dim)
							xtmp <- cbind(xtmp,eval(parse(text=RA[[k]][[jj]]$specs$term[i]),envir=data))
						if(nrow(xtmp)!=length(RA[[k]]$fac))
							xtmp <- xtmp[RA[[k]]$fac,]
						}
					else
						{
						fout <- qhelp2(drawsc[[k]][[j]],ZRA[[k]]$tildeZ,response,eta)
						xtmp <- RA[[k]][[jj]]$rwx
						RA[[k]][[jj]]$names <- paste(RA[[k]][[jj]]$rwnames,":rw",sep="")
						}

					if(icheck)
						{
						fout <- fout[iZ,]
						if(is.matrix(xtmp))
							xtmp <- xtmp[iZ,]
						else
							xtmp <- xtmp[iZ]
						}

					fout <- cbind(xtmp,fout)
					fout <- fout[order(fout[,1]),]

					if(class(RA[[k]][[jj]]) == "mrf.smooth")
						tmpnam <- RA[[k]][[jj]]$names[1]
					else
						tmpnam <- RA[[k]][[jj]]$names
					colnames(fout) <- c(tmpnam,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
					matrownames[j] <- RA[[k]][[jj]]$term
					redrawsc <- redraw(Z[[k]][[j]]$RQ,drawsc[[k]][[j]])
					gamma <- sfout(redrawsc,type=2)
					gcn <- paste(RA[[k]][[jj]]$label,":",1:nrow(gamma),sep="")
					rownames(gamma) <- gcn
					attr(fout,"smooth.coef") <- gamma
					attr(fout,"s") <- Z[[k]][[j]]$s
					attr(fout,"rq") <- Z[[k]][[j]]$RQ
					attr(fout,"smooth.hyper") <- sfout(lambda[[k]][[j]])
					attr(fout,"smooth.hyper.draws") <- lambda[[k]][[j]]
					attr(fout,"smooth.variance") <- sfout(S2o/lambda[[k]][[j]])
					attr(fout,"smooth.variance.draws") <- S2o/lambda[[k]][[j]]
					attr(fout,"smooth.ceffect") <- mean(cRSS[[k]][[j]])
					attr(fout,"smooth.ceffect.draws") <- cRSS[[k]][[j]]
					attr(fout,"smooth.coef.draws") <- drawsc[[k]][[j]]
					attr(fout,"smooth.coef.mean") <- mdrawsc[[k]][[j]]
					attr(fout,"smooth.coef.draws.utr") <- redrawsc

					if(!is.null(RA[[k]][[jj]]$rwcheck))
						sm.specs <- RA[[k]][[jj]]$rwspecs
					else
						sm.specs <- RA[[k]][[jj]]$specs
					attr(fout,"smooth.specs") <- sm.specs
					attr(fout,"term.type") <- "smooth"
					if(!is.null(RA[[k]][[jj]]$rwcheck))
						class(fout) <- paste(RA[[k]][[jj]]$rwtype,".gibbs",sep="")
					else
						{
						if(class(sm.specs)!="mrf.smooth.spec")
							class(fout) <- "sm.gibbs"
						else
							class(fout) <- "mrf.gibbs"
						}
					ra.fout[[mapind[RA[[k]]$indspec[2,] == j]]] <- fout
					smooth.mat[k,] <- attr(fout,"smooth.hyper")
					}
				rownames(smooth.mat) <- matrownames
				colnames(smooth.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
				attr(ra.fout,"smooth.mat") <- smooth.mat
				smooth.mat.tmp <- rbind(smooth.mat.tmp,smooth.mat)
				}
			if(RRl[k] > 0)
				{
				pmed <- psd <- pmean <- pqu2p5 <- pqu10 <- pqu90 <- pqu97p5 <- xnam <- rep(0,RRl[k])
				redrawsl <- redraw(desRl[[k]]$RQ,drawTZRl[[k]])
				lid <- 1:length(RA[[k]]$indspec[1,])
				for(j in 1:RRl[k])
					{
					lidok <- lid[RA[[k]]$indspec[1,] == j]
					pmed[j] <- median(redrawsl[j,])
					psd[j] <- sd(redrawsl[j,])
					pmean[j] <- mean(redrawsl[j,])
					pqu2p5[j] <- quantile(redrawsl[j,], probs = 0.025)
					pqu97p5[j] <- quantile(redrawsl[j,], probs = 0.975)
					pqu10[j] <- quantile(redrawsl[j,], probs = 0.1)
					pqu90[j] <- quantile(redrawsl[j,], probs = 0.9)
					xnam[j] <- RA[[k]][[lidok]]$term
					}
				betas <- cbind(pmean,psd,pqu2p5,pqu10,pmed,pqu90,pqu97p5)
				attr(betas,"linear.coef.draws") <- drawTZRl[[k]]
				attr(betas,"linear.coef.mean") <- drawTZRlM[[k]]
				attr(betas,"dimnames") <- list(xnam,c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5"))

				class(betas) <- "gibbsfit"

				for(j in 1:RRl[k])
					{
					lidok <- lid[RA[[k]]$indspec[1,] == j]
					fpmed <- (design[[k]][,j]*pmed[j])[RA[[k]]$fac]
					fpmean <- (design[[k]][,j]*pmean[j])[RA[[k]]$fac]
					fpqu2p5 <- (design[[k]][,j]*pqu2p5[j])[RA[[k]]$fac]
					fpqu10 <- (design[[k]][,j]*pqu10[j])[RA[[k]]$fac]
					fpqu90 <- (design[[k]][,j]*pqu90[j])[RA[[k]]$fac]
					fpqu97p5 <- (design[[k]][,j]*pqu97p5[j])[RA[[k]]$fac]
					xtmp <- RA[[k]][[lidok]]$x[RA[[k]]$fac]

					if(icheck)
						{
						fpmed <- fpmed[iZ]
						fpmean <- fpmean[iZ]
						fpqu2p5 <- fpqu2p5[iZ]
						fpqu10 <- fpqu10[iZ]
						fpqu90 <- fpqu90[iZ]
						fpqu97p5 <- fpqu97p5[iZ]
						xtmp <- xtmp[iZ]
						}
					partial.resid <- response - (eta - fpmed)
					pcat80 <- (fpqu10 < 0 & fpqu90 < 0)*(-1) + (fpqu10 <= 0 & fpqu90 >= 0)*0 + (fpqu10 > 0 & fpqu90 > 0)*1
					pcat95 <- (fpqu2p5 < 0 & fpqu97p5 < 0)*(-1) + (fpqu2p5 <= 0 & fpqu97p5 >= 0)*0 + (fpqu2p5 > 0 & fpqu97p5 > 0)*1
					fout <- cbind(xtmp,fpmean,fpqu2p5,fpqu10,fpmed,fpqu90,fpqu97p5,partial.resid,pcat95,pcat80)
					fout <- fout[order(fout[,1]),]
					colnames(fout) <- c(RA[[k]][[lidok]]$term,
							    "pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
					attr(fout,"linear.coef") <- betas[j,]
					attr(fout,"linear.ceffect") <- mean(cRLL[[k]][j,])
					attr(fout,"term.type") <- "linear"
					attr(fout,"linear.coef.draws") <- drawTZRl[[k]][j,]
					attr(fout,"linear.coef.draws.utr") <- redrawsl[j,]
					attr(fout,"linear.coef.mean") <- drawTZRlM[[k]][j,]
					class(fout) <- "linear.gibbs"

					ra.fout[[mapind[RA[[k]]$indspec[1,] == j]]] <- fout
					}
				attr(ra.fout,"lin.mat") <- betas
				betastmp <- rbind(betastmp,betas)
				}
			if(RRi[k] > 0)
				{
				if(which > 2)
					facind <- RA[[k]]$fac[iZ]
				else
					facind <- RA[[k]]$fac
				iNS <- extract.random(RRi[k],ins[[k]]$RR,ins[[k]]$RRl,ins[[k]]$RRs,ins[[k]]$RRi,
                                              	      ins[[k]]$fitR,ins[[k]]$RA,ins[[k]]$lambdaR,ins[[k]]$tildeZra,
                                              	      ins[[k]]$drawTildeGammaR,ins[[k]]$drawTildeGammaRM,ins[[k]]$dPR,response,eta,
                                              	      ins[[k]]$drawB,ins[[k]]$drawBM,ins[[k]]$S2R,ins[[k]]$ins,
                                              	      ins[[k]]$drawTZRl,ins[[k]]$drawTZRlM,ins[[k]]$desRl,ins[[k]]$design,
                                              	      ins[[k]]$metaR,ins[[k]]$newR,TRUE,facind,ins[[k]]$cRMS,ins[[k]]$cRSS,ins[[k]]$cRLL,
                                                      S2R,ins[[k]]$Z,which,data,ins[[k]]$stetaR)
				for(j in 1:RRi[k])
					{
					ra.fout[[mapind[RA[[k]]$indspec[3,] == j]]] <- iNS[[j]]
					betastmp <- rbind(betastmp,attr(iNS[[j]],"lin.mat"))
					smooth.mat.tmp <- rbind(smooth.mat.tmp,attr(iNS[[j]],"smooth.mat"))
					varmat <- rbind(varmat,attr(iNS[[j]],"var.mat"))
					}
				}
			refout[[k]] <- list(terms=ra.fout,effects=effects)
			attr(refout[[k]],"lin.mat") <- betastmp
			attr(refout[[k]],"smooth.mat") <- smooth.mat.tmp
			}
		else
			refout[[k]] <- list(terms=NULL,effects=effects)
		attr(refout[[k]],"term.type") <- "random"
		attr(refout[[k]],"var.mat") <- varmat
		class(refout[[k]]) <- "gibbs"
		}
	return(refout)	
	}

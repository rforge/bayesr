r.setup <- function(R, terms, ra, this, pframe, N, wcheck, 
                    hyperasigma, keepers, hyperatau, res, 
                    data, W, dcheck, speed) 
	{
   	newR <- newRsm <- fRtmp <- OTzRa <- tildeZra <- ttildeZra <- RA <- ZZ <- Z <- indspec <- 
      	fr <- IR <- tau2R <- drawBetaR <- lambdaR <- drawTildeGammaR <- drawTildeGammaRM <- fitR <- 
	drawB <- drawBM <- etaR <- metaR <- stetaR <- dPR <- tau2R <- hyperaTau2newR <- desRl <- IXRl <- tZRl <- ttZRl <- drawTZRl <- drawTZRlM <- 
     	newBR <- Ni <- ins <- termsi <- check <- la <- center <- frl <- design <- rind <- 
	cRM <- cRS <- cRMSS <- cRSS <- cRML <- cRLL <- rlm <- DRAWERS <- DRAWERSM <- nilor <- 
	DRAWERS1 <- DRAWERSM1 <- nilor1 <- DRAWERS1L <- DRAWERSM1L <- nilor1L <- nilore <- vector("list",R)
    	BR <- haRnew <- hbRnew <- RR <- hspec <- dPRl <- RRl <- RRs <- RRi <- cRL <- bcheck <- rcc <- rep(0, R)
  	rwcheck <- rep(FALSE,R)

      	R <- as.integer(R)
      	fR <- matrix(0, N, R)
      	f <- rep(0,N)
      	sigma2R <- rep(1e-04, R)
      	S2R <- matrix(0, R, keepers)
	cR <- rep(0,R)
	cRMS <- matrix(0, R, keepers)

    	for(k in 1:R) 
		{
		RA[[k]] <- smooth.construct(eval(parse(text = terms[ran[k]])),data)
		# RA[[k]] <- eval(parse(text = terms[ra[k]]), envir=data)

		obsz <- nrow(RA[[k]]$Z)
        	#if(N != obsz) 
		#	{
            	#	stoptxt <- paste("Variable length differing from N in term '", terms[ran[k]], "'!", sep = "")
            	#	stop(stoptxt, call. = FALSE)
        	#	}
        	if(wcheck) 
            		RA[[k]]$Z <- W*RA[[k]]$Z

        	BR[k] <- as.integer(RA[[k]]$n)
		DRAWERS[[k]] <- rep(0,BR[k])
		DRAWERSM[[k]] <- rep(0,BR[k])
		nilor[[k]] <- rep(0,BR[k])
		nilore[[k]] <- rep(0,BR[k])
        	haRnew[k] <- RA[[k]]$haRnew
		if(is.null(RA[[k]]$byfacid))
			{
			Z[[k]] <- list(tildeZ=RA[[k]]$Z,s=rep(1,BR[k]),RQ=diag(1,BR[k]))
			ZZ[[k]] <- diag(crossprod(RA[[k]]$Z))
			}
		else
			{
			bcheck[k] <- 1
			Z[[k]] <- trans.design(RA[[k]]$Z,RA[[k]]$K)
			ZZ[[k]] <- rep(1,BR[k])
			}
        	RR[k] <- RA[[k]]$length
		rind[[k]] <- RA[[k]]$fac
		if(speed)
          		rind[[k]] <- as.integer(rind[[k]] - 1)
		cRM[[k]] <- RA[[k]]$mx

        	if(RA[[k]]$hspec)
			{
			indspec[[k]] <- RA[[k]]$indspec
			RRl[k] <- as.integer(max(indspec[[k]][1,]))
			RRs[k] <- as.integer(max(indspec[[k]][2,]))
			RRi[k] <- as.integer(max(indspec[[k]][3,]))
			hspec[k] <- 1

			if(RRs[k] > 0)
				{
            			newR[[k]] <- newRsm[[k]] <- OTzRa[[k]] <- tildeZra[[k]] <- ttildeZra[[k]] <- drawTildeGammaR[[k]] <- 
				drawTildeGammaRM[[k]] <- IR[[k]] <- fr[[k]] <- cRMSS[[k]] <- cRSS[[k]] <- lambdaR[[k]] <- la[[k]] <- 
				DRAWERS1[[k]] <- DRAWERSM1[[k]] <- nilor1[[k]] <- fRtmp[[k]] <- vector("list",RRs[k])
            			newBR[[k]] <- dPR[[k]] <- center[[k]] <- cRS[[k]] <- rep(0, RRs[k])
				tau2R[[k]] <- hyperaTau2newR[[k]] <- rep(1e-04, RRs[k])
				}
			if(RRi[k] > 0)
				termsi[[k]] <- check[[k]] <- rep(0, RRi[k])
                	for(j in 1:RR[k]) 
				{
				if(RA[[k]][[j]]$check)
					{
					ii <- indspec[[k]][3,j]
					check[[k]][ii] <- 1
					termsi[[k]][ii] <- RA[[k]][[j]]$call
					Ni[[k]] <- as.integer(nrow(RA[[k]][[j]]$Z))
					}
				else
					{
                  			if(RA[[k]][[j]]$spec == 2) 
						{
						ii <- as.integer(indspec[[k]][2,j])

						if(is.null(RA[[k]][[j]]$xt$constr))
							{
                  					tildeZra[[k]][[ii]] <- RA[[k]][[j]]$tildeZ
							tildeZra[[k]][[ii]]$cs <- as.integer(0)
							}
						if(!is.null(RA[[k]][[j]]$xt$constr) && speed)
							{
							tildeZra[[k]][[ii]] <- list(tildeZ=RA[[k]][[j]]$X,
                                                                                    s=RA[[k]][[j]]$S[[1]],
                                                                                    RQ=diag(1,ncol(RA[[k]][[j]]$S[[1]])),
                                                                                    cs=as.integer(RA[[k]][[j]]$xt$constr),
                                                                                    ZZ=crossprod(RA[[k]][[j]]$X))
							}
						OTzRa[[k]][[ii]] <- tildeZra[[k]][[ii]]$tildeZ
						tildeZra[[k]][[ii]]$tildeZ <- tr(OTzRa[[k]][[ii]],RA[[k]][[j]]$ind)
						fRtmp[[k]][[ii]] <- rep(0,nrow(tildeZra[[k]][[ii]]$tildeZ))
                  				ttildeZra[[k]][[ii]] <- t(tildeZra[[k]][[ii]]$tildeZ)
						newR[[k]][[ii]] <- as.integer(RA[[k]][[j]]$ind)
						if(speed)
                					newRsm[[k]][[ii]] <- as.integer(newR[[k]][[ii]] - 1)
                  				dPR[[k]][ii] <- RA[[k]][[j]]$dPR
						DRAWERS1[[k]][[ii]] <- rep(0,dPR[[k]][ii])
						DRAWERSM1[[k]][[ii]] <- rep(0,dPR[[k]][ii])
						nilor1[[k]][[ii]] <- rep(0,dPR[[k]][ii])
						if(!is.matrix(tildeZra[[k]][[ii]]$tildeZ))
							tildeZra[[k]][[ii]]$tildeZ <- matrix(tildeZra[[k]][[ii]]$tildeZ,1,1)
						newBR[[k]][ii] <- as.integer(nrow(tildeZra[[k]][[ii]]$tildeZ))

						mrfcheck <- paste(strsplit(RA[[k]][[j]]$term,"")[[1]][1:3],collapse="")
						if(mrfcheck == "mrf")
							cRMSS[[k]][[ii]] <- list(a=c(RA[[k]][[j]]$mx),b=TRUE,c=tildeZra[[k]][[ii]]$RQ)
						else
							{
							if(is.null(RA[[k]][[j]]$xt$constr))
								cRMSS[[k]][[ii]] <- list(a=c(trans.design(RA[[k]][[j]]$mx)$tildeZ),b=FALSE)
							if(!is.null(RA[[k]][[j]]$xt$constr) && speed)
								cRMSS[[k]][[ii]] <- list(a=c(RA[[k]][[j]]$mx),b=FALSE)
							}
                    				fr[[k]][[ii]] <- RA[[k]][[j]]$fr
                    				drawTildeGammaR[[k]][[ii]] <- RA[[k]][[j]]$drawTildeGammaR
                    				drawTildeGammaRM[[k]][[ii]] <- RA[[k]][[j]]$drawTildeGammaR
						cRSS[[k]][[ii]] <- rep(0, keepers)
                  				IR[[k]][[ii]] <- RA[[k]][[j]]$IR
                    				hyperaTau2newR[[k]][ii] <- RA[[k]][[j]]$hyperaTau2newR

						if(!is.null(RA[[k]][[j]]$rwcheck))
							rwcheck[k] <- TRUE
						if(!is.null(RA[[k]][[j]]$sp))
							{
							if(is.null(RA[[k]][[j]]$rwcheck))
								la[[k]][[ii]] <- list(a=FALSE,b=RA[[k]][[j]]$sp,c=FALSE)
							else
								la[[k]][[ii]] <- list(a=FALSE,b=RA[[k]][[j]]$sp,c=TRUE)
							}
						else
							{
							if(is.null(RA[[k]][[j]]$rwcheck))
								la[[k]][[ii]] <- list(a=TRUE,b=RA[[k]][[j]]$sp,c=FALSE)
							else
								la[[k]][[ii]] <- list(a=TRUE,b=RA[[k]][[j]]$sp,c=TRUE)
							}
						if(is.null(RA[[k]][[j]]$xt$center))
							center[[k]][ii] <- 1
						else
							center[[k]][ii] <- RA[[k]][[j]]$xt$center
                      				lambdaR[[k]][[ii]] <- RA[[k]][[j]]$Rlambda
                 	 			}
					if(RA[[k]][[j]]$spec == 1)
						{
						ii <- as.integer(indspec[[k]][1,j])
						design[[k]] <- cbind(design[[k]],RA[[k]][[j]]$x)
						rlm[[k]] <- c(rlm[[k]],RA[[k]][[j]]$x[approxm(RA[[k]][[j]]$x)])
						}
					}
				}
			if(RRs[k] > 0)
        			{
				dPR[[k]] <- as.integer(dPR[[k]])
				center[[k]] <- as.integer(center[[k]])
				}
			if(RRl[k] > 0)
				{
				desRl[[k]] <- trans.design(design[[k]])
				tZRl[[k]] <- desRl[[k]]$tildeZ
				ttZRl[[k]] <- t(tZRl[[k]])
				dPRl[k] <- as.integer(ncol(design[[k]]))
				DRAWERS1L[[k]] <- rep(0,dPRl[k])
				DRAWERSM1L[[k]] <- rep(0,dPRl[k])
				nilor1L[[k]] <- rep(0,dPRl[k])
				IXRl[[k]] <- rep(1,dPRl[k])
				cRML[[k]] <- list(a=desRl[[k]]$RQ,b=rlm[[k]],c=length(rlm[[k]]),d=rep(1,length(rlm[[k]])))
				drawTZRl[[k]] <- matrix(0,dPRl[k],keepers)
				drawTZRlM[[k]] <- matrix(0,dPRl[k],keepers)
				frl[[k]] <- rep(0,BR[k])
				cRLL[[k]] <- matrix(0,nrow=length(rlm[[k]]),ncol=keepers)
				}
        		}

            	drawB[[k]] <- RA[[k]]$drawB
            	drawBM[[k]] <- RA[[k]]$drawB
            	stetaR[[k]] <- RA[[k]]$drawB
            	drawBetaR[[k]] <- RA[[k]]$drawBetaR
            	etaR[[k]] <- RA[[k]]$etaR
		metaR[[k]] <- RA[[k]]$etaR
		rcc[k] <- RA[[k]]$center
		}

	for(k in 1:R)
		if(RRi[k] > 0)
			{
			rai <- 1:RRi[k]
			ins[[k]] <- r.setup(RRi[k], termsi[[k]], rai, this, pframe, Ni[[k]], wcheck, 
                                	    hyperasigma, keepers, hyperatau, res, RA[[k]]$usedata, W, dcheck, speed) 
			}

	return(list(fR=fR,RA=RA,sigma2R=sigma2R,ZZ=ZZ,metaR=metaR,newR=newR,cRM=cRM,
		      	Z=Z,fr=fr,BR=as.integer(BR),drawBetaR=drawBetaR,center=center,newBR=newBR,
           		tau2R=tau2R,tildeZra=tildeZra,dPR=dPR,ttildeZra=ttildeZra,cRMS=cRMS,cRMSS=cRMSS,cRS=cRS,cRSS=cRSS,
                  	IR=IR,hyperaTau2newR=hyperaTau2newR,hbRnew=hbRnew,OTzRa=OTzRa, 
                  	haRnew=haRnew,lambdaR=lambdaR,drawTildeGammaR=drawTildeGammaR, 
                 	drawB=drawB,S2R=S2R,indspec=indspec,cR=cR,cRL=cRL,cRML=cRML,cRLL=cRLL, 
                  	res=res,RR=RR,etaR=etaR,ins=ins,design=design,rind=rind,rwcheck=rwcheck,
                  	hspec=as.integer(hspec),check=check,R=R,la=la,RRs=as.integer(RRs),RRl=as.integer(RRl),RRi=as.integer(RRi),
                  	desRl=desRl,dPRl=as.integer(dPRl),IXRl=IXRl,tZRl=tZRl,ttZRl=ttZRl,drawTZRl=drawTZRl,frl=frl,
                  	drawTZRlM=drawTZRlM,drawTildeGammaRM=drawTildeGammaRM,drawBM=drawBM,
                  	DRAWERS=DRAWERS,DRAWERSM=DRAWERSM,nilor=nilor,DRAWERS1=DRAWERS1,DRAWERSM1=DRAWERSM1,
                  	nilor1=nilor1,DRAWERS1L=DRAWERS1L,DRAWERSM1L=DRAWERSM1L,nilor1L=nilor1L,f=f,fRtmp=fRtmp,newRsm=newRsm,nilore=nilore,
			bcheck=as.integer(bcheck),stetaR=stetaR,rcc=as.integer(rcc)))
	}

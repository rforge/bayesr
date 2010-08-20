gibbs <- function(formula, data, weights, family = "gaussian", iter = 1200, burnin = 200, thinning = 1, 
		  hyperasigma = 1e-04, hyperbsigma = 1e-04, hyperatau = 1e-04, hyperbtau = 1e-04, 
		  trace = TRUE, dots = FALSE, speed = TRUE, ...)
	{
	res <- list()
	res$call <- match.call()
	res$formula <- formula <- as.formula(formula)
	varNames <- all.vars(formula)
	pframe <- parent.frame()
	lvn <- length(varNames)
	tf <- terms(formula, specials = c("s","te","m","r"))
	itcpt <- attr(tf,"intercept")

	terms <- attr(tf, "term.labels")
	vars <- as.character(attr(tf, "variables"))
	vars <- vars[vars != "list"]
	res$vars <- vars
	res$terms <- terms
	res$tf <- tf
	lt <- length(terms)
	
	if(attr(tf, "response") < 1)
		stop("No response variable specified!")

	if(missing(data))
		{
		dcheck <- FALSE
		data <- .GlobalEnv
		}
	else
		dcheck <- TRUE
	this <- .GlobalEnv

	response <- eval(parse(text = vars[1]), envir=data)

	res$N <- N <- length(response)
	if(N < 2)
		stop("Not enough observations!")
	L <- K <- R <- M <- 0

	if(burnin > iter)
		{
		warning("Burnin exceeds number of iterations, burnin is set to default!")
		burnin <- 0
		}
	if(thinning > (iter-burnin) || thinning == (iter-burnin))
		if(thinning != 1)
			{
			warning("Thinning argument is specified wrong, thinning is set to default!")
			thinning <- 1
			}
	if(thinning <= 0)
		{
		warning("Thinning argument is set below zero, thinning is set to default!")
		thinning <- 1
		}
	keepers <- round((iter-burnin)/thinning)
	res$keepers <- keepers

	wcheck <- !missing(weights)
	if(wcheck)
		{
		if(length(weights) != N)
			stop("Lenght of weights is differing from number of observations!")
		W <- 1/sqrt(weights)
		response <- W*response
		}
	else 
		W <- rep(1,N)

	cat("Setting up model...\n")
	tnames <- NULL
	smstuff <- list(NULL)
	mustuff <- list(NULL)
	rastuff <- list(NULL)
	ALL <- 0
	if(lt > 0)
		{
		ind <- c(1:lt)
		sm <- attr(tf,"specials")$s-1
		te <- attr(tf,"specials")$te-1
		mul <- attr(tf,"specials")$m-1
		ran <- attr(tf,"specials")$r-1
		lin <- ind[!ind %in% c(sm,te,mul,ran)]
		smooth <- ind[!ind %in% c(mul,ran,lin)]

		lsm <- length(sm)
		lte <- length(te)
		M <- length(mul)
		R <- length(ran)
		L <- length(lin)
		K <- lsm+lte
		KS <- MS <- RS <- mfun <- mfunid <- NULL
		ALL <- M+R+L+K
		byfun <- FALSE

		if(K > 0)
			{
			# smooth function setup
			tildeZ <- OldK <- cK <- ttildeZ <- new <- newsm <- cmcheck <- sm.specs <- vector("list", K)
			Rank <- cme <- lbdcheck <- NULL
			jj <- 1
			k <- minsm <- 0
			kk <- K
			while(jj < (K+1))
				{
				mfuncheck <- FALSE
				smtmp <- eval(parse(text = terms[smooth[jj-minsm]]), envir=data)
				if(class(smtmp)!="list")
					smtmp <- stewrap(smtmp,data,terms[smooth[jj-minsm]])	
				if(!attr(smtmp,"byfun"))
					smtmp <- list(smtmp)
				else
					{
					if(attr(smtmp,"type") == "m")
						{
						terms[smooth[jj]] <- smtmp
						if(M > 0)
							mul <- c(mul,smooth[jj])
						else
							mul <- smooth[jj]
						mfuncheck <- TRUE
						smooth <- smooth[smooth!=smooth[jj]]
						kk <- kk - 1
						jj <- jj + 1
						minsm <- minsm + 1
						}
					}
				if(!mfuncheck)
					{
					nk <- length(smtmp)
					kk <- kk + nk - 1
					for(j in 1:nk)
						{
						if(j > 1)
							{
							byfun <- TRUE
							smooth <- c(smooth,length(terms)+1)
							terms <- c(terms,smtmp[[j]]$spec$label)
							}
						k <- k + 1
						if(N!=length(smtmp[[j]]$ind))
							{
							stoptxt <- paste("Variable length differing from N in term '",terms[smooth[k]],"'!",sep="")
							stop(stoptxt, call.=FALSE)
							}
						if(wcheck)
							smtmp[[j]]$X <- weights*smtmp[[j]]$X

						if(smtmp[[j]]$center)
							cme <- c(cme,1)
						else
							cme <- c(cme,0)
						if(is.null(smtmp[[j]]$constr))
							{
							tildeZ[[k]] <- trans.design(smtmp[[j]]$X,smtmp[[j]]$S[[1]])
							tildeZ[[k]]$cs <- as.integer(0)
							OldK[[k]] <- tildeZ[[k]]$tildeZ
							cK[[k]] <- c(trans.design(smtmp[[j]]$mx)$tildeZ)
							}
						else
							{
							tildeZ[[k]] <- list(tildeZ=smtmp[[j]]$X,s=smtmp[[j]]$S[[1]],RQ=diag(1,ncol(smtmp[[j]]$S[[1]])),
                                                                    	    cs=as.integer(smtmp[[j]]$constr),ZZ=crossprod(smtmp[[j]]$X))
							OldK[[k]] <- smtmp[[j]]$X
							cK[[k]] <- smtmp[[j]]$mx
							}
						tildeZ[[k]]$tildeZ <- tr(tildeZ[[k]]$tildeZ,smtmp[[j]]$ind)
						tildeZ[[k]]$lambda <- smtmp[[j]]$sp
						if(!is.null(tildeZ[[k]]$lambda))
							lbdcheck <- c(lbdcheck,1)
						else
							lbdcheck <- c(lbdcheck,0)
						ttildeZ[[k]] <- t(tildeZ[[k]]$tildeZ)
						new[[k]] <- smtmp[[j]]$ind
						if(speed)
            						newsm[[k]] <- new[[k]] - 1
						if(class(smtmp[[j]]) == "mrf.smooth")
							cmcheck[[k]] <- list(a=TRUE,b=tildeZ[[k]]$RQ)
						else
							cmcheck[[k]] <- list(a=FALSE)
						Rank <- c(Rank,qr(smtmp[[j]]$S[[1]])$rank)
						sm.specs[[k]] <- smtmp[[j]]$spec
						attr(sm.specs[[k]],"smooth.construct") <- smtmp[[j]]
						}
					jj <- jj + 1
					}
				}
			K <- kk
			if(K > 0)
				{
				dP <- Nnew <- rep(0, (K))
				I <- vector("list", K)

				fps <- matrix(0,N,K)
				csK <- matrix(0,K,keepers)
				sce <- rep(0,K)
			
				drawTildeGamma <- mDTG <- drawsmC <- vector("list", K)				
				tau2 <- rep(0.0001,K)
				lambda <- matrix(0,K,keepers)	
				for(k in 1:K)
					{
					dP[k] <- ncol(tildeZ[[k]]$tildeZ)
					Nnew[k] <- nrow(tildeZ[[k]]$tildeZ)
					I[[k]] <- rep(1,dP[k])
					drawTildeGamma[[k]] <- matrix(0,dP[k],keepers)
					mDTG[[k]] <- matrix(0,dP[k],keepers)
					drawsmC[[k]] <- rep(0,dP[k])
					}
				hyperaTau2new <- hyperatau + (Rank)/2
			
				if(speed)
          				{
          				smstuff <- list(fps,tildeZ,tau2,OldK,ttildeZ,hyperbtau,hyperaTau2new,
                          				lambda,as.integer(dP),drawTildeGamma,mDTG,
                          				as.integer(cme),cmcheck,as.numeric(sce),cK,csK,newsm,as.integer(lbdcheck))
          				rm(fps,tau2,OldK,ttildeZ,hyperaTau2new,lbdcheck,
             				lambda,drawTildeGamma,mDTG,cme,cmcheck,sce,cK,csK,newsm)
          				invisible(gc())
          				}
				if(byfun)
					res$terms <- terms
				}
			M <- length(mul)
			}
		if(M > 0)
			{
			# some multilevel fun
			tildeZmu <- ttildeZmu <- MU <- muNames <- IM <- fpsM <- flM <- tau2M <- hyperaTau2newM <- vector("list", M)
			dPM <- RankM <- fitM <- ind <- cM <- cMM <- csM <- mult.specs <- drawsmCmu <- vector("list", M)
			nl <- vector(length=M)
			nmu <- nsm <- nlm <- 0
			hypident <- rep(FALSE,M)

			fM <- matrix(0,N,M)

			spec <- vector(length=M)
			lambdaM <- drawTildeGammaM <- FM <- mlambdacheck <- mlcheck <- mcentering <- muz <- vector("list", M)
			
			for(k in 1:M)
				{
				mutmp <- eval(parse(text = terms[mul[k]]), envir=data)
				hypident[k] <- mutmp$hypident
				if(N!=mutmp$N)
					stop("Variable lengths differing in m() term")
				faclev <- levels(as.factor(mutmp$fac))
				spec[k] <- mutmp$smcheck
				nl[k] <- mutmp$nl
				ind[[k]] <- mult.specs[[k]] <- vector("list", nl[k])				

				if(spec[k] == 2)
					{
					nsm <- nsm + 1
					dPM[[k]] <- as.integer(rep(0,nl[k]))
					if(speed)
              					mlambdacheck[[k]] <- as.numeric(rep(0,nl[k]))
					else
              					mlambdacheck[[k]] <- rep(NA,nl[k])
          				mcentering[[k]] <- rep(NA,nl[k])
					mlcheck[[k]] <- rep(0,nl[k])
					RankM[[k]] <- vector(length = nl[k])
					lambdaM[[k]] <- matrix(0,nl[k],keepers)
					if(hypident[k])
						tau2M[[k]] <- 0.0001
					else	
						tau2M[[k]] <- rep(0.0001,(nl[k]))
					csM[[k]] <- matrix(0,nl[k],keepers)
					cM[[k]] <- rep(0,nl[k])
					FM[[k]] <- drawTildeGammaM[[k]] <- vector("list", nl[k])						
					tildeZmu[[k]] <- ttildeZmu[[k]] <- IM[[k]] <- vector("list", nl[k])
					fpsM[[k]] <- cMM[[k]] <- drawsmCmu[[k]] <- vector("list", nl[k])
					tnames <- NULL

					for(j in 1:nl[k])
						{
						ind[[k]][[j]] <- as.integer(mutmp[[j]]$ind)
						if(wcheck)
							mutmp[[j]]$X <- W[ind[[k]][[j]]]*mutmp[[j]]$X
						if(is.null(mutmp[[j]]$xt$constr))
							{
							tildeZmu[[k]][[j]] <- trans.design(mutmp[[j]]$X,mutmp[[j]]$S[[1]])
							tildeZmu[[k]][[j]]$cs <- as.integer(0)
							}
						else
							tildeZmu[[k]][[j]] <- list(tildeZ=mutmp[[j]]$X,s=mutmp[[j]]$S[[1]],
                                                                                   RQ=diag(1,ncol(mutmp[[j]]$S[[1]])),
                                                                                   cs=as.integer(mutmp[[j]]$xt$constr),
                                                                                   ZZ=crossprod(mutmp[[j]]$X))
						ttildeZmu[[k]][[j]] <- t(tildeZmu[[k]][[j]]$tildeZ) 

						if(mutmp[[j]]$ismrf)
							{
							cMM[[k]][[j]] <- list(a=mutmp[[j]]$ismrf,b=c(mutmp[[j]]$mx),c=tildeZmu[[k]][[j]]$RQ)
							smx <- mutmp[[j]]$x
							mutmp$names <- mutmp$term[1]
							}
						else
							{
							if(is.null(mutmp[[j]]$constr))
								cMM[[k]][[j]] <- list(a=mutmp[[j]]$ismrf,b=c(trans.design(mutmp[[j]]$mx)$tildeZ))
							else
								cMM[[k]][[j]] <- list(a=mutmp[[j]]$ismrf,b=mutmp[[j]]$mx)
							smx <- mutmp[[j]]$x
							}
						if(!is.null(mutmp[[j]]$lambda))
              						{
							mlambdacheck[[k]][j] <- mutmp[[j]]$lambda[j]
                  					mlcheck[[k]][j] <- 1
							}             
						mcentering[[k]][j] <- mutmp[[j]]$center
						dPM[[k]][j] <- as.integer(ncol(tildeZmu[[k]][[j]]$tildeZ))
						RankM[[k]][j] <- qr(mutmp[[j]]$S[[1]])$rank				
						IM[[k]][[j]] <- rep(1,dPM[[k]][j])
						drawsmCmu[[k]][[j]] <- rep(0,dPM[[k]][j])
						nm <- mutmp[[j]]$n
						nmu <- c(nmu,nm)
						fpsM[[k]][[j]] <- rep(0,nm)
						if(j > 1)
							tildeZmu[[k]][[j]]$intcpt <- rep(1,nm)
						drawTildeGammaM[[k]][[j]] <- matrix(0,dPM[[k]][j],keepers)

						if(length(mutmp$term) == 1)
							{
							xnam <- paste(mutmp$term,":",faclev[j],sep="")
							tnames <- c(tnames,xnam)
							}
						else
							{
							xnam1 <- paste(mutmp$term[1],":",faclev[j],sep="")
							xnam2 <- paste(mutmp$term[2],":",faclev[j],sep="")
							xnam <- c(xnam1,xnam2)
							xnams <- paste(mutmp$term[1],",",mutmp$term[1],":",faclev[j],sep="")
							tnames <- rbind(tnames,xnams)
							}
						mult.specs[[k]][[j]] <- mutmp[[j]]$specs
						}
					if(speed)
              					mlambdacheck[[k]] <- as.numeric(mlambdacheck[[k]])
          				mlcheck[[k]] <- as.integer(mlcheck[[k]])  
					attr(mult.specs[[k]],"m.term") <- terms[mul[k]]
					attr(mult.specs[[k]],"tnames") <- tnames
					attr(mult.specs[[k]],"data") <- mutmp$data
					hyperaTau2newM[[k]] <- hyperatau + (RankM[[k]])/2
					if(hypident[k])
						hyperaTau2newM[[k]] <- sum(hyperaTau2newM[[k]])
					}
				if(spec[k] == 1)
					{
					nlm <- nlm + 1
					dPM[[k]] <- rep(0,nl[k])	
					tildeZmu[[k]] <- ttildeZmu[[k]] <- IM[[k]] <- vector("list", nl[k])
					flM[[k]] <- vector("list", nl[k])
					FM[[k]] <- drawTildeGammaM[[k]] <- vector("list", nl[k])
					sigma2M <- 0.0001
					csM[[k]] <- vector("list",nl[k])
					tnames <- NULL

					for(j in 1:nl[k])
						{
						ind[[k]][[j]] <- as.integer(mutmp[[j]]$ind)
						if(wcheck)
							mutmp[[j]]$X <- W[ind[[k]][[j]]]*mutmp[[j]]$X
						tildeZmu[[k]][[j]] <- trans.design(mutmp[[j]]$X,mutmp[[j]]$S[[1]])
						ttildeZmu[[k]][[j]] <- t(tildeZmu[[k]][[j]]$tildeZ)
						dPM[[k]][j] <- ncol(tildeZmu[[k]][[j]]$tildeZ)
						IM[[k]][[j]] <- rep(1,dPM[[k]][j])
						nm <- mutmp[[j]]$n
						nmu <- c(nmu,nm)
						flM[[k]][[j]] <- rep(0,nm)
						cMM[[k]][[j]] <- list(a=tildeZmu[[k]][[j]]$RQ,b=rep(1,length(mutmp[[j]]$mx)),
                                                                      c=length(mutmp[[j]]$mx),d=mutmp[[j]]$mx)
						drawTildeGammaM[[k]][[j]] <- matrix(0,dPM[[k]][j],keepers)
						cM[[k]] <- rep(0,nl[k])
						csM[[k]][[j]] <- matrix(0,length(mutmp[[j]]$mx),keepers)

						xnam <- paste(mutmp$specs$term,".",faclev[j],sep="")
						mult.specs[[k]][[j]] <- list(x=mutmp[[j]]$X,tnames=xnam,
					                     		     term=terms[mul[k]],type="linear")
						tnames <- c(tnames,xnam)
						}
					attr(mult.specs[[k]],"tnames") <- tnames
					}
				}
      
			drawTGMM <- drawTildeGammaM
			if(speed)
          			{
          			mustuff <- list(fM,as.integer(spec),hypident,tau2M,as.integer(nl),mlambdacheck,tildeZmu,
                          			ind,dPM,cM,cMM,fpsM,mcentering,hyperaTau2newM,lambdaM,
                          			drawTildeGammaM,drawTGMM,csM,ttildeZmu,mlcheck)
          			rm(fM,hypident,tau2M,mlambdacheck,cM,cMM,fpsM,mcentering,hyperaTau2newM,lambdaM,
             		           drawTildeGammaM,drawTGMM,csM,ttildeZmu,mlcheck)   
          			invisible(gc())                                         
          			}        
			}
		if(R > 0)
			{
			# random effects setup 
   			newR <- fRtmp <- OTzRa <- tildeZra <- ttildeZra <- RA <- ZZ <- Z <- indspec <- 
      			fr <- IR <- tau2R <- drawBetaR <- lambdaR <- drawTildeGammaR <- drawTildeGammaRM <- fitR <- 
			drawB <- drawBM <- etaR <- metaR <- stetaR <- dPR <- tau2R <- hyperaTau2newR <- desRl <- IXRl <- tZRl <- ttZRl <- drawTZRl <- drawTZRlM <- 
     		 	newBR <- Ni <- ins <- termsi <- check <- la <- center <- frl <- design <- rind <- 
			cRM <- cRS <- cRMSS <- cRSS <- cRML <- cRLL <- rlm <- newRsm <- vector("list",R)
    			BR <- haRnew <- hbRnew <- RR <- hspec <- dPRl <- RRl <- RRs <- RRi <- cRL <- bcheck <- rcc <- rep(0, R)
			rwcheck <- rep(FALSE,R)
          
            		R <- as.integer(R)
        		fR <- matrix(0, N, R)
        		sigma2R <- rep(1e-04, R)
        		S2R <- matrix(0, R, keepers)
			cR <- rep(0,R)
			cRMS <- matrix(0, R, keepers)

    			for(k in 1:R) 
				{
				RA[[k]] <- eval(parse(text = terms[ran[k]]), envir=data)
				obsz <- nrow(RA[[k]]$Z)
        			if(N != obsz) 
					{
            				stoptxt <- paste("Variable length differing from N in term '", terms[ran[k]], "'!", sep = "")
            				stop(stoptxt, call. = FALSE)
        				}
        			if(wcheck) 
            				RA[[k]]$Z <- W*RA[[k]]$Z
        			BR[k] <- RA[[k]]$n
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
					RRl[k] <- max(indspec[[k]][1,])
					RRs[k] <- max(indspec[[k]][2,])
					RRi[k] <- max(indspec[[k]][3,])
					hspec[k] <- 1

					if(RRs[k] > 0)
						{
            					newR[[k]] <- newRsm[[k]]<- OTzRa[[k]] <- tildeZra[[k]] <- ttildeZra[[k]] <- drawTildeGammaR[[k]] <- 
						drawTildeGammaRM[[k]] <- IR[[k]] <- fr[[k]] <- cRMSS[[k]] <- cRSS[[k]] <- lambdaR[[k]] <- la[[k]] <- 
						fRtmp[[k]] <- vector("list",RRs[k])
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
							Ni[[k]] <- nrow(RA[[k]][[j]]$Z)
							}
						else
							{
                  					if(RA[[k]][[j]]$spec == 2) 
								{
								ii <- indspec[[k]][2,j]
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
								newR[[k]][[ii]] <- RA[[k]][[j]]$ind
                						if(speed)
                    							newRsm[[k]][[ii]] <- as.integer(newR[[k]][[ii]] - 1)
                  						dPR[[k]][ii] <- RA[[k]][[j]]$dPR
								if(!is.matrix(tildeZra[[k]][[ii]]$tildeZ))
									tildeZra[[k]][[ii]]$tildeZ <- matrix(tildeZra[[k]][[ii]]$tildeZ,1,1)
								newBR[[k]][ii] <- nrow(tildeZra[[k]][[ii]]$tildeZ)

								mrfcheck <- paste(strsplit(RA[[k]][[j]]$term,"")[[1]][1:3],collapse="")
								if(class(RA[[k]][[j]]) == "mrf.smooth")
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
								if(!is.null(RA[[k]][[j]]$lambda))
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
								ii <- indspec[[k]][1,j]
								design[[k]] <- cbind(design[[k]],RA[[k]][[j]]$x)
								rlm[[k]] <- c(rlm[[k]],RA[[k]][[j]]$x[approxm(RA[[k]][[j]]$x)])
								}
							}
						}
					if(RRl[k] > 0)
						{
						desRl[[k]] <- trans.design(design[[k]])
						tZRl[[k]] <- desRl[[k]]$tildeZ
						ttZRl[[k]] <- t(tZRl[[k]])
						dPRl[k] <- ncol(design[[k]])
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
			
			if(speed)
          			{
          			rastuff <- list(fR,hspec,rwcheck,sigma2R,Z,etaR,BR,drawBetaR,cR,cRM,RRi,ins,
                          			metaR,RRs,fr,la,tau2R,tildeZra,OTzRa,dPR,ttildeZra,center,
                          			hyperaTau2newR,cRS,cRMSS,lambdaR,drawTildeGammaR,drawTildeGammaRM,
                          			cRSS,newRsm,RRl,frl,tZRl,ttZRl,cRL,cRML,drawTZRl,drawTZRlM,cRLL,
                          			haRnew,drawB,drawBM,S2R,cRMS,dPRl,fRtmp,rind,ZZ,bcheck,stetaR,as.integer(rcc))
                          
          			rm(fR,hspec,rwcheck,sigma2R,etaR,BR,drawBetaR,cR,cRM,ins,
            			   metaR,fr,la,tau2R,OTzRa,dPR,ttildeZra,center,stetaR,
             			   hyperaTau2newR,cRS,cRMSS,lambdaR,drawTildeGammaR,drawTildeGammaRM,
             			   frl,tZRl,ttZRl,cRL,cRML,drawTZRl,drawTZRlM,cRLL,
             			   haRnew,drawB,drawBM,S2R,cRMS,dPRl,fRtmp,newRsm,ZZ,bcheck,rcc)
          			invisible(gc())
          			}
			}

		tnames <- NULL
		if(L > 0)
			{
			# setup linear stuff
			LIN <- linNames <- linposition <- vector("list", L)
			types <- NULL
			kk <- 1
			for(k in 1:L)
				{
				tmp <- eval(parse(text = paste("model.matrix(~",terms[lin[k]],")")), envir=data)[,2]
				#tmp <- eval(parse(text = terms[lin[k]]), envir=data)

				if(N!=length(tmp))
					stop("Variable lengths differing in linear terms!",call.=FALSE)
				for(i in 1:lvn)
					if(any(grep(varNames[i], terms[lin[k]])))
						linNames[[k]] <- terms[lin[k]]
				if(is.factor(tmp))
					{
					nt <- levels(tmp)
					nlinlev <- nlevels(tmp)
					tmp <- diag(nlinlev)[tmp,]
					colnames(tmp) <- nt
					cs <- colSums(tmp)
					ind.l <- 1:length(cs)
					which <- ind.l[cs==max(cs)]
					if(length(which)>1)
						which <- which[1]
					ind.l <- ind.l[ind.l!=which]
					tmp <- tmp[,ind.l]
					if(is.matrix(tmp))
						nt <- colnames(tmp)
					else
						nt <- nt[ind.l]
					for(j in 1:length(nt))
						nt[j] <- paste(linNames[[k]],":",nt[j],sep="")
					tnames <- c(tnames,nt)
					types <- c(types,rep("factor",length(nt)))
					LIN[[k]] <- tmp
					linposition[[k]] <- c(1:length(ind.l))+kk-1
					attr(linposition[[k]],"factorcheck") <- TRUE
					kk <- kk + length(ind.l)
					}
				else
					{
					tnames <- c(tnames,linNames[[k]])
					types <- c(types,"numeric")
					LIN[[k]] <- tmp
					linposition[[k]] <- kk
					attr(linposition[[k]],"factorcheck") <- FALSE
					kk <- kk + 1
					}
				}
			}
		}
	if(ALL < 1 & itcpt < 1)
		stop("Nothing to do for gibbs()!")
	if(ALL == 1)
		ALL <- 2
	res$K <- K
	res$L <- L
	res$R <- R
	res$M <- M

	# Xlin <- matrix(0,N,L)	
	Xlin <- NULL
	if(itcpt > 0)				
		Xdesign <- matrix(1,N,1)
	else
		Xdesign <- matrix(0,N,1)
	Lorig <- L	
	if(L > 0)
		{
		for(k in 1:L)
			Xlin <- cbind(Xlin,LIN[[k]])
			#Xlin[,k] <- LIN[[k]]
		Ltmp <- ncol(Xlin)
		if(Ltmp > L)
			{	
			tmp1 <- tmp2 <- vector("list",Ltmp)
			for(k in 1:Ltmp)
				{
				tmp1[[k]] <- Xlin[,k]
				tmp <- tnames[k]
				attr(tmp,"type") <- types[k]
				tmp2[[k]] <- tmp
				}
			res$L <- L <- Ltmp
			LIN <- tmp1
			}
		}
    	Xdesign <- cbind(Xdesign,Xlin)
	if(sum(Xdesign)==0)
		itcpt <- 0

	#if(M > 0)
	#	{
	#	for(k in 1:M)
	#		if(spec[k] == 2)
	#			{
	#			Xdesign <- cbind(Xdesign,muz[[k]]$Z[,2:nl[k]])
	#			tnames <- c(tnames,attr(mult.specs[[k]],"tnames")[2:nl[k]])
	#			}
	#	}
	if(wcheck)
		Xdesign <- W*Xdesign

	bm <- rep(0,ncol(Xdesign))
	bm[1] <- 1
	if(ncol(Xdesign) > 1)
		{
		for(j in 2:ncol(Xdesign))
			{
			tmp <- approxm(Xdesign[,j])
			bm[j] <- Xdesign[tmp,j]
			}
		}
	LL <- length(bm)
	if(itcpt < 1 & L < 1)
		tildeX <- list(tildeZ=Xdesign,RQ=1,s=1)
	else
		tildeX <- c(trans.design(Xdesign))
	if(itcpt < 1)
		{
		if(is.matrix(tildeX$RQ))
			tildeX$RQ[1,1] <- 0.0
		else
			tildeX$RQ <- 0.0
		}
	ttildeX <- t(tildeX$tildeZ)
	tbrq <- t(tildeX$RQ)
	nclin <- ncol(Xdesign)
	Ixx <- rep(1,nclin)

	# Initialize vectors/matrices
	flin <- rep(0,N)
	ce <- 0
	bs <- 0
	cL <- matrix(0,(LL-1),keepers)
	CE <- rep(0,keepers)
	drawTildeBeta <- matrix(0,nclin,keepers)
	drawTBML <- matrix(0,nclin,keepers)
	S2 <- rep(0,keepers)
	drawTbeta <- rep(0,nclin)
	sigma2 <- 0.0001 
	eta <- rep(0,N)
	meta <- rep(0,N)
	diceta <- 0

	# Hyperparameters for gamma distributions
	hyperaSigma2new <- hyperasigma + N/2

	# Data augmentation preparation
	LO <- HI <- NULL
	if(family == "binomial")
		{
		origresponse <- as.integer(response)
		lo <- c(-Inf, 0)
		hi <- c(0, Inf)
		LO <- lo[response + 1]
		HI <- hi[response + 1]

		prob <- TRUE
		sigma2 <- 1
		}
	else
		{
		prob  <- FALSE
		origresponse <- response
		# sdy <- sd(response)
		# sdy2 <- sdy^2
		# response <- response*sdy
		}

	itrthin <- 1
	steps <- 1000
	thinsteps <- floor(seq(burnin,iter,length=(keepers)))
	thinsteps <- thinsteps[thinsteps!=0]
	maxmoniter <- 100

	#######################################################################################################################
      	# Sampler                                                                                                             #
	#######################################################################################################################

	ptm <- proc.time()
	if(speed)
		{
		out <- sampling(N,
                                iter,
                    		thinsteps,
                    		response, 
                    		sigma2,                   
                    		K,
                    		eta,
                    		flin,
                    		Ixx,
                    		tildeX,
                    		ttildeX,
                    		hyperbsigma,
                    		hyperaSigma2new,
                    		S2,
                    		drawTildeBeta,
                    		drawTBML,
                    		LL,
                    		tbrq,
                    		bm,
                    		cL,
                    		CE,
                    		smstuff,
                    		mustuff,
                    		M,
                    		rastuff,
                    		R,
                    		prob,
                    		LO,
                    		HI,
				itcpt,
				ALL)
		now <- proc.time()
    		rm(mustuff,smstuff,rastuff)
    		invisible(gc()) 
    		# smstuff
    		if(K > 0)
        		{
        		lambda <- out$outsm$lambda
        		drawTildeGamma <- out$outsm$drawtg
        		mDTG <- out$outsm$mdrawtg
        		csK <- out$outsm$csK
        		}
    		if(M > 0)
        		{
        		lambdaM <- out$outmu$lambdaM
        		drawTildeGammaM <- out$outmu$drawTildeGammaM
        		drawTGMM <- out$outmu$drawTGMM
        		csM <- out$outmu$csM
       		}
    		if(R > 0)
        		{
        		drawB <- out$outra$drawB
        		drawBM <- out$outra$drawBM
        		S2R <- out$outra$S2R
        		cRMS <- out$outra$cRMS   
        		drawTildeGammaR <- out$outra$drawTildeGammaR 
        		drawTildeGammaRM <- out$outra$drawTildeGammaRM
        		lambdaR <- out$outra$lambdaR
        		cRSS <- out$outra$cRSS  
        		drawTZRl <- out$outra$drawTZRl
        		drawTZRlM <- out$outra$drawTZRlM  
        		cRLL <- out$outra$cRLL 
        		metaR <- out$outra$metaR 
        		ins <- out$outra$ins
			stetaR <- out$outra$stetaR
        		}
		S2 <- out$S2
		drawTildeBeta <- out$lcoef
		drawTBML <- out$mlcoef
		drawTbeta <- out$ldraws
		meta <- out$meta
		diceta <- out$diceta
		itrthin <- out$itrthin + 1
		cL <- out$cL
		CE <- out$CE
		}
	else{
  	cat("Starting\n")
	# Iteration
	for(i in 1:iter)
		{
		# Data augmentation
		sigma <- sigma2^0.5

		if(prob)
			response <- rtruncated(N,LO,HI,pnorm,qnorm,eta,sigma)
	
		# Sampling for nonlinear functions
		if(K > 0)
			for(k in 1:K)
				{
				eta <- eta - fps[,k]

				if(is.null(tildeZ[[k]]$lambda))
					lambd <- sigma2/tau2[k]
				else
					lambd <- tildeZ[[k]]$lambda
			
				if(tildeZ[[k]]$cs > 0)
					{
					help <- tildeZ[[k]]$ZZ + lambd*tildeZ[[k]]$s
					mutmp <- crossprod(tildeZ[[k]]$tildeZ,(response - eta))

					for(j in 1:20)
						{
						for(d in 1:dP[k])
							{
							mcsmu <- 0
							if(d > 1)
								mcsmu <- sum(drawsmC[[k]][1:(d-1)]*help[d,1:(d-1)])
							if(d < dP[k])
								mcsmu <- mcsmu + sum(drawsmC[[k]][(d+1):dP[k]]*help[d,(d+1):dP[k]])

							mcsmu <- (mutmp[d] - mcsmu)/help[d,d]

							if(tildeZ[[k]]$cs==1)
								{
								if(d==1)
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=-Inf,high=drawsmC[[k]][2])
								if(d==dP[k])
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmC[[k]][d-1],high=Inf)
								if(d > 1 && d < dP[k])
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmC[[k]][d-1],high=drawsmC[[k]][d+1])
								}
							else
								{
								if(d==1)
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmC[[k]][2],high=Inf)
								if(d==dP[k])
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=-Inf,high=drawsmC[[k]][d-1])
								if(d > 1 && d < dP[k])
									drawsmC[[k]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmC[[k]][d+1],high=drawsmC[[k]][d-1])
								}
							}
						}

					draws <- drawsmC[[k]] 
					sums <- crossprod(draws,tildeZ[[k]]$s)%*%draws
					}
				else
					{
					help <- 1/(1+lambd*tildeZ[[k]]$s)
					mutmp <- help*crossprod(OldK[[k]],(response - eta))
					draws <- rnorm(dP[k],mutmp,(help*sigma2)^0.5)
					sums <- sum(draws^2*tildeZ[[k]]$s)
					}

				hyperbTau2new <- hyperbtau + sums/2
				tau2[k] <- 1/(rgamma2(1,hyperaTau2new[k],hyperbTau2new)) # 1/(rgamma(1,shape=hyperaTau2new[k],rate=hyperbTau2new))

				# Evaluate center effect
				ce <- ce - sce[k]
				if(cmcheck[[k]]$a)
					sce[k] <- cK[[k]]%*%(cmcheck[[k]]$b%*%draws)
				else
					sce[k] <- cK[[k]]%*%draws
				ce <- ce + sce[k]

#id <- 1:length(draws)
#drawstmp <- draws
#draws <- rtr2(help,draws)
#matplot(x=id,y=cbind(drawstmp,draws,c(drawstmp-mean(drawstmp))),type="l",lty=1,col=c(1,2,3))
#abline(h=sum(draws),col=2,lty=2)
#abline(h=sum(c(drawstmp-mean(drawstmp))),col=3,lty=2)
#Sys.sleep(0.2)
				ftmp <- crossprod(ttildeZ[[k]],draws)

				# Centering
				if(cme[k] == 1)
					{
					# draws <- rtr(Sigma2,draws)
					ftmp <- ftmp - mean(ftmp)
					# ftmp <- rtr2(rep(1,length(ftmp)),ftmp)
					# ftmp <- mycenter(ftmp,Nnew[k])
					}
				fps[,k] <- ftmp[new[[k]]]
				eta <- eta + fps[,k]	

				# Storing
				if(i%in%thinsteps)
					{
					lambda[k,itrthin] <- tau2[k] #lambd
					drawTildeGamma[[k]][,itrthin] <- draws
					mDTG[[k]][,itrthin] <- mutmp
					if(ALL > 1)
						csK[k,itrthin] <- ce - sce[k]
					else
						csK[k,itrthin] <- ce
					}
				}

		# Sampling for multilevel effects
		if(M > 0)
			for(k in 1:M)
				{
				eta <- eta - fM[,k]
				r <- response - eta

				# Smooth multilevel effects
				if(spec[k] == 2)
					{
					if(hypident[k])
						{
						sum2 <- 0
						lambd <- sigma2/tau2M[[k]][1]
						}
					for(j in 1:nl[k])
						{
						if(!hypident[k])
							{
							if(is.na(mlambdacheck[[k]][j]))
								lambd <- sigma2/tau2M[[k]][j] 
							else
								lambd <- mlambdacheck[[k]][j]
							}

						if(tildeZmu[[k]][[j]]$cs > 0)
							{
							help <- tildeZmu[[k]][[j]]$ZZ + lambd*tildeZmu[[k]][[j]]$s
							mutmp <- crossprod(tildeZmu[[k]][[j]]$tildeZ,r[ind[[k]][[j]]])

							for(jj in 1:20)
								{
								for(d in 1:dPM[[k]][j])
									{
									mcsmu <- 0
									if(d > 1)
										mcsmu <- sum(drawsmCmu[[k]][[j]][1:(d-1)]*help[d,1:(d-1)])
									if(d < dPM[[k]][j])
										mcsmu <- mcsmu + sum(drawsmCmu[[k]][[j]][(d+1):dPM[[k]][j]]*help[d,(d+1):dPM[[k]][j]])

									mcsmu <- (mutmp[d] - mcsmu)/help[d,d]

									if(tildeZmu[[k]][[j]]$cs==1)
										{
										if(d==1)
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=-Inf,high=drawsmCmu[[k]][[j]][2])
										if(d==dPM[[k]][j])
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmCmu[[k]][[j]][d-1],high=Inf)
										if(d > 1 && d < dPM[[k]][j])
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmCmu[[k]][[j]][d-1],high=drawsmCmu[[k]][[j]][d+1])
										}
									else
										{
										if(d==1)
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmCmu[[k]][[j]][2],high=Inf)
										if(d==dPM[[k]][j])
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=-Inf,high=drawsmCmu[[k]][[j]][d-1])
										if(d > 1 && d < dPM[[k]][j])
											drawsmCmu[[k]][[j]][d] <- rtvnorm2(1,mean=mcsmu,sd=sqrt(sigma2/help[d,d]),low=drawsmCmu[[k]][[j]][d+1],high=drawsmCmu[[k]][[j]][d-1])
										}
									}
								}
							draws <- drawsmCmu[[k]][[j]]
							sums <- crossprod(draws,tildeZmu[[k]][[j]]$s)%*%draws
							}
						else
							{
							help <- 1/(1+lambd*tildeZmu[[k]][[j]]$s)

							mutmp <- help*crossprod(tildeZmu[[k]][[j]]$tildeZ,r[ind[[k]][[j]]]) 
							draws <- rnorm(dPM[[k]][j],mutmp,(help*sigma2)^0.5)
							sums <- sum(draws^2*tildeZmu[[k]][[j]]$s)
							}

						# Draws for tau2M
						if(!hypident[k])
							{
							hyperbTau2newM <- hyperbtau + sums/2
							tau2M[[k]][j] <- 1/(rgamma2(1,hyperaTau2newM[[k]][j],hyperbTau2newM))
							}
						else
							sum2 <- sum2 + sums

						# Evaluate center effect
						ce <- ce - cM[[k]][j]
						if(cMM[[k]][[j]]$a)
							cM[[k]][j] <- cMM[[k]][[j]]$b%*%(cMM[[k]][[j]]$c%*%draws)
						else
							cM[[k]][j] <- cMM[[k]][[j]]$b%*%draws
						ce <- ce + cM[[k]][j]
					
						fpsM[[k]][[j]] <- crossprod(ttildeZmu[[k]][[j]],draws)
						# Centering
						if(mcentering[[k]][j])
							{
							#draws <- rtr(Sigma2M,draws)
							fpsM[[k]][[j]] <- fpsM[[k]][[j]] - mean(fpsM[[k]][[j]])
							}

						fM[ind[[k]][[j]],k] <- fpsM[[k]][[j]]

						# Storing
						if(i%in%thinsteps)
							{
							if(hypident[k])
								lambdaM[[k]][j,itrthin] <- tau2M[[k]][1]
							else
								lambdaM[[k]][j,itrthin] <- tau2M[[k]][j]
							drawTildeGammaM[[k]][[j]][,itrthin] <- draws
							drawTGMM[[k]][[j]][,itrthin] <- mutmp
							if(ALL > 1)
								csM[[k]][j,itrthin] <- ce - cM[[k]][j]
							else
								csM[[k]][j,itrthin] <- ce
							}
						}
					if(hypident[k])
						{
						hyperbTau2newM <- hyperbtau + sum2/2
						tau2M[[k]][1] <- 1/(rgamma2(1,hyperaTau2newM[[k]],hyperbTau2newM))
						}
					}

				# Linear multilevel effects
				if(spec[k] == 1)
					for(j in 1:nl[k])
						{
						mutmp <- crossprod(tildeZmu[[k]][[j]]$tildeZ,r[ind[[k]][[j]]]) 
						draws <- rnorm(dPM[[k]][j],mutmp,(sigma2^0.5)*IM[[k]][[j]])
						fM[ind[[k]][[j]],k] <- crossprod(ttildeZmu[[k]][[j]],draws)

						# Evaluate center effect
						ce <- ce - cM[[k]][j]
						tmp <- (cMM[[k]][[j]]$a%*%draws)*cMM[[k]][[j]]$d

						cM[[k]][j] <- cMM[[k]][[j]]$b%*%tmp
						ce <- ce + cM[[k]][j]

						# Storing
						if(i%in%thinsteps)
							{
							drawTildeGammaM[[k]][[j]][,itrthin] <- draws
							drawTGMM[[k]][[j]][,itrthin] <- mutmp
							for(jj in 1:cMM[[k]][[j]]$c)
								{
								if(ALL > 1)
									csM[[k]][[j]][jj,itrthin] <- ce - tmp[jj,]
								else
									csM[[k]][[j]][jj,itrthin] <- ce
								}
							}
						}
				eta <- eta + fM[,k]
				}

		# Sampling for random effects
		if(R > 0)
			for(k in 1:R)
				{
				eta <- eta - fR[,k]

				if(hspec[k] > 0)
					{
					# 1st stage
					if(rwcheck[k])
						{
						lambdR <- (sigma2/sigma2R[k])
						SIGMA <- 1/(ZZ[[k]] + lambdR*Z[[k]]$s)
						}
					else
						{
						lambdR <- (sigma2/sigma2R[k])*Z[[k]]$s
						SIGMA <- 1/(ZZ[[k]] + lambdR)
						}

					MU <- SIGMA*(crossprod(Z[[k]]$tildeZ,(response - eta)) + lambdR*etaR[[k]])
					drawers <- rnorm(BR[k],MU,(sigma2*SIGMA)^0.5)
					drawBetaR[[k]] <- Z[[k]]$RQ%*%drawers

					ce <- ce - cR[k]
					cR[k] <- cRM[[k]]%*%drawBetaR[[k]]
					ce <- ce + cR[k]

					# 3rd stage
					if(RRi[k] > 0)
						{
						for(kk in 1:ins[[k]]$R)
							{
							etaR[[k]] <- etaR[[k]] - ins[[k]]$fR[,kk]
					
							if(ins[[k]]$hspec[kk] > 0)
								{
								# Stage 1
								if(ins[[k]]$rwcheck[kk])
									{
									lambdRi <- sigma2R[k]/ins[[k]]$sigma2R[kk]
									SIGMAi <- 1/(ins[[k]]$ZZ[[kk]] + lambdRi*ins[[k]]$Z[[kk]]$s)
									}
								else
									{
									lambdRi <- sigma2R[k]/ins[[k]]$sigma2R[kk]*ins[[k]]$Z[[kk]]$s
									SIGMAi <- 1/(ins[[k]]$ZZ[[kk]] + lambdRi)
									}

								MUi <- SIGMAi*(crossprod(ins[[k]]$Z[[kk]]$tildeZ,(drawBetaR[[k]] - etaR[[k]])) + lambdRi*ins[[k]]$etaR[[kk]])		
								drawersi <- rnorm(ins[[k]]$BR[kk],MUi,(sigma2R[k]*SIGMAi)^0.5) 
								ins[[k]]$drawBetaR[[kk]] <- ins[[k]]$Z[[kk]]$RQ%*%drawersi

								ce <- ce - ins[[k]]$cR[kk]
								ins[[k]]$cR[kk] <- ins[[k]]$cRM[[kk]]%*%ins[[k]]$drawBetaR[[kk]]
								ce <- ce + ins[[k]]$cR[kk]

								if(ins[[k]]$RRs[kk] > 0)
									{
									for(jj in 1:ins[[k]]$RRs[kk])
										{
										ins[[k]]$etaR[[kk]] <- ins[[k]]$etaR[[kk]] - ins[[k]]$fr[[kk]][[jj]]

										if(ins[[k]]$la[[kk]][[jj]]$a)
											lambdRi <- ins[[k]]$sigma2R[kk]/ins[[k]]$tau2R[[kk]][jj]
										else
											lambdRi <- ins[[k]]$la[[kk]][[jj]]$b
	
										helpRi <- 1/(1+lambdRi*ins[[k]]$tildeZra[[kk]][[jj]]$s)

										mutmp <- helpRi*crossprod(ins[[k]]$OTzRa[[kk]][[jj]],ins[[k]]$drawBetaR[[kk]] - ins[[k]]$etaR[[kk]])
										draws <- rnorm(ins[[k]]$dPR[[kk]][jj],mutmp,(helpRi*ins[[k]]$sigma2R[kk])^0.5)

										ftmpri <- crossprod(ins[[k]]$ttildeZra[[kk]][[jj]],draws)
										if(ins[[k]]$center[[kk]][jj] > 0)
											{
											# ftmpri <- mycenter(ftmpri,ins[[k]]$newBR[[kk]][jj])
											ftmpri <- ftmpri - mean(ftmpri)
											}
											
										sums <- sum(draws^2*ins[[k]]$tildeZra[[kk]][[jj]]$s)
										hyperbTau2newRi <- hyperbtau + sums/2
										ins[[k]]$tau2R[[kk]][jj] <- 1/(rgamma2(1,ins[[k]]$hyperaTau2newR[[kk]][jj],hyperbTau2newRi))

										ce <- ce - ins[[k]]$cRS[[kk]][jj]
										if(ins[[k]]$cRMSS[[kk]][[jj]]$b)
											ins[[k]]$cRS[[kk]][jj] <- ins[[k]]$cRMSS[[kk]][[jj]]$a(ins[[k]]$cRMSS[[kk]][[jj]]$c%*%draws)
										else
											ins[[k]]$cRS[[kk]][jj] <- ins[[k]]$cRMSS[[kk]][[jj]]$a%*%draws
										ce <- ce + ins[[k]]$cRS[[kk]][jj]

										if(i%in%thinsteps)
											{
											ins[[k]]$lambdaR[[kk]][[jj]][itrthin] <- ins[[k]]$tau2R[[kk]][jj]
											ins[[k]]$drawTildeGammaR[[kk]][[jj]][,itrthin] <- draws
											ins[[k]]$drawTildeGammaRM[[kk]][[jj]][,itrthin] <- mutmp
											if(ALL > 1)
												ins[[k]]$cRSS[[kk]][[jj]][itrthin] <- ce - ins[[k]]$cRS[[kk]][jj]
											else
												ins[[k]]$cRSS[[kk]][[jj]][itrthin] <- ce
											}
										ins[[k]]$fr[[kk]][[jj]] <- ftmpri[ins[[k]]$newR[[kk]][[jj]]]
										ins[[k]]$etaR[[kk]] <- ins[[k]]$etaR[[kk]] + ins[[k]]$fr[[kk]][[jj]]
										}
									}

								# Linear
								if(ins[[k]]$RRl[kk] > 0)
									{
									ins[[k]]$etaR[[kk]] <- ins[[k]]$etaR[[kk]] - ins[[k]]$frl[[kk]]
									mutmp <- crossprod(ins[[k]]$tZRl[[kk]],ins[[k]]$drawBetaR[[kk]] - ins[[k]]$etaR[[kk]])
									draws <- rnorm(ins[[k]]$dPRl[kk],mutmp,(ins[[k]]$sigma2R[kk]*ins[[k]]$IXRl[[kk]])^0.5)
									ins[[k]]$frl[[kk]] <- crossprod(ins[[k]]$ttZRl[[kk]],draws)
									ins[[k]]$etaR[[kk]] <- ins[[k]]$etaR[[kk]] + ins[[k]]$frl[[kk]]

									ce <- ce - ins[[k]]$cRL[[kk]]
									tmp <- (ins[[k]]$cRML[[kk]]$a%*%draws)*ins[[k]]$cRML[[kk]]$b
									ins[[k]]$cRL[[kk]] <- ins[[k]]$cRML[[kk]]$d%*%tmp
									ce <- ce + ins[[k]]$cRL[[kk]]

									# Storing
									if(i%in%thinsteps)
										{
										ins[[k]]$drawTZRl[[kk]][,itrthin] <- draws
										ins[[k]]$drawTZRlM[[kk]][,itrthin] <- mutmp
										ins[[k]]$cRLL[[kk]][,itrthin] <- tmp
										}
									}

								# Draws for random sigma2
								hbRnew <- hyperbsigma + crossprod(ins[[k]]$drawBetaR[[kk]] - ins[[k]]$etaR[[kk]])/2
								ins[[k]]$sigma2R[kk] <- 1/(rgamma2(1,ins[[k]]$haRnew[kk],hbRnew))

								# Storing
								if(i%in%thinsteps)
									{		
									ins[[k]]$drawB[[kk]][,itrthin] <- drawersi # - ins[[k]]$etaR[[kk]] # ins[[k]]$drawBetaR[[kk]]
									ins[[k]]$drawBM[[kk]][,itrthin] <- MUi
									ins[[k]]$stetaR[[kk]][,itrthin] <- drawersi - ins[[k]]$etaR[[kk]]
									ins[[k]]$S2R[kk,itrthin] <- ins[[k]]$sigma2R[kk]
									if(ALL > 1)
										ins[[k]]$cRMS[kk,itrthin] <- ce - ins[[k]]$cR[kk]
									else
										ins[[k]]$cRMS[kk,itrthin] <- ce 
									}
								}
							else
								{
								SIGMAi <- 1/(ins[[k]]$ZZ[[kk]] + ins[[k]]$sigma2R[kk]/sigma2R[k]*ins[[k]]$Z[[kk]]$s)
								MUi <- SIGMAi*crossprod(ins[[k]]$Z[[kk]]$tildeZ,(drawBetaR[[k]] - etaR[[k]]))

								SIGMAi <- (sigma2R[k]*SIGMAi)^0.5
								drawersi <- rnorm(ins[[k]]$BR[kk],MUi,SIGMAi)
								ins[[k]]$drawBetaR[[kk]] <- ins[[k]]$Z[[kk]]$RQ%*%drawersi

								ce <- ce - ins[[k]]$cR[kk]
								ins[[k]]$cR[kk] <- ins[[k]]$cRM[[kk]]%*%ins[[k]]$drawBetaR[[kk]]
								ce <- ce + ins[[k]]$cR[kk]
				
								# Draws for random sigma2
								hbRnew <- hyperbsigma + crossprod((((drawersi^2)*ins[[k]]$Z[[kk]]$s)),ins[[k]]$ZZ[[kk]])/2
								ins[[k]]$sigma2R[kk] <- 1/(rgamma2(1,ins[[k]]$haRnew[kk],hbRnew))

								# Storing
								if(i%in%thinsteps)
									{
									ins[[k]]$drawB[[kk]][,itrthin] <- drawersi # ins[[k]]$drawBetaR[[kk]]
									ins[[k]]$drawBM[[kk]][,itrthin] <- MUi
									ins[[k]]$S2R[kk,itrthin] <- ins[[k]]$sigma2R[kk]
									if(ALL > 1)
										ins[[k]]$cRMS[kk,itrthin] <- ce - ins[[k]]$cR[kk]	
									else
										ins[[k]]$cRMS[kk,itrthin] <- ce
									}
								}

							if(ins[[k]]$bcheck[kk] > 0)
								ins[[k]]$fR[,kk] <- ins[[k]]$Z[[kk]]$tildeZ%*%drawersi
							else
								ins[[k]]$fR[,kk] <- ins[[k]]$drawBetaR[[kk]][ins[[k]]$rind[[kk]]]

							if(ins[[k]]$rcc[kk] > 0)
								ins[[k]]$fR[,kk] <- ins[[k]]$fR[,kk] - mean(ins[[k]]$fR[,kk])
							etaR[[k]] <- etaR[[k]] + ins[[k]]$fR[,kk]
							metaR[[k]] <- metaR[[k]] + etaR[[k]]
							}
						 }
					# 2nd stage
					# Smooth
					if(RRs[k] > 0)
						{
						for(j in 1:RRs[k])
							{
							etaR[[k]] <- etaR[[k]] - fr[[k]][[j]]

							if(la[[k]][[j]]$a)
								lambdR <- sigma2R[k]/tau2R[[k]][j]
							else
								lambdR <- la[[k]][[j]]$b

							helpR <- 1/(1+lambdR*tildeZra[[k]][[j]]$s)
							mutmp <- helpR*crossprod(OTzRa[[k]][[j]],(drawBetaR[[k]] - etaR[[k]]))
							draws <- rnorm(dPR[[k]][j],mutmp,(helpR*sigma2R[k])^0.5)

							# draws <- rtr(sigma2GammaR,draws)
              						ftmpr <- crossprod(ttildeZra[[k]][[j]],draws)

							if(center[[k]][j] > 0)

								{
								ftmpr <- ftmpr - mean(ftmpr)
								# ftmpr <- mycenter(ftmpr,newBR[[k]][j])
								}
							# Draws for random tau2
							# sum <- crossprod(((draws^2)*tildeZra[[k]][[j]]$s),IR[[k]][[j]])
							sums <- sum(draws^2*tildeZra[[k]][[j]]$s)
							hyperbTau2newR <- hyperbtau + sums/2
							tau2R[[k]][j] <- 1/(rgamma2(1,hyperaTau2newR[[k]][j],hyperbTau2newR)) 

							ce <- ce - cRS[[k]][j]
							if(cRMSS[[k]][[j]]$b)
								cRS[[k]][j] <- cRMSS[[k]][[j]]$a%*%(cRMSS[[k]][[j]]$c%*%draws)
							else
								cRS[[k]][j] <- cRMSS[[k]][[j]]$a%*%draws
							ce <- ce + cRS[[k]][j]

							# Storing
							if(i%in%thinsteps)
								{
								lambdaR[[k]][[j]][itrthin] <- tau2R[[k]][j]
								drawTildeGammaR[[k]][[j]][,itrthin] <- draws
								drawTildeGammaRM[[k]][[j]][,itrthin] <- mutmp
								if(ALL > 1)
									cRSS[[k]][[j]][itrthin] <- ce - cRS[[k]][j]
								else
									cRSS[[k]][[j]][itrthin] <- ce
								}
							fr[[k]][[j]] <- ftmpr[newR[[k]][[j]]]
							etaR[[k]] <- etaR[[k]] + fr[[k]][[j]]
							metaR[[k]] <- metaR[[k]] + etaR[[k]]
							}	
						}
					# Linear
					if(RRl[k] > 0)
						{
						etaR[[k]] <- etaR[[k]] - frl[[k]]
						mutmp <- crossprod(tZRl[[k]],(drawBetaR[[k]] - etaR[[k]]))
						draws <- rnorm(dPRl[k],mutmp,(sigma2R[k]^0.5*IXRl[[k]]))
						frl[[k]] <- crossprod(ttZRl[[k]],draws)
						etaR[[k]] <- etaR[[k]] + frl[[k]]
						metaR[[k]] <- metaR[[k]] + etaR[[k]]

						ce <- ce - cRL[[k]]
						tmp <- (cRML[[k]]$a%*%draws)*cRML[[k]]$b
						cRL[[k]] <- cRML[[k]]$d%*%tmp
						ce <- ce + cRL[[k]]

						# Storing
						if(i%in%thinsteps)
							{
							drawTZRl[[k]][,itrthin] <- draws
							drawTZRlM[[k]][,itrthin] <- mutmp
							if(ALL > 1)
								cRLL[[k]][,itrthin] <- ce - tmp
							else
								cRLL[[k]][,itrthin] <- ce
							}
						}

					# Draws for random sigma2
					hbRnew <- hyperbsigma + crossprod(drawBetaR[[k]] - etaR[[k]])/2
					sigma2R[k] <- 1/(rgamma2(1,haRnew[k],hbRnew))

					# Storing
					if(i%in%thinsteps)
						{		
						drawB[[k]][,itrthin] <- drawers # - etaR[[k]] # drawBetaR[[k]]
						drawBM[[k]][,itrthin] <- MU
						stetaR[[k]][,itrthin] <- drawers - etaR[[k]]
						S2R[k,itrthin] <- sigma2R[k]
						if(ALL > 1)
							cRMS[k,itrthin] <- ce - cR[k]
						else
							cRMS[k,itrthin] <- ce
						}
					}
				else
					{
					SIGMA <- 1/(ZZ[[k]] + (sigma2/sigma2R[k])*Z[[k]]$s)
					MU <- SIGMA*crossprod(Z[[k]]$tildeZ,(response - eta))

					drawers <- rnorm(BR[k],MU,(sigma2*SIGMA)^0.5)
					drawBetaR[[k]] <- Z[[k]]$RQ%*%drawers

					ce <- ce - cR[k]
					cR[k] <- cRM[[k]]%*%drawBetaR[[k]]
					ce <- ce + cR[k]
				
					# Draws for random sigma2
					sums <- sum((drawers^2)*Z[[k]]$s)
					hbRnew <- hyperbsigma + sums/2
					sigma2R[k] <- 1/(rgamma2(1,haRnew[k],hbRnew))

					# Storing
					if(i%in%thinsteps)
						{
						drawB[[k]][,itrthin] <- drawers # drawBetaR[[k]]
						drawBM[[k]][,itrthin] <- MU
						S2R[k,itrthin] <- sigma2R[k]
						if(ALL > 1)
							cRMS[k,itrthin] <- ce - cR[k]
						else
							cRMS[k,itrthin] <- ce
						}
					}

				# Get new predictor
				if(bcheck[k] > 0)
					fR[,k] <- Z[[k]]$tildeZ%*%drawers
				else
					fR[,k] <- drawBetaR[[k]][rind[[k]]]

				if(rcc[k] > 0)
					fR[,k] <- fR[,k] - mean(fR[,k])
				eta <- eta + fR[,k]	
				}
			
		# Sampling for linear predictor
		if(itcpt > 0 || L > 0)
			{
			eta <- eta - flin 
			mutmp <- crossprod(tildeX$tildeZ,(response - eta))
			draws <- rnorm(nclin,mutmp,(sigma*Ixx))

			# Evaluate center effects
			ce <- ce - bs
			bb <- crossprod(tbrq,draws)*bm
			bs <- sum(bb)
			ce <- ce + bs

			flin <- crossprod(ttildeX,draws)
			eta <- eta + flin
			}

		# Draws for sigma2
		if(!prob)
			{
			hyperbSigma2new <- hyperbsigma + crossprod(response - eta)/2
			sigma2 <- 1/(rgamma2(1,hyperaSigma2new,hyperbSigma2new))
			}

		# Storing and iterthin update
		if(i%in%thinsteps)
			{
			# storeETA[,itrthin] <- eta
			meta <- meta + eta
			diceta <- diceta + dev(response,eta,sigma2)
			S2[itrthin] <- sigma2
			if(itcpt > 0 || L > 0)
				{
				drawTildeBeta[,itrthin] <- draws
				drawTBML[,itrthin] <- mutmp
				
				if(L > 0)
					for(k in 2:LL)
						cL[(k-1),itrthin] <- ce - bb[k,1]
				}
			CE[itrthin] <- ce
			itrthin <- itrthin + 1
			}

		# Iterations
		if(trace && i == steps)
			{
			if(dots)
				cat(".")
			else
				cat("MCMC iteration:",i,"\n")
			steps <- steps + 1000
			}
		if(i == 100)
			{
			time <- proc.time() - ptm
			time <- as.vector(time[3]*iter/100)
			cat("Approximate sampling time is:",time,"sec\n")
			}
		}
	now <- proc.time()
	cat("MCMC iteration:",iter,"\n")
	}
	
	if(dots)
		cat("\n")
	samptime <- now - ptm
	res$samptime <- samptime
	cat("Total sampling time was:",samptime[3],"sec\n")
	#######################################################################################################################
        # End Sampler                                                                                                         #
	#######################################################################################################################

	# Only if data augmentation
	response <- origresponse

	# Generate output
	cat("Writing output...\n")

	meta <- meta/(itrthin-1)
	diceta <- diceta/(itrthin-1)
	ce <- mean(CE)

	# Extract fitted smooth functions and coefficients
	res$fout <- list()

	# Get predictor
	eta <- c(meta)
	if(K > 0)
		{
		tmpind <- rep(list(1:N),K)
		tmp <- extract(drawTildeGamma,lambda,
                               response,eta,tildeZ,csK,K,
                               TRUE,new,mDTG,S2,sm.specs,tmpind,data)
		attr(res$fout,"smooth.mat") <- tmp[[2]]
		res$fout[smooth] <- tmp[[1]]
		}
	if(M > 0)
		{
		tmp <- extract.mult(M,spec,nlm,nsm,nl,csM,FM,lambdaM,
                                    tildeZmu,response,eta,drawTildeGammaM,
                                    fitM,ind,dPM,S2,mult.specs,drawTGMM,data,itcpt)
		res$fout[mul] <- tmp[[1]]
		attr(res$fout,"lin.mu.mat") <- tmp[[2]]
		attr(res$fout,"smooth.mu.mat") <- tmp[[3]]
		}
	if(R > 0)
		{
        	res$fout[ran] <- extract.random(R, RR, RRl, RRs, RRi, fitR, RA, 
            			                lambdaR, tildeZra, drawTildeGammaR, drawTildeGammaRM,
            			                dPR, response, eta, drawB, drawBM, S2R, ins, drawTZRl, drawTZRlM, 
					        desRl, design, metaR, newR, FALSE, NULL, 
					        cRMS,cRSS,cRLL,S2,Z,1,data,stetaR)
		attr(res$fout,"ranpos") <- ran
		}
	# Extract linear effects
	linear.fout <- extract.linear(L,nclin,drawTildeBeta,drawTBML,tildeX,LIN,response,eta,cL,tnames,1,TRUE,itcpt)
	if(Lorig > 0)
		for(k in 1:Lorig)
			{
			if(attr(linposition[[k]],"factorcheck"))
				{
				tmp <- list(linear.fout$fits[linposition[[k]]])
				attr(tmp,"factorcheck") <- "factor"
				class(tmp) <- "linear.gibbs"
				}
			else
				{
				tmp <- linear.fout$fits[[linposition[[k]]]]
				attr(tmp,"factorcheck") <- "nonfactor"
				}
			attr(tmp,"term.type") <- "linear"
			res$fout[[lin[k]]] <- tmp
			}
	attr(res$fout,"lin.mat") <- linear.fout$beta
	res$coefficients <- getcoefs(res$fout,linear.fout$beta)
	res$residuals <- as.vector(response - eta)
	res$response <- response
	s2pmed <- median(S2)
	s2pmean <- mean(S2)
	s2pqu2p5 <- quantile(S2, probs = 0.025)
	s2pqu97p5 <- quantile(S2, probs = 0.975)
	s2pqu10 <- quantile(S2, probs = 0.1)
	s2pqu90 <- quantile(S2, probs = 0.9)
	s2psd <- sd(S2)
	res$fitted <- eta
	s2f <- cbind(s2pmean,s2psd,s2pqu2p5,s2pqu10,s2pmed,s2pqu90,s2pqu97p5)
	colnames(s2f) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	rownames(s2f) <- "Sigma2"
	attr(res$fitted,"variance") <- s2f
	attr(res$fitted,"variance.draws") <- S2
	Dm <- dev(response,meta,s2pmean)
	pd <- diceta - Dm
	res$DIC <- list(DIC = (diceta + pd), pd = pd)

	res$iterations <- iter
	res$burnin <- burnin
	res$thinning <- thinning
	
	class(res) <- "gibbs"
	return(res)
	}

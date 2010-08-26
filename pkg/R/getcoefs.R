getcoefs <- function(fout,beta)
	{
	coefs <- beta
	coefdim <- dim(coefs)
	coefdim[2] <- coefdim[2] - 1
	coefn <- rownames(coefs)
	coefs <- c(coefs[,1],coefs[,3:7])
	dim(coefs) <- coefdim
	rownames(coefs) <- coefn
	ok <- FALSE
	
	if(!is.null(attr(fout,"lin.mu.mat")))
		{
		ok <- TRUE
		tmp <- attr(fout,"lin.mu.mat")
		tmpdim <- dim(tmp)
		tmpdim[2] <- tmpdim[2] - 1
		tmpn <- rownames(tmp)
		tmp <- c(tmp[,1],tmp[,3:7])
		dim(tmp) <- tmpdim
		rownames(tmp) <- tmpn
		coefs <- rbind(coefs,tmp)
		}
	if(length(fout) > 0)
		{
		ok <- TRUE
		for(j in 1:length(fout))
			{
			if(attr(fout[[j]],"term.type") == "smooth")
				coefs <- rbind(coefs,attr(fout[[j]],"smooth.coef"))
			if(attr(fout[[j]],"term.type") == "m")
				{
				lx <- length(fout[[j]])
				if(!inherits(fout[[j]][[1]],"linear.gibbs"))
					for(i in 1:lx)
						coefs <- rbind(coefs,attr(fout[[j]][[i]],"smooth.coef"))
				}
			if(attr(fout[[j]],"term.type") == "random")
				{
				eff <- fout[[j]]$effects
				effn <- rownames(eff)
				deff <- dim(eff)
				eff <- c(eff[,2],eff[,4:8])
				deff[2] <- deff[2] - 4
				dim(eff) <- deff
				rownames(eff) <- effn
				coefs <- rbind(coefs,eff)
				if(!is.null(fout[[j]]$terms))
					{
					if(!is.null(attr(fout[[j]],"lin.mat")))
						{
						tmp <- attr(fout[[j]],"lin.mat")
						tmpdim <- dim(tmp)
						tmpdim[2] <- tmpdim[2] - 1
						tmpn <- rownames(tmp)
						tmp <- c(tmp[,1],tmp[,3:7])
						dim(tmp) <- tmpdim
						rownames(tmp) <- tmpn
						coefs <- rbind(coefs,tmp)
						}
					for(i in 1:length(fout[[j]]$terms))
						{
						if(attr(fout[[j]]$terms[[i]],"term.type") == "smooth")
							coefs <- rbind(coefs,attr(fout[[j]]$terms[[i]],"smooth.coef"))
						if(attr(fout[[j]]$terms[[i]],"term.type") == "random")
							{
							eff <- fout[[j]]$terms[[i]]$effects
							effn <- rownames(eff)
							deff <- dim(eff)
							eff <- c(eff[,2],eff[,4:8])
							deff[2] <- deff[2] - 4
							dim(eff) <- deff
							rownames(eff) <- effn
							coefs <- rbind(coefs,eff)
							if(!is.null(fout[[j]]$terms[[i]]$terms))
								{
								if(!is.null(attr(fout[[j]]$terms[[i]],"lin.mat")))
									{
									tmp <- attr(fout[[j]]$terms[[i]],"lin.mat")
									tmpdim <- dim(tmp)
									tmpdim[2] <- tmpdim[2] - 1
									tmpn <- rownames(tmp)
									tmp <- c(tmp[,1],tmp[,3:7])
									dim(tmp) <- tmpdim
									rownames(tmp) <- tmpn
									coefs <- rbind(coefs,tmp)
									}
								for(jj in 1:length(fout[[j]]$terms[[i]]$terms))
									{
									if(attr(fout[[j]]$terms[[i]]$terms[[jj]],"term.type") == "smooth")
										coefs <- rbind(coefs,attr(fout[[j]]$terms[[i]]$terms[[jj]],"smooth.coef"))
									}
								}
							}
						}
					}
				}
			}
		}
	if(ok)
		{
		if(all(coefs[1,]==0))
			coefs <- coefs[2:nrow(coefs),]
		}
	colnames(coefs) <- c("pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	return(coefs)
	}

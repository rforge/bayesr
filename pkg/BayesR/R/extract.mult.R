extract.mult <- function(M,spec,nlm,nsm,nl,csM,FM,lambda,Z,response,eta,
                         drawsc,fitM,ind,dPM,S2,mult.specs,mcdraws,data,itcpt)
	{
	fout <- vector("list",length=M)
	sm.par.mat <- lin.par.mat <- NULL

	for(k in 1:M)
		{
		if(spec[k] == 2)
			{
			new <- vector("list", nl[k])
			for(j in 1:nl[k])
				new[[k]][[j]] <- 1:nrow(Z[[k]][[j]]$tildeZ)
			tmp <- extract(drawsc[[k]],lambda[[k]],response,eta,Z[[k]],
                                       csM[[k]],nl[k],TRUE,new[[k]],mcdraws[[k]],S2,mult.specs[[k]],ind[[k]],attr(mult.specs[[k]],"data"))
			fout[[k]] <- tmp[[1]]
			par.mat <- tmp[[2]]
			sm.par.mat <- rbind(sm.par.mat,par.mat)
      			}
		if(spec[k] == 1)
			{
			par.mat <- matrix(0,nl[k],7)
			for(j in 1:nl[k])
				{
				if(is.null(ncol(mult.specs[[k]][[j]]$x)))
					xtmp <- list(mult.specs[[k]][[j]]$x)
				else
					{
					xtmp <- vector("list",length=ncol(mult.specs[[k]][[j]]$x))
					for(jj in 1:ncol(mult.specs[[k]][[j]]$x))
						xtmp[[jj]] <- mult.specs[[k]][[j]]$x[,j]	
					}
				tmp <- extract.linear(dPM[[k]][j],dPM[[k]][j],drawsc[[k]][[j]],
                                              	      mcdraws[[k]][[j]],Z[[k]][[j]],xtmp,
                                              	      response[ind[[k]][[j]]],eta[ind[[k]][[j]]],csM[[k]][[j]],
                                              	      mult.specs[[k]][[j]]$tnames,0,FALSE,1)
				fout[[k]][[j]] <- tmp$fits[[1]]
				par.mat[j,] <- tmp$beta[1,]
				}
			rownames(par.mat) <- attr(mult.specs[[k]],"tnames")
			colnames(par.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
			lin.par.mat <- rbind(lin.par.mat,par.mat)
			}
		attr(fout[[k]],"term.type") <- "m"
		attr(fout[[k]],"par.mat") <- par.mat
		attr(fout[[k]],"m.term") <- attr(mult.specs[[k]],"m.term")
		class(fout[[k]]) <- "gibbs"
		}

	return(list(fout,lin.par.mat,sm.par.mat))
	}

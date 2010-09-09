plot.gibbs <- function(x, which = 1, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
		       colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, grid = 30, theta = 40, phi = 40, image = FALSE,
		       map = NULL, names = FALSE, values = NULL, range = NULL, pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, dgts = 2, lpos = c(0.2,0),
                       cex = 1, ask = FALSE, byplots = FALSE, p3d = FALSE, pcat = NULL, border = NULL, ...)
	{
	if(!inherits(x,"gibbs"))
		stop("Argument x is not a gibbs object!")
	if(diagnostics == 4)
		{
		bdraws <- attr(attr(x$fout,"lin.mat"),"linear.coef.draws")[1,]
		coefplot(bdraws,xlab,ylab,main,"(Intercept)",acf,FALSE,...)
		}
	if(diagnostics == 1)
		{
		yes <- TRUE
		residuals <- residuals(x)
		eta <- fitted(x)

           	dens <- density(residuals)
		ymaxdens <- max(dens$y)

		hst <- hist(residuals, plot=FALSE)
		yhistmax <- max(hst$density)

		if(yhistmax < ymaxdens)
			ylimits <- c(0,ymaxdens)
		else
			ylimits <- c(0,yhistmax)		

		xlab <- "residuals"
		ylab <- "density"

		if(is.null(main))
			par(mfrow=c(1,2))
		if(!is.null(main))
			par(mfrow=c(1,2),oma=c(0,0,3,0))

		hist(residuals,freq=FALSE,xlab=xlab,ylab=ylab,ylim=ylimits,main=main[1],cex.main=1,font.main=1,...)
		lines(dens)
		box()
		ylab <- xlab
		xlab <- "fitted values"
		plot(eta,residuals,xlab=xlab,ylab=ylab,main=main[2],cex.main=1,font.main=1,cex=cex,...)
		abline(h=0)
		if(!is.null(main) && (length(main)>2))
			mtext(main[3],side=3,line=0,outer=TRUE,font=2,cex=1.5)
		}
	if(diagnostics == 2 || diagnostics == 3 || diagnostics == FALSE)
		{
		if(is.null(attr(x,"term.type")))
			{
			if(length(x$terms) < which[1] || which[1] == 0)
				stop("Argument which is specified wrong, nothing to plot!")
			X <- x$fout[[which[1]]]
			}
		else 
			{
			X <- x
			if(attr(x,"term.type") == "random")
				which <- 1
			}
		type <- attr(X,"term.type")
		if(type == "smooth")
			{
			if(inherits(X,"sm.gibbs"))
				{
				plot.sm.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			              	      colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				              cex,ask,border,nrc,pal,...)
				}
			if(inherits(X,"mrf.gibbs"))
				{
				plot.mrf.gibbs(X,resid,map,names,values,colored,col,lwdc,lwdconf,range,pal,legend, 
			               	       scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos, 
				               const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
				}
			}
		if(type == "m")
			{
			lx <- length(X)
			if(length(which)>1)
				{
				start <- which[2]
				end <- which[2]
				lx <- 1
				}
			else
				{
				start <- 1
				end <- lx
				}
			if(inherits(X[[1]],"mrf.gibbs"))
				{
				if(!is.null(map))
					{
					if(ask)
						par(ask=TRUE)
					else
						setmfrow(lx)
					}
				}
			else
				{
				if(lx > 1)
					{
					if(ask)
						par(ask=TRUE)
					else
						setmfrow(lx)
					}
				}
			if(!is.null(xlab))
				if(length(xlab) != lx)
					xlab <- rep(xlab,lx)
			if(!is.null(ylab))
				if(length(ylab) != lx)
					ylab <- rep(ylab,lx)
			if(!is.null(zlab))
				if(length(zlab) != lx)
					zlab <- rep(zlab,lx)
			if(inherits(X[[1]],"sm.gibbs"))
				{
				for(j in start:end)
					{
					plot.sm.gibbs(X[[j]],resid,jit,const,diagnostics,acf,xlab[j],ylab[j],zlab[j],main[j], 
			              	              colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				                      cex,ask,border,nrc,pal,...)
					}
				}
			if(inherits(X[[1]],"mrf.gibbs"))
				{
				for(j in start:end)
					{
				        plot.mrf.gibbs(X[[j]],resid,map,names,values,colored,col,lwdc,lwdconf,range,pal,legend, 
			               	               scale,nrc,xlab[j],ylab[j],zlab[j],main[j],dgts,cex,lpos, 
				                       const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
					}
				}
			if(inherits(X[[1]],"linear.gibbs"))
				{
				for(j in start:end)
					{
					plot.linear.gibbs(X[[j]],resid,jit,const,diagnostics,acf,xlab[j],ylab[j], 
                             	 			  main[j],colored,col,lwdc,lwdconf,range,
                              			          cex,ask,...)
					}
				}
			}
		if(type == "random")
			{
			effplot <- FALSE
			if(!is.list(X))
				X <- list(effects=X,terms=NULL)
			if(!is.null(X$terms))
				{
				if(length(which) > 1)
					{
					if(length(which) == 2)
						{
						if(which[2] == 0)
							which[2] <- 1
						X <- X$terms
						if(which[2] > length(X))
							stop("Argument which is specified wrong, nothing to plot!")
						X <- X[[which[2]]]
						if(attr(X,"term.type") == "smooth")
							{
							if(inherits(X,"sm.gibbs"))
								{
								plot.sm.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			              	      				      colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				              				      cex,ask,border,nrc,pal,...)
								}
							if(inherits(X,"mrf.gibbs"))
								{
								plot.mrf.gibbs(X,resid,map,names,values,colored,col,lwdc,lwdconf,range,
									       pal,legend,scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos, 
				               				       const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
								}
							}
						if(attr(X,"term.type") == "linear")
							{
							plot.linear.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab, 
                              					          main,colored,col,lwdc,lwdconf,range,
                              					          cex,ask,...)
							}	
						if(attr(X,"term.type") == "random" && !byplots)
							{
							X <- X$effects
							plot.random.gibbs(X,resid,const,diagnostics,acf,xlab,ylab,zlab, 
                              					          main,colored,col,lwdc,lwdconf,map,names, 
                              					          range,pal,legend,scale,nrc,dgts,lpos,
                              					          cex,ask,p3d,theta,phi,pcat,border,...)
							}
						if(attr(X,"term.type") == "random" && byplots)
							{
							byyes1 <- byyes2 <- FALSE
							if(!is.null(attr(X,"byplots")))
								byyes2 <- TRUE
							if(!byyes2)
								if(!is.null(attr(X$effects,"byplots")))
									byyes1 <- TRUE
							if(byyes1 || byyes2)
								{
								if(byyes1)
									byp <- attr(X$effects,"byplots")
								if(byyes2)
									byp <- attr(X,"byplots")
								lbyp <- length(byp)*length(byp[[1]])
								mfrowpr <- as.integer(par()$mfrow)
								if(mfrowpr[1] == 1 && mfrowpr[2] == 1)
									setmfrow(lbyp)
								# needs to be checked!
								#if(length(which)<2)
								#	{
									start <- 1
									end <- length(byp)		
								#	}
								#else
								#	{
								#	if(length(byp)<which[2])
								#		start <- end <- 1
								#	else
								#		start <- end <- which[2]
								#	}
								for(j in start:end)
									{
									for(jj in 1:length(byp[[j]]))
										{
										if(attr(byp[[j]][[jj]],"term.type") == "smooth")
											{
											if(inherits(byp[[j]][[jj]],"sm.gibbs"))
												{
												plot.sm.gibbs(byp[[j]][[jj]],resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			              	      									colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				              									cex,ask,border,nrc,pal,...)
												}
											if(inherits(byp[[j]][[jj]],"mrf.gibbs"))
												{
												plot.mrf.gibbs(byp[[j]][[jj]],resid,map,names,values,colored,col,
										       				lwdc,lwdconf,range,pal,legend, 
			               	       					       				scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos, 
				               					       				const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
												}
											}
										if(attr(byp[[j]][[jj]],"term.type") == "linear")
											{
											plot.linear.gibbs(byp[[j]][[jj]],resid,jit,const,diagnostics,acf,xlab,ylab, 
                              						          			main,colored,col,lwdc,lwdconf,range,
                              						          			cex,ask,...)
											}
										}
									}
								effplot <- FALSE
								}		
							else
								effplot <- TRUE
							}
						}

					if(length(which) == 3)
						{
						if(which[2] == 0)
							stop("Argument which is specified wrong, nothing to plot!")
						if(which[3] == 0)
							which[3] <- 1
						X <- X$terms
						X <- X[[which[2]]]
						X <- X$terms
						if(which[3] > length(X))
							stop("Argument which is specified wrong, nothing to plot!")
						X <- X[[which[3]]]
						if(attr(X,"term.type") == "smooth")
							{
							if(inherits(X,"sm.gibbs"))
								{
								plot.sm.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			              	      				      colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				              				      cex,ask,border,nrc,pal,...)
								}
							if(inherits(X,"mrf.gibbs"))
								{
								plot.mrf.gibbs(X,resid,map,names,values,colored,col,lwdc,lwdconf,range,pal,legend, 
			               	       				       scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos, 
				               				       const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
								}
							}
						if(attr(X,"term.type") == "linear")
							{
							plot.linear.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab, 
                              						  main,colored,col,lwdc,lwdconf,range,
                              						  cex,ask,...)
							}
						if(attr(X,"term.type") == "random" && !byplots)
							{
							which <- which[3:length(which)]
							plot.gibbs(X,which,resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			     					   colored,col,lwdc,lwdconf,grid,theta,phi,image,map,names,values,range,pal,legend, 
								   scale,nrc,dgts,lpos,cex,ask,byplots,p3d,pcat,border, ...)
							}
						if(attr(X,"term.type") == "random" && byplots)
							{
							byyes1 <- byyes2 <- FALSE
							if(!is.null(attr(X,"byplots")))
								byyes2 <- TRUE
							if(!byyes2)
								if(!is.null(attr(X$effects,"byplots")))
									byyes1 <- TRUE
							if(byyes1 || byyes2)
								{
								if(byyes1)
									byp <- attr(X$effects,"byplots")
								if(byyes2)
									byp <- attr(X,"byplots")
								lbyp <- length(byp)*length(byp[[1]])
								mfrowpr <- as.integer(par()$mfrow)
								if(mfrowpr[1] == 1 && mfrowpr[2] == 1)
									setmfrow(lbyp)
								# needs to be checked!
								#if(length(which)<2)
								#	{
									start <- 1
									end <- length(byp)		
								#	}
								#else
								#	{
								#	if(length(byp)<which[2])
								#		start <- end <- 1
								#	else
								#		start <- end <- which[2]
								#	}
								for(j in start:end)
									{
									for(jj in 1:length(byp[[j]]))
										{
										if(attr(byp[[j]][[jj]],"term.type") == "smooth")
											{
											if(inherits(byp[[j]][[jj]],"sm.gibbs"))
												{
												plot.sm.gibbs(byp[[j]][[jj]],resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			              	      									colored,col,lwdc,lwdconf,grid,theta,phi,image, 
				              									cex,ask,border,nrc,pal,...)
												}
											if(inherits(byp[[j]][[jj]],"mrf.gibbs"))
												{
												plot.mrf.gibbs(byp[[j]][[jj]],resid,map,names,values,colored,col,
										       				lwdc,lwdconf,range,pal,legend, 
			               	       					       				scale,nrc,xlab,ylab,zlab,main,dgts,cex,lpos, 
				               					       				const,diagnostics,acf,p3d,theta,phi,pcat,border,...)
												}
											}
										if(attr(byp[[j]][[jj]],"term.type") == "linear")
											{
											plot.linear.gibbs(byp[[j]][[jj]],resid,jit,const,diagnostics,acf,xlab,ylab, 
                              						          			main,colored,col,lwdc,lwdconf,range,
                              						          			cex,ask,...)
											}
										}
									}
								effplot <- FALSE
								}		
							else
								effplot <- TRUE
							}
						}
					if(length(which) > 3)
						{
						if(which[4] == 0)
							which[4] <- 1
						X <- X$terms
						X <- X[[which[2]]]
						X <- X$terms
						X <- X[[which[3]]]
						X <- list(fout=X$terms,terms=rep("NA",length(X$terms)))
						class(X) <- "gibbs"
						which <- which[4:length(which)]
						plot.gibbs(X,which,resid,jit,const,diagnostics,acf,xlab,ylab,zlab,main, 
			     				   colored,col,lwdc,lwdconf,grid,theta,phi,image,map,names,values,range,pal,legend, 
							   scale,nrc,dgts,lpos,cex,ask,byplots,p3d,pcat,border, ...)
						return(invisible(NULL))
						}
					}
				else
					effplot <- TRUE
				}
			if(length(which)==1 && !byplots)
				effplot <- TRUE
			if(effplot)
				{
				X <- X$effects
				plot.random.gibbs(X,resid,const,diagnostics,acf,xlab,ylab,zlab, 
                              			  main,colored,col,lwdc,lwdconf,map,names, 
                              			  range,pal,legend,scale,nrc,dgts,lpos,
                              			  cex,ask,p3d,theta,phi,pcat,border,...)
				}
			}
		if(type == "linear")
			{
			plot.linear.gibbs(X,resid,jit,const,diagnostics,acf,xlab,ylab, 
                             	 	  main,colored,col,lwdc,lwdconf,range,cex,ask,...)
			}
		}
	return(invisible(NULL))
	}

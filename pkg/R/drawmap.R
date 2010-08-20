drawmap <- function(map, names = FALSE, values = NULL, colored = FALSE, range = NULL, pal = "rainbow",  
		    legend = TRUE, scale = 0.2, nrc = 100, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
		    dgts = 4, cex = 1, lpos = c(0.2,0), p3d = FALSE, theta = 0, phi = 25, pcat = NULL, border = NULL, 
		    density = NULL, resid = FALSE, ...)
	{
	if(missing(map))
		stop("Map object is missing!")
	if(!is.list(map))
		stop("Argument map must be a list() of polygons!")

	args <- as.list(substitute(list(...)))[-1]

	nm <- length(map)
	nam <- names(map)

	maxx <- minx <- maxy <- miny <- maxxe <- minxe <- maxye <- minye <- rep(0,nm)
	for(k in 1:nm)
		{
		maxxe[k] <- maxx[k] <- max(map[[k]][,1], na.rm=TRUE)
		minxe[k] <- minx[k] <- min(map[[k]][,1], na.rm=TRUE)
		maxye[k] <- maxy[k] <- max(map[[k]][,2], na.rm=TRUE)
		minye[k] <- miny[k] <- min(map[[k]][,2], na.rm=TRUE)	
		}
	maxx <- max(maxx)
	minx <- min(minx)
	maxy <- max(maxy)
	miny <- min(miny)

	if(is.null(values))
		legend <- FALSE

	if(legend)
		{
		ry <- abs(maxy - miny)
		rx <- abs(maxx - minx)
		maxy <- maxy + scale*ry
		miny <- miny - scale*ry
		maxx <- maxx + scale*rx
		minx <- minx - scale*rx
		}

	Y <- c(miny,maxy)
	X <- c(minx,maxx)
	
	foutcheck <- racheck <- pcatcheck <- FALSE
	if(!is.null(values))
		{
		if(inherits(values,"mrf.gibbs"))
			{
			foutcheck <- TRUE
			fout <- values
			if(!is.null(pcat))
				{
				if(pcat == 80)
					{
					values <- unique(cbind(values[,1],values[,10]))
					args$pcat80 <- NULL
					}
				if(pcat == 95)
					{
					values <- unique(cbind(values[,1],values[,9]))
					args$pcat95 <- NULL
					}
				pcatcheck <- TRUE
				nrc <- 3
				}
			else
				values <- unique(values[,1:2])
			}
		if(inherits(values,"random.gibbs"))
			{
			foutcheck <- TRUE
			fout <- values
			if(!is.null(pcat))
				{
				if(pcat == 80)
					{
					values <- unique(cbind(values[,1],values[,10]))
					args$pcat80 <- NULL
					}
				if(pcat == 95)
					{
					values <- unique(cbind(values[,1],values[,9]))
					args$pcat95 <- NULL
					}
				pcatcheck <- TRUE
				nrc <- 3
				}
			else
				values <- unique(values[,1:2])
			residr <- attr(fout,"random.partial.resid")
			racheck <- TRUE
			}
		identify <- as.character(values[,1])
		values <- values[,2]
		minv <- min(values)
		maxv <- max(values)
		maxvt <- round(maxv, digits=dgts)
		minvt <- round(minv, digits=dgts)

		#if(nrc > length(unique(values)))
		#	nrc <- length(unique(values))
		}
	rangevals <- NULL
	if(!is.null(range))
		{
		if(length(range) != 2)
			stop("Argument range must be a vector with a lower and upper limit!")
		if(range[1] > range[2])
			range <- c(range[2],range[1])

		values[values < range[1]] <- range[1]
		values[values > range[2]] <- range[2]

		minv <- min(c(values,range[1]))
		maxv <- max(c(values,range[2]))
		# minvt <- paste("<", round(minv, digits=dgts))
		# maxvt <- paste(">", round(maxv, digits=dgts))
		minvt <- paste(round(minv, digits=dgts))
		maxvt <- paste(round(maxv, digits=dgts))

		rangevals <- range
		}
	if(!is.null(values) && !colored)
		{
		if(any(values < 0))
			values <- (values + abs(min(values)))
		else
			values <- values - min(values)
		if(pcatcheck)
			valcut <- values
		else
			valcut <- cut(c(values,rangevals), nrc, include.lowest=TRUE)
		valcol <- as.vector(valcut, mode="numeric")
		mvc <- max(valcol) + 1
		valcol <- mvc - valcol
		# mvc <- mvc - 1

		colors <- gray(valcol/mvc)
		lgc <- gray(nrc:1/nrc)
		}
	if(!is.null(values) && colored)
		{
		if(length(which(pal == "rainbow")) != 0)
			colors <- rainbow(nrc)
		if(length(which(pal == "heat")) != 0)
			colors <- heat.colors(nrc)
		if(length(which(pal == "terrain")) != 0)
			colors <- terrain.colors(nrc)
		if(length(which(pal == "topo")) != 0)
			colors <- topo.colors(nrc)
		if(length(which(pal == "cm")) != 0)
			colors <- cm.colors(nrc)
		if(length(which(pal == "hcl")) != 0)
			{
			if(is.null(args$hcl.par))
				hcl.par <- list(h=c(130,25), c=100, l=c(90,70))
			else
				hcl.par <- args$hcl.par
			args$hcl.par <- NULL
                	colors <- colorspace::diverge_hcl(nrc,h=hcl.par$h,c=hcl.par$c,l=hcl.par$l)
			}
		if(length(which(pal == c("rainbow","heat","terrain","topo","cm","hcl"))) == 0)
			{
			warning("Argument pal specified wrong, set to default!")
			colors <- rainbow(nrc)
			}
		colors <- rev(colors)
		lgc <- rev(colors)
		if(pcatcheck)
			valcut <- c(values)
		else		
			valcut <- cut(c(values,rangevals), nrc)	
		valcol <- as.vector(valcut, mode="numeric") 
		mvc <- max(valcol)+1
		valcol <- mvc - valcol
		if(is.null(args$swap))
			swap <- FALSE
		else
			swap <- TRUE
		args$swap <- NULL
		if(swap)
			{
			valcol <- rev(valcol)
			lgc <- rev(lgc)
			}
		colors <- colors[valcol]
		}
	if(is.null(xlab))
		xlabel <- ""
	else
		xlabel <- xlab
	if(is.null(ylab))
		ylabel <- ""
	else
		ylabel <- ylab

	if(p3d)
		{
		if(foutcheck)
			{
			if(racheck)
				zmat <- matrix(seq(min(c(unlist(residr),fout[,2])),max(c(unlist(residr),fout[,2])),length=4),2,2)
			else
				zmat <- matrix(seq(min(c(fout[,8],fout[,2])),max(c(fout[,8],fout[,2])),length=4),2,2)
			}
		else
			zmat <- matrix(runif(4),2,2)
		mzmat <- mean(zmat)

		args$y <- substitute(Y)
		args$x <- substitute(X)
		args$z <- substitute(zmat)
		args$xlab <- xlab
		args$ylab <- ylab
		args$zlab <- zlab
		args$main <- main
		args$col <- NA
		args$border <- NA
		if(is.null(args$box))
			args$box <- FALSE		
		if(is.null(args$zlim))
			args$zlim <- c(min(zmat),max(zmat))
		if(is.null(args$axes))
			args$axes <- FALSE
		if(is.null(args$expand))
			{
			if(foutcheck)
				args$expand <- 1
			else
				args$expand <- 0
			}
		args$theta <- theta
		args$phi <- phi
		pmat <- do.call("persp",args)
		if(foutcheck)
			foutlevels <- as.integer(levels(as.factor(fout[,1])))
		if(foutcheck & resid)
			{
			cme <- 1
			idind <- 1:length(identify)
			for(k in 1:nm)
				{
				if(nam[k]%in%identify)
					{
					pos <- centroidpos(map[[k]])
					if(racheck)
						e <- residr[[idind[identify==nam[k]]]]
					else
						e <- fout[fout[,1]==foutlevels[cme],8]
					e <- e - (mean(e) - mzmat)
					vals <- values[cme] -(mean(e) - mzmat)
					check <- e<mzmat
					if(any(check==TRUE))
						e <- e[check]
					else
						e <- mzmat
					pos1 <- rep(pos[1],length(e))
					pos2 <- rep(pos[2],length(e))
					t3d <- trans3d(pos1,pos2,e,pmat)
					
					points(t3d$x,t3d$y,cex=cex)
					#if(vals < mzmat)
						#{
						#t3d <- trans3d(pos[1],pos[2],vals,pmat)
						#points(t3d$x,t3d$y,cex=cex,col=colors[cme],pch=16)
						#t3d <- trans3d(pos[1],pos[2],vals + 0.05*(max(zmat)-min(zmat)),pmat)
						#text(t3d$x,t3d$y,round(values[cme],dgts),cex=cex)
						#}
					cme <- cme + 1
					}
				}
			}
		cme <- 1
		for(k in 1:nm)
			{
			t3d <- trans3d(map[[k]][,1],map[[k]][,2],mzmat,pmat)
			if(is.null(values))
				polygon(x=t3d$x,y=t3d$y,lwd=0.3,border=border,density=density)
			else
				{
				if(nam[k]%in%identify)
					{
					polygon(x=t3d$x,y=t3d$y,col=colors[cme],border=border)
					cme <- cme + 1
					}
				else
					polygon(x=t3d$x,y=t3d$y,border=border,density=density)	
				}
			if(names)
				{
				pos <- centroidpos(map[[k]])
				if(is.null(nam))
					txt <- paste(k)
				else
					txt <- nam[k]
				t3d <- trans3d(pos[1],pos[2],mzmat,pmat)
				text(t3d$x,t3d$y,txt,cex=cex)
				}
			}
		if(foutcheck & resid)
			{
			cme <- 1
			for(k in 1:nm)
				{
				if(nam[k]%in%identify)
					{
					pos <- centroidpos(map[[k]])
					if(racheck)
						e <- residr[[idind[identify==nam[k]]]]
					else
						e <- fout[fout[,1]==foutlevels[cme],8]
					e <- e - (mean(e) - mzmat)
					vals <- values[cme] -(mean(e) - mzmat)
					check <- e>=mzmat
					if(any(check==TRUE))
						e <- e[check]
					else
						e <- mzmat
					pos1 <- rep(pos[1],length(e))
					pos2 <- rep(pos[2],length(e))
					t3d <- trans3d(pos1,pos2,e,pmat)
					points(t3d$x,t3d$y,cex=cex)
					#if(vals < mzmat)
						#{
						#t3d <- trans3d(pos[1],pos[2],vals,pmat)
						#points(t3d$x,t3d$y,cex=cex,col=colors[cme],pch=16)
						#t3d <- trans3d(pos[1],pos[2],vals + 0.05*(max(zmat)-min(zmat)),pmat)
						#text(t3d$x,t3d$y,round(values[cme],dgts),cex=cex)
						#}
					cme <- cme + 1
					}
				}
			}
		}
	else
		{
		args$x <- X
		args$y <- Y
		args$ylab <- ylabel
		args$xlab <- xlabel
		args$axes <- FALSE
		args$type <- "n"
		args$main <- main
		if(!is.null(args$add))
			{
			if(args$add)
				par(new = TRUE)
			args$add <- NULL
			}
		do.call("plot",args)
		# plot(Y~X, xlab=xlabel, ylab=ylabel, axes=FALSE, type="n", main=main,...)
		if(!is.null(values))
			{
			cme <- 1
			for(k in 1:nm)
				{
				if(nam[k]%in%identify)	
					{
					polygon(map[[k]][,1],map[[k]][,2],col=colors[cme],border=border)
					cme <- cme + 1
					}
				else
					polygon(map[[k]][,1],map[[k]][,2],border=border,density=density)

				if(names)
					centroidtext(map[[k]],nam,k,cex)
				}
			}
		}
	if(legend && !is.null(values))
		{
		ylo <- yro <- Y[1] + (0.7 * (Y[2] - Y[1]))/12 + lpos[2]*(Y[2] - Y[1])
               	ylu <- yru <- Y[1] + (0.3 * (Y[2] - Y[1]))/12 + lpos[2]*(Y[2] - Y[1])
            	tylu <- tyru <- Y[1] + lpos[2]*(Y[2] - Y[1])
            	xlu <- xlo <- X[1] + lpos[1] * (X[2] - X[1])
            	xru <- xro <- X[1] + 0.3 * (X[2] - X[1]) + lpos[1] * (X[2] - X[1])
            	step <- (xru - xlu)/nrc
		for (i in 0:(nrc-1)) 
			{
			if(p3d)
				{
				t3d <- trans3d(c(xlo + step * i, xlo + step * (i+1), 
                  		xlu + step * (i+1), xlu + step * i),c(ylo, yro, yru, ylu),mzmat,pmat)
				polygon(t3d$x, t3d$y, col=lgc[i+1], border=lgc[i+1])
				}
			else
				{
                		polygon(c(xlo + step * i, xlo + step * (i+1), 
                  		xlu + step * (i+1), xlu + step * i), c(ylo, 
                  		yro, yru, ylu), col=lgc[i+1], border=lgc[i+1])
				}
            	}
		if(p3d)
			{
			t3d <- trans3d(c(xlo, xro, xru, xlu, xlo),c(ylo, yro, yru, ylu, ylo),mzmat,pmat)
			lines(t3d$x, t3d$y, col = "black")
			t3d <- trans3d(xlu + 0.5 * step,tylu,mzmat,pmat)
			text(t3d$x, t3d$y, minvt, cex=cex)	
			t3d <- trans3d(xru - 0.5 * step,tyru,mzmat,pmat)		
            		text(t3d$x, t3d$y, maxvt, cex=cex)

			if(minv < 0 && maxv > 0) 
				{
				if(pcatcheck)
					help <- c(-1,0,1)
				else
					help <- cut(c(0, minv, maxv), nrc)
               		help <- as.vector(help, mode = "numeric")
				t3d <- trans3d(xlu + step * (help[1] - 0.5),tylu,mzmat,pmat)
               		text(t3d$x, t3d$y, "0",cex=cex)
				}
			else
				{
				if(pcatcheck)
					help <- c(-1,0,1)
				else
					help <- cut(c(minv, maxv), nrc)
               		help <- as.vector(help, mode = "numeric")
				t3d <- trans3d(xlu + (nrc*step)*0.5,tylu,mzmat,pmat)
               		text(t3d$x, t3d$y, paste(round(((maxv-minv)/2+minv), digits=dgts)),cex=cex)
				}

			}
		else
			{
			lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, ylu, ylo), col = "black")
			text(xlu + 0.5 * step, tylu, minvt, cex=cex)			
            		text(xru - 0.5 * step, tyru, maxvt, cex=cex)

			if(minv < 0 && maxv > 0) 
				{
				if(pcatcheck)
					help <- cut(c(0, -1, 1), nrc)
				else
					help <- cut(c(0, minv, maxv), nrc)
               		help <- as.vector(help, mode = "numeric")
               		text(xlu + step * (help[1] - 0.5), tylu, "0",cex=cex)
				}
			else
				{
				if(pcatcheck)
					help <- cut(c(0, -1, 1), nrc)
				else
					help <- cut(c(minv, maxv), nrc)
               			help <- as.vector(help, mode = "numeric")
               			text(xlu + (nrc*step)*0.5, tylu, paste(round(((maxv-minv)/2+minv), digits=dgts)),cex=cex)
				}
			}
		}
	else
		{
		if(!p3d)
			{
			for(k in 1:nm)
				{
				polygon(map[[k]][,1],map[[k]][,2],border=border,density=density)
				if(names)
					centroidtext(map[[k]],nam,k,cex)
				}
			}
		}
	return(invisible(list(X=X,Y=Y)))
	}


plot.bnd <- function(map, names = FALSE, values = NULL, colored = FALSE, range = NULL, pal = "rainbow",  
		     legend = TRUE, scale = 0.2, nrc = 100, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
		     dgts = 4, cex = 1, lpos = c(0.2,0), p3d = FALSE, theta = 0, phi = 25, pcat = NULL, border = NULL, 
		     density = NULL, resid = FALSE, ...)
	{
	drawmap(map,names,values,colored,range,pal,legend,scale,nrc,xlab,ylab,zlab,main,
		dgts,cex,lpos,p3d,theta,phi,pcat,border,density,resid,...)
	return(invisible(NULL))
	}

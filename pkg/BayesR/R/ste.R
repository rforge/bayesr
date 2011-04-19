ste <- function(x,z, by = NULL, degree = c(3,3), knots = c(20,20), orderpenalty = c(2,2), lambda = NULL, 
		   center = TRUE, dgts = TRUE, map = NULL, ind = NULL)
	{
	indmap <- NULL
	if(!is.null(map))
		{
		cents <- matrix(0,nrow=length(map),ncol=2)
		#binds <- NULL
		for(i in 1:length(map))
			{
			cents[i,] <- centroidpos(map[[i]])
			#binds <- rbind(binds,map[[i]])
			}
		ind <- as.factor(ind)
		Z<-diag(nlevels(ind))[ind,]
		x <- as.vector(Z%*%cents[,1])
		z <- as.vector(Z%*%cents[,2])
		#x[x==min(x)] <- min(binds[,1])
		#x[x==max(x)] <- max(binds[,1])
		#y[y==min(y)] <- min(binds[,2])
		#y[y==max(y)] <- max(binds[,2])
		indmap <- ind
		}

	if(missing(x) || missing(z))
		stop("Variable missing in ste()!")
	else
		{
		call <- match.call(expand.dots=FALSE)
		initcall <- all.names(call)
		callsm <- as.character(as.expression(sys.call(0)))
		callsm <- paste(strsplit(callsm," ","")[[1]],collapse="")
		names <- initcall[initcall!="ste"]
		if(!is.null(map))
			names <- c(paste(names[length(names)],".xco",sep=""),paste(names[length(names)],".yco",sep=""))
		if(any("$"%in%names))
			names <- c(paste(names[2],"$",names[3],sep=""),paste(names[5],"$",names[6],sep=""))
		if(dgts != FALSE || dgts != F)
			{
			m <- cbind(x,z)
			xn <- get.unique(m,digits=dgts)
			ind <- xn$ind
			iind <- xn$iind
			xu <- xn$xu[,1]
			zu <- xn$xu[,2]
			}
		else
			{
			xu <- x
			zu <- z
			ind <- iind <- 1:length(xu)
			}
		n <- length(x)
		m <- length(z)
		mx <- approxm(x)
		mz <- approxm(z)
		if(n != m)
			stop("Variable length differing in ste()!")
		else
			{
			if(length(degree)==1)
				degree <- c(degree,degree)
			if(length(knots)==1)
				knots <- c(knots,knots)
			if(length(orderpenalty)==1)
				orderpenalty <- c(orderpenalty,orderpenalty)
			if((orderpenalty[1] > 0) && (orderpenalty[2] > 0))
				{
				ps1 <- sps(x=x,degree=degree[1],knots=knots[1],orderpenalty=orderpenalty[1])
				ps2 <- sps(x=z,degree=degree[2],knots=knots[2],orderpenalty=orderpenalty[2])

				I1 <- diag(1, nrow(ps1$K), ncol(ps1$K))
				I2 <- diag(1, nrow(ps2$K), ncol(ps2$K))
    
     		B <- matrix(0, n, 0)
				bmx <- ps1$basis[mx,]
				bmz <- ps2$basis[mz,]
				mxz <- matrix(0,1,0)
        for (k in 1:ncol(ps1$basis)) 
					{
					B <- cbind(B, ps1$basis[,k]*ps2$basis)
					mxz <- c(mxz,bmx[k]*bmz)
					}
				K <- kronecker(I2,ps1$K)+kronecker(ps2$K,I1)
			
				if(!is.null(by))
					{
					if(is.factor(by))
						{
						if(length(levels(by)) < 3)
							by <- (by==as.integer(levels(by)[2]))*1
						else
							stop("Too many levels in ste() argument by, see function mu()!")
						}
					B <- B*by
					}

				return(list(basis=B,by=by,K=K,x=x,z=z,degree=degree,knots=knots,
                                            orderpenalty=orderpenalty,names=names,lambda=lambda,mx=mxz,center=center,
				            ind=ind,xu=xu,zu=zu,iind=iind,type="te",callsm=callsm,indmap=indmap))
				}
			else
				{
				ps1 <- sps(x=x,degree=degree[1],knots=knots[1],orderpenalty=orderpenalty[1])
				ps2 <- sps(x=z,degree=degree[2],knots=knots[2],orderpenalty=orderpenalty[2])
    
     		B <- matrix(0, n, 0)
				bmx <- ps1$basis[mx,]
				bmz <- ps2$basis[mz,]
				mxz <- matrix(0,1,0)
        for (k in 1:ncol(ps1$basis)) 
					{
					B <- cbind(B, ps1$basis[,k]*ps2$basis)
					mxz <- c(mxz,bmx[k]*bmz)
					}
				K <- NULL

				if(!is.null(by))
					{
					if(is.factor(by))
						{
						if(length(levels(by)) < 3)
							by <- (by==as.integer(levels(by)[2]))*1
						else
							stop("Too many levels in ste() argument by, see function mu()!")
						}
					B <- B*by
					}
				
				return(list(basis=B,by=by,K=K,x=x,z=z,degree=degree,knots=knots,
					    orderpenalty=orderpenalty,names=names,mx=mxz,center=center,
				            ind=ind,xu=xu,zu=zu,iind=iind,type="te",callsm=callsm,indmap=indmap))
				}
			}
		}
	}


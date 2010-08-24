ma <- function(x,z, by = NULL, knots=20, c=11.75, kappa=2, loc=NULL, lambda = NULL, center = TRUE, dgts = TRUE, map = NULL, ind = NULL)
	{
	if(!is.null(map))
		{
		cents <- matrix(0,nrow=length(map),ncol=2)
		#binds <- NULL
		for(i in 1:length(map))
			{
			cents[i,] <- centroidpos(map[[i]])
			#binds <- rbind(binds,map[[i]])
			}
		#ind <- as.factor(ind)
		#Z<-diag(nlevels(ind))[ind,]
		#x <- as.vector(Z%*%cents[,1])
		#z <- as.vector(Z%*%cents[,2])
		x <- as.vector(cents[ind,1])
		z <- as.vector(cents[ind,2])
		#x[x==min(x)] <- min(binds[,1])
		#x[x==max(x)] <- max(binds[,1])
		#y[y==min(y)] <- min(binds[,2])
		#y[y==max(y)] <- max(binds[,2])
		}

	if(missing(x) || length(x) == 0)
		stop("No variable defined in ma()!")
	
	call <- match.call(expand.dots=FALSE)
	initcall <- all.names(call)
	callsm <- as.character(as.expression(sys.call(0)))
	callsm <- paste(strsplit(callsm," ","")[[1]],collapse="")
	names <- initcall[initcall!="ma"]
	indmap <- NULL
	if(!is.null(map))
		{
		names <- c(paste(names[length(names)],".xco",sep=""),paste(names[length(names)],".yco",sep=""))
		indmap <- ind
		}
	mx <- approxm(x)
	
	if(missing(z))
		{
		if(dgts != FALSE || dgts != F)
			{
			xn <- get.unique(x,digits=dgts)
			ind <- xn$ind
			iind <- xn$iind
			xu <- xn$xu
			}
		else
			{
			xu <- x
			ind <- iind <- 1:length(xu)
			}
		if(is.null(knots) && is.null(loc))
			{
			B <- rdist(x)
			phi <- max(B)/c
			basis <- Matern(B, range=phi, scale=c, smoothness=kappa)
			dim <- 1
			
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function mu()!")
					}
				basis <- basis*by
				}

			return(list(basis=basis,by=by,K=basis,dim=dim,x=x,kappa=kappa,names=names,lambda=lambda,c=c,mx=basis[mx,],center=center,
				    ind=ind,xu=xu,iind=iind,knots=knots,loc=loc,type="ma",callsm=callsm))
			}
		if(!is.null(knots) && is.null(loc))
			{
			n <- length(x)
			if(knots < 0 || knots >= n)
				stop("Inappropriate number of knots selected!")
			
			get <- collect(1)
			if(is.null(get$response))
				{
				B <- rdist(x)
				phi <- max(B)/c
				basis <- Matern(B, range=phi, scale=c, smoothness=kappa)
				K <- basis
				kn <- knots
				warning("No y variable available for finding optimal design points, used all points instead!")
				}
			else
				{		
				yx <- cbind(get$response,x)
				kn <- cover.design(R=yx, nd=knots)
				kn <- sort(kn[,2])	
				# kn <- seq(min(x), max(x), length=knots)
				B <- rdist(x,kn)
				K <- rdist(kn)

				phi1 <- max(B)/c
				phi2 <- max(K)/c
				basis <- Matern(B, range=phi1, scale=c, smoothness=kappa)
				K <- Matern(K, range=phi2, scale=c, smoothness=kappa)
				}
			dim <- 1
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function m()!")
					}
				basis <- basis*by
				}

			return(list(basis=basis,by=by,K=K,loc=kn,dim=dim,x=x,kappa=kappa,names=names,lambda=lambda,c=c,mx=basis[mx,],center=center,
				    ind=ind,xu=xu,iind=iind,knots=knots,type="ma",callsm=callsm))
			}
		if(!is.null(loc))
			{
			loc <- as.vector(loc)
			B <- rdist(x,loc)
			K <- rdist(loc)

			phi1 <- max(B)/c
			phi2 <- max(K)/c
			basis <- Matern(B, range=phi1, scale=c, smoothness=kappa)
			K <- Matern(K, range=phi2, scale=c, smoothness=kappa)
			dim <- 1
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function m()!")
					}
				basis <- basis*by
				}

			return(list(basis=basis,by=by,K=K,loc=loc,dim=dim,x=x,kappa=kappa,names=names,lambda=lambda,c=c,mx=basis[mx,],center=center,
				    ind=ind,xu=xu,iind=iind,knots=knots,type="ma",callsm=callsm))
			}
		}
	else
		{
		n <- length(x)
		m <- length(z)
		if(n != m)
			stop("Variable lengths differing")
		if(dgts != FALSE || dgts != F)
			{
			X <- cbind(x,z)
			xn <- get.unique(X,digits=dgts)
			ind <- xn$ind
			iind <- xn$iind
			xu <- xn$xu[,1]
			zu <- xn$xu[,2]
			X <- cbind(xu,zu)
			}
		else
			{
			xu <- x
			zu <- z
			X <- cbind(x,z)
			xn <- list(xu=unique(X))
			ind <- iind <- 1:n
			}
		if(is.null(knots) && is.null(loc))
			{
			basis <- matrix(0,nrow=nrow(X),ncol=nrow(X))
			for(j in 1:nrow(X))
				{
				knl <- matrix(X[j,],nrow=1,ncol=2)
				tmpdist <- rdist(X,knl)
				phi <- max(tmpdist)/c
				basis[,j] <- Matern(tmpdist, range=phi, scale=c, smoothness=kappa)
				}	
			dim <- 2
			mxz <- mcenter(x,z,c,kappa,X)
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function mu()!")
					}
				basis <- basis*by
				}
			return(list(basis=basis,by=by,K=basis,dim=dim,x=x,z=z,kappa=kappa,names=names,lambda=lambda,c=c,mx=mxz,center=center,
				    ind=ind,xu=xu,zu=zu,iind=iind,loc=loc,knots=knots,type="ma",callsm=callsm,indmap=indmap))
			}
		if(!is.null(knots) && is.null(loc))
			{
			if(knots < 0 || knots >= n)
				stop("Inappropriate number of knots selected!")
			if(knots == n)
				{
				warning("Number of knots is too large in ma(), set to n-1!")
				knots <- knots - 1
				}
			options(warn = -1)
			options(verbose = TRUE)
		
			kn <- cover.design(R=xn$xu, nd=knots)
			kn <- kn[,1:2]
			
			K <- rdist(kn,kn)
			phi2 <- max(K)/c

			# phi <- max(di)/c
			basis <- matrix(0,nrow=nrow(X),ncol=knots)
			for(j in 1:knots)
				{
				knl <- matrix(kn[j,],nrow=1,ncol=2)
				tmpdist <- rdist(X,knl)
				phi <- max(tmpdist)/c
				basis[,j] <- Matern(tmpdist, range=phi, scale=c, smoothness=kappa)
				}	
			K <- Matern(K, range=phi2, scale=c, smoothness=kappa)

			dim <- 2
			mxz <- mcenter(x,z,c,kappa,kn)

			options(warn = 1)
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function mu()!")
					}
				basis <- basis*by
				}
			return(list(basis=basis,by=by,K=K,loc=kn,dim=dim,x=x,z=z,kappa=kappa,names=names,lambda=lambda,c=c,mx=mxz,center=center,
				    ind=ind,xu=xu,zu=zu,iind=iind,knots=knots,type="ma",callsm=callsm,indmap=indmap))
			}
		if(!is.null(loc))
			{
			if(!is.matrix(loc))
				stop("Argument loc must be a matrix!")

			K <- rdist(loc,loc)
			phi2 <- max(K)/c
			basis <- matrix(0,nrow=nrow(X),ncol=nrow(loc))
			for(j in 1:nrow(loc))
				{
				knl <- matrix(loc[j,],nrow=1,ncol=2)
				tmpdist <- rdist(X,knl)
				phi <- max(tmpdist)/c
				basis[,j] <- Matern(tmpdist, range=phi, scale=c, smoothness=kappa)
				}	
			K <- Matern(K, range=phi2, scale=c, smoothness=kappa)

			dim <- 2
			mxz <- mcenter(x,z,c,kappa,loc)
			if(!is.null(by))
				{
				if(is.factor(by))
					{
					if(length(levels(by)) < 3)
						by <- (by==as.integer(levels(by)[2]))*1
					else
						stop("Too many levels in ma() argument by, see function mu()!")
					}
				basis <- basis*by
				}

			return(list(basis=basis,by=by,K=K,loc=loc,dim=dim,x=x,z=z,kappa=kappa,names=names,lambda=lambda,c=c,mx=mxz,center=center,
				    ind=ind,xu=xu,zu=zu,iind=iind,knots=knots,type="ma",callsm=callsm,indmap=indmap))
			}
		}
	}

mrf <- function(ind, neighbor, by = NULL, lambda = NULL, center = TRUE)
	{
	if(missing(ind))
		stop("Location variable missing!")
	if(missing(neighbor) && !is.null(lambda))
		stop("Neighbor information matrix missing!")

	call <- match.call(expand.dots=FALSE)
	initcall <- all.names(call)
	callsm <- as.character(as.expression(sys.call(0)))
	callsm <- paste(strsplit(callsm," ","")[[1]],collapse="")
	names <- initcall[initcall!="mrf"]

	oind <- ind

	if(!is.factor(ind))
		ind <- as.factor(ind)
	nl <- nlevels(ind)
	Z <- diag(nl)[ind,]

	if(!missing(neighbor))
		if(nl < ncol(neighbor))
			{
			warning("Need to adjust neighbormatrix, more neighbors than different levels in index vector argument in mrf()!")
			nn <- colnames(neighbor)
			nnu <- !duplicated(nn)
			neighbor <- neighbor[nnu,nnu]
			if(ncol(neighbor)!=nl)
				stop("Neighbormatrix in mrf() is not compatible with index vector, cannot create appropriate design matrix, please check for uniqueness of regionnames!")
			}

	if(missing(neighbor) && is.null(lambda))
		neighbor <- matrix(0,nl,nl)
	else
		{
		neighbor <- neighbor*(-1)
		diag(neighbor) <- abs(rowSums(neighbor))
		}
	if(!is.null(by))
		{
		if(is.factor(by))
			{
			if(length(levels(by)) < 3)
				by <- (by==as.integer(levels(by)[2]))*1
			else
				stop("Too many levels in mrf() argument by, see function mu()!")
			}
		Z <- Z*by
		}
	
	IU <- get.unique(as.integer(ind),digits=0)
	ind <- IU$ind
	iind <- IU$iind
	iu <- IU$xu

	mx <- incidentm(Z)

	return(list(basis=Z,by=by,K=neighbor,names=names,lambda=lambda,mrfind=ind,mx=mx,center=center,
		    ind=ind,iu=iu,iind=iind,type="mrf",callsm=callsm,x=oind))
  	}

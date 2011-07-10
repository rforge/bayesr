sps <- function(x, by = NULL, degree = 3, knots = 20, orderpenalty = 2, lambda = NULL, center = TRUE, dgts = TRUE, constr = NULL)
	{
	if(is.null(x))
		stop("No variable specified in ps()")
	if(degree < 0)
		stop("There is no negative degree in ps()!")
	if(knots < 2)
		stop("Need at least 2 knots for basis function construction in ps()")
	if(orderpenalty < 0)
		stop("There is no negative penalty in ps()!")
	constrdo <- NULL
	if(!is.null(constr))
		if(constr!=1 && constr!=2)
			{
			warning("Wrong specified argument constr in ps(), set to default!")
			constr <- NULL
			}
	call <- match.call(expand.dots=FALSE)
	initcall <- all.names(call)
	callsm <- as.character(as.expression(sys.call(0)))
	callsm <- paste(strsplit(callsm," ","")[[1]],collapse="")

	names <- initcall[initcall!="ps"]		
	if(any("$"%in%names))
		names <- paste(names[2],names[1],names[3],sep="")
	if(is.null(by))
		{
		byname <- NULL
		mx <- approxm(x)
		}
	else
		{
		byname <- names[2]
		}
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
		ind <- iind <- 1:length(x)
		}

	if(is.null(by))
		{
		if(orderpenalty == 0)
			{
			basis <- bspline(x, degree, knots)
			K <- NULL

			return(list(basis=basis,K=K,x=x,by=by,names=names,lambda=lambda,mx=basis[mx,],center=center,ind=ind,xu=xu,iind=iind,
					degree=degree,knots=knots,orderpenalty=orderpenalty,type="ps",callsm=callsm,
					byname=byname,constr=constr))
			}
		else
			{
			basis <- bspline(x, degree, knots)
			K <- penalty(orderpenalty,basis)

			return(list(basis=basis,K=K,x=x,by=by,names=names,lambda=lambda,mx=basis[mx,],center=center,ind=ind,xu=xu,iind=iind,
					degree=degree,knots=knots,orderpenalty=orderpenalty,type="ps",callsm=callsm,
					byname=byname,constr=constr))
			}
		}
	if(!is.null(by))
		{
		if(is.factor(by))
			{
			if(length(levels(by)) < 3)
				by <- (by==as.integer(levels(by)[2]))*1
			else
				stop("Too many levels in ps() argument by, see function mu()!")
			}
		# if(all(-1 < by | by < 1))
		#	{
		#	cat("yo!\n")
		#	by <- by + 1
		#	}
		mx <- approxm(x*by)
		if(orderpenalty == 0)
			{
			basis <- bspline(x, degree, knots)
			K <- NULL
			}
		else
			{
			basis <- bspline(x, degree, knots)
			K <- penalty(orderpenalty,basis)
			}

		# by <- diag(t(U)%*%diag(by)%*%U)

		if(is.factor(by))
			basis <- basis*(c(by)-1)
		else
			basis <- basis*by
		return(list(basis=basis,K=K,x=x,by=by,names=names,lambda=lambda,mx=basis[mx,],center=center,ind=ind,xu=xu,iind=iind,
				degree=degree,knots=knots,orderpenalty=orderpenalty,type="ps",callsm=callsm,
				byname=byname,constr=constr))
		}
	}

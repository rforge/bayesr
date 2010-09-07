stewrap <- function(object,data,call)
	{
	ocall <- call
	geosm <- FALSE
	byfun <- FALSE
	if(!is.null(object$xt$geo))
		geosm <- object$xt$geo
	if(object$by!="NA")
		{
		by <- eval(parse(text=object$by),envir=data)
		if(is.factor(by))
			{
			Xok <- byfun(by,object,data,call)
			byfun <- TRUE
			}
		}
	if(object$dim < 2 & !geosm)
		{
		Xtmp <- eval(parse(text = object$term[1]),envir=data)
		if(!is.null(object$xt$dgts))
			{
			xn <- get.unique(Xtmp,digits=object$xt$dgts)
			ind <- xn$ind
			iind <- xn$iind
			xu <- xn$xu
			}
		else
			{
			xu <- Xtmp
			ind <- iind <- 1:length(Xtmp)
			}
		if(object$by!="NA")
			{
			if(!byfun)
				mx <- approxm(Xtmp*by)
			}
		else
			mx <- approxm(Xtmp)
		if(!byfun)
			{
			Xok <- smooth.construct(object,data,NULL)
			mx <- Xok$X[mx,]
			}
		}
	if(object$dim > 1 || geosm)
		{
		m <- NULL
		if(!is.null(object$mp))
			if(object$mp)
				object$mp <- FALSE
		if(!geosm)
			{
			for(i in 1:object$dim)
				m <- cbind(m,eval(parse(text = object$term[i]),envir=data))
			if(!is.null(object$xt$dgts))
				{
				xn <- get.unique(m,digits=object$xt$dgts)
				ind <- xn$ind
				iind <- xn$iind
				}
			else
				ind <- iind <- 1:nrow(m)
			datlist <- list()
			bmx <- "NA"
			if(object$by!="NA")
				{
				if(!byfun)
					{
					bmx <- eval(parse(text=paste(object$by,"[approxm(",object$by,")]")),envir=data)
					datlist$bmx <- c(bmx,rnorm(100,mean=bmx,sd=0.00001))
					}
				}
			for(i in 1:object$dim)
				{
				what <- eval(parse(text=object$term[i]),envir=data)
				eval(parse(text=paste("datlist$mx",i,"<-what[approxm(what)]",sep="")))
				txt <- paste("datlist$mx",i,"<-c(datlist$mx",i,",rnorm(100,mean=datlist$mx",i,",sd=2))",sep="")
				eval(parse(text=txt))
				call <- sub(object$term[i],paste("mx",i,sep=""),call)
				}
			call <- sub(object$by,"bmx",call)
			Xok <- smooth.construct(object,data,NULL)
			Xmx <- eval(parse(text=paste("smooth.construct(",call,",datlist,NULL)")),envir=datlist)
			mx <- Xmx$X[1,] ### needs to be checked !!!!!!!!!
			}
		else
			{
			m <- list()
			for(i in 1:object$dim)
				m[[i]] <- eval(parse(text = object$term[i]),envir=data)
			ind <- iind <- 1:length(m[[1]])
			Xok <- smooth.construct(object,data,NULL)
			}
		}
	if(!byfun)
		{
		if(is.null(object$xt$center))
			object$xt$center <- TRUE
		if(is.null(Xok$mx))
			Xok$mx <- mx
		Xok$ind <- ind
		Xok$center <- object$xt$center
		Xok$iind <- iind
		Xok$constr <- object$xt$constr
		Xok$specs <- object
		Xok$specs$label <- object$label # paste(strsplit(ocall," ","")[[1]],collapse="")
		Xok$specs$xt <- delete.NULLs(Xok$xt)
		Xok$specs$m <- Xok$m
		Xok$specs$call <- ocall
		Xok$specs$dim <- Xok$dim
		Xok$specs$term <- Xok$term
		Xok$specs$data <- data
		attr(Xok,"byfun") <- FALSE
		}

	return(Xok)
	}


# byfun
byfun <- function(by,object,data,call)
	{
	if(!is.null(object$xt$type) & (object$xt$type=="m"))
		{
		bs <-  strsplit(class(object),"",".")[[1]]
		bsc <- NULL
		ok <- TRUE
		for(i in 1:length(bs))
			{
			if(bs[i]==".")
				ok <- FALSE
			if(ok)
				bsc <- c(bsc,bs[i])
			}
		bsc <- paste(bsc,collapse="")
		mterm <- strsplit(object$label,""," ")[[1]]
		mterm <- paste(mterm[1:(length(mterm)-1)],collapse="")
		mterm <- paste(mterm,",k=",object$bs.dim,sep="")
		mterm <- paste(mterm,",m=",object$p.order,sep="")
		xtn <- names(object$xt)
		xtnspecs <- "list("
		for(i in 1:length(xtn))
			if(xtn[i]!="hypident" && xtn[i]!="center" && xtn[i]!="type")
				{
				if(i < length(xtn))
					xtnspecs <- paste(xtnspecs,xtn[i],"=",object$xt[xtn[i]],",",sep="")
				else
					xtnspecs <- paste(xtnspecs,xtn[i],"=",object$xt[xtn[i]],sep="")
				}
		xtnspecs <- paste(xtnspecs,")",sep="")
		mterm <- paste(mterm,",xt=",xtnspecs,sep="")
		mterm <- paste(mterm,",bs=\"",bsc,"\")",sep="")
		center <- object$xt$center
		if(is.null(center))
			center <- TRUE
		hypi<- object$xt$hypident
		if(is.null(hypi))
			hypi <- FALSE
		if(is.null(object$sp))
			Xok <- paste("m(",object$by,",",mterm,",sp=NULL,center=",center,",hypident=",hypi,")",sep="")
		else
			Xok <- paste("m(",object$by,",",mterm,",sp=",object$sp,",center=",center,",hypident=",hypi,")",sep="")
		attr(Xok,"type") <- "m"
		}
	else
		{
		x <- NULL
		for(i in 1:object$dim)
			x <- cbind(x,eval(parse(text=object$term[i]),envir=data))
		nt <- levels(by)
		nl <- nlevels(by)
		bym <- diag(nl)[by,]
		byn <- object$by
		colnames(bym) <- nt
		out <- object
		out$by <- "NA"
		oterm <- object$term
		Xok <- list()
		for(i in 1:length(nt))
			{
			dat <- list()
			out$term <- paste(object$term,".",nt[i],sep="")
			ncall <- call
			label <- out$label
			for(j in 1:length(out$term))
				{
				ncall <- sub(oterm[j],out$term[j],ncall)
				label <- sub(oterm[j],out$term[j],label)
				if(is.matrix(bym))
					{
					eval(parse(text=paste("dat$",out$term[j],"<- x[,j]*bym[,i]",sep="")))
					eval(parse(text=paste("dat$byokcheck <- (bym[,i]==1)",sep="")))
					}
				else
					{
					eval(parse(text=paste("dat$",out$term[j],"<- x[,j]*bym",sep="")))
					eval(parse(text=paste("dat$byokcheck <- (bym==1)",sep="")))
					}
				}
			out$label <- label
			Xok[[i]] <- stewrap(out,dat,ncall)
			}
		attr(Xok,"type") <- "nom"
		}
	attr(Xok,"byfun") <- TRUE

	return(Xok)
	}


# create matern smooth construct
smooth.construct.ma.smooth.spec<-function(object,data,knots)
	{ 
	m <- object$p.order
  	if(is.na(m[1])) 
		m <- c(11.75,2)
 	if(object$bs.dim<0) 
		object$bs.dim <- 20
	geosm <- FALSE
  	loc <- object$xt$loc
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	x <- eval(parse(text=object$term[1]),envir=data)
	if(object$dim<2)
		{
		get <- collect(1)
		if(is.null(get$response))
			loc <- seq(min(x),max(x),length=object$bs.dim)
		X <- ma(x,by=by,knots=object$bs.dim,c=m[1],kappa=m[2], 
   	   	        loc=loc,map=object$xt$map,ind=object$xt$ind)
		}
	else
		{
		eval(parse(text=paste(object$term[1],"<- eval(parse(text=object$term[1]),envir=data)")))
		eval(parse(text=paste(object$term[2],"<- eval(parse(text=object$term[2]),envir=data)")))
		if(!is.null(object$xt$geo))
			geosm <- object$xt$geo
		if(!geosm)
			X <- eval(parse(text=paste("ma(",object$term[1],",",object$term[2],",by=by,knots=object$bs.dim,c=m[1],kappa=m[2],loc=loc)")))
		else
			X <- eval(parse(text=paste("ma(by=by,knots=object$bs.dim,c=m[1],kappa=m[2],map=",object$term[2],",ind=",object$term[1],")")))
		}
  	object$X <- X$basis # the finished model matrix
  	object$rank <- qr(X$K)$rank  # penalty rank
	object$S <- list(X$K)
  	object$null.space.dim <- ncol(X$basis)  # dim. of unpenalized space

  	# store "ma" specific stuff ...
	if(geosm)
		{
		object$term <- X$names
		txt <- paste("list(",X$names[1],"=X$x,",X$names[2],"=X$z)",sep="")
		tmp <- eval(parse(text=txt))
		tmp$id <- X$indmap
		assign("geospline.centroid.data",tmp,envir=.GlobalEnv)
		}
	else 
		geosm <- NULL
	xt <- object$xt
	if(is.null(xt))
 		xt<-list(loc=X$loc,geo=geosm)
	else
		{
		xt$loc <- X$loc
		xt$geo <- geosm
		}
	object$xt <- xt
	object$m<-m
	object$mx <- X$mx
  	object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  	class(object)<-"ma.smooth"  # Give object a class

  	return(object)
	}

Predict.matrix.ma.smooth<-function(object,data)
	{ 
  	m <- object$m;     # spline order (3=cubic)
  	loc<-object$xt$loc    # knot locations
	x <- eval(parse(text=object$term[1]),envir=data)
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	if(object$dim<2)
		{
		X <- ma(x,by=by,knots=object$bs.dim,c=m[1],kappa=m[2], 
   	   	        loc=loc,map=object$xt$map,ind=object$xt$ind)
		}
	else
		{
		z <- eval(parse(text=object$term[2]),envir=data)
		X <- ma(x,z,by=by,knots=object$bs.dim,c=m[1],kappa=m[2], 
   	   	        loc=loc,map=object$xt$map,ind=object$xt$ind)
		}
	return(X$basis)
	}


# create own tensor smooth construct
smooth.construct.ste.smooth.spec<-function(object,data,knots)
	{ 
	m <- object$p.order
  	if(is.na(m[1])) 
		m <- list(degree=c(3,3),orderpenalty=c(2,2))
	if(length(m$degree)<2)
		m$degree <- rep(m$degree,2)
	if(length(m$orderpenalty)<2)
		m$orderpenalty <- rep(m$orderpenalty,2)
 	if(object$bs.dim<0) 
		object$bs.dim <- c(20,20)
	if(length(object$bs.dim)<2)
		object$bs.dim <- rep(object$bs.dim,2)
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	eval(parse(text=paste(object$term[1],"<- eval(parse(text=object$term[1]),envir=data)")))
	eval(parse(text=paste(object$term[2],"<- eval(parse(text=object$term[2]),envir=data)")))
	geosm <- FALSE
	if(!is.null(object$xt$geo))
		geosm <- object$xt$geo
	if(!geosm)
		X <- eval(parse(text=paste("ste(",object$term[1],",",object$term[2],",by=by,degree=m$degree,knots=object$bs.dim,orderpenalty=m$orderpenalty)")))
	else
		X <- eval(parse(text=paste("ste(by=by,degree=m$degree,knots=object$bs.dim,orderpenalty=m$orderpenalty,map=",object$term[2],",ind=",object$term[1],")")))
  	object$X <- X$basis # the finished model matrix
  	object$rank<-qr(X$K)$rank  # penalty rank
	object$S <- list(X$K)
  	object$null.space.dim <- ncol(X$basis)  # dim. of unpenalized space
  	# store "ste" specific stuff ...
	if(object$dim<2)
		object$dim <- 2
	if(geosm)
		{
		object$term <- X$names
		tmp <- eval(parse(text=paste("list(",X$names[1],"=X$x,",X$names[2],"=X$z)",sep="")))
		tmp$id <- X$indmap
		assign("geospline.centroid.data",tmp,envir=.GlobalEnv)
		}
	object$m <- m	
	object$mx <- X$mx
	object$x <- cbind(X$x,X$z) # !!!!!!! check for geospline
  	object$df<-ncol(object$X)     # maximum DoF (if unconstrained)
  	class(object)<-"ste.smooth"  # Give object a class
  	return(object)
	}

Predict.matrix.ste.smooth<-function(object,data)
	{ 
  	m <- object$m;     # spline order (3=cubic)
  	loc <- object$xt$loc    # knot locations
	x <- eval(parse(text=object$term[1]),envir=data)
	z <- eval(parse(text=object$term[2]),envir=data)
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	X <- ste(x,z,by=by,degree=m$degree,knots=object$bs.dim,orderpenalty=m$orderpenalty, 
		 map=object$xt$map,ind=object$xt$ind)
	return(X$basis)
	}


# create mrf smooth construct
smooth.construct.mrf.smooth.spec<-function(object,data,knots)
	{ 
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by[1]))
		by <- NULL
	ind <- eval(parse(text=object$term),envir=data)
	neighbor <- object$xt$neighbor
	if(is.null(neighbor))
		neighbor <- object$xt[[1]]
	X <- mrf(ind=ind,neighbor=neighbor,by=by)
  	object$X <- X$basis # the finished model matrix
  	object$rank <- qr(X$K)$rank  # penalty rank
	object$S <- list(X$K)
  	object$null.space.dim <- ncol(X$basis)  # dim. of unpenalized space
  	object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
	object$mx <- X$mx
  	class(object)<-"mrf.smooth"  # Give object a class
  	return(object)
	}

Predict.matrix.mrf.smooth<-function(object,data)
	{ 
	ind <- eval(parse(text=object$term),envir=data)
	neighbor <- object$xt$neighbor
	if(is.null(neighbor))
		neighbor <- object$xt[[1]]
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	X <- mrf(ind=ind,neighbor=neighbor,by=by)
	return(X$basis)
	}


# create my ps smooth construct
smooth.construct.sps.smooth.spec<-function(object,data,knots)
	{ 
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	x <- eval(parse(text=object$term),envir=data)
	m <- object$p.order
  	if(is.na(m)) 
		m <- list(degree=3,orderpenalty=2)
	if(object$bs.dim<0)
		object$bs.dim <- 20
	X <- sps(x,by=by,degree=m[[1]],knots=object$bs.dim,orderpenalty= m[[2]])
  	object$X <- X$basis # the finished model matrix
  	object$rank <- qr(X$K)$rank  # penalty rank
	object$S <- list(X$K)
  	object$null.space.dim <- ncol(X$basis)  # dim. of unpenalized space
  	object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
	object$mx <- X$mx
	object$m <- m
  	class(object)<-"sps.smooth"  # Give object a class
  	return(object)
	}

Predict.matrix.sps.smooth<-function(object,data)
	{ 
	x <- eval(parse(text=object$term),envir=data)
	by <- eval(parse(text=object$by),envir=data)
	if(is.na(by))
		by <- NULL
	X <- sps(x,by=by,degree=object$m[[1]],knots=object$bs.dim,orderpenalty=object$m[[2]])
	return(X$basis)
	}

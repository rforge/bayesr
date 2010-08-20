predict.gibbs <- function(object, which = 1, newdata = NULL, ...)
	{
	if(is.null(which))
		return(fitted(object))
	if(length(object$terms) < which[1] || which[1] == 0)
		stop("Argument which is specified wrong, nothing to predict!")
	call <- match.call()
	if(!is.null(newdata))
		{
		ncall <- all.names(call)
		nname <- ncall[length(ncall)]
		}
	X <- object$fout[[which[1]]]
	type <- attr(X,"term.type")
	if(type == "smooth")
		{
		if(inherits(X,"ps.gibbs"))
			{
			if(!is.null(newdata))
				pred <- pspredict(X,newdata,nname)
			else
				return(X)
			}
		if(inherits(X,"te.gibbs"))
			{
			if(!is.null(newdata))
				pred <- tepredict(X,newdata,nname)
			else
				return(X)
			}
		if(inherits(X,"ma.gibbs"))
			{
			if(!is.null(newdata))
				pred <- mapredict(X,newdata,nname)
			else
				return(X)
			}
		if(inherits(X,"mrf.gibbs"))
			return(X)
		}
	if(type == "linear")
		{
		if(is.null(attr(X,"factorcheck")))
			attr(X,"factorcheck") <- "nonfactor"
		if(attr(X,"factorcheck") == "nonfactor")
			{
			if(!is.null(newdata))
				pred <- linearpredict(X,newdata,nname)
			else
				return(X)
			}
		else
			return(X)
		}
	if(type == "mu")
		{
		if(!is.null(newdata))
			{
			nl <- length(X)
			pred <- vector("list",nl)

			if(inherits(X[[1]],"ps.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- pspredict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			if(inherits(X[[1]],"te.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- tepredict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			if(inherits(X[[1]],"ma.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- mapredict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			if(inherits(X[[1]],"mrf.gibbs"))
				return(X)
			if(inherits(X[[1]],"linear.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- linearpredict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			attr(pred,"term.type") <- "mu"
			attr(pred,"par.mat") <- attr(X,"par.mat")
			attr(pred,"mu.term") <- attr(X,"mu.term")
		      class(pred) <- "gibbs"
			}
		else
			return(X)
		}
	if(type == "random")
		{
		effplot <- FALSE
		if(!is.null(X$terms))
			{
			if(length(which) > 1)
				{
				if(length(which) == 2)
					{
					if(which[2] == 0)
						which[2] <- 1
					s2r <- attr(X$effects,"random.variance.draws")
					X <- X$terms
					if(which[2] > length(X))
						stop("Argument which is specified wrong, nothing to predict!")
					X <- X[[which[2]]]

					if(attr(X,"term.type") == "smooth")
						{
						if(inherits(X,"ps.gibbs"))
							{
							if(!is.null(newdata))
								pred <- pspredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"te.gibbs"))
							{
							if(!is.null(newdata))
								pred <- tepredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"ma.gibbs"))
							{
							if(!is.null(newdata))
								pred <- mapredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"mrf.gibbs"))
							return(X)
						}
					if(attr(X,"term.type") == "linear")
						{
						if(is.null(attr(X,"factorcheck")))
							attr(X,"factorcheck") <- "nonfactor"
						if(attr(X,"factorcheck") == "nonfactor")
							{
							if(!is.null(newdata))
								pred <- linearpredict(X,newdata,nname)
							else
								return(X)
							}
						else
							return(X)
						}	
					if(attr(X,"term.type") == "random")
						return(X)
					}
				if(length(which) == 3)
					{
					if(which[2] == 0)
						stop("Argument which is specified wrong, nothing to predict!")
					if(which[3] == 0)
						which[3] <- 1
					X <- X$terms
					X <- X[[which[2]]]
					s2r <- attr(X$effects,"random.variance.draws")
					X <- X$terms
					if(which[3] > length(X))
						stop("Argument which is specified wrong, nothing to do!")
					X <- X[[which[3]]]
					if(attr(X,"term.type") == "smooth")
						{
						if(inherits(X,"ps.gibbs"))
							{
							if(!is.null(newdata))
								pred <- pspredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"te.gibbs"))
							{
							if(!is.null(newdata))
								pred <- tepredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"ma.gibbs"))
							{
							if(!is.null(newdata))
								pred <- mapredict(X,newdata,nname)
							else
								return(X)
							}
						if(inherits(X,"mrf.gibbs"))
							return(X)
						}
					if(attr(X,"term.type") == "linear")
						{
						if(is.null(attr(X,"factorcheck")))
							attr(X,"factorcheck") <- "nonfactor"
						if(attr(X,"factorcheck") == "nonfactor")
							{
							if(!is.null(newdata))
								pred <- linearpredict(X,newdata,nname)
							else
								return(X)
							}
						else
							return(X)
						}	
					}
				}
			else
				effplot <- TRUE
			}
		if(length(which) == 1)
			effplot <- TRUE
		if(effplot)
			return(X)
		}
	return(invisible(pred))
	}


pspredict <- function(x,newdata,name)
	{
	specs <- attr(x,"smooth.specs")
	coefs <- attr(x,"smooth.coef.draws")
    	minx <- min(x[,1]) - 0.001
    	maxx <- max(x[,1]) + 0.001
    	step <- (maxx - minx)/(specs$knots - 1)
    	k <- seq(minx - specs$degree * step, maxx + specs$degree * step, by = step)

	bas <- spline.des(knots = k, newdata, 
                        ord = (specs$degree + 1), 
                        outer.ok = TRUE)$design
	nd <- ncol(coefs)
	nf <- matrix(0,length(newdata),nd)
	rq <- attr(x,"rq")
	for(i in 1:nd)
		nf[,i] <- bas%*%(rq%*%coefs[,i])
	pred <- matrix(0,length(newdata),10)
	pred[,1] <- newdata
	pred[,2] <- apply(nf,1,mean)
	pred[,3] <- apply(nf,1,quantile,probs=0.025)
	pred[,4] <- apply(nf,1,quantile,probs=0.1)
	pred[,5] <- apply(nf,1,quantile,probs=0.5)
	pred[,6] <- apply(nf,1,quantile,probs=0.9)
	pred[,7] <- apply(nf,1,quantile,probs=0.975)
      pred[,9] <- (pred[,4] < 0 & pred[,6] < 0) * (-1) + (pred[,4] <= 0 & pred[,6] >= 0) * 0 + (pred[,4] > 0 & pred[,6] > 0) * 1
      pred[,10] <- (pred[,3] < 0 & pred[,7] < 0) * (-1) + (pred[,3] <= 0 & pred[,7] >= 0) * 0 + (pred[,3] > 0 & pred[,7] > 0) * 1
	pred <- pred[order(pred[,1]),]
	cnams <- colnames(x)
	cnams[1] <- name
	colnames(pred) <- cnams

	attr(pred,"smooth.specs") <- attr(x,"smooth.specs")
	attr(pred,"smooth.specs")$names <- name
	attr(pred,"smooth.coef") <- attr(x,"smooth.coef")
	attr(pred,"smooth.prediction") <- TRUE
	attr(pred,"smooth.edf") <- attr(x,"smooth.edf")
	attr(pred,"smooth.ceffect") <- attr(x,"smooth.ceffect")
	attr(pred,"smooth.ceffect.draws") <- attr(x,"smooth.ceffect.draws")
	attr(pred,"smooth.coef.draws") <- attr(x,"smooth.coef.draws")
	attr(pred,"smooth.variance.draws") <- attr(x,"smooth.variance.draws")
	attr(pred,"smooth.coef.mean") <- attr(x,"smooth.coef.mean")

	class(pred) <- "ps.gibbs"
	return(pred)
	}

tepredict <- function(x,newdata,name)
	{
	if(!is.matrix(newdata))
		stop("Argument newdata is not a matrix!")
	if(ncol(newdata) != 2)
		stop("Wrong number of columns in newdata!")

	specs <- attr(x,"smooth.specs")
	coefs <- attr(x,"smooth.coef.draws")

    	minx <- min(x[,1]) - 0.001
    	maxx <- max(x[,1]) + 0.001
    	step1 <- (maxx - minx)/(specs$knots[1] - 1)
    	k1 <- seq(minx - specs$degree[1] * step1, maxx + specs$degree[1] * step1, by = step1)

	bas1 <- spline.des(knots = k1, newdata[,1], 
                         ord = (specs$degree[1] + 1), 
                         outer.ok = TRUE)$design

    	minz <- min(x[,2]) - 0.001
    	maxz <- max(x[,2]) + 0.001
    	step2 <- (maxz - minz)/(specs$knots[2] - 1)
    	k2 <- seq(minz - specs$degree[2] * step2, maxz + specs$degree[2] * step2, by = step2)

	bas2 <- spline.des(knots = k2, newdata[,2], 
                         ord = (specs$degree[2] + 1), 
                         outer.ok = TRUE)$design

     	bas <- matrix(0, nrow(newdata), 0)
     	for (k in 1:ncol(bas1)) 
		bas <- cbind(bas, bas1[,k]*bas2)

	nd <- ncol(coefs)
	nf <- matrix(0,nrow(newdata),nd)
	rq <- attr(x,"rq")

	for(i in 1:nd)
		nf[,i] <- bas%*%(rq%*%coefs[,i])

	pred <- matrix(0,nrow(newdata),11)
	pred[,1] <- newdata[,1]
	pred[,2] <- newdata[,2]
	pred[,3] <- apply(nf,1,mean)
	pred[,4] <- apply(nf,1,quantile,probs=0.025)
	pred[,5] <- apply(nf,1,quantile,probs=0.1)
	pred[,6] <- apply(nf,1,quantile,probs=0.5)
	pred[,7] <- apply(nf,1,quantile,probs=0.9)
	pred[,8] <- apply(nf,1,quantile,probs=0.975)
      pred[,10] <- (pred[,5] < 0 & pred[,7] < 0) * (-1) + (pred[,5] <= 0 & pred[,7] >= 0) * 0 + (pred[,5] > 0 & pred[,7] > 0) * 1
      pred[,11] <- (pred[,4] < 0 & pred[,8] < 0) * (-1) + (pred[,4] <= 0 & pred[,8] >= 0) * 0 + (pred[,4] > 0 & pred[,8] > 0) * 1
	pred <- pred[order(pred[,1]),]
	cnams <- colnames(x)
	if(is.null(colnames(newdata)))
		{
		cnams[1] <- paste(name,":",cnams[1],sep="")
		cnams[2] <- paste(name,":",cnams[2],sep="")
		}
	else
		{
		cnams[1] <- colnames(newdata)[1]
		cnams[2] <- colnames(newdata)[2]
		}
	colnames(pred) <- cnams
	attr(pred,"smooth.specs") <- attr(x,"smooth.specs")
	attr(pred,"smooth.specs")$names <- cnams[1:2]
	attr(pred,"smooth.coef") <- attr(x,"smooth.coef")
	attr(pred,"smooth.prediction") <- TRUE
	attr(pred,"knots") <- list(k1,k2)
	attr(pred,"smooth.edf") <- attr(x,"smooth.edf")
	attr(pred,"smooth.ceffect") <- attr(x,"smooth.ceffect")
	attr(pred,"smooth.ceffect.draws") <- attr(x,"smooth.ceffect.draws")
	attr(pred,"smooth.coef.draws") <- attr(x,"smooth.coef.draws")
	attr(pred,"smooth.variance.draws") <- attr(x,"smooth.variance.draws")
	attr(pred,"smooth.coef.mean") <- attr(x,"smooth.coef.mean")

	class(pred) <- "te.gibbs"
	return(pred)
	}

mapredict <- function(x,newdata,name)
	{
	specs <- attr(x,"smooth.specs")
	coefs <- attr(x,"smooth.coef.draws")
	if(specs$dim == 1)
		{		
		bas <- ma(newdata,c=specs$c,kappa=specs$kappa,loc=specs$loc,dgts=FALSE)$basis
            #matplot(x=newdata[order(newdata)],y=bas[order(newdata),],type="l")
		nd <- ncol(coefs)
		nf <- matrix(0,length(newdata),nd)
		rq <- attr(x,"rq")
		for(i in 1:nd)
			nf[,i] <- bas%*%(rq%*%coefs[,i])
		pred <- matrix(0,length(newdata),10)
		pred[,1] <- newdata
		pred[,2] <- apply(nf,1,mean)
		pred[,3] <- apply(nf,1,quantile,probs=0.025)
		pred[,4] <- apply(nf,1,quantile,probs=0.1)
		pred[,5] <- apply(nf,1,quantile,probs=0.5)
		pred[,6] <- apply(nf,1,quantile,probs=0.9)
		pred[,7] <- apply(nf,1,quantile,probs=0.975)
      	pred[,9] <- (pred[,4] < 0 & pred[,6] < 0) * (-1) + (pred[,4] <= 0 & pred[,6] >= 0) * 0 + (pred[,4] > 0 & pred[,6] > 0) * 1
      	pred[,10] <- (pred[,3] < 0 & pred[,7] < 0) * (-1) + (pred[,3] <= 0 & pred[,7] >= 0) * 0 + (pred[,3] > 0 & pred[,7] > 0) * 1
		pred <- pred[order(pred[,1]),]
		cnams <- colnames(x)
		cnams[1] <- name
		colnames(pred) <- cnams
		attr(pred,"smooth.specs") <- attr(x,"smooth.specs")
		attr(pred,"smooth.specs")$names <- name
		attr(pred,"smooth.specs")$dim <- 1
		attr(pred,"smooth.coef") <- attr(x,"smooth.coef")
		attr(pred,"smooth.prediction") <- TRUE
		attr(pred,"smooth.edf") <- attr(x,"smooth.edf")
		attr(pred,"smooth.ceffect") <- attr(x,"smooth.ceffect")
		attr(pred,"smooth.ceffect.draws") <- attr(x,"smooth.ceffect.draws")
		attr(pred,"smooth.coef.draws") <- attr(x,"smooth.coef.draws")
		attr(pred,"smooth.variance.draws") <- attr(x,"smooth.variance.draws")
		attr(pred,"smooth.coef.mean") <- attr(x,"smooth.coef.mean")

		class(pred) <- "ma.gibbs"
		}
	else
		{
		if(!is.matrix(newdata))
			stop("Argument newdata is not a matrix!")
		if(ncol(newdata) != 2)
			stop("Wrong number of columns in newdata!")
		bas <- ma(newdata[,1],newdata[,2],c=specs$c,kappa=specs$kappa,loc=specs$loc,dgts=FALSE)$basis
		nd <- ncol(coefs)
		nf <- matrix(0,nrow(newdata),nd)
		rq <- attr(x,"rq")

		for(i in 1:nd)
			nf[,i] <- bas%*%(rq%*%coefs[,i])

		pred <- matrix(0,nrow(newdata),11)
		pred[,1] <- newdata[,1]
		pred[,2] <- newdata[,2]
		pred[,3] <- apply(nf,1,mean)
		pred[,4] <- apply(nf,1,quantile,probs=0.025)
		pred[,5] <- apply(nf,1,quantile,probs=0.1)
		pred[,6] <- apply(nf,1,quantile,probs=0.5)
		pred[,7] <- apply(nf,1,quantile,probs=0.9)
		pred[,8] <- apply(nf,1,quantile,probs=0.975)
      	pred[,10] <- (pred[,5] < 0 & pred[,7] < 0) * (-1) + (pred[,5] <= 0 & pred[,7] >= 0) * 0 + (pred[,5] > 0 & pred[,7] > 0) * 1
      	pred[,11] <- (pred[,4] < 0 & pred[,8] < 0) * (-1) + (pred[,4] <= 0 & pred[,8] >= 0) * 0 + (pred[,4] > 0 & pred[,8] > 0) * 1
		pred <- pred[order(pred[,1]),]
		cnams <- colnames(x)
		if(is.null(colnames(newdata)))
			{
			cnams[1] <- paste(name,":",cnams[1],sep="")
			cnams[2] <- paste(name,":",cnams[2],sep="")
			}
		else
			{
			cnams[1] <- colnames(newdata)[1]
			cnams[2] <- colnames(newdata)[2]
			}
		colnames(pred) <- cnams
		attr(pred,"smooth.specs") <- attr(x,"smooth.specs")
		attr(pred,"smooth.specs")$names <- cnams[1:2]
		attr(pred,"smooth.coef") <- attr(x,"smooth.coef")
		attr(pred,"smooth.prediction") <- TRUE
		attr(pred,"smooth.edf") <- attr(x,"smooth.edf")
		attr(pred,"smooth.ceffect") <- attr(x,"smooth.ceffect")
		attr(pred,"smooth.ceffect.draws") <- attr(x,"smooth.ceffect.draws")
		attr(pred,"smooth.coef.draws") <- attr(x,"smooth.coef.draws")
		attr(pred,"smooth.variance.draws") <- attr(x,"smooth.variance.draws")
		attr(pred,"smooth.coef.mean") <- attr(x,"smooth.coef.mean")

		class(pred) <- "ma.gibbs"
		return(pred)
		}
	return(pred)
	}

linearpredict <- function(x,newdata,name)
	{
	coefs <- attr(x,"linear.coef")
	pred <- matrix(0,length(newdata),10)
	pred[,1] <- newdata
	pred[,2] <- newdata*coefs[1]
	pred[,3] <- newdata*coefs[3]
	pred[,4] <- newdata*coefs[4]
	pred[,5] <- newdata*coefs[5]
	pred[,6] <- newdata*coefs[6]
	pred[,7] <- newdata*coefs[7]
      pred[,9] <- (pred[,4] < 0 & pred[,6] < 0) * (-1) + (pred[,4] <= 0 & pred[,6] >= 0) * 0 + (pred[,4] > 0 & pred[,6] > 0) * 1
      pred[,10] <- (pred[,3] < 0 & pred[,7] < 0) * (-1) + (pred[,3] <= 0 & pred[,7] >= 0) * 0 + (pred[,3] > 0 & pred[,7] > 0) * 1
	pred <- pred[order(pred[,1]),]
	cnams <- colnames(x)
	cnams[1] <- name
	colnames(pred) <- cnams
	attr(pred,"smooth.specs") <- list(names=name)
	attr(pred,"linear.coef.draws") <- attr(x,"linear.coef.draws")
	attr(pred,"linear.ceffect") <- attr(x,"linear.ceffect") 
	attr(pred,"linear.coef") <- attr(x,"linear.coef")
	attr(pred,"linear.coef.mean") <- attr(pred,"linear.coef.mean")
	attr(pred,"linear.prediction") <- TRUE

	class(pred) <- "linear.gibbs"
	return(pred)
	}


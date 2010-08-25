predict.gibbs <- function(object, which = 1, newdata = NULL, ...)
	{
	if(is.null(which))
		return(fitted(object))
	if(length(object$terms) < which[1] || which[1] == 0)
		stop("Argument which is specified wrong, nothing to predict!")
	call <- match.call()
	ok <- FALSE
	if(!is.null(newdata))
		{
		if((!is.list(newdata)) && (!is.data.frame(newdata)))
			stop("Argument newdata must be a list or data.frame()!")
		else
			ok <- TRUE
		}
	X <- object$fout[[which[1]]]
	type <- attr(X,"term.type")
	if(type == "smooth")
		{
		if(ok)
			{
			if(inherits(X,"sm.gibbs"))
				pred <- sm.predict(X,newdata)
			else
				return(X)
			}
		else
			return(X)
		}
	if(type == "linear")
		{
		if(is.null(attr(X,"factorcheck")))
			attr(X,"factorcheck") <- "nonfactor"
		if(attr(X,"factorcheck") == "nonfactor")
			{
			if(ok)
				pred <- lin.predict(X,newdata)
			else
				return(X)
			}
		else
			return(X)
		}
	if(type == "m")
		{
		if(ok)
			{
			nl <- length(X)
			pred <- vector("list",nl)
			if(inherits(X[[1]],"sm.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- sm.predict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			if(inherits(X[[1]],"mrf.gibbs"))
				return(X)
			if(inherits(X[[1]],"linear.gibbs"))
				{
				for(j in 1:nl)
					pred[[j]] <- lin.predict(X[[j]],newdata,paste(nname,":",j,sep=""))
				}
			attr(pred,"term.type") <- "m"
			attr(pred,"par.mat") <- attr(X,"par.mat")
			attr(pred,"m.term") <- attr(X,"m.term")
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
						if(inherits(X,"sm.gibbs"))
							{
							if(ok)
								pred <- sm.predict(X,newdata)
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
							if(ok)
								pred <- lin.predict(X,newdata)
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
						if(inherits(X,"sm.gibbs"))
							{
							if(ok)
								pred <- sm.predict(X,newdata)
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
							if(ok)
								pred <- lin.predict(X,newdata)
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


sm.predict <- function(x,newdata)
	{
	specs <- attr(x,"smooth.specs")
	object <- attr(specs,"smooth.construct")

	X <- Predict.matrix(object,newdata)
	coefs <- attr(x,"smooth.coef.draws.utr") 
	fitted <- qhelp(coefs,X)
	
	fout <- NULL
	for(i in 1:length(specs$term))
	fout <- cbind(fout,newdata[specs$term[i]][[1]])
	k <- ncol(fout)
	fout <- cbind(fout,fitted$mean)
	fout <- cbind(fout,fitted$q2)
	fout <- cbind(fout,fitted$q21)
	fout <- cbind(fout,fitted$median)
	fout <- cbind(fout,fitted$q11)
	fout <- cbind(fout,fitted$q1)
	fout <- cbind(fout,0)

        fout <- cbind(fout, (fout[,k+3] < 0 & fout[,k+5] < 0) * (-1) + (fout[,k+3] <= 0 & fout[,k+5] >= 0) * 0 + (fout[,k+3] > 0 & fout[,k+5] > 0) * 1)
        fout <- cbind(fout, (fout[,k+2] < 0 & fout[,k+6] < 0) * (-1) + (fout[,k+2] <= 0 & fout[,k+6] >= 0) * 0 + (fout[,k+2] > 0 & fout[,k+6] > 0) * 1)
	fout <- fout[order(fout[,1]),]

	colnames(fout) <- colnames(x)

	attr(fout,"smooth.variance") <- attr(x,"smooth.variance")
	attr(fout,"smooth.variance.draws") <- attr(x,"smooth.variance.draws")
	attr(fout,"smooth.hyper") <- attr(x,"smooth.hyper")
	attr(fout,"smooth.hyper.draws") <- attr(x,"smooth.hyper.draws")
	attr(fout,"smooth.ceffect") <- attr(x,"smooth.ceffect")
	attr(fout,"smooth.ceffect.draws") <- attr(x,"smooth.ceffect.draws")
	attr(fout,"smooth.coef.draws") <- attr(x,"smooth.coef.draws")
	attr(fout,"smooth.coef.mean") <- attr(x,"smooth.coef.mean")
	attr(fout,"smooth.coef.draws.utr") <- attr(x,"smooth.coef.draws.utr")
	attr(fout,"smooth.specs") <- attr(x,"smooth.specs")
	attr(fout,"term.type") <- "smooth"
	class(fout) <- "sm.gibbs"

	return(fout)
	}


lin.predict <- function(x,newdata)
	{
	name <- colnames(x)[1]
	if(!name%in%names(newdata))
		stop("no data to predict at")
	newdata <- newdata[name][[1]]
	coefs <- attr(x,"linear.coef")
	fout <- matrix(0,length(newdata),10)
	fout[,1] <- newdata
	fout[,2] <- newdata*coefs[1]
	fout[,3] <- newdata*coefs[3]
	fout[,4] <- newdata*coefs[4]
	fout[,5] <- newdata*coefs[5]
	fout[,6] <- newdata*coefs[6]
	fout[,7] <- newdata*coefs[7]
      	fout[,9] <- (fout[,4] < 0 & fout[,6] < 0) * (-1) + (fout[,4] <= 0 & fout[,6] >= 0) * 0 + (fout[,4] > 0 & fout[,6] > 0) * 1
      	fout[,10] <- (fout[,3] < 0 & fout[,7] < 0) * (-1) + (fout[,3] <= 0 & fout[,7] >= 0) * 0 + (fout[,3] > 0 & fout[,7] > 0) * 1
	fout <- fout[order(fout[,1]),]
	colnames(fout) <- colnames(x)
	attr(fout,"smooth.specs") <- list(names=name)
	attr(fout,"linear.coef.draws") <- attr(x,"linear.coef.draws")
	attr(fout,"linear.ceffect") <- attr(x,"linear.ceffect") 
	attr(fout,"linear.coef") <- attr(x,"linear.coef")
	attr(fout,"linear.coef.mean") <- attr(fout,"linear.coef.mean")
	class(fout) <- "linear.gibbs"
	return(fout)
	}

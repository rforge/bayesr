cprobs <- function(x, which = 1)
	{
	if(length(x$terms) < which[1] || which[1] == 0)
		stop("Argument which is specified wrong, nothing to plot!")
	X <- x$fout[[which[1]]]
	type <- attr(X,"term.type")
	here <- FALSE
	s2 <- attr(x$fitted,"variance.draws")
	nd <- length(s2)
	if(type == "smooth")
		{
		term <- attr(X,"smooth.specs")$term
		lambda <- attr(X,"smooth.hyper.draws")
		beta <- attr(X,"smooth.coef.draws")
		mbeta <- attr(X,"smooth.coef.mean")
		s <- attr(X,"s")
		prob <- spcprob(s2,lambda,beta,mbeta,s,nd)
		}
	if(type == "linear")
		{
		if(is.null(attr(X,"factorcheck")))
			attr(X,"factorcheck") <- "nonfactor"
		if(attr(X,"factorcheck") == "nonfactor")
			{
			term <- colnames(X)[1]
			lambda <- rep(1,nd)
			beta <- attr(X,"linear.coef.draws")
			mbeta <- attr(X,"linear.coef.mean")
			s <- 0
			prob <- spcprob(s2,lambda,beta,mbeta,s,nd)
			}
		else
			{
			X <- X[[1]]
			nl <- length(X)
			prob <- rep(0,nl)
			fac <- rep("",nl)
			term <- strsplit(colnames(X[[1]])[1],":")[[1]][1]
			here <- TRUE
			for(j in 1:nl)
				{
				fac[j] <- colnames(X[[j]])[1]
				lambda <- rep(1,nd)
				beta <- attr(X[[j]],"linear.coef.draws")
				mbeta <- attr(X[[j]],"linear.coef.mean")
				s <- 0
				prob[j] <- spcprob(s2,lambda,beta,mbeta,s,nd)
				}
			}
		}
	if(type == "mu")
		{
		nl <- length(X)
		prob <- rep(0,nl)
		fac <- rep("",nl)
		term <- strsplit(colnames(X[[1]])[1],":")[[1]][1]
		here <- TRUE
		if(inherits(X[[1]],"linear.gibbs"))
			{
			for(j in 1:nl)
				{
				fac[j]<- colnames(X[[j]])[1]
				lambda <- rep(1,nd)
				beta <- attr(X[[j]],"linear.coef.draws")
				mbeta <- attr(X[[j]],"linear.coef.mean")
				s <- 0
				prob[j] <- spcprob(s2,lambda,beta,mbeta,s,nd)
				}
			}
		else
			{
			for(j in 1:nl)
				{
				fac[j] <- attr(X[[j]],"smooth.specs")$term
				lambda <- attr(X[[j]],"smooth.hyper.draws")
				beta <- attr(X[[j]],"smooth.coef.draws")
				mbeta <- attr(X[[j]],"smooth.coef.mean")
				s <- attr(X[[j]],"s")
				prob[j] <- spcprob(s2,lambda,beta,mbeta,s,nd)
				}
			}
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
						stop("Argument which is specified wrong, nothing to do!")
					X <- X[[which[2]]]
					if(attr(X,"term.type") == "smooth")
						{
						term <- attr(X,"smooth.specs")$term
						lambda <- attr(X,"smooth.hyper.draws")
						beta <- attr(X,"smooth.coef.draws")
						mbeta <- attr(X,"smooth.coef.mean")
						s <- attr(X,"s")
						prob <- spcprob(s2r,lambda,beta,mbeta,s,nd)
						}
					if(attr(X,"term.type") == "linear")
						{
						term <- colnames(X)[1]
						lambda <- rep(1,nd)
						beta <- attr(X,"linear.coef.draws")
						mbeta <- attr(X,"linear.coef.mean")
						s <- 0
						prob <- spcprob(s2r,lambda,beta,mbeta,s,nd)
						}	
					if(attr(X,"term.type") == "random")
						{
						term <- colnames(X$effects)[1]
						X <- X$effects
						beta <- attr(X,"random.coefs.draws")
						mbeta <- attr(X,"random.coefs.mean")
						s <- attr(X,"s")
						lambda <- rep(1,nd)
						prob <- spcprob(s2r,lambda,beta,mbeta,s,nd)
						}
					}
				if(length(which) == 3)
					{
					if(which[2] == 0)
						stop("Argument which is specified wrong, nothing to do!")
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
						term <- attr(X,"smooth.specs")$term
						lambda <- attr(X,"smooth.hyper.draws")
						beta <- attr(X,"smooth.coef.draws")
						mbeta <- attr(X,"smooth.coef.mean")
						s <- attr(X,"s")
						prob <- spcprob(s2r,lambda,beta,mbeta,s,nd)
						}
					if(attr(X,"term.type") == "linear")
						{
						term <- colnames(X)[1]
						lambda <- rep(1,nd)
						beta <- attr(X,"linear.coef.draws")
						mbeta <- attr(X,"linear.coef.mean")
						s <- 0
						prob <- spcprob(s2r,lambda,beta,mbeta,s,nd)
						}	
					}
				}
			else
				effplot <- TRUE
			}
		if(length(which) == 1)
			effplot <- TRUE
		if(effplot)
			{
			term <- colnames(X$effects)[1]
			X <- X$effects
			beta <- attr(X,"random.coefs.draws")
			mbeta <- attr(X,"random.coefs.mean")
			s <- attr(X,"s")
			lambda <- rep(1,nd)
			prob <- spcprob(s2,lambda,beta,mbeta,s,nd)
			}
		}
	if(here)
		{
		cat("Contour probabilities of term ",term,":","\n",sep="")
		for(j in 1:nl)
			cat(fac[j],": ",prob[j],"\n",sep="")
		}
	else
		cat("Contour probability of term ",term,": ",prob,"\n",sep="")
	return(invisible(prob))

	}
qhelp <- function(draws,basis)
	{
	n <- nrow(basis)
	k <- ncol(basis)
	iter <- ncol(draws)
	np1 <- iter*0.025
	np2 <- iter*0.5
	np3 <- iter*0.975
	np4 <- iter*0.1
	np5 <- iter*0.9

	out <- .Call("qhelp",
			as.numeric(draws),
			as.numeric(basis),
			as.integer(iter),
			as.integer(n),
			as.integer(k),
			as.numeric(np1),
			as.numeric(np2),
			as.numeric(np3),
			as.numeric(np4),
			as.numeric(np5))
	}


redraw <- function(rq, draws)
	{
	out <- .Call("redraw",rq,draws)
	}


qhelp2 <- function(draws,basis,response=NULL,eta=NULL,ind=NULL)
	{
	n <- nrow(basis)
	k <- ncol(basis)
	iter <- ncol(draws)
	np1 <- iter*0.025
	np2 <- iter*0.5
	np3 <- iter*0.975
	np4 <- iter*0.1
	np5 <- iter*0.9

	if(is.null(response) && is.null(eta) && is.null(ind))
		{
		response <- rep(0,n)
		eta <- rep(0,n)
		ind <- 1:n
		}
	if(!is.null(response) && !is.null(eta) && is.null(ind))
		ind <- 1:n

	out <- .Call("qhelp2",
			as.numeric(draws),
			as.numeric(basis),
			as.integer(iter),
			as.integer(n),
			as.integer(k),
			as.numeric(np1),
			as.numeric(np2),
			as.numeric(np3),
			as.numeric(np4),
			as.numeric(np5),
			as.numeric(response),
			as.numeric(eta),
			as.integer(ind))
	colnames(out) <- c("pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")

	return(out)
	}


sfout <- function(x)
	{
	if(!is.matrix(x))
		x <- matrix(x,nrow=1)
	else
		{
		if(ncol(x)<2)
			x <- matrix(x,nrow=1)
		}
	x <- cbind(apply(x,1,mean),apply(x,1,sd),t(apply(x,1,quantile,probs=c(0.025,0.1,0.5,0.9,0.975))))
	colnames(x) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	rownames(x) <- NULL

	return(x)
	}

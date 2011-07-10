tr <- function(x,ind)
	{
	yes <- idcheck <- FALSE
	ind <- as.integer(ind)
	if(!is.matrix(x))
		{
		if(is.vector(x))
			{
			yes <- TRUE
			if(identical(x,ind))
				idcheck <- TRUE
			if(is.factor(x))
				{
				x <- as.integer(x)
				if(identical(x,ind))
					idcheck <- TRUE
				}
			x <- matrix(x,length(x),1)
			}
		else
			stop("x is not a vector or matrix!")
		}
	if(is.matrix(x) && ncol(x) < 2)
		yes <- TRUE
	n <- nrow(x)
	k <- ncol(x)

	if(idcheck)
		{
		out <- unique(ind)
		}
	else
		{
		uind <- as.integer(levels(as.factor(ind)))
		nl <- length(uind)
		umat <- matrix(0,nl,k)
		ind2 <- rep(0,n)
		out <- .Call("tr",
		as.numeric(x),
		as.numeric(umat),
		as.integer(ind),
		as.numeric(ind2),
		as.integer(uind),
		as.integer(n),
		as.integer(k),
		as.integer(nl))
		}
	if(yes)
		return(as.vector(out))
	else
		return(matrix(out,nl,k))
	}
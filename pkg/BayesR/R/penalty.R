penalty <- function(orderpenalty,basis)
	{
	D <- diag(rep(1,dim(basis)[2]))
  	for(i in 1:orderpenalty)
    		D <- diff(D)
  	K <- crossprod(D,D)

	return(K)
	}


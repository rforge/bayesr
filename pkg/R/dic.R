dic <- function(response,eta,meta,sigma2,msigma2)
	{
	n <- length(sigma2)
	D <- rep(0,n)
	for(i in 1:n) 
		D[i] <- dev(response,eta[,i],sigma2[i])
	mD <- mean(D)
	Dm <- dev(response,meta,msigma2)
	pd <- mD - Dm

	return(list(DIC=(mD + pd),pd=pd))
	}
mcenter <- function(x, z, c, kappa, kn)
	{
	mdx <- approxm(x)
	mdz <- approxm(z)

	X <- cbind(mdx,mdz)
	B <- rdist(X,kn)
	phi <- max(B)/c
	mx <- Matern(B, range=phi, scale=c, smoothness=kappa)

	return(mx)
	}
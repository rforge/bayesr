incidentm <- function(x)
	{
	cs <- colSums(x)
	mx <- rep(0,ncol(x))
	check <- cs==max(cs)
	mx[check] <- 1
	if(sum(mx) > 1)
		{
		ind <- 1:ncol(x)
		mx <- rep(0,ncol(x))
		ind <- ind[check]
		mx[ind[1]] <- 1
		}
	return(mx)
	}

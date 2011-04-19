incidentm2 <- function(x)
	{
	tx <- table(x)
	ind <- 1:length(tx)
	what <- ind[tx==max(tx)]
	out <- rep(0,length(tx))
	out[what[1]] <- 1
	return(out)
	}
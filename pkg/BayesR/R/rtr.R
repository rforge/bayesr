rtr <- function(help,coefs)
	{
	U <- 1/sum(help)*help
	return(coefs - U*sum(coefs))		
	}


rtr2 <- function(help,coefs)
	{
	out <- .Call("rtr2",
	             as.numeric(help),
	             as.numeric(coefs))	
	return(out)	
	}

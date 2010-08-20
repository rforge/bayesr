gcv <- function(object)
	{
	e <- residuals(object)
	n <- length(e)
	pd <- object$DIC$pd
	tmp <- (e/(1-pd/n))^2
	out <- sum(tmp)
	return(out/n)
	}
fitted.gibbs <- function(object, type = 1)
	{
	if(type < 2)
		return(object$fitted)
	else
		{
		f <- object$fout
		class(f) <- "gibbsfit"
		return(f)
		}
	}

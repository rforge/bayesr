cpoint <- function(a,b)
	{
	out <- .Call("cpoint",a,b)
	return(out)
	}
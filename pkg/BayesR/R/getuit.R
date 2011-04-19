getuit <- function(x, xu, n, m)
	{
	check <- is.matrix(x)
	if(check)
		ch <- 2
	else
		ch <- 1
	.Call("getuit",
	       as.numeric(x),
	       as.numeric(xu),	
	       as.integer(n),
	       as.integer(m),
		 as.integer(ch))
	}
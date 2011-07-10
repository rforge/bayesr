ctrdist <- function(map)
	{
	n <- as.integer(length(map))
	out <- .Call("cdist",
			 map,
			 n,
			 numeric(n*n))
	return(out)
	}
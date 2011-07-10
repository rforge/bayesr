centroidpos <- function(polygon)
	{
	p <- polygon
	np <- (nrow(p) - 1)

	if(is.na(p[1,1]))
		{
		p <- p[2:(np+1),]
		np <- np -1
		}
	if((p[1,1] != p[(np+1),1]) || (p[1,2] != p[(np+1),2]))
		p[(np+1),] <- p[1,]

	out <- cpos(p,np)

	return(out)
	}


centroids <- function(map)
	{
	n <- length(map)
	cp <- matrix(0,n,2)
	for(i in 1:n)
		cp[i,] <- centroidpos(map[[i]])
	return(cp)
	}

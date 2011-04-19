rw <- function(x,z)
	{
	if(missing(z))
		{
		what <- paste(match.call(expand.dots=FALSE)[2])
		return(list(x=x,name=what))
		}
	else
		{
		what1 <- paste(match.call(expand.dots=FALSE)[2])
		what2 <- paste(match.call(expand.dots=FALSE)[3])
		return(list(x=cbind(x,z),name=c(what1,what2)))
		}
	}
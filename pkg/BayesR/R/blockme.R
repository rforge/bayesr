blockme <- function(x,z)
	{
	if(is.null(x))
		return(z)
	else
		{
		nx <- nrow(x)
		nz <- nrow(z)
		cx <- ncol(x)
		cz <- ncol(z)
		m <- matrix(0,nx,cz)
		out <- cbind(x,m)
		m <- matrix(0,nz,cx)
		out <- rbind(out,cbind(m,z))
		return(out)
		}
	}
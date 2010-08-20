approxm <- function(x)
	{
	if(all(x==1))
		mx <- 1
	else
		{
		ox <- x[x!=0]
		if(all(ox==1))
			{
			ind <- 1:length(x)
			ind <- ind[x!=0]
			mx <- ind[1]
			}
		else
			{
			mx <- 1:length(x)
			mxm <- abs(x - mean(x))
			mx <- mx[mxm==min(mxm)]
			if(length(mx) > 1)
				{
				mx <- 1:length(x)
				mxm <- abs(x - mean(x))
				check <- mxm!=mean(x)
				mx <- mx[check]
				mxm <- mxm[check]
				mx <- mx[mxm==min(mxm)]
				}
			}
		}
	return(mx[1])
	}

setmfrow <- function(np,n = FALSE)
	{
	if(!n)
		{
		if(np == 1)
			par(mfrow=c(1,1),oma=c(5,5,5,5))
		if(np == 2)
			par(mfrow=c(1,2),oma=c(5,5,5,5))
		if(np == 3)
			par(mfrow=c(1,3),oma=c(5,5,5,5))
		if(np == 4)
			par(mfrow=c(2,2),oma=c(5,5,5,5))
		if(np == 5 || np == 6)
			par(mfrow=c(2,3),oma=c(5,5,5,5))
		if(np > 6)
			{
			par(mfrow=c(3,3),oma=c(5,5,5,5))
			par(ask=TRUE)
			}
		}
	else
		{
		if(np == 1)
			par(mfrow=c(1,1))
		if(np == 2)
			par(mfrow=c(1,2))
		if(np == 3)
			par(mfrow=c(1,3))
		if(np == 4)
			par(mfrow=c(2,2))
		if(np == 5 || np == 6)
			par(mfrow=c(2,3))
		if(np > 6)
			par(mfrow=c(round((np+np*0.1)/2.5),3))
		}
	}

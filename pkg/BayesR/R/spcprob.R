spcprob <- function(s2,lambda,beta,mbeta,s,n)
	{
	np <- n*0.5
	if(n%%2 == 0)
		np <- trunc((np + np + 1)/2) - 1
	else
		np <- trunc(np) 
	k <- length(s)
	.Call("cprobs",
		as.numeric(s2),
		as.numeric(lambda),
		as.numeric(beta),
		as.numeric(mbeta),
		as.numeric(s),
		as.integer(np),
		as.integer(n),
		as.integer(k))
	}
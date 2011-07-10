rtvnorm <- function(n, mean = 0, sd = 1, low = -Inf, high = Inf)
	{
	if(length(mean) != n)
		mean <- rep(mean,n)
	if(length(sd) != n)
		sd <- rep(sd[1],n)
	if(length(low) != n)
		low <- rep(low,n)
	if(length(high) != n)
		high <- rep(high,n)
	.Call("rtvnorm",
		as.integer(n),
		as.numeric(mean),
		as.numeric(sd),
		as.numeric(low),
		as.numeric(high))
	}


rtvnorm2 <- function(n, mean = 0, sd = 1, low = -Inf, high = Inf)
	{
	if(length(mean) != n)
		mean <- rep(mean,n)
	if(length(sd) != n)
		sd <- rep(sd[1],n)
	if(length(low) != n)
		low <- rep(low,n)
	if(length(high) != n)
		high <- rep(high,n)
	.Call("rtvnorm2",
		as.integer(n),
		as.numeric(mean),
		as.numeric(sd),
		as.numeric(low),
		as.numeric(high))
	}


rgamma2 <- function(n, a, b)
	{
	.Call("rgamma2",
		as.integer(n),
		as.numeric(a),
		as.numeric(b))
	}

bspline <- function(x, degree = 3, knots = 20)
	{
	minx <- min(x)-0.001
	maxx <- max(x)+0.001
	step <- (maxx-minx)/(knots-1)
	k <- seq(minx-degree*step,maxx+degree*step,by=step)

	basis <- spline.des(knots=k,x,ord=(degree+1),outer.ok=TRUE)$design

	return(basis)
	}


dev <- function(response,eta,sigma2)
	{
	# D <- -2*sum(dnorm(response, mean = eta, sd = sqrt(sigma2), log = TRUE))
	D <- sum(log(2*pi*sigma2)+1/sigma2*(response-eta)^2)

	return(D)
	}
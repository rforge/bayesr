print.gibbs <- function(x,...)
	{
	object <- x
	cat("Call:\n")
	print(object$call)
	cat("-------------------------------------------------------------------\n")
	cat("DIC =", object$DIC$DIC, " pd =", object$DIC$pd, " n =", object$N, "\n")
	cat("Iterations =", object$iterations, " burnin =", object$burnin, " thinning =", object$thinning, "\n")
	}
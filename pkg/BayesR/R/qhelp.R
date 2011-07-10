qhelp <- function(draws,basis)
	{
	n <- nrow(basis)
	k <- ncol(basis)
	iter <- ncol(draws)
	np1 <- iter*0.025
	np2 <- iter*0.5
	np3 <- iter*0.975
	np4 <- iter*0.1
	np5 <- iter*0.9

	out <- .Call("qhelp",
			as.numeric(draws),
			as.numeric(basis),
			as.integer(iter),
			as.integer(n),
			as.integer(k),
			as.numeric(np1),
			as.numeric(np2),
			as.numeric(np3),
			as.numeric(np4),
			as.numeric(np5))
	}


redraw <- function(rq, draws)
	{
	out <- .Call("redraw",rq,draws)
	}

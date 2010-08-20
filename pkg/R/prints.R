print.sm.gibbs <- function(x,...)
	{
	out <- x[,1:ncol(x)]
	rnam <- rownames(x)
	cnam <- colnames(x)
	rownames(out) <- rnam
	colnames(out) <- cnam
	print(out,...)
	}

print.mrf.gibbs <- function(x,...)
	{
	out <- x[,1:ncol(x)]
	rnam <- rownames(x)
	cnam <- colnames(x)
	rownames(out) <- rnam
	colnames(out) <- cnam
	print(out,...)
	}

print.random.gibbs <- function(x,...)
	{
	out <- x[,1:ncol(x)]
	rnam <- rownames(x)
	cnam <- colnames(x)
	rownames(out) <- rnam
	colnames(out) <- cnam
	print(out,...)
	}
	
print.linear.gibbs <- function(x,...)
	{
	out <- x[,1:ncol(x)]
	rnam <- rownames(x)
	cnam <- colnames(x)
	rownames(out) <- rnam
	colnames(out) <- cnam
	print(out,...)
	}

print.gibbsfit <- function(x,...)
	{
	n <- length(x)
	cat(paste("list() with",n,"fitted gibbs() term object(s), to get results for the first term just type fitted(gibbsobject)[[1]], also plot this term with plot(fitted(gibbsobject)[[1]])\n"))
	}

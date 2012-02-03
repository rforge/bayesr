cpos <- function(p,np)
{
	o <- c(0,0)
	.Call("cpos",
		as.numeric(p),
		as.integer(np),
		as.numeric(o))
}

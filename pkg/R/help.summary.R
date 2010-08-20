help.summary <- function(tmp,r.smoothhyp,r.coeflin,s2r)
	{
	if(!is.null(tmp$terms))
		for(j in 1:length(tmp$terms))
			if(attr(tmp$terms[[j]],"term.type") == "random")
				{
				r.smoothhyp <- rbind(r.smoothhyp,attr(tmp$terms[[j]],"smooth.mat"))
				r.coeflin <- rbind(r.coeflin,attr(tmp$terms[[j]],"lin.mat"))
				s2r <- rbind(s2r,attr(tmp$terms[[j]]$effects,"random.variance"))
				}
	}
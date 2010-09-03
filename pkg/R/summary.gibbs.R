summary.gibbs <- function(object, digits = 4,...)
	{
	K <- object$K
	M <- object$M
	R <- object$R
	smoothhyp <- smoothmuhyp <- smoothrahyp <- coeflinmu <- ra <- NULL
	coeflin <- round(attr(object$fout,"lin.mat"),digits)

	if(K > 0)
		smoothhyp <- round(attr(object$fout,"smooth.mat"),digits)
	if(M > 0)
		{
		if(!is.null(attr(object$fout,"smooth.mu.mat")))
			smoothmuhyp <- round(attr(object$fout,"smooth.mu.mat"),digits)
		if(!is.null(attr(object$fout,"lin.mu.mat")))
			coeflinmu <- round(attr(object$fout,"lin.mu.mat"),digits)
		}
	if(R > 0)
		{
		ranpos <- attr(object$fout,"ranpos")
		r.coeflin <- r.smoothhyp <- s2ro <- lr <- sr <- NULL
		for(k in 1:R)
			{
			tmp <- object$fout[[ranpos[k]]]
			s2r <- attr(tmp,"var.mat")
			r.smoothhyp <- attr(tmp,"smooth.mat")
			r.coeflin <- attr(tmp,"lin.mat")
			type <- attr(tmp,"term.type")
			
			if(!is.null(tmp$terms))
				for(j in 1:length(tmp$terms))
					if(attr(tmp$terms[[j]],"term.type") == "random")
						{
						r.smoothhyp <- rbind(r.smoothhyp,attr(tmp$terms[[j]],"smooth.mat"))
						r.coeflin <- rbind(r.coeflin,attr(tmp$terms[[j]],"lin.mat"))
						}
						
			if(!is.null(r.smoothhyp))
				r.smoothhyp <- unique(round(r.smoothhyp,digits))
			if(!is.null(r.coeflin))
				r.coeflin <- unique(round(r.coeflin,digits))
			s2r <- round(s2r,digits)
			s2ro <- rbind(s2ro,s2r)
			lr <- rbind(lr,r.coeflin)
			sr <- rbind(sr,r.smoothhyp)
			}
		ra <- list(smoothhyp=sr,coeflin=lr,s2r=s2ro)
		}
	S2 <- round(attr(fitted(object),"variance"),digits)
	object$DIC$DIC <- round(object$DIC$DIC,digits)
	object$DIC$pd <- round(object$DIC$pd,digits)
	method <- "MCMC"
	if(!is.null(attr(object,"method")))
		method <- attr(object,"method")
	res <- list(call=object$call,K=K,M=M,R=R,coeflin=coeflin,smoothhyp=smoothhyp,
                    smoothmuhyp=smoothmuhyp,smoothrahyp=smoothrahyp,S2=S2,method=method,
                    DIC=object$DIC,burnin=object$burnin,thinning=object$thinning,iter=object$iterations,
                    coeflinmu=coeflinmu,N=object$N,ra=ra,digits=digits,samptime=object$samptime[3])
	class(res) <- "summary.gibbs"
	res
	}

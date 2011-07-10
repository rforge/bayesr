trans.design <- function(basis,K)
	{
	nc <- ncol(basis)

	if((missing(K) || is.null(K) || K == 0) && (nc == 1 || is.null(nc)))
		{
		XX <- crossprod(basis,basis)
		r <- chol(XX)
		rinv <- solve(r)
		rr <- crossprod(rinv,rinv)
		SVD <- svd(rr)
		Q <- SVD$u
		tildeX <- basis%*%rinv%*%Q
		RQ <- rinv%*%Q

		return(list(tildeZ=tildeX,s=0,RQ=RQ))
		}
	if(nc > 1 || !is.null(nc))
		{
		if(missing(K) || is.null(K) || K == 0)
			{
			bb <- crossprod(basis,basis)
			r <- try(chol(bb), silent=TRUE)
			if(inherits(r, "try-error"))
				{
				bb <- bb + diag(diag(bb) + 0.00001)
				r <- chol(bb)
				}
			rinv <- solve(r)
			u <- crossprod(rinv,rinv)
			SVD <- svd(u)
			Q <- SVD$u
			tildeZ <- basis%*%rinv%*%Q
			s <- rep(0,ncol(basis))
			RQ <- rinv%*%Q

			return(list(tildeZ=tildeZ,s=s,RQ=RQ))
			}
		else
			{
			bb <- crossprod(basis,basis)
			r <- try(chol(bb), silent=TRUE)
			if(inherits(r, "try-error"))
				{
				bb <- bb + diag(diag(bb) + 0.00001)
				r <- chol(bb)
				}
			rinv <- solve(r)

			u <- crossprod(rinv,K)%*%rinv
			SVD <- svd(u)
			Q <- SVD$u
			tildeZ <- basis%*%rinv%*%Q
			s <- SVD$d
			RQ <- rinv%*%Q

			return(list(tildeZ=tildeZ,s=s,RQ=RQ))
			}
		}
	}
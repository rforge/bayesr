gibbsupdate <- function(Z,tZ,sigma2,tau2,O,response,eta,dP,
		        ce,sce,cmcheck,cK,new,cme,hypa,hypb)
	{
	if(is.null(Z$lambda))
		lambd <- (sigma2/tau2)
	else
		lambd <- Z$lambda
			
	help <- 1/(1+lambd*Z$s)
		
	mutmp <- help*crossprod(O,(response - eta))
	draws <- rnorm(dP,mutmp,(help*sigma2)^0.5)

	# Evaluate center effect
	ce <- ce - sce
	if(cmcheck$a)
		sce <- cK%*%(cmcheck$b%*%draws)
	else
		sce <- cK%*%draws
	ce <- ce + sce

	ftmp <- crossprod(tZ,draws)

	# Centering
	if(cme == 1)
		{
		# draws <- rtr(Sigma2,draws)
		ftmp <- ftmp - mean(ftmp)
		# ftmp <- mycenter(ftmp,Nnew[k])
		}
	fs <- ftmp[new]

	sums <- sum(draws^2*Z$s)
	hypb <- hypb + sums/2
	tau2 <- 1/(rgamma(1,shape=hypa,rate=hypb)) 

	return(list(mutmp=mutmp,draws=draws,sce=sce,ce=ce,tau2=tau2,fs=fs))
	}

findModel <- function(formula, data, weights, family = "normal", iter = 1200, burnin = 200, thinning = 1, 
			    hyperasigma = 1e-04, hyperbsigma = 1e-04, hyperatau = 1e-04, hyperbtau = 1e-04, 
			    trace = FALSE, dots = FALSE, epsilon = 0.1, ...)
	{
	model <- gibbs(formula,data,weights,family,iter,burnin,thinning, 
			   hyperasigma,hyperbsigma,hyperatau,hyperbtau, 
			   trace,dots)
	DIC <- dicone <- model$DIC$DIC
	diff <- 10
	form <- formone <- model$formula
	terms <- model$terms

	while(diff > epsilon)
		{
		lt <- length(terms)
		if(lt == 1)
			return(list(formula=form,DIC=DIC))
		ind <- 1:(lt+1)
		indo <- 1:lt
		forms <- vector("list",length=(lt+1))
		dic <- rep(0,(lt+1))
		dic[(lt+1)] <- dicone
		forms[[lt+1]] <- formone
		for(i in 1:lt)
			{
			tmp <- terms[indo!=i]
			ltmp <- length(tmp)
			this <- rep("",ltmp)
			if(ltmp > 1)
				for(j in 1:ltmp)
					{
					if(j < ltmp)
						this[j] <- paste(tmp[j],"+",sep="")
					else
						this[j] <- tmp[j]
					}
			else
				this <- tmp
			this <- paste(this,collapse="")
			chf <- as.character(form)
			chf <- paste(chf[2],chf[1],chf[3],sep="")
			txt <- paste("form <- update(",chf,",~",this,")",sep="")
			eval(parse(text = txt))
			model <- gibbs(form,data,weights,family,iter,burnin,thinning, 
			   		   hyperasigma,hyperbsigma,hyperatau,hyperbtau, 
			   		   trace,dots)
			forms[[i]] <- form
			dic[i] <- model$DIC$DIC
			}
		take <- ind[dic==min(dic)]
		diff <- DIC - dic[take]
		DIC <- dic[take]
		form <- forms[[take]]
		}

	return(list(formula=form,DIC=DIC))
	}
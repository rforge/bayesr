gibbs2 <- function(formula,family="gaussian",data=NULL,
		   method="MCMC",iter=1200,burnin=200,thinning=1,eps=0.0001)
	{
	res <- list()
	tf <- terms(formula, specials = c("s","te","m","r"))
	terms <- attr(tf, "term.labels")
	res$terms <- terms
	lt <- length(terms)
	gp <- interpret.gam(formula)
	if(is.null(data))
		data <- parent.frame()
	dat <- model.matrix(gp$fake.formula,data)
	if(attr(tf,"intercept")==1)
		{
		vars <- colnames(dat)
		vars <- vars[2:length(vars)]
		if(ncol(dat)<3)
			dat <- matrix(dat[,2:ncol(dat)],ncol=1)
		}
	else
		vars <- colnames(dat)
	resp <- eval(parse(text=gp$response),envir=data)
	dat <- cbind(resp,dat)
	colnames(dat) <- c(gp$response,vars)
	
	path<-paste(.libPaths()[1],"/BayesR/bayesxtemp",sep="")
	dir.create(path)
	patho <- getwd()
	setwd(path)
	write.table(dat,"data.raw",col.names=TRUE,row.names=FALSE,quote=FALSE)
	addterm <- NULL
	T <- list()
	if(lt>0)
		{
		ind <- c(1:lt)
		sm <- attr(tf,"specials")$s-1
		te <- attr(tf,"specials")$te-1
		lin <- ind[!ind %in% c(sm,te)]
		smooth <- ind[!ind %in% lin]
		L <- length(lin)
		K <- length(sm)+length(te)
		if(K>0)
			{
			for(k in 1:K)
				{
				smtmp <- eval(parse(text = terms[smooth[k]]))
				T[[smooth[k]]] <- bayesx.setup(smtmp)
				addterm <- c(addterm,T[[smooth[k]]]$term)
				}
			}
		if(L>0)
			for(k in 1:L)
				{
				addterm <- c(addterm,terms[lin[k]])
				T[[lin[k]]] <- list(term=terms[lin[k]],dim=1,spec="lin.spec")
				}
		}
	if(method=="MCMC")
		fullformula <- paste("b.regress",gp$response,"=",addterm,", family =",family,"iterations =",iter,"burnin =",burnin,"step =",thinning,"predict using d")
	else
		fullformula <- paste("b.regress",gp$response,"=",addterm,", family =",family,"eps =",eps,"predict using d")
	# now its time for bayesx
	outfile <- "bayesx.prg"
	cat("% usefile bayesx.prg\n",file=outfile,append=FALSE)
	cat("dataset d\n",file=outfile,append=TRUE)
	cat(paste("d.infile using ",path,"/data.raw\n",sep=""),file=outfile,append=TRUE)

	if(method=="MCMC")
		cat("bayesreg b\n",file=outfile,append=TRUE)
	else
		cat("remlreg b\n",file=outfile,append=TRUE)
	cat(paste("b.outfile = ",path,"/estim\n",sep=""),file=outfile,append=TRUE)

	cat(fullformula,"\n",file=outfile,append=TRUE)
	cat(paste("b.outfile = ",path,"/samples\n",sep=""),file=outfile,append=TRUE)
	cat("b.getsample\n",file=outfile,append=TRUE)

	# estimate with BayesX
	cat("Starting:\n")
	ptm <- proc.time()
	system(paste("BayesX ",path,"/bayesx.prg",sep=""))
	now <- proc.time()
	samptime <- now - ptm
	res$samptime <- samptime
	cat("Total sampling time was:",samptime[3],"sec\n")

	res$fout <- get.bayesx.output(T,path,resp)
	res$iterations <- iter
	res$burnin <- burnin
	res$thinning <- thinning
	# system(paste("rm -Rf",path))
	setwd(patho)
	class(res) <- "gibbs"
	return(res)
	}


bayesx.setup <- function(object)
	{
	if(is.na(object$p.order))
		object$p.order <- 0
	if(object$bs.dim<0)
		object$bs.dim <- 8
	if(is.null(object$m))
		object$m <- 2
	else
		{
		if(length(object$m)>1)
			{
			if(object$m[2] > 2)
				{
				warning("Order of the difference penalty not supported by BayesX, set to 2!")
				object$m <- 2
				}
			}
		else
			{
			if(object$m > 2)
				{
				warning("Order of the difference penalty not supported by BayesX, set to 2!")
				object$m <- 2
				}
			}
		}
	if(class(object)=="tp.smooth.spec")
		stop("BayesX does not support smooths of class tp.smooth.spec!") 
	if(class(object)=="ps.smooth.spec" && object$dim<2)
		{
		nrknots <- object$bs.dim - object$p.order[1]
		termo <- object$term
		term <- paste(termo,"(psplinerw",object$m[1],",nrknots=",nrknots,")",sep="")
		}
	return(list(term=term,dim=object$dim,spec=class(object),to=termo))
	}


get.bayesx.output <- function(terms,path,response)
	{
	setwd(path)
	fout <- vector("list",length=length(terms))
	smooth.mat <- NULL
	eta <- read.table("estim_predictmean.raw",header=TRUE)$linpred
	for(k in 1:length(terms))
		{
		te <- terms[[k]]
		if(te$dim<2)
			{
			if(te$spec=="ps.smooth.spec")
				{
				what <- paste("estim_f_",te$to,"_pspline.res",sep="")
				fout[[k]] <- as.matrix(read.table(what,header=TRUE))
				name <- colnames(fout[[k]])[2]
				fout[[k]] <- fout[[k]][,2:ncol(fout[[k]])]
				resids <- response - eta + fout[[k]][,2]
				fout[[k]] <- cbind(fout[[k]][,1:(ncol(fout[[k]])-2)],  resids,  fout[[k]][,(ncol(fout[[k]])-1):ncol(fout[[k]])])
				fout[[k]] <- fout[[k]][order(fout[[k]][,1]),]
				colnames(fout[[k]]) <- c(name,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat80","pcat90")
				attr(fout[[k]],"smooth.specs") <- list(dim=1)
				attr(fout[[k]],"term.type") <- "smooth"
				class(fout[[k]]) <- "sm.gibbs"
				what <- paste("estim_f_",te$to,"_pspline_lambda.res",sep="")
				tmp <- as.matrix(read.table(what,header=TRUE))
				smooth.mat <- rbind(smooth.mat,tmp)
				}
			}
		else
			{
			fout <- 1
			}
		}
	attr(fout,"smooth.mat") <- smooth.mat
	lin.mat <- read.table("estim_FixedEffects1.res",header=TRUE)
	linnames <- names(lin.mat)
	
	lin.mat <- lin.mat[2:ncol(lin.mat)]
	attr(fout,"lin.mat") <- lin.mat
	attr(fout,"variance") <- as.matrix(read.table("estim_scale.res",header=TRUE))
	return(fout)
	}

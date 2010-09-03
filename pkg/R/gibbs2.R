gibbs2 <- function(formula,family="gaussian",data=NULL,
		   method="MCMC",iter=1200,burnin=200,thinning=1,eps=0.0001,
		   path=NULL,rm=TRUE)
	{
	res <- list()
	res$call <- match.call()
	res$formula <- formula <- as.formula(formula)
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
		else
			dat <- dat[,2:ncol(dat)]
		}
	else
		vars <- colnames(dat)
	resp <- eval(parse(text=gp$response),envir=data)
	dat <- cbind(resp,dat)
	colnames(dat) <- c(gp$response,vars)
	if(is.null(path))	
		path <- paste(.libPaths()[1],"/BayesR/bayesxtemp",sep="")
	dir.create(path)
	patho <- getwd()
	setwd(path)
	write.table(dat,"data.raw",col.names=TRUE,row.names=FALSE,quote=FALSE)
	res$N <- nrow(dat)
	addterm <- NULL
	T <- list()
	res$M <- res$R <- K <- L <- M <- R <- 0
	if(lt>0)
		{
		ind <- c(1:lt)
		sm <- attr(tf,"specials")$s-1
		te <- attr(tf,"specials")$te-1
		lin <- ind[!ind %in% c(sm,te)]
		smooth <- ind[!ind %in% lin]
		res$L <- L <- length(lin)
		res$K <- K <- length(sm)+length(te)
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
			{
			datlin <- NULL
			for(k in 1:L)
				{
				addterm <- c(addterm,terms[lin[k]])
				T[[lin[k]]] <- list(term=terms[lin[k]],dim=1,spec="lin.spec")
				}
			}
		}
	addterm <- paste(addterm,collapse=" + ")
	if(method=="MCMC")
		fullformula <- paste("b.regress",gp$response,"=",addterm,", family =",family,"iterations =",iter,"burnin =",burnin,"step =",thinning,"predict using d")
	else
		fullformula <- paste("b.regress",gp$response,"=",addterm,", family =",family,"eps =",eps,"using d")
	# now its time for bayesx
	outfile <- "bayesx.prg"
	cat(paste("% usefile ",path,"/bayesx.prg\n",sep=""),file=outfile,append=FALSE)
	cat("dataset d\n",file=outfile,append=TRUE)
	cat(paste("d.infile using ",path,"/data.raw\n",sep=""),file=outfile,append=TRUE)

	if(method=="MCMC")
		cat("bayesreg b\n",file=outfile,append=TRUE)
	else
		cat("remlreg b\n",file=outfile,append=TRUE)
	cat(paste("b.outfile = ",path,"/estim\n",sep=""),file=outfile,append=TRUE)

	cat(fullformula,"\n",file=outfile,append=TRUE)
	if(method=="MCMC")
		{
		cat(paste("b.outfile = ",path,"/samples\n",sep=""),file=outfile,append=TRUE)
		cat("b.getsample\n",file=outfile,append=TRUE)
		}
	# estimate with BayesX
	cat("Starting:\n")
	ptm <- proc.time()
	system(paste("BayesX ",path,"/bayesx.prg",sep=""))
	now <- proc.time()
	samptime <- now - ptm
	res$samptime <- samptime
	cat("Total sampling time was:",samptime[3],"sec\n")

	res$fitted <- list()
	res$fout <- get.bayesx.output(T,path,resp,as.data.frame(dat),method)
	attr(res$fitted,"variance") <- attr(res$fout,"variance")
	res$DIC <- attr(res$fout,"DIC")
	attr(res$fout,"DIC") <- NULL
	res$iterations <- iter
	res$burnin <- burnin
	res$thinning <- thinning
	if(rm)
		system(paste("rm -Rf",path))
	setwd(patho)
	attr(res,"method") <- method
	class(res) <- "gibbs"

	return(res)
	}


bayesx.setup <- function(object)
	{
	if(class(object)=="tensor.smooth.spec")
		{
		termo <- object$term
		object$p.order <- object$margin[[1]]$p.order
		object$bs.dim <- (object$margin[[1]]$bs.dim*2) + 1
		if(is.na(object$p.order[1]))
			object$p.order <- c(3,2)
		if(object$bs.dim<0)
			object$bs.dim <- 8
		else
			{
			if(length(object$p.order)>1)
				{
				if(object$p.order[2] > 2)
					{
					warning("Order of the difference penalty not supported by BayesX, set to 2!")
					object$p.order <- c(object$p.order[1],2)
					}
				}
			}
		nrknots <- object$bs.dim - object$p.order[1]
		term <- paste(termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
		}
	else
		{
		if(is.na(object$p.order[1]))
			object$p.order <- c(3,2)
		if(object$bs.dim<0)
			object$bs.dim <- 8
		else
			{
			if(length(object$p.order)>1)
				{
				if(object$p.order[2] > 2)
					{
					warning("Order of the difference penalty not supported by BayesX, set to 2!")
					object$p.order <- c(object$p.order[1],2)
					}
				}
			}
		if(class(object)=="tp.smooth.spec")
			stop("BayesX does not support smooths of class tp.smooth.spec!") 
		if(class(object)=="ps.smooth.spec" && object$dim<2)
			{
			nrknots <- object$bs.dim - object$p.order[1]
			termo <- object$term
			term <- paste(termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
			}
		}

	return(list(term=term,dim=object$dim,spec=class(object),to=termo,label=object$label))
	}


get.bayesx.output <- function(terms,path,response,data,method)
	{
	setwd(path)
	fout <- vector("list",length=length(terms))
	smooth.mat <- NULL
	if(method=="REML")
		{
		eta <- read.table("estim_predict.raw",header=TRUE)$eta
		lin.mat <- as.matrix(read.table("estim_FixedEffects.res",header=TRUE))
		linnames <- as.character(lin.mat[,2])
		lin.mat <- lin.mat[,3:(ncol(lin.mat)-3)]
		}
	else
		{
		eta <- read.table("estim_predictmean.raw",header=TRUE)$linpred
		lin.mat <- as.matrix(read.table("estim_FixedEffects1.res",header=TRUE))
		linnames <- as.character(lin.mat[,2])
		lin.mat <- lin.mat[,3:(ncol(lin.mat)-2)]
		}

	mode(lin.mat) <- "numeric"
	if(!is.matrix(lin.mat))
		lin.mat <- matrix(lin.mat,nrow=1)
	if(linnames[1]=="const")
		linnames[1] <- "(Intercept)"
	rownames(lin.mat) <- linnames
	if(method=="MCMC")
		colnames(lin.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	else
		colnames(lin.mat) <- c("pmode","ci95lower","ci80lower","std","ci80upper","ci95upper")
	for(k in 1:length(terms))
		{
		te <- terms[[k]]
		if(te$dim<2)
			{
			if(te$spec=="ps.smooth.spec")
				{
				dat <- as.numeric(as.matrix(data[te$to]))
				ind <- get.unique(dat)$ind[order(dat)]
				what <- paste("estim_f_",te$to,"_pspline.res",sep="")
				fout[[k]] <- as.matrix(read.table(what,header=TRUE))
				name <- colnames(fout[[k]])[2]
				fout[[k]] <- fout[[k]][,2:ncol(fout[[k]])]
				fout[[k]] <- fout[[k]][order(fout[[k]][,1]),]
				fout[[k]] <- fout[[k]][ind,]
				resids <- response[order(dat)] - eta[order(dat)] + fout[[k]][,2]
				fout[[k]] <- cbind(fout[[k]][,1:(ncol(fout[[k]])-2)],  resids,  fout[[k]][,(ncol(fout[[k]])-1):ncol(fout[[k]])])
				if(method=="MCMC")
					colnames(fout[[k]]) <- c(name,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
				else
					colnames(fout[[k]]) <- c(name,"pmode","ci95lower","ci80lower","std","ci80upper","ci95upper","partial.resid","pcat95","pcat80")
				attr(fout[[k]],"smooth.specs") <- list(dim=1)
				attr(fout[[k]],"term.type") <- "smooth"
				class(fout[[k]]) <- "sm.gibbs"
				what <- paste("estim_f_",te$to,"_pspline_var.res",sep="")
				tmp <- as.matrix(read.table(what,header=TRUE))
				if(method=="MCMC")
					tmp <- matrix(tmp[,1:(ncol(tmp)-2)],nrow=1)
				rownames(tmp) <- te$label
				smooth.mat <- rbind(smooth.mat,tmp)
				}
			if(te$spec=="lin.spec")
				{
				x <- unlist(data[te$term])
				coefs <- lin.mat[linnames==te$term,]
				if(method=="MCMC")
					{
					fpqu2p5 <- x*coefs[3]
					fpqu10 <- x*coefs[4]
					fpqu90 <- x*coefs[6]
					fpqu97p5 <- x*coefs[7]
					fpmean <- x*coefs[1]
					fpmed <- x*coefs[5]
					}
				else
					{
					fpqu2p5 <- x*coefs[2]
					fpqu10 <- x*coefs[3]
					fpqu90 <- x*coefs[5]
					fpqu97p5 <- x*coefs[6]
					fpmean <- x*coefs[1]
					fpmed <- coefs[4]
					}
				e <- response - (eta - fpmean)
				pcat80 <- (fpqu10 < 0 & fpqu90 < 0)*(-1) + (fpqu10 <= 0 & fpqu90 >= 0)*0 + (fpqu10 > 0 & fpqu90 > 0)*1
				pcat95 <- (fpqu2p5 < 0 & fpqu97p5 < 0)*(-1) + (fpqu2p5 <= 0 & fpqu97p5 >= 0)*0 + (fpqu2p5 > 0 & fpqu97p5 > 0)*1
				fout[[k]] <- cbind(x,fpmean,fpqu2p5,fpqu10,fpmed,fpqu90,fpqu97p5,e,pcat95,pcat80)
				fout[[k]] <- fout[[k]][order(fout[[k]][,1]),]
				if(method=="MCMC")
					colnames(fout[[k]]) <- c(te$term,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
				else
					colnames(fout[[k]]) <- c(te$term,"pmode","ci95lower","ci80lower","std","ci80upper","ci95upper","partial.resid","pcat95","pcat80")
				rownames(fout[[k]]) <- NULL
				attr(fout[[k]],"linear.coef") <- coefs
				# attr(fout,"linear.ceffect") <- tmpceffect[k]
				# attr(fout,"linear.coef.draws") <- draws[k+plus,]
				# attr(fout,"linear.coef.draws.utr") <- redraws[k+plus,]
				# attr(fout,"linear.coef.mean") <- mdraws[k+plus,]
				attr(fout[[k]],"term.type") <- "linear"
				class(fout[[k]]) <- "linear.gibbs"
				}
			}
		else
			{
			if(te$spec=="tensor.smooth.spec")
				{
				dat <- data[te$to]
				dat <- as.matrix(dat,ncol=ncol(dat))
				ind <- get.unique(dat)$ind[order(dat[,1])]
				what <- paste("estim_f_",te$to[1],"_",te$to[2],"_pspline.res",sep="")
				fout[[k]] <- as.matrix(read.table(what,header=TRUE))
				fout[[k]] <- fout[[k]][order(fout[[k]][,1]),]
				fout[[k]] <- fout[[k]][ind,]
				names <- colnames(fout[[k]])[2:3]
				fout[[k]] <- fout[[k]][,2:ncol(fout[[k]])]
				resids <- response - eta + fout[[k]][,3]
				fout[[k]] <- cbind(fout[[k]][,1:(ncol(fout[[k]])-2)],  resids,  fout[[k]][,(ncol(fout[[k]])-1):ncol(fout[[k]])])
				if(method=="MCMC")
					colnames(fout[[k]]) <- c(names,"pmean","pqu2p5","pqu10","pmed","pqu90","pqu97p5","partial.resid","pcat95","pcat80")
				else
					colnames(fout[[k]]) <- c(names,"pmode","ci95lower","ci80lower","std","ci80upper","ci95upper","partial.resid","pcat95","pcat80")
				attr(fout[[k]],"smooth.specs") <- list(dim=2)
				attr(fout[[k]],"term.type") <- "smooth"
				class(fout[[k]]) <- "sm.gibbs"
				what <- paste("estim_f_",te$to[1],"_",te$to[2],"_pspline_var.res",sep="")
				tmp <- as.matrix(read.table(what,header=TRUE))
				if(method=="MCMC")
					tmp <- matrix(tmp[,1:(ncol(tmp)-2)],nrow=1)
				rownames(tmp) <- te$label
				smooth.mat <- rbind(smooth.mat,tmp)
				}
			}
		}
	if(method=="MCMC")
		colnames(smooth.mat) <- c("pmean","psd","pqu2p5","pqu10","pmed","pqu90","pqu97p5")
	else
		colnames(smooth.mat) <- c("variance","smoothpar","df","stopped") 
	attr(fout,"smooth.mat") <- smooth.mat
	attr(fout,"lin.mat") <- lin.mat
	attr(fout,"variance") <- as.matrix(read.table("estim_scale.res",header=TRUE))
	rownames(attr(fout,"variance")) <- "Sigma2"
	if(method=="MCMC")
		{
		Dm <- dev(response,eta,attr(fout,"variance")[,1])
		diceta <- read.table("estim_deviance_sample.raw",header=TRUE)$unstandardized_deviance
		diceta <- sum(diceta)/length(diceta)
		pd <- diceta - Dm
		attr(fout,"DIC") <- list(DIC = (diceta + pd), pd = pd)
		}
	else
		{
		mf <- read.table("estim_modelfit.raw",header=TRUE)
		attr(fout,"DIC") <- list(DIC = mf$bic, pd = mf$df)
		}

	return(fout)
	}

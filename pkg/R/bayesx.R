write.bayesx.input <- function(formula, family="gaussian", data=NULL, method="MCMC",
                               iter=1200, burnin=200, thinning=1, file=NULL)
	{
	specs <- list()
	specs$family <- family
	specs$method <- method
	specs$iter <- iter
	specs$burnin <- burnin
	specs$thinning <- thinning
	data.file=NULL
	if(is.null(data))
		data <- parent.frame()
	else
		{
		if(is.character(data))
			{
			data.file <- data
			data <- read.table(data,header=TRUE)
			}
		if(!is.data.frame(data))
			data <- as.data.frame(data)
		}
	if(is.null(file))
		file <- paste(tempdir(),"/bayesx",sep="")
	dir.create(file,showWarnings=FALSE)
	tf <- terms(formula, specials = c("s","te","m","r"))
	terms <- attr(tf, "term.labels")
	lt <- length(terms)
	gp <- interpret.gam(formula)
	fake.formula <- gp$fake.formula
	fake.labels <- attr(terms(fake.formula),"term.labels")
	racheck <- FALSE
	byra <- rep("NA",length(fake.labels))
	for(k in 1:length(fake.labels))
		{
		splits <- strsplit(fake.labels[k],""," ")[[1]][1:2]
		if(length(splits)>1)
			{
			splits <- paste(splits[1:2],collapse="")
			if(splits=="r(")
				{
				terms[terms==fake.labels[k]] <- "99999999111"
				raterm <- eval(parse(text=fake.labels[k]),envir=data)
				if(!is.null(raterm$byterms))
					byra[k] <- raterm$byterms
				fake.labels[k] <- raterm$facname
				terms[terms=="99999999111"] <- fake.labels[k]	
				racheck <- TRUE
				}
			}
		}
	if(racheck)
		{
		for(k in 1:length(byra))
			if(byra[k]!="NA")
				fake.labels <- c(fake.labels,byra[k])
		fake.formula <- as.character(fake.formula)
		fake.formula <- as.formula(paste(fake.formula[2],"~",paste("+",fake.labels,collapse=""),collapse=""))
		}
	specs$file <- file
	fileo <- getwd()
	setwd(file)
	if(is.null(data.file))
		{
		dat <- model.matrix(fake.formula,data)
		if(attr(tf,"intercept")==1 || racheck)
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
		write.table(dat,"bayesx.data.raw",col.names=TRUE,row.names=FALSE,quote=FALSE)
		specs$dat.file <- file
		data.file.check <- FALSE
		}
	else
		{
		data.file.check <- TRUE
		specs$dat.file <- data.file
		}
	outfile <- "bayesx.input.prg"
	cat(paste("% usefile ",file,"/bayesx.input.prg\n",sep=""),file=outfile,append=FALSE)
	specs$N <- nrow(dat)
	addterm <- NULL
	T <- list()
	specs$L <- specs$K <- specs$M <- specs$R <- 0
	L <- K <- M <- R <- 0
	if(lt>0)
		{
		ind <- c(1:lt)
		sm <- attr(tf,"specials")$s-1
		te <- attr(tf,"specials")$te-1
		ran <- attr(tf,"specials")$r-1
		lin <- ind[!ind %in% c(sm,te,ran)]
		smooth <- ind[!ind %in% c(lin,ran)]
		specs$L <- L <- length(lin)
		specs$K <- K <- length(sm)+length(te)
		specs$R <- R <- length(ran)
		if(K>0)
			{
			for(k in 1:K)
				{
				smtmp <- eval(parse(text = terms[smooth[k]]))
				T[[smooth[k]]] <- bayesx.setup(smtmp,file,outfile)
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
		if(R>0)
			{
			for(k in 1:R)
				{
				if(byra[k]=="NA")
					addterm <- c(addterm,paste(terms[ran[k]],"(random)",sep=""))
				else
					addterm <- c(addterm,paste(byra[k],"*",terms[ran[k]],"(random)",sep=""))
				T[[ran[k]]] <- list(term=terms[ran[k]],name=terms[ran[k]],
                                                    dim=1,class="random.spec",label=terms[ran[k]],byra=byra[k])
				}
			}
		}
	addterm <- paste(addterm,collapse=" + ")
	if(method=="MCMC")
		{
		fullformula <- paste("b.regress ",gp$response,"=",addterm,
                                     ", family=",family," iterations=",iter,
                                     " burnin=",burnin," step=",thinning,
                                     " predict using d",sep="")
		}
	else
		{
		fullformula <- paste("b.regress ",gp$response,"=",addterm,
                                     ", family=",family," eps=",eps," using d",sep="")
		}
	cat("dataset d\n",file=outfile,append=TRUE)
	if(data.file.check)
		cat(paste("d.infile using ",data.file,"\n",sep=""),file=outfile,append=TRUE)
	else
		cat(paste("d.infile using ",file,"/bayesx.data.raw\n",sep=""),file=outfile,append=TRUE)
	if(method=="MCMC")
		cat("bayesreg b\n",file=outfile,append=TRUE)
	else
		cat("remlreg b\n",file=outfile,append=TRUE)
	cat(paste("b.outfile = ",file,"/model\n",sep=""),file=outfile,append=TRUE)
	cat(fullformula,"\n",file=outfile,append=TRUE)
	if(method=="MCMC")
		cat("b.getsample\n",file=outfile,append=TRUE)
	bayesx.prg <- paste(paste(readLines(paste(file,"/bayesx.input.prg",sep=""),n=-1),collapse=" \n")," \n",sep="")
	specs$prg <- bayesx.prg
	specs$terms <- T
	setwd(fileo)

	return(invisible(specs))
	}


bayesx <- function(formula, family="gaussian", data=NULL, method="MCMC",
                   iter=1200, burnin=200, thinning=1, file=NULL)
	{
	res <- list()
	res$call <- match.call()
	res$formula <- formula <- as.formula(formula)
	res$bayesx.setup <- write.bayesx.input(formula,family,data,method,
                                               iter,burnin,thinning,file)
	# now estimate with BayesX
	res$bayesx.run <- run.bayesx(res$bayesx.setup$file)
	res$fout <- read.bayesx.output(res$bayesx.setup$file,method)
	res$fitted <- attr(res$fout,"fitted")
	attr(res$fout,"fitted") <- NULL
	class(res) <- "bayesx"
	return(res)
	}


run.bayesx <- function(file, prg.name="bayesx.input.prg")
	{
	fileo <- getwd()
	setwd(file)
	cat("Starting:\n")
	ptm <- proc.time()
	ok <- 0	
	# ok <- try(system(paste("BayesX ",file,"/",prg.name,sep="")))
	now <- proc.time()
	samptime <- now - ptm
	samptime <- samptime
	cat("Total run time was:",samptime[3],"sec\n")
	setwd(fileo)

	return(list(ok=ok,samptime=samptime[3]))
	}


read.bayesx.output <- function(file,method="MCMC")
	{
	mcheck <- FALSE
	if(method=="MCMC")
		mcheck <- TRUE
	fout <- list()
	fileo <- getwd()
	setwd(file)
	files <- list.files()
	usefiles <- file.specs <- eta <- NULL
	for(i in 1:length(files))
		{
		split <- strsplit(files[i],""," ")[[1]]
		n <- length(split)
		splitcheck <- paste(split[(n-2):n],collapse="")
		if(splitcheck=="res")
			{
			usefiles <- c(usefiles,files[i])
			file.specs <- c(file.specs,paste(split[(n-10):n],collapse=""))
			}
		if(mcheck)
			{
			if(n > 14)
				{
				get.eta <- paste(split[(n-14):n],collapse="")
				if(get.eta=="predictmean.raw")
					eta <- read.table(files[i],header=TRUE)
				}
			}
		}
	k <- 1
	etacheck <- !is.null(eta)
	smooth.mat <- var.mat <- s2mat <- NULL
	if(etacheck)
		{
		eta.names <- names(eta)
		response <- eta[eta.names[1]]
		if(mcheck)
			pred <- eta$linpred
		}
	for(i in 1:length(usefiles))
		{
		if(file.specs[i]=="pspline.res" || file.specs[i]=="spatial.res")
			{
			fout[[k]] <- read.bayesx.res(usefiles[i],etacheck,eta,response,pred,mcheck,FALSE)
			if(file.specs[i]=="pspline.res")
				class(fout[[k]]) <- "sm.gibbs"
			if(file.specs[i]=="spatial.res")
				{
				class(fout[[k]]) <- "mrf.gibbs"
				attr(fout[[k]],"coef") <- unique(fout[[2]][,c(2,3,4,5,6,7)])
				}
			attr(fout[[k]],"term.type") <- "smooth"
			smooth.mat <- rbind(smooth.mat,attr(fout[[k]],"variance"))
			k <- k + 1
			}
		if(file.specs[i]=="_random.res")
			{
			fout[[k]] <- read.bayesx.res(usefiles[i],etacheck,eta,response,pred,mcheck,TRUE)
			attr(fout[[k]],"term.type") <- "random"
			var.mat <- rbind(var.mat,attr(fout[[k]],"variance"))
			class(fout[[k]]) <- "random.gibbs"
			k <- k + 1
			}
		}
	if(any(grep("FixedEffects",files)))
		{
		draws <- NULL
		which <- grep("FixedEffects",files)
		for(i in 1:length(which))
			{
			if(any(grep("sample",files[which[i]])))
				{
				draws <- read.table(files[which[i]],header=TRUE)
				draws$intnr <- NULL
				draws <- t(as.matrix(draws))
				}
			else
				lin.mat <- read.table(files[which[i]],header=TRUE)
			}
		lin.mat$paramnr <- NULL
		varnames <- as.character(lin.mat$varname)
		if(varnames[1]=="const")
			varnames[1] <- "(Intercept)"
		lin.mat$varname <- NULL
		lin.mat <- as.matrix(lin.mat)
		cnames <- colnames(lin.mat)
		lin.mat <- lin.mat[,1:(ncol(lin.mat)-2)]
		colnames(lin.mat) <- cnames[1:(length(cnames)-2)]
		rownames(lin.mat) <- varnames
		attr(lin.mat,"coef.draws.utr") <- draws
		attr(fout,"lin.mat") <- lin.mat
		if(etacheck)
			{
			for(i in 1:length(varnames))
				{
				if(varnames[i]!="(Intercept)")
					{
					x <- unlist(eta[varnames[i]])
					xid <- rownames(lin.mat)==varnames[i]
					coefso <- as.numeric(lin.mat[xid,])
					coefs <- coefso[c(1,3,4,5,6,7)]
					ftmp <- NULL
					for(j in 1:length(coefs))
						ftmp <- cbind(ftmp,coefs[j]*x)
					pcat80 <- (ftmp[,3] < 0 & ftmp[,5] < 0)*(-1) + (ftmp[,3] <= 0 & ftmp[,5] >= 0)*0 + (ftmp[,3] > 0 & ftmp[,5] > 0)*1
					pcat95 <- (ftmp[,2] < 0 & ftmp[,6] < 0)*(-1) + (ftmp[,2] <= 0 & ftmp[,6] >= 0)*0 + (ftmp[,2] > 0 & ftmp[,6] > 0)*1
					ftmp <- cbind(x,ftmp,pcat95,pcat80)
					rownames(ftmp) <- NULL
					partial.resid <- response - pred + ftmp[,2]
					ftmp <- as.matrix(cbind(ftmp,partial.resid))
					colnames(ftmp) <- c(varnames[i],colnames(lin.mat)[c(1,3,4,5,6,7)],"pcat95","pcat80","partial.resid")
					class(ftmp) <- "linear.gibbs"
					attr(ftmp,"coef") <- coefso
					if(!is.null(draws))
						attr(ftmp,"coef.draws.utr") <- draws[xid,]
					fout[[k]] <- ftmp
					}
				}
			}
		}
	if(etacheck)
		fitted <- as.numeric(unlist(pred))
	else
		fitted <- "NA"
	if(any(grep("scale",files)))
		{
		which <- grep("scale",files)
		which2 <- any(grep("scale_sample",files))
		if(which2)
			which <- which[which!=grep("scale_sample",files)]
		if(length(which)<2)
			{
			s2mat <- read.table(files[which],header=TRUE)
			s2nam <- names(s2mat)
			s2mat <- matrix(s2mat,nrow=1)
			colnames(s2mat) <- s2nam
			rownames(s2mat) <- "Sigma2"
			attr(fitted,"variance") <- s2mat
			}
		if(which2)
			{
			which2 <- grep("scale_sample",files)
			if(length(which2)<2)
				{
				samples <- read.table(files[which2],header=TRUE)
				samples$intnr <- NULL
				attr(fitted,"variance.draws") <- as.numeric(unlist(samples))
				}
			}
		}
	attr(fout,"fitted") <- fitted
	attr(fout,"lin.mat") <- lin.mat
	attr(fout,"smooth.mat") <- smooth.mat

	setwd(fileo)
	
	return(fout)
	}


read.bayesx.res <- function(file,etacheck,eta,response,pred,mcheck,racheck)
	{
	fout <- read.table(file,header=TRUE)
	if(ncol(fout)>10)
		dim <- 2
	else
		dim <- 1
	fout$intnr <- NULL
	name <- colnames(fout)[1]
	raoke <- FALSE
	if(etacheck)
		{
		dat <- as.numeric(as.matrix(eta[name]))
		ind <- get.unique(dat)$ind[order(dat)]
		if(racheck)
			{
			partial.resid <- list()
			lev <- unlist(fout[name])
			levels <- levels(as.factor(lev))
			lev <- as.integer(lev[ind])
			e <- unlist(response - pred + fout$pmean[ind]) 
			for(j in 1:length(levels))
				partial.resid[[j]] <- e[lev==levels[j]]
			raoke <- TRUE
			}
		else
			{
			fout <- fout[ind,]
			if(mcheck)
				partial.resid <- response - pred + fout$pmean
			fout <- cbind(fout,partial.resid)
			names(fout)[length(names(fout))] <- "partial.resid"
			}
		}
	fout <- as.matrix(fout)
	rownames(fout) <- NULL
	attr(fout,"specs") <- list(dim=dim,term=name,label=name)
	if(raoke)
		attr(fout,"partial.resid") <- partial.resid
	if(mcheck)
		{
		filetmpl <- strsplit(file,""," ")[[1]]
		filetmpl <- paste(filetmpl[1:(length(filetmpl)-4)],collapse="")
		fsup <- c("_var.res","_variance_sample.raw","_lambda.res","_sample.raw")
		fsupspec <- c("variance","variance.draws","hyper","coef.draws.utr")
		fsupcheck <- list.files()
		for(k in 1:length(fsup))
			{
			ftxt <- paste(filetmpl,fsup[k],sep="",collapse="")
			if(ftxt%in%fsupcheck)
				{
				fget <- read.table(ftxt,header=TRUE)
				fget$intnr <- NULL
				fget <- as.matrix(fget)
				if(nrow(fget)!=2 && ncol(fget)!=9)
					fget <- t(fget)
				else
					fget <- fget[,1:(ncol(fget)-2)]
				if(!is.matrix(fget))
					{
					dimnam <- attr(fget,"names")
					fget <- matrix(fget,nrow=1)
					colnames(fget) <- dimnam
					}
				if(fsupspec[k]=="variance")
					rownames(fget) <- name
				attr(fout,fsupspec[k]) <- fget
				}
			else
				attr(fout,fsupspec[k]) <- 0
			}
		}

	return(fout)
	}


bayesx.setup <- function(object,path,outfile)
	{
	if(class(object)=="tp.smooth.spec")
		stop("BayesX does not support smooths of class tp.smooth.spec!") 
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
		nrknots <- object$bs.dim - object$p.order[1] + 1
		if(object$by=="NA")
			term <- paste(termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
		else
			term <- paste(object$by,"*",termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
		}
	if(class(object)=="ps.smooth.spec" && object$dim<2)
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
		nrknots <- object$bs.dim - object$p.order[1] + 1
		termo <- object$term
		constr <- object$xt$constr
		if(is.null(constr))
			{
			if(object$by=="NA")
				term <- paste(termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
			else
				term <- paste(object$by,"*",termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
			}
		else
			{
			if(constr<2)
				{
				if(object$by=="NA")
					term <- paste(termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],",monotone=increasing)",sep="")
				else
					term <- paste(object$by,"*",termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],",monotone=increasing)",sep="")
				}
			else
				{
				if(object$by=="NA")
					term <- paste(termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],",monotone=decreasing)",sep="")
				else
					term <- paste(object$by,"*",termo,"(psplinerw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],",monotone=decreasing)",sep="")
				}
			}
		}
	if(class(object)=="mrf.smooth.spec")
		{
		termo <- object$term
		if(object$by=="NA")
			term <- paste(object$term[1],"(spatial,map=bayesxmap)",sep="")	
		else
			term <- paste(object$by,"*",object$term[1],"(spatial,map=bayesxmap)",sep="")	
		class(object$xt$map) <- "bnd" 
		BayesX::write.bnd(object$xt$map,paste(path,"/bayesxmap.bnd",sep=""),replace=TRUE)
		cat("map bayesxmap\n",file=outfile,append=TRUE)
		cat("bayesxmap.infile using",paste(path,"/bayesxmap.bnd\n",sep=""),file=outfile,append=TRUE)
		}

	return(list(term=term,dim=object$dim,class=class(object),name=termo,label=object$label,by=object$by))
	}

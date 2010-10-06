parse.bayesx.input <- function(formula, family = "gaussian", data = NULL, 
                               weights = NULL, method = "MCMC", outfile = NULL, offset = NULL,
                               iter = 1200, burnin = 200, maxint = NULL, step = 1, aresp = 1, bresp = 0.005, 
			       distopt = "nb", reference = NULL, zipdistopt = "NA", 
                               eps = 1e-05, lowerlim = 0.001, maxit = 400, maxchange = 1e+06,
                               leftint = NULL, lefttrunc = NULL, state = NULL,...)
	{
	p <- list()
	if(is.null(data))
		data <- parent.frame()
	p$formula <- formula
	p$family <- family
	p$call <- match.call()
	p$terms <- terms(formula, specials = c("s","te","r"),keep.order=TRUE)
	p$data <- data
	p$weights <- weights
	p$method <- method
	p$outfile <- outfile
	p$offset <- offset
	p$iter <- iter
	p$burnin <- burnin
	p$maxint <- maxint
	p$step <- step
	p$aresp <- aresp
	p$bresp <- bresp
	p$distopt <- distopt
	p$reference <- reference
	p$zipdistopt <- zipdistopt
	p$eps <- eps
	p$lowerlim <- lowerlim
	p$maxit <- maxit
	p$maxchange <- maxchange
	p$leftint <- leftint
	p$lefttrunc <- lefttrunc
	p$state <- state

	return(p)
	}


write.bayesx.input <- function(object = NULL)
	{
	if(is.null(object))
		stop("nothing to do!")
	
	data.file=NULL
	if(is.null(object$data))
		object$data <- parent.frame()
	else
		{
		if(is.character(object$data))
			data.file <- object$data
		}
	if(is.null(object$outfile))
		object$outfile <- paste(tempdir(),"/bayesx",sep="")
	dir.create(object$outfile,showWarnings=FALSE)

	gp <- interpret.gam(object$formula)
	fake.labels <- attr(terms(gp$fake.formula), "term.labels")
	rafoo <- function(ra,data,outfile,k)
		{
		file.name.h <- NULL
		if(!is.null(ra$ins))
			{
			file.name <- file.name.h <- paste(outfile,"/bayesx.ra.data",k,".raw",sep="")	
			for(j in 1:length(ra$ins))
				file.name.h <- c(file.name.h,rafoo(ra$ins[[j]],data,outfile,as.integer(paste(k,j,sep=""))))
			fake.formula <- ra$fake.formula	
			dat <- model.matrix(fake.formula,data)
			dat.names <- colnames(dat)
			dat <- as.data.frame(dat)
			if(dat.names[1]=="(Intercept)")	
				dat <- dat[dat.names[2:ncol(dat)]]	
			write.table(dat,file.name,col.names=TRUE,row.names=FALSE,quote=FALSE)
			}
		return(file.name.h)
		}
	rah <- FALSE
	ra.file.names.h <- list()
	rk <- 1
	for(k in 1:length(fake.labels))
		{
		st <- strsplit(fake.labels[k],""," ")[[1]]
		if(length(st) > 1)
			{
			ft <- paste(st[1:2],collapse="")
			if(ft=="r(")
				{
				ra <- eval(parse(text=fake.labels[k]))
				fake.labels[k] <- ra$term
				if(!is.null(ra$ins))
					{
					rah <- TRUE
					ra.file.names.h[[rk]] <- rafoo(ra,object$data,object$outfile,k)
					rk <- rk + 1
					}
				}
			}
		}
	if(!rah)
		ra.file.names.h <- NULL
	if(is.null(data.file))
		{
		if(length(fake.labels)>1)
			for(k in 2:length(fake.labels))
				fake.labels[k] <- paste("+",fake.labels[k],sep="")	
		fake.formula <- paste("~-1+",paste(fake.labels,collapse=""),sep="",collapse="")
		dat <- model.matrix(as.formula(fake.formula),object$data)
		vars <- colnames(dat)
		resp <- eval(parse(text=gp$response),envir=object$data)
		dat <- cbind(resp,dat)
		colnames(dat) <- c(gp$response,vars)
		data.file <- paste(object$outfile,"/bayesx.data.raw",sep="")
		write.table(dat,data.file,col.names=TRUE,row.names=FALSE,quote=FALSE)
		}
	add.terms <- attr(object$terms,"term.labels")
	whatis <- rep("lin",length(add.terms))
	bt <- NULL
	sk <- attr(object$terms,"specials")$s
	tk <- attr(object$terms,"specials")$te
	rk <- attr(object$terms,"specials")$r
	if(!is.null(sk))
		sk <- sk - 1
	if(!is.null(tk))
		tk <- tk - 1
	if(!is.null(rk))
		rk <- rk - 1
	whatis[sk] <- "s"
	whatis[tk] <- "te"
	whatis[rk] <- "r"
	for(k in 1:length(add.terms))
		{
		if(whatis[k]=="lin")
			bt <- c(bt,add.terms[k])
		else
			{
			st <- eval(parse(text=add.terms[k]))
			bt <- c(bt,bayesx.construct(st))	
			}	
		}
	wd <- getwd()
	setwd(object$outfile)
	prg.file <- "bayesx.input.prg"
	cat(paste("% usefile ",object$outfile,"/bayesx.input.prg\n",sep=""),file=prg.file,append=FALSE)

	bt <- paste(bt,collapse=" + ")
	if(object$method=="MCMC")
		{
		fullformula <- paste("b.regress ",gp$response," = ",bt,
                                     ", family=",object$family," iterations=",object$iter,
                                     " burnin=",object$burnin," step=",object$step,
                                     " predict using d",sep="")
		}
	else
		{
		fullformula <- paste("b.regress ",gp$response," = ",bt,
                                     ", family=",object$family," eps=",object$eps," using d",sep="")
		}
	cat("dataset d\n",file=prg.file,append=TRUE)
	cat(paste("d.infile using ",data.file,"\n",sep=""),file=prg.file,append=TRUE)

	if(object$method=="MCMC")
		cat("bayesreg b\n",file=prg.file,append=TRUE)
	else
		cat("remlreg b\n",file=prg.file,append=TRUE)
	cat(paste("b.outfile = ",object$outfile,"/model\n",sep=""),file=prg.file,append=TRUE)
	cat(fullformula,"\n",file=prg.file,append=TRUE)
	if(object$method=="MCMC")
		cat("b.getsample\n",file=prg.file,append=TRUE)
	bayesx.prg <- paste(paste(readLines(paste(object$outfile,"/bayesx.input.prg",sep=""),n=-1),collapse=" \n")," \n",sep="")
	setwd(wd)

	return(invisible(list(prg=bayesx.prg,file.dir=object$outfile)))
	}


bayesx.construct <- function(object){ 
	UseMethod("bayesx.construct")}


bayesx.construct.ps.smooth.spec <- function(object)
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
				warning("order of the difference penalty not supported by BayesX, set to 2!")
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

	return(term)
	}


bayesx.construct.tensor.smooth.spec <- function(object)
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
				warning("order of the difference penalty not supported by BayesX, set to 2!")
				object$p.order <- c(object$p.order[1],2)
				}
			}
		}
	nrknots <- object$bs.dim - object$p.order[1] + 1
	if(object$by=="NA")
		term <- paste(termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")
	else
		term <- paste(object$by,"*",termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],",nrknots=",nrknots,",degree=",object$p.order[1],")",sep="")

	return(term)
	}


bayesx.construct.ra.smooth.spec <- function(object)
	{
	term <- object$term
	term <- paste(term,"(random)",sep="")

	return(term)
	}


bayesx.construct.rw1.smooth.spec <- function(object)
	{
	term <- object$term
	term <- paste(term,"(rw1)",sep="")

	return(term)
	}


bayesx.construct.rw2.smooth.spec <- function(object)
	{
	term <- object$term
	term <- paste(term,"(rw1)",sep="")

	return(term)
	}


r <- function(id, method = NULL, by = NA, xt = NULL)
	{
    	term <- deparse(substitute(id), backtick = TRUE, width.cutoff = 500)
    	by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
	ins <- fake.formula <- intcpt <- NULL
    	if(by.var == ".") 
        	stop("by=. not allowed")
    	if(term == ".") 
        	stop("r(.) not yet supported.")
	call <- match.call()
    	label <- paste("r(",term)
	if(is.null(method) && by.var=="NA")
		label <- paste(label,")",collapse="")
	if(!is.null(method) && by.var=="NA")
		{
		mlabel <- as.character(call[3])
		split <- strsplit(mlabel,""," ")[[1]]
		if(split[1]!="~")
			{
			label <- paste(label,",~",mlabel,")",collapse="")
			mf <- terms.formula(as.formula(paste("~",mlabel,collapse="")),specials=c("s","te","r"))
			}
		else
			{
			label <- paste(label,",",mlabel,")",collapse="")
			mf <- terms.formula(as.formula(mlabel),specials=c("s","te","r"))
			}
		intcpt <- attr(mf,"intercept")
		mterms <- attr(mf, "term.labels")
		for(k in 1:length(mterms))
			{
			ins[[k]] <- eval(parse(text=mterms[k]))
			if(!is.list(ins[[k]]))
				{
				ins[[k]] <-list(term=mterms[k],label=mterms[k])
				class(ins[[k]]) <- "lin.smooth.spec"
				}
			fake.formula <- c(fake.formula,ins[[k]]$term)
			}
		if(length(fake.formula)>1)
			for(k in 2:length(fake.formula))
				fake.formula[k] <- paste("+",fake.formula[k])
		if(intcpt > 0)
			fake.formula <- as.formula(paste("~",paste(fake.formula,collapse=""),collapse=""))
		else
			fake.formula <- as.formula(paste("~-1+",paste(fake.formula,collapse=""),collapse=""))
		}
	if(is.null(method) && by.var!="NA")
		label <- paste(label,",by=",by.var,")",collapse="")
	label <- paste(strsplit(paste(as.expression(label))," ","")[[1]],collapse="")
	ret <- list(term=term,label=label,by=by.var,xt=xt,ins=ins,fake.formula=fake.formula,call=call)
	class(ret) <- "ra.smooth.spec"

	ret 
	}


bayesx <- function(formula, family = "gaussian", data = NULL, 
                   weights = NULL, method = "MCMC", outfile = NULL, offset = NULL,
                   iter = 1200, burnin = 200, maxint = NULL, step = 1, aresp = 1, bresp = 0.005, 
		   distopt = "nb", reference = NULL, zipdistopt = "NA", 
                   eps = 1e-05, lowerlim = 0.001, maxit = 400, maxchange = 1e+06,
                   leftint = NULL, lefttrunc = NULL, state = NULL,...)
	{
	res <- list()
	res$call <- match.call()
	res$formula <- formula <- as.formula(formula)

	# setup files for bayesx
	res$bayesx.setup <- parse.bayesx.input(formula,family,data, 
                                               weights,method,outfile,offset,
                                               iter,burnin,maxint,step,aresp,bresp, 
		                               distopt,reference,zipdistopt, 
                                               eps,lowerlim,maxit,maxchange,
                                               leftint,lefttrunc,state,...)
	res$bayesx.prg <- write.bayesx.input(res$bayesx.setup)

	# now estimate with BayesX
	res$bayesx.run <- run.bayesx(res$bayesx.prg$file.dir)

	# get the output
	#res$fout <- term.order(res$bayesx.setup$term.labels,
        #                       read.bayesx.output(res$bayesx.setup$file,method))

	#res$terms <- res$bayesx.setup$term.labels
	#res$fitted <- attr(res$fout,"fitted")
	#attr(res$fout,"fitted") <- NULL
	#res$residuals <- attr(res$fout,"residuals")
	#attr(res$fout,"residuals") <- NULL
	#if(method=="MCMC")
		#{
		#res$DIC<- attr(res$fout,"DIC")
		#res$N <- res$bayesx.setup$N
		#res$iter <- res$bayesx.setup$iter
		#res$burnin <- res$bayesx.setup$burnin
		#res$thinning <- res$bayesx.setup$thinning
		#res$samptime <- res$bayesx.run$samptime
		#attr(res$fout,"DIC") <- NULL
		#}

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
	ok <- try(system(paste("BayesX ",file,"/",prg.name,sep="")))
	now <- proc.time()
	samptime <- now - ptm
	samptime <- samptime
	cat("Total run time was:",samptime[3],"sec\n")
	setwd(fileo)

	return(list(ok=ok,samptime=samptime))
	}


term.order <- function(terms,fin)
	{
	k <- length(terms)
	if(k==length(fin))
		{
		fattr <- attributes(fin)
		id <- rep(0,k)
		for(i in 1:k)
			for(j in 1:k)
				{
				name <- colnames(fin[[j]])[1]
				if(any(grep(name,terms[i])))
					id[i] <- j
				}
		fin <- fin[id]
		attributes(fin) <- fattr
		}

	return(fin)
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
		else
			{
			if(n > 11)
				{
				get.eta <- paste(split[(n-10):n],collapse="")
				if(get.eta=="predict.raw")
					eta <- read.table(files[i],header=TRUE)
				}
			}
		}
	k <- 1
	etacheck <- !is.null(eta)
	smooth.mat <- var.mat <- s2mat <- NULL
	if(etacheck)
		{
		if(mcheck)
			{
			eta.names <- names(eta)
			pred <- eta$linpred
			response <- eta[eta.names[1]]
			}
		else
			{
			pred <- eta$eta
			response <- 0
			}
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
			if(length(attr(fout[[k]],"variance"))>1)
				{
				smooth.mat <- rbind(smooth.mat,attr(fout[[k]],"variance"))
				k <- k + 1
				}
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
		if(!is.matrix(lin.mat))
			lin.mat <- matrix(lin.mat,nrow=1)
		colnames(lin.mat) <- cnames[1:(length(cnames)-2)]
		rownames(lin.mat) <- varnames
		attr(lin.mat,"coef.draws.utr") <- draws
		attr(fout,"lin.mat") <- lin.mat
		if(etacheck & mcheck)
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
					attr(ftmp,"term.type") <- "linear"
					fout[[k]] <- ftmp
					}
				}
			}
		}
	if(etacheck)
		{
		fitted <- as.numeric(unlist(pred))
		residuals <- response - fitted
		}
	else
		{		
		fitted <- "NA"
		residuals <- "NA"
		}
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
			s2mat <- matrix(unlist(s2mat),nrow=1)
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
	if(any(grep("deviance",files)))
		{
		which <- grep("deviance",files)
		devi <- read.table(files[which],header=TRUE)
		pd <- devi$unstandardized_deviance[length(devi$unstandardized_deviance)-1]
		DIC <- devi$unstandardized_deviance[length(devi$unstandardized_deviance)]
		attr(fout,"DIC") <- list(DIC=DIC,pd=pd)
		}
	attr(fout,"fitted") <- fitted
	attr(fout,"residuals") <- unlist(residuals)
	attr(fout,"lin.mat") <- lin.mat
	attr(fout,"smooth.mat") <- smooth.mat
	attr(fout,"var.mat") <- var.mat

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
	if(etacheck && mcheck)
		{
		dat <- as.numeric(as.matrix(eta[name]))
		ind <- get.unique(dat)$ind[order(dat)]
		if(racheck)
			{
			partial.resid <- list()
			lev <- unlist(fout[name])
			levels <- levels(as.factor(lev))
			lev <- as.integer(lev[ind])
			if(mcheck)
				e <- unlist(response - pred + fout$pmean[ind]) 
			else
				e <- unlist(response - pred + fout$pmode[ind]) 
			for(j in 1:length(levels))
				partial.resid[[j]] <- e[lev==levels[j]]

			raoke <- TRUE
			}
		else
			{
			fout <- fout[ind,]
			if(mcheck)
				partial.resid <- response - pred + fout$pmean
			else
				partial.resid <- response - pred + fout$pmode
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


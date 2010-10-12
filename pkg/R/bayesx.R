parse.bayesx.input <- function(formula, family = "gaussian", data = NULL, 
                               weights = NULL, method = "MCMC", outfile = NULL, offset = NULL,
                               iter = 1200, burnin = 200, maxint = 150, step = 1, aresp = 1, bresp = 0.005, 
			       distopt = "nb", reference = NULL, zipdistopt = "NA", 
                               begin = NULL, level1 = NULL, level2 = NULL,
                               eps = 1e-05, lowerlim = 0.001, maxit = 400, maxchange = 1e+06,
                               leftint = NULL, lefttrunc = NULL, state = NULL, ...)
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
	p$begin <- begin
	p$level1 <- level1
	p$level2 <- level2
	p$eps <- eps
	p$lowerlim <- lowerlim
	p$maxit <- maxit
	p$maxchange <- maxchange
	p$leftint <- leftint
	p$lefttrunc <- lefttrunc
	p$state <- state
	class(p) <- "bayesx.input"

	return(p)
	}


write.bayesx.input <- function(object = NULL)
	{
	if(is.null(object))
		stop("nothing to do!")
	if(class(object)!="bayesx.input")
		stop("object must be of class bayesx.input")
	
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
		for(i in 1:length(vars))
			{
			if(any(grep(":",vars[i])))
				{
				split <- strsplit(vars[i],":","")[[1]]
				vars[i] <- paste(split[1],"_",split[2],sep="")
				}
			}
		if(!is.null(object$offset))
			{
			dat <- cbind(dat,object$offset)
			off <- "offset"
			vars <- c(vars,off)
			}
		if(!is.null(object$weights))
			{
			dat <- cbind(dat,object$weights)
			wght <- "weights"
			vars <- c(vars,wght)
			}
		resp <- eval(parse(text=gp$response),envir=object$data)
		dat <- cbind(resp,dat)
		colnames(dat) <- c(gp$response,vars)
		data.file <- paste(object$outfile,"/bayesx.data.raw",sep="")
		write.table(dat,data.file,col.names=TRUE,row.names=FALSE,quote=FALSE)
		}

	wd <- getwd()
	setwd(object$outfile)
	prg.file <- "bayesx.input.prg"
	cat(paste("% usefile ",object$outfile,"/bayesx.input.prg\n",sep=""),file=prg.file,append=FALSE)

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
			{
			fcheck <- eval(parse(text=add.terms[k]),envir=object$data)
			if(is.factor(fcheck))
				{
				fcheck <- as.formula(paste("~",add.terms[k]))
				fcheck <- model.matrix(fcheck,data=object$data) 
				fcheck <- colnames(fcheck)[2:ncol(fcheck)]
				bt<-c(bt,fcheck)
				}
			else
				{
				if(any(grep(":",add.terms[k])))
					{
					split <- strsplit(add.terms[k],":","")[[1]]
					add.terms[k] <- paste(split[1],"_",split[2],sep="")
					}
				bt <- c(bt,add.terms[k])
				}
			}
		else
			{
			st <- eval(parse(text=add.terms[k]))
			bt <- c(bt,bayesx.construct(st,object$outfile,prg.file,object$data))	
			}	
		}
	if(!is.null(object$offset))
		bt <- c(bt,paste(off,"(offset)",sep=""))
	bt <- paste(bt,collapse=" + ")
	if(object$method=="MCMC")
		{
		fullformula <- paste("b.regress ",gp$response," = ",bt,
                                     ", family=",object$family," iterations=",object$iter,
                                     " burnin=",object$burnin," step=",object$step, 
                                     " maxint=",object$maxint," aresp=",object$aresp,
                                     " bresp=",object$bresp,sep="")
		if(!is.null(object$weights))
			fullformula <- paste(fullformula," weightvar=",wght,sep="")
		if(object$family=="nbinomial")
			fullformula <- paste(fullformula," distopt=",object$distopt,sep="")
		if(object$family=="multinomial")
			if(!is.null(object$reference))
				fullformula <- paste(fullformula," reference=",object$reference,sep="")
		if(object$family=="zip")
			if(object$distopt!="NA")
				fullformula <- paste(fullformula," distopt=",object$distopt,sep="")
		if(!is.null(object$begin))
			fullformula <- paste(fullformula," begin=",object$begin,sep="")
		if(!is.null(object$level1))
			fullformula <- paste(fullformula," level1=",object$level1,sep="")
		if(!is.null(object$level2))
			fullformula <- paste(fullformula," level2=",object$level2,sep="")
                fullformula <- paste(fullformula,"predict using d")
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


bayesx.construct <- function(object,file,prg,data){ 
	UseMethod("bayesx.construct")}


bayesx.construct.ps.smooth.spec <- function(object,file,prg,data)
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
	xt <- object$xt
	term <- paste(termo,"(psplinerw",object$p.order[2],",nrknots=",
                      nrknots,",degree=",object$p.order[1],sep="")
	if(object$by!="NA")
		term <- paste(object$by,"*",term,sep="")
	if(!is.null(xt$constr))
		{
		if(xt$constr<2)
			term <- paste(term,",monotone=increasing",sep="")
		else
			term <- paste(term,",monotone=decreasing",sep="")
		}
	if(!is.null(xt$a))
		term <- paste(term,",a=",xt$a,sep="")
	if(!is.null(xt$b))
		term <- paste(term,",b=",xt$a,sep="")
	if(!is.null(xt$contourprob))
		term <- paste(term,",contourprob=",xt$contourprob,sep="")
	if(!is.null(xt$derivative))
		term <- paste(term,",derivative",sep="")
	if(!is.null(xt$gridsize))
		term <- paste(term,",gridsize=",xt$gridsize,sep="")
	if(!is.null(xt$I))
		term <- paste(term,",I=",xt$I,sep="")
	if(!is.null(xt$lambda))
		term <- paste(term,",lambda=",xt$lambda,sep="")
	term <- paste(term,")",sep="")

	return(term)
	}


bayesx.construct.tensor.smooth.spec <- function(object,file,prg,data)
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
	xt <- object$xt
	term <- paste(termo[1],"*",termo[2],"(pspline2dimrw",object$p.order[2],
                      ",nrknots=",nrknots,",degree=",object$p.order[1],sep="")
	if(object$by!="NA")
		term <- paste(object$by,"*",term,sep="")
	if(!is.null(xt$a))
		term <- paste(term,",a=",xt$a,sep="")
	if(!is.null(xt$b))
		term <- paste(term,",b=",xt$a,sep="")
	if(!is.null(xt$contourprob))
		term <- paste(term,",contourprob=",xt$contourprob,sep="")
	if(!is.null(xt$derivative))
		term <- paste(term,",derivative",sep="")
	if(!is.null(xt$gridsize))
		term <- paste(term,",gridsize=",xt$gridsize,sep="")
	if(!is.null(xt$I))
		term <- paste(term,",I=",xt$I,sep="")
	if(!is.null(xt$lambda))
		term <- paste(term,",lambda=",xt$lambda,sep="")
	term <- paste(term,")",sep="")

	return(term)
	}


bayesx.construct.ra.smooth.spec <- function(object,file,prg,data)
	{
	term <- object$term
	term <- paste(term,"(random)",sep="")

	return(term)
	}


bayesx.construct.rw1.smooth.spec <- function(object,file,prg,data)
	{
	term <- object$term
	term <- paste(term,"(rw1)",sep="")

	return(term)
	}


bayesx.construct.rw2.smooth.spec <- function(object,file,prg,data)
	{
	term <- object$term
	term <- paste(term,"(rw1)",sep="")

	return(term)
	}


bayesx.construct.mrf.smooth.spec <- function(object,file,prg,data)
	{
	map <- object$xt$map
	if(is.null(map))
		{
		map <- object$xt[[1]]
		if(is.null(map))
			{
			map <- object$xt
			if(is.null(map))
				stop("need to supply a map object in argument xt!")
			if(!is.list(map) && !inherits(map,"bnd"))
				stop("need to supply a map object in argument xt!")
			}
		}
	if(!inherits(map,"bnd"))
		if(is.list(map))
			class(map) <- "bnd"
	bndfile <- paste(file,"/bayesx.map.bnd",sep="")
	prgfile <- paste(file,"/",prg,sep="")
	write.bnd(map=map,file=bndfile,replace=TRUE)
	cat("map bayesxmap\n",file=prgfile,append=TRUE)
	txt <- paste("bayesxmap.infile using",paste(file,"/bayesx.map.bnd",sep=""),"\n")
	cat(txt,file=prgfile,append=TRUE)
	term <- object$term
	term <- paste(term,"(spatial,map=bayesxmap)",sep="")

	return(term)
	}


fcheck <- function(labels,data)
	{
	labels <- as.list(labels)
	for(k in 1:length(labels))
		{
		tmp <- eval(parse(text=labels[[k]]),envir=data)
		if(is.factor(tmp))
			{
			fnam <- labels[[k]]
			tmp <- as.formula(paste("~",labels[[k]]))
			tmp <- model.matrix(tmp,envir=data)
			tmp <- colnames(tmp)[2:ncol(tmp)]
			labels[[k]] <- tmp
			attr(labels[[k]],"is.factor") <- TRUE
			attr(labels[[k]],"term") <- fnam
			}
		else
			{
			attr(labels[[k]],"is.factor") <- FALSE
			if(is.list(tmp))
				attr(labels[[k]],"term") <- tmp$term
			else
				attr(labels[[k]],"term") <- labels[[k]]
			}
		}
	return(labels)
	}


bayesx <- function(formula, family = "gaussian", data = NULL, 
                   weights = NULL, method = "MCMC", outfile = NULL, offset = NULL,
                   iter = 1200, burnin = 200, maxint = 150, step = 1, aresp = 1, bresp = 0.005, 
                   distopt = "nb", reference = NULL, zipdistopt = "NA", 
                   begin = NULL, level1 = NULL, level2 = NULL,
                   eps = 1e-05, lowerlim = 0.001, maxit = 400, maxchange = 1e+06,
                   leftint = NULL, lefttrunc = NULL, state = NULL, ...)
	{
	res <- list()
	res$call <- match.call()
	res$formula <- formula <- as.formula(formula)
	if(is.character(data))
		data <- read.table(data,header=TRUE)

	# setup files for bayesx
	res$bayesx.setup <- parse.bayesx.input(formula,family,data, 
                                               weights,method,outfile,offset,
                                               iter,burnin,maxint,step,aresp,bresp, 
		                               distopt,reference,zipdistopt, 
					       begin,level1,level2,
                                               eps,lowerlim,maxit,maxchange,
                                               leftint,lefttrunc,state,...)
	res$bayesx.prg <- write.bayesx.input(res$bayesx.setup)

	labels <- attr(res$bayesx.setup$terms,"term.labels")
	llabels <- fcheck(labels,data)

	# now estimate with BayesX
	res$bayesx.run <- run.bayesx(res$bayesx.prg$file.dir)

	# get the output
	res$fout <- term.order(llabels,read.bayesx.output(res$bayesx.prg$file.dir,method,labels))

	# maybe remove output folder
	if(!is.character(data))
		{
		wd <- getwd()
		setwd(res$bayesx.prg$file.dir)
		files <- list.files()
		for(k in 1:length(files))
			file.remove(files[k])
		setwd(wd)
		file.remove(res$bayesx.prg$file.dir)
		}

	res$terms <- res$bayesx.setup$term.labels
	res$fitted <- attr(res$fout,"fitted")
	attr(res$fout,"fitted") <- NULL
	res$residuals <- attr(res$fout,"residuals")
	attr(res$fout,"residuals") <- NULL
	if(method=="MCMC")
		{
		res$DIC<- attr(res$fout,"DIC")
		res$N <- length(res$residuals)
		res$samptime <- res$bayesx.run$samptime
		attr(res$fout,"DIC") <- NULL
		}

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
	m <- length(fin)
	fattr <- attributes(fin)
	id <- rep(0,k)
	for(i in 1:k)
		{
		tc <- attr(terms[[i]],"term")
		for(j in 1:m)
			{
			if(length(tc)<2)
				{
				if(ncol(fin[[j]])<11)
					{
					name <- colnames(fin[[j]])[2]
					if(tc==name)
						id[i] <- j
					}
				}
			else
				{
				if(ncol(fin[[j]])>10)
					{
					name <- colnames(fin[[j]])[2:3]
					if(tc[1]==name[1] && tc[2]==name[2])
						id[i] <- j
					}
				}
			}
		}
	fin <- fin[id]
	attributes(fin) <- fattr

	return(fin)
	}


read.bayesx.output <- function(file,method="MCMC",labels=NULL)
	{
	if(!is.null(labels))
		{
		nl <- length(labels)
		labs <- rep("",nl)
		for(k in 1:nl)
			{
			labtmp <- strsplit(labels[k],""," ")[[1]]
			if(length(labtmp)>3)
				{
				check <- paste(labtmp[1:2],sep="",collapse="")
				if(check%in%c("s(","r("))
					labtmp <- eval(parse(text=paste(labtmp,sep="",collapse="")))$label
				else
					{
					check <- paste(labtmp[1:3],sep="",collapse="")
					if(check=="te(")
						labtmp <- eval(parse(text=paste(labtmp,sep="",collapse="")))$label
					}
				}
			labels[k] <- paste(labtmp,sep="",collapse="")
			}
		}
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
			fout[[k]] <- read.bayesx.res(usefiles[i],etacheck,eta,response,pred,mcheck,FALSE,labels)
			if(file.specs[i]=="pspline.res")
				attr(fout[[k]],"type") <- "bayesx.pspline"
			if(file.specs[i]=="spatial.res")
				{
				attr(fout[[k]],"type") <- "bayesx.mrf"
				attr(fout[[k]],"coef") <- unique(fout[[k]][,c(2,3,4,5,6,7)])
				}
			if(length(attr(fout[[k]],"variance"))>1)
				{
				smooth.mat <- rbind(smooth.mat,attr(fout[[k]],"variance"))
				k <- k + 1
				}
			}
		if(file.specs[i]=="_random.res")
			{
			fout[[k]] <- read.bayesx.res(usefiles[i],etacheck,eta,response,pred,mcheck,TRUE,labels)
			attr(fout[[k]],"type") <- "bayesx.random"
			var.mat <- rbind(var.mat,attr(fout[[k]],"variance"))
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
					partial.resid <- cbind(x,partial.resid)
					colnames(partial.resid) <- c(varnames[i],"partial.resid")
					ftmp <- as.matrix(ftmp)
					colnames(ftmp) <- c(varnames[i],colnames(lin.mat)[c(1,3,4,5,6,7)],"pcat95","pcat80")
					ftmp <- ftmp[order(ftmp[,2]),]
					intnr <- 1:nrow(ftmp)
					ftmp <- cbind(intnr,ftmp)
					fout[[k]] <- as.data.frame(ftmp)
					attr(fout[[k]],"coef") <- coefso
					if(!is.null(draws))
						attr(fout[[k]],"coef.draws.utr") <- draws[xid,]
					attr(fout[[k]],"type") <- "bayesx.linear"
					attr(fout[[k]],"partial.resid") <- partial.resid
					k <- k + 1
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


read.bayesx.res <- function(file,etacheck,eta,response,pred,mcheck,racheck,labels)
	{
	fout <- read.table(file,header=TRUE)
	if(ncol(fout)>10)
		{
		dim <- 2
		name <- colnames(fout)[2:3]
		}
	else
		{
		dim <- 1
		name <- colnames(fout)[2]
		}
	raoke <- FALSE
	if(etacheck && mcheck)
		{
		dat <- as.numeric(as.matrix(eta[name[1]]))
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
			if(mcheck)
				partial.resid <- response - pred + fout$pmean[ind]
			else
				partial.resid <- response - pred + fout$pmode[ind]
			if(dim>1)	
				{
				partial.resid <- cbind(fout[ind,2],fout[ind,3],partial.resid)
				colnames(partial.resid) <- c(colnames(fout)[2:3],"partial.resid")
				}
			else
				{
				partial.resid <- cbind(fout[ind,2],partial.resid)
				colnames(partial.resid) <- c(colnames(fout)[2],"partial.resid")
				}
			}
		}
	fout <- as.matrix(fout)
	fout <- fout[order(fout[,2]),]
	rownames(fout) <- NULL
	fout <- as.data.frame(fout)
	if(!is.null(labels))
		{
		for(k in 1:length(labels))
			if(any(grep(name[1],labels[k])))
				{
				if(dim>1)
					{
					if(any(grep(name[2],labels[k])))
						lname <- labels[k]
					}
				else
					lname <- labels[k]
				}
		}
	else
		lname <- name
	if(dim>1 && is.null(labels))
		lname <- name <- paste(lname[1],"_",lname[2],sep="")
	attr(fout,"specs") <- list(dim=dim,term=name,label=lname)
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
					rownames(fget) <- lname
				attr(fout,fsupspec[k]) <- fget
				}
			else
				attr(fout,fsupspec[k]) <- 0
			}
		}

	return(fout)
	}


plot.bayesx <- function(x, which = 1, resid = FALSE, map = NULL, xa = 2, y = c(3, 4, 5, 7, 8), z = 4, ylim = NULL, 
                        lty = c(1, 2, 3, 2, 3), cols = rep(1, length(y)), month, year, step = 12, 
                        xlab, ylab, mode = 1, ticktype = "detailed", expand = 1, d = 1, theta = 40, phi = 40, 
		        regionvar=2, plotvar=3, limits, mcols="hcl", nrcolors=100, 
                        swapcolors=FALSE, pcat=FALSE, hcl.par=list(h=c(130,25), c=100, l=c(90,70)), 
                        hsv.par=list(s=1, v=1), legend=TRUE, drawnames=FALSE, cex.names=0.7, 
                        cex.legend=0.7, mar.min=2, density=15, ...)
	{
	if(class(x)!="bayesx")
		stop("argument x must be a bayesx object")
	X <- x$fout[[which]]
	if(!is.null(attr(X,"type")))
		{
		if(attr(X,"type")=="bayesx.pspline")
			{
			if(ncol(X)>10)
				BayesX::plotsurf(X,xa,y[1],z,mode,ticktype,expand,d,theta,phi,...)
			else
				{
				if(resid)
					{
					resids <- attr(X,"partial.resid")
					ylim <- range(resids[,2])
					}
				BayesX::plotnonp(X,xa,y,ylim,lty,cols,month,year,step,xlab,ylab,...)
				if(resid)
					points(x=resids[,1],y=resids[,2])
				}
			}
		if(attr(X,"type")=="bayesx.linear")
			{
			if(resid)
				{
				resids <- attr(X,"partial.resid")
				ylim <- range(resids[,2])
				}
			BayesX::plotnonp(X,xa,y,ylim,lty,cols,month,year,step,xlab,ylab,...)	
			if(resid)
				points(x=resids[,1],y=resids[,2])	
			}
		if(attr(X,"type")=="bayesx.random")
			{
			if(x$bayesx.setup$method=="MCMC")
				{
				id <- X[,2]
				eff <- as.matrix(cbind(X$pqu2p5,X$pqu10,X$pmed,X$pqu90,X$pqu97p5))
				grp <- peff <- NULL
				for(i in 1:length(id))
					{
					grp <- c(grp,rep(id[i],ncol(eff)))
					peff <- c(peff,eff[i,])
					}
				if(missing(xlab))
					xlab <- NULL
				if(missing(ylab))
					ylab <- NULL
				boxplot(peff~grp,...)
				if(!is.null(xlab))
					mtext(text=xlab,side=1,line=3)
				if(!is.null(ylab))
					mtext(text=ylab,side=2,line=3)
				}
			}

		if(attr(X,"type")=="bayesx.mrf")
			{     
			if(is.null(map))
				stop("a map object should be supplied!")
			if(!inherits(map,"bnd"))
				{
				class(map) <- "bnd"
				attr(map,"height2width") <- 1
				}
			if(is.null(attr(map,"height2width")))
				attr(map,"height2width") <- 1
     			BayesX::drawmap(data=X,map=map,regionvar=regionvar,plotvar=plotvar,limits=limits,
                                  cols=mcols,nrcolors=nrcolors,swapcolors=swapcolors,pcat=pcat,hcl.par=hcl.par,
                                  hsv.par=hsv.par,legend=legend,drawnames=drawnames,cex.names=cex.names,
                                  cex.legend=cex.legend,mar.min=mar.min,density=density,...)
			}
		}
	} 


print.bayesx <- function(x,...)
	{
	cat("Call:\n")
	print(x$call)
	cat("-------------------------------------------------------------------\n")
	cat("DIC =", x$DIC$DIC, " pd =", x$DIC$pd, " n =", x$N, " samptime =",paste(round(x$bayesx.run$samptime[3],2),"sec",sep=""), "\n")
	cat("Iterations =", x$bayesx.setup$iter, " burnin =", x$bayesx.setup$burnin, " step =", x$bayesx.setup$step, "\n")
	}


summary.bayesx <- function(object, digits = 4,...)
	{
	res <- list()
	res$call <- object$call
	res$DIC <- object$DIC
	res$setup <- object$bayesx.setup
	if(!is.null(attr(object$fout,"lin.mat")))
		res$lin.mat <- round(attr(object$fout,"lin.mat"),digits)
	if(!is.null(attr(object$fout,"smooth.mat")))
		res$smooth.mat <- round(attr(object$fout,"smooth.mat"),digits)
	if(!is.null(attr(object$fout,"var.mat")))
		res$var.mat <- round(attr(object$fout,"var.mat"),digits)
	res$S2 <-  round(attr(object$fitted,"variance"),digits)
	res$N <- object$N
	res$DIC <- object$DIC
	res$samptime <- round(object$bayesx.run$samptime[3],2)
	res$iter <- object$bayesx.setup$iter
	res$burn <- object$bayesx.setup$burn
	res$step <- object$bayesx.setup$step
	res$method <- object$bayesx.setup$method

	class(res) <- "summary.bayesx"
	res
	}


print.summary.bayesx <- function(x,...)
	{
	cat("Call:\n")
	print(x$call)
	tcheck <- rownames(x$lin.mat)
	tcheck2 <- as.character(x$lin.mat)
	if(!is.null(x$smooth.mat))
		{
		tcheck <- c(tcheck,rownames(x$smooth.mat))
		tcheck2 <- c(tcheck2,as.character(x$smooth.mat))
		}
	if(is.null(x$smooth.mat) && is.null(x$lin.mat))
		{
		tcheck <- c(tcheck,rownames(x$S2))
		tcheck2 <- c(tcheck2,as.character(x$S2))
		}
	ltcheck <- length(tcheck)
	maxsize <- rep(0,ltcheck)
	for(j in 1:ltcheck)
		maxsize[j] <- length(strsplit(tcheck[j],"","")[[1]])
	maxsize <- max(maxsize)

	ltcheck2 <- length(tcheck2)
	maxsize2 <- rep(0,ltcheck2)
	for(j in 1:ltcheck2)
		maxsize2[j] <- length(strsplit(tcheck2[j],"","")[[1]])
	maxsize2 <- max(maxsize2)

	if(maxsize2 > 7)
		maxsize <- maxsize + maxsize2*7 + 7	
	else
		maxsize <- maxsize + maxsize2*7 + 7
	liner <- paste(rep("-",maxsize),collapse="")

	fc <- TRUE
	if(nrow(x$lin.mat) < 2)
		{
		if(all(x$lin.mat[1,]==0))
			fc <- FALSE
		}
	else
		{
		if(all(x$lin.mat[1,]==0))
			{
			m <- ncol(x$lin.mat)
			nc <- colnames(x$lin.mat)
			nr <- rownames(x$lin.mat)[2:nrow(x$lin.mat)]
			x$lin.mat <- matrix(x$lin.mat[2:nrow(x$lin.mat),],ncol=m)
			colnames(x$lin.mat) <- nc
			rownames(x$lin.mat) <- nr
			}
		}
	if(fc || (!is.null(smooth.mat)))
		{
		cat(liner,"\n")
		cat("Fixed effects estimation results:\n")
		cat("---\n")
		}
	if(fc)
		{
		cat("Parametric Coefficients:\n")
		printCoefmat(x$lin.mat)
		}
	
	if(!is.null(x$smooth.mat))
		{
		if(fc)
			cat("-\n")
		cat("Smooth terms variances:\n")
		ls <- ncol(x$smooth.mat)
		terms <- colnames(x$smooth.mat)	
		rn <- rownames(x$smooth.mat)
		#for(j in 1:ls)
			#{
			#cat(terms[j],"\n")
			#newmat <- matrix(x$smoothhyp[j,],nrow=1,ncol=ncol(x$smoothhyp))
			#rownames(newmat) <- ""
			#colnames(newmat) <- rn
			#printCoefmat(newmat)
			#}
		printCoefmat(x$smooth.mat)
		}
	cat(liner,"\n")

	if(!is.null(x$var.mat))
		{
		cat("Random effects estimation results:\n")
		cat("---\n")

		if(!is.null(x$ra$coeflin))
			{
			cat("Parametric Coefficients:\n")
			printCoefmat(x$ra$coeflin)
			cat("---\n")
			}
		if(!is.null(x$ra$smoothhyp))
			{
			cat("Smooth terms variances:\n")
			printCoefmat(x$ra$smoothhyp)
			cat("---\n")
			}
		cat("Random variance estimation results:\n")
		cat("---\n")
		printCoefmat(x$var.mat)		
		cat(liner,"\n")
		}
	cat("Global variance estimation results:\n")
	printCoefmat(x$S2)
	cat(liner,"\n")
	if(x$method=="MCMC")
		{
		cat("DIC = ",x$DIC$DIC,"  pd = ",x$DIC$pd,"  n = ",x$N,"  samptime = ",x$samptime,"sec\n", sep="")
		cat("Iterations =",x$iter," burnin =",x$burn," step =",x$step,"\n")
		}
	else
		cat("BIC = ",x$DIC$DIC,"  df = ",x$DIC$pd,"  n = ",x$N,"  estimtime = ",x$samptime,"sec\n", sep="")
	}

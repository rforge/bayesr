r <- function(fac, method = NULL, by = NULL, xt = NULL)
	{
	missme <- missing(fac)
	if(missme && is.null(by))
		stop("Level variable missing in ra()!")

	call <- match.call(expand.dots=FALSE)
	get <- collect(2)

	allvars <- all.vars(call)
	if(missing(fac))
		{
		fact <- fac <- 1
		only <- TRUE
		}
	else
		{
		fact <- allvars[1]
		only <- FALSE
		}
	byvars <- ""
	if(!is.null(by))
		{
		if(only && is.null(method))
			initcall2 <- call[2]
		if(only && !is.null(method))
			initcall2 <- call[3]
		if(!only && !is.null(method))
			initcall2 <- call[4]
		if(!only && is.null(method))
			initcall2 <- call[3]
		if(is.null(eval(parse(text=as.character(initcall2)))))
			initcall2 <- call[3]

		if(!any(grep("~",initcall2)))
			form <- as.formula(paste("~",initcall2))
		else
			form <- as.formula(paste(initcall2))
		byvars <- all.vars(form)
		}

	pr <- sys.parents()
	isthere <- TRUE
	for(k in 1:length(pr))
		{
		ev <- sys.frame(pr[k])
		gf <- ls(envir=ev)
		if("usedata" %in% gf)
			{
			isthere <- FALSE
			usedata <- get("usedata",ev)
			}
		}
	if(isthere)
		{
		data <- NULL
		for(k in 1:length(pr))
			{
			ev <- sys.frame(pr[k])
			gf <- ls(envir=ev)
			if("dcheck" %in% gf)
				{
				dcheck <- get("dcheck",ev)
				if(dcheck)
					data <- get("data", envir = ev)					
				}
			}
		thisisitev <- new.env()
		usedata <- list()
		if(is.null(data))
			data <- .GlobalEnv
		if(!only)
			{
			fac <- tmp <- eval(parse(text = fact),envir=data)
			if(is.numeric(fac))
				fac <- tmp <- as.factor(round(fac))
			eval(parse(text = paste("usedata$",fact," <- tmp")))
			lef <- length(fac)
			startme <- 2
			if(is.null(method))
				startme <- 1
			}
		else
			lef <- startme <- 1
		for(j in startme:length(allvars))
			{
			if(!is.na(allvars[j]))
				{
				tmp <- eval(parse(text=allvars[j]),envir=data)
				if(length(tmp) == lef && !allvars[j]%in%byvars && j!=1)
					tmp <- tr(tmp,c(fac))
				eval(parse(text = paste("usedata$",allvars[j],"<-tmp",sep="")))
				}
			}
		if(only)
			{
			fac <- rep(1,length(usedata[[1]]))
			usedata$nofac <- fac
			fact <- paste(strsplit(paste(as.expression(call))," ","")[[1]],collapse="")
			}
		assign("usedata",usedata,envir=thisisitev)
		}
	data <- usedata
	if(!is.factor(fac))
		fac <- as.factor(round(fac))
	nl <- nlevels(fac)
	Z <- diag(nl)[fac,]
	if(is.vector(Z))
		Z <- matrix(1,nrow=length(Z),ncol=1)
	bystructure <- byterms <- byrownid <- byfacid <- rintc <- NULL
	K <- diag(1,nl)
	if(!is.null(by))
		{
		K <- NULL
		bystructure <- list()
		TF <- terms.formula(form)
		terms <- byterms <- attr(TF, "term.labels")
		last <- 1
		termstruct <- matrix(0,length(terms),2)
		rownames(termstruct) <- terms
		colnames(termstruct) <- c("start","end")
		tmpa <- specs <- vector("list",length=length(terms))
		for(d in 1:length(terms))
			{
			tmpa[[d]] <- eval(parse(text = terms[d]))
			if(is.list(tmpa[[d]]))
				{
				tmpa[[d]] <- stewrap(tmpa[[d]],data,terms[d])
				if(class(tmpa[[d]]) == "mrf.smooth")
					{
					smx <- tmpa[[d]]$mrfind
					tmpa[[d]]$names <- tmpa[[d]]$names[1]
					}
				else
					{
					smx <- NULL
					for(i in 1:tmpa[[d]]$specs$dim)
						smx <- cbind(smx,eval(parse(text=tmpa[[d]]$specs$term[i]),envir=data))
					}
				specs[[d]] <- list(specs=tmpa[[d]]$specs,x=smx,term=terms[d],bytype="smooth")
				}
			else
				specs[[d]] <- list(bytype="linear",x=tmpa[[d]],term=terms[d])
			}
		rintc <- FALSE
		if(strsplit(as.character(initcall2),"","")[[1]][1] == "1")
			{
			rintc <- TRUE
			last <- 0
			out <- NULL
			termstruct <- rbind(c(0,0),termstruct)
			rownames(termstruct) <- c("(Intercept)",terms)
			
			for(j in 1:ncol(Z))
				{
				byfacid <- c(byfacid,j)
				byrownid <- c(byrownid,paste("(Intercept):",j,sep=""))
				check <- Z[,j] != 0
				out <- cbind(out,Z[,j])
				last <- last + 1
				termstruct[1,] <- c(last,last)
				K <- matrix(0,1,1)
				for(d in 1:length(terms))
					{
					if(is.list(tmpa[[d]]))
						{
						K <- blockme(K,tmpa[[d]]$S[[1]])
						tmp <- tmpa[[d]]$X
						nct <- ncol(tmp)
						byfacid <- c(byfacid,rep(j,nct))
						txt <- paste(strsplit(sub(tmpa[[d]]$specs$term,paste(tmpa[[d]]$specs$term,":", j,sep=""),
                                                 	     terms[d])," ","")[[1]],collapse="")
						rownames(termstruct)[d + 1] <- txt
						byrownid <- c(byrownid,paste(txt,":",1:nct,sep=""))
						termstruct[d+1,] <- c(last + 1,last + nct)
						last <- last + ncol(tmp)
						}
					else
						{
						K <- cbind(0,blockme(K,diag(1,1)))
						txt <- paste(terms[d],":",j,sep="") 
						byrownid <- c(byrownid,txt)
						byfacid <- c(byfacid,j)
						tmp <- tmpa[[d]]
						last <- last + 1
						termstruct[d+1,] <- c(last,last)
						rownames(termstruct)[d + 1] <- txt
						}
					tmp <- Z[,j]*tmp
					out <- cbind(out,tmp)
					}
				attr(termstruct,"used") <- check
				bystructure[[j]] <- termstruct
				}
			}		
		else
			{
			out <- NULL
			for(j in 1:ncol(Z))
				{
				check <- Z[,j] != 0
				for(d in 1:length(terms))
					{
					if(is.list(tmpa[[d]]))
						{
						K <- blockme(K,tmpa[[d]]$S[[1]])
						tmp <- tmpa[[d]]$X
						nct <- ncol(tmp)
						byfacid <- c(byfacid,rep(j,nct))
						withthis <- paste(tmpa[[d]]$specs$term,":", j,sep="")
						txt <- terms[d]
						for(jj in 1:length(withthis))
							{
							txt <- paste(strsplit(sub(tmpa[[d]]$specs$term[jj],withthis[jj],
                                                                     txt)," ","")[[1]],collapse="")
							}
						rownames(termstruct)[d] <- txt
						byrownid <- c(byrownid,paste(txt,":",1:nct,sep=""))
						termstruct[d,] <- c(last,last + nct - 1)
						last <- last + ncol(tmp)
						}
					else
						{
						K <- blockme(K,diag(1,1))
						txt <- paste(terms[d],":",j,sep="")
						byrownid <- c(byrownid,txt)
						byfacid <- c(byfacid,j)
						tmp <- tmpa[[d]]
						termstruct[d,] <- c(last,last)
						rownames(termstruct)[d] <- txt
						last <- last + 1
						}
					tmp <- Z[,j]*tmp
					out <- cbind(out,tmp)
					}
				attr(termstruct,"used") <- check
				bystructure[[j]] <- termstruct
				}
			}
		attr(bystructure,"specs") <- specs
		colnames(out) <- rep("",ncol(out))
		Z <- out
		nl <- ncol(Z)
		}
	if(!is.null(method))
		{
		if(missme && !is.null(by))
			initcall <- call[2]
		else
			initcall <- call[3]
		if(!any(grep("~",initcall)))
			form <- as.formula(paste("~",initcall))
		else
			form <- as.formula(paste(initcall))
		vars <- all.vars(form)

		F <- length(fac)

		TF <- terms(form, specials = "r")
		terms <- attr(TF, "term.labels")
		if(attr(TF,"intercept") > 0)
			terms <- c("1",terms)
		specials <- attr(TF,"specials")$r
		lt <- length(terms)
		speccheck <- rep(FALSE,lt)
		if(attr(TF,"intercept") > 0)
			speccheck[specials+1] <- TRUE 
		else
			speccheck[specials] <- TRUE
		indspec <- matrix(0,3,lt)
		itr1 <- itr2 <- itr3 <- 0

		RA <- vector("list",lt)
		RA$hspec <- FALSE

		for(k in 1:lt)
			{
			oterm <- terms[k]
			nf <- as.formula(paste("~",terms[k]))

			rwcheck <- paste(strsplit(terms[k],"","")[[1]][1:3],collapse="")
			rwcheck <- rwcheck == "rw("
			if(rwcheck)
				vs <- vsout <- terms[k]
			else
				vs <- vsout <- all.vars(nf)				
			
			matchtxt2 <- c("s(","te(")
			matchtxt <- c("s(","te(")

			yes <- here <- FALSE
      			whatis <- 0
			spec <- 1

			if(speccheck[k])
				{
				itr3 <- itr3 + 1
				indspec[3,k] <- itr3
				RA[[k]] <- eval(parse(text = terms[k]),envir=data)
				}
			else
				{
				splits <- strsplit(terms[k],"","")[[1]]
				sp2 <- paste(splits[1:2],collapse="")
				sp3 <- paste(splits[1:3],collapse="")

				for(p in 1:2)
					{
					if(matchtxt2[p]==sp2)
						{
						take <- matchtxt[p]
						yes <- TRUE
						spec <- 2
						}	
					else
						{
						if(matchtxt2[p]==sp3)
							{
							take <- matchtxt[p]
							yes <- TRUE
              						whatis <- 1
							spec <- 2
							}	
						}
					}
				if(rwcheck)
					spec <- 4
				if(any("F"%in%vs))
					vs <- vs[vs!="F"]
				if(any("T"%in%vs))
					vs <- vs[vs!="F"]
				lvs <- length(vs)

				if(any(grep("[$]",terms[k])))
					{
					if(yes)
						{	
						if(lvs == 2)
							{				
							split <- unlist(strsplit(terms[k],","))
							split <- unlist(strsplit(split[1],take))[2]
							vs <- sub("[(]","",split)
							vs <- sub("[)]","",vs)
							lvs <- 1
							}
						if(lvs == 3)
							{
							split1 <- unlist(strsplit(terms[k],","))
							split2 <- unlist(strsplit(split1[1],take))[2]
							vs <- split2[1]
							vs <- sub("[(]","",vs)
							vs <- sub("[)]","",vs)
							vs2 <- sub(" ","",split1[2])
							vs <- c(vs,vs2)
							lvs <- 2
							}
						}
					else
						{
						vs <- paste(vs[1],"$",vs[2],sep="")
						lvs <- 1
						}
					}
				if(spec != 4 & oterm!="1")
					{
					obs <- rep(0,lvs)
					VS <- vector("list",lvs)
					for(j in 1:lvs)
						{
						VS[[j]] <- eval(parse(text = vs[j]),envir=data)
						obs[j] <- length(VS[[j]])
						}
					if(lvs > 1)
						if(obs[1] != obs[2])
							if(!any(grep("mrf",terms[k])))
								{
								stoptxt <- paste("Length of variable ",vs[1]," is differing from ",vs[2],"!",sep="")
								stop(stoptxt)
								}
					L <- obs[1]

					if(L > F)
						{
						text <- paste("Too many observations in term ",terms[k],"!",sep="")
						stop(text)
						}	
					if(L == F)
						{
						here <- TRUE
						if(yes)
							{
							takenow <- unlist(strsplit(take,""))
							takenow <- paste(takenow[1:(length(takenow)-1)],collapse="")

              						if(whatis==1)
                						split <- vs
              						else
                						{
                						split <- unlist(strsplit(terms[k],takenow))[2]
								split <- unlist(strsplit(split,""))
								split <- paste(split[2:(length(split)-1)],collapse="")
               	 						if(j > 1)
									split <- unlist(strsplit(split,","))
								for(j in 1:length(split))
									split[j] <- paste(strsplit(split[j]," ")[[1]],sep="",collapse="")
                						}
							for(j in 1:lvs)
								{
								if(j==2)
                    							{
                    							if(whatis==0)
                      								{
                      								var <- tr(c(VS[[j]]),c(fac))
                      								check <- getc <- paste("new",k,j,sep="")
                      								eval(parse(text=paste(check,"<- var")))
                      								split[j] <- getc
                      								}
                    							}
                						else
                    							{
                    							var <- tr(c(VS[[j]]),c(fac))
                    							check <- getc <- paste("new",k,j,sep="")
                    							eval(parse(text=paste(check,"<- var")))
                    							split[j] <- getc
                    							}
								if(j < 2)
									{
                  							if(whatis==0 || take == "mrf(")
                    								split[j] <- paste("(",split[j],sep="")
									}
								}
							scheck <- FALSE
							if(length(split) > 1)
								scheck <- TRUE
							split <- paste(split,collapse=",")
							terms[k] <- paste(takenow,split,")",sep="")
							}
						else
							{
							check <- getc <- paste("new",k,sep="")
							var <- tr(c(VS[[1]]),c(fac))
							eval(parse(text = paste(check," <- var",sep="")))
							terms[k] <- getc
							}
						}
					}
				if(spec == 2)
					{
					RA[[k]] <- stewrap(eval(parse(text = terms[k])),data,terms[k])
					itr2 <- itr2 + 1
					indspec[2,k] <- itr2
					# test 
					vs <- RA[[k]]$spec$term
					}
				if(spec == 1)
					{
					if(terms[k]=="1")
						{
						RA[[k]]$X <- rep(1,nl)
						RA[[k]]$term <- oterm <- "(Intercept)"
						RA[[k]]$mx <- 1
						}
					else
						{
						RA[[k]]$X <- c(eval(parse(text = terms[k]), envir=data))
						RA[[k]]$term <- terms[k]
						RA[[k]]$mx <- RA[[k]]$x[approxm(RA[[k]]$x)]
						}
					RA[[k]]$x <- c(RA[[k]]$X)
					itr1 <- itr1 + 1
					indspec[1,k] <- itr1 		
					}
				if(spec == 4)
					{
					if(is.null(by))
						stop("Need argument by to use method rw()!")
					what <- eval(parse(text = terms[k]))
					RA[[k]]$X <- diag(1,nl)
					RA[[k]]$S <- list(K)
					K <- RA[[k]]$X
					RA[[k]]$x <- 1:nl 
					RA[[k]]$mx <- rep(1,nl)
					RA[[k]]$ind <- 1:nl
					RA[[k]]$rwcheck <- 1
					RA[[k]]$center <- TRUE
					RA[[k]]$rwx <- what$x
					RA[[k]]$rwnames <- what$name
					RA[[k]]$rwtype <- tmpa[[1]]$type
					RA[[k]]$type <- "rw"
					tmp <- paste(tmpa[[1]]$names,":rw",sep="")
					RA[[k]]$rwspecs <- list(by=tmpa[[1]]$by,loc=tmpa[[1]]$loc,dim=tmpa[[1]]$dim,kappa=tmpa[[1]]$kappa,names=tmp,
                                               		        lambda=tmpa[[1]]$lambda,c=tmpa[[1]]$c,center=tmpa[[1]]$center,knots=tmpa[[1]]$knots,
								degree=tmpa[[1]]$degree,term=tmpa[[1]]$callsm,orderpenalty=tmpa[[1]]$orderpenalty)
					itr2 <- itr2 + 1
					indspec[2,k] <- itr2
					spec <- 2
					vs <- paste(vs,":",1:nl,sep="")
					}

				RA[[k]]$names <- vs
				RA[[k]]$spec <- spec
				RA[[k]]$term <- oterm
				RA[[k]]$hspec <- FALSE
				RA[[k]]$check <- FALSE

				if(k > 1)
					{
					if(!is.null(RA[[k]]$x) && !is.null(RA[[k-1]]$x))
						{
						xtmp1 <- is.matrix(RA[[k]]$x)
						xtmp2 <- is.matrix(RA[[k-1]]$x)
						if(xtmp1)
							xtmp1 <- nrow(RA[[k]]$x)
						else
							xtmp1 <- length(RA[[k]]$x)
						if(xtmp2)
							xtmp2 <- nrow(RA[[k-1]]$x)
						else
							xtmp2 <- length(RA[[k-1]]$x)
						if(xtmp1 != xtmp2)
							if(!any(grep("mrf",terms[k])))
								{
								stoptxt <- paste("Length of variable in term ",
                                                                                 RA[[k]]$term," is differing from ",
										 RA[[k-1]]$term,"!",sep="")
								stop(stoptxt)
								}
						}
					}
				if(get$gcheck)
					{
					RA[[k]]$tildeZ <- trans.design(RA[[k]]$X,RA[[k]]$S[[1]])
					RA[[k]]$ttildeZ <- t(RA[[k]]$tildeZ$tildeZ)
					RA[[k]]$dPR <- ncol(RA[[k]]$tildeZ$tildeZ)
					RA[[k]]$fr <- rep(0,nl)
					RA[[k]]$drawTildeGammaR <- matrix(0,RA[[k]]$dPR,get$keepers)
					RA[[k]]$IR <- rep(1,RA[[k]]$dPR)

					if(spec == 2)
						{
						RA[[k]]$Rlambda <- rep(0,get$keepers)
						RA[[k]]$fitR <- matrix(0,nl,get$keepers)
						RA[[k]]$tau2R <- 0.0001
						RA[[k]]$hyperaTau2newR <- get$hyperatau + (qr(RA[[k]]$S[[1]])$rank)/2
						}
					}
				}
			}
		RA$hspec <- TRUE
		RA$length <- lt
		RA$indspec <- indspec
		}
	if(is.null(method))
		{
		call <- match.call(expand.dots=FALSE)
		names <- all.names(call)
		names <- names[names!="r"] ### check
		
		RA <- vector("list",1)		
		RA[[1]]$names <- names
		RA[[1]]$term <- terms
		RA$length <- 1
		RA$hspec <- FALSE
		}

	nl <- ncol(Z)
	RA$Z <- Z

	RA$n <- nl
	RA$call <- paste(strsplit(paste(as.expression(call))," ","")[[1]],collapse="")
	RA$check <- TRUE
	if(only)
		RA$fac <- RA$facorig <- 1:ncol(Z)
	else
		{
		RA$fac <- as.integer(fac)
		RA$facorig <- fac
		}
	RA$mx <- as.numeric(incidentm2(RA$fac))
	RA$facname <- fact
	RA$byterms <- byterms
	RA$bywhat <- bystructure
	RA$K <- K
	RA$type <- "ra"
	RA$byfac <- byrownid
	RA$byfacid <- byfacid
	RA$rintc <- rintc
	RA$center <- 0
	if(!is.null(xt$center))
		if(xt$center)
			RA$center <- 1

	if(get$gcheck)
		{
		RA$ZZ <- as.numeric(diag(crossprod(Z,Z)))
		RA$tZ <- t(Z)
		RA$BR <- as.integer(nl)
		RA$drawBetaR <- as.numeric(rep(0,RA$n))
		RA$fR <- as.numeric(rep(0,nrow(Z)))
		RA$sigma2 <- 0.0001
		RA$drawB <- matrix(0,nl,get$keepers)
		RA$etaR <- as.numeric(rep(0,nl))
		RA$haRnew <- as.numeric(get$hyperasigma + qr(K)$rank/2)
		RA$S2R <- as.numeric(rep(0,get$keepers))
		RA$usedata <- usedata
		}

	return(RA)
	}

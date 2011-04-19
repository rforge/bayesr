m <- function(fac, method = NULL, sp = NULL, center = TRUE, hypident = FALSE)
	{
	if(missing(fac))
		stop("Level variable missing in m()!")

	call <- match.call(expand.dots=FALSE)

	if(!is.null(method))
		{
		initcall <- call[3]
		if(strsplit(as.character(initcall),"","")[[1]][1] == "~")
			{
			tmp <- strsplit(as.character(initcall),"","")[[1]]
			initcall <- paste(tmp[2:length(tmp)],collapse="")
			initcall <- parse(text=initcall)
			tmpcall <- paste(strsplit(as.character(initcall),"")[[1]][1:3],collapse="")
			names <- all.names(initcall)
			F <- L <- 1
			}
		else
			{
			tmpcall <- paste(strsplit(as.character(initcall),"")[[1]][1:3],collapse="")
			names <- all.names(initcall)
			if(mode(initcall) == "call")
				initcall <- paste(initcall)
			object <- eval(parse(text = initcall), envir=sys.frame(-1))
			data <- collect(2)$data
			is.sm <- FALSE
			if(is.list(object))
				{
				X <- smooth.construct(object,data,NULL)
				L <- nrow(X$X)
				is.sm <- TRUE
				}
			else
				L <- length(object)
			if(L == 0)
				stop("Variable of length zero in method in m()!")
			F <- length(fac)
			if((L < F || L > F) && method != 1)
				stop("Variable lengths differing in m()!")
			}
		if(L == F)
			{
			dat <- list()
		    	if(is.sm)
				{
				nx <- length(fac)
				if(!is.factor(fac))
					fac <- as.factor(fac)
				nl <- nlevels(fac)
				x <- matrix(NA,nrow=nx,ncol=object$dim)
				for(i in 1:object$dim)
					x[,i] <- eval(parse(text = object$term[i]), envir=data)
				IND <- 1:nx
				if(length(center) == 1)
					center <- rep(center,nl)
				levtxt <- levels(fac)
				levfac <- as.integer(levtxt)
				checkfac <- fac
				storage.mode(checkfac) <- "integer"

				bycheck <- FALSE; mby <- MBY <- NULL
				if(object$by!="NA")
					{
					bycheck <- TRUE
					mby <- eval(parse(text = object$by), envir=data)
					MBY <- object$by
					}
				x <- as.data.frame(cbind(x,mby))
				xtxt <- c(object$term,MBY)
				names(x) <- xtxt
				MU <- list()
				object <- eval(parse(text=initcall))
				ismrf <- (class(object) == "mrf.smooth.spec")
				mlab <- object$label
				
				for(j in 1:nl)
					{
					initcall2 <- initcall
					ind <- IND[checkfac==levfac[j]]
					X <- as.data.frame(x[ind,])
					xlabel <- paste(xtxt,".",j,sep="")
					names(X) <- xlabel
					for(jj in 1:length(xlabel))
						{
						eval(parse(text=paste("dat$",xlabel[jj],"<-X[,jj]",sep="")))
						initcall2 <- sub(xtxt[jj],xlabel[jj],initcall2)
						}
					if(bycheck)
						{
						object$by <- xlabel[object$dim+1]
						initcall2 <- sub(xtxt[length(xtxt)],xlabel[length(xlabel)],initcall2)
						}
					object$term <- xlabel[1:object$dim]		
					if(!is.null(object$margin))
						for(jj in 1:object$dim)
							{
							object$margin[[jj]]$term <- xlabel[jj]
							if(bycheck)
								object$margin[[jj]]$by <- xlabel[length(xlabel)]
							}
					object$label <- paste(mlab,".",levtxt[j],sep="")
					MU[[j]] <- stewrap(object,X,initcall2)	
					MU[[j]]$ind <- ind
					MU[[j]]$n <- nrow(X)
					MU[[j]]$lambda <- sp[j]
					MU[[j]]$center <- center[j]
					MU[[j]]$ismrf <- ismrf
					}
				MU$nl <- nl
				MU$smcheck <- 2
				MU$term <- object$term
				MU$fac <- fac
				MU$x <- x
				MU$hypident <- hypident
				MU$N <- nx
				MU$data <- dat

				return(MU)
				}
			else
				{
				f <- as.formula(paste("~",initcall))
				cm <- charmatch("1",as.character(initcall))
				v <- all.vars(f)
			
				lv <- length(v)
				if(lv > 1)
					{
					var <- list()
					for(k in 1:lv)
						var[[k]] <- eval(parse(text = v[k]), envir=sys.frame(-1))
					}

				nf <- length(fac)
				IND <- 1:nf
				if(!is.factor(fac))
					fac <- as.factor(fac)
				nl <- nlevels(fac)
				# Z <- diag(nl)[fac,]
				levfac <- as.integer(levels(fac))
				checkfac <- fac
				storage.mode(checkfac) <- "integer"

				MU <- list()
				if(is.na(cm) && length(v) == 1)
					{
					X <- list()
					for(j in 1:nl)
						{
						#X <- Z[,j]*method
						#ind <- IND[X!=0]
						ind <- IND[checkfac==levfac[j]]
						tmp <- method[ind]
						MU[[j]] <- list(X=tmp, S=list(NULL))
						MU[[j]]$ind <- ind
						MU[[j]]$n <- length(X)	
						MU[[j]]$mx <- tmp[approxm(tmp)]
						}
					MU$name <- v
					MU$ok <- TRUE
					}
				if(is.na(cm) && length(v) > 1)
					{
					X <- list()
					for(j in 1:nl)
						{
						#ind <- IND[Z[,j]!=0]
						ind <- IND[checkfac==levfac[j]]
						nv <- length(ind)
						mat <- matrix(0,nv,lv)
						mx <- rep(0,lv)
						for(k in 1:lv)
							{
							# V <- Z[,j]*var[[k]]
							# tmp <- V[V!=0]
							tmp <- var[[k]][ind]
							mx[k] <- tmp[approxm(tmp)]
							mat[,k] <- tmp
							}
						MU[[j]] <- list(X=mat, S=list(NULL))
						MU[[j]]$ind <- ind
						MU[[j]]$n <- nv	
						MU[[j]]$mx <- mx	
						}
					MU$name <- v
					}
				if(!is.na(cm))
					stop("Invalid method specified in m()!")

				# MU$trans <- diag(1,nf) bzw. nx
				# MU$Z <- diag(1,nf)
				MU$nl <- nl
				MU$smcheck <- 1
				MU$fac <- fac
				MU$hypident <- hypident
				MU$N <- length(fac)
				MU$specs <- list(term=v,call=v)

				return(MU)
				}
			}
		}
	if(is.null(method))
		stop("No method specified in m()!")
	}

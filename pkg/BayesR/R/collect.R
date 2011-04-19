collect <- function(frame)
	{
	gcheck <- FALSE
	ev <- sys.frame(frame)
	pr <- sys.parents()

	for(k in 1:length(pr))
		{
		ev <- sys.frame(pr[k])
		gf <- ls(envir=ev)

		if("keepers" %in% gf)
			{
			gcheck <- TRUE
			keepers <- get("keepers", envir=ev)
			response <- try(get("response", envir=ev),silent=TRUE)
			hyperasigma <- get("hyperasigma", envir=ev)
			hyperatau <- get("hyperatau", envir=ev)
			dcheck <- get("dcheck", envir=ev)
			if(dcheck)
				this <- get("this", envir=ev)
			else
				this <- data <- NULL
			data <- get("data", envir=ev)
			pframe <- get("pframe", envir=ev)

			return(list(gcheck=gcheck,keepers=keepers,hyperasigma=hyperasigma,dcheck=dcheck,
                         hyperatau=hyperatau,this=this,data=data,response=response))
			}
		}
	return(list(gcheck=gcheck))
	}
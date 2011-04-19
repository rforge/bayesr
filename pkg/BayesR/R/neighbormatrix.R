neighbormatrix <- function(map, type = 1, scale = NULL)
	{
	if(missing(map))
		stop("Map is missing!")
	if(!is.list(map))
		warning("Map object is not of list() type!")
	if(type < 1 || type > 3)
		warning("Wrong type selected, type is set to default!")
	N <- length(map)
	regionnames <- names(map)
	adjmatrix <- matrix(0,N,N)
	colnames(adjmatrix) <- regionnames
	rownames(adjmatrix) <- regionnames

	if(type == 2)
		{
		B <- matrix(0,nrow=N,ncol=2)
		for(i in 1:N)
			B[i,] <- centroidpos(map[[i]])
		# B <- ctrdist(map)
		# b <- B$distance
		b <- rdist(B)
		b <- b/max(b)
		if(is.null(scale))
			scale <- median(b)
		adjmatrix[b < scale] <- 1
		}
	if(type == 3)
		{
		B <- matrix(0,nrow=N,ncol=2)
		for(i in 1:N)
			B[i,] <- centroidpos(map[[i]])

		# B <- ctrdist(map)
		# cent <- B$centroids
		dtr <- delaunayn(B)

		for(i in 1:(N+1))
			{
			w1 <- dtr[i,1]
			w2 <- dtr[i,2]
			w3 <- dtr[i,3]
			adjmatrix[w1,w2] <- 1
			adjmatrix[w1,w3] <- 1
			adjmatrix[w2,w1] <- 1
			adjmatrix[w3,w1] <- 1
			adjmatrix[w2,w3] <- 1
			adjmatrix[w3,w2] <- 1
			}
		}
	if(type == 1)
		{
		for(i in 1:N)
			{
			for(j in 1:N)
				{
				if((i!=j) && (j > i))
					{
					max01x <- max(map[[i]][,1], na.rm=TRUE)
					max02x <- max(map[[j]][,1], na.rm=TRUE)
					min01x <- min(map[[i]][,1], na.rm=TRUE)
					min02x <- min(map[[j]][,1], na.rm=TRUE)
					max01y <- max(map[[i]][,2], na.rm=TRUE)
					max02y <- max(map[[j]][,2], na.rm=TRUE)
					min01y <- min(map[[i]][,2], na.rm=TRUE)
					min02y <- min(map[[j]][,2], na.rm=TRUE)

					if(max01y > max02y)
						{
						check01 <- min01y <= max02y
						check03 <- (max01y > max02y) && (min01y < min02y)
						}
					if(max01y < max02y)
						{
						check01 <- min02y <= max01y
						check03 <- (max02y > max01y) && (min02y < min01y)
						}
					if(max01x > max02x)
						{
						check02 <- min01x <= max02x
						check04 <- (max01x > max02x) && (min01x < min02x)
						}
					if(max01x < max02x)
						{
						check02 <- min02x <= max01x
						check04 <- (max02x > max01x) && (min02x < min01x)
						}

					if(check03 && check04)
						adjmatrix[i,j] <- 1
					if((check01 && check02))
						{
						#adjmatrix[i,j] <- cpoint(map[[i]],map[[j]])

						check1 <- is.element(map[[i]][,1],map[[j]][,1])
						next1 <- map[[i]][check1,]

						check2 <- is.element(map[[j]][,1],map[[i]][,1])
						next2 <- map[[j]][check2,]

						check <- matrix(is.element(next1,next2)*1,ncol=2)

						if(any(rowSums(check) > 1))
							adjmatrix[i,j] <- 1

						# needs to be corrected...
						#else
						#	{
						#	l1 <- nrow(map[[i]])
						#	l2 <- nrow(map[[j]])
				
						#	slope01 <- intc01 <- pcheck01 <- checkem <- rep(0,l1)
						#	slope02 <- intc02 <- pcheck02 <- rep(0,l2)

						#	for(k in 1:l1)
						#		{
						#		if(k < l1)
						#			{
						#			slope01[k] <- (map[[i]][k+1,2]-map[[i]][k,2])/(map[[i]][k+1,1]-map[[i]][k,1])
						#			intc01[k] <- map[[i]][k,2] - slope01[k]*map[[i]][k,1]
						#			pcheck01[k] <- map[[i]][k,2] + map[[i]][k,1]
						#			}
						#		if(k == l1)
						#			{
						#			slope01[k] <- (map[[i]][k,2]-map[[i]][1,2])/(map[[i]][k,1]-map[[i]][1,1])
						#			intc01[k] <- map[[i]][k,2] - slope01[k]*map[[i]][k,1]
						#			pcheck01[k] <- map[[i]][k,2] + map[[i]][k,1]
						#			}
						#		}
						#	for(k in 1:l2)
						#		{
						#		if(k < l2)
						#			{
						#			slope02[k] <- (map[[j]][k+1,2]-map[[j]][k,2])/(map[[j]][k+1,1]-map[[j]][k,1])
						#			intc02[k] <- map[[j]][k,2] - slope02[k]*map[[j]][k,1]
						#			pcheck02[k] <- map[[j]][k,2] + map[[j]][k,1]
						#			}
						#		if(k == l2)
						#			{
						#			slope02[k] <- (map[[j]][k,2]-map[[j]][1,2])/(map[[j]][k,1]-map[[j]][1,1])
						#			intc02[k] <- map[[j]][k,2] - slope02[k]*map[[j]][k,1]
						#			pcheck02[k] <- map[[j]][k,2] + map[[j]][k,1]
						#			}
						#		}

						#	checkem <- slope01 %in% slope02
						#	if(any(checkem))
						#		{
						#		slope01 <- slope01[checkem]
						#		slope02 <- slope02[checkem]
						#		intc01 <- intc01[checkem]
						#		intc02 <- intc02[checkem]
						#		px01 <- map[[i]][,1][checkem]
						#		px02 <- map[[j]][,1][checkem]
						#		py01 <- map[[i]][,2][checkem]
						#		py02 <- map[[j]][,2][checkem]

						#		s1 <- length(slope01)
						#		for(k in 1:s1)
						#			{
						#			ident <- intc01[k] + slope01[k]*px02
						#			if(any(ident %in% py01))
						#				adjmatrix[i,j] <- 1
						#			}
						#		}
						#	if(any(pcheck01 %in% pcheck02))
						#		adjmatrix[i,j] <- 1
						}
					}
				adjmatrix[j,i] <- adjmatrix[i,j]
				}
			}
		}
	diag(adjmatrix) <- 0
	
	return(adjmatrix)
	}
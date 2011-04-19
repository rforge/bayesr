readbndfile <- function(path = path, name = name)
	{
	options(warn = -1)
	data.raw <- scan(path, what = list("", ""), sep = ",", quote = "")
	data.numeric <- list(as.numeric(data.raw[[1]]), as.numeric(data.raw[[2
		]]))
	anzkreise <- sum(is.na(data.numeric[[1]])) - sum(data.raw[[1]] == 
		"is.in")
	cat("Note: map consists of", anzkreise, "polygons\n")
	cat("Reading map ...\n")
	map <- list()
	i <- 1
	for(k in 1:anzkreise) {
		j <- 1
		npoints <- data.numeric[[2]][i]
		if(is.na(data.numeric[[1]][i + 1]) && is.na(data.numeric[[2]][i + 1])) {
			npoints <- npoints + 1
		}
		elem1 <- c()
		elem2 <- c()
		for(j in 1:npoints) {
			elem1 <- c(elem1, data.numeric[[1]][i + j])
			elem2 <- c(elem2, data.numeric[[2]][i + j])
		}
		map[[k]] <- matrix(c(elem1, elem2), ncol = 2)
		names(map)[k] <- substring(data.raw[[1]][i], 2, nchar(data.raw[[
			1]][i]) - 1)
		i <- i + npoints + 1
	}
	if(sum(is.na(as.numeric(names(map)))) == 0) {
		map <- map[order(as.numeric(names(map)))]
		cat("Note: regions sorted by number\n")
	}
	else {
		map <- map[order(names(map))]
		cat("Note: regions sorted by name\n")
	}
	count <- 1
	if(length(map) > 1) {
		for(i in 2:anzkreise) {
			if(names(map)[i] != names(map)[i - 1])
				count <- count + 1
		}
	}
	cat("Note: map consists of", count, "regions\n")
	assign(name, map, pos=1)
	options(warn = 1)
	return(invisible())
	}
centroidtext <- function(polygon, nam, counter, cex,...)
	{
	pos <- centroidpos(polygon)
				
	if(is.null(nam))
		txt <- paste(counter)
	else
		txt <- nam[counter]

	text(pos[1], pos[2], txt, cex=cex,...)
	}
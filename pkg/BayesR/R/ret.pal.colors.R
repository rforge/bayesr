ret.pal.colors <- function(nrc,pal,z)
	{
	surf.colors <- function(x, col = terrain.colors(20)) 
		{
		if(is.na(col[1]))
			{
			hgt <- 0.25 * (z[-ncol(x),-ncol(x)] + z[-1,-ncol(x)] + z[-ncol(x),-1] + z[-1,-1])
        		hgt <- (hgt - min(hgt))/ (max(hgt) - min(hgt))
			colors <- gray(1 - hgt)
			}
		else
			{
  			hgt <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] + x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)])/4
  			colors = col[cut(hgt, breaks = length(col), include.lowest = T)]
			}

  		return(colors)
		}
	if(is.null(nrc))
		nrc <- length(unique(values))
	if(pal == "rainbow")
		colors <- surf.colors(z,col=rainbow(nrc))
	if(pal == "heat")
		colors <- surf.colors(z,col=heat.colors(nrc))
	if(pal == "terrain")
		colors <- surf.colors(z,col=terrain.colors(nrc))
	if(pal == "topo")
		colors <- surf.colors(z,col=topo.colors(nrc))
	if(pal == "cm")
		colors <- surf.colors(z,col=cm.colors(nrc))
	if(pal == "gray")
		colors <- surf.colors(z,col=NA)
	if(!any(pal%in%c("rainbow","heat","terrain","topo","cm","gray")))
		{
		warning("Argument pal specified wrong, set to default!")
		colors <- heat.colors(nrc)
		}
	return(colors)
	}



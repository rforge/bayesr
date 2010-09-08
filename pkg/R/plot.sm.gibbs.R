plot.sm.gibbs <- function(x, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
			  colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, grid = 30, theta = 40, phi = 40, image = FALSE, 
	                  cex = 1, ask = FALSE, border = NULL, nrc = 100, pal = NULL, ...)
	{
	if(diagnostics == FALSE) 
		{
		if(attr(x,"specs")$dim == 1)
			smoothone(x,resid,jit,const,diagnostics,acf,xlab,ylab,main,colored,col,lwdc,lwdconf,cex,...)
		else
			smoothtwo(x,grid,xlab,ylab,zlab,main,theta,phi,image,const,col,colored,resid,cex,border,nrc,pal,...)
		}
	else
		{
		if(diagnostics == 2)
			coefplot(attr(x,"variance.draws"),xlab,ylab,main,attr(x,"specs")$label,acf,TRUE,...)
		else
			coefplot(attr(x,"coef.draws.utr"),xlab,ylab,main,attr(x,"specs")$label,acf,FALSE,...)
		}
	}

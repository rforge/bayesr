plote <- function(model, which = 1, resid = FALSE, main = "(a)", xlab = "x", ylab = "Effect of x", file= "D:\\graphic.eps", 
                width = 4.4, height = 4, fonts = c("serif"),
		    jit = TRUE, const = FALSE, diagnostics = FALSE,  acf = FALSE, zlab = NULL, 
                colored = TRUE, col = c("azure3","azure4"), lwdc = 1, lwdconf = 0, grid = 30, theta = 40, 
                phi = 40, image = FALSE, map = NULL, names = FALSE, values = NULL, range = NULL, 
                pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, dgts = 2, lpos = c(0.2,0),
                cex = 1, ask = FALSE, byplots = FALSE, p3d = FALSE, pcat = NULL, border = NULL, subcex = 1.2, maincex = 1.2, 
		    xlabo = "", ylabo = "", maino = "", ...)
	{
	postscript(file=file,horizontal=FALSE,width=width,height=height,fonts=fonts)

	if(inherits(model,"gibbs"))
		{
		plot(x=model,which=which,resid=resid,xlab=xlabo,ylab=ylabo,main=maino,
           		jit=jit, const=const, diagnostics=diagnostics,  acf=acf, zlab=zlab, 
           		colored=colored, col=col, lwdc=lwdc, lwdconf=lwdconf, grid=grid, theta=theta, 
           		phi=phi, image=image, map=map, names=names, values=values, range=range, 
          	 	pal=pal, legend=legend, scale=scale, nrc=nrc, dgts=dgts, lpos=lpos,
           		cex=cex, ask=ask, byplots=byplots, p3d=p3d, pcat=pcat, border=border, ...)
		}
	else
		{
		plot(x=model,resid=resid,xlab=xlabo,ylab=ylabo,main=maino,
           		jit=jit, const=const, diagnostics=diagnostics,  acf=acf, zlab=zlab, 
           		colored=colored, col=col, lwdc=lwdc, lwdconf=lwdconf, grid=grid, theta=theta, 
           		phi=phi, image=image, map=map, names=names, values=values, range=range, 
          	 	pal=pal, legend=legend, scale=scale, nrc=nrc, dgts=dgts, lpos=lpos,
           		cex=cex, ask=ask, byplots=byplots, p3d=p3d, pcat=pcat, border=border, ...)
		}

	mtext(side=3,text=main,line=1,family=fonts,cex=maincex)
	mtext(side=2,text=ylab,line=2,family=fonts,cex=subcex)
	mtext(side=1,text=xlab,line=2,family=fonts,cex=subcex)

	# legend("topright",legend=c("estimate","true function"),lty=c(1,1),col=c(1,2),cex=0.6)

	dev.off()
	}
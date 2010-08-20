write2dofile <- function(x,file="D:\\tmp",append=FALSE)
	{
	nameofx <- colnames(x)[1]
	tmp <- strsplit(nameofx,":","")[[1]]
	if(length(tmp)>1)
		{
		nameofx <- paste(tmp[1],"_",tmp[2],sep="")
		colnames(x)[1] <- nameofx
		}
	x<-x[!duplicated(x[,1]),]
	x<-cbind(1:nrow(x),x[,1:7],x[,9:10])
	colnames(x)[1]<-"intnr"
	x<-as.data.frame(x)
	write.table(x,row.names=FALSE,quote=FALSE,file=paste(file,"_",nameofx,".res",sep=""))

	doout <- paste(file,".do",sep="")
      write("set scheme s1mono \n\nclear",file=doout,append=append)
	
	write(paste("infile intnr",nameofx,"pmean pqu2p5 pqu10 pmed pqu90 pqu97p5 pcat95 pcat80 using",
                   file=paste(file,"_",nameofx,".res",sep="")),file=doout,append=TRUE)
	write("drop in 1",file=doout,append=TRUE)

	xrange <- c(min(x[,2]),max(x[,2]))
	xstep <- round((xrange[2]-xrange[1])/10,digits=2)
	xrange <- round(xrange,digits=2)

	yrange <- c(min(x[,3]),max(x[,3]))
	ystep <- round((yrange[2]-yrange[1])/5,digits=2)
	yrange <- round(yrange,digits=2)

	write(paste("graph twoway rarea pqu2p5 pqu97p5 ",nameofx,", bcolor(gs13) || rarea pqu10 pqu90 ",nameofx,", bcolor(gs10) || /*
 */ scatter pmean ",nameofx,", c(l) m(i) clpattern(l) clcolor(gs0) /*
 */ ytitle(\"Effect\") subtitle(\"Effect of ",nameofx,"\") xsize() xtitle(\"",nameofx,"\") xlab() ylab() legend(off)",sep=""),
		file=doout,append=TRUE)
	write(paste("graph export ",file,"_",nameofx,".eps, replace",sep=""),file=doout,append=TRUE)

	return(invisible(NULL))
	}
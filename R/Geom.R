Cindex <- function(x,thresh=NULL,connect.method="C") {
   if(!is.null(thresh)) x[x < thresh] <- 0
   NP <- sum(colSums(x==0,na.rm=TRUE),na.rm=TRUE)
   x[x==0] <- NA
   x <- as.im(x)
   x <- connected(x,method=connect.method)
   NC <- max(as.numeric(x$v),na.rm=TRUE)
   return(1 - (NC-1)/(sqrt(NP)+NC))
} # end of 'Cindex' function.

Sindex <- function(x,thresh=NULL,loc=NULL) {
   if(!is.null(thresh)) x[x < thresh] <- 0
   n <- sum(colSums(x>0,na.rm=TRUE),na.rm=TRUE)
   n2 <- sqrt(n)
   if(floor(n2)==n2) Pmin <- 4*sqrt(n)
   else Pmin <- 2*(floor(2*n2)+1)
   if(is.null(loc)) {
	xdim <- dim(x)
	loc <- cbind(rep(1:xdim[1],xdim[2]),rep(1:xdim[2],each=xdim[1]))
   }
   id <- c(x)!=0
   corners <- apply(loc[id,],2,range,finite=TRUE)
   P <- 2*(diff(corners[,1])+1) + 2*(diff(corners[,2])+1)
   return(list(Sindex=Pmin/P, Pmin=Pmin, P=P))
} # end of 'Sindex' function.

Aindex <- function(x,thresh=NULL,dx=1,dy=1) {
   if(is.null(thresh)) thresh <- 1e-8 
   x <- as.im(x)
   x <- solutionset(x>=thresh)
   ch <- convexhull(x)
   A <- sum(colSums(x$m,na.rm=TRUE),na.rm=TRUE)*dx*dy
   Aconvex <- area.owin(ch)*dx*dy
   Aindex <- A/Aconvex
   return(list(Aindex=Aindex,A=A,Aconvex=Aconvex,dx=dx,dy=dy))
} # end of 'Aindex' function.

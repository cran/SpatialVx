GeoBoxPlot <- function(x, areas, ...) {
   ##
   ## Function to create the geographic boxplots of Willmott et al. (2007).
   ##
   ## Arguments:
   ##
   ## 'x' numeric vector of data to be box-plotted.
   ## 'areas' numeric vector of same length as 'x' giving the areas corresponding to each value of 'x'.
   ## '...' optional arguments to the 'boxplot' function.  The argument 'plot' is not allowed.
   ##
   ## Value: A plot is produced.  The values of the statistics going into the
   ##	boxplot's are returned invisibly.
   ##
   if( any( is.na(x)) | any(is.na( areas))) warning("GeoBoxPlot: missing values are not handled.")
   out <- boxplot( x, plot=FALSE, ...)
   ox <- order( x)
   a <- areas[ ox]
   a.frac <- a/sum( a)
   a2 <- cumsum( a.frac)
   x <- x[ox]
   n <- length( x)
   out$stats[,1] <- c( min(x), x[ min( (1:n)[ a2>0.25])], x[ min( (1:n)[ a2>0.5])], x[ min( (1:n)[ a2>0.75])], max( x))
   out <- bxp(out, ...)
   invisible( out)
}

KernelGradFUN <- function(x,ktype="LoG", nx=10, ny=12, sigma=1) return(kernel2dsmooth(x,kernel.type=ktype, nx=nx, ny=ny, sigma=sigma))

S1 <- function(X,Y,gradFUN="KernelGradFUN", ...) {
   Xgrad <- do.call(gradFUN, c(list(X),list(...)))
   Ygrad <- do.call(gradFUN, c(list(Y),list(...)))
   denom <- sum(colSums(pmax(abs(Xgrad),abs(Ygrad),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE)
   numer <- sum(colSums(abs(Ygrad - Xgrad),na.rm=TRUE),na.rm=TRUE)
   return(100*numer/denom)
}

ACC <- function(X,Y,Xclim=NULL,Yclim=NULL) {
   if(!is.null(Xclim)) X <- X - Xclim
   if(!is.null(Yclim)) Y <- Y - Yclim
   return(cor(c(X),c(Y)))
}

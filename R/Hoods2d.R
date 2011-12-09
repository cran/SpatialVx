hoods2d <-
function( obj, which.methods = c("mincvr", "multi.event", "fuzzy", "joint", "fss", "pragmatic"), verbose = FALSE) {
   if( verbose) begin.time <- Sys.time() 
   thresholds <- obj$thresholds
   q <- dim( thresholds)[1]
   levels <- obj$levels
   l <- length( levels)
   Y <- get( obj$Fcst.name)
   X  <- get( obj$Vx.name)
   xdim <- obj$xdim
   binmat <- matrix(0, xdim[1], xdim[2])
   outmat <- matrix( NA, l, q)
   out <- hoods2dSetUpLists( which.methods=which.methods, mat=outmat)
   out$prep.object <- as.character( substitute( obj))
   bigN <- obj$Nxy
   sub <- obj$subset
   if( verbose) cat("Looping through thresholds.\n")
   for( threshold in 1:q) {
      	if( any( c("mincvr", "mincvr", "multi.event", "fuzzy", "joint", "fss", "pragmatic") %in% which.methods)) {
	   if( verbose) cat("\n", "Setting up binary objects for thresholds = ", thresholds[ threshold,], "\n")
	   Ix <- Iy <- binmat
	   Ix[ X >= thresholds[threshold,2]] <- 1
	   Iy[ Y >= thresholds[threshold,1]] <- 1
        } # end of if find 'Ix' and 'Iy' stmt.
	if( "fss" %in% which.methods) {
	   if( is.null( sub)) f0 <- mean( Ix, na.rm=TRUE)
	   else f0 <- mean( c(Ix)[ subset], na.rm=TRUE)
	   if( threshold==1) {
		out$fss$fss.uniform <- 0.5 + f0/2 
		out$fss$fss.random  <- f0
	   } else {
		out$fss$fss.uniform <- c( out$fss$fss.uniform, 0.5 + f0/2)
		out$fss$fss.random  <- c( out$fss$fss.random, f0)
	   } # end of if else 'threshold' is 1 stmts.
	} # end of if 'fss' method.
	if( verbose) cat( "Looping through levels.\n")
        for( level in 1:l) {
	   if( verbose) cat("Neighborhood length = ", levels[ level], "\n")
	   # levelW <- kernel2dsmooth( X, kernel.type="boxcar", n=levels[level], xdim=xdim, setup=TRUE)
	   levelW <- do.call( obj$smooth.fun, c(list( x=X, lambda=levels[level], W=NULL, setup=TRUE), obj$smooth.params))
	   if( any( c( "mincvr", "multi.event", "fuzzy", "joint", "pragmatic", "fss") %in% which.methods)) {
		if( any( c( "mincvr", "multi.event", "fuzzy", "joint", "fss") %in% which.methods)) 
			sPx <- do.call( obj$smooth.fun, c(list( x=Ix, lambda=levels[level], W=levelW), obj$smooth.params))
			# sPx <- kernel2dsmooth( Ix, kernel.type="boxcar", n=levels[level], W=levelW, xdim=xdim)
		sPy <- do.call( obj$smooth.fun, c(list(x=Iy, lambda=levels[level], W=levelW), obj$smooth.params))
		# sPy <- kernel2dsmooth( Iy, kernel.type="boxcar", n=levels[level], W=levelW, xdim=xdim)
	   } # end of if any 'mincvr', 'multi', 'fuzzy', 'joint' or 'pragmatic' stmts.
	   if( any( c( "mincvr", "multi.event") %in% which.methods)) {
                sIx <- sIy <- binmat
                sIx[ sPx >= obj$Pe[level]] <- 1
		sIy[ sPy >= obj$Pe[level]] <- 1
		if( "mincvr" %in% which.methods) {
		   tmp <- MinCvg2dfun( sIy=sIy, sIx=sIx, subset=sub)
		   out$mincvr$pod[level,threshold] <- tmp$pod
		   out$mincvr$far[level,threshold] <- tmp$far
		   out$mincvr$ets[level,threshold] <- tmp$ets
		} # end of 'mincvr' part.
		if( "multi.event" %in% which.methods) {
		   tmp <- multicon2dfun( sIy=sIy, Ix=Ix, subset=sub)
		   out$multi.event$pod[level,threshold] <- tmp$pod
		   out$multi.event$f[level,threshold] <- tmp$f
		   out$multi.event$hk[level,threshold] <- tmp$hk
		} # end of 'multi' stmts.
           } # end of 'mincvr'/'multi' stmts.
	   if( any( c("fuzzy", "joint") %in% which.methods)) {
		tmp <- fuzzyjoint2dfun( sPy=sPy, sPx=sPx, subset=sub)
		if( "fuzzy" %in% which.methods) {
		   out$fuzzy$pod[level,threshold] <- tmp$fuzzy$pod
		   out$fuzzy$far[level,threshold] <- tmp$fuzzy$far
		   out$fuzzy$ets[level,threshold] <- tmp$fuzzy$ets
		} # end of if 'fuzzy' stmt.
		if( "joint" %in% which.methods) {
		   out$joint$pod[level,threshold] <- tmp$joint$pod
		   out$joint$far[level,threshold] <- tmp$joint$far
		   out$joint$ets[level,threshold] <- tmp$joint$ets
		} # end of if 'joint' methods.
	   } # end of if 'fuzzy/joint' stmts.
	   if( "fss" %in% which.methods) out$fss$values[level,threshold] <- fss2dfun(sPy=sPy, sPx=sPx, subset=sub)
	   if( "pragmatic" %in% which.methods) {
		tmp <- pragmatic2dfun( sPy=sPy, Ix=Ix, subset=sub)
		out$pragmatic$bs[level,threshold] <- tmp$bs
		out$pragmatic$bss[level,threshold] <- tmp$bss
	   }
	} # end of for 'level' loop.
   } # end of for 'threshold' loop.
   if( verbose) print( Sys.time() - begin.time)
   class( out) <- "hoods2d"
   return( out)
} # end of 'hoods2d' function.

hoods2dPrep <-
function( Fcst.name, Vx.name, thresholds=NULL, Pe=NULL, levels=NULL, max.n=NULL, subset=NULL,
                                loc=NULL, qs=NULL, field.type="", units=NULL, smooth.fun="hoods2dsmooth", smooth.params=NULL, mapit=FALSE) {
   out <- list()
   out$Fcst.name <- Fcst.name
   out$Vx.name <- Vx.name
   Fcst <- get( Fcst.name)
   Obs <- get( Vx.name)
   dimF <- dim( Fcst)
   out$xdim <- dimF
   out$map <- mapit
   if( any(dimF != dim( Obs))) stop("fss2d: dim of Obs must be the same as dim of Fcst.")
   Nxy <- prod( dimF)
   out$Nxy <- Nxy
   if( is.null( levels)) {
        if( is.null( max.n)) max.n <- 2*max(dimF)-1
        else {
           if( max.n %%2 == 0) max.n <- max.n-1
           if( max.n > 2*max( dimF)-1) stop(paste("fss2d: max.n must be less than 2N-1, where N is ", max(dimF), sep=""))
        } # end of if else 'max.n' stmts.
        if( max.n < 1) stop("fss2d: max.n must be a positive integer.")
        levels <- seq(1,max.n,2)
   } # end of if no 'levels' given.
   out$levels <- levels
   out$max.n <- max.n
   if( is.null( thresholds)) {
	thresholds <- cbind( quantile( c(Fcst), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95)),
                                quantile(c(Obs), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95)))
	qs <- as.character( c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95))
   } else if( !is.matrix( thresholds)) thresholds <- cbind( thresholds, thresholds)
   out$thresholds <- thresholds
   out$smooth.fun <- smooth.fun
   out$smooth.params <- smooth.params
   if( is.null( Pe)) out$Pe <- 1/(levels^2)
   if( length( Pe) == 1) out$Pe <- rep( Pe, length( levels))
   else  out$Pe <- Pe
   out$subset <- subset
   out$loc <- loc
   out$qs <- qs
   out$units <- units
   class( out) <- "hoods2dPrep"
   return( out)
} # end of 'hoods2dPrep' function.

plot.hoods2dPrep <- function(x, ...) {
   X <- get(x$Vx.name)
   Y <- get(x$Fcst.name)
   zl <- range(c(c(X),c(Y)),finite=TRUE)
   par(mfrow=c(1,2),mar=c(4.1,2.1,4.1,2.1))
   image(X, col=c("gray", tim.colors(64)), axes=FALSE, zlim=zl, main=paste(x$Vx.name, "\n", x$field.type, " (", x$units, ")", sep=""))
   if(!is.null(x$thresholds)) contour(X, levels=x$thresholds[,2], add=TRUE, lwd=0.5, lty=2, col="gray")
   if(x$map) {
	par(usr=c(range(x$loc[,1],finite=TRUE),range(x$loc[,2], finite=TRUE)))
	map(add=TRUE)
	map(add=TRUE,database="state",lty=2)
   }
   image(Y, col=c("gray", tim.colors(64)), axes=FALSE, zlim=zl, main=paste(x$Fcst.name, "\n", x$field.type, " (", x$units, ")", sep=""))
   if(!is.null(x$thresholds)) contour(Y, levels=x$thresholds[,1], add=TRUE, lwd=0.5, lty=2, col="gray")
   if(x$map) {
        par(usr=c(range(x$loc[,1],finite=TRUE),range(x$loc[,2], finite=TRUE)))
        map(add=TRUE)
	map(add=TRUE,database="state",lty=2)
   }
   image.plot(Y, col=c("gray", tim.colors(64)), zlim=zl, legend.only=TRUE, legend.mar=2.1, horizontal=TRUE)
   # X <- as.im(X)
   # hist(X, col="darkblue", main="", probability=TRUE, ylim=c(0,1), xlab="Grid value")
   # Y <- as.im(Y)
   # hist(Y, col="darkblue", main="", probability=TRUE, ylim=c(0,1), xlab="Grid value")
   invisible()
} # end of 'plot.hoods2dPrep' function.

hist.hoods2dPrep <- function(x, ...) {
   X <- get(x$Vx.name)
   Y <- get(x$Fcst.name)
   u <- x$thresholds
   udim <- dim(u)
   udim[1] <- udim[1]+1
   par(mfrow=udim)
   for(i in 1:udim[1]) {
	if(i==1) {
	   m1 <- paste(x$Vx.name, "\n", "All ", x$field.type, " values", sep="")
	   m2 <- paste(x$Fcst.name, "\n", "All ", x$field.type, " values", sep="")
	   tmpX <- c(X)
	   tmpY <- c(Y)
	} else {
	   m1 <- paste("Only values >= ", u[i-1,2], " ", x$units, sep="")
	   m2 <- paste("Only values >= ", u[i-1,1], " ", x$units, sep="")
	   tmpX <- c(X[X>=u[i-1,2]])
	   tmpY <- c(Y[Y>=u[i-1,1]])
	}
	if(i==udim[1]) x1 <- paste(x$field.type, " (", x$units, ")", sep="")
	else x1 <- ""
	hist(tmpX, main=m1, xlab=x1, ...)
	hist(tmpY, main=m2, xlab=x1, ...)
   } # end of for 'i' loop.
   invisible()
} # end of 'hist.hoods2dPrep' function.

plot.hoods2d <-
function(x, ...) {
   a <- get( x$prep.object)
   mets <- names( x)
   if( "mincvr" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$mincvr$ets, args=a, main=paste("Min. Coverage: Gilbert Skill Score (GSS)", sep="")))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$mincvr$far, args=a, main="Min. Coverage: False Alarm Ratio"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$mincvr$pod, args=a, main="Min. Coverage: Hit Rate"))
   } # end of if 'mincvr' stmts.
   if( "multi.event" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$multi.event$pod, args=a, main=paste("Multi-Event Contingency Table: Hit Rate", sep="")))
	a$ylab <- "False Alarm Rate"
        try( hoods2dPlot( x$multi.event$f, args=a, main="Multi-Event Contingency Table: False Alarm Rate"))
	a$ylab <- "Hanssen-Kuipers Score"
	try(hoods2dPlot( x$multi.event$hk, args=a, main="Multi-Event Contingency Table: HK"))
   } # end of if 'multi.event' stmts.
   if( "fuzzy" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$fuzzy$ets, args=a, main="Fuzzy Logic: Gilbert Skill Score (GSS)"))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$fuzzy$far, args=a, main="Fuzzy Logic: False Alarm Ratio (FAR)"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$fuzzy$pod, args=a, main="Fuzzy Logic: Hit Rate"))
   } # end of if 'fuzzy' stmts.
   if( "joint" %in% mets) {
	par( mfrow=c(4,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$joint$ets, args=a, main="Joint Probability: Gilbert Skill Score (GSS)"))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$joint$far, args=a, main="Joint Probability: False Alarm Ratio (FAR)"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$joint$pod, args=a, main="Joint Probability: Hit Rate"))
   } # end of if 'joint' stmts.
   if( "fss" %in% mets) {
	look <- a
	look$values <- x$fss$values
	look$fss.random <- x$fss$fss.random
	look$fss.uniform <- x$fss$fss.uniform
	try( fss2dPlot( look, main="Fractions Skill Score (FSS)"))
   } # end of if 'fss' stmts.
   if( "pragmatic" %in% mets) {
	par( mfrow=c(2,2), bg="beige")
	a$ylab <- "Brier Score"
	try( hoods2dPlot( x$pragmatic$bs, args=a, main="Pragmatic: Brier Score (BS)"))
	a$ylab <- "Brier Skill Score"
	try( hoods2dPlot( x$pragmatic$bss, args=a, main="Pragmatic: Brier Skill Score (BSS)"))
   } # end of if 'pragmatic' stmts.
   invisible()
} # end of 'plot.hoods2d' function.

fss2dPlot <-
function(x, ...) {
   ##
   ## Function to make useful plots for the output from 'fss2d'.
   ##
   ## Arguments:
   ##
   ## 'x' a list object of class "fss2d".
   ## '...' optional arguments to the 'image' and 'image.plot' functions.
   ##
   ## Value: No value is returned, but a two-panel plot is created: the first is an image plot (quilt plot) showing the FSS scores
   ##	for each threshold/level combination.  The second is a line plot with a different colored line for each threshold.  The second
   ##	plot also shows the FSS for a random forecast (dotted line) and based on a uniform distribution (dotted line).
   ##
   # image plot
   odim <- dim( x$values)
   q <- odim[2]
   l <- odim[1]
   if( is.null( odim)) stop("fss2d: values must be a matrix.")
   if( q == 1 & l == 1) stop("fss2d: values must be a matrix with at least one dimension > 1.")
   par( mfrow=c(1,2), mar=c(5.1, 4.1, 3.1, 4.1), bg="beige")
   # image plot
   if( is.null( x$qs) & all( x$thresholds[,1] == x$thresholds[,2])) {
        a1.labels <- as.character( x$thresholds[,1])
        xlb <- paste("Threshold (", x$units, ")", sep="")
   } else if( !is.null( x$qs)) {
        a1.labels <- x$qs
        xlb <- "Threshold (quantile)"
   } else {
        a1.labels <- as.character( 1:dim( x$thresholds)[1])
        xlb <- "Threshold index number"
   }
   image( t( x$values), xaxt="n", yaxt="n", xlab=xlb, ylab="Neighborhood size (grid squares)", col=c("grey", tim.colors(64)), ...)
   text( x=seq(0,1,,q)[ rep(1:q,l)], y=seq(0,1,,l)[ rep(1:l,each=q)], labels=round( t( x$values), digits=2))
   axis( 1, at=seq(0,1,,q), labels=a1.labels)
   if( any( x$thresholds[,1] != x$thresholds[,2])) warning("fss2dPlot: thresholds differ for the two fields, labels only for the first")
   axis( 2, at=seq(0,1,,length(x$levels)), labels=x$levels)
   # lines( seq( 0, 1,,q), x$fss.random, lty=1, col="darkorange")
   # lines( seq( 0, 1,,q), x$fss.uniform, lty=2, col="darkred")
   image.plot( x$values, legend.only=TRUE, col=c("grey", tim.colors(64)), ...)

   # line plot
   look <- x$values
   matplot( look, ylim = c(0,1), ylab = "FSS", xlab = "Neighborhood size (grid squares)", type = "l", lty = 1, axes = FALSE , lwd = 2)
   abline(h=c(x$fss.uniform), col=1:q, lty=2)
   abline(h=c(x$fss.random), col=1:q, lty=3)
   axis(2)
   box()
   axis(1, at = 1:l, labels = x$levels)
   grid()
   legend("topleft", legend = a1.labels, col = 1:q, title = xlb, inset = 0.02, lwd = 2 )
   invisible()
} # end of 'fss2dPlot' function.

Fourier2d <-
function(x, bigdim=NULL, kdim=NULL) {
   ##
   ## Function to compute the Fast Fourier Transform (FFT) of a field 'x'.
   ##
   ## Arguments:
   ##
   ## 'x' 'n X m' numeric matrix.
   ## 'bigdim' numeric vector of length 2 giving the optimized dimension for the FFT.  If NULL, then 'kdim' must be supplied,
   ##	and it will be determined.
   ## 'kdim' dimension of the kernel matrix before expansion.  Only used if 'bigdim' is not supplied.
   ##
   ## Details: This function creates a matrix of zeros with twice the dimension of the input field, 'x'.
   ##	'x' is then put in the upper left corner of the matrix of zeros, missing values are set to zero,
   ##	and the FFT is calculated for the large matrix.  This is what is returned.  The output can then be
   ##	used by the 'kernel2dsmooth' function in order to possibly save an FFT in computing time at some point.
   ##
   ## See also: 'fft', 'kernel2dsmooth', 'hoods2d'
   ##
   ## Value: '2n X 2m' numeric matrix giving the Fourier transformed field.
   ##
   xdim <- dim( x)
   if( is.null( bigdim)) {
	if(is.null(kdim)) stop("Fourier2d: one of bigdim or kdim must be supplied.")
	bigdim <- xdim + kdim -1
	if( bigdim[1] <= 1024) bigdim[1] <- 2^ceiling(log2(bigdim[1]))
        else bigdim[1] <- ceiling( bigdim[1]/512)*512
        if( bigdim[2] <= 1024) bigdim[2] <- 2^ceiling(log2(bigdim[2]))
        else bigdim[2] <- ceiling( bigdim[2]/512)*512
   }
   out <- matrix( 0, bigdim[1], bigdim[2])
   out[1:xdim[1],1:xdim[2]] <- x
   out[ is.na( out)] <- 0
   return( fft( out))
} # end of 'Fourier2d' function.

vxstats <-
function(Fcst, Obs, which.stats=c("bias", "ts", "ets", "pod", "far", "f", "hk", "bcts", "bcets", "mse"), subset=NULL) {
   ##
   ## Function to calculate various traditional verification statistics for a gridded
   ## verification set.
   ##
   ## Arguments:
   ##
   ## 'Fcst', 'Obs' 'k X m' logical or numeric matrices of forecast and observed values, resp.
   ## 'which.stats' character vector telling which verification statistics should be computed.
   ## 'subset' numeric vector indicating a subset of points over which to calculate the statistics.  If NULL, then the entire
   ##	fields are used.
   ##
   ## Details:
   ##
   ## The possible statistics that can be computed, as determined by 'which.stats' are:
   ##
   ##	"bias" the number of forecast events divided by the number of observed events (sometimes called frequency bias).
   ##	"ts" threat score, given by hits/(hits + misses + false alarms)
   ##	"ets" equitable threat score, given by (hits - hits.random)/(hits + misses + false alarms - hits.random), where
   ##		'hits.random' is the number of observed events times the number of forecast events divided by the total
   ##		number of forecasts.
   ##	"pod" probability of detecting an observed event (aka, hit rate).  It is given by hits/(hits + misses).
   ##	"far" false alarm ratio, given by (false alarms)/(hits + false alarms).
   ##	"f" false alarm rate (aka probability of false detection) is given by (false alarms)/(correct rejections + false alarms).
   ##	"hk" Hansen-Kuipers Score is given by the difference between the hit rate ("pod") and the false alarm rate ("f").
   ##	"mse" mean square error (not a contingency table statistic, but can be used with binary fields).  This is the only
   ##		statistic that can be calculated here that does not require binary fields for 'Fcst' and 'Obs'.
   ##
   ## Warnings: It is up to the user to provide the appropriate type of fields for the given statistics
   ##	to be computed.  For example, they must be binary for all types of 'which.stats' except "mse".
   ##
   ## Value: list object with component names the same as 'which.stats' giving the single numeric
   ##	value of each statistic for the two fields.
   ##
   out <- list()
   xdim <- dim( Fcst)
   if( "mse" %in% which.stats) {
        if( is.null( subset)) {
         Nxy <- sum( colSums( !is.na( Fcst) & !is.na( Obs), na.rm=TRUE), na.rm=TRUE)
         out$mse <- sum( colSums( (Fcst - Obs)^2, na.rm=TRUE), na.rm=TRUE)/Nxy
        } else {
           out$mse <- mean( (c(Fcst)[ subset] - c(Obs)[subset])^2, na.rm=TRUE)
        }
   } # end of if do MSE.
   if( any( is.element(c("bias", "ts", "ets", "pod", "far", "f", "hk", "bcts", "bcets"), which.stats))) {
	if( !is.logical( Fcst)) Fcst <- as.logical( Fcst)
	if( !is.logical( Obs)) Obs <- as.logical( Obs)
	if( is.null( subset)) {
	   hits <- sum( colSums( matrix( Fcst & Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   miss <- sum( colSums( matrix( !Fcst & Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   fa   <- sum( colSums( matrix( Fcst & !Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   if( any( c("ets", "f", "hk") %in% which.stats)) cn <- sum( colSums( matrix( !Fcst & !Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	} else {
	   hits <- sum( c(Fcst)[subset] & c(Obs)[subset], na.rm=TRUE)
           miss <- sum( !c(Fcst)[subset] & c(Obs)[subset], na.rm=TRUE)
           fa   <- sum( c(Fcst)[subset] & !c(Obs)[subset], na.rm=TRUE)
           if( any( c("ets", "f", "hk") %in% which.stats)) cn <- sum( !c(Fcst)[subset] & !c(Obs)[subset], na.rm=TRUE)
	}
	if( "bias" %in% which.stats) {
	   if( (hits + fa == 0) & (hits + miss == 0)) out$bias <- 1
	   else if( hits + miss == 0) out$bias <- (hits + fa)/(1e-8)
	   else out$bias <- (hits + fa)/(hits + miss)
 	}
	if( "ts" %in% which.stats) {
	   if( hits == 0) out$ts <- 0
	   else out$ts <- hits/(hits + miss + fa)
	}
	if( "ets" %in% which.stats) {
	   if( (hits + miss == 0) | (hits + fa == 0)) hits.random <- 0
	   else hits.random <- (hits + miss)*(hits + fa)/(hits + miss + fa + cn)
	   if( hits + miss + fa == 0) out$ets <- 0
	   else out$ets <- (hits - hits.random)/(hits + miss + fa - hits.random)
	}
	if( any( c("pod", "hk") %in% which.stats)) {
	   if( hits + miss == 0) pod <- 0
	   else pod <- hits/(hits + miss)
	   if( "pod" %in% which.stats) out$pod <- pod
	} 
	if( "far" %in% which.stats) {
	   if( hits + fa == 0) out$far <- 0
	   else out$far <- fa/(hits + fa)
	}
	if( any( c("f", "hk") %in% which.stats)) {
	   if( cn + fa == 0) f <- 0
	   else f <- fa/(cn + fa)
	   if( "f" %in% which.stats) out$f <- f
	   if( "hk" %in% which.stats) out$hk <- pod - f
	}
	if(any(is.element(c("bcts", "bcets"),which.stats))) {
	   nF <- hits + fa
	   nO <- hits + miss
	   lf <- log(nO/miss)
	   Ha <- nO - (fa/lf)*W((nO/fa)*lf)
	   if(is.element("bcts",which.stats)) out$bcts <- Ha/(2*nO-hits)
	   if(is.element("bcets",which.stats)) out$bcets <- (Ha - (nO^2)/(hits+miss+fa+cn))/(2*nO-Ha-(nO^2)/(hits+miss+fa+cn))
	}
    } # end of if any contingency table scores stmts.
   class( out) <- "vxstats"
   return( out)
} # end of 'vxstats' function.

hoods2dSetUpLists <-
function( which.methods, mat) {
   out <- list()
   if( "mincvr" %in% which.methods) {
      out$mincvr <- list()
      out$mincvr$pod <- out$mincvr$far <- out$mincvr$ets <- mat
   } # end of if calculate minimum coverage stmts.
   if( "multi.event" %in% which.methods) {
      out$multi.event <- list()
      out$multi.event$pod <- out$multi.event$f <- out$multi.event$hk <- mat
   } # end of if calculate "multi.event" stmts.
   if( "fuzzy" %in% which.methods) {
      out$fuzzy <- list()
      out$fuzzy$pod <- out$fuzzy$far <- out$fuzzy$ets <- mat
   } # end of if do fuzzy logic method.
   if( "joint" %in% which.methods) {
      out$joint <- list()
      out$joint$pod <- out$joint$far <- out$joint$ets <- mat
   } # end of if do joint prob method.
   if( "fss" %in% which.methods) {
      out$fss <- list()
      out$fss$values <- mat
      out$fss$fss.random <- out$fss$fss.uniform <- numeric(0)
   } # end of if do FSS method.
   if( "pragmatic" %in% which.methods) {
      out$pragmatic <- list()
      out$pragmatic$bs <- out$pragmatic$bss <- mat
   } # end of do pragmatic method.
   return( out)
} # end of 'hoods2dSetUpLists' function.

MinCvg2dfun <-
function( sIy, sIx, subset=NULL) {
   ##
   ## Function to calculate the minimum coverage neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy', 'sIx' (optional) 'k X m' binary forecast and observed matrices, resp.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##	all of the points are used.
   ##
   ## Value: list object with components: "pod", "far" and "ets".
   ##
   out <- vxstats( sIy, sIx, which.stats=c("pod", "far", "ets"), subset=subset)
   return( out)
} # end of 'MinCvg2dfun' function.

multicon2dfun <-
function(sIy, Ix, subset=NULL) {
   ##
   ## Function to calculate the multi-event contingency table neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy' 'k X m' binary smoothed forecast matrix.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components "pod" and "f" and "hk".
   ##
   out <- vxstats( sIy, Ix, which.stats=c("pod", "f", "hk"), subset=subset)
   return( out)
} # end of 'multicon2dfun' function.

fuzzyjoint2dfun <-
function( sPy, sPx, subset=NULL) {
   ##
   ## Function to calculate the fuzzy logic and joint probability neighborhood methods.
   ##
   ## Arguments:
   ##
   ## 'sPy', 'sPx' smoothed 'k X m' forecast and observed matrices, resp., output from 'kernel2dsmooth'.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components: "fuzzy" and "joint", each of which are themselves
   ##	list objects each with components: "pod", "far" and "ets".
   ##
   out <- list()
   vxfun <- function(n11, n01, n10, n00) {
	pod <- n11/(n11 + n01)
	far <- n10/(n11 + n10)
	hits.random <- (n11 + n01)*(n11 + n10)/(n11 + n01 + n10 + n00)
	ets <- (n11 - hits.random)/(n11 + n01 + n10 - hits.random)
	return( list(pod=pod, far=far, ets=ets))
   } # end of internal 'vxfun' function.
   if( is.null( subset)) {
      hits <- sum( colSums( sPx < sPy, na.rm=TRUE), na.rm=TRUE)
      miss <- sum( colSums( sPx < (1-sPy), na.rm=TRUE), na.rm=TRUE)
      fa   <- sum( colSums( (1-sPx) < sPy, na.rm=TRUE), na.rm=TRUE)
      cn   <- sum( colSums( (1-sPx) < (1-sPy), na.rm=TRUE), na.rm=TRUE)
      out$fuzzy <- vxfun( hits, miss, fa, cn)
      hits <- sum( colSums( sPx*sPy, na.rm=TRUE), na.rm=TRUE)
      miss <- sum( colSums( sPx*(1-sPy), na.rm=TRUE), na.rm=TRUE)
      fa   <- sum( colSums( (1-sPx)*sPy, na.rm=TRUE), na.rm=TRUE)
      cn   <- sum( colSums( (1-sPx)*(1-sPy), na.rm=TRUE), na.rm=TRUE)
      out$joint <- vxfun( hits, miss, fa, cn)
   } else {
      hits <- sum( c(sPx)[subset] < c(sPy)[subset], na.rm=TRUE)
      miss <- sum( c( sPx)[subset] < (1-c(sPy)[subset]), na.rm=TRUE)
      fa   <- sum( (1-c(sPx)[subset]) < c(sPy)[subset], na.rm=TRUE)
      cn   <- sum( (1-c(sPx)[subset]) < (1-c(sPy)[subset]), na.rm=TRUE)
      out$fuzzy <- vxfun( hits, miss, fa, cn)
      hits <- sum( (c(sPx)[subset])*(c(sPy)[subset]), na.rm=TRUE)
      miss <- sum( (c( sPx)[subset])*(1-c(sPy)[subset]), na.rm=TRUE) 
      fa   <- sum( (1-c(sPx)[subset])*c(sPy)[subset], na.rm=TRUE)
      cn   <- sum( (1-c(sPx)[subset])*(1-c(sPy)[subset]), na.rm=TRUE)
      out$joint <- vxfun( hits, miss, fa, cn)
   } # end of if else 'subset' stmt.
   return( out)
} # end of 'fuzzyjoint2dfun' function.

pragmatic2dfun <-
function( sPy, Ix, mIx=NULL, subset=NULL) {
   ##
   ## Function to calculate the pragmatic neighborhood statistics.
   ##
   ## Argmunets:
   ##
   ## 'sPy' 'k X m' matrix of smoothed forecast event frequencies.
   ## 'Ix' 'k X m' binary matrix of observed events.
   ## 'mIx' single numeric giving the base rate.  If NULL, it will be computed here.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components 'bs' and 'bss' giving the Brier and Brier Skill Scores.
   ##
   bs <- vxstats( sPy, Ix, which.stats="mse", subset=subset)$mse
   if( is.null( subset)) {
      if( is.null( mIx)) mIx <- mean( Ix, na.rm=TRUE)
      denom <- sum( colSums( (mIx - Ix)^2, na.rm=TRUE), na.rm=TRUE)/sum( !is.na( Ix), na.rm=TRUE)
   } else {
      if( is.null( mIx)) mIx <- mean( c(Ix)[subset], na.rm=TRUE)
      denom <- mean( (mIx - c(Ix)[subset])^2, na.rm=TRUE)
   }
   bss <- 1 - bs/denom
   return( list( bs=bs, bss=bss))
} # end of 'pragmatic2dfun' function.

upscale2d <- function(object, thresholds=NULL, verbose=FALSE) {
   out <- list()
   if( is.null( thresholds)) thresholds <- object$thresholds
   if( !is.matrix( thresholds)) thresholds <- cbind( thresholds, thresholds)
   q <- dim( thresholds)[1]
   levels <- object$levels
   l <- length( levels)
   out$l <- l
   out$q <- q
   out$thresholds <- thresholds
   out$levels <- levels
   out$rmse <- numeric(l)+NA
   out$bias <- out$ts <- out$ets <- matrix( NA, l, q)
   X <- get( object$Vx.name)
   Y <- get( object$Fcst.name)
   xdim <- dim( X)
   sub <- object$subset
   for( level in 1:l) {
      levelW <- kernel2dsmooth( Y, kernel.type="boxcar", n=levels[level], setup=TRUE)
      sYy <- kernel2dsmooth( Y, W=levelW, xdim=xdim)
      sYx <- kernel2dsmooth( X, W=levelW, xdim=xdim)
      out$rmse[level] <- upscale2dfun(sYy=sYy, sYx=sYx, which.stats="rmse", subset=sub)$rmse
      if( verbose) cat("\n", "RMSE for neighborhood length = ", levels[level], " is ", out$rmse[level], "\n")
      for( threshold in 1:q) {
         tmp <- upscale2dfun(sYy=sYy, sYx=sYx, threshold=c(thresholds[threshold,]), which.stats=c("bias", "ts", "ets"), subset=sub)
         out$bias[level,threshold] <- tmp$bias
         out$ts[level,threshold] <- tmp$ts
         out$ets[level,threshold] <- tmp$ets
         if( verbose) {
	   cat("Bias for thresholds: ", thresholds[threshold], " is ", tmp$bias, "\n")
	   cat("Threat score is ", tmp$ts, "\n")
	   cat("GSS is ", tmp$ets, "\n")
	 }
       } # end of for 'threshold' loop.
   } # end of for 'level' loop.
   class( out) <- "upscale2d"
   return( out)
} # end of 'upscale2d' function.

upscale2dfun <-
function( sYy, sYx, threshold=NULL, which.stats=c("rmse", "bias", "ts", "ets"), subset=NULL) {
   ##
   ## Function to calculate the upscaling neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy', 'sIx' 'k X m' upscaled forecast and observed matrices (e.g., as output from 'kernel2dsmooth'), resp.
   ## 'threshold' (optional) numeric of length 2 giving the value over which to calculate the statistics.  If NULL,
   ##	only RMSE is calculated.  The forecast threshold is first, and the observbed one second.
   ## 'which.stats' character vector telling which statistics to compute.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components "rmse", and if 'threshold' is not NULL, "bias", "gss", and "ts".
   ##
   out <- list()
   if( "rmse" %in% which.stats) out$rmse <- sqrt( vxstats( sYy, sYx, which.stats="mse", subset=subset)$mse)
   if( !is.null( threshold)) {
	xdim <- dim( sYy)
	binmat <- matrix(0, xdim[1], xdim[2]) 
	sIx <- sIy <- binmat
	sIx[ sYx >= threshold[ 2]] <- 1
	sIy[ sIy >= threshold[ 1]] <- 1
	tmp <- vxstats( sIy, sIx, which.stats=c("bias", "ts", "ets")[ c("bias", "ts", "ets") %in% which.stats], subset=subset)
	if( "bias" %in% which.stats) out$bias <- tmp$bias
	if( "ts" %in% which.stats) out$ts <- tmp$ts
	if( "ets" %in% which.stats) out$ets <- tmp$ets
   } # end of if 'threshold' stmts.
   return( out)
} # end of 'upscale2dfun'.

hoods2dPlot <-
function(x, args, ...) {
   require( fields)
   odim <- dim( x)
   q <- odim[2]
   l <- odim[1]
   if( is.null( odim)) stop("hoods2dPlot: values must be a matrix.")
   if( q == 1 & l == 1) stop("hoods2dPlot: values must be a matrix with at least one dimension > 1.")
   if( is.null( args$qs) & all( args$thresholds[,1] == args$thresholds[,2])) {
        a1.labels <- as.character( args$thresholds[,1])
        xlb <- paste("Threshold (", args$units, ")", sep="")
   } else if( !is.null( args$qs)) {
        a1.labels <- args$qs
        xlb <- "Threshold (quantile)"
   } else {
        a1.labels <- as.character( 1:dim( args$thresholds)[1])
        xlb <- "Threshold index number"
   }

   image( t( x), xaxt="n", yaxt="n", xlab=xlb, ylab="Neighborhood size (grid squares)", col=c("grey", heat.colors(12)), ...)
   text( x=seq(0,1,,q)[ rep(1:q,l)], y=seq(0,1,,l)[ rep(1:l,each=q)], labels=round( t( x), digits=2))
   axis( 1, at=seq(0,1,,q), labels=a1.labels)
   if( any( args$thresholds[,1] != args$thresholds[,2])) warning("hoods2dPlot: thresholds differ for the two fields, labels only for the first")
   axis( 2, at=seq(0,1,,length(args$levels)), labels=args$levels)
   image.plot( t( x), legend.only=TRUE, col=c("grey", heat.colors(12)), ...)

   matplot( x, ylab=args$ylab, xlab="Neighborhood size (grid squares)", type = "l", lty = 1, axes = FALSE , lwd = 2)
   axis(2)
   box()
   axis(1, at = 1:l, labels = args$levels)
   grid()
   legend("topleft", legend = a1.labels, col = 1:q, title = xlb, inset = 0.02, lwd = 2 )
   invisible()
} # end of 'hoods2dPlot' function.

upscale2dPlot <-
function(object, args, ...) {
   par( mfrow=c(4,2), bg="beige")
   try( hoods2dPlot( object$ets, args=args, main="Upscaling: Gilbert Skill Score (GSS)", ...))
   try( hoods2dPlot( object$ts, args=args, main="Upscaling: Threat Score (TS)", ...))
   try( hoods2dPlot( object$bias, args=args, main="Upscaling: Bias", ...))
   try( plot( args$levels, object$rmse, type="b", xlab="Neighborhood Length (grid squares)", ylab="RMSE", main="Upscaling: RMSE", col="darkblue"))
} # end of 'upscale2dPlot' function.

plot.upscale2d <- function(x, ...) {
   upscale2dPlot( object=x, args=list( levels=x$levels, thresholds=x$thresholds, units=x$units, qs=x$qs), ...)
   invisible()
} # end of 'plot.upscale2d' function.

fss2dfun <- function(sPy, sPx, subset=NULL, verbose=FALSE) {
   ## 
   ## Function to calculate the FSS neighborhood statistics.  Not to be confused
   ## with 'fss2d', which does the same thing, but this function begins with the smoothed
   ## fields, and only calculates the final score.  To be used internally by 'fss2d'.
   ## 
   ## Arguments:
   ## 
   ## 'sPy', 'sPx' smoothed 'k X m' forecast and observed matrices, resp., output from 'kernel2dsmooth'.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: single numeric giving the FSS value.
   ##
   id1 <- !is.na( sPy)
   id2 <- !is.na( sPx)
   if( verbose) cat("Finding the numbers of non-missing values for each field.\n")
   if( is.null( subset)) {
      N1 <- sum( colSums( id1, na.rm=TRUE), na.rm=TRUE)
      N2 <- sum( colSums( id2, na.rm=TRUE), na.rm=TRUE)
      Nxy <- sum( colSums( id1 & id2, na.rm=TRUE), na.rm=TRUE)
      if( verbose) cat(Nxy, " total number of non-missing points.\n")
      num <- sum( colSums( (sPy - sPx)^2, na.rm=TRUE), na.rm=TRUE)/Nxy
      if( verbose) cat("MSE = ", num, "\n")
      denom <- sum( colSums( sPy^2, na.rm=TRUE), na.rm=TRUE)/N1 + sum( colSums( sPx^2, na.rm=TRUE), na.rm=TRUE)/N2
   } else {
      num <- mean( (c(sPy)[subset]-c(sPx)[subset])^2, na.rm=TRUE)
      denom <- mean( (c(sPy)[subset])^2+(c(sPx)[subset])^2, na.rm=TRUE)
   }
   if( verbose) cat("Reference MSE = ", denom, "\n")
   return(1-num/denom)
} # end of 'fss2dfun' function.

pphindcast2d <- function(obj, which.score="ets", verbose=FALSE, ...) {
   Fcst <- get( obj$Fcst.name)
   Obs <- get( obj$Vx.name)
   xdim <- obj$xdim
   Nxy <- obj$Nxy
   subset <- obj$subset
   thresholds <- obj$thresholds
   levels <- obj$levels
   q <- dim( thresholds)[1]
   l <- length( levels)
   outmat <- Pthresh <- matrix( NA, l, q)
   out <- list()
   out$which.score <- which.score
   findthresh <- function( p, Ix, sPx, binmat, which.score, subset=NULL) {
	sIx <- binmat
	sIx[ sPx >= p] <- 1
	return( -vxstats( sIx, Ix, which.stats=which.score, subset=subset)[[which.score]])
   } # end of 'findthresh' internal function.
   binmat <- matrix(0, xdim[1], xdim[2])
   for( threshold in 1:q) {
	Ix <- Iy <- binmat
	Ix[ Obs >= thresholds[threshold,2]] <- 1
	Iy[ Fcst >= thresholds[threshold,1]] <- 1
	for( level in 1:l) {
	   Wlvl <- kernel2dsmooth( Ix, kernel.type="boxcar", n=levels[ level], setup=TRUE)
	   sPy <- kernel2dsmooth( Iy, kernel.type="boxcar", n=levels[level], W=Wlvl, xdim=xdim, Nxy=Nxy)
	   sPx <- kernel2dsmooth( Ix, kernel.type="boxcar", n=levels[level], W=Wlvl, xdim=xdim, Nxy=Nxy)
   	   Pthresh[ level, threshold] <- optim( 0, findthresh, Ix=Ix, sPx=sPx, binmat=binmat, which.score=which.score, subset=subset, 
							lower=0, upper=1, method="L-BFGS-B", ...)$par
	   if( verbose) cat("Pthresh = ", Pthresh[level,threshold], " for obs threshold = ", thresholds[ threshold,2], " and level = ", levels[level], "\n")
	   sIy <- binmat
	   sIy[ sPy >= Pthresh[ level, threshold]] <- 1
	   outmat[ level, threshold] <- vxstats( sIy, Ix, which.stats=which.score)[[which.score]]
	} # end of for 'level' loop.
   } # end of for 'threshold' loop.
   out$values <- outmat
   out$Pthresh <- Pthresh
   return( out)
} # end of 'pphindcast' function.

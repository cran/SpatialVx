locmeasures2dPrep <- function(object, k=NULL, alpha=0.1, bdconst=NULL, p=2) {

   out <- object
   attr(out, "alpha") <- alpha
   if(is.null(bdconst)) attr(out, "bdconst") <- Inf
   else attr(out, "bdconst") <- bdconst
   attr(out, "p") <- p
   attr(out, "k") <- k

   return(out)

} # end of 'locmetric2dPrep' function.

locmeasures2d <- function(object, which.stats=c("bdelta", "haus", "qdmapdiff", "med", "msd", "ph", "fom"),
        distfun="distmapfun", distfun.params=NULL, k=NULL, alpha=0.1, bdconst=NULL, p=2, ...) {

    UseMethod("locmeasures2d", object)

} # end of 'locmeasures2d' function.

locmeasures2d.default <- function(object, which.stats=c("bdelta", "haus", "qdmapdiff", "med", "msd", "ph", "fom"),
        distfun="distmapfun", distfun.params=NULL, k=NULL, alpha=0.1, bdconst=NULL, p=2, ..., Y, thresholds=NULL) {

    args <- list(...)
    if(missing(Y) && is.element("Xhat", names(args))) Y <- args$Xhat
    if( !all( dim( Y ) == dim( object ) ) ) {

	stopmsg <- paste("locmeasures2d: dim of observed field (", dim( object ), 
		") must equal dim of forecast field (", dim( Y ), ")", sep = "" )
	stop( stopmsg )

    }

    obj <- make.SpatialVx(X=object, Xhat=Y, thresholds=thresholds)

    out <- locmeasures2d.SpatialVx(object=obj, which.stats=which.stats, distfun=distfun,
						k=k, alpha=alpha, bdconst=bdconst, p=p, ...)
    return(out)

} # end of 'locmeasures2d.default' function.

locmeasures2d.SpatialVx <- function(object, which.stats=c("bdelta", "haus", "qdmapdiff", "med", "msd", "ph", "fom"),
	distfun="distmapfun", distfun.params=NULL, k=NULL, alpha=0.1, bdconst=NULL, p=2, ..., time.point = 1, obs = 1, model = 1 ) {

    object <- locmeasures2dPrep(object=object, k=k, alpha=alpha, bdconst=bdconst, p=p)
    a <- attributes(object)

    thresholds <- a$thresholds
    q <- dim( thresholds$X )[ 1 ]

    if(is.null(a$k) && "qdmapdiff" %in% which.stats) stop("locmeasures2d: must supply k in call to locmeasures2dPrep to do qdmapdiff method.")

    if(!is.null(a$k)) nk <- length(a$k)
    else nk <- 1

    if(!is.null(a$p)) np <- length(a$p)
    else np <- 1

    if(is.null(a$alpha)) nalpha <- 1
    else nalpha <- length(a$alpha)

    out <- LocListSetup(a=a, which.stats=which.stats, nthresh=q, np=np, nk=nk, nalpha=nalpha)
    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    ## Begin: Get the data sets
    dat <- datagrabber( object, time.point = time.point, obs = obs, model = model )

    Obs <- dat$X
    Fcst <- dat$Xhat

    mainname <- a$data.name
    vxname <- a$obs.name[ obs ]
    attr(out, "data.name") <- c(mainname, vxname, a$model.name[ model ] )
    attr(out, "thresholds") <- list( X = thresholds$X[, obs ], Xhat = thresholds$Xhat[ , model ] )
    ## End: Get the data sets

    xdim <- a$xdim
    x.id <- im(Obs)
    y.id <- im(Fcst)

    for(threshold in 1:q) {

	Ix <- solutionset( x.id >= thresholds$X[threshold, obs ] )
	Iy <- solutionset( y.id >= thresholds$Xhat[threshold, model ] )

	if( "bdelta" %in% which.stats) for(p in 1:np) { 
			tmpDelta <- try(deltametric(Iy,Ix,p=a$p[p], c=a$bdconst, ...))
			if(class(tmpDelta) != "try-error") out$bdelta[p, threshold] <- tmpDelta
			} # end of for 'p' loop.

	if( "haus" %in% which.stats) out$haus[threshold] <- deltametric(Iy,Ix,p=Inf,c=Inf, ...)

	if( "qdmapdiff" %in% which.stats) for( k in 1:nk) out$qdmapdiff[k,threshold] <- locperf(X=Ix, Y=Iy, which.stats="qdmapdiff", k=a$k[k],
											distfun=distfun, distfun.params)$qdmapdiff
	if( w.id <- any( c("med", "msd") %in% which.stats)) {

	   tmp <- locperf( X = Ix, Y = Iy, which.stats = c( "med", "msd" )[ w.id ], distfun = distfun, distfun.params )
	   if( "med" %in% which.stats) {

		out$medMiss[ threshold ] <- tmp$medMiss
		out$medFalseAlarm[ threshold ] <- tmp$medFalseAlarm

	   }

	   if( "msd" %in% which.stats) {

		out$msdMiss[ threshold ] <- tmp$msdMiss
		out$msdFalseAlarm[ threshold ] <- tmp$msdFalseAlarm

	    }

	} # end of if do partial and/or modified Hausdorff metric.

	if( "fom" %in% which.stats) {

	   for( i in 1:nalpha) {

		out$fom[i,threshold] <- locperf(X=Ix, Y=Iy, which.stats="fom", alpha=a$alpha[i], distfun=distfun, distfun.params)$fom

	   } # end of for 'a' loop.

	}

   } # end of for 'threshold' loop.

   class(out) <- "locmeasures2d"
   return(out)

} # end of 'locmeasures.SpatialVx' function.

metrV <- function(x, ...) {
    UseMethod("metrV", x)
} # end of 'metrV' function.

metrV.SpatialVx <- function(x, time.point=1, obs = 1, model=1, lam1=0.5, lam2=0.5, distfun="distmapfun", verbose=FALSE, ...) {

    a <- attributes(x)
    out <- list()
    attributes(out) <- a

    u <- a$thresholds
    thresholds <- cbind( u$X[, obs ], u$Xhat[, model ] )

    mainname <- a$data.name
    vxname <- a$obs.name[ obs ]

    if(length(model) == 1) {

	dat <- datagrabber(x, time.point=time.point, obs = obs, model = model )
	X <- dat$X
        Xhat <- dat$Xhat

	out <- metrV.default(x=X, xhat=Xhat, thresholds=thresholds, lam1=lam1, lam2=lam2, distfun=distfun,
				a=a, verbose=verbose, ...)

    } else {

	if(length(model) != 2) stop("metrV.SpatialVx: invalid model argument.  Must have length 1 or 2.")

	dat <- datagrabber(x, time.point=time.point, obs = obs, model = model[ 1 ] )
	dat2 <- datagrabber(x, time.point=time.point, obs = obs, model=model[ 2 ] )

	X <- dat$X
        Xhat <- dat$Xhat
        Xhat2 <- dat2$Xhat

	out <- metrV.default(x=X, xhat=Xhat, xhat2=Xhat2, thresholds=thresholds, lam1=lam1, lam2=lam2, distfun=distfun,
                                a=a, verbose=verbose, ...)

    } # end of if else 1 or 2 models stmts.

    attr(out, "data.name") <- c(mainname, vxname, a$model.name[ model ] )
    attr(out, "thresholds") <- thresholds
    a <- attributes(out)

    return(out)
} # end of 'metrV.SpatialVx' function.

metrV.default <- function(x, xhat, xhat2=NULL, thresholds, lam1=0.5, lam2=0.5, distfun="distmapfun", a=NULL, verbose=FALSE, ...) {

    if(!is.matrix(x) || !is.matrix(xhat)) stop("metrV.default: invalid x and/or xhat argument.")
    if((!missing(xhat2) || !is.null(xhat2)) && !is.matrix(xhat2)) stop("metrV.default: xhat2 must be NULL or a matrix.")

    M1im <- im(xhat)
    Oim <- im(x)
    if(is.m2 <- !is.null(xhat2)) M2im <- im(xhat2)
    
    q <- dim(thresholds)[1]
    out <- list()

    if(!is.null(a)) attributes(out) <- a
    else attr(out, "thresholds") <- thresholds

    out$OvsM1 <- matrix(NA, q, 3)
    colnames( out$OvsM1) <- c("distOV", "distob", "metrV")

    if(is.m2) {
	out$OvsM2 <- out$OvsM1
	out$M1vsM2 <- out$OvsM1
    }

    for(threshold in 1:q) {
	Ix <- solutionset(Oim >= thresholds[threshold,1])
	Im1 <- solutionset(M1im >= thresholds[threshold,2])
	if(is.m2) Im2 <- solutionset(M2im >= thresholds[threshold,3]) 

	out$OvsM1[threshold,1] <- sqrt(sum(colSums((Ix$m - Im1$m)^2, na.rm=TRUE), na.rm=TRUE))
	out$OvsM1[threshold,2] <- distob( Ix, Im1, distfun=distfun, ...)
	out$OvsM1[threshold,3] <- lam1*out$OvsM1[threshold,1] + lam2*out$OvsM1[threshold,2]

	if(verbose) cat("O vs M1 distOV for thresholds: ", thresholds[threshold,], " = ", out$OvsM1[threshold,], "\n")

	if(is.m2) {
	   out$OvsM2[threshold,1] <- sqrt(sum(colSums((Ix$m - Im2$m)^2, na.rm=TRUE), na.rm=TRUE))
           out$OvsM2[threshold,2] <- distob(Ix, Im2, distfun=distfun, ...)
           out$OvsM2[threshold,3] <- lam1*out$OvsM2[threshold,1] + lam2*out$OvsM2[threshold,2]

	   if(verbose) cat("O vs M2 distOV for thresholds: ", thresholds[threshold,], " = ", out$OvsM2[threshold,], "\n")

	   out$M1vsM2[threshold,1] <- sqrt(sum(colSums((Im1$m - Im2$m)^2, na.rm=TRUE), na.rm=TRUE))
	   out$M1vsM2[threshold,2] <- abs( out$OvsM1[threshold,2] - out$OvsM2[threshold,2]) 
	   out$M1vsM2[threshold,3] <- lam1*out$M1vsM2[threshold,1] + lam2*out$M1vsM2[threshold,2] 

	   if(verbose) cat("M1 vs M2 distOV for thresholds: ", thresholds[threshold,], " = ", out$M1vsM2[threshold,], "\n")
	}
    } # end of for 'threshold' loop.
    class( out) <- "metrV"
    return( out)
} # end of 'metrV.default' function.

print.metrV <- function(x, ...) {

    a <- attributes(x)

    if(!is.null(a$msg)) print(a$msg)
    if(!is.null(a$data.name)) {
	cat("\n", "Comparison for:\n")
	print(a$data.name)
    }

    cat("\n", "Thresholds applied:\n")
    if(!is.null(a$qs)) print(a$qs)
    else print(a$thresholds)

    cat("\n", "Observed field vs Forecast 1:\n")
    print(x$OvsM1)

    if(!is.null(x$OvsM2)) {
	cat("\n", "Observed field vs Forecast 2:\n")
	print(x$OvsM2)
	cat("\n", "Forecast 1 vs Forecast 2:\n")
	print(x$M1vsM2)
    }

    invisible()
} # end of 'print.metrV' function.

distob <- function(X,Y, distfun="distmapfun", ...) {
   # X and Y are objects output from 'solutionset'.
   #
   # If both X and Y contain events, then this function returns the mean error distance with
   # respect to X (i.e., averaging over pixels in X the distances d(x,Y)).
   #
   # Otherwise, it returns zero if neither X nor Y contain any events, and max( dim(X))
   # if only one contains no events.
   nX <- sum(colSums(as.matrix( X ), na.rm=TRUE),na.rm=TRUE)
   nY <- sum(colSums(as.matrix( Y ), na.rm=TRUE), na.rm=TRUE)
   if( nX==0 & nY==0) return(0)
   else if( nX==0 | nY==0) return( max( dim( as.matrix( X ) ), na.rm=TRUE) )
   else out <- locperf( X=X, Y=Y, which.stats = "med", distfun=distfun, ...)$medMiss
   return(out)
}


distmapfun <- function(x, ...) {

   return( distmap(x, ...)$v)

} # end of 'distmapfun' function.

locperf <- function(X, Y, which.stats = c("qdmapdiff", "med", "msd", "ph", "fom", "minsep"), alpha=0.1, k=4, distfun="distmapfun", a=NULL, ...) { 

   out = LocListSetup(a=a, which.stats, nthresh=1)

   bb = boundingbox(as.rectangle(X), as.rectangle(Y))
   X = rebound( X, bb )
   Y = rebound( Y, bb )
   # dY = distmap(Y, ...)
   dY = do.call(distfun, list(x=Y, ...))

   # dX <- distmap(X, ...)
   dX = do.call(distfun, list(x = X, ...))

    if(!any( as.matrix( Y ))) {

	dY = unique(c(dY))
	dYcheck = FALSE

    } else dYcheck = TRUE 

    if(!any( as.matrix( X ))) dXcheck = FALSE
    else dXcheck = TRUE

   if(any(c("med", "msd", "fom", "minsep", "ph" ) %in% which.stats)) {

	if(dYcheck) Z = dY[ as.matrix( X ) ]
	else if(any(as.matrix( X ))) {

	    Z = as.matrix( X ) + NA
	    Z[ as.matrix( X ) ] <- dY

	} else Z = 0

	if( dXcheck ) Zother = dX[ as.matrix( Y ) ]
	else if( any( as.matrix( Y ) ) ) {

	    Zother = as.matrix( Y ) + NA
	    Zother[ as.matrix( Y ) ] <- dX

	} else Zother = 0

    }

   if(any(c("msd", "fom") %in% which.stats)) {

	Z2      = Z^2
	Zother2 = Zother^2

	if( "fom" %in% which.stats) {

	    if(dYcheck && dXcheck) N <- max( sum( colSums( as.matrix( X ), na.rm=TRUE), na.rm=TRUE),
						sum( colSums( as.matrix( Y ), na.rm=TRUE), na.rm=TRUE), na.rm=TRUE)
	    else if(dYcheck && !dXcheck) N <- sum( colSums( as.matrix( Y ), na.rm=TRUE), na.rm=TRUE)
	    else if(!dYcheck && dXcheck) N <- sum( colSums( as.matrix( X ), na.rm=TRUE), na.rm=TRUE)
	    else N <- 1e16

	} # end of if do fom stmts.

   }

  if( "ph" %in% which.stats ) {

    if( k < 1 || 
	(k <= sum( colSums( as.matrix( X ), na.rm = TRUE ), na.rm = TRUE ) || 
	 k <= sum( colSums( as.matrix( Y ), na.rm = TRUE ), na.rm = TRUE ) ) ) {

	    if( dXcheck || dYcheck ) {

	        if( k >= 1 ) out$ph <- max( c( sort( Z, decreasing = TRUE )[ k ], sort( Zother, decreasing = TRUE )[ k ] ), na.rm = TRUE )
	        else if( k >= 0 && k < 1 ) out$ph <- max( quantile( Z, probs = k ), quantile( Zother, probs = k ), na.rm = TRUE )
	        else out$ph <- NA

	    } else out$ph = max( sort( Z, decreasing = TRUE )[ 1 ], sort( Zother, decreasing = TRUE )[ 1 ], na.rm = TRUE )

       } else {

	out$ph <- NA

	} # end of if else make sure there are enough events to calculate the k-th highest value.
   }

   if( "qdmapdiff" %in% which.stats ) {

	diffXY = sort( c(abs(dX - dY)), decreasing=TRUE)

	if( "qdmapdiff" %in% which.stats) {

	    if(dXcheck || dYcheck) {

	   	if( k >= 1) out$qdmapdiff = diffXY[k]
	   	else if( k >= 0 && k < 1) out$qdmapdiff = quantile( diffXY, probs=k)
	   	else out$qdmapdiff = NA

	    } else out$qdmapdiff <- diffXY[1]

	}

   } # end of if do 'qdmapdiff' stmts.

   if( "med" %in% which.stats) {

	out$medMiss <- mean( Z, na.rm=TRUE)
	out$medFalseAlarm <- mean( Zother, na.rm = TRUE )

   }

   if( "msd" %in% which.stats) {

	out$msdMiss <- mean( Z2, na.rm=TRUE )
	out$msdFalseAlarm <- mean( Zother2, na.rm = TRUE )

   }

   if( "fom" %in% which.stats) out$fom <- sum( 1 / ( 1 + alpha * Z2 ), na.rm=TRUE ) / N
   if( "minsep" %in% which.stats) out$minsep <- min( c(Z), na.rm=TRUE )

   return( out)

} # end of 'locperf' function.

LocListSetup <- function(a, which.stats= c("bdelta", "haus", "qdmapdiff", "med", "msd", "ph", "fom", "minsep"),
			    nthresh=1, np=1, nk=1, nalpha=1) {

   out <- list()
   attributes(out) <- a

   q <- nthresh
   outvec <- numeric(q)+NA
   if( "bdelta" %in% which.stats ) out$bdelta <- matrix( NA, np, q)
   if( "haus" %in% which.stats ) out$haus <- outvec
   if( "qdmapdiff" %in% which.stats) out$qdmapdiff <- matrix( NA, nk, q)
   if( "med" %in% which.stats ) out$medMiss <- out$medFalseAlarm <- outvec
   if( "msd" %in% which.stats ) out$msdMiss <- out$msdFalseAlarm <- outvec
   if( "ph" %in% which.stats ) out$ph <- matrix( NA, nk, q )
   if( "fom" %in% which.stats ) out$fom <- matrix( NA, nalpha, q)
   if( "minsep" %in% which.stats) out$minsep <- outvec
   return( out)
} # end of 'LocListSetup' function.

summary.locmeasures2d <- function(object, ...) {
   x <- attributes(object)
   u <- x$thresholds
   if(!is.null( x$qs)) lu <- x$qs
   else if( all( u[,1] == u[,2])) lu <- as.character( u[,1])
   else lu <- paste( "mod thresh ", u[,1], ", vx thresh ", u[,2], sep="")
   k <- x$k
   p <- x$p
   a <- x$alpha
   print(x$msg)
   cat("\n", "Comparison for:\n")
   print(x$data.name)
   cat("Threshold(s) is (are):\n")
   print(lu)
   if( !is.null( object$bdelta)) {
	y <- object$bdelta
	rownames(y) <- paste( "p = ", p, "; ", sep="")
	colnames(y) <-  lu
	cat("Baddeley Delta Metric with c = ", x$bdconst, "\n")
	print(y)
   }
   if( !is.null( object$haus)) {
	y <- object$haus
	y <- matrix( y, nrow=1)
	colnames(y) <- lu
	cat("\n", "Hausdorff distance\n")
	print(y)
   }
   if( !is.null( object$qdmapdiff)) {
	y <- object$qdmapdiff
	rownames( y) <- paste("k = ", as.character( k), "; ", sep="")
	colnames( y) <- lu
	cat("\n", "Quantile (if k in (0,1) or k-th highest (if k = 1, 2, ...) difference in distance maps.\n")
	print( y)
   }
   if( !is.null( object$medMiss)) {
        y <- object$medMiss
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean error distance (Miss)\n")
        print( y)
   }
   if( !is.null( object$medFalseAlarm)) {
        y <- object$medFalseAlarm
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean error distance (FalseAlarm)\n")
        print( y)
   }
   if( !is.null( object$msdMiss)) {
	y <- object$msdMiss
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean square error distance (Miss)\n")
        print( y)
   }
   if( !is.null( object$msdFalseAlarm )) {
        y <- object$msdFalseAlarm 
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean square error distance (FalseAlarm)\n")
        print( y)
   }

   if( !is.null( object$qdmapdiff)) {
        y <- object$ph
        rownames( y) <- paste("k = ", as.character( k), "; ", sep="")
        colnames( y) <- lu
        cat("\n", "Partial Hausdorff distance\n")
        print( y)
   }

   if( !is.null( object$fom)) {
        y <- object$fom
	rownames( y) <- paste("alpha = ", x$alpha, "; ", sep="")
	colnames( y) <- lu
        cat("\n", "Pratt\'s figure of merit (FOM)\n")
        print( y)
   }
   invisible()
} # end of 'summary.locmeasures2d' function.

print.locmeasures2d <- function(x, ...) {
    a <- attributes(x)
    print(a$msg)
    cat("Comparison for: \n")
    print(a$data.name)

    print(summary(x, ...))

    invisible()
} # end of 'print.locmeasures2d' function.

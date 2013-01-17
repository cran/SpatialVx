locmeasures2dPrep <- function( Vx.name, Fcst.name, thresholds=NULL, k=NULL, alpha=0.1, bdconst=NULL, p=2, loc=NULL, qs=NULL, units=NULL) {
   out <- list()
   data.name <- c(Vx.name, Fcst.name)
   names(data.name) <- c("verification","forecast")
   out$data.name <- data.name
   Fcst <- get(Fcst.name)
   Obs <- get(Vx.name)
   if (is.null(thresholds)) {
        thresholds <- cbind(quantile(c(Fcst), probs = c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95)),
			    quantile(c(Obs), probs = c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95)))
        qs <- as.character(c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95))
   }
   else if (!is.matrix(thresholds)) thresholds <- cbind(thresholds, thresholds)
   out$thresholds <- thresholds
   out$qs <- qs
   out$alpha <- alpha
   out$loc <- loc
   out$units <- units
   out$xdim <- dim( Fcst)
   out$Nxy <- prod( out$xdim)
   if( is.null( bdconst)) out$bdconst <- Inf
   else out$bdconst <- bdconst
   out$p <- p
   out$k <- k
   return( out)
} # end of 'locmetric2dPrep' function.

locmeasures2d <- function( object, which.stats=c("baddeley", "hausdorff", "ph", "mhd", "med", "msd", "fom"),
	distfun="distmapfun", distfun.params=NULL, ...) {
   thresholds <- object$thresholds
   q <- dim( thresholds)[1]
   if( is.null(object$k) & "ph" %in% which.stats) stop("locmeasures2d: must supply k in call to locmeasures2dPrep to do ph method.")
   if( !is.null(object$k)) nk <- length( object$k)
   else nk <- 1
   if( !is.null(object$p)) np <- length( object$p)
   else np <- 1
   if( is.null(object$alpha)) nalpha <- 1
   else nalpha <- length( object$alpha)
   out <- LocListSetup(which.stats=which.stats, nthresh=q, np=np, nk=nk, nalpha=nalpha)
   out$prep.obj <- as.character(substitute(object))
   Obs <- get( object$data.name[1])
   Fcst <- get( object$data.name[2])
   xdim <- object$xdim
   x.id <- im( Obs)
   y.id <- im( Fcst)
   for( threshold in 1:q) {
	Ix <- solutionset(x.id>=thresholds[threshold,2])
	Iy <- solutionset(y.id>=thresholds[threshold,1])
	if( "baddeley" %in% which.stats) for( p in 1:np) { 
			tmpDelta <- try(deltametric(Iy,Ix,p=object$p[p], c=object$bdconst, ...))
			if(class(tmpDelta) != "try-error") out$baddeley[p, threshold] <- tmpDelta
			} # end of for 'p' loop.
	if( "hausdorff" %in% which.stats) out$hausdorff[threshold] <- deltametric(Iy,Ix,p=Inf,c=Inf, ...)
	if( "ph" %in% which.stats) for( k in 1:nk) out$ph[k,threshold] <- locperf(X=Ix, Y=Iy, which.stats="ph", k=object$k[k],
												distfun=distfun, distfun.params)$ph
	if( w.id <- any( c("mhd", "med", "msd") %in% which.stats)) {
	   tmp <- locperf(X=Ix, Y=Iy, which.stats=c("ph","mhd", "med", "msd")[w.id], distfun=distfun, distfun.params)
	   if( "mhd" %in% which.stats) out$mhd[threshold] <- tmp$mhd
	   if( "med" %in% which.stats) out$med[threshold] <- tmp$med
	   if( "msd" %in% which.stats) out$msd[threshold] <- tmp$msd
	} # end of if do partial and/or modified Hausdorff metric.
	if( "fom" %in% which.stats) {
	   for( i in 1:nalpha) {
		out$fom[i,threshold] <- locperf(X=Ix, Y=Iy, which.stats="fom", alpha=object$alpha[i], distfun=distfun, distfun.params)$fom
	   } # end of for 'a' loop.
	}
   } # end of for 'threshold' loop.
   class( out) <- "locmeasures2d"
   return( out)
} # end of 'locmeasures function.

metrV <- function(object1, object2=NULL, lam1=0.5, lam2=0.5, distfun="distmapfun", verbose=FALSE, ...) {
   ## This works, and is fast, but does not match results in Zhu et al.
   M1 <- get( object1$data.name[2])
   M1im <- im( M1)
   O <- get( object1$data.name[1])
   Oim <- im( O)
   if( is.m2 <- !is.null( object2)) {
	M2 <- get( object2$data.name[2])
	M2im <- im( M2)
   }
   thresholds <- object1$thresholds
   q <- dim( thresholds)[1]
   out <- list()
   out$prep.obj1 <- as.character( substitute( object1))
   out$OvsM1 <- matrix( NA, q, 3)
   colnames( out$OvsM1) <- c("distOV", "distob", "metrV")
   if( is.m2) {
	out$prep.obj2 <- as.character( substitute( object2))
	out$OvsM2 <- out$OvsM1
	out$M1vsM2 <- out$OvsM1
   }
   for( threshold in 1:q) {
	Ix <- solutionset( Oim >= thresholds[threshold,2])
	Im1 <- solutionset( M1im >= thresholds[threshold,1])
	if(is.m2) Im2 <- solutionset( M2im >= thresholds[threshold,1]) 
	out$OvsM1[threshold,1] <- sqrt( sum( colSums( (Ix$m - Im1$m)^2, na.rm=TRUE), na.rm=TRUE))
	out$OvsM1[threshold,2] <- distob( Ix, Im1, distfun=distfun, ...)
	out$OvsM1[threshold,3] <- lam1*out$OvsM1[threshold,1] + lam2*out$OvsM1[threshold,2]
	if( verbose) cat("O vs M1 distOV for thresholds: ", thresholds[threshold,], " = ", out$OvsM1[threshold,], "\n")
	if( is.m2) {
	   out$OvsM2[threshold,1] <- sqrt( sum( colSums( (Ix$m - Im2$m)^2, na.rm=TRUE), na.rm=TRUE))
           out$OvsM2[threshold,2] <- distob( Ix, Im2, distfun=distfun, ...)
           out$OvsM2[threshold,3] <- lam1*out$OvsM2[threshold,1] + lam2*out$OvsM2[threshold,2]
	   if( verbose) cat("O vs M2 distOV for thresholds: ", thresholds[threshold,], " = ", out$OvsM2[threshold,], "\n")
	   out$M1vsM2[threshold,1] <- sqrt( sum( colSums( (Im1$m - Im2$m)^2, na.rm=TRUE), na.rm=TRUE))
	   out$M1vsM2[threshold,2] <- abs( out$OvsM1[threshold,2] - out$OvsM2[threshold,2]) 
	   out$M1vsM2[threshold,3] <- lam1*out$M1vsM2[threshold,1] + lam2*out$M1vsM2[threshold,2] 
	   if( verbose) cat("M1 vs M2 distOV for thresholds: ", thresholds[threshold,], " = ", out$M1vsM2[threshold,], "\n")
	}
   } # end of for 'threshold' loop.
   class( out) <- "metrV"
   return( out)
} # end of 'metrV' function.

distob <- function(X,Y, distfun="distmapfun", ...) {
   # X and Y are objects output from 'solutionset'.
   #
   # If both X and Y contain events, then this function returns the mean error distance with
   # respect to X (i.e., averaging over pixels in X the distances d(x,Y)).
   #
   # Otherwise, it returns zero if neither X nor Y contain any events, and max( dim(X))
   # if only one contains no events.
   nX <- sum(colSums(X$m, na.rm=TRUE),na.rm=TRUE)
   nY <- sum(colSums(Y$m, na.rm=TRUE), na.rm=TRUE)
   if( nX==0 & nY==0) return(0)
   else if( nX==0 | nY==0) return(max(dim(X$m), na.rm=TRUE))
   else out <- locperf( X=X, Y=Y, which.stats="med", distfun=distfun, ...)$med
   return(out)
}

distmapfun <- function(x, ...) {
   return( distmap(x, ...)$v)
} # end of 'distmapfun' function.

locperf <- function(X,Y, which.stats=c("ph", "mhd", "med", "msd", "fom", "minsep"), alpha=0.1, k=4, distfun="distmapfun", ...) {
   out <- LocListSetup(which.stats, nthresh=1)
   bb <- bounding.box(as.rectangle(X), as.rectangle(Y))
   X <- rebound(X, bb)
   Y <- rebound(Y, bb)
   # dY <- distmap(Y, ...)$v
   dY <- do.call(distfun, list(x=Y, ...))
   if( any( c("med", "msd", "fom", "minsep") %in% which.stats)) Z <- dY[X$m]
   if( any( c("msd", "fom") %in% which.stats)) {
	Z2 <- Z^2
	if( "fom" %in% which.stats) N <- max( sum( colSums(X$m, na.rm=TRUE), na.rm=TRUE), sum( colSums(Y$m, na.rm=TRUE), na.rm=TRUE), na.rm=TRUE)
   }
   if( any( c("ph", "mhd") %in% which.stats)) {
	# dX <- distmap(X, ...)$v
	dX <- do.call(distfun, list(x=X, ...))
	diffXY <- sort( c( abs(dX - dY)), decreasing=TRUE)
	if( "ph" %in% which.stats) {
	   if( k >= 1) out$ph <- diffXY[k]
	   else if( k >= 0 & k < 1) out$ph <- quantile( diffXY, probs=k)
	   else out$ph <- NA
	}
	if( "mhd" %in% which.stats) {
	   V <- dX[Y$m]
	   out$mhd <- max( c( mean(Z, na.rm=TRUE), mean(V, na.rm=TRUE), na.rm=TRUE))
	}
   } # end of if do partial- and/or modified- Hausdorff stmts.
   if( "med" %in% which.stats) out$med <- mean( Z, na.rm=TRUE)
   if( "msd" %in% which.stats) out$msd <- mean( Z2, na.rm=TRUE)
   if( "fom" %in% which.stats) out$fom <- sum( 1/(1+alpha*Z2), na.rm=TRUE)/N
   if( "minsep" %in% which.stats) out$minsep <- min( c(Z), na.rm=TRUE)
   return( out)
}

LocListSetup <- function(which.stats= c("baddeley", "hausdorff", "ph", "mhd", "med", "msd", "fom", "minsep"), nthresh=1, np=1, nk=1, nalpha=1) {
   out <- list()
   q <- nthresh
   outvec <- numeric(q)+NA
   if( "baddeley" %in% which.stats) out$baddeley <- matrix( NA, np, q)
   if( "hausdorff" %in% which.stats) out$hausdorff <- outvec
   if( "ph" %in% which.stats) out$ph <- matrix( NA, nk, q)
   if( "mhd" %in% which.stats) out$mhd <- outvec
   if( "med" %in% which.stats) out$med <- outvec
   if( "msd" %in% which.stats) out$msd <- outvec
   if( "fom" %in% which.stats) out$fom <- matrix( NA, nalpha, q)
   if( "minsep" %in% which.stats) out$minsep <- outvec
   return( out)
} # end of 'LocListSetup' function.

summary.locmeasures2d <- function(object, ...) {
   x <- get( object$prep.obj)
   u <- x$thresholds
   if( !is.null( x$qs)) lu <- x$qs
   else if( all( u[,1] == u[,2])) lu <- as.character( u[,1])
   else lu <- paste( "mod thresh ", u[,1], ", vx thresh ", u[,2], sep="")
   k <- x$k
   p <- x$p
   a <- x$alpha
   cat("Comparison between model ", x$data.name[2], " and verification ", x$data.name[1], ":\n")
   cat("Threshold(s) is (are):\n")
   print( lu)
   if( !is.null( object$baddeley)) {
	y <- object$baddeley
	rownames( y) <- paste( "p = ", as.character( p), "; ", sep="")
	colnames( y) <-  lu
	cat("Baddeley Delta Metric with c = ", x$bdconst, "\n")
	print( y)
   }
   if( !is.null( object$hausdorff)) {
	y <- object$hausdorff
	y <- matrix( y, nrow=1)
	colnames( y) <- lu
	cat("\n", "Hausdorff distance\n")
	print( y)
   }
   if( !is.null( object$ph)) {
	y <- object$ph
	rownames( y) <- paste("k = ", as.character( k), "; ", sep="")
	colnames( y) <- lu
	cat("\n", "Partial Hausdorff distance\n")
	print( y)
   }
   if( !is.null( object$mhd)) {
	y <- object$mhd
	y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "modified Hausdorff distance\n")
        print( y)
   }
   if( !is.null( object$med)) {
        y <- object$med
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean error distance\n")
        print( y)
   }
   if( !is.null( object$msd)) {
	y <- object$msd
        y <- matrix( y, nrow=1)
        colnames( y) <- lu
        cat("\n", "Mean square error distance\n")
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
}

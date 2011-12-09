FeatureSuite <- function(prep.object, fun.object, verbose=FALSE) {
   if( verbose) begin.time <- Sys.time()
   Y <- get( prep.object$Fcst.name)
   X  <- get( prep.object$Vx.name)
   xdim <- prep.object$xdim
   Yim <- im( Y)
   Xim <- im( X)
   out <- list()
   if( verbose) cat("\n", "Calling feature identification function ", fun.object$identify$fun, "\n")
   feats <- do.call(fun.object$identify$fun, c( list(object=prep.object), fun.object$identify$args))
   if( !is.null(fun.object$merge)) {
      if( verbose) cat("\n", "Merging features within fields with ", fun.object$merge$fun, "\n")
      feats2 <- do.call( fun.object$merge$fun, c( list( x=feats, object=prep.object), fun.object$merge$args))
   } else feats2 <- NULL
   if( !is.null(fun.object$match)) {
      if( verbose) cat("\n", "Matching features between fields with ", fun.object$match$fun, "\n")
      matches <- do.call( fun.object$match$fun, c( list( x=feats, y=feats2, object=prep.object), fun.object$match$args))
   } else matches <- NULL
   if( verbose) cat("\n", "Analyzing features with ", fun.object$analysis$fun, "\n")
   out <- do.call( fun.object$analysis$fun, c( list( x=feats, y=feats2, matches=matches, object=prep.object), fun.object$analysis$args))
   out$prep.obj <- as.character(substitute( prep.object))
   out$fun.obj <- as.character(substitute( fun.object))
   if( verbose) print( Sys.time() - begin.time)
   class( out) <- "FeatureSuite"
   return( out)
} # end of 'FeatureSuite' function.

FeatureFunPrep <- function(identfun=NULL, identfun.args=NULL,
			mergefun=NULL, mergefun.args=NULL,
			matchfun=NULL, matchfun.args=NULL,
			analysisfun=NULL, analysisfun.args=NULL) {
   out <- list()
   if( !is.null( identfun)) out$identify <- list( fun=identfun, args=identfun.args)
   if( !is.null( mergefun)) out$merge <- list( fun=mergefun, args=mergefun.args)
   if( !is.null( matchfun)) out$match <- list( fun=matchfun, args=matchfun.args)
   if( !is.null( analysisfun)) out$analysis <- list( fun=analysisfun, args=analysisfun.args)
   return( out)
} # end of 'FeatureFunPrep' function.

FeatureSuitePrep <- function(Fcst.name, Vx.name, loc=NULL, units=NULL) {
   out <- list()
   out$Fcst.name <- Fcst.name
   out$Vx.name <- Vx.name
   Y <- get( Fcst.name)
   X <- get( Vx.name)
   xdim <- dim( X)
   out$xdim <- xdim
   out$loc <- loc
   out$units <- units
   class( out) <- "FeaturesSuitePrep"
   return( out)
} # end of 'FeatureSuitePrep' function.

convthresh <- function(object, smoothfun="disk2dsmooth", smoothpar=1, smoothfunargs=NULL, thresh=1e-8, idfun="disjointer", zero.down=FALSE, ...) {
   X <- get(object$Vx.name)
   Y <- get(object$Fcst.name)
   xdim <- dim(X)
   Xsm <- do.call(smoothfun, c(list( x=X, lambda=smoothpar), smoothfunargs))
   Ysm <- do.call(smoothfun, c(list( x=Y, lambda=smoothpar), smoothfunargs))
   if( zero.down) {
	Xsm[ Xsm < 0] <- 0
	Xsm <- zapsmall( Xsm)
	Ysm[ Ysm < 0] <- 0
	Ysm <- zapsmall( Ysm)
   }
   if( is.null( thresh)) thresh <- c(quantile(c(Xsm), 0.75), quantile(c(Ysm), 0.75))
   else if( length( thresh) == 1) thresh <- c(thresh, thresh)
   sIx <- sIy <- matrix(0, xdim[1], xdim[2])
   sIx[ Xsm >= thresh[2]] <- 1
   sIy[ Ysm >= thresh[1]] <- 1
   X.feats <- do.call(idfun, c( list( x=sIx), list(...)))
   Y.feats <- do.call(idfun, c( list( x=sIy), list(...)))
   Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])
   if(!is.null(X.feats)) for( i in 1:length( X.feats)) Xlab[ X.feats[[i]][["m"]]] <- i
   else X.feats <- NULL
   if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[ Y.feats[[j]][["m"]]] <- j
   else Y.feats <- NULL
   if(is.null(X.feats)) warning("convthresh: No values above threshold in verification field.")
   if(is.null(Y.feats)) warning("convthresh: No values above threshold in forecast field.")
   return( list( X.feats=X.feats, Y.feats=Y.feats, X.labeled=Xlab, Y.labeled=Ylab))
} # end of 'convthresh' function.

disjointer <- function(x, method="C") {
   x[ x==0] <- NA
   if(any(!is.na(x))) {
	out <- as.im( x)
   	out <- connected( X=out, method=method)
   	out <- tiles(tess(image=out))
   } else out <- NULL
   return( out)
} # end of 'disjointer' function.

salIDfun <- function(object, fac=0.06666667, q=0.95, wash.out=NULL, thresh=NULL, idfun="disjointer", ...) {
   X <- get(object$Vx.name)
   Y <- get(object$Fcst.name)
   xdim <- object$xdim
   Ix <- Iy <- matrix(0, xdim[1], xdim[2])
   if(is.null(thresh)) {
	if(is.null(wash.out)){
	   thresh <- quantile(c(Y), probs=q)
	   thresh <- c(thresh, quantile(c(X), probs=q))
	} else {
	   thresh <- quantile(c(Y[Y>=wash.out]), probs=q)
	   thresh <- c(thresh, quantile(c(X[X>=wash.out]), probs=q))
	}
	thresh <- thresh*fac
   } else if(length(thresh)==1) thresh <- c(thresh, thresh)
   Ix[X >= thresh[2]] <- 1
   Iy[Y >= thresh[1]] <- 1
   X.feats <- do.call(idfun, c(list(x=Ix), list(...)))
   Y.feats <- do.call(idfun, c(list(x=Iy), list(...)))
   Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])
   if(!is.null(X.feats)) for( i in 1:length( X.feats)) Xlab[ X.feats[[i]][["m"]]] <- i
   else X.feats <- NULL
   if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[ Y.feats[[j]][["m"]]] <- j
   else Y.feats <- NULL
   if(is.null(X.feats)) warning("salIDfun: No values above threshold in verification field.")
   if(is.null(Y.feats)) warning("salIDfun: No values above threshold in forecast field.")
   return( list( X.feats=X.feats, Y.feats=Y.feats, X.labeled=Xlab, Y.labeled=Ylab))
} # end of 'salIDfun' function.

saller <- function(x, object, y=NULL, matches=NULL, d=NULL) {
   out <- list()
   if(is.null(y)) tmp <- x
   else tmp <- y
   y <- tmp$Y.feats
   x <- tmp$X.feats
   binX <- im(tmp$X.labeled)
   binX <- solutionset(binX > 0)
   binY <- im(tmp$Y.labeled)
   binY <- solutionset(binY > 0)

   # Amplitude
   X <- get(object$Vx.name, pos=".GlobalEnv")
   Y <- get(object$Fcst.name, pos=".GlobalEnv")
   DomRmod <- mean(Y,na.rm=TRUE)
   DomRobs <- mean(X,na.rm=TRUE)
   A <- 2*(DomRmod - DomRobs)/(DomRmod + DomRobs)
   out$A <- A

   # Location
   if(is.null(d)) d <- max(object$xdim, na.rm=TRUE)
   num <- centdist(binY,binX)
   L1 <- num/d
   intRamt <- function(id,x) return(sum(x[id$m],na.rm=TRUE))
   RnMod <- as.numeric(unlist(lapply(y, intRamt, x=Y)))
   RnObs <- as.numeric(unlist(lapply(x, intRamt, x=X)))
   xRmodN <- as.numeric(unlist(lapply(y, centdist, y=binY)))
   xRobsN <- as.numeric(unlist(lapply(x, centdist, y=binX)))
   RmodSum <- sum( RnMod, na.rm=TRUE)
   RobsSum <- sum( RnObs, na.rm=TRUE)
   rmod <- sum( RnMod*xRmodN, na.rm=TRUE)/RmodSum
   robs <- sum( RnObs*xRobsN, na.rm=TRUE)/RobsSum
   L2 <- 2*abs(rmod - robs)/d
   out$L <- L1 + L2

   # Structure
   Rmaxer <- function(id, x) return(max(x[id$m], na.rm=TRUE))
   RnMaxMod <- as.numeric(unlist(lapply(y, Rmaxer, x=Y)))
   RnMaxObs <- as.numeric(unlist(lapply(x, Rmaxer, x=X)))
   VmodN <- RnMod/RnMaxMod
   VobsN <- RnObs/RnMaxObs
   Vmod <- sum(RnMod*VmodN, na.rm=TRUE)/RmodSum
   Vobs <- sum(RnObs*VobsN, na.rm=TRUE)/RobsSum
   out$S <- 2*(Vmod - Vobs)/(Vmod + Vobs)
   return(out)
} # end of 'saller' function.

centdist <- function(x,y) {
   xcen <- centroid.owin(x)
   ycen <- centroid.owin(y)
   return(sqrt((xcen$x - ycen$x)^2 + (xcen$y - ycen$y)^2))
} # end of 'centdist' function.

deltamm <- function(x, y = NULL, object=NULL, max.delta=Inf, verbose=FALSE, ...) {
   if( verbose) begin.time <- Sys.time()
   out <- list() 

   if(is.null(y)) {
      Y <- x$Y.feats
      X <- x$X.feats
   } else {
      Y <- y$Y.feats
      X <- y$X.feats
   }
   xdim <- dim( Y[[1]][["m"]])

   m <- length( Y)
   n <- length( X)
   mn <- m*n
   Upsilon <- Psi <- Ksi <- matrix(NA,m,n)

   if( verbose) cat("\n", "Finding Upsilon matrix with delta between each individual object across fields.\n")
   for( i in 1:m) {
	if(verbose) cat(i, "\n")
	for( j in 1:n) {
	   if(verbose) cat(j, " ")
	   Upsilon[i,j] <- deltametric(Y[[i]], X[[j]], ...)
	   if(verbose & (j==n)) cat("\n", "\n")
	} # end of for 'j' loop.
   } # end of for 'i' loop.
	
   if( verbose) cat("\n", "Upsilon found.\n")
   jk <- t( apply( Upsilon, 1, order))
   il <- apply( Upsilon, 2, order)

   if( verbose) cat("\n", "Finding Psi matrix.\n")
   for( i in 1:m) {
	if(verbose) cat(i, "\n")
	k <- jk[i,]
	k <- unique(k)
	for( j in 1:length(k)) {
	   if(verbose) cat(j, " ")
           if( j==1) Xtmp <- X[[k[j]]]
           else Xtmp <- union.owin( Xtmp, X[[k[j]]])
           Psi[i,j] <- deltametric( Y[[i]], Xtmp, ...)
	   if(verbose & (j==n)) cat("\n", "\n")
        } # end of for 'j' loop.
   } # end of for 'i' loop.
   if( verbose) cat("\n", "Psi matrix found.\n")

   if( verbose) cat("\n", "Finding Ksi matrix.\n")
   for(j in 1:n) {
	if(verbose) cat(j, "\n")
	k <- il[,j]
	k <- unique(k)
        for( i in 1:length(k)) {
	   if(verbose) cat(i, " ")
           if( i==1) Ytmp <- Y[[k[i]]]
           else Ytmp <- union.owin( Ytmp, Y[[k[i]]])
           Ksi[i,j] <- deltametric( Ytmp, X[[j]], ...)
	   if(verbose & (i==m)) cat("\n", "\n")
        } # end of for 'i' loop.
   } # end of for 'j' loop.
   if( verbose) cat("\n", "Ksi matrix found.\n")

   mo <- 1:m
   vo <- 1:n
   movo <- cbind( rep(mo,n), rep(vo,each=m))
   movo <- rbind( movo, movo, movo)

   hooklist <- list()
   # Set up first mn components of the list giving the individual objects for each field.
   for(i in 1:mn) {
	tmphook <- list()
	tmphook$f <- movo[i,1]
	tmphook$vx  <- movo[i,2]
	hooklist[[i]] <- tmphook
   } # end of for 'i' loop.

   # Set up hooklist for the Psi matrix containing mergings of forecast field.
   for(i in (mn+1):(2*mn)) {
	tmphook <- list()
	fi <- movo[i,1]
	vi <- movo[i,2]
	tmphook$f <- fi
	tmphook$vx <- jk[fi,1:vi]
	hooklist[[i]] <- tmphook
   } # end of for 'i' loop.

   # Set up hooklist for the Ksi matrix containing mergings of verification field.
   for(i in (2*mn+1):(3*mn)) {
        tmphook <- list()
        fi <- movo[i,1]
        vi <- movo[i,2]
        tmphook$f <- il[1:fi,vi]
        tmphook$vx <- vi
        hooklist[[i]] <- tmphook
   } # end of for 'i' loop.

   bigQ <- array( NA, dim=c(m,n,3))
   Indy <- 1:(3*mn)
   bigQ[,,1] <- Upsilon
   bigQ[,,2] <- Psi
   bigQ[,,3] <- Ksi
   out$Q <- bigQ
   bigQ <- c(bigQ)
   iter <- 1
   mm <- list()
   if( any( bigQ <= max.delta)) {
	   if( verbose) cat("\n", "Looping through Q= {Upsilon,Psi,Ksi} to find estimated best merges and matches.\n")
	   while((length(mo)>0) & (length(vo)>0) & (iter <= 3*mn)) {
	      if(verbose) cat(iter, " ")
	      minQ <- min(bigQ[Indy], na.rm=TRUE)
	      if(is.na(minQ)) {
		warning("deltamm: minimum Q is NA.  Stopping loop now.")
		break
	      }
	      if(minQ > max.delta) break
	      else {
	         id <- Indy[bigQ[Indy]==minQ]
	         if(length(id)>1) id <- id[1]
		 else if(length(id)==0) break
		 fobj <- hooklist[[id]][["f"]]
		 vxobj <- hooklist[[id]][["vx"]]
	         mm[[iter]] <- cbind(fobj, vxobj)
		 colnames(mm[[iter]]) <- NULL
		 findusedIndy <- numeric(0)
	         for(i in Indy) {
			if(any(is.element(hooklist[[i]][["f"]],fobj)) | any(is.element(hooklist[[i]][["vx"]],vxobj))) findusedIndy <- c(findusedIndy, i)
		 } # end of for 'i' loop.
		 findusedIndy <- unique(findusedIndy)
		 Indy <- Indy[!is.element(Indy,findusedIndy)]
		 # bigQ <- bigQ[!is.element(Indy,findusedIndy)]
		 fobj <- unique(fobj)
		 vxobj <- unique(vxobj)
		 mo <- mo[!is.element(mo,fobj)]
		 vo <- vo[!is.element(vo,vxobj)]
	      } # end of if else break bc minQ too big stmts.
	      iter <- iter+1
	   } # end of while mo and vo both still have values loop.
   } # end of if any delta's are small enough stmts.

   if(verbose) cat("\n", "Finished merging and matching.  Finishing up bookeeping.\n")
   ol <- length( mm)
   frunmatched <- mo
   vxunmatched <- vo

   # put things back in the same form as they were entered in case another pass is desired.
   out$X.feats <- out$Y.feats <- list()
   for( i in 1:ol) {
	mo <- unique( mm[[i]][,1])
	vo <- unique( mm[[i]][,2])
	lmo <- length( mo)
	lvo <- length( vo)
	out$Y.feats[[i]] <- Y[[mo[1]]]
	if( lmo > 1) for( j in 2:lmo) out$Y.feats[[i]] <- union.owin( Y[[mo[j]]], out$Y.feats[[i]])
	out$X.feats[[i]] <- X[[vo[1]]]
	if( lvo > 1) for( j in 2:lvo) out$X.feats[[i]] <- union.owin( X[[vo[j]]], out$X.feats[[i]])
   } # end of for 'i' loop.
 
   if( length(frunmatched)>0) {
	newfrunmatched <- (ol+1):(ol+length( frunmatched))
	for( i in newfrunmatched) out$Y.feats[[i]] <- Y[[i]]
   } else newfrunmatched <- numeric(0)

   if( length(vxunmatched)>0) {
	newvxunmatched <- (ol+1):(ol+length( vxunmatched))
	for( i in newvxunmatched) out$X.feats[[i]] <- X[[i]]
   } else newvxunmatched <- numeric(0)
   out$mm.old.labels <- list( mm=matrix( unlist( mm), ncol=2, byrow=TRUE), unmatched=list(fcst=frunmatched, vx=vxunmatched))
   out$mm.new.labels <- list( mm=cbind(1:ol,1:ol), unmatched=list(fcst=newfrunmatched, vx=newvxunmatched))

   Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])
   for( i in 1:length( out$Y.feats)) Ylab[ out$Y.feats[[i]][["m"]]] <- i
   for( j in 1:length( out$X.feats)) Xlab[ out$X.feats[[j]][["m"]]] <- j
   out$X.labeled <- Xlab
   out$Y.labeled <- Ylab
   if(!is.null(object)) {
	out$names <- list(Fcst=object$Fcst.name, Vx=object$Vx.name)
	out$prep <- as.character(substitute(object))
   }
   class(out) <- "deltamm"

   if( verbose) print( Sys.time() - begin.time)
   return( out)
} # end of 'deltamm' function. 

plot.deltamm <- function(x,...) {
   z1 <- x$X.labeled
   z2 <- x$Y.labeled
   if(!is.null(x$names)) {
	z1.name <- paste(x$names$Vx, " Feature Field", sep="")
	z2.name <- paste(x$names$Fcst, " Feature Field", sep="")
   } else {
	z1.name <- "Verification Feature Field"
	z2.name <- "Forecast Feature Field"
   }
   numObjs <- max( c(c(z1),c(z2)), na.rm=TRUE)
   numMatched <- dim( x$mm.new.labels$mm)[1]
   if(length(x$mm.new.labels$unmatched$vx)>0) z1[z1>numMatched] <- numMatched+1
   if(length(x$mm.new.labels$unmatched$fcst)>0) z2[z2>numMatched] <- numMatched+1

   N <- max( c(length(x$X.feats), length(x$Y.feats)), na.rm=TRUE)
   # z1 <- im(z1)
   # z2 <- im(z2)
   par( mfrow=c(1,2))
   image.plot(z1, col=c("grey", rainbow(N)), zlim=c(0,N))
   image.plot(z2, col=c("grey", rainbow(N)), zlim=c(0,N))
   # plot(z1,col=colourmap(c(0,2:(numMatched+2)), range=c(0,numMatched+1)), zlim=c(0,numMatched+1),main=z1.name, ...)
   # plot(z2,col=colourmap(c(0,2:(numMatched+2)), range=c(0,numMatched+1)), zlim=c(0,numMatched+1),main=z2.name, ...)
   invisible()
} # end of 'plot.deltamm' function.

centmatch <- function(x,y=NULL,object=NULL, criteria=1, const=14, verbose=FALSE) {
   if(is.null(y)) {
	xdim <- dim(x$X.labeled)
	Y <- x$Y.feats
	X <- x$X.feats
   } else {
	xdim <- dim(y$X.labeled)
	X <- y$X.feats
	Y <- y$Y.feats
   }
   m <- length(Y)
   n <- length(X)
   if(criteria != 3) {
      Ax <- numeric(n)
      Ay <- numeric(m)
   }

   Dcent  <- matrix(NA,m,n)

   if(verbose) {
	if(criteria != 3) cat("\n", "Looping through each feature in each field to find the centroid differences.\n")
	else cat("\n", "Looping through each feature in each field to find the areas and centroid differences.\n")
   }
   for(i in 1:m) {
	if(verbose) cat(i, "\n")
	if(criteria != 3) {
	   tmpy <- FeatureProps(x=Y[[i]], which.props=c("centroid", "area"))
	   Ay[i] <- sqrt(tmpy$area)
	} else tmpy <- FeatureProps(x=Y[[i]], which.props="centroid")
	ycen <- tmpy$centroid
	for(j in 1:n) {
	   if(verbose) cat(j, " ")
	   if(criteria != 3) {
		tmpx <- FeatureProps(x=X[[j]], which.props=c("centroid", "area"))
	   	Ax[j] <- sqrt(tmpx$area)
	   } else tmpx <- FeatureProps(x=X[[j]], which.props="centroid")
	   xcen <- tmpx$centroid
	   Dcent[i,j] <- sqrt((xcen$x - ycen$x)^2 + (xcen$y - ycen$y)^2)
	} # end of for 'j' loop.
	if(verbose) cat("\n")
   } # end of for 'i' loop.
   if(criteria != 3) {
	Ay <- matrix( rep(Ay,n), m, n)
   	Ax <- matrix( rep(Ax,m), m, n, byrow=TRUE)
   }
   if(criteria == 1) Dcomp <- Ay + Ax
   else if(criteria == 2) Dcomp <- (Ax + Ay)/2
   else if(criteria == 3) Dcomp <- matrix(const,m,n)
   else stop("centmatch: criteria must be 1, 2 or 3.")
   Dcomp <- Dcent < Dcomp

   fmatches <- matrix(NA,m,2)
   if(verbose) cat("Final loops to determine matches.\n")
   for(i in 1:m) {
	if(verbose) cat(i, "\n")
	id <- Dcomp[i,]
	z <- Dcent[i,]
	if(sum(id, na.rm=TRUE) > 1) hold <- (1:n)[z==min(z,na.rm=TRUE)]
	else if(any(id)) hold <- (1:n)[id]
	else hold <- NA
	fmatches[i,] <- c(i, hold)
   } # end of for 'i' loop.
   out <- list()
   out$mm.old.labels <- list()
   mm <- fmatches[!is.na(fmatches[,2]),]
   out$mm.old.labels$mm <- mm
   out$mm.old.labels$unmatched <- list(fcst=(1:m)[is.na(fmatches[,2])], vx=(1:n)[!is.element(1:n,fmatches[,2])])
   ol <- dim(fmatches[!is.na(fmatches[,2]),])[1]
   out$mm.new.labels <- list()
   if(length(ol)>0) {
	out$mm.new.labels$mm <- cbind(1:ol,1:ol)
        if(m>ol) funmatched <- (1:m)[ol:m]
   	else funmatched <- numeric(0)
   	if(n>ol) vxunmatched <- (1:n)[ol:n]
   	else vxunmatched <- numeric(0)
	X.feats <- Y.feats <- list()
        if(ol >= 1) {
          for(i in 1:ol) {
            X.feats[[i]] <- X[[mm[i,2]]]
            Y.feats[[i]] <- Y[[mm[i,1]]]
          } # end of for 'i' loop.
     	} # end of if any matches stmts.
	unm <- length(funmatched)
   	unvx <- length(vxunmatched)
   	if(unm>0) for(i in (ol+1):(ol+unm-1)) Y.feats[[i]] <- Y[[funmatched[i]]]
   	if(unvx>0) for(i in (ol+1):(ol+unvx-1)) X.feats[[i]] <- X[[vxunmatched[i]]]
   } else {
	out$mm.new.labels$mm <- NULL
	X.feats <- X
	Y.feats <- Y
	funmatched <- 1:m
	vxunmatched <- 1:n
   }
   out$mm.new.labels$unmatched <- list(fcst=funmatched, vx=vxunmatched)
   out$X.feats <- X.feats
   out$Y.feats <- Y.feats
   Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])
   for(i in 1:length(X.feats)) Xlab[X.feats[[i]][["m"]]] <- i
   for(i in 1:length(Y.feats)) Ylab[Y.feats[[i]][["m"]]] <- i
   out$X.labeled <- Xlab
   out$Y.labeled <- Ylab
   class(out) <- "centmatch"
   return(out)
} # end of 'centmatch' function.

FeatureProps <- function(x, Im=NULL, which.props=c("centroid", "area", "axis", "intensity"), areafac=1, q=c(0.25, 0.90), ...) {
   out <- list()
   if(is.element("centroid", which.props)) out$centroid <- centroid.owin(x)
   if(is.element("area", which.props)) out$area <- sum(colSums(x$m, na.rm=TRUE), na.rm=TRUE)*areafac
   if(is.element("axis", which.props)) out$axis <- FeatureAxis(x=x,fac=areafac, ...)
   if(is.element("intensity", which.props)) {
	ivec <- matrix(length(q), nrow=1)
	colnames(ivec) <- as.character(q)
	ivec[1,] <- quantile(c(Im[x$m]), probs=q)
	out$intensity <- ivec
   }
   return(out)
} # end of 'FeatureProps' function.

FeatureComps <- function(Y, X, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area", "bdelta", "haus", "ph", "mhd", "med", "msd", "fom", "minsep"),
		sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {

   id1 <- is.element(c("cent.dist", "angle.diff", "area.ratio", "int.area"), which.comps)
   if(any(id1)) {
	list1 <- character(0)
	if( is.element("cent.dist", which.comps)) list1 <- c(list1, "centroid")
	if( any(is.element(c("area.ratio","int.area"), which.comps))) list1 <- c(list1, "area")
	if( is.element("angle.diff", which.comps)) list1 <- c(list1, "axis")
   }
   id2 <- is.element(c("ph", "mhd", "med", "msd", "fom", "minsep"), which.comps)
   if(any(id2)) list2 <- c("ph", "mhd", "med", "msd", "fom", "minsep")[id2]

   if(any(id1)) {
	Xsingle.props <- FeatureProps(x=X, which.props=list1, areafac=sizefac^2)
	Ysingle.props <- FeatureProps(x=Y, which.props=list1, areafac=sizefac^2)
   }

   if(any(id2)) out <- locperf(X=X, Y=Y, which.stats=list2, alpha=alpha, k=k, distfun=distfun, ...)
   else out <- list()

   if(is.element("cent.dist", which.comps)) {
	Xcent.x <- Xsingle.props$centroid$x
	Xcent.y <- Xsingle.props$centroid$y
	Ycent.x <- Ysingle.props$centroid$x
	Ycent.y <- Ysingle.props$centroid$y
	out$cent.dist <- sqrt( (Ycent.x - Xcent.x)^2 + (Ycent.y - Xcent.y)^2)*sizefac
   } # end of if do centroid distance stmt.

   if(is.element("angle.diff", which.comps)) {
	phiX <- Xsingle.props$axis$OrientationAngle$MajorAxis*pi/180
	phiY <- Ysingle.props$axis$OrientationAngle$MajorAxis*pi/180
	out$angle.diff <- abs(atan2(sin(phiX-phiY),cos(phiX-phiY))*180/pi)
   } # end of if do angle difference stmts.

   if(any(is.element(c("area.ratio","int.area"), which.comps))) {
	Xa <- Xsingle.props$area
        Ya <- Ysingle.props$area
   } # end of if do any area calculations stmt.

   if(is.element("area.ratio", which.comps)) out$area.ratio <- min(Xa,Ya)/max(Xa,Ya)

   if(is.element("int.area", which.comps)) {
	denom <- (Xa + Ya)/2
	XY <- intersect.owin(X,Y)
	XYa <- FeatureProps(XY, which.props="area", areafac=sizefac^2)$area
	out$int.area <- XYa/denom
   } # end of if do area ratio stmt.

   if(is.element("bdelta", which.comps)) out$bdelta <- deltametric(X,Y, p=p, c=c)
   if(is.element("haus", which.comps)) out$haus <- deltametric(X,Y,p=Inf,c=Inf)
   return(out)
} # end of 'FeatureComps' function.

FeatureAxis <- function(x, fac=1, flipit=FALSE, twixt=FALSE) {
   out <- list()
   # out$feature.name <- as.character(substitute(x))
   if( flipit) x <- flipxy(x)
   out$z <- x
   ch <- convexhull(x)
   out$chull <- ch
   pts <- unname(cbind(ch$bdry[[1]][["x"]], ch$bdry[[1]][["y"]]))
   out$pts <- pts
   axfit <- sma(y~x, data.frame(x=pts[,1],y=pts[,2]))
   out$axfit <- axfit
   axis.x <- c(axfit$from[[1]], axfit$to[[1]])
   a <- axfit$coef[[1]][1,1]
   b <- axfit$coef[[1]][2,1]
   axis.y <- a + b*axis.x

   if(any(c(is.na(axis.x),is.na(axis.y)))) return(NULL)
   axwin <- owin(xrange=range(axis.x), yrange=range(axis.y))
   MajorAxis <- as.psp(data.frame(x0=axis.x[1], y0=axis.y[1], x1=axis.x[2], y1=axis.y[2]),window=axwin)

   theta <- angles.psp(MajorAxis)
   if((0 <= theta) & (theta <= pi/2)) theta2 <- pi/2 - theta
   else theta2 <- 3*pi/2 - theta
   tmp <- rotate(ch,theta2)
   tmp <- bounding.box(tmp)
   l <- tmp$xrange[2] - tmp$xrange[1]
   theta <- theta*180/pi
   if(twixt) {
	if((theta > 90) & (theta <= 270)) theta <- theta - 180
   	else if((theta > 270) & (theta <= 360)) theta <- theta - 360
   } # end of if force theta to be between +/- 90 degrees.

   MidPoint <- midpoints.psp(MajorAxis)

   r <- lengths.psp(MajorAxis)*fac

   phi <- angles.psp(rotate(MajorAxis,pi/2))

   MinorAxis <- as.psp(data.frame(xmid=MidPoint$x, ymid=MidPoint$y, length=l/fac, angle=phi),window=axwin)
   phi <- phi*180/pi

   out$MajorAxis <- MajorAxis
   out$MinorAxis <- MinorAxis
   out$OrientationAngle <- list(MajorAxis=theta, MinorAxis=phi)
   out$aspect.ratio <- l/r
   out$MidPoint <- MidPoint
   out$lengths <- list(MajorAxis=r, MinorAxis=l)
   out$sma.fit <- axfit
   class(out) <- "FeatureAxis"
   return(out)
} # end of 'FeatureAxis' function.

plot.FeatureAxis <- function(x, ...) {
   args <- list(...)
   par(bg="beige")
   plot( x$z, col="darkblue", main="", ...)
   plot( x$chull, add=TRUE)
   plot( x$MajorAxis, add=TRUE, col="yellow", lwd=1.5)
   plot( x$MajorAxis, add=TRUE, lty=2, lwd=1.5)
   plot( x$MinorAxis, add=TRUE, col="yellow", lwd=1.5)
   plot( x$MinorAxis, add=TRUE, lty=2, lwd=1.5)
   plot( x$MidPoint, add=TRUE, col="darkorange")
   invisible()
} # end of 'plot.FeatureAxis' function.

FeatureMatchAnalyzer <- function(x, y=NULL, matches=NULL, object=NULL,
		which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area", "bdelta", "haus", "ph", "mhd", "med", "msd", "fom", "minsep"),
		sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {
   if(!is.null(matches)) obj <- matches
   else if(!is.null(y)) obj <- y
   else obj <- x
   Yfeats <- obj$Y.feats
   Xfeats <- obj$X.feats
   if(is.null(obj$mm.new.labels$mm)) stop("FeatureAnalyzer: This function requires matches!")
   n <- dim(obj$mm.new.labels$mm)[1]
   out <- list()
   for(i in 1:n) out[[i]] <- FeatureComps(Y=Yfeats[[i]], X=Xfeats[[i]], which.comps=which.comps, sizefac=sizefac, alpha=alpha, k=k, p=p, c=c, distfun=distfun, ...)
   return(out)
} # end of 'FeatureAnalyzer' function.

bearing <- function( point1, point2, deg=TRUE, aty="compass") {
   if( is.null( dim( point1)) & length( point1)==2) point1 <- matrix( point1, 1, 2)
   if( is.null( dim( point2)) & length( point2)==2) point2 <- matrix( point2, 1, 2)
   
   if( deg) {
      # convert latitudes to radians
      point1[,2] <- point1[,2]*pi/180
      point2[,2] <- point2[,2]*pi/180
   } 
   # compute difference in longitude
   # dlon <- (point1[,1] - point2[,1])
   dlon <- (point1[,1] - point2[,1])
   if( deg) dlon <- dlon*pi/180
   S <- cos( point2[,2])*sin( dlon)
   Cval <- cos( point1[,2])*sin( point2[,2]) - sin( point1[,2])*cos( point2[,2])*cos( dlon)
   out <- atan2( S, Cval)
   
   # convert to degrees
   if( deg) out <- out*180/pi
   
   # out[out < 0 ] <-  out[out < 0] + 360
   # out[ out > 360] <- NA
   
   # print( out)
   
   if(aty == "radial") {
       # convert to polar coordinate type angle
       out[ (out >= 0)  &  (out <= 45)]  <- abs(out[ (out >= 0)  &  (out <= 45)] - 45) + 45
       out[ out > 45] <- 90 - out[ out > 45]
       out[ out < 0] <- out[ out < 0] + 360
      } # end of if aty is "radial" stmts.
   return( out)
} # end of 'bearing' function.

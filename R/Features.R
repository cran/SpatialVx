FeatureSuite <- function(object, fun.object, time.point=1, model=1, verbose=FALSE) {
    if( verbose) begin.time <- Sys.time()

    a <- attributes(object)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    xdim <- a$xdim
    Yim <- im(Xhat)
    Xim <- im(X)
    out <- list()
    attributes(out) <- a

    if(verbose) cat("\n", "Calling feature identification function ", fun.object$identify$fun, "\n")
    feats <- do.call(fun.object$identify$fun, c(list(object=object), fun.object$identify$args))
    out$features <- feats

    if( !is.null(fun.object$merge)) {

      if( verbose) cat("\n", "Merging features within fields with ", fun.object$merge$fun, "\n")
      feats2 <- do.call(fun.object$merge$fun, c( list(x=feats, object=object), fun.object$merge$args))
      out$merges <- feats2

    } else feats2 <- NULL

    if(!is.null(fun.object$match)) {

      if( verbose) cat("\n", "Matching features between fields with ", fun.object$match$fun, "\n")
      matches <- do.call(fun.object$match$fun, c(list(x=feats, y=feats2, object=object), fun.object$match$args))
      out$matches <- matches

    } else matches <- NULL

    if( verbose) cat("\n", "Analyzing features with ", fun.object$analysis$fun, "\n")
    out$results <- do.call(fun.object$analysis$fun,
			c(list(x=feats, y=feats2, matches=matches, object=object), fun.object$analysis$args))
    out$fun.obj <- fun.object

    dn <- a$data.name
    if(length(dn) == a$nforecast + 2) dn <- dn[-(1:2)]
    else dn <- dn[-1]

    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    if(length(a$data.name) == a$nforecast + 2) attr(out, "data.name") <- c(a$data.name[1:2], dn[model.num])
    else attr(out, "data.name") <- c(a$data.name[1], dn[model.num])

    if( verbose) print(Sys.time() - begin.time)
    attr(out, "time.point") <- time.point
    attr(out, "model") <- model
    class(out) <- "FeatureSuite"
    return(out)
} # end of 'FeatureSuite' function.

summary.FeatureSuite <- function(object,...) {

    args <- list(...)
    if(is.null(args$silent)) silent <- FALSE
    else silent <- args$silent
    interest <- args$interest
    con <- args$con
    FUN <- object$fun.obj
    PREP <- attributes(object)
    
    if(!silent) {
	if(length(PREP$data.name) == 3) {
	    vxname <- PREP$data.name[2]
	    fcstname <- PREP$data.name[3]
	} else {
	    vxname <- PREP$data.name[1]
            fcstname <- PREP$data.name[2]
	}

        cat("\n", PREP$msg, "\n")
        cat("\n", "Features identified by function: ", FUN$identify$fun, "\n")
        cat("\n", "Number of Features identified in ", vxname, " = ", length(object$features$X.feats), "\n")
        cat("\n", "Number of Features identified in ", fcstname, " = ", length(object$features$Y.feats), "\n")
        cat("\n", "\n")
        if(!is.null(object$merges)) {
	    cat("Objects merged using function: ", FUN$merge$fun, "\n")
	    cat("\n", "Number of Features identified in ", vxname, " after merging = ", length(object$merges$X.feats), "\n")
	    cat("\n", "Number of Features identified in ", fcstname, " after merging = ", length(object$merges$Y.feats), "\n")
        }
        if(!is.null(object$matches)) {
	    cat("\n", "Objects matched via function: ", FUN$match$fun, "\n")
	    cat("\n", "Number of matched objects = ", dim(object$matches$mm.new.labels$mm)[1], "\n")
        }
    } # end of if '!silent' stmts.

    PREP$names <- NULL
    out <- list()
    attributes(out) <- PREP
    out$feature.summary <- summary(object$features, silent=silent)
    out$analysis.summary <- summary(object$results, silent=silent, interest=interest, con=con)
    invisible(out)
} # end of 'summary.FeatureSuite' function.

print.FeatureSuite <- function(x, ...) {
    a <- attributes(x)
    if(!is.null(a$msg)) cat(a$msg)
    if(!is.null(a$data.name)) {
	cat("\n", "Features for:\n")
	print(a$data.name)
	dn <- a$data.name
	if(length(dn)==3) dn <- dn[-1]
        cat(length(x$features$X.feats), " ", dn[1], " features found.\n")
	cat(length(x$features$Y.feats), " ", dn[2], " features found.\n")
    } else {
	cat(length(x$features$X.feats), " verification features found.\n")
        cat(length(x$features$Y.feats), " forecast features found.\n")
    }
    invisible()
} # end of 'print.FeatureSuite' function.

plot.FeatureSuite <- function(x, ...) {
   FUN <- x$fun.obj
   PREP <- attributes(x)

    if(length(PREP$data.name) == 3) {
        vxname <- PREP$data.name[2]
        fcstname <- PREP$data.name[3]
    } else {
        vxname <- PREP$data.name[1]
        fcstname <- PREP$data.name[2]
    }

   # par(mfrow=c(2,2))
   # zl <- range(c(c(X),c(Y)),finite=TRUE)
   # image(X, main=PREP$data.name[1], col=c("grey",tim.colors(64)), zlim=zl, axes=FALSE)
   # image(Y, main=PREP$data.name[2], col=c("grey",tim.colors(64)), zlim=zl, axes=FALSE)
   # image.plot(X, col=c("grey",tim.colors(64)), zlim=zl, legend.only=TRUE)
   if(!is.null(x$matches)) plot(x$matches, no.label=TRUE, plot.set=TRUE)
   else if(!is.null(x$merges)) plot(x$merges)
   else plot(x$features)
   plot(x$results)
} # end of 'plot.FeatureSuite' function.

FeatureFunPrep <- function(identfun=NULL, identfun.args=NULL,
			mergefun=NULL, mergefun.args=NULL,
			matchfun=NULL, matchfun.args=NULL,
			analysisfun=NULL, analysisfun.args=NULL) {
   out <- list()
   if(!is.null(identfun)) out$identify <- list(fun=identfun, args=identfun.args)
   if(!is.null(mergefun)) out$merge <- list(fun=mergefun, args=mergefun.args)
   if(!is.null(matchfun)) out$match <- list(fun=matchfun, args=matchfun.args)
   if(!is.null(analysisfun)) out$analysis <- list(fun=analysisfun, args=analysisfun.args)
   return( out)
} # end of 'FeatureFunPrep' function.

# FeatureSuitePrep <- function(Vx.name, Fcst.name, loc=NULL, units=NULL) {
#    out <- list()
#    data.name <- c(Vx.name, Fcst.name)
#    names(data.name) <- c("verification","forecast")
#    out$data.name <- data.name
#    Y <- get( Fcst.name)
#    X <- get( Vx.name)
#    xdim <- dim( X)
#    out$xdim <- xdim
#    out$loc <- loc
#    out$units <- units
#    class( out) <- "FeaturesSuitePrep"
#    return( out)
# } # end of 'FeatureSuitePrep' function.

convthresh <- function(object, smoothfun="disk2dsmooth", smoothpar=1, smoothfunargs=NULL, thresh=1e-8, idfun="disjointer", zero.down=FALSE, time.point=1, model=1, ...) {
   
    a <- attributes(object)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat
    ## End: Get the data sets

    xdim <- a$xdim

    Xsm <- do.call(smoothfun, c(list(x=X, lambda=smoothpar), smoothfunargs))
    Ysm <- do.call(smoothfun, c(list(x=Y, lambda=smoothpar), smoothfunargs))

    if( zero.down) {
	Xsm[ Xsm < 0] <- 0
	Xsm <- zapsmall( Xsm)
	Ysm[ Ysm < 0] <- 0
	Ysm <- zapsmall( Ysm)
    }

    if( is.null( thresh)) thresh <- c(quantile(c(Xsm), 0.75), quantile(c(Ysm), 0.75))
    else if( length( thresh) == 1) thresh <- c(thresh, thresh)

    sIx <- sIy <- matrix(0, xdim[1], xdim[2])
    sIx[ Xsm >= thresh[1]] <- 1
    sIy[ Ysm >= thresh[2]] <- 1

    X.feats <- do.call(idfun, c( list( x=sIx), list(...)))
    Y.feats <- do.call(idfun, c( list( x=sIy), list(...)))

    Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])

    if(!is.null(X.feats)) for( i in 1:length( X.feats)) Xlab[ X.feats[[i]][["m"]]] <- i
    else X.feats <- NULL

    if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[ Y.feats[[j]][["m"]]] <- j
    else Y.feats <- NULL

    # if(is.null(X.feats)) warning("convthresh: No values above threshold in verification field.")
    # if(is.null(Y.feats)) warning("convthresh: No values above threshold in forecast field.")

    out <- list()
    attributes(out) <- a

    out$X <- X
    out$Xhat <- Y
    out$X.feats <- X.feats
    out$Y.feats <- Y.feats
    out$X.labeled <- Xlab
    out$Y.labeled <- Ylab
    out$identifier.function <- "convthresh"
    out$identifier.label <- "Convolution Threshold"

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
	dn <- a$data.name[-(1:2)]
	vxname <- a$data.name[1:2]
    } else {
	dn <- a$data.name[-1]
	vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    class(out) <- "features"
    return(out)
} # end of 'convthresh' function.

threshsizer <- function(object, thresh=1e-8, Ncontig=50, idfun="disjointer", time.point=1, model=1, ...) {

    a <- attributes(object)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat
    ## End: Get the data sets

    xdim <- a$xdim

    if(length(thresh)==1) thresh <- c(thresh, thresh)

    sIx <- sIy <- matrix(0, xdim[1], xdim[2])
    sIx[ X >= thresh[1]] <- 1
    sIy[ Y >= thresh[2]] <- 1

    X.feats0 <- do.call(idfun, c( list( x=sIx), list(...)))
    Y.feats0 <- do.call(idfun, c( list( x=sIy), list(...)))

    if(length(X.feats0) == 0) X.feats <- NULL
    else X.feats <- list()

    if(length(Y.feats0) == 0) Y.feats <- NULL
    else Y.feats <- list()

    Nfun <- function(Obj) return(sum(colSums(Obj[["m"]], na.rm=TRUE), na.rm=TRUE))
    Xnums <- c(unlist(lapply(X.feats0, Nfun)))
    Ynums <- c(unlist(lapply(Y.feats0, Nfun)))

    Xnums <- Xnums >= Ncontig
    Ynums <- Ynums >= Ncontig

    Xj <- (1:length(Xnums))[Xnums]
    Yj <- (1:length(Ynums))[Ynums]

    if(length(X.feats0) > 0) for(i in 1:length(Xj)) X.feats[[i]] <- X.feats0[[ Xj[i] ]]
    if(length(Y.feats0) > 0) for(i in 1:length(Yj)) Y.feats[[i]] <- Y.feats0[[ Yj[i] ]]

    Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])

    if(!is.null(X.feats)) for( i in 1:length(X.feats)) Xlab[ X.feats[[i]][["m"]]] <- i
    else X.feats <- NULL

    if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[ Y.feats[[j]][["m"]]] <- j
    else Y.feats <- NULL

    # if(is.null(X.feats)) warning("convthresh: No values above threshold in verification field.")
    # if(is.null(Y.feats)) warning("convthresh: No values above threshold in forecast field.")

    out <- list()
    attributes(out) <- a

    out$X <- X
    out$Xhat <- Y
    out$X.feats <- X.feats
    out$Y.feats <- Y.feats
    out$X.labeled <- Xlab
    out$Y.labeled <- Ylab
    out$identifier.function <- "threshsizer"
    out$identifier.label <- paste(Ncontig, " Connected pts > ", thresh, sep="")

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    class(out) <- "features"
    return(out)
} # end of 'threshsizer' function.

summary.features <- function(object,...) {
    args <- list(...)
    a <- attributes(object)

    if(is.null(args$silent)) silent <- FALSE
    else silent <- args$silent

    X <- object$X.feats
    Y <- object$Y.feats

    out <- list()
    b <- a
    b$names <- NULL
    attributes(out) <- b

    n <- length(X)
    m <- length(Y)

    if(!is.null(object$Xhat)) {
	Im1 <- object$X
	Im2 <- object$Xhat
	holdX <- matrix(NA, n, 7)
        holdY <- matrix(NA, m, 7)
        colnames(holdX) <- colnames(holdY) <- c("centroidX", "centroidY", "area", "OrientationAngle", "AspectRatio",
                                                "Intensity0.25", "Intensity0.9")
	wpr <- c("centroid", "area", "axis", "intensity")
	do.int <- TRUE
    } else {
	Im1 <- Im2 <- NULL
	holdX <- matrix(NA, n, 5)
        holdY <- matrix(NA, m, 5)
        colnames(holdX) <- colnames(holdY) <- c("centroidX", "centroidY", "area", "OrientationAngle", "AspectRatio")
	wpr <- c("centroid", "area", "axis")
	do.int <- FALSE
    }

    for(i in 1:n) {

      tmp <- FeatureProps(X[[i]], Im=Im1, which.props=wpr)
      holdX[i,1:2] <- c(tmp$centroid$x, tmp$centroid$y)
      holdX[i,3] <- tmp$area
      if(!is.null(c(tmp$axis$OrientationAngle$MajorAxis))) holdX[i,4] <- c(tmp$axis$OrientationAngle$MajorAxis)
      if(!is.null(tmp$axis$aspect.ratio)) holdX[i,5] <- tmp$axis$aspect.ratio
      if(do.int) holdX[i,6:7] <- c(tmp$intensity)

    }

    for(i in 1:m) {

      tmp <- FeatureProps(Y[[i]], Im=Im2, which.props=wpr)
      holdY[i,1:2] <- c(tmp$centroid$x, tmp$centroid$y)
      holdY[i,3] <- tmp$area
      if(!is.null(c(tmp$axis$OrientationAngle$MajorAxis))) holdY[i,4] <- c(tmp$axis$OrientationAngle$MajorAxis)
      if(!is.null(tmp$axis$aspect.ratio)) holdY[i,5] <- tmp$axis$aspect.ratio
      if(do.int) holdY[i,6:7] <- c(tmp$intensity)

    }

    if(!silent) {
	cat("\n", "Verification field (", object$data.name[1], ") feature properties:\n")
	print(holdX)
	cat("\n", "Forecast field (", object$data.name[2], ") feature properties:\n")
	print(holdY)
    }
    out$X <- holdX
    out$Y <- holdY
    class(out) <- "summary.features"
    invisible(out)
} # end of 'summary.features' function.

print.features <- function(x, ...) {

    a <- attributes(x)

    print(a$msg)
    print(a$data.name)

    if(length(a$data.name) == 3) {
	vxname <- a$data.name[2]
	fcstname <- a$data.name[3]
    } else {
	vxname <- a$data.name[1]
	fcstname <- a$data.name[2]
    }

    print(x$identifier.function)
    print(x$identifier.label)

    cat("\n", length(x$X.feats), " ", vxname, " objects identified.\n")
    cat(length(x$Y.feats), " ", fcstname, " objects identified.\n")
    invisible()
} # end of 'print.features' function.

plot.features <- function(x, ..., set.pw=FALSE, loc.byrow=TRUE, horizontal=TRUE) {

    if(!is.logical(set.pw)) {
	if(is.numeric(set.pw) && length(set.pw) == 2) par(mfrow=set.pw, oma=c(0,0,2,0))
	else stop("plot.features: invalid set.pw argument.")
    } else if(!set.pw) par(oma=c(0,0,2,0))
    else if(set.pw) par(mfrow=c(1,2), oma=c(0,0,2,0))

    a <- attributes(x)

    if(length(a$data.name) == 3) {
        vxname <- a$data.name[2]
        fcstname <- a$data.name[3]
    } else {
        vxname <- a$data.name[1]
        fcstname <- a$data.name[2]
    }

    X <- x$X.labeled
    Y <- x$Y.labeled
    m <- max(c(c(X),c(Y)),na.rm=TRUE)
    zl <- c(0,m)

    domap <- a$map
    proj <- a$projection
    xd <- a$xdim

    if(domap) {
	locr <- apply(a$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))
    }

    if(proj) loc <- list(x=matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    
    if(domap) {

	if(proj) {

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$y)
	    axis(2, at=ax$y, labels=ax$y)
	    poly.image(loc$x, loc$y, X, col=c("white", rainbow(m)), add=TRUE, xaxt="n", yaxt="n", zlim=zl)
	    map(add=TRUE, lwd=1.5)
	    map(add=TRUE, database="state")
	    title(vxname)

	    map(xlim=locr[,1], ylim=locr[,2], type = "n")
	    axis(1, at=ax$x, labels=ax$y)
	    axis(2, at=ax$y, labels=ax$y)
            poly.image(loc$x, loc$y, Y, col=c("white", rainbow(m)), add=TRUE, xaxt="n", yaxt="n", zlim=zl)
            map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
	    title(fcstname)

	} else {

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$y)
	    axis(2, at=ax$y, labels=ax$y)
	    image(as.image(X, nx=xd[1], ny=xd[2], x=a$loc), col=c("white", rainbow(m)), add=TRUE,
		xaxt="n", yaxt="n", zlim=zl)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
	    title(vxname)

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$y)
	    axis(2, at=ax$y, labels=ax$y)
	    image(as.image(Y, nx=xd[1], ny=xd[2], x=a$loc), col=c("white", rainbow(m)), add=TRUE,
		xaxt="n", yaxt="n", zlim=zl)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")

	    title(fcstname)

	}
    } else {

	if(proj) {

	    poly.image(loc$x, loc$y, X, col=c("white", rainbow(m)), main=vxname, xaxt="n", yaxt="n", zlim=zl)
	    poly.image(loc$x, loc$y, Y, col=c("white", rainbow(m)), main=fcstname, xaxt="n", yaxt="n", zlim=zl)

	} else {

            image(X, col=c("white", rainbow(m)), main=vxname, xaxt="n", yaxt="n", zlim=zl)
            image(Y, col=c("white", rainbow(m)), main=fcstname, xaxt="n", yaxt="n", zlim=zl)

	}
    }

    image.plot(Y, col=c("white", rainbow(m)), zlim=zl, legend.only=TRUE, horizontal=horizontal)

    title("")
    mtext(a$msg, line=0.05, outer=TRUE)

    invisible()

} # end of 'plot.features' function.

plot.summary.features <- function(x, ...) {

   X <- x$X
   Y <- x$Y

   n <- dim(X)[1]
   m <- dim(Y)[1]
   k <- dim(X)[2]

   goodX <- sum(!is.na(X[,"OrientationAngle"]),na.rm=TRUE) > 2
   goodY <- sum(!is.na(Y[,"OrientationAngle"]),na.rm=TRUE) > 2

   if(k == 5) par(mfrow=c(2,2))
   else par(mfrow=c(4,2))
   
   hist(X[,"area"], col="darkblue", breaks="FD", main="Verification \nFeature Area", xlab="area")
   hist(Y[,"area"], col="darkorange", breaks="FD", main="Forecast \nFeature Area", xlab="area")

   if(goodX & goodY) {

	plot(X[,"OrientationAngle"], X[,"AspectRatio"], pch=19, col="darkblue",
	    xlab="Major Axis Orientation Angle", ylab="Aspect Ratio")

   	points(Y[,"OrientationAngle"], Y[,"AspectRatio"], pch=19, col="darkorange")

   	legend("topright", legend=c("Verification", "Forecast"), pch=19, col=c("darkblue", "darkorange"), bty="n")
   }

   if(k==7) {
	hist(X[,"Intensity0.25"], col="darkblue", breaks="FD", main="Verification Feature \nIntensity (lower quartile)",
	    xlab="intensity")
	hist(X[,"Intensity0.9"], col="darkblue", breaks="FD", main="Verification Feature \nIntensity (0.9 quartile)",
	    xlab="intensity")
	hist(Y[,"Intensity0.25"], col="darkblue", breaks="FD", main="Forecast Feature \nIntensity (lower quartile)",
	    xlab="intensity")
        hist(Y[,"Intensity0.9"], col="darkblue", breaks="FD", main="Forecast Feature \nIntensity (0.9 quartile)",
	    xlab="intensity")
   }

    xl <- range(c(X[,"centroidX"],Y[,"centroidX"]), finite=TRUE)
    yl <- range(c(X[,"centroidY"],Y[,"centroidY"]), finite=TRUE)

    plot(X[,"centroidX"], X[,"centroidY"], pch=19, col="darkblue", xlim=xl, ylim=yl, xlab="", ylab="",
	main="Feature centroids")
    points(Y[,"centroidX"], Y[,"centroidY"], pch=19, col="darkorange")
    legend("topright", legend=c("Verification", "Forecast"), pch=19, col=c("darkblue", "darkorange"))

   invisible()
} # end of 'plot.summary.features' function.

disjointer <- function(x, method="C") {
   x[ x==0] <- NA
   if(any(!is.na(x))) {
	out <- as.im( x)
   	out <- connected(X=out, method=method)
   	out <- tiles(tess(image=out))
   } else out <- NULL
   return( out)
} # end of 'disjointer' function.

threshfac <- function(object, fac=0.06666667, q=0.95, wash.out=NULL, thresh=NULL, idfun="disjointer",
    time.point=1, model=1, ...) {

    a <- attributes(object)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat
    ## End: Get the data sets

    xdim <- a$xdim
    Ix <- Iy <- matrix(0, xdim[1], xdim[2])

    if(is.null(thresh)) {

	if(is.null(wash.out)){

	   thresh <- quantile(c(X), probs=q)
	   thresh <- c(thresh, quantile(c(Y), probs=q))

	} else {

	   thresh <- quantile(c(X[X >= wash.out]), probs=q)
	   thresh <- c(thresh, quantile(c(Y[Y >= wash.out]), probs=q))

	}
	thresh <- thresh * fac

    } else if(length(thresh)==1) thresh <- c(thresh, thresh)

    Ix[X >= thresh[1]] <- 1
    Iy[Y >= thresh[2]] <- 1

    X.feats <- do.call(idfun, c(list(x=Ix), list(...)))
    Y.feats <- do.call(idfun, c(list(x=Iy), list(...)))

    Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])

    if(!is.null(X.feats)) for( i in 1:length( X.feats)) Xlab[X.feats[[i]][["m"]]] <- i
    else X.feats <- NULL

    if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[Y.feats[[j]][["m"]]] <- j
    else Y.feats <- NULL

    # if(is.null(X.feats)) warning("threshfac: No values above threshold in verification field.")
    # if(is.null(Y.feats)) warning("threshfac: No values above threshold in forecast field.")

    out <- list()
    attributes(out) <- a

    out$X.feats <- X.feats
    out$Y.feats <- Y.feats
    out$X.labeled <- Xlab
    out$Y.labeled <- Ylab
    out$identifier.function <- "threshfac"
    out$identifier.label <- "Threshold"

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    class(out) <- "features"
    return(out)
} # end of 'threshfac' function.

saller <- function(object, x, y=NULL, matches=NULL, d=NULL, time.point=1, model=1) {

    out <- list()
    a <- attributes(object)
    if(!is.null(a$names)) a$names <- NULL
    attributes(out) <- a
    xdim <- a$xdim 

    if(is.null(y)) tmp <- x
    else tmp <- y
    y <- tmp$Y.feats
    x <- tmp$X.feats

    binX <- im(tmp$X.labeled)
    binX <- solutionset(binX > 0)
    binY <- im(tmp$Y.labeled)
    binY <- solutionset(binY > 0)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat
    ## End: Get the data sets

    # Amplitude
    DomRmod <- mean(Y,na.rm=TRUE)
    DomRobs <- mean(X,na.rm=TRUE)
    A <- 2*(DomRmod - DomRobs)/(DomRmod + DomRobs)
    out$A <- A

    # Location
    if(is.null(d)) d <- max(a$xdim, na.rm=TRUE)
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

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model
 
    attr(out, "data.name") <- c(vxname, dn[model.num])

    class(out) <- "saller"
    return(out)
} # end of 'saller' function.

print.saller <- function(x, ...) {
    a <- attributes(x)
    cat(a$msg, "\n")
    print(a$data.name)
    cat("\n\n")
    b <- c(x$S, x$A, x$L)
    names(b) <- c("S", "A", "L")
    print(b)
    invisible(b)
} # end of 'print.saller' function.

summary.saller <- function(object,...) {
   # args <- list(...)
    a <- attributes(object)
    print(a$msg)
    print(a$data.name)

    cat("\n", "Structure Component (S): ", object$S, "\n")
    cat("\n", "Amplitude Component (A): ", object$A, "\n")
    cat("\n", "Location Component (L): ", object$L, "\n")
    invisible()
} # end of 'summary.saller' function.

plot.saller <- function(x, ...) {
   invisible() 
} # end of 'plot.saller' function.

centdist <- function(x,y) {
   xcen <- centroid.owin(x)
   ycen <- centroid.owin(y)
   return(sqrt((xcen$x - ycen$x)^2 + (xcen$y - ycen$y)^2))
} # end of 'centdist' function.

deltamm <- function(object, x, y = NULL, max.delta=Inf, time.point=1, model=1, verbose=FALSE, ...) {

    if( verbose) begin.time <- Sys.time()
    out <- list()
    if(!missing(object)) {
	a <- attributes(object)
	attributes(out) <- a
    }

    if(is.null(y)) {
      # No merging has been performed already.
      Y <- x$Y.feats
      X <- x$X.feats
    } else {
      # For 'FeatureSuite' wrapper function, this is if merging is first performed.
      # Included here in case it is desired to run through this method twice to get
      # better results.
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

    if(!missing(object)) {

        attr(out, "time.point") <- time.point
        attr(out, "model") <- model

        if(length(a$data.name) == a$nforecast + 2) {
            dn <- a$data.name[-(1:2)]
            vxname <- a$data.name[1:2]
        } else {
            dn <- a$data.name[-1]
            vxname <- a$data.name[1]
        }
        if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
        else model.num <- model

        attr(out, "data.name") <- c(vxname, dn[model.num])

    } # end of if '!missing(object)' stmts.

    class(out) <- "matched"

    if( verbose) print( Sys.time() - begin.time)
    return( out)
} # end of 'deltamm' function. 

plot.matched <- function(x,..., set.pw=FALSE, horizontal=TRUE, loc.byrow=TRUE) {

    a <- attributes(x)

    args <- list(...)

    z1 <- x$X.labeled
    z2 <- x$Y.labeled

    if(!is.null(a$data.name)) {

	dn <- a$data.name

	if(length(dn) == 3) {

	    vxname <- dn[2]
	    fcstname <- dn[3]

	} else {

	    vxname <- dn[1]
	    fcstname <- dn[2]

	}

	z1.name <- paste(vxname, "\nFeature Field", sep="")
	z2.name <- paste(fcstname, "\nFeature Field", sep="")

    } else {

	    z1.name <- "Verification\nFeature Field"
	    z2.name <- "Forecast\nFeature Field"

    } # end of if '!is.null(a$data.name)' stmts.

    numMatched <- dim( x$mm.new.labels$mm)[1]

    if(length(x$mm.new.labels$unmatched$vx)>0) z1[z1>numMatched] <- numMatched+1
    if(length(x$mm.new.labels$unmatched$fcst)>0) z2[z2>numMatched] <- numMatched+1 
    N <- max(c(length(x$X.feats), length(x$Y.feats)), na.rm=TRUE)

    if(!is.logical(set.pw) && is.numeric(set.pw)) {

        if(length(set.pw) != 2) stop("plot.matched: invalid set.pw argument.")

	par(mfrow=set.pw, oma=c(0,0,2,0))
    } else if(set.pw) par(mfrow=c(1,2), oma=c(0,0,2,0))
    else if(!is.null(a$msg)) par(oma=c(0,0,2,0))

    if(is.null(a$projection)) proj <- FALSE
    else proj <- a$projection

    if(is.null(a$map)) domap <- FALSE
    else domap <- a$map

    if(!is.null(a$xdim)) xd <- a$xdim
    else xd <- dim(z1)

    if(proj) loc <- list(x=matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    
    if(domap) {
	locr <- apply(a$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))
    }

    if(domap) {

	if(proj) {

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    poly.image(loc$x, loc$y, z1, add=TRUE, col=c("white", rainbow(N)), zlim=c(0,N))
	    map(add=TRUE, lwd=1.5)
	    map(add=TRUE, database="state")
	    title(z1.name)

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    poly.image(loc$x, loc$y, z2, add=TRUE, col=c("white", rainbow(N)), zlim=c(0,N))
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
            title(z2.name)

	} else {

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    image(as.image(z1, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE), col=c("white", rainbow(N)),
		zlim=c(0,N), add=TRUE)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
            title(z1.name)

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    image(as.image(z2, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE), col=c("white", rainbow(N)),
                zlim=c(0,N), add=TRUE)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
            title(z2.name)

	}
    } else {

	if(proj) {

	    poly.image(loc$x, loc$y, z1, add=TRUE, col=c("white", rainbow(N)), zlim=c(0,N))
	    title(z1.name)

	    poly.image(loc$x, loc$y, z2, add=TRUE, col=c("white", rainbow(N)), zlim=c(0,N))
	    title(z2.name)

	} else {

            image(z1, col=c("white", rainbow(N)), zlim=c(0,N), main=z1.name)
            image(z2, col=c("white", rainbow(N)), zlim=c(0,N), main=z2.name)

	}
    }
    image.plot(z1, col=c("white", rainbow(N)), zlim=c(0,N), legend.only=TRUE, horizontal=horizontal)

    if(!is.null(a$msg)) {
	title("")
	mtext(a$msg, line=0.05, outer=TRUE)
    }

    invisible()
} # end of 'plot.matched' function.

# interester <- function(object, which.comps=c("cent.dist","angle.diff","area.ratio","int.area","bdelta","haus","ph","mhd","med","msd","fom","minsep"),
# 			interest=1, sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", angle.thresh=0.8, verbose=FALSE, ...) {
#    out <- list()
#    if(length(interest)==1) interest <- rep(interest, length(which.comps))
#    if(is.null(names(interest))) names(interest) <- which.comps
#    if(verbose) begin.tiid <- Sys.time()
#    if(is.null(object$Y.feats)) {
# 	y <- object$X.feats
# 	out$comparison.type <- "within"
#    } else {
# 	y <- object$Y.feats
# 	out$comparison.type <- "between"
#    }
#    x <- object$X.feats
#    m <- length(y)
#    n <- length(x)
#    res <- matrix(NA, n, m)
#    if(verbose) cat("\n", "Finding interest for ", m*n, " pairs of features.\n")
#    if(is.element("bearing",which.comps)) {
# 	which.comps <- which.comps[which.comps != "bearing"]
# 	warning("interester: bearing is not computed for interest calculations.")
#    }
#    if(is.element("angle.diff",which.comps)) {
# 	angle.conf.x <- numeric(n)
# 	angle.conf.y <- numeric(m)
# 	for( i in 1:n) {
# 	   ar <- FeatureAxis(x[[i]])$aspect.ratio
# 	   if(ar > angle.thresh) angle.conf.x[i] <- 1
# 	} # end of for 'i' loop.
# 	for( j in 1:m) {
# 	   ar <- FeatureAxis(x[[j]])$aspect.ratio
#            if(ar > angle.thresh) angle.conf.y[j] <- 1
# 	} # end of for 'j' loop.
#    }
#    for(i in 1:n) {
# 	if(verbose) cat(i, "\n")
# 	for(j in 1:m) {
# 	   if(verbose) cat(j, " ")
# 	   int <- interest
# 	   if(is.element("angle.diff",which.comps)) int["angle.diff"] <- int["angle.diff"]*angle.conf.x[i]*angle.conf.y[j]
# 	   tmp <- FeatureComps(X=x[[i]], Y=y[[j]], which.comps=which.comps, sizefac=sizefac, alpha = alpha, k = k, p = p, c = c, distfun = distfun, ...)
# 	   res[i,j] <- unlist(tmp)[which.comps]*int[which.comps]
# 	} # end of for 'j' loop.
#    } # end of for 'i' loop.
#    if(verbose) cat("\n")
#    if(verbose) print(Sys.time() - begin.tiid)
#    return(res)
# } # end of 'interester' function.

centmatch <- function(object, x, y=NULL, criteria=1, const=14, time.point=1, model=1, verbose=FALSE) {

    out <- list()
    if(!missing(object)) {
	a <- attributes(object)
	attributes(out) <- a
    }

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

    if(!missing(object)) {
        attr(out, "time.point") <- time.point
        attr(out, "model") <- model

        if(length(a$data.name) == a$nforecast + 2) {
            dn <- a$data.name[-(1:2)]
            vxname <- a$data.name[1:2]
        } else {
            dn <- a$data.name[-1]
            vxname <- a$data.name[1]
        }
        if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
        else model.num <- model

        attr(out, "data.name") <- c(vxname, dn[model.num])
    }

    class(out) <- "matched"
    return(out)
} # end of 'centmatch' function.

FeatureProps <- function(x, Im=NULL, which.props=c("centroid", "area", "axis", "intensity"), areafac=1, q=c(0.25, 0.90), ...) {
    out <- list()
    if(is.element("centroid", which.props)) out$centroid <- centroid.owin(x)
    if(is.element("area", which.props)) out$area <- sum(colSums(x$m, na.rm=TRUE), na.rm=TRUE)*areafac
    if(is.element("axis", which.props)) out$axis <- FeatureAxis(x=x,fac=areafac, ...)
    if(is.element("intensity", which.props)) {
	ivec <- matrix(NA, ncol=length(q), nrow=1)
	colnames(ivec) <- as.character(q)
	ivec[1,] <- quantile(c(Im[x$m]), probs=q)
	out$intensity <- ivec
    }
    return(out)
} # end of 'FeatureProps' function.

FeatureComps <- function(Y, X, which.comps=c("cent.dist", "angle.diff", "area.ratio",
					     "int.area", "bdelta", "haus", "ph", "mhd",
					     "med", "msd", "fom", "minsep", "bearing"),
		sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", deg=TRUE, aty="compass", ...) {

   id1 <- is.element(c("cent.dist", "angle.diff", "area.ratio", "int.area", "bearing"), which.comps)
   if(any(id1)) {
	list1 <- character(0)
	if( any(is.element(c("cent.dist","bearing"), which.comps))) list1 <- c(list1, "centroid")
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

   if(is.element("bearing", which.comps)) out$bearing <- bearing(cbind(Ysingle.props$centroid$x,Ysingle.props$centroid$y),
									cbind(Xsingle.props$centroid$x,Xsingle.props$centroid$y),
								deg=deg, aty=aty)

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

summary.FeatureAxis <- function(object, ...) {
   cat("\n", "Mid-point of Axis is at:\n")
   print(paste("(", object$MidPoint$x, ", ", object$MidPoint$y, ")", sep=""))
   cat("\n", "Major Axis length = ", object$lengths$MajorAxis, "\n")
   cat("\n", "Major Axis Angle = ", object$OrientationAngle$MajorAxis, " degrees\n")
   cat("\n", "Minor Axis length = ", object$lengths$MinorAxis, "\n")
   cat("\n", "Minor Axis Angle = ", object$OrientationAngle$MinorAxis, " degrees\n")
   cat("\n", "Aspect ratio = ", object$aspect.ratio, "\n")
   cat("\n\n", "sma fit summary (see help file for function sma from package smatr)\n\n")
   print(summary(object$sma.fit))
   invisible()
} # end of 'summary.FeatureAxis' function.

FeatureMatchAnalyzer <- function(x, y=NULL, matches=NULL, object=NULL,
    which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
		    "bdelta", "haus", "ph", "mhd", "med", "msd", "fom", "minsep", "bearing"),
    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {

    if(!is.null(matches)) obj <- matches
    else if(!is.null(y)) obj <- y
    else obj <- x

    Yfeats <- obj$Y.feats
    Xfeats <- obj$X.feats

    if(is.null(obj$mm.new.labels$mm)) {
	# stop("FeatureAnalyzer: This function requires matches!")
	out <- "No matches found"
    } else {
        n <- dim(obj$mm.new.labels$mm)[1]

        out <- list()
        for(i in 1:n) out[[i]] <- FeatureComps(Y=Yfeats[[i]], X=Xfeats[[i]], which.comps=which.comps,
					sizefac=sizefac, alpha=alpha, k=k, p=p, c=c, distfun=distfun, ...)
    }

    class(out) <- "FeatureMatchAnalyzer"
    return(out)
} # end of 'FeatureMatchAnalyzer' function.

print.FeatureMatchAnalyzer <- function(x, ...) {
    if(is.list(x)) {
        n <- length(x)
        hold <- x[[1]]
        m <- length(hold)
        cn <- names(hold)

        for(i in 1:n) {
	    l <- unlist(lapply(x[[i]], length))
	    if(any(l==0)) x[[i]][[(1:m)[l==0]]] <- NA
        } # end of for 'i' loop.

        out <- c(unlist(x))
        attributes(out) <- NULL
        out <- matrix(out, n, m, byrow=TRUE)
        colnames(out) <- cn
        print(out)
        invisible(out)
    } else {
	attributes(x) <- NULL
	print(c(x))
	invisible()
    }
} # end of 'print.FeatureMatchAnalyzer' function.

summary.FeatureMatchAnalyzer <- function(object, ...) {

    if(is.list(object)) {

        args <- list(...)
        if(is.null(args$silent)) silent <- FALSE
        else silent <- args$silent
        n <- length(object)
        m <- length(object[[1]])
        res <- matrix(NA, n, m)
        colnames(res) <- names(object[[1]])

        for( i in 1:m) {

	    if(!silent) {

   	        if(names(object[[1]])[i] == "ph") cat("\n", "Partial Hausdorff Distance:\n")
   	        else if(names(object[[1]])[i] == "mhd") cat("\n", "Modified Hausdorff Distance:\n")
   	        else if(names(object[[1]])[i] == "med") cat("\n", "Mean Error Distance:\n")
   	        else if(names(object[[1]])[i] == "msd") cat("\n", "Mean Square Error Distance:\n")
   	        else if(names(object[[1]])[i] == "fom") cat("\n", "Pratt\'s Figure of Merit:\n")
   	        else if(names(object[[1]])[i] == "minsep") cat("\n", "Minimum Separation Distance:\n")
   	        else if(names(object[[1]])[i] == "cent.dist") cat("\n", "Centroid Distance:\n")
   	        else if(names(object[[1]])[i] == "angle.diff") cat("\n", "Angle Difference:\n")
   	        else if(names(object[[1]])[i] == "area.ratio") cat("\n", "Area Ratio:\n")
   	        else if(names(object[[1]])[i] == "int.area") cat("\n", "Intersection Area:\n")
   	        else if(names(object[[1]])[i] == "bdelta") cat("\n", "Baddeley\'s Delta Metric:\n")
   	        else if(names(object[[1]])[i] == "haus") cat("\n", "Hausdorff Distance:\n")
   	        else if(names(object[[1]])[i] == "bearing") cat("\n", "Bearing:\n")

	    } # end of if '!silent' stmt.

	    hold <- numeric(n)+NA
	    for(j in 1:n) if(length(object[[j]][[i]])>0) hold[j] <- object[[j]][[i]]
	    res[,i] <- hold
	    if(!silent) print(hold)
       } # end of for 'i' loop.

   if(!is.null(args$interest)) {

	int <- args$interest
	a <- matrix(int, n, m, byrow=TRUE)

	if(!is.null(args$con)) {
	   con <- args$con
	   con <- match.fun(con)
	   a <- con(res, a, which.comps=names(object[[1]]))
	}

	out <- list()
	out$match.properties <- res
	b <- rowSums(a*res, na.rm=TRUE)
	out$object.interest <- b 
	out$interest.values <- int
	if(!silent) {
	   cat("\n", "Matched object interest values.\n")
	   print(b)
	}
   } else out <- res
       return(invisible(out))
    } else print(object)
    invisible()
} # end of 'summary.FeatureMatchAnalyzer' function.

plot.FeatureMatchAnalyzer <- function(x, ..., type=c("all", "ph", "mhd", "med", "msd", "fom", "minsep",
					"cent.dist", "angle.diff", "area.ratio", "int.area", "bearing",
					"bdelta", "haus"), set.pw=FALSE) {

    if(is.list(x)) {
        type <- tolower(type)
        type <- match.arg(type)

        a <- attributes(x)
        y <- summary(x, silent=TRUE)
        n <- dim(y)[2]

	if(is.null(n)) {
	    n <- 1
	    y <- matrix(y, ncol=1)
	}

        if(type == "all") {
            if((n > 1) && set.pw) {
    	        if(n==2) par(mfrow=c(1,2),oma=c(0,0,2,0))
    	        else if(n==3 | n==4) par(mfrow=c(2,2),oma=c(0,0,2,0))
    	        else if(n==5 | n==6) par(mfrow=c(2,3),oma=c(0,0,2,0))
    	        else if(n==7 | n==8) par(mfrow=c(2,4),oma=c(0,0,2,0))
    	        else if(n==9) par(mfrow=c(3,3),oma=c(0,0,2,0))
    	        else if(n==10 | n==11 | n==12) par(mfrow=c(3,4),oma=c(0,0,2,0))
    	        else par(mfrow=c(4,4),oma=c(0,0,2,0))
            } else par(oma=c(0,0,2,0))
    
            for(i in 1:n) {
                if(is.element(colnames(y)[i], c("angle.diff","bearing"))) {
                    if(colnames(y)[i] == "angle.diff") t1 <- "Angle Difference"
                    else t1 <- "Bearing from Xhat to X\n obj. centroids (ref. = north)"
                    circ.plot(rad(y[,i]), main=t1, shrink=1.5)
    	        } else {
                    if(colnames(y)[i] == "ph") t1 <- "Partial Hausdorff \nDistance"
                    else if(colnames(y)[i] == "mhd") t1 <- "Modified Hausdorff \nDistance"
                    else if(colnames(y)[i] == "med") t1 <- "Mean Error Distance"
                    else if(colnames(y)[i] == "msd") t1 <- "Mean Square Error \nDistance"
                    else if(colnames(y)[i] == "fom") t1 <- "Pratt\'s Figure \nof Merit"
                    else if(colnames(y)[i] == "minsep") t1 <- "Minimum Separation \nDistance"
                    else if(colnames(y)[i] == "cent.dist") t1 <- "Centroid Distance"
                    else if(colnames(y)[i] == "area.ratio") t1 <- "Area Ratio"
                    else if(colnames(y)[i] == "int.area") t1 <- "Intersection Area"
                    else if(colnames(y)[i] == "bdelta") t1 <- "Baddeley\'s Delta Metric"
                    else if(colnames(y)[i] == "haus") t1 <- "Hausdorff Distance"
                    barplot(y[,i], main=t1, ...)
    	        }
            } # end of for 'i' loop.
        } else {
	    if(type=="ph") t1 <- "Partial Hausdorff \nDistance"
            else if(type == "mhd") t1 <- "Modified Hausdorff \nDistance"
            else if(type == "med") t1 <- "Mean Error Distance"
            else if(type == "msd") t1 <- "Mean Square Error \nDistance"
            else if(type == "fom") t1 <- "Pratt\'s Figure \nof Merit"
            else if(type == "minsep") t1 <- "Minimum Separation \nDistance"
            else if(type == "cent.dist") t1 <- "Centroid Distance"
	    else if(type == "angle.diff") t1 <- "Angle Difference"
            else if(type == "area.ratio") t1 <- "Area Ratio"
            else if(type == "int.area") t1 <- "Intersection Area"
	    else if(type == "bearing") t1 <- "Bearing from Xhat to X\n obj. centroids (ref. = north)"
            else if(type == "bdelta") t1 <- "Baddeley\'s Delta Metric"
            else if(type == "haus") t1 <- "Hausdorff Distance"

	    i <- (1:ncol(y))[colnames(y) == type]

	    if(is.element(type, c("angle.diff", "bearing"))) circ.plot(rad(y[,i]), main=t1, shrink=1.5)
	    else barplot(y[,i], main=t1, ...)
        }

            if(!is.null(a$msg) && set.pw) mtext(a$msg, line=0.05, outer=TRUE)
    } # end of if 'is.list(x)' stmts. 
    invisible()

} # end of 'plot.FeatureMatchAnalyzer' function.

bearing <- function(point1, point2, deg=TRUE, aty="compass") {

   if(is.null( dim( point1)) & length( point1)==2) point1 <- matrix(point1, 1, 2)
   if(is.null( dim( point2)) & length( point2)==2) point2 <- matrix(point2, 1, 2)
   
   if(deg) {

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

# compositer <- function(x, ...) {
#     UseMethod("compositer", x)
# } # end of 'compositer' function.
# 
# compositer.SpatialVx <- function(x, ..., time.point=NULL, model=1, identfun="threshsizer", verbose=FALSE) {
# 
#     if(verbose) begin.tiid <- Sys.time()
# 
#     theCall <- match.call()
# 
#     a <- attributes(x)
# 
#     ## Begin: Get the data sets
# #     if(!is.null(time.point)) {
# #         if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
# #         else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
# #         else if(!missing(model)) dat <- datagrabber(object, model=model)
# #         else dat <- datagrabber(object)
# # 
# #         X <- dat$X
# #         Xhat <- dat$Xhat
# #     } else {
# # 	X <- x[[1]]
# # 	if(is.list(x[[2]])) {
# # 	    if(is.function(model)) Xhat <- do.call(model, c(list(x[[2]]), list(...)))
# # 	    else if(length(model) == 1) Xhat <- x[[2]][[model]]
# # 	    else stop("compositer.SpatialVx: invalid model argument.")
# # 	} else Xhat <- x[[2]]
# #     }
#     ## End: Get the data sets
# 
#     ## Internal function to apply compositing strategy to each
#     ## forecast model and to verification field(s).
#     cfun <- function(Z, ..., ifun, time.point, model) {
# 
# 	o <- do.call(ifun, c(list(object=Z, time.point=time.point, model=model), list(...)))
# 
# 	ff <- function(x, obj) {
# 	    y <- c(obj)
# 	    z <- as.logical(c(x$m))
# 	    y[ !z ] <- 0
# 	    return(y)
# 	} # end of internal-internal 'ff' function.
# 
# 	X <- lapply(o$X.feats, ff, obj=o$X)
# 	Xhat <- lapply(o$Y.feats, ff, obj=o$Xhat)
# 
# 	nobj1 <- length(X)
# 	nobj2 <- length(Xhat)
# 
# 	out <- list()
# 	out$X <- matrix(unlist(X), nrow=nobj1, byrow=TRUE)
# 	out$Xhat <- matrix(unlist(Xhat), nrow=nobj2, byrow=TRUE)
# 
# 	return(out)
#     } # end of internal 'cfun' function.
# 
#     out <- list()
#     attributes(out) <- a
# 
#     loop.through.time <- (length(a$xdim) > 2) && (is.null(time.point) || length(time.point) > 1)
#     if(loop.through.time) {
# 	if(is.character(time.point)) tlab <- time.point
#         else tlab <- NULL
#     }    
# 
#     if(is.null(time.point) && length(a$xdim) == 3) time.point <- 1:a$xdim[3]
#     else if(is.null(time.point) && length(a$xdim) == 2) time.point <- 1
#     if(!is.numeric(time.point)) time.point <- (1:a$xdim[3])[time.point == a$time]
# 
#     if(!loop.through.time) res <- cfun(Z=x, ..., time.point=time.point, model=model, ifun=identfun)
#     else {
# 
# 	nvx <- dim(X)[3]
# 	res <- list()
# 	res$X <- numeric(0)
# 	res$Xhat <- numeric(0)
# 
# 	if(verbose) cat("\n", "Looping through time.\n")
# 	for(i in 1:nvx) {
# 	    if(verbose) cat(time.point[i], "\n")
# 	    tmp <- cfun(x, ..., fun=identfun, time.point=time.point[i], model=model)
# 	    X <- tmp$X
# 	    Xhat <- tmp$Xhat
# 
# 	    n1 <- dim(X)
# 	    if(is.null(n1)) n1 <- 1
# 
# 	    n2 <- dim(Xhat)
# 	    if(is.null(n2)) n2 <- 1
# 
# 	    if(!is.null(X)) rownames(X) <- paste(tlab[i], ".obj", 1:n1, sep="") 
# 	    if(!is.null(Xhat)) rownames(Xhat) <- paste(tlab[i], ".obj", 1:n2, sep="")
# 
# 	    res$X <- rbind(res$X, X)
# 	    res$Xhat <- rbind(res$Xhat, Xhat)
# 
# 	} # end of for 'i' loop.
# 
#     } # end of if else more than one time point stmts.
# 
#     # TO DO: Should now have matrices with rows corresponding to unique objects in each
#     # field.  Need to check this.  Next, need to center each object onto a new relative
#     # grid.  Will not worry about the relative grid size here.  Will allow for that in
#     # subsequent functions.
#     # out$X.composites <- res$X
#     # out$Xhat.composites <- res$Xhat
# 
#     attr(out, "time.point") <- time.point
#     attr(out, "model") <- model
#     attr(out, "call") <- theCall
# 
#     attr(out, "identifier.function") <- identfun
# 
#     if(verbose) print(Sys.time() - begin.tiid)
# 
#     class(out) <- "compositer"
#     return(out)
# } # end of 'compositer.SpatialVx' function.
# 
# compositer.default <- function(x, ..., xhat, identfun="threshsizer") {
#     tmp <- make.SpatialVx(x, xhat, data.name=c(deparse(substitute(x)), deparse(substitute(xhat))))
#     res <- compositer.SpatialVx(tmp, ..., identfun=identfun)
#     return(res)
# } # end of 'compositer.default' function.

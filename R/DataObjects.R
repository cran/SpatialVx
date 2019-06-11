make.SpatialVx <- function(X, Xhat, thresholds=NULL, loc=NULL, projection=FALSE, subset=NULL, time.vals = NULL,
			    reg.grid = TRUE, map = FALSE, loc.byrow = FALSE, field.type="", units="", data.name = "",
			    obs.name = "X", model.name = "Xhat",
			    q=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95), qs=NULL) {

    if( is.list( X ) ) {

	nobs <- length( X )
	xdim <- dim( X[[ 1 ]] )

    } else {

	nobs <- 1
	xdim <- dim( X )

    }

    if(is.list(Xhat)) {

	nforecast <- length(Xhat)
	ydim <- dim(Xhat[[1]])

    } else {

	nforecast <- 1
	ydim <- dim( Xhat )

    }

    if(all(xdim != ydim)) {

	stopmsg <- paste("make.SpatialVx: dim of X (", xdim[ 1 ], " by ", xdim[ 2 ],
		") must be the same as dim of (each component of) Xhat (", ydim[ 1 ], " by ",
		ydim[ 2 ], ")", sep = "" )

	stop( stopmsg )
	# stop("make.SpatialVx: dim of X must be the same as dim of (each component of) Xhat")

    } # end of stop bc dims are not equal stmt.

    out <- list(X, Xhat)

    if( is.null( thresholds ) ) {

	if( nobs > 1 ) {

	    othresh <- numeric(0)
	    for( i in 1:nobs ) othresh <- cbind( othresh, quantile( c( X[[ i ]] ), probs = q, na.rm = TRUE ) )

	} else othresh <- quantile( c( X ), probs = q, na.rm = TRUE )

	if( nforecast > 1 ) {

	    fthresh <- numeric(0)
	
	    for( i in 1:nforecast ) cbind( fthresh, quantile( c( Xhat[[ i ]] ), probs = q, na.rm = TRUE ) )

	} else fthresh <- quantile( c( Xhat ), probs = q, na.rm = TRUE )

	qs <- as.character(q)
        # qs <- paste( as.character(q), "quantile" )
	# rownames( othresh ) <- qs
	# rownames( fthresh ) <- qs

	thresholds <- list( X = othresh, Xhat = fthresh )

    } else if( is.list( thresholds ) ) {

	threshnames <- names( thresholds )
	if( !( "X" %in% threshnames ) || !( "Xhat" %in% threshnames ) ) {

	    stop( "make.SpatialVx: invalid thresholds argument.  List must have components named X and Xhat." )

	} 

	odim <- dim( thresholds$X )
 	fdim <- dim( thresholds$Xhat )

	if( odim[ 1 ] != fdim[ 1 ] ) stop( "make.SpatialVx: invalid thresholds argument.  X must have same number of thresholds (rows) as Xhat." )

	nth <- odim[ 1 ]

	if( is.null( odim ) ) thresholds$X <- matrix( thresholds$X, length( thresholds$X ), nobs )
	else if( odim[ 2 ] == 1 && nobs > 1 ) thresholds$X <- matrix( thresholds$X, odim[ 1 ], nobs )
	else if( odim[ 2 ] > 1 && odim[ 2 ] != nobs ) stop( "make.SpatialVx: invalid thresholds argument." )

	if( is.null( fdim ) ) thresholds$Xhat <- matrix( thresholds$Xhat, length( thresholds$Xhat ), nforecast )
        else if( fdim[ 2 ] == 1 && nforecast > 1 ) thresholds$X <- matrix( thresholds$Xhat, fdim[ 1 ], nforecast )
        else if( fdim[ 2 ] > 1 && fdim[ 2 ] != nobs ) stop( "make.SpatialVx: invalid thresholds argument." )

	qs <- paste("threshold", 1:nth )

    } else if( is.vector(thresholds) ) {

        qs <- as.character(thresholds)
        # thresholds <- cbind(thresholds, thresholds)
        nth <- length( thresholds )
        thresholds <- list( X = matrix( thresholds, nth, nobs ),
                            Xhat = matrix( thresholds, nth, nforecast ) )

    } else stop( "make.SpatialVx: invalid thresholds argument.  Must be a vector or list." )

    # udim <- dim(thresholds)
    # if( udim[2] == 2) colnames(thresholds) <- c("X", "Xhat")
    # else if( udim[2] > 2) colnames(thresholds) <- c("X", paste("Xhat", 1:(udim[2] - 1), sep=""))

    # rownames(thresholds) <- qs

    if( is.null( loc ) ) loc <- cbind( rep( 1:xdim[ 1 ], xdim[ 2 ] ), rep( 1:xdim[ 2 ], each = xdim[ 1 ] ) )

    if(is.null(time.vals)) {

	if(length(xdim) == 3) time.vals <- 1:xdim[3]
	else time.vals <- 1

    }

    if( !missing( data.name ) ) msg <- data.name
    else msg <- ""

    if(field.type != "" && units != "") msg <- paste(msg, "\n", field.type, " (", units, ")", sep="")
    else if(field.type != "") msg <- paste(msg, "\n", field.type, sep="")
    else if(units != "") msg <- paste(msg, "\n", " (", units, ")", sep="")


    class(out) <- "SpatialVx"
    attr(out, "xdim") <- xdim
    # attr(out, "udim") <- udim
    attr(out, "time") <- time.vals
    attr(out, "thresholds") <- thresholds
    attr(out, "loc") <- loc
    attr(out, "loc.byrow") <- loc.byrow
    attr(out, "subset") <- subset
    attr(out, "data.name") <- data.name
    attr(out, "obs.name" ) <- obs.name
    attr(out, "model.name" ) <- model.name
    attr(out, "nobs" ) <- nobs
    attr(out, "nforecast") <- nforecast
    attr(out, "field.type") <- field.type
    attr(out, "units") <- units
    attr(out, "projection") <- projection
    attr(out, "reg.grid") <- reg.grid
    attr(out, "map") <- map
    attr(out, "qs") <- qs
    attr(out, "msg") <- msg

    return(out)

} # end of 'make.SpatialVx' function.

print.SpatialVx <- function(x, ...) {
    tmp <- attributes(x)
    cat("SpatialVx data object: ", tmp$data.name, "\n")
    cat(tmp$msg, "\n")

    cat("\n\nObservations: ", tmp$obs.name, "\n\n" )
    cat("\n\nModels: ", tmp$model.name, "\n\n" )
    # if(tmp$field.type != "") cat("Field type: ", tmp$field.type)
    # if(tmp$units != "") cat(" (", tmp$units, ")", "\n")

    cat("Field dimensions: ", tmp$xdim, "\n")

    if(tmp$nforecast > 1) cat(tmp$nforecast, " forecasts to be evaluated/verified.\n")
    else cat("1 forecast to be evaluated/verified.\n")

    u <- tmp$thresholds
    cat("\nThresholds\n")
    print(u)

    if(!is.null(tmp$levels)) {
        cat("\n", "Neighborhood levels to be applied for neighborhood methods:\n")
        print(tmp$levels)
        cat("\n", "Smoothing function to be applied for neighborhood methods:\n")
        print(tmp$smooth.fun)
        if(!is.null(tmp$smooth.params)) {
            cat("\n", "Smoothing function parameters:\n")
            print(tmp$smooth.params)
        }
    }

    invisible()

} # end of 'print.SpatialVx' function.

summary.SpatialVx <- function(object, ...) {

    out <- list()

    x <- object
    print(x)
    tmp <- attributes(x)

    print( out$X <- summary(x[[1]]))

    if(tmp$nforecast > 1) print(out$Xhat <- lapply(x[[2]], summary))
    else print(out$Xhat <- summary(x[[2]]))

    invisible(out)
} # end of 'summary.SpatialVx' function.

plot.SpatialVx <- function( x, ..., time.point = 1, obs = 1, model = 1, col = c( "gray", tim.colors( 64 ) ), zlim, mfrow = c(1, 2) ) {

    a <- attributes( x )

    # if( a$map ) class( x ) <- "SpatialVxMap"
    # else class( x ) <- "SpatialVxNoMap"

    if( a$map ) {

	if( !missing( zlim ) ) plot.SpatialVxMap( x = x, ..., time.point = time.point, obs = obs, model = model, col = col, zlim = zlim, mfrow = mfrow )
	else plot.SpatialVxMap( x = x, ..., time.point = time.point, obs = obs, model = model, col = col, mfrow = mfrow )

    } else {

	if( !missing( zlim ) ) plot.SpatialVxNoMap( x = x, ..., time.point = time.point, obs = obs, model = model, col = col, zlim = zlim, mfrow = mfrow )
	else plot.SpatialVxNoMap(  x = x, ..., time.point = time.point, obs = obs, model = model, col = col, mfrow = mfrow )

    }

    # UseMethod( "plot", x )

}

plot.SpatialVxMap <- function( x, ..., time.point = 1, obs = 1, model = 1, col = c( "gray", tim.colors( 64 ) ), zlim, mfrow = c(1, 2) ) {

    if( !is.null( mfrow ) ) {

        op <- par()

        par( mfrow = mfrow, oma = c(0, 0, 2, 0) )

    }

    a <- attributes( x )

    loc <- a$loc
    xdim <- a$xdim

    l <- list( x = matrix( loc[, 1 ], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ),
		y = matrix( loc[, 2 ], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ) )

    lr <- apply( loc, 2, range, finite = TRUE )

    dat <- datagrabber( x, ..., time.point = time.point, obs = obs, model = model )
    X <- dat$X
    Xhat <- dat$Xhat

    if( missing( zlim ) ) zlim <- range( c( c(X), c(Xhat) ), finite = TRUE )

    if( is.function( time.point ) ) tp <- as.character( substitute( time.point ) )
    else tp <- time.point

    msg <- a$data.name
    msgX <- paste( a$obs.name[ obs ], "\n", "(", a$field.type, ", ", a$units, ", time = ", tp, ")", sep = "" )
    msgXhat <- paste( a$model.name[ model ], "\n", "(", a$field.type, ", ", a$units, ", time = ", tp, ")", sep = "" )

    map( xlim = lr[,1], ylim = lr[,2], type = "n" )
    title( msgX )
    poly.image( l$x, l$y, X, col = col, zlim = zlim, add = TRUE, ... )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    map.axes()

    map( xlim = lr[,1], ylim = lr[,2], type = "n" )
    title( msgXhat )
    poly.image( l$x, l$y, Xhat, col = col, zlim = zlim, add = TRUE, ... )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    map.axes()

    image.plot( X, zlim = zlim, col = col, legend.only = TRUE, ... )

    if( !is.null( mfrow ) ) {

	title("")
	mtext( msg, line = 0.5, outer = TRUE )

        par( mfrow = op$mfrow )

    }

    invisible()

} # end of 'plot.SpatialVxMap' function.

plot.SpatialVxNoMap <- function( x, ..., time.point = 1, obs = 1, model = 1, col = c( "gray", tim.colors( 64 ) ), zlim, mfrow = c(1, 2) ) {

    if( !is.null( mfrow ) ) {

        op <- par()

        par( mfrow = mfrow, oma = c(0, 0, 2, 0) )

    }

    a <- attributes( x )

    dat <- datagrabber( x, ..., time.point = time.point, obs = obs, model = model )
    X <- dat$X
    Xhat <- dat$Xhat

    if( missing( zlim ) ) zlim <- range( c( c(X), c(Xhat) ), finite = TRUE )

    if( is.function( time.point ) ) tp <- as.character( substitute( time.point ) )
    else tp <- time.point

    msg <- a$data.name
    msgX <- paste( a$obs.name[ obs ], "\n", "(", a$field.type, ", ", a$units, ", time = ", tp, ")", sep = "" )
    msgXhat <- paste( a$model.name[ model ], "\n", "(", a$field.type, ", ", a$units, ", time = ", tp, ")", sep = "" )

    image( X, col = col, zlim = zlim, main = msgX, ... )
    image( Xhat, col = col, zlim = zlim, main = msgXhat, ... )
    image.plot( X, zlim = zlim, col = col, legend.only = TRUE, ... )

    if( !is.null( mfrow ) ) {

	title("")
        mtext( msg, line = 0.5, outer = TRUE )

        par( mfrow = op$mfrow )

    }

    invisible()

} # end of 'plot.SpatialVxNoMap' function.

datagrabber.SpatialVx <- function(x, ..., time.point = 1, obs = 1, model = 1) {

    tmp <- attributes(x)

    if(!missing(time.point)) {

	if( !is.function( time.point ) ) {

            if(length(time.point) > 1) stop("datagrabber: length of time.point must be one.")
            tiid <- tmp$time
            Nt <- length(tiid)
            if(!is.numeric(time.point)) time.point <- (1:Nt)[ tiid == time.point ]
            if(is.na(time.point)) stop("datagrabber: invalid time.point argument.")

	}

    }

    if( !is.numeric( obs ) ) stop( "datagrabber: invalid obs argument.  Must be numeric." )
    if( !is.numeric( model ) ) stop( "datagrabber: invalid model argument.  Must be numeric." )

    xdim <- tmp$xdim
    nobs <- tmp$nobs
    nf <- tmp$nforecast

    Vx <- x[[1]]
    if( nobs > 1 ) Vx <- Vx[[ obs ]]
    
    Fcst <- x[[2]]
    if(nf > 1) Fcst <- Fcst[[ model ]]

    if(length(xdim) == 3) {

	if( is.function( time.point ) ) {

	    afun <- match.fun( time.point )

	    Vx <- apply( Vx, 1:2, afun, ... )
	    Fcst <- apply( Fcst, 1:2, afun, ... )

	} else {

	    Vx <- Vx[,,time.point]
	    Fcst <- Fcst[,,time.point]

	} # end of if take a single time point or a function over time.

    } # if more than one time point present stmts.

    out <- list( X = Vx, Xhat = Fcst )

    return( out ) 

} # end of 'datagrabber.SpatialVx' function.

datagrabber.features <- function(x, ...) {

    return(list(X = x$X, Xhat = x$Xhat))

} # end of 'data.grabber.features' function.

datagrabber.matched <- function(x, ...) {

    return(list(X = x$X, Xhat = x$Xhat))

} # end of 'data.grabber.matched' function.

thresholder <- function( x, type = c( "binary", "replace.below" ), th, rule = ">=", replace.with = 0, ... ) {

    if( !is.element( rule, c( ">=", ">", "<", "<=" ) ) ) stop("thresholder: invalid rule argument.")

    UseMethod( "thresholder" )

} # end of 'thresholder' function.

thresholder.default <- function( x, type = c( "binary", "replace.below" ), th, rule = ">=", replace.with = 0, ... ) {

    type <- match.arg( type )

    if( rule == ">=" ) id <- x >= th
    else if( rule == ">" ) id <- x > th
    else if( rule == "<=" ) id <- x <= th
    else if( rule == "<" ) id <- x < th

    xdim <- dim( x )

    if( is.null( xdim ) ) out <- numeric( length( x ) ) + replace.with
    else out <- matrix( replace.with, xdim[ 1 ], xdim[ 2 ] )

    if( type == "binary" ) out[ id ] <- 1 
    else out[ id ] <- x[ id ]

    return( out )

} # end of 'thresholder.default' function.

thresholder.SpatialVx <- function( x, type = c( "binary", "replace.below" ), th, rule = ">=", replace.with = 0, ...,
    time.point = 1, obs = 1, model = 1 ) {

    type <- match.arg( type )

    dat <- datagrabber( x, time.point = time.point, obs = obs, model = model )

    X <- dat$X
    Xhat <- dat$Xhat

    a <- attributes( x )

    u <- a$thresholds

    if( is.list( u ) ) {

	u1 <- u$X[ th, obs ]
	u2 <- u$Xhat[ th, model ]

    } else {

	u1 <- u2 <- u[ th ]

    }

    out1 <- thresholder( X, type = type, th = u1, rule = rule, replace.with = replace.with )
    out2 <- thresholder( Xhat, type = type, th = u2, rule = rule, replace.with = replace.with )

    return( list( X = out1, Xhat = out2 ) )

} # end of 'thresholder.SpatialVx' function.

hist.SpatialVx <- function( x, ..., time.point = 1, obs = 1, model = 1, threshold.num = NULL ) {

   tmp <- attributes( x )

    if( is.null( threshold.num ) ) {

	dat <- datagrabber( x, time.point = time.point, obs = obs, model = model )
	qs <- ""

    } else {

	dat <- thresholder( x, type = "replace.below", th = threshold.num, time.point = time.point, obs = obs, model = model )
	qs <- paste( ">= ", tmp$qs[ threshold.num ], sep = "" )

    }

   X <- dat$X
   Xhat <- dat$Xhat

   dn <- tmp$data.name
    vxname <- tmp$obs.name[ obs ]
    fcname <- tmp$model.name[ model ]
   nobs <- tmp$nobs
   nf <- tmp$nforecast

    if( !missing( obs ) ) if( length( obs ) > 1 ) stop("hist.SpatialVx: length of obs argument must be one." )

   if(!missing(model)) if(length(model) > 1) stop("hist.SpatialVx: length of model argument must be one.")

    h1 <- hist( X, main = dn, xlab = paste( vxname, "(", tmp$field.type, ", ", tmp$units, ")", sep = "" ), ... )
    h2 <- hist( Xhat, main = dn, xlab = paste( fcname, "(", tmp$field.type, ", ", tmp$units, ")", sep = "" ), ... )

   invisible( list( X = h1, Xhat = h2 ) )

} # end of 'hist.SpatialVx' function.

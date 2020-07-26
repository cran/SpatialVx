wavePurifyVx <- function( x, climate = NULL, which.stats = c("bias",
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL,
    rule = ">=", return.fields = FALSE, verbose = FALSE, ... ) {

    UseMethod( "wavePurifyVx" )

} # end of 'wavePurifyVx' function.

wavePurifyVx.SpatialVx <- function( x, climate = NULL, which.stats = c("bias",
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL,
    rule = ">=", return.fields = FALSE, verbose = FALSE, ...,
     time.point = 1, obs = 1, model = 1 ) {

    a <- attributes( x )

    dat <- datagrabber( x, time.point = time.point, obs = obs, model = model )

    u <- a$thresholds
    u$X <- u$X[, obs, drop = FALSE ]
    u$Xhat <- u$Xhat[, model, drop = FALSE ]

    thresholds <- u

    res <- wavePurifyVx( x = dat, climate = climate, which.stats = which.stats,
		thresholds = thresholds, return.fields = return.fields,
		verbose = verbose, ... )

    attr( res, "xdim" ) <- a$xdim
    attr( res, "time" ) <- a$time[ time.point ]
    attr( res, "loc" )  <- a$loc
    attr( res, "loc.byrow" ) <- a$loc.byrow
    attr( res, "data.name" ) <- a$data.name
    attr( res, "obs.name" ) <- a$obs.name[ obs ]
    attr( res, "model.name" ) <- a$model.name[ model ]
    attr( res, "field.type" ) <- a$field.type
    attr( res, "units" ) <- a$units
    attr( res, "projection" ) <- a$projection # maybe antiquated.
    attr( res, "reg.grid" ) <- a$reg.grid # perhaps also antiquated.
    attr( res, "map" ) <- a$map
    if( is.null( thresholds ) ) attr( res, "qs" ) <- a$qs
    attr( res, "msg" ) <- a$msg

    return( res )

} # end of 'wavePurifyVx.SpatialVx' function.

wavePurifyVx.default <- function( x, climate = NULL, which.stats = c("bias", 
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL, 
    rule = ">=", return.fields = FALSE, verbose = FALSE, ... ) {

    if( verbose ) begin.time <- Sys.time()

    if( verbose ) print( thresholds )

    out <- list()

    if( is.list( x ) && all( is.element( c( "X", "Xhat" ), names( x ) ) ) ) {

        X <- x$X
        Xhat <- x$Xhat

    } else if( is.list( x ) && length( x ) == 2 ) {

	X <- x[[ 1 ]]
	Xhat <- x[[ 2 ]]

    } else if( is.array( x ) && dim( x )[ 3 ] == 2 ) {

	X <- x[,, 1 ]
	Xhat <- x[,, 2 ]

    } else stop( "wavePurifyVx: invalid x argument." )

    attributes(out) <- list(...)

    xdim <- dim( X )

    if( !all( xdim == dim( Xhat ) ) ) stop( "wavePurifyVx: dimensions of verification and modelfield must be the same." )

    if( is.null( thresholds ) && any( is.element( c( "bias", "ts", "ets", "pod", "far", "f", "hk" ), which.stats ) ) ) {

        thresholds <- cbind( quantile( c( Xhat ), probs = c( 0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95 ) ),
		quantile( c( X ), probs = c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95) ) )

        attr( out, "qs" ) <- c( "0", "0.1", "0.25", "0.33", "0.5", "0.66", "0.75", "0.9", "0.95" )

    } else {

	if( !is.list( thresholds ) && is.vector( thresholds ) ) thresholds <- list( X = cbind( thresholds ), Xhat = cbind( thresholds ) )

	if( !is.list( thresholds ) ) stop( "wavePurifyVx: invalid thresholds argument.  Must be a list object." )
	if( !all( is.element( c( "X", "Xhat" ), names( thresholds ) ) ) )
	    stop( "wavePurifyVx: invalid thresholds argument.  Must have X and Xhat components." )

        thresholds <- cbind( thresholds$Xhat, thresholds$X )

        colnames( thresholds ) <- c( "Forecast", "Verification" )

        out$thresholds <- thresholds 

	# attr( out, "qs" ) <- 1:dim( thresholds )[ 1 ]

    } # end of if thresholds are null or not stmts.

    attr( out, "thresholds" ) <- thresholds
    
    args <- list(...)

    if( isnt.J <- is.null( args$J ) ) {

        J <- floor( log2( min( xdim ) ) )
        args$J <- J

    } else J <- args$J

    out$args <- args

    if( all( floor( log2( xdim ) ) == ceiling( log2( xdim ) ) ) ) dyadic <- TRUE
    else dyadic <- FALSE

    if (verbose) cat("\n", "Denoising the fields.\n")

    if( dyadic ) {

        if (isnt.J) {

            Z <- denoise.dwt.2d( x = X, J = J, ... )
            Y <- denoise.dwt.2d( x = Xhat, J = J, ... )

            if ( !is.null( climate ) ) Climate <- denoise.dwt.2d(x = climate, J = J, ...)

        } else {

            Z <- denoise.dwt.2d( x = X, ... )
            Y <- denoise.dwt.2d( x = Xhat, ... )

            if( !is.null( climate ) ) Climate <- denoise.dwt.2d( x = climate, ... )

        }

    } else {

        if( isnt.J ) {

            Z <- denoise.modwt.2d( x = X, J = J, ... )
            Y <- denoise.modwt.2d( x = Xhat, J = J, ... )

            if( !is.null( climate ) ) Climate <- denoise.modwt.2d( x = climate, J = J, ... )

        } else {

            Z <- denoise.modwt.2d( x = X, ... )
            Y <- denoise.modwt.2d( x = Xhat, ... )

            if( !is.null( climate ) ) Climate <- denoise.modwt.2d( x = climate, ... )

        }

    }

    if( return.fields ) {

        out$X <- X
        out$Xhat <- Xhat

        out$X2 <- Z
        out$Y2 <- Y

        if( !is.null( climate ) ) {

            out$Climate <- climate
            out$Climate2 <- Climate

        }

    }

    if( !is.null( thresholds ) ) q <- dim( thresholds )[ 1 ]
    else q <- 1

    if( is.element( "bias", which.stats ) ) out$bias <- numeric(q) + NA
    if( is.element( "ts", which.stats ) ) out$ts <- numeric(q) + NA
    if( is.element( "ets", which.stats ) ) out$ets <- numeric(q) + NA
    if( is.element( "pod", which.stats ) ) out$pod <- numeric(q) + NA
    if( is.element( "far", which.stats ) ) out$far <- numeric(q) + NA
    if( is.element( "f", which.stats ) ) out$f <- numeric(q) + NA
    if( is.element( "hk", which.stats ) ) out$hk <- numeric(q) + NA
    if( is.element( "mse", which.stats ) ) out$mse <- numeric(q) + NA
    if( !is.null( climate ) ) out$acc <- numeric(q) + NA

    if( verbose ) {

        cat("\n", "Looping through thresholds = \n")
        print( thresholds )
        cat("\n")

    }

    for( threshold in 1:q ) {

        if( verbose ) cat( threshold, " " )

        if( is.element( "mse", which.stats ) ) {

            if( !is.null( thresholds ) ) {

		X2 <- thresholder( x = Z, type = "replace.below", th = thresholds[ threshold, 2 ], rule = rule )
		Y2 <- thresholder( x = Y, type = "replace.below", th = thresholds[ threshold, 1 ], rule = rule )

                out$mse[ threshold ] <- vxstats( Y2, X2, which.stats = "mse" )$mse

            } else out$mse[ threshold ] <- vxstats( Y, Z, which.stats = "mse" )$mse

        }

        if( any( is.element( c( "bias", "ts", "ets", "pod", "far", "f", "hk" ), which.stats ) ) ) {

	    Xbin <- thresholder( x = Z, type = "binary", th = thresholds[ threshold, 2 ], rule = rule )
	    Ybin <- thresholder( x = Y, type = "binary", th = thresholds[ threshold, 1 ], rule = rule )

            if( threshold == 1 ) {

                dostats <- which.stats
                if( is.element( "mse", dostats ) ) dostats <- dostats[ dostats != "mse" ]

            }

            tmp <- vxstats( Ybin, Xbin, which.stats = dostats )

            if( is.element( "bias", which.stats ) ) out$bias[ threshold ] <- tmp$bias
            if( is.element( "ts", which.stats ) )   out$ts[ threshold ]   <- tmp$ts
            if( is.element( "ets", which.stats ) )  out$ets[ threshold ]  <- tmp$ets
            if( is.element( "pod", which.stats ) )  out$pod[ threshold ]  <- tmp$pod
            if( is.element( "far", which.stats ) )  out$far[ threshold ]  <- tmp$far
            if( is.element( "f", which.stats ) )    out$f[ threshold ]    <- tmp$f
            if( is.element( "hk", which.stats ) )   out$hk[ threshold ]   <- tmp$hk

        }

        if( !is.null( climate ) ) {

            if( !is.null( thresholds ) ) {

		X2 <- thresholder( x = Z, type = "replace.below", th = thresholds[ threshold, 2 ], rule = rule )
                Y2 <- thresholder( x = Y, type = "replace.below", th = thresholds[ threshold, 1 ], rule = rule )
		Clim <- thresholder( x = Climate, type = "replace.below", th = thresholds[ threshold, 1 ], rule = rule )

                denom <- sqrt( sum( colSums( ( Y2 - Clim )^2, na.rm = TRUE ), 
                  na.rm = TRUE ) ) * sqrt( sum( colSums( ( X2 - Clim )^2, 
                  na.rm = TRUE ), na.rm = TRUE ) )

                numer <- sum( diag( t( Y2 - Clim ) %*% ( X2 - Clim ) ), na.rm = TRUE )

            } else {

                denom <- sqrt( sum( colSums( ( Y - Climate )^2, na.rm = TRUE ), 
                  na.rm = TRUE ) ) * sqrt( sum( colSums( ( Z - Climate )^2, 
                  na.rm = TRUE ), na.rm = TRUE ) )

                numer <- sum( diag( t( Y - Climate ) %*% ( Z - Climate ) ), na.rm = TRUE )

            }

            out$acc[ threshold ] <- numer / denom

        }

    }

    class(out) <- "wavePurifyVx"

    if( verbose ) print( Sys.time() - begin.time )

    return( out )

} # end of 'wavePurifyVx.default' function.

summary.wavePurifyVx <- function( object, ... ) {

   object.attr <- attributes( object )

   if( !is.null( object.attr$thresholds ) ) {

        u <- Thresholds <- object.attr$thresholds
        q <- dim( u )[ 1 ]

        cat("\n", "Thresholds applied are:\n")
        print( Thresholds )

           if( is.null( object.attr$qs ) ) {

                if( all( u[, 1 ] == u[, 2 ] ) ) ulab <- as.character( u[, 1 ] )
                else ulab <- as.character( 1:q )

           } else ulab <- object.attr$qs

           if( !is.null( object$bias ) ) {

                bias <- matrix( object$bias, nrow = 1 )
                colnames( bias ) <- ulab
                cat("\n", "Frequency Bias: \n")
                print( bias )

           }

           if(!is.null(object$ts)) {

                ts <- matrix( object$ts, nrow = 1 )
                colnames( ts ) <- ulab
                cat("\n", "Threat Score: \n")
                print( ts )

           }

           if( !is.null( object$ets ) )  {

                ets <- matrix( object$ets, nrow = 1 )
                colnames( ets ) <- ulab
                cat("\n", "Gilbert Skill Score: \n")
                print( ets )

           }

           if( !is.null( object$pod ) )  {

                pod <- matrix( object$pod, nrow = 1 )
                colnames( pod ) <- ulab
                cat("\n", "Probability of Detecting an Event (POD): \n")
                print( pod )

           }

           if( !is.null( object$far ) ) {

                far <- matrix( object$far, nrow = 1 )
                colnames( far ) <- ulab
                cat("\n", "False Alarm Ratio: \n")
                print( far )

           }

           if( !is.null( object$f ) ) {

                f <- matrix( object$f, nrow = 1 )
                colnames( f ) <- ulab
                cat("\n", "False Alarm Rate: \n")
                print( f )

           }

           if( !is.null( object$hk ) ) {

                hk <- matrix( object$hk, nrow = 1 )
                colnames( hk ) <- ulab
                cat("\n", "Hanssen-Kuipers Score: \n")
                print( hk )

           }

           if( !is.null( object$mse ) ) {

                mse <- matrix( object$mse, nrow = 1 )
                colnames( mse ) <- ulab
                cat("\n", "MSE: \n")
                print( mse )

           }

           if( !is.null( object$acc ) ) {

                acc <- matrix( object$acc, nrow = 1 )
                colnames( acc ) <- ulab
                cat("\n", "ACC: \n")
                print( acc )

           }

   } else {

        if( !is.null( object$mse ) ) {

                mse <- matrix( object$mse, nrow = 1 )
                cat("\n", "MSE: \n")
                print( mse )

           }

           if( !is.null( object$acc ) ) {

                acc <- matrix( object$acc, nrow = 1 )
                cat("\n", "ACC: \n")
                print( acc )

           }

   } # end of whether or not thresholds applied stmts.

   invisible()

} # end of 'summary.wavePurifyVx' function.

plot.wavePurifyVx <- function(x, ..., col = c( "gray", tim.colors(64) ), zlim, mfrow, horizontal = TRUE, type = c("stats", "fields") ) {

    a <- attributes( x )

    type <- tolower( type )
    type <- match.arg( type )

    if( type == "stats" ) {

	if( !missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxStats( x = x, ..., col = col, zlim = zlim, mfrow = mfrow, horizontal = horizontal, type = type  )
	else if( missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxStats( x = x, ..., col = col, mfrow = mfrow, horizontal = horizontal, type = type  )
	else if( !missing( zlim ) && missing( mfrow ) ) plot.wavePurifyVxStats( x = x, ..., col = col, zlim = zlim, horizontal = horizontal, type = type  )
	else  plot.wavePurifyVxStats( x = x, ..., col = col, horizontal = horizontal, type = type  )

	# class( x ) <- "wavePurifyVxStats"

    } else if( type == "fields" && a$map ) {

	if( !missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxFieldsMap( x = x, ..., col = col, zlim = zlim, mfrow = mfrow, horizontal = horizontal, type = type  )
        else if( missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxFieldsMap( x = x, ..., col = col, mfrow = mfrow, horizontal = horizontal, type = type  )
        else if( !missing( zlim ) && missing( mfrow ) ) plot.wavePurifyVxFieldsMap( x = x, ..., col = col, zlim = zlim, horizontal = horizontal, type = type  )
        else  plot.wavePurifyVxFieldsMap( x = x, ..., col = col, horizontal = horizontal, type = type  )

	# class( x ) <- "wavePurifyVxFieldsMap"

    } else if( type == "fields" && !a$map ) {

	if( !missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxFieldsNoMap( x = x, ..., col = col, zlim = zlim, mfrow = mfrow, horizontal = horizontal, type = type  )
        else if( missing( zlim ) && !missing( mfrow ) ) plot.wavePurifyVxFieldsNoMap( x = x, ..., col = col, mfrow = mfrow, horizontal = horizontal, type = type  )
        else if( !missing( zlim ) && missing( mfrow ) ) plot.wavePurifyVxFieldsNoMap( x = x, ..., col = col, zlim = zlim, horizontal = horizontal, type = type  )
        else  plot.wavePurifyVxFieldsNoMap( x = x, ..., col = col, horizontal = horizontal, type = type  )

	# class( x ) <- "wavePurifyVxFieldsNoMap"

    }

    # UseMethod( "plot", x )

} # end of 'plot.wavePurifyVx' function.

plot.wavePurifyVxFieldsMap <- function( x, ..., col = c( "gray", tim.colors(64) ), zlim, mfrow, horizontal = TRUE ) {

    op <- par()

    a <- attributes( x )

    loc <- a$loc

    xdim <- a$xdim

    l <- list( x = matrix( loc[, 1 ], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ),
		y = matrix( loc[, 2 ], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ) ) 

    lr <- apply( loc, 2, range, finite = TRUE )

    ax <- list( x = pretty( round( loc[, 1], digits = 2 ) ), 
                y = pretty( round( loc[, 2], digits = 2 ) ) )

    X <- x$X
    Xhat <- x$Xhat

    # Have the fields been de-noised and returned with the object?
    denX <- !is.null( x$X2 )

    if( denX ) {

	X2 <- x$X2
	Y2 <- x$Y2

	if( !is.null( x$Climate ) ) {

            Clim <- x$Climate
            Clim2 <- x$Climate2

            if( missing( zlim ) ) zlim <- range( c( c( X ), c( Xhat ), c( X2 ), c( Y2 ), c( Clim ), c( Clim2 ) ), finite = TRUE )

	    if( missing( mfrow ) ) mfrow <- c( 2, 3 )

        } else {

	    if( missing( zlim ) ) zlim <- range( c( c( X ), c( Xhat ), c( X2 ), c( Y2 ) ), finite = TRUE )

	    if( missing( mfrow ) ) mfrow <- c( 2, 2 )

	} # end of if climatology is used or not stmts.

	par( mfrow = mfrow, oma = c( 0, 0, 2, 0 ) )


    } else if( !denX ) stop("plot.wavePurifyVx: return.fields must be TRUE in call to wavePurifyVx to use type fields.")
    # end of if else 'denX' stmts.

    # Image of Verification field (X)
    map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

    # axis( 1, at = ax$x, labels = ax$x )
    # axis( 2, at = ax$y, labels = ax$y )

    if( !is.null( a$obs.name ) ) title( a$obs.name )
    else title( "X" )

    poly.image( l$x, l$y, X, col = col, zlim = zlim, add = TRUE, ... )

    map( add = TRUE, lwd = 1.5 )
    map( add = TRUE, database = "state" )
    map.axes()

    # Image of Model field (Xhat)
    map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

    # axis( 1, at = ax$x, labels = ax$x )
    # axis( 2, at = ax$y, labels = ax$y )

    if( !is.null( a$model.name ) ) title( a$model.name )
    else title( "Xhat" )

    poly.image( l$x, l$y, Xhat, col = col, zlim = zlim, add = TRUE, ... )

    map( add = TRUE, lwd = 1.5 )
    map( add = TRUE, database = "state" )
    map.axes()

    # Cllimate field (if present)
    if (!is.null(x$Climate)) {

	 map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

        # axis( 1, at = ax$x, labels = ax$x )
        # axis( 2, at = ax$y, labels = ax$y )

        title( "Climatology" )

        poly.image( l$x, l$y, Clim, col = col, zlim = zlim, add = TRUE, ... )

        map( add = TRUE, lwd = 1.5 )
        map( add = TRUE, database = "state" )
	map.axes()

    } # end of if climatology present stmts.

    # Image of de-noised Verification field (X2)
    map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

    # axis( 1, at = ax$x, labels = ax$x )
    # axis( 2, at = ax$y, labels = ax$y )

    if( !is.null( a$obs.name ) ) title( paste( a$obs.name, " (de-noised)", sep = "" ) )
    else title( "X (de-noised)" )

    poly.image( l$x, l$y, X2, col = col, zlim = zlim, add = TRUE, ... )

    map( add = TRUE, lwd = 1.5 )
    map( add = TRUE, database = "state" )
    map.axes()

    # Image of de-noised Model field (Y2)
    map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

    # axis( 1, at = ax$x, labels = ax$x )
    # axis( 2, at = ax$y, labels = ax$y )

    if( !is.null( a$model.name ) ) title( paste( a$model.name, " (de-noised)", sep = "" ) )
    else title( "Xhat (de-noised)" )

    poly.image( l$x, l$y, Y2, col = col, zlim = zlim, add = TRUE, ... )

    map( add = TRUE, lwd = 1.5 )
    map( add = TRUE, database = "state" )
    map.axes()

     # Cllimate field (if present)
    if (!is.null(x$Climate)) {

         map( xlim = lr[, 1 ], ylim = lr[, 2 ], type = "n" )

        # axis( 1, at = ax$x, labels = ax$x )
        # axis( 2, at = ax$y, labels = ax$y )

        title( "Climatology (de-noised)" )

        poly.image( l$x, l$y, Clim2, col = col, zlim = zlim, add = TRUE, ... )

        map( add = TRUE, lwd = 1.5 )
        map( add = TRUE, database = "state" )
	map.axes()

    } # end of if climatology present stmts.

    image.plot( X, legend.only = TRUE, horizontal = horizontal, col = col, zlim = zlim, ... )

    mtext( a$msg, line = 0.05, outer = TRUE )

    par( mfrow = op$mfrow, oma = op$oma )

    invisible()

} # end of 'plot.wavePurifyVxFieldsMap' function.

plot.wavePurifyVxFieldsNoMap <- function( x, ..., col = c( "gray", tim.colors(64) ), zlim, mfrow, horizontal = TRUE ) {

    op <- par()

    a <- attributes( x )

    X <- x$X
    Xhat <- x$Xhat

    denX <- !is.null( x$X2 )

    # Have the fields been de-noised and returned with the object?
    if( denX ) {

        X2 <- x$X2
        Y2 <- x$Y2

        if( !is.null( x$Climate ) ) {

            Clim <- x$Climate
            Clim2 <- x$Climate2

            if( missing( zlim ) ) zlim <- range( c( c( X ), c( Xhat ), c( X2 ), c( Y2 ), c( Clim ), c( Clim2 ) ), finite = TRUE )

            if( missing( mfrow ) ) mfrow <- c( 2, 3 )

        } else {

            if( missing( zlim ) ) zlim <- range( c( c( X ), c( Xhat ), c( X2 ), c( Y2 ) ), finite = TRUE )

            if( missing( mfrow ) ) mfrow <- c( 2, 2 )

        } # end of if climatology is used or not stmts.

        par( mfrow = mfrow, oma = c( 0, 0, 2, 0 ) )


    } else stop("plot.wavePurifyVx: return.fields must be TRUE in call to wavePurifyVx to use type fields.")

    # Image of Verification field (X)
    image( X, col = col, zlim = zlim  )
    if( !is.null( a$obs.name ) ) title( a$obs.name )
    else title( "X" )

    # Image of MOdel field (Xhat)
    image( Xhat, col = col, zlim = zlim )
    if( !is.null( a$model.name ) ) title( a$model.name )
    else title( "Xhat" )

    # Image of climatology (if present)
    if( !is.null( x$Climate ) ) {

        image( Clim, col = col, zlim = zlim )
        title( "Climatology" )

    } # end of if climatology present stmt.

    # Image of de-noised Verification field (X2)
    image( X2, col = col, zlim = zlim )
    if( !is.null( a$obs.name ) ) title( paste( a$obs.name, " (de-noised)", sep = "" ) )
    else title( "X (de-noised)" )

    # Image of de-noised Model field (Y2)
    image( Y2, col = col, zlim = zlim )
    if( !is.null( a$model.name ) ) title( paste( a$model.name, " (de-noised)", sep = "" ) )
    else title( "Xhat (de-noised)" )

     # Image of de-noised climatology (if present)
    if( !is.null( x$Climate ) ) {

        image( Clim2, col = col, zlim = zlim )
        title( "Climatology (de-noised)" )

    } # end of if climatology present stmt.

    image.plot( X, legend.only = TRUE, horizontal = horizontal, col = col, zlim = zlim )

    mtext( a$msg, line = 0.05, outer = TRUE )

    par( mfrow = op$mfrow, oma = op$oma )

    invisible()

} # end of 'plot.wavePurifyVxFieldsNoMap' function.

plot.wavePurifyVxStats <- function( x, ..., mfrow ) {

    op <- par()

    a <- attributes( x )
    nm <- names( x )

    if( missing( mfrow ) ) {

        if( any( is.element( c( "bias", "ts", "ets", "pod", "far", "hk" ), nm ) ) ) {

            if( all( is.element( c( "mse", "acc" ), nm ) ) ) mfrow <- c( 2, 2 )
	    else if( any( is.element( c( "mse", "acc" ), nm ) ) ) mfrow <- c( 1, 2 )
            else mfrow <- c( 1, 1 )

        } else {

	    if( all( is.element( c( "mse", "acc" ), nm ) ) ) mfrow <- c( 1, 2 )
	    else mfrow <- c( 1, 1 )

        } # if else any of the first set of types stmts.

    }

    par( mfrow = mfrow, oma = c( 0, 0, 2, 0 ) )

     if( !is.null( a$thresholds ) ) {

        if( is.list( a$thresholds ) && "X" %in% names( a$thresholds ) ) q <- dim( a$thresholds$X )[ 1 ]
	else if( is.matrix( a$thresholds ) ) q <- dim( a$thresholds )[ 1 ]
	else stop( "plot.wavePurifyVxStats: sorry, but there is a problem with the thresholds.  My fault!" )

        if (q >= 2) {

            if( is.null( a$units ) ) unt <- NULL
            else unt <- paste("(", a$units, ")", sep = "")

            if( any( is.element( c( "bias", "ts", "ets", "pod", "far", "hk" ), nm ) ) ) {

                if( is.element( "hk", nm ) ) yl <- c( -1, 1 )
                else if( is.element( "ets", nm ) ) yl <- c( -1/3, 1 )
                else yl <- c( 0, 1 )

                plot( 1:q, runif( q, -1, 1 ), type = "n", xaxt = "n", 
                  ylim = yl, xlab = paste( "Threshold ", unt, sep = "" ), ylab = "" )

                legnames <- c( "TS", "GSS", "POD", "FAR", "F", "HK",
				"Bias" )[ is.element( c( "ts", "ets",
				"pod", "far", "f", "hk", "bias" ), nm ) ]

                nl <- length( legnames )

                if( nl == 7 ) cl <- c( 1:6, 8 )
                else cl <- 1:nl

                legend( "topright", legend = legnames, lty = cl, col = cl, bty = "n")

                i <- 1

                if( is.element( "ts", nm ) ) {

                  lines( 1:q, x$ts, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if( is.element( "ets", nm ) ) {

                  lines( 1:q, x$ets, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if( is.element( "pod", nm ) ) {

                  lines( 1:q, x$pod, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if (is.element("far", nm)) {

                  lines( 1:q, x$far, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if (is.element("f", nm)) {

                  lines( 1:q, x$f, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if (is.element("hk", nm)) {

                  lines( 1:q, x$hk, col = i, lty = i, lwd = 1.5 )
                  i <- i + 1

                }

                if (is.element("bias", nm)) {

                  if (i == 7) i <- 8

                  brg <- range( x$bias, finite = TRUE )

                  if( length( unique( brg ) ) == 1 ) brg <- c( brg[ 1 ] - 0.01, brg[ 2 ] + 0.01 )

                  par(usr = c(1, q, brg))
                  lines(1:q, x$bias, col = i, lty = i, lwd = 1.5)
                  bax <- pretty(sort(unique(x$bias)))
                  axis(4, at = bax, labels = bax, col = i)

                }

            } # end of if any non-MSE statistics included stmts.

            if (is.element("mse", nm)) 
                plot(1:q, x$mse, type = "l", col = "darkblue", 
                  xlab = paste("Threshold ", unt, sep = ""), 
                  ylab = "MSE")
            if (is.element("acc", nm)) 
                plot(1:q, x$acc, type = "l", col = "darkblue", 
                  xlab = paste("Threshold ", unt, sep = ""), 
                  ylab = "ACC")
        }
    }

    mtext( a$msg, line = 0.05, outer = TRUE )

    par( mfrow = op$mfrow, oma = op$oma )

    invisible()

} # end of 'plot.wavePurifyVxStats' function.

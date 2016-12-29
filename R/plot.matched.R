plot.matched <- function( x, mfrow = c(1, 2), ... ) {

    a <- attributes( x )
    if( a$map ) class( x ) <- "matchedMap"
    else class( x ) <- "matchedNoMap"

    UseMethod( "plot", x )

} # end of 'plot.matched' function.

plot.matchedMap <- function( x, mfrow = c(1, 2), ... ) {

    a <- attributes( x )

    if( !is.null( mfrow ) ) {

	op <- par()
	par( mfrow = mfrow, oma = c(0, 0, 2, 0) )
	msg <- paste( a$data.name, "\n(", a$field.type, ", ", a$units, ", time = ", a$time.point, ")", sep = "" )

    }

    if( !( "deltamm" %in% x$match.type[ 1 ] ) && ( is.null( x$MergeForced ) || !x$MergeForced ) )
	warning("plot.matched: If MergeForce is not run, then colors will not identify matches.")

    matches <- x$matches

    xdim <- a$xdim

    Aun <- x$unmatched$X
    Bun <- x$unmatched$Xhat
    nX <- length( Aun )
    nY <- length( Bun )

    X <- x$X.labeled
    Xhat <- x$Y.labeled

    if ( any( dim( matches ) == 0 ) ) {

        nm <- 0
        X[ X > 0 ] <- 1

        Xhat <- x$Y.labeled
        Xhat[ Xhat > 0 ] <- 1
        icolX <- icolXhat <- icol <- c("white", "gray")

    } else {

        nm <- dim( matches )[ 1 ]
	icolX <- icolXhat <- icol <- c( "white", rainbow( nm ) )

	if( nX > 0 ) {

	    X[ X > nm ] <- nm + 1
	    icolX <- c( icolX, "gray" )

	}

	if( nY > 0 ) {

	    Xhat[ Xhat > nm ] <- nm + 1
	    icolXhat <- c( icolXhat, "gray" )

	}

	if( nX > 0 || nY > 0 ) icol <- c( icol, "gray" )

    }

    if( nX > nY ) Z <- X
    else Z <- Xhat

    loc <- a$loc

    l <- list( x = matrix( loc[, 1], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ),
		y = matrix( loc[, 2], xdim[ 1 ], xdim[ 2 ], byrow = a$loc.byrow ) )

    lr <- apply( loc, 2, range, finite = TRUE )

    map( xlim = lr[, 1], ylim = lr[, 2], type = "n" )
    title( a$obs.name )
    poly.image( l$x, l$y, X, col = icolX, add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )

    map( xlim = lr[, 1], ylim = lr[, 2], type = "n" )
    title( a$model.name )
    poly.image( l$x, l$y, Xhat, col = icolX, add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )

    image.plot( Z, col = icol, legend.only = TRUE, ... )

    if( !is.null( mfrow ) ) {

	title("")
	mtext( msg, line = 0.5, outer = TRUE )
	par( mfrow = op$mfrow, oma = op$oma )

    }

    invisible()

} # end of 'plot.matchedMap' function.

plot.matchedNoMap <- function( x, mfrow = c(1, 2), ... ) {

    a <- attributes( x )

    if( !is.null( mfrow ) ) {

        op <- par()
        par( mfrow = mfrow, oma = c(0, 0, 2, 0) )
        msg <- paste( a$data.name, "\n(", a$field.type, ", ", a$units, ", time = ", a$time.point, ")", sep = "" )

    }

    if( !( "deltamm" %in% x$match.type[ 1 ] ) && ( is.null( x$MergeForced ) || !x$MergeForced ) )
        warning("plot.matched: If MergeForce is not run, then colors will not identify matches.")

    matches <- x$matches

    xdim <- dim( x$X.feats )

    Aun <- x$unmatched$X
    Bun <- x$unmatched$Xhat
    nX <- length( Aun )
    nY <- length( Bun )

    X <- x$X.labeled
    Xhat <- x$Y.labeled

    if ( any( dim( matches ) == 0 ) ) {

        nm <- 0
        X[ X > 0 ] <- 1

        Xhat <- x$Y.labeled
        Xhat[ Xhat > 0 ] <- 1
        icolX <- icolXhat <- icol <- c("white", "gray")

    } else {

        nm <- dim( matches )[ 1 ]
        icolX <- icolXhat <- icol <- c( "white", rainbow( nm ) )

        if( nX > 0 ) {

            X[ X > nm ] <- nm + 1
            icolX <- c( icolX, "gray" )

        }

        if( nY > 0 ) {

            Xhat[ Xhat > nm ] <- nm + 1
            icolXhat <- c( icolXhat, "gray" )

        }

        if( nX > 0 || nY > 0 ) icol <- c( icol, "gray" )

    }

    if( nX > nY ) Z <- X
    else Z <- Xhat

    image( X, col = icolX )
    title( a$obs.name )
    image( Xhat, col = icolXhat )
    title( a$model.name )
    image.plot( Z, col = icol, legend.only = TRUE, ... )

    if( !is.null( mfrow ) ) {

	title("")
        mtext( msg, line = 0.5, outer = TRUE )
        par( mfrow = op$mfrow, oma = op$oma )

    }

    invisible()

} # end of 'plot.matchedNoMap' function.

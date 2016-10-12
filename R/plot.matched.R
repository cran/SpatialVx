plot.matched <- function( x, ... ) {

    a <- attributes( x )

    if( is.null( x$MergeForced ) || !x$MergeForced ) warning("plot.matched: If MergeForce is not run, then colors will not identify matches.")

    matches <- x$matches

    xdim <- dim( x$X.feats )

    Aun <- x$unmatched$X
    Bun <- x$unmatched$Xhat
    nX <- length( Aun )
    nY <- length( Bun )

    X <- x$X.labeled
    Xhat <- x$Y.labeled

    if (any(dim(matches) == 0)) {

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

    isLocList <- is.list( a$loc )

    if( isLocList && !( c("x", "y") %in% names( a$loc ) ) ) {

	a$loc <- cbind( rep( 1:xdim[ 1 ], xdim[ 2 ] ),
		    rep( 1:xdim[ 2 ], each = xdim[ 1 ] ) )

	isLocList <- FALSE

    } # end of if loc is not what is expected or hoped for.

    if( isLocList ) {

        if( a$map ) {

	    lr <- lapply( a$loc, range, finite = TRUE )
	    map( xlim = lr$x, ylim = lr$y, type = "n" )
	    poly.image( a$loc$x, a$loc$y, X, col = icolX, add = TRUE )
	    map( add = TRUE )
            map( database = "state", add = TRUE )

	    lr <- lapply( a$loc, range, finite = TRUE )
            map( xlim = lr$x, ylim = lr$y, type = "n" )
            poly.image( a$loc$x, a$loc$y, Xhat, col = icolX, add = TRUE )
            map( add = TRUE )
            map( database = "state", add = TRUE )

	} else {

	    poly.image( a$loc$x, a$loc$y, X, col = icolX ) 
            poly.image( a$loc$x, a$loc$y, Xhat, col = icolX )

	}

	image.plot( Z, col = icol, legend.only = TRUE, ... )

    } else if( is.matrix( a$loc ) && dim( a$loc )[ 2 ] == 2 ) {

	if( a$map ) {

	    lr <- apply( a$loc, 2, range, finite = TRUE )
	    map( xlim = lr[,1], ylim = lr[,2], type = "n" )
	    image( X, col = X, add = TRUE )
	    map( add = TRUE )
            map( database = "state", add = TRUE )

	    map( xlim = lr[,1], ylim = lr[,2], type = "n" )
	    image( Xhat, col = icolXhat, add = TRUE )
	    map( add = TRUE )
	    map( database = "state", add = TRUE )

	} else {

	    image( X, col = icolX )
	    image( Xhat, col = icolXhat )
	}

        image.plot( Z, col = icol, legend.only = TRUE, ... )

    } # end of if else use projection mapping or not stmts.

    invisible()

} # end of 'plot.matched' function.

iwarper <- function( x0, x1, nc = 4, labcol = "magenta", col = c( "gray", tim.colors(64) ), zlim, cex = 2, alwd = 1.25, ... ) {

    par( mfrow = c(2,2) )

    xdim <- dim( x0 )

    if( missing( x0 ) ) stop( "iwarper: must supply an x0 matrix." )
    if( missing( x1 ) ) stop( "iwarper: must supply a x1 matrix." )

    if( !all( xdim == dim( x1 ) ) ) stop( "iwarper: x0 and x1 must have same dimension." )

    loc <- cbind( rep( 1:xdim[ 1 ], xdim[ 2 ] ), rep( 1:xdim[ 2 ], each = xdim[ 1 ] ) )

    if( missing( zlim ) ) zl <- range( c( c(x0), c(x1) ), finite = TRUE )
    else zl <- zlim

    l <- list( x = matrix( loc[,1], xdim[ 1 ], xdim[ 2 ] ),
		y = matrix( loc[,2], xdim[ 1 ], xdim[ 2 ] ) )


    # image.plot( l$x, l$y, x0, main = "0-energy image", col = col, zlim = zl, ... )
    poly.image( l$x, l$y, x0, main =  "0-energy image", col = col, zlim = zl, ... )
    contour( x = seq(min(loc[,1]), max(loc[,1]),, nrow(x1) ), y = seq( min(loc[,2]), max(loc[,2]),, ncol(x1) ),  x1, add = TRUE, col = "lightgray", lty = 2 )

    cat( "\nChoose ", nc, " 0-energy control locations.\n")

    p0 <- p1 <- numeric( 0 )

    for( i in 1:nc ) {

	hold0 <- locator( 1 )
	p0 <- rbind( p0, unlist( hold0 ) )

	# points( p0[i,1], p0[i,2], pch = as.character( i ), cex = cex, col = labcol )
	text( p0[i,1], p0[i,2], labels = as.character( i ), cex = cex, col = labcol )

    } # end of for 'i' loop.

    # image.plot( l$x, l$y, x1, main = "1-energy image", col = col, zlim = zl, ... )
    poly.image( l$x, l$y, x1, main = "1-energy image", col = col, zlim = zl, ... )
    contour( x = seq(min(loc[,1]), max(loc[,1]),, nrow(x0) ), y = seq( min(loc[,2]), max(loc[,2]),, ncol(x0) ), x0, add = TRUE, col = "lightgray", lty = 2 )

    cat( "\nChoose ", nc, " 1-energy control locations.\n")

    for( i in 1:nc ) {

        hold1 <- locator( 1 )
        p1 <- rbind( p1, unlist( hold1 ) )

        text( p1[i,1], p1[i,2], labels = as.character( i ), cex = cex, col = labcol )

    } # end of for 'i' loop.

    image.plot( x1, legend.only = TRUE, col = col, zlim = zl, ... )

    # image.plot( l$x, l$y, x1 - x0, main = "Original Error Field\nx1 - x0", col = tim.colors(64), ... )

    theta <- warpTpsMatrices( p0 = p0, s = loc )
    wpts <- warpTps( p1 = p1, B = theta$B )

    x1.def <- Fint2d( x1, Ws = wpts, s = loc )

    plot( c(0,1), type = "n", xaxt = "n", yaxt= "n", xlab = "", ylab = "" )
    # image.plot( l$x, l$y, x1.def, main = "Deformed 1-energy image", col = col, zlim = zl, ... )
    poly.image( l$x, l$y, x1.def, main = "Deformed 1-energy image", col = col, zlim = zl, ... )
    image.plot( x1.def, legend.only = TRUE, col = col, zlim = zl, ... )

    arrows(p1[,1], p1[,2], p0[,1], p0[,2], col="magenta", length=0.2, lwd = alwd )

    # image.plot( l$x, l$y, x1.def - x0, main = "Resulting Error Field\nx1.def - x0", col = tim.colors(64), ... )

    args <- list( ... )

    out <- list()

    out$Im0 <- x0
    out$Im1 <- x1
    out$Im1.def <- x1.def

    out$p0 <- p0
    out$p1 <- p1

    out$warped.locations <- wpts
    out$s <- loc

    if( !is.null( args$imethod ) ) out$imethod <- args$imethod

    out$theta <- theta

    par( mfrow = c(1,1) )

    class( out ) <- "iwarped"

    return( out )

} # end of 'iwarper' function.

summary.iwarped <- function( object, ... ) {

    out <- list( ... )

    rmse0 <- sqrt( mean( (object$Im1 - object$Im0)^2, na.rm = TRUE ) )
    rmse1 <- sqrt( mean( (object$Im1.def - object$Im0)^2, na.rm = TRUE ) )
    redux <- ( (rmse0 - rmse1) / rmse0 ) * 100

    k <- nrow( object$p0 )

    Be <- object$theta$iL[ 1:k, 1:k ]
    p1 <- object$p1

    minBe <- sum( diag( t( p1 ) %*% Be %*% p1 ) )

    res <- c( rmse0, rmse1, redux, minBe )
    names( res ) <- c( "RMSE0", "RMSE1", "% RMSE Reduction", "Bending Energy" )

    silent <- out$silent
    if( is.null( out$silent ) ) silent <- FALSE

    if( !silent ) {

        a <- attributes( object )
        print( paste( a$data.name, collapse = " vs ") )

        print( res )

    }

    out$result <- res

    invisible( out )

} # end of 'summary.iwarped' function.

plot.iwarped <- function( x, col = c( "gray", tim.colors(64) ), alwd = 1.5, ... ) {

    par( mfrow = c(2,3) )

    image.plot( x$Im0, main = "0-energy field", col = col, ... )
    image.plot( x$Im1, main = "1-energy field", col = col, ... )
    image.plot( x$Im1 - x$Im0, main = "Error Field", col = tim.colors(64), ... )

    xdim <- dim( x$Im0 ) 
    h <- matrix( sqrt( ( x$warped.locations[,1] - x$s[,1] )^2 + ( x$warped.locations[,2] - x$s[,2] )^2 ), xdim[ 1 ], xdim[ 2 ] )
    image.plot( matrix( x$s[,1], xdim[ 1 ], xdim[ 2 ] ),
		matrix( x$s[,2], xdim[ 1 ], xdim[ 2 ] ),
		h, main = "Distance travelled", col = c( "gray", tim.colors(64) ), ... )
    arrows(x$p1[,1], x$p1[,2], x$p0[,1], x$p0[,2], col="magenta", length=0.2, lwd = alwd )

    image.plot( x$Im1.def, main = "Deformed 1-energy field", col = col, ... )
    image.plot( x$Im1.def - x$Im0, main = "Error Field\n(after warping)", col = tim.colors(64), ... )

    par( mfrow = c(1,1) )

    invisible()

} # end of 'plot.iwarped' function.

# see warper.R for the original function, which needs some work.
# this function relies on 'equal.axis'.

plot.warped.lonlat <- function( x, col = c( "gray", tim.colors(64) ), alwd = 1.5, ..., lonlat = NULL ) {

    if( is.null( lonlat ) ) stop( "plot.warped.lonlat: using the wrong function." )

    if( is.list( lonlat ) ) xy = equal.axis( x = c( lonlat$x ), y = c( lonlat$y ) )
    else if( is.matrix( lonlat ) && dim( lonlat )[ 2 ] == 2 ) {

	xdim <- dim( x$Im0 )
	xy = equal.axis( x = lonlat[, 1], y = lonlat[, 2] )
	lonlat <- list( x = matrix( lonlat[, 1], xdim[ 1 ], xdim[ 2 ] ),
			y = matrix( lonlat[, 2], xdim[ 1 ], xdim[ 2 ] ) ) 

    } else stop( "plot.warped.lonlat: invalid lonlat argument." ) 

    zl1 = range( c( c( x$Im0 ), c( x$Im1 ), c( x$Im1.def ) ), finite = TRUE )
    zl2 = range( c( c( x$Im1 - x$Im0 ), c( x$Im1.def - x$Im0 ) ), finite = TRUE )

    par( mfrow = c(2,3), mar = c(6.1, 2.1, 6.1, 0.1) )

    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "0-energy field", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, x$Im0, col = col, ..., zlim = zl1, add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    image.plot( x$Im0, legend.only = TRUE, horizontal = TRUE, col = col, zlim = zl1 )

    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "1-energy field", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, x$Im1, col = col, ..., zlim = zl1, add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    image.plot( x$Im1, legend.only = TRUE, horizontal = TRUE, col = col, zlim = zl1 )

    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "Error Field", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, x$Im1 - x$Im0, col = tim.colors(64), zlim = zl2, ..., add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    image.plot( x$Im1 - x$Im0, legend.only = TRUE, horizontal = TRUE, col = col, zlim = zl2 )

    xdim <- dim( x$Im0 )
    h <- matrix( sqrt( ( x$warped.locations[, 1] - x$s[, 1] )^2 + ( x$warped.locations[, 2] - x$s[, 2] )^2 ),
        xdim[ 1 ], xdim[ 2 ] )
    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "Distance Travelled", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, h, col = c( "gray", tim.colors(64) ), ..., add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )

    fitX = lm( xPrime ~ x, data = data.frame( x = x$s[,1], xPrime = c(lonlat$x) ) )
    x1   = predict( fitX, newdata = list( x = x$p0[,1] ) )
    x2   = predict( fitX, newdata = list( x = x$p1[,1] ) )

    fitY = lm( yPrime ~ y, data = data.frame( y = x$s[,2], yPrime = c(lonlat$y) ) )
    y1   = predict( fitY, newdata = list( y = x$p0[,2] ) )
    y2   = predict( fitY, newdata = list( y = x$p1[,2] ) )

    arrows(x1, x2, y1, y2, col="magenta", length=0.2, lwd = alwd )
    image.plot( h, legend.only = TRUE, horizontal = TRUE, col = c( "gray", tim.colors(64) ) )

    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "Deformed 1-energy field", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, x$Im1.def, col = col, ..., zlim = zl1, add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    image.plot( x$Im1.def, legend.only = TRUE, horizontal = TRUE, col = col, zlim = zl1 )

    plot( xy$x, xy$y, type = "n", xlab = "", ylab = "", main = "Error Field\n(after warping)", cex.main = 1.25 )
    poly.image( lonlat$x, lonlat$y, x$Im1.def - x$Im0, col = tim.colors(64), zlim = zl2, ..., add = TRUE )
    map( add = TRUE )
    map( database = "state", add = TRUE )
    image.plot( lonlat$x, lonlat$y, x$Im1.def - x$Im0, legend.only = TRUE, horizontal = TRUE,
	col = tim.colors(64), zlim = zl2 )

    par( mfrow = c(1,1) )

    invisible()

} # end of 'plot.warped.lonlat' function.

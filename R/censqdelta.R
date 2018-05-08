censqdelta <- function( x, y, N, const = Inf, p = 2, ... ) {

    # op <- par()

    if( !all( is.element( unique( x ), c( 0, 1 ) ) ) ) x <- thresholder( x, th = 0, rule = ">" )
    if( !all( is.element( unique( y ), c( 0, 1 ) ) ) ) y <- thresholder( y, th = 0, rule = ">" )

    # par( mfrow = c( 2, 2 ) )

    # image( x, col = c( "white", "darkblue" ) )
    # contour( y, add = TRUE )
    # image( y, col = c( "white", "darkblue" ) )
    # contour( x, add = TRUE )

    # dx <- as.im( x )
    # dy <- as.im( y )

    # dx <- solutionset( dx > 0 )
    # dy <- solutionset( dy > 0 )

    # dx <- distmap( dx )
    # if( is.finite( const ) ) dx$v[ dx$v > const ] <- const
    # dy <- distmap( dy )
    # if( is.finite( const ) ) dy$v[ dy$v > const ] <- const

    # image( x + y, col = c( "white", "lightblue", "darkblue" ) )
    # title( "Binary A and B superimposed" )
    # image.plot( abs( as.matrix( dx ) - as.matrix( dy ) ), main = "Abs. Diff. in Distance Maps \nd(x, A) - d(x, B)" )

    xdim <- dim( x )
    loc <- cbind( rep( 1:xdim[ 1 ], xdim[ 2 ] ), rep( 1:xdim[ 2 ], each = xdim[ 1 ] ) )

    if( missing( N ) ) N <- max( xdim )

    if( N %% 2 == 0 ) N <- N + 1

    # xcen <- imomenter( x, loc = loc )$centroid
    # ycen <- imomenter( y, loc = loc )$centroid

    xy <- x
    xy[ x == 0 & y > 0 ] <- 1

    xycen <- imomenter( xy )$centroid

    bigDloc <- cbind( rep( 1:N, N ), rep( 1:N, each = N ) )

    bigDcen <- rep( ( (N - 1) / 2 ) + 1, 2 )

    cendiff <- bigDcen - xycen

    xloc <- loc[ c( as.logical( x ) ), ]
    yloc <- loc[ c( as.logical( y ) ), ]

    X <- Y <- matrix( 0, N, N )

    idX <- xloc + matrix( cendiff, dim( xloc )[ 1 ], 2, byrow = TRUE )
    idY <- yloc + matrix( cendiff, dim( yloc )[ 1 ], 2, byrow = TRUE )

    goodIDx <- idX >= 1 & idX <= N
    goodIDy <- idY >= 1 & idY <= N

    if( !all( goodIDx ) ) {

	warning( "censqdelta: centering pushes observed data outside of new domain.  Removing some data.  Maybe choose larger N?" )
	if( !any( goodIDx ) ) stop( "censqdelta: No observed data remains after centering in this domain." )
	idX <- idX[ goodIDx[,1] & goodIDx[,2], ]

    }

    if( !all( goodIDy ) ) {

	warning( "censqdelta: centering pushes forecast data outside of new domain.  Maybe choose larger N?" )
	if( !any( goodIDy ) ) stop( "censqdelta: No forecast data remains after centering in this domain." )
        idY <- idY[ goodIDy[,1] & goodIDy[,2], ]

    }

    X[ idX ] <- 1
    Y[ idY ] <- 1

    # dX <- as.im( X )
    # dY <- as.im( Y )

    # dX <- solutionset( dX > 0 )
    # dY <- solutionset( dY > 0 )

    # dX <- distmap( dX )
    # if( is.finite( const ) ) dX$v[ dX$v > const ] <- const
    # dY <- distmap( dY )
    # if( is.finite( const ) ) dY$v[ dY$v > const ] <- const

    # image( X + Y, col = c( "white", "lightblue", "darkblue" ) )
    # title( "Binary A and B superimposed\n centered on new square grid" )
    # image.plot( abs( as.matrix( dX ) - as.matrix( dY ) ), main = "Abs. Diff. in Distance Maps of\nnew grid: d(x, A) - d(x, B)" )
    
    # contour( Y, add = TRUE ) 
    # image( Y, col = c( "white", "darkblue" ) )
    # contour( X, add = TRUE )

    X <- as.im( X )
    Y <- as.im( Y )

    X <- solutionset( X > 0 )
    Y <- solutionset( Y > 0 )

    # par( mfrow = op$mfrow )

    return( deltametric( X, Y, p = p, c = const ) )

} # end of 'modBadder' function.

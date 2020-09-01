warpTps <- function( p1, B, ... ) {

    return( B %*% p1 )

} # end of 'warpTps' function.

warpTpsMatrices <- function( p0, s, ... ) {

    # 'p0' is a matrix of control locations in d >= 2 dimensions.
    # 's' is a matrix of all locations in the image in d dimensions.

    p0dim <- dim( p0 )
    nc <- p0dim[ 1 ]
    d  <- p0dim[ 2 ]

    # Radial basis functions between all locations and control.
    Ls <- warpingRadialBasis( s, p0 )
    Ls <- cbind( Ls, 1, s )

    # Radial basis functions between control locations.
    L <- warpingRadialBasis( p0 )
    L <- cbind( L, 1, p0 )

    # Calculating the warping matrices.
    L <- rbind( L, c( rep( 1, nc ), 0, 0, 0 ), cbind( t( p0 ), matrix( 0, 2, 3 ) ) )
    iL <- solve( L )

    B <- Ls %*% ( iL %*% rbind( diag( nc ), matrix( 0, 3, nc ) ) )

    return( list( L = L, iL = iL, B = B ) )

} # end of 'warpTpsMatrices' function.

warpingRadialBasis <- function( x, y = x ) {

    d <- dim( x )[ 2 ]

    if( dim( y )[ 2 ] != d ) stop( "warpingRadialBasis: invalid y argument." )

    h <- rdist( x, y )

    if( d %% 2 == 0 ) {

        L <- h^d * log( h )
        L[ h == 0 ] <- 0

    } else L <- h^d

    return( L )

} # end of 'warpingRadialBasis' function.



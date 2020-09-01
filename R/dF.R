dF <- function( x ) {

    xdim <- dim( x )
    Kx <- Ky <- matrix( 0, 5, 5 )
    Kx[ 3, 1 ] <- -1
    Ky[ 1, 3 ] <- -1
    Kx[ 3, 3 ] <- 1
    Ky[ 3, 3 ] <- 1

    Kx <- 0.5 * Kx
    Ky <- 0.5 * Ky

    dFx <- kernel2dsmooth( x, K = Kx )
    dFy <- kernel2dsmooth( x, K = Ky )

    return( list( dFx = dFx, dFy = dFy ) )

} # end of 'dF' function.

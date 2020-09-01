prinwarp <- function( x, ... ) {

    if( !is.element( class( x ), c( "warped", "iwarped" ) ) ) stop( "prinwarp: invalid x argument. Must be of class warped or iwarped." )
    if( is.null( x$theta ) ) stop( "prinwarp: x is missing theta component." )
    else if( is.null( x$theta$iL ) ) stop( "prinwarp: theta component of x must have an iL matrix." )

    out <- list()

    theCall <- match.call()

    data.name <- deparse( substitute( x ) )

    k <- nrow( x$p0 )

    Be <- x$theta$iL[ 1:k, 1:k ]

    G <- eigen( Be )

    out$Be <- Be

    out$eigenBe <- G

    out$id <- 1:( k - 3 )

    prinw <- G$vectors %*% t( x$theta$B )
    Lam <- matrix( G$values, ncol = 1 )
    hold <- t( x$p1 ) %*% Lam
    partialw <- matrix( c( hold ), 2, k ) %*% G$vectors %*% prinw

    out$principal.warps <- prinw

    out$partial.warps <- partialw

    out$partial.warp.scores <- t( x$p1 ) %*% G$vectors

    attr( out, "data.name" ) <- data.name
    attr( out, "call" ) <- theCall

    class( out ) <- "prinwarped"

    return( out )

} # end of 'prinwarp' function.

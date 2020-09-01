warpTpsConverter <- function( B.new, B.old, p1.old ) {

    #
    # Convert an old set of 1-energy control points to a
    # new set to correspond with a new set of 0-energy
    # control points for Tps warp functions.
    #

    BnBn <- t( B.new ) %*% B.new
    BnBo <- t( B.new ) %*% ( B.old %*% p1.old )

    p1.new <- solve( BnBn, BnBo )

    return( p1.new )

} # end of 'warpTpsConverter' function.

equal.axis <- function(x, y) {

	pin <- par()$pin

	xr <- x[ 2. ] - x[ 1. ]
	yr <- y[ 2. ] - y[ 1. ]

	my <- mean( y )
	mx <- mean( x )

	r2 <- xr / ( pin[ 1. ] )
	r1 <- yr / ( pin[ 2. ] )

	if(r1 < r2) {

		yr <- (yr * r2) / r1
		y[ 1. ] <- my - yr / 2.
		y[ 2. ] <- my + yr / 2.

	} else {

		xr <- (xr * r1) / r2
		x[ 1. ] <- mx - xr / 2.
		x[ 2. ] <- mx + xr / 2.

	}

	return( list( x = x, y = y ) )

} # end of 'equal.axis' function.


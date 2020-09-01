Q <- function( p, p0, s, Im0, Im1, imethod, B, beta, Cmat = NULL, do.gr = FALSE, dF = NULL, ... ) {

    # 'p' are the 1-energy (p1) column-stacked locations: nc * d
    #    numeric vector (although, for now, d = 2).  Last value is
    #    the nuisance parameter, sigma_epsilon.
    # 'p0' are the 0-energy control locations: nc X d matrix.
    # 's' are the N X d full set of locations (N >= nc).
    # 'B' the Tps "B" matrix (calculated by 'warpTpsMatrices').
    # 'beta' single numeric giving the penalty term
    #    (chosen a priori by user).
    # 'Im0' k X m (where k * m = N) matrix giving the 0-energy image.
    # 'Im1' k X m matrix giving the 1-energy image.
    # 'imethod' character string naming the interpolation method.
    # 'Cmat' nc X nc precision matrix for penalty function.
    # 'do.gr' logical, should the gradients also be returned.
    # '...' not used.

    np <- length( p )
    nc <- nrow( p0 )

    # p1 <- cbind( p[ 1:nc ], p[ (nc + 1):( np  - 1 ) ] )
    p1 <- cbind( p[ 1:nc ], p[ (nc + 1):np ] )

    # Technically, there is a nuisance parameter, sigma, but it is not really
    # needed if not using a true likelihood (cf. 'Qgauss').  Leaving it here
    # for posterity in case it can be useful later.
    sigma <- sigma2 <- 1

    # Intensity part of the loss function.
    wpts <- warpTps( p1 = p1, B = B )
    Im1.def <- Fint2d( Im1, Ws = wpts, s = s, method = imethod, derivs = do.gr )

    if( do.gr ) {

	# Calculate components of the gradient function.

        # Gradient associated with the interpolation part.
        dX <- Im1.def$dx
        dY <- Im1.def$dy

	# Convert Im1.def back to its form for calculating the
	# objective function.
        Im1.def <- Im1.def$xy

	# Gradient associated with the pixels themselves.
        # This part is done before warping!
        if( is.null( dF ) ) dF <- dF( Im1 )

    } # end of if 'do.gr' stmt.

    ImD <- Im1.def - Im0
    ImD2 <- ImD^2

    if( !do.gr ) {

	ll <- sum( colSums( ImD2 ) ) / sigma2
	if( is.na( ll ) ) ll <- sum( colSums( ImD2, na.rm = TRUE ), na.rm = TRUE ) / sigma2

    } # end of if calculate the objective function stmt.

    # Penalty part of the loss function.
    xdelta <- p1[,1, drop = FALSE] - p0[,1, drop = FALSE]
    ydelta <- p1[,2, drop = FALSE] - p0[,2, drop = FALSE]

    if( do.gr ) {

	dLx = dF$dFx * ( ImD / ( sigma2 ) ) * dX
	dLy = dF$dFy * ( ImD / ( sigma2 ) ) * dY
	dL = c( t( B ) %*% cbind( c( dLx ), c( dLy ) ) )

    } # end of if 'do.gr' stmt.

    if( beta > 0 ) {

	if( is.null( Cmat ) ) stop( "Q: must specify Cmat argument when beta > 0." )

	if( !do.gr ) {

            pen <- beta * ( t(xdelta) %*% Cmat %*% xdelta + 
		    t(ydelta) %*% Cmat %*% ydelta ) / 2 # + nc * log( beta )

	} else dpen <- beta * c( t(xdelta) %*% Cmat, t(ydelta) %*% Cmat )

    } else if( beta == 0 ) pen <- dpen <- 0
    else stop("Q: beta must be non-negative.")

    if( !do.gr ) res <- ll + pen
    else res <- dL + dpen

    return( res )

} # end of 'Q' function.

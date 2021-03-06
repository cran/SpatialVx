Qnonzero <- function( p, p0, s, Im0, Im1, imethod, B, beta, Cmat = NULL, do.gr = FALSE, dF = NULL, binary = FALSE, ... ) {

    # 'p' are log( sigma_e ), the 1-energy (p1) column-stacked locations: nc * d
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

    # if( ( np %% 2 ) != 0 ) stop( "Qgauss: invalid p argument.  Must have odd length as first value is for nuisance parameter, sigma." )

    nc <- nrow( p0 )

    p1 <- matrix( p, ncol = 2 )
    # p <- p[ -1 ]
    # p1 <- cbind( p[ 1:nc ], p[ (nc + 1):(np - 1) ] )

    # Intensity part of the loss function.
    wpts <- warpTps( p1 = p1, B = B )
    Im1.def <- Fint2d( Im1, Ws = wpts, s = s, method = imethod, derivs = do.gr )

    if( do.gr ) {

	# Calculate components of the gradient function.

	# Gradient associated with the interpolation part.
	dX <- Im1.def$dx
	dY <- Im1.def$dy

	# Now convert Im1.def back to what it would be for
	# calculating the objective function.
	Im1.def <- Im1.def$xy

	# Gradient associated with the pixels themselves.
        # This part is done before warping!
	if( is.null( dF ) ) dF <- dF( Im1 )

    } # end of if 'do.gr' stmt.

    id <- ( c( Im0 ) > 0 ) & ( c( Im1.def ) > 0 )
    N <- sum( id, na.rm = TRUE )

    ImD <- Im1.def - Im0
    ImD2 <- ImD^2

    if( binary ) sigma2 <- sigma <- 1
    else {

	sigma2 <- var( c( ImD )[ id ], na.rm = TRUE )
	if( sigma2 < 0.01 ) sigma2 <- 0.01
	sigma <- sqrt( sigma2 )

    }

    if( !do.gr ) {

	kappa <- ( N / 2 ) * log( 2 * pi * sigma ) # part of likelihood not involving control points.
        ll <- sum( c( ImD2 )[ id ], na.rm = TRUE ) / ( 2 * sigma2 ) + kappa

    } # end of if calculating the objective function stmt.

    # Penalty part of the loss function.
    xdelta <- p1[, 1, drop = FALSE ] - p0[, 1, drop = FALSE ]
    ydelta <- p1[, 2, drop = FALSE ] - p0[, 2, drop = FALSE ]

    if( do.gr ) {

	dLx <- dF$dFx * ( ImD / ( sigma2 ) ) * dX
	dLy <- dF$dFy * ( ImD / ( sigma2 ) ) * dY
	dL <- c( t( B )[, id, drop = FALSE ] %*% cbind( c( dLx )[ id ], c( dLy )[ id ] ) )

    } # end of if 'do.gr' stmt.

    if( beta > 0 ) {

	if( is.null( Cmat ) ) stop( "Qnonzero: must specify Cmat argument when beta > 0." )

        if( !do.gr ) {

	    pen <- beta * ( t(xdelta) %*% Cmat %*% xdelta + 
		    t(ydelta) %*% Cmat %*% ydelta ) / 2

	} else dpen <- beta * c( t(xdelta) %*% Cmat, t(ydelta) %*% Cmat )

    } else if( beta == 0 ) pen <- dpen <- 0
    else stop("Qgauss: beta must be non-negative.")

    if( !do.gr ) res <- ll + pen
    else res <- dL + dpen

    return( res )

} # end of 'Qnonzero' function.

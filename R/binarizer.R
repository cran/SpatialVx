##
## Need to re-write 'thresholder' to allow for X > u1 & X <= u2
##
binarizer <- function( X, Xhat, threshold = NULL, rule = c( ">", "<", ">=", "<=", "<>", "><", "=<>", "<>=", "=<>=",
						   "=><", "><=", "=><=" ), value = c( "matrix", "owin" ), ... ) {

	xdim <- dim( X )
	if( !all( xdim == dim( Xhat ) ) ) stop( "binarizer: X and Xhat must have same dimensions." )

	rule <- match.arg( rule )
	value <- match.arg( value )

	if( missing( threshold ) || is.null( threshold ) ) {

		if( any( class( X ) == "owin" ) ) Z <- X
		else if( any( class( X ) == "im" ) ) Z <- solutionset( X > 0 )
		else stop( "binarizer: invalid X argument" )

		if( any( class( Xhat ) == "owin" ) ) Zhat <- Xhat
		else if( any( class( Xhat ) == "im" ) ) Zhat <- solutionset( Zhat > 0 )
		else stop( "binarizer: invalid Xhat argument" )

		return( list( Z, Zhat ) )

	} else {

		dbl <- rule %in% c( "<>", "><", "=<>", "<>=", "=<>=", "=><", "><=", "=><=" )

		if( is.matrix( threshold ) ) {

			dth <- dim( threshold )
		        good <- is.element( dth[ 1 ], 1:2 ) & is.element( dth[ 2 ], 1:2 )
			if( !good ) stop( "gfommer: invalid threshold argument." )

			if( any( dth == 1 ) ) threshold <- c( threshold )
			else if( dbl ) {

				u11 <- threshold[1,1]; u12 <- threshold[2,1]
				# u1 <- c( threshold[1,1], threshold[1,2] )
				u21 <- threshold[2,1]; u22 <- threshold[2,2]
				# u2 <- c( threshold[1,2], threshold[2,2] )

			} else if( !dbl & all( dth == 2 ) ) stop( "binarizer: too many threshold values." )

		} # end of if threshold is a matrix stmts.

		if( is.vector( threshold ) ) {
			
			nu <- length( threshold )
			if( nu == 1 & dbl ) stop( "binarizer: only one threshold with a double rule argument." )
			else if( nu == 1 ) { u11 <- threshold; u12 <-  threshold }
			else if( nu == 2 & !dbl ) { u11 <- threshold[ 1 ]; u12 <- threshold[ 2 ] }
			else if( nu == 2 & dbl ) {
			
				u11 <- threshold[ 1 ]; u12 <- threshold[ 2 ]	
				# u1 <- c( threshold[ 1 ], threshold[ 1 ] )
				u21 <- threshold[ 1 ]; u22 <- threshold[ 2 ]
				# u2 <- c( threshold[ 2 ], threshold[ 2 ] )

			} # end of if else stmt.

		} 

	} # end of if else missing threshold stmts.

	# if( class( X ) == "matrix" ) Z <- im( X )
	if( any( class( X ) == "owin" ) ) X <- as.matrix( X )
	# else if( class( X ) != "owin" ) stop( "binarizer: invalid X argument.  Must be a matrix or owin object." )

	# if( class( Xhat ) == "matrix" ) Zhat <- im( Xhat )
	if( any( class( Xhat ) == "owin" ) ) Xhat <- as.matrix( Xhat )
	# else if( class( Xhat ) != "owin" ) stop( "binarizer: invalid Xhat argument.  Must be a matrix or owin object." )

	if( rule == ">" ) {

		Z <- matrix( as.numeric( X > u11 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat > u12 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z > u11 )
		# Zhat <- solutionset( Z > u12 )

	} else if( rule == "<" ) {

		Z <- matrix( as.numeric( X < u11 ), xdim[ 1 ], xdim[ 2 ] )
	        Zhat <- matrix( as.numeric( Xhat < u12 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z < u11 )
		# Zhat <- solutionset( Z < u12 )

	} else if( rule == ">=" ) {

		Z <- matrix( as.numeric( X >= u11 ), xdim[ 1 ], xdim[ 2 ] )
	        Zhat <- matrix( as.numeric( Xhat >= u12 ), xdim[ 1 ], xdim[ 2 ] )
		# Z    <- solutionset( Z >= u11 )
		# Zhat <- solutionset( Z >= u12 )

	} else if( rule == "<=" ) {

		Z <- matrix( as.numeric( X <= u11 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat <= u12 ), xdim[ 1 ], xdim[ 2 ] )
		# Z    <- solutionset( Z <= u11 )
		# Zhat    <- solutionset( Zhat <= u12 )

	} else if( rule == "<>" ) {

		Z <- matrix( as.numeric( X < u11 & X > u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat < u21 & Xhat > u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z    <- solutionset( Z < u11 & Z > u12 )
		# Zhat    <- solutionset( Zhat < u21 & Zhat > u22 )
	       
	} else if( rule == "><" ) {

		Z <- matrix( as.numeric( X > u11 & X < u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat > u21 & Xhat < u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z    <- solutionset( Z > u11 & Z < u12 )
		# Zhat    <- solutionset( Zhat > u21 & Zhat < u22 )
	       
	} else if( rule == "=<>" ) {

		Z <- matrix( as.numeric( X <= u11 & X > u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat <= u21 & Xhat > u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z    <- solutionset( Z <= u11 & Z > u21 )
		# Zhat <- solutionset( Zhat <= u21 & Zhat > u22 )
	       
	} else if( rule == "<>=" ) {

		Z <- matrix( as.numeric( X < u11 & X >= u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat < u21 & Xhat >= u22 ), xdim[ 1 ], xdim[ 2 ] ) 
		# Z <- solutionset( Z < u11 & Z >= u12 )
		# Zhat <- solutionset( Zhat < u21 & Zhat >= u22 )

	} else if( rule == "=<>=" ) {

		Z <- matrix( as.numeric( X <= u11 & X >= u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat <= u21 & Xhat >= u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z <= u11 & Z >= u12 )
		# Zhat <- solutionset( Zhat <= u21 & Zhat >= u22 )
	       
	} else if( rule == "=><" ) {

		Z <- matrix( as.numeric( X >= u11 & X < u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat >= u21 & Xhat < u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z >= u11 & Z < u12 )
		# Zhat <- solutionset( Zhat >= u21 & Zhat < u22 )
	       
	} else if( rule == "><=" ) {

		Z <- matrix( as.numeric( X > u11 & X <= u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat > u21 & Xhat <= u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z > u11 & Z <= u12 )
		# Zhat <- solutionset( Zhat > u21 & Zhat <= u22 )
	       
	} else if( rule == "=><=" ) {

		Z <- matrix( as.numeric( X >= u11 & X <= u12 ), xdim[ 1 ], xdim[ 2 ] )
		Zhat <- matrix( as.numeric( Xhat >= u21 & Xhat <= u22 ), xdim[ 1 ], xdim[ 2 ] )
		# Z <- solutionset( Z >= u11 & Z <= u12 )
		# Zhat <- solutionset( Zhat >= u21 & Zhat <= u22 )

	} # end of if else rule stmts.

	Z <- im( Z )
	Z <- solutionset( Z > 0 )
	Zhat <- im( Zhat )
	Zhat <- solutionset( Zhat > 0 )

	if( value == "matrix" ) {

		Z <- as.matrix( Z )
		Zhat <- as.matrix( Zhat )

	} # end of if 'value' stmt.

	return( list( Z, Zhat ) )

} # end of 'binarizer' function.

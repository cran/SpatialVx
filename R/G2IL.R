G2IL <- function( X, Xhat, threshold, beta, alpha = 0, rule = ">", ... ) {

	hold <- binarizer( X = X, Xhat = Xhat, threshold = threshold, rule = rule, value = "owin" )
	Z <- hold[[ 1 ]]
	Zhat <- hold[[ 2 ]]
	G    <- Gbeta( X = Z, Xhat = Zhat, beta = beta, alpha = alpha, ... )

	IA <- as.matrix( Z ) == 1
	IB <- as.matrix( Zhat ) == 1

	nA <- sum( IA )
	nB <- sum( IB )
	nAB <- sum( IA & IB )
	N <- prod( dim( X ) )

	if( nA > 0 & nB > 0 ) {
		
		tmp <- qqplot( c( X[ IA ] ), c( X[ IB ] ), plot.it = FALSE )
		mae <- mean( abs( tmp$x - tmp$y ), na.rm = TRUE )
		rmse <- sqrt( mean( ( tmp$x - tmp$y )^2, na.rm = TRUE ) )

	} else if( nA == nB ) mae <- rmse <- 0
	else if( nA > 0 ) mae <- rmse <- max( abs( X ), na.rm = TRUE )
	else if( nB > 0 ) mae <- rmse <- max( abs( Xhat ), na.rm = TRUE )

	out <- ubalancer( ( nA + nB - 2 * nAB ) * ( 1 + mae ), alpha = alpha, beta = beta )

	a <- attributes( G )
	nomen <- names( a$components )
	comps <- c( a$components, rmse )
	names( comps ) <- c( nomen, "RMSE (sorted)" )

	attr( out, "components" ) <- comps
	attr( out, "beta" ) <- beta
	attr( out, "alpha" ) <- alpha
	attr( out, "threshold" ) <- c( threshold, rule )
	attr( out, "data.name" ) <- c( deparse( substitute( X ) ), deparse( substitute( Xhat ) ) )
	# attr( out, "weights" ) <- c( w, 1 - w )

	class( out ) <- "G2IL"

	return( out )

} # end of 'GbetaIL' function.

print.G2IL <- function( x, ... ) {

	a <- attributes( x )
        cat( "Observation (A) = ", a$data.name[ 1 ], "\n" )
	cat( "Model (B) = ", a$data.name[ 2 ], "\n" )
	cat( "Threshold and rule: ")
	cat( as.character( a$rule ), a$threshold, "\n" )
	if( a$alpha != 0 ) cat( "alpha = ", a$alpha, "\n" )
	cat( "beta = ", attributes( x )$beta, "\n" )
	cat( "G2betaIL(A,B) = ", c( x ), "\n" )
	cat( "\n\n", "Component parts and asymmetric Gbeta:\n\n" )
	print( a$components )

} # end of 'print.GbetaIL' function.

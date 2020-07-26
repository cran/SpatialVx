GbetaIL <- function( X, Xhat, threshold, beta, alpha = 0, rule = ">", w = 0.5, ... ) {

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

	if( missing( beta ) ) beta <- N^2 / 2

	if( nA > 0 & nB > 0 ) {
		
		tmp <- qqplot( c( X[ IA ] ), c( X[ IB ] ), plot.it = FALSE )
		rho <- cor( tmp$x, tmp$y )

	} else if( nA == nB ) rho <- 1
	else if( nA > 0 ) rho <- 1 - ( nA / N )
	else if( nB > 0 ) rho <- 1 - ( nB / N )

	rho <- max( rho, 0 )

	out <- w * G + (1 - w) * rho

	a <- attributes( G )
	nomen <- names( a$components )
	comps <- c( a$components, rho )
	names( comps ) <- c( nomen, "theta (sorted)" ) # In the paper, I call rho 'theta'.

	attr( out, "components" ) <- comps
	attr( out, "beta" ) <- beta
	attr( out, "alpha" ) <- alpha
	attr( out, "threshold" ) <- c( threshold, rule )
	attr( out, "data.name" ) <- c( deparse( substitute( X ) ), deparse( substitute( Xhat ) ) )
	attr( out, "weights" ) <- c( w, 1 - w )

	class( out ) <- "GbetaIL"

	return( out )

} # end of 'GbetaIL' function.

print.GbetaIL <- function( x, ... ) {

	a <- attributes( x )
        cat( "Observation (A) = ", a$data.name[ 1 ], "\n" )
	cat( "Model (B) = ", a$data.name[ 2 ], "\n" )
	cat( "Threshold and rule: ")
	cat( as.character( a$rule ), a$threshold, "\n" )
	if( a$alpha != 0 ) cat( "alpha = ", a$alpha, "\n" )
	cat( "beta = ", attributes( x )$beta, "\n" )
	cat( "GbetaIL(A,B) = ", c( x ), "\n" )
	cat( "\n\n", "Component parts and asymmetric Gbeta:\n\n" )
	print( a$components )

} # end of 'print.GbetaIL' function.

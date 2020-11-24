Gbeta <- function( X, Xhat, threshold, beta, alpha = 0, rule = ">", ... ) {

	hold <- binarizer( X = X, Xhat = Xhat, threshold = threshold, rule = rule, value = "owin", ... )
	Z <- hold[[ 1 ]]
	Zhat <- hold[[ 2 ]]

	dA <- distmap( Z )
	dB <- distmap( Zhat )

	IA <- as.matrix( Z )
	IB <- as.matrix( Zhat )
	if( !all( dim( IA ) == dim( IB ) ) ) warning( "Gbeta: X and Xhat have different dimensions." )

	N <- prod( dim( IA ) )

	if( missing( beta ) ) beta <- N^2 / 2

	nA <- sum( IA, na.rm = TRUE )
	nB <- sum( IB, na.rm = TRUE )
	nAB <- sum( IA == 1 & IB == 1, na.rm = TRUE )
	# pA <- nA / N
	# pB <- nB / N

	term1 <- ( nA + nB - 2 * nAB )
	term1b <- term1 / sqrt( beta - alpha )
	term1c <- ( nB - nAB ) / sqrt( beta - alpha )
	term1d <- ( nA - nAB ) / sqrt( beta - alpha )
	# term1 <- min( c( nA, nB ) ) / max( c( 1, nA, nB ) ) # Could be zero

	medAB <- sum( dA * IB, na.rm = TRUE ) / max( c( 1, nB ) ) # Could be zero
	medBA <- sum( dB * IA, na.rm = TRUE ) / max( c( 1, nA ) ) # Could be zero

	# term2 <- max( c( 1e-16, min( c( medAB, medBA ) ) ) ) / max( c( 1e-16, medAB, medBA ) )
	term2 <- medAB * nB
	term2b <- term2 / sqrt( beta - alpha )
	term3 <- medBA * nA
	term3b <- term3 / sqrt( beta - alpha )
	term4 <- ( term2 + term3 ) / sqrt( beta - alpha )

	x <- term1b * term4
	# Ind <- as.numeric( x <= beta + 1 / alpha )

	# out <- wobbler( x = x, alpha = alpha, beta = beta )
	# outAB <- ubalancer( x = ( nB - nAB ) * term2, alpha = alpha, beta = beta )
	# outBA <- ubalancer( x = ( nA - nAB ) * term3, alpha = alpha, beta = beta )
	# out <- ubalancer( x = x, alpha = alpha, beta = beta )
	outAB <- pmax( 1 - ( term1c * term2b - alpha ), 0 )
	outBA <- pmax( 1 - ( term1d * term3b - alpha ), 0 )
	out <- pmax( 1 - ( x - alpha ), 0 )

	res <- c( nA, nB, nAB, term1, medAB, medBA, term2, term3, outAB, outBA )

	names( res ) <- c( "nA", "nB", "nAB", "nA + nB - 2nAB", "medAB", "medBA", "medAB * nB",
			 "medBA * nA", "asymGbetaAB", "asymGbetaBA" )

	attr( out, "components" ) <- res
	attr( out, "beta" ) <- beta
	attr( out, "alpha" ) <- alpha
	if( !missing( threshold ) ) attr( out, "threshold" ) <- c( threshold, rule )
	attr( out, "data.name" ) <- c( deparse( substitute( X ) ), deparse( substitute( Xhat ) ) )
	class( out ) <- "Gbeta"

	return( out )

} # end of 'Gbeta' function.

print.Gbeta <- function( x, ... ) {

	a <- attributes( x )
	cat( "Observation (A) = ", a$data.name[ 1 ], "\n" )
	cat( "Model (B) = ", a$data.name[ 2 ], "\n" )
	if( !is.null( a$threshold ) ) {
		
		cat( "Threshold and rule: ")
        	cat( a$threshold[ 1 ], "(", a$threshold[ 2 ], ")\n" )

	}
	if( a$alpha != 0 ) cat( "alpha = ", a$alpha, "\n" )
	cat( "beta = ", attributes( x )$beta, "\n" )
	cat( "Gbeta(A,B) = ", c( x ), "\n" )
	cat( "\n\n", "Component parts and asymmetric Gbeta:\n\n" )
	print( a$components )

} # end of 'print.Gbeta' function.

# wobbler <- function( x, alpha, beta, eta = 2, plot = FALSE, ..., ylab = "GPM", xlab = "distance index" ) {
# 
# 	res <- exp( -pmax( 1 + alpha * ( x - beta ), 0 )^eta )
# 
# 	if( plot ) plot( x, res, type = "l", ylim = c(0,1), xlab = xlab, ylab = ylab, ... )
# 
# 	return( res )
# 
# } # end of 'wobbler' function.

ubalancer <- function( x, alpha, beta, plot = FALSE, ..., ylab = "GPM", xlab = "distance index" ) {

	res <- pmax( 1 - ( x - alpha ) / ( beta - alpha ), 0 )

	if( plot ) plot( x, res, type = "l", ylim = c(0,1), xlab = xlab, ylab = ylab, ... )

	return( res )

} # end of 'ubalancer' function.

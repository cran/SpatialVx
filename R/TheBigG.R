TheBigG <- function( X, Xhat, threshold, rule = ">", ... ) {

	hold <- binarizer( X = X, Xhat = Xhat, threshold = threshold, rule = rule, value = "owin", ... )
	Z <- hold[[ 1 ]]
	Zhat <- hold[[ 2 ]]

	dA <- distmap( Z )
	dB <- distmap( Zhat )

	IA <- as.matrix( Z )
	IB <- as.matrix( Zhat )
	if( !all( dim( IA ) == dim( IB ) ) ) warning( "TheBigG: X and Xhat have different dimensions." )

	N <- prod( dim( IA ) )

	nA <- sum( IA, na.rm = TRUE )
	nB <- sum( IB, na.rm = TRUE )
	nAB <- sum( IA == 1 & IB == 1, na.rm = TRUE )

	term1 <- nA + nB - 2 * nAB

	medAB <- sum( dA * IB, na.rm = TRUE ) / max( c( 1, nB ) ) # Could be zero
	medBA <- sum( dB * IA, na.rm = TRUE ) / max( c( 1, nA ) ) # Could be zero

	term2 <- medAB * nB
	term3 <- medBA * nA
	term4 <- term2 + term3

	x <- term1 * term4

	out <- x^( 1 / 3 )
	outAB <- ( nB - nAB ) * term2
	outBA <- ( nA - nAB ) * term3
	res <- c( nA, nB, nAB, term1, medAB, medBA, term2, term3, outAB, outBA )

	names( res ) <- c( "nA", "nB", "nAB", "nA + nB - 2nAB", "medAB", "medBA", "medAB * nB",
	 		 "medBA * nA", "asymG.AB", "asymG.BA" )

	attr( out, "components" ) <- res
	if( !missing( threshold ) ) attr( out, "threshold" ) <- c( threshold, rule )
	attr( out, "data.name" ) <- c( deparse( substitute( X ) ), deparse( substitute( Xhat ) ) )
	class( out ) <- "TheBigG"

	return( out )

} # end of 'TheBigG' function.

print.TheBigG <- function( x, ... ) {

	a <- attributes( x )
	cat( "Observation (A) = ", a$data.name[ 1 ], "\n" )
	cat( "Model (B) = ", a$data.name[ 2 ], "\n" )
	if( !is.null( a$threshold ) ) {
		
		cat( "Threshold and rule: ")
        	cat( a$threshold[ 1 ], "(", a$threshold[ 2 ], ")\n" )

	}
	cat( "G(A,B) = ", c( x ), "\n" )
	cat( "\n\n", "Component parts and asymmetric G:\n\n" )
	print( a$components )

	invisible()

} # end of 'print.TheBigG' function.

warper <- function( Im0, Im1, p0, init, s, imethod = "bicubic",
    lossfun = "Q", lossfun.args = list( beta = 0, Cmat = NULL ),
    grlossfun = "defaultQ", lower, upper, verbose = FALSE, ... ) {

    begin.tiid <- Sys.time()

    theCall <- match.call()

    data.name <- c( deparse( substitute( Im0 ) ), deparse( substitute( Im1 ) ) )

    dF <- NULL

    if( verbose ) cat( "\nInitializing Tps Warp Matrices.\n" )

    nc <- nrow( p0 )
    xdim <- dim( Im0 )

    theta <- warpTpsMatrices( p0 = p0, s = s )

    if( ( lossfun == "Q" || lossfun == "Qgauss" || lossfun == "Qnonzero" ) && lossfun.args$beta > 0 && is.null( lossfun.args$Cmat ) ) {

	if( verbose ) {

	    cat( "\nUsing penalized Gaussian loss function with Cmat unspecified.\n" )
	    cat("Using the bending energy matrix for Cmat (the penalty precision matrix).\n" ) 

	} # end of if verbose stmt.

	lossfun.args$Cmat <- theta$iL[ 1:nc, 1:nc ]

    } # end of if using "Q" for the loss function, but Cmat is not specified stmts.

    if( is.element( lossfun, c( "Q", "Qgauss", "Qnonzero" ) ) && missing( lower ) ) {

	lower <- 1
	if( verbose ) cat("\n", "Using ", lower, " for lower in call to nlminb optimization routine.\n" )

    } # end of if 'lossfun' is Q/Qgauss and missing 'lower' stmt.

    if(  is.element( lossfun, c( "Q", "Qgauss", "Qnonzero" ) ) && missing( upper ) ) {

        if( verbose ) cat("\n", "Using domain dimension information for upper in call to nlminb optimization routine.\n" )
	upper <- c( rep( xdim[ 1 ], nc ), rep( xdim[ 2 ], nc ) )

    } # end of if 'lossfun' is Q and upper stmt.

    WarpObjFun <- function( x, loc0, loc, Bmat, im0, im1, imet, lossf, lossfargs, dF ) {

	c( do.call( lossf, c( list( p = x, p0 = p0, s = loc, B = Bmat,
		    Im0 = im0, Im1 = im1, imethod = imet ), lossfargs ) ) )

    } # end of internal warping objective function.

    if( is.element( lossfun, c( "Q", "Qgauss", "Qnonzero" ) ) && ( !is.null( grlossfun ) &&
	(grlossfun == "defaultQ" || grlossfun == "defaultQgauss" || grlossfun == "defaultQnonzero" ) ) ) {

 	if( verbose ) cat( "\n", "Using default gradient function for ", lossfun, " loss.\n" )

	dF <- dF( Im1 )

 	 grfn <- function( x, loc0, loc, Bmat, im0, im1, imet, lossf, lossfargs, dF ) {

             c( do.call( lossf, c( list( p = x, p0 = p0, s = loc, B = Bmat,
                     Im0 = im0, Im1 = im1, imethod = imet, dF = dF, do.gr = TRUE ), lossfargs ) ) )

         } # end of internal warping objective function.

     } else if( !is.null( grlossfun ) ) {

 	grfn <- match.fun( grlossfun )

     } else grfn <- NULL
    # end of if use default gradient for Q stmts.

    if( missing( init ) ) {

	init <- c( p0 )

    } else {

	if( verbose ) cat( "\nOptimizing over objective function.\n" )

	# fit <- try( optimx( par = init, fn = WarpObjFun, gr = grfn, loc0 = p0, loc = s, Bmat = theta$B,
	# 	method = omethod, im0 = Im0, im1 = Im1, imet = imethod, lossf = lossfun,
	# 	lossfargs = lossfun.args, lower = lower, upper = upper, ... ), silent = !verbose ) 

	fit <- try( nlminb( start = init, objective = WarpObjFun, gradient = grfn, loc0 = p0, loc = s, Bmat = theta$B,
		im0 = Im0, im1 = Im1, imet = imethod, lossf = lossfun,
		lossfargs = lossfun.args, dF = dF, lower = lower, upper = upper, ... ), silent = !verbose )

    } # end of if else 'init' stmt.

    if( class( fit ) == "try-error" ) {

	if( verbose ) cat( "\nSorry, but something failed to work.\nReturning error object." )
	return( fit )

    }

    # p1 <- matrix( summary( fit, order = value )[1, 1:(nc * 2) ], ncol = 2 )
    p1 <- matrix( fit$par, nc, 2 )

    wpts <- warpTps( p1 = p1, B = theta$B )

    out <- list()

    out$Im0 <- Im0
    out$Im1 <- Im1
    out$Im1.def <- Fint2d( Im1, Ws = wpts, s = s )
    out$p0 <- p0
    out$p1 <- p1
    out$sigma <- sd( c( out$Im1.def ) )
    out$warped.locations <- wpts
    out$init <- init
    out$s <- s
    out$imethod <- imethod
    out$lossfun <- lossfun
    out$lossfun.args <- lossfun.args
    out$theta <- theta
    out$arguments <- list( ... )
    out$fit <- fit

    out$proc.time <- Sys.time() - begin.tiid
    if( verbose ) print( out$proc.time )

    attr( out, "data.name" ) <- data.name
    attr( out, "call" ) <- theCall

    class( out ) <- "warped"

    return( out )

} # end of 'warper' function.

plot.warped <- function( x, col = c( "gray", tim.colors(64) ), alwd = 1.5, zlim, ... ) {

    par( mfrow = c(2,3) )

    if( missing( zlim ) ) zlim <- range( c( c( x$Im0 ), c( x$Im1 ), c( x$Im1.def ) ), finite = TRUE )

    sq <- seq( 0, 1,, dim( x$s )[ 1 ] )

    image( x$Im0, main = "0-energy field", col = col, zlim = zlim, ... )
    # points( sq[ x$p0[, 1 ] ], sq[ x$p0[, 2 ] ], pch = 19, col = "yellow" )
    image( x$Im1, main = "1-energy field", col = col, zlim = zlim, ... )
    # points( sq[ x$p1[, 1 ] ], sq[ x$p1[, 2 ] ], pch = 19, col = "yellow" )
    image.plot( x$Im0, legend.only = TRUE, col = col, zlim = zlim, ... )
    image.plot( x$Im1 - x$Im0, main = "Error Field", col = tim.colors(64), ... )

    xdim <- dim( x$Im0 ) 
    h <- matrix( sqrt( ( x$warped.locations[,1] - x$s[,1] )^2 + ( x$warped.locations[,2] - x$s[,2] )^2 ),
	xdim[ 1 ], xdim[ 2 ] )
    image.plot( matrix( x$s[,1], xdim[ 1 ], xdim[ 2 ] ),
                matrix( x$s[,2], xdim[ 1 ], xdim[ 2 ] ),
                h, main = "Distance travelled", col = c( "gray", tim.colors(64) ), ... )
    arrows(x$p1[,1], x$p1[,2], x$p0[,1], x$p0[,2], col="magenta", length=0.2, lwd = alwd ) 

    image.plot( x$Im1.def, main = "Deformed 1-energy field", col = col, ... )
    image.plot( x$Im1.def - x$Im0, main = "Error Field\n(after warping)", col = tim.colors(64), ... )
    
    par( mfrow = c(1,1) )

    invisible()

} # end of 'plot.warped' function.
# 
print.warped <- function( x, ... ) {

    a <- attributes( x )

    print( a$data.name )
    print( a$call )

    y <- summary( x, silent = TRUE )

    print( y )

    invisible( y )

} # end of 'print.warped' function.

summary.warped <- function( object, ... ) {

    out <- list( ... )

    rmse0 <- sqrt( mean( (object$Im1 - object$Im0)^2, na.rm = TRUE ) )
    rmse1 <- sqrt( mean( (object$Im1.def - object$Im0)^2, na.rm = TRUE ) )
    redux <- ( (rmse0 - rmse1) / rmse0 ) * 100

    k <- nrow( object$p0 )

    Be <- object$theta$iL[ 1:k, 1:k ]
    p1 <- object$p1

    # Affine part
    A <- object$theta$iL[ (k + 1):(k + 3), 1:k ] %*% p1
    xpart <- A[ , 1 ]
    ypart <- A[ , 2 ]

    minBe <- sum( diag( t( p1 ) %*% Be %*% p1 ) )

    res <- c( rmse0, rmse1, redux, xpart, ypart, minBe )
    names( res ) <- c( "RMSE0", "RMSE1", "% RMSE Reduction", "intercept(x)", "a_x(x)", "a_x(y)", 
			 "intercept(y)", "a_y(x)", "a_y(y)", "Bending Energy" )

    silent <- out$silent
    if( is.null( out$silent ) ) silent <- FALSE 

    if( !silent ) {

	a <- attributes( object )
	print( paste( a$data.name, collapse = " vs ") )

	print( res )

    } 

    out$result <- res

    invisible( out )    

} # end of 'summary.warped' function.



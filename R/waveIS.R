waveIS <- function( x, th = NULL, J = NULL, wavelet.type = "haar", levels = NULL, max.n = NULL,
    smooth.fun = "hoods2dsmooth", smooth.params = NULL, rule = ">=", verbose = FALSE, ... ) {

    UseMethod( "waveIS", x )

}

waveIS.SpatialVx <- function( x, th = NULL, J = NULL, wavelet.type = "haar",
    levels = NULL, max.n = NULL, smooth.fun = "hoods2dsmooth",
    smooth.params = NULL, rule = ">=", verbose = FALSE, ...,
    time.point = 1, obs = 1, model = 1 ) {

    a <- attributes( x )

    if( is.null( th ) ) th <- list( X = cbind( a$thresholds$X[, obs ] ), Xhat = cbind( a$thresholds$Xhat[, model ] ) )

    dat <- datagrabber( x, time.point = time.point, obs = obs, model = model )

    obj <- list( X = dat$X, Xhat = dat$Xhat )

    out <- waveIS( x = obj, th = th, J = J, wavelet.type = wavelet.type,
	    levels = levels, max.n = max.n, smooth.fun = smooth.fun,
	    smooth.params = smooth.params, rule = rule, verbose = verbose, ... )

    nomen <- names( out )

    attributes( out ) <- a
    names( out ) <- nomen
    attr( out, "time" ) <- a$time[ time.point ]
    attr( out, "obs.name" ) <- a$obs.name[ obs ]
    attr( out, "model.name" ) <- a$model.name[ model ]
    attr( out, "rule" ) <- rule
    attr( out, "levels" ) <- levels
    attr( out, "class" ) <- "waveIS"

    return( out )

} # end of 'waveIS.SpatialVx' function.

waveIS.default <- function ( x, th = NULL, J = NULL, wavelet.type = "haar", levels = NULL, max.n = NULL,
    smooth.fun = "hoods2dsmooth", smooth.params = NULL, rule = ">=", verbose = FALSE, ... ) {

    if (verbose) begin.time <- Sys.time()

    if( is.null( th ) ) stop( "waveIS: invalid th argument.  Must specify thresholds." ) 
    if( !is.list( th ) ) stop( "waveIS: invalid th argument.  Must be a list with components X and Xhat." )
    if( !all( c( "X", "Xhat" ) %in% names( th ) ) ) stop( "waveIS: invalid th argument.  Must have components X and Xhat." )

    if( is.list( x ) && all( is.element( c( "X", "Xhat" ), names( x ) ) ) ) obj <- list( x$X, x$Xhat )
    else if( is.list( x ) && length( x ) == 2 ) obj <- x 
    else if( is.array( x ) && dim( x )[ 3 ] == 2 )  obj <- list( x[,, 1], x[,, 2] )
    else stop( "waveIS: invalid x argument." )

    xdim <- dim( obj[[ 1 ]] )
    if( !all( xdim == dim( obj[[ 2 ]] ) ) ) stop( "waveIS: invalid x argument.  Dimensions of two fields must be the same." )

    attr( obj, "xdim" ) <- xdim

    object <- hoods2dPrep( object = obj, levels = levels, max.n = max.n,
        smooth.fun = smooth.fun, smooth.params = smooth.params)

    out <- list()

    thresholds <- th
    q <- dim( thresholds$X )[ 1 ]

    X <- obj[[ 1 ]]
    Y <- obj[[ 2 ]]

    bigN <- prod( xdim )
    Jtry <- log2( xdim )

    if( all( floor( Jtry ) == ceiling( Jtry ) ) ) dyadic <- TRUE
    else dyadic <- FALSE

    if( is.null( J ) ) {

        if (dyadic) J <- min( Jtry, na.rm = TRUE )
        else J <- 4

    }

    if (dyadic) wmeth <- "DWT"
    else wmeth <- "MODWT"

    out$wave.method <- wmeth

    eX <- eY <- MSE <- SS <- matrix(NA, nrow = J, ncol = q)

    Bias <- mserand <- numeric(q) + NA

    if (verbose) cat("\n", "Looping through thresholds = \n")

    for (threshold in 1:q) {

        if (verbose) cat("\nThreshold ", threshold, "\n")

        Xbin <- thresholder( obj[[ 1 ]], type = "binary",
	    th = th$X[ threshold ], rule = rule )

	Ybin <- thresholder( obj[[ 2 ]], type = "binary",
	    th = th$Xhat[ threshold ], rule = rule )

        s <- sum(colSums(Xbin, na.rm = TRUE), na.rm = TRUE)
        B <- sum(colSums(Ybin, na.rm = TRUE), na.rm = TRUE) / s
        s <- s / bigN

        Bias[threshold] <- B

        MSE.random <- B * s * (1 - s) + s * (1 - B * s)

        mserand[threshold] <- MSE.random

        if (dyadic) {

            wv.X <- dwt.2d(Xbin, wf = wavelet.type, J = J)
            wv.Y <- dwt.2d(Ybin, wf = wavelet.type, J = J)
            wv.diff <- dwt.2d(Ybin - Xbin, wf = wavelet.type, J = J)

        } else {

            wv.X <- modwt.2d(Xbin, wf = wavelet.type, J = J)
            wv.Y <- modwt.2d(Ybin, wf = wavelet.type, J = J)
            wv.diff <- modwt.2d(Ybin - Xbin, wf = wavelet.type, J = J)

        }

        if( threshold == 1 ) lnames <- makeWaveNames( J = J )

        eX[, threshold ] <- energizer( wv.X, lnames = lnames, J = J )
        eY[, threshold ] <- energizer( wv.Y, lnames = lnames, J = J )

        for (j in 1:J) {

            hold <- detailer(wv.diff, level = j, which.space = "wavelet", 
                trans.type = wmeth, lnames = lnames, J = J)

            N <- length( hold[ !is.na( hold ) ] )

            mse.ui <- sum( colSums( hold, na.rm = TRUE ), na.rm = TRUE ) / N
            MSE[j, threshold ] <- mse.ui
            SS[j, threshold ] <- 1 - ( mse.ui * (J + 1) ) / MSE.random

        }

    }

    if (verbose) cat("\n", "Finished computing DWTs for each threshold.\n")

    out$J <- J
    out$EnVx <- eX
    out$EnFcst <- eY
    out$MSE <- MSE
    out$Bias <- Bias
    out$SS <- SS
    out$MSE.random <- mserand

    attr( out, "rule" ) <- rule
    attr( out, "levels" ) <- levels
    attr( out, "thresholds" ) <- th
    attr( out, "qs" ) <- 1:( dim( th$X )[ 1 ] )
    attr( out, "xdim" ) <- xdim
    attr( out, "map" ) <- FALSE

    class(out) <- "waveIS"

    return(out)

} # end of 'waveIS.default' function.

plot.waveIS <- function( x, main1 = "X", main2 = "Y",
    which.plots = c( "all", "mse", "ss", "energy" ),
    level.label = NULL, ... ) {

    op <- par()

    a <- attributes( x )

    x <- summary( x, silent = TRUE )
 
    J <- x$J

    xdim <- a$xdim

    thresholds <- a$qs

    q <- length(thresholds)

    dyadic <- all( floor( log2( xdim ) ) == ceiling( log2( xdim ) ) )

    if(dyadic) levels <- paste("2^", J-(1:J), sep="")
    else {

       levels <- 1:J
       if(is.null(level.label)) level.label <- "Level"

    } # end of if else 'dyadic stmts.

    if( any( is.element( c( "all", "mse" ), which.plots ) ) ) {

        par( mfrow = c( 2, 1 ), mar = c( 5.1, 4.1, 4.1, 5.1 ) )

        image( x$MSE, col = c( "grey", tim.colors(64) ), main = paste( "MSE (", main1, " vs ", main2, ")", sep="" ),
		xlab = level.label, ylab = paste( "threshold (", a$units, ")", sep=""), axes = FALSE )

        axis( 1, at = seq( 0, 1,, J), labels = levels )
        axis( 2, at = seq( 0, 1,, q), labels = thresholds )

        image.plot( x$MSE, col = c( "grey", tim.colors(64) ), legend.only = TRUE )

        MSEperc <- x$MSEperc

        image( MSEperc, col = c( "grey", tim.colors(64) ), zlim = c( 0, 100 ), main = "MSE %", xlab = level.label,
                                ylab = paste( "threshold (", a$units, ")", sep=""), axes = FALSE )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot(MSEperc, col=c("grey",tim.colors(64)), zlim=c(0,100),legend.only=TRUE)

    } # end of if do MSE stmts.

    if( any( is.element( c( "all", "ss" ), which.plots ) ) ) {

        par( mfrow = c( 1, 1 ), mar = c( 5.1, 4.1, 4.1, 5.1 ) )

        image( x$SS, col = c( "grey", tim.colors(64) ),
	    main = paste( "IS Skill Score (", main1, " vs ", main2, ")", sep=""), xlab = level.label,
            ylab = paste( "threshold (", a$units, ")", sep = "" ), axes = FALSE )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$SS, col = c( "grey", tim.colors(64) ), legend.only = TRUE )

   } # end of if do SS stmts.

    if( any( is.element( c( "all", "energy" ), which.plots ) ) ) {

        zl <- range( c( c( x$EnVx ), c( x$EnFcst ) ), finite = TRUE )

        par( mfrow = c(3, 2), mar = c( 5.1, 4.1, 4.1, 5.1 ) )

        image( x$EnVx, col = c( "grey", tim.colors(64) ), main = paste( main1, " Energy", sep="" ),
	    xlab = level.label, ylab = paste( "threshold (", a$units, ")", sep = "" ), axes = FALSE,
	    zlim = zl )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$EnVx, col = c( "grey", tim.colors(64) ), legend.only = TRUE, zlim = zl )

        image( x$EnFcst, col = c( "grey", tim.colors(64) ), zlim = zl, main = paste( main2, " Energy", sep="" ),
	    xlab = level.label, ylab = paste( "threshold (", a$units, ")", sep="" ), axes = FALSE )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$EnFcst, col = c( "grey", tim.colors(64) ), zlim = zl, legend.only = TRUE )

        image( x$EnVx.perc, col = c( "grey", tim.colors(64) ), main = paste( main1, " Energy %", sep="" ),
	    xlab = level.label, ylab = paste( "threshold (", a$units, ")", sep="" ), axes = FALSE, zlim = c(0, 100) )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$EnVx.perc, col = c( "grey", tim.colors(64) ), legend.only = TRUE, zlim = c( 0, 100 ) )

        image( x$EnFcst.perc, col =c( "grey", tim.colors(64) ), zlim = c(0, 100), main = paste( main2, " Energy %", sep=""),
	    xlab = level.label, ylab = paste( "threshold (", a$units, ")", sep="" ), axes = FALSE )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$EnFcst.perc, col = c( "grey", tim.colors(64) ), zlim = c(0, 100), legend.only = TRUE )

        image( x$EnRelDiff, col = c( "grey", tim.colors(64) ),
	    main = paste( "Energy Relative Difference (", main1, " vs ", main2, ")", sep="" ), xlab = level.label,
	    ylab = paste( "threshold (", a$units, ")", sep="" ), axes = FALSE, zlim = c(-1, 1) )

        axis( 1, at = seq( 0, 1,, J ), labels = levels )
        axis( 2, at = seq( 0, 1,, q ), labels = thresholds )

        image.plot( x$EnRelDiff, col = c( "grey", tim.colors(64) ), zlim = c( -1, 1 ), legend.only = TRUE )

   } # end of if plot energy results.

    par( mfrow = op$mfrow, mar = op$mar )

   invisible()

} # end of 'plot.waveIS' function.

summary.waveIS <- function( object, ... ) {

   x <- object
   a <- attributes( object )

   args <- list( ... )

   if( !is.null( args$silent ) ) silent <- args$silent
   else silent <- FALSE

   J <- x$J

   MSE <- x$MSE
   colnames( MSE ) <- a$qs

   MSEu <- apply( MSE, 2, sum, na.rm = TRUE )

   MSEperc <- t( t( MSE ) / MSEu ) * 100

   SSu <- matrix( 1 - MSEu / x$MSE.random, nrow = 1 )
   SS <- x$SS
   colnames( SS ) <- a$qs

   EnVx <- x$EnVx
   colnames( EnVx ) <- a$qs

   EnFcst <- x$EnFcst
   colnames( EnFcst ) <- a$qs

   EnVx.u <- apply( EnVx, 2, sum, na.rm = TRUE )

   EnFcst.u  <- apply( EnFcst, 2, sum, na.rm = TRUE )

   EnVx.perc <- t( t( EnVx ) / EnVx.u ) * 100
   EnFcst.perc <- t( t( EnFcst ) / EnFcst.u ) * 100

   EnRelDiff <- ( EnFcst - EnVx ) / ( EnFcst + EnVx )

   Bias <- matrix( x$Bias, nrow = 1 )
   colnames(Bias) <- a$qs

   xdim <- a$xdim

   dyadic <- all( floor( log2( xdim ) ) == ceiling( log2( xdim ) ) )

   if( dyadic ) level.names <- paste( "2^", J - ( 1:J ), sep = "" )
   else level.names <- as.character( 1:J )

   rownames( MSE ) <- rownames( EnVx ) <- rownames( SS ) <- rownames( MSEperc ) <-
        rownames( EnFcst ) <- rownames( EnVx.perc ) <- rownames( EnFcst.perc ) <-
        rownames( EnRelDiff ) <- level.names

    if(!silent) {

        cat("\n", x$Fcst.name, " model field compared with verification field ", x$Vx.name, "\n")

        cat("MSE by threshold (columns) and level (rows).\n")

        print(round(MSE, digits=4))

        cat("\n", "MSE% by threshold (columns) and level (rows).\n")

        print(round(MSEperc, digits=0))

        cat("\n", "Field skill score by threshold.\n")

        print(round(SSu, digits=4))

        cat("\n", "IS Skill Score\n")

        print(round(SS, digits=4))

        cat("\n", "Verification energy\n")

        print(round(EnVx, digits=4))

        cat("\n", "% Verification energy\n")

        print(round(EnVx.perc,digits=0))

        cat("\n", "Forecast energy\n")

        print(round(EnFcst, digits=4))

        cat("\n", "% Forecast energy\n")

        print(round(EnFcst.perc, digits=0))

        cat("\n", "Energy Relative Difference\n")

        print(round(EnRelDiff, digits=4))

        cat("\n", "Frequency Bias\n")

        print(round(Bias, digits=4))

   } # end of if '!silent' stmts.

   invisible( c( x, list( MSEu = MSEu, MSEperc = MSEperc, SSu = SSu,
	EnVx.u = EnVx.u, EnFcst.u = EnFcst.u, EnVx.perc = EnVx.perc,
	EnFcst.perc = EnFcst.perc, EnRelDiff = EnRelDiff ) ) )

} # end of 'summary.waveIS' function.


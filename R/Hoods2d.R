kernel2d <-
function( x, n, W=NULL, K=NULL, X=NULL, xdim=NULL, Nxy=NULL, setup=FALSE, verbose=FALSE) {
   ##
   ## Function to compute the boxcar kernel convolution smooth on a 2d image field.
   ##
   ## Arguments:
   ##
   ## 'x' 'k X r' matrix giving the grid point (pixel) values to be smoothed.
   ## 'n' numeric giving the neighborhood length.  If 'n' is not an integer, it will be converted by 'floor'.  If it is even,
   ##	one will be subtracted from it.  If it is less than 1, the function will error.  If it is equal to 1, 'x' will be
   ##	returned as is (no smoothing applied) in order to ensure consistency when neighborhood lengths of 1 are desired as, e.g.,
   ##	part of an analysis of several threshold/level combinations.
   ## 'W' (optional) '2k X 2r' matrix of Fourier transformed kernel weights.  If NULL, these will be found, and
   ##	if 'setup' is TRUE, these are what is returned, rather than a smoothed field.
   ## 'K' (optional) '2k X 2r' matrix giving the kernel to be applied to the smoothing (e.g., if one wanted to
   ##	apply a different kernel other than the boxcar kernel, then these would be passed in here.  If
   ##	NULL, then the boxcar kernel is applied (i.e., nearest neighbors with a neighborhood length of 'n'.
   ## 'X' (optional) '2k X 2r' matrix giving the Fourier transformed expanded image of 'x'.  This is for when
   ##	multiple calls to this program are required (e.g., for different neighborhood lengths).  If supplied, it
   ##	will save one FFT in computation.  If 'W' is also supplied, then it will save two FFT's.
   ## 'xdim' (optional) numeric vector of length 2 giving the dimensions of 'x'.  If NULL, these are found here.
   ## 'Nxy' (optional) numeric giving the total number of grid points in 'x'.  If NULL, it will be calculated here.
   ## 'setup' logical, should just the Fourier transformed smoothing weights be calculated and returned?  If TRUE,
   ##	then 'W' should be NULL, or you will not get what you think.
   ## 'verbose' logical, should progress information be printed to the screen?
   ##
   ## Details:
   ##
   ## This 2-d spatial kernel smoother applies the smoother of Roberts and Lean (2008) to a spatial field.  Specifically,
   ## If X is a matrix of grid points, then the returned field, denoted by Ebert (2008) as <X>_s, is a smoothed field such
   ## that the value at each grid point '<X>_s[i,j]' is given by:
   ##
   ##	<X>_s[i,j] = sum_k sum_l X[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2]*K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2],
   ##
   ## where n is the neighborhood length, k,l = 1, ..., n, and K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2] (default) is
   ## constant, and equal to 1/n^2.
   ##
   ## In order to be fast, loops are avoided completely.  Instead, the convolution theorem is applied with a Fast Fourier
   ## Transform (FFT).  If the weights 'W' are supplied, then you will save one FFT in computation time.  See Romanyuk (2005)
   ## for more information on this approach.
   ##
   ## The convolution theorem says that the Fourier transform of a convolution between two functions f and g is equal to the
   ## product of the Fourier transformed functions.  That is, if F denotes the Fourier transform, and * the convolution
   ## operator, F( f*g ) = k F(f)F(g), where 'k' is a scaling factor.  The neighborhood smooth is given by a convolution
   ## between the field and a boxcar kernel (i.e., a square around a point with constant value 1/n^2).  Because of the FFT,
   ## this enables a fast way to compute this convoultion.
   ##
   ## In order to zero-pad the field, and perform a cyclic convolution, it is necessary to expand the field, 'x', and re-arrange
   ## the boxcar kernel (or else it will not be centered on the points).
   ## 
   ## References:
   ##
   ## Ebert E. E., 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework. Meteorol. Appl.,
   ##   15:51--64. DOI: 10.1002/met.25
   ##   Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf
   ##
   ## Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations
   ##   from high-resolution forecasts of convective events.  Mon. Wea. Rev., 136:78--97.
   ##   DOI: 10.1175/2007MWR2123.1.
   ##
   ## Romanyuk, Y. A., 2005: "Discrete Signal Transformations."  Basics of Digital Signal Processing.  Part I.  Chapter 3.
   ##
   ## Value: If 'setup' is FALSE, then the smoothed field is returned, which is a numeric matrix with the same dimension as 'x'.
   ##	Otherwise, a '2k X 2r' matrix of Fourier transformed kernel weights is returned for future use with this function.
   ##
   n <- floor(n)
   if( n%%2 == 0) n <- n-1
   if( n==1) return( x)
   if( n < 1) stop("kernel2d: n must be an odd whole number >= 1.")
   m <- (n-1)/2
   n2 <- n^2
   if( is.null( xdim)) xdim <- dim( x)
   if( is.null( Nxy)) Nxy <- prod( xdim)
   out <- matrix( 0, 2*xdim[1], 2*xdim[2])
   out[1:xdim[1],1:xdim[2]] <- x
   out[ is.na( out)] <- 0
   if( is.null(W)) {
	if( is.null( K)) {
	   if( verbose) cat("constucting the kernel matrix.\n")
	   Ktmp <- matrix( 1/n2, n, n)
	   Ktop <- cbind( Ktmp[(m+1):n, (m+1):n], matrix(0, n-m, 2*xdim[2]-n), Ktmp[(m+1):n, 1:m])
	   Kmid <- matrix(0, 2*xdim[1]-n, 2*xdim[2])
	   if( m == 1) Kbot <- c( Ktmp[1:m,(m+1):n], rep(0, 2*xdim[2]-n), Ktmp[1:m,1:m])
	   else Kbot <- cbind( Ktmp[1:m,(m+1):n], matrix(0, m, 2*xdim[2]-n), Ktmp[1:m,1:m])
	   K <- rbind( Ktop, Kmid, Kbot)
	   if( verbose) cat("The kernel matrix has been constructed.  It has dimension", dim( K), "\n")
   	}
	if( verbose) cat("Finding the FFT of the kernel matrix.\n")
   	W <- fft( K)/(4*Nxy)
	if( verbose) cat("FFT of kernel matrix found.\n")
	if( setup) return( W)
   } # end of if 'is.null(W)' stmts.
   if( verbose) cat("Performing the convolution.\n")
   if( !is.null( X)) out <- Re( fft( X*W, inverse=TRUE))[1:xdim[1], 1:xdim[2]]
   else out <- Re( fft( fft( out)*W, inverse=TRUE))[1:xdim[1], 1:xdim[2]]
   if( verbose) cat("The convolution has been carried out.\n")
   return( zapsmall(out))
} # end of 'kernel2d' function.

hoods2d <-
function( obj, which.methods = c("upscaling", "mincvr", "multi.event", "fuzzy", "joint", "fss", "pragmatic"), verbose = FALSE) {
   ##
   ## Function to calculate several of the neighborhood smoothing methods (see Gilleland et al., 2009) for spatial forecast verification
   ## described in Ebert (2008).
   ##
   ## Arguments:
   ##
   ## 'obj' a list object output from the 'hoods2dPrep' function.
   ## 'which.methods' character vector stating which methods are to be executed.  Default is for the entire list to be executed.  See Details
   ##	section for specific option information.
   ##
   ## Details:
   ##
   ## This function uses an object from the function 'hoods2dPrep' that includes some of the options utilized by this function, including
   ## the thresholds and neighborhood lengths (levels) to be used.
   ##
   ##  The neighborhood methods (cf. Ebert 2008; Gilleland et al., 2009) apply a smoothing filter to either the raw forecast
   ##	(and possibly also the observed) field(s) or to the binary counterpart(s) determined by thresholding.  The specific smoothing filter
   ##	applied for these methods could be of any type, but those described in Ebert (2008) are generally taken to be "neighborhood" filters.
   ##	In some circles, this could be referred to as convolution filters with a boxcar kernel.  Because the smoothing filter can be represented
   ##	this way, it is possible to use the convolution theorem with the Fast Fourier Transform (FFT) to perform the neighborhood smoothing
   ##	very quickly. The particular approach used here "zero pads" the field, and replaces all missing values with zero as well.  If any missing
   ##	values are introduced after the convolution, they are removed.
   ##
   ## If zero-padding is undesireable, then two options are available.  1. Give a subset to the 'hoods2dPrep' function (e.g., some tile within the
   ## domain) so that the final statistics are calculated only on this subset.  2. Extrapolate the fields before applying this function (and 'hoods2dPrep').
   ## In the case of 2, you might want to also give it the subset (e.g., to give it only the original un-extrapolated fields).
   ##
   ## To simplify the notation for the descriptions of the specific methods employed here, the notation of Ebert (2008) is adopted.  That is, if
   ## a method uses neighborhood smoothed observations (NO), then the observed field is given to be <X>s, or the associated binary field, by <Ix>s.  
   ## Otherwise, if the observation field is not smoothed (denoted by SO in Ebert, 2008), then simply X or Ix are used.  Similarly, for the forecast
   ## field, <Y>s or <Iy>s are used for neighborhood smoothed forecast fields (NF).  If it is the indicator fields that are smoothed (i.e., the original
   ## fields are first thresholded), then the resulting fields are indicated by <Px>s and <Py>s, resp.  Below, NO-NF indicates that a neighborhood
   ## smoothed observed field (<Yx>s, <Ix>s, or <Px>s) is compared with a neighborhood smoothed forecast field, and SO-NF indicates that the observed
   ## field is not smoothed.
   ## 
   ## Options for 'which.methods' include:
   ##           "upscaling": (NO-NF) the rmse is calculated between <X>s and <Y>s, and indicator fields <Ix>s and <Iy>s are created using the 'thresholds'
   ##                           given by the argument 'obj' applied to <X>s and <Y>s.  Bias, threat score (ts) and equitable threat score (ets, which is
   ##				also known as Gilbert Skill Score, GSS) are calculated on <Ix>s and <Iy>s.
   ##           "mincvr": (NO-NF) Compares <Ix>s and <Iy>s by thresholding the neighborhood smoothed fields <Px>s and <Py>s (i.e., smoothed versions of
   ##				Ix and Iy) to obtain <Ix>s and <Iy>s.  Indicator fields <Ix>s and <Iy>s are created by thresholding <Px>s and <Py>s by
   ##				frequency threshold 'Pe' given by the 'obj' argument.  Scores calculated between <Ix>s and <Iy>s include:
   ##				probability of detecting an event (pod, also known as the hit rate), false alarm ratio (far) and ets.
   ##           "multi.event": (SO-NF) Multi-event Contingency Table, which compares the binary observed field Ix against the smoothed forecast indicator
   ##                           field, <Iy>s, which is determined similarly as for "mincvr" (i.e., using Pe as a threshold on <Py>s).  The hit rate and
   ##                           false alarm rate (F) are calculated.
   ##           "fuzzy": (NO-NF) The fuzzy logic approach compares <Px>s to <Py>s by creating a new contingency table where hits = sum_i min(<Px>s_i,<Py>s_i),
   ##				misses = sum_i min(<Px>s_i,1-<Py>s_i), false alarms = sum_i min(1-<Px>s_i,<Py>s_i), and
   ##				correct negatives = sum_i min(1-<Px>s_i,1-<Py>s_i).
   ##		"joint": (NO-NF) Similar to "fuzzy" above, but hits  = sum_i prod(<Px>s_i,<Py>s_i), misses = sum_i prod(<Px>s_i,1-<Py>s_i), false alarms =
   ##				sum_i prod(1-<Px>s_i,<Py>s_i), and correct negatives = sum_i prod(1-<Px>s_i,1-<Py>s_i).
   ##		"fss": (NO-NF) Compares <Px>s and <Py>s directly using a Fractions Brier and Fractions Skill Score (FBS and FSS, resp.), where FBS is
   ##				the mean square difference between <Px>s and <Py>s, and the FSS is one minus the FBS divided by a reference MSE given by
   ##				the sum of the sum of squares of <Px>s and <Py>s individually, divided by the total.
   ##		"pragmatic": (SO-NF) Compares Ix with <Py>s, calculating the Brier and Brier Skill Score (BS and BSS, resp.), where the reference forecast
   ##				used for the BSS is taken to be the mean square error between the base rate and Ix.
   ##
   ## References:
   ##
   ## Ebert EE, 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.
   ##   Meteorol. Appl., 15:51--64. DOI: 10.1002/met.25 Available at
   ##   http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf
   ##
   ## Gilleland, E., D. Ahijevych, B.G. Brown, B. Casati, and E.E. Ebert, 2009. Intercomparison of Spatial
   ##   Forecast Verification Methods. Wea. Forecasting, 24, 1416--1430, DOI: 10.1175/2009WAF2222269.1.
   ##
   ## See Also: 'fss2dPlot' (can be called on an object output from this function simply by 'plot'), 'kernel2d', 'hoods2dPrep', 'fft'.
   ##
   ## Value: A list object of class 'hoods2d' with components determined by 'which.methods' (named by the elements of 'which.methods').  Each
   ##	component is itself a list object containing relevant components to the given method.  For example, hit rate is abbreviated pod here, and if this
   ##	is an output for a method, then there will be a component named pod (all lower case).  The Gilbert Skill Score is abbreviated 'ets' (equitable threat
   ##	score; again all lower case here).
   ##

   if( verbose) begin.time <- Sys.time() 
   thresholds <- obj$thresholds
   q <- dim( thresholds)[1]
   levels <- obj$levels
   l <- length( levels)
   Y <- get( obj$Fcst.name)
   X  <- get( obj$Obs.name)
   xdim <- obj$xdim
   binmat <- matrix(0, xdim[1], xdim[2])
   outmat <- matrix( NA, l, q)
   out <- hoods2dSetUpLists( which.methods=which.methods, mat=outmat)
   out$prep.object <- as.character( substitute( obj))
   bigN <- obj$Nxy
   sub <- obj$subset
   if( "upscaling" %in% which.methods) {
           if( verbose) cat("Finding FFT of forecast and observed fields for upscaling method.\n")
           fftY <- Fourier2d( Y, xdim=xdim)
           fftX <- Fourier2d( X, xdim=xdim)
        } # end of if 'upscaling' stmts.
   if( verbose) cat("Looping through thresholds.\n")
   for( threshold in 1:q) {
      	if( any( c("mincvr", "mincvr", "multi.event", "fuzzy", "joint", "fss", "pragmatic") %in% which.methods)) {
	   if( verbose) cat("\n", "Setting up binary objects for thresholds = ", thresholds[ threshold,], "\n")
	   Ix <- Iy <- binmat
	   Ix[ X >= thresholds[threshold,2]] <- 1
	   Iy[ Y >= thresholds[threshold,1]] <- 1
        } # end of if find 'Ix' and 'Iy' stmt.
	if( "fss" %in% which.methods) {
	   if( is.null( sub)) f0 <- mean( Ix, na.rm=TRUE)
	   else f0 <- mean( c(Ix)[ subset], na.rm=TRUE)
	   if( threshold==1) {
		out$fss$fss.uniform <- 0.5 + f0/2 
		out$fss$fss.random  <- f0
	   } else {
		out$fss$fss.uniform <- c( out$fss$fss.uniform, 0.5 + f0/2)
		out$fss$fss.random  <- c( out$fss$fss.random, f0)
	   } # end of if else 'threshold' is 1 stmts.
	} # end of if 'fss' method.
	if( verbose) cat( "Looping through levels.\n")
        for( level in 1:l) {
	   if( verbose) cat("Neighborhood length = ", levels[ level], "\n")
	   levelW <- kernel2d( X, n=levels[level], xdim=xdim, Nxy=bigN, setup=TRUE)
	   if( any( c( "mincvr", "multi.event", "fuzzy", "joint", "pragmatic", "fss") %in% which.methods)) {
		if( any( c( "mincvr", "multi.event", "fuzzy", "joint", "fss") %in% which.methods)) sPx <- kernel2d( Ix, levels[level], W=levelW, xdim=xdim, Nxy=bigN)
		sPy <- kernel2d( Iy, levels[level], W=levelW, xdim=xdim, Nxy=bigN)
	   } # end of if any 'mincvr', 'multi', 'fuzzy', 'joint' or 'pragmatic' stmts.
	   if( any( c( "mincvr", "multi.event") %in% which.methods)) {
                sIx <- sIy <- binmat
                sIx[ sPx >= obj$Pe[level]] <- 1
		sIy[ sPy >= obj$Pe[level]] <- 1
		if( "mincvr" %in% which.methods) {
		   tmp <- MinCvg2dfun( sIy=sIy, sIx=sIx, subset=sub)
		   out$mincvr$pod[level,threshold] <- tmp$pod
		   out$mincvr$far[level,threshold] <- tmp$far
		   out$mincvr$ets[level,threshold] <- tmp$ets
		} # end of 'mincvr' part.
		if( "multi.event" %in% which.methods) {
		   tmp <- multicon2dfun( sIy=sIy, Ix=Ix, subset=sub)
		   out$multi.event$pod[level,threshold] <- tmp$pod
		   out$multi.event$f[level,threshold] <- tmp$f
		   out$multi.event$hk[level,threshold] <- tmp$hk
		} # end of 'multi' stmts.
           } # end of 'mincvr'/'multi' stmts.
	   if( any( c("fuzzy", "joint") %in% which.methods)) {
		tmp <- fuzzyjoint2dfun( sPy=sPy, sPx=sPx, subset=sub)
		if( "fuzzy" %in% which.methods) {
		   out$fuzzy$pod[level,threshold] <- tmp$fuzzy$pod
		   out$fuzzy$far[level,threshold] <- tmp$fuzzy$far
		   out$fuzzy$ets[level,threshold] <- tmp$fuzzy$ets
		} # end of if 'fuzzy' stmt.
		if( "joint" %in% which.methods) {
		   out$joint$pod[level,threshold] <- tmp$joint$pod
		   out$joint$far[level,threshold] <- tmp$joint$far
		   out$joint$ets[level,threshold] <- tmp$joint$ets
		} # end of if 'joint' methods.
	   } # end of if 'fuzzy/joint' stmts.
	   if( "fss" %in% which.methods) out$fss$values[level,threshold] <- fss2dfun(sPy=sPy, sPx=sPx, subset=sub)
	   if( "pragmatic" %in% which.methods) {
		tmp <- pragmatic2dfun( sPy=sPy, Ix=Ix, subset=sub)
		out$pragmatic$bs[level,threshold] <- tmp$bs
		out$pragmatic$bss[level,threshold] <- tmp$bss
	   }
	   if( "upscaling" %in% which.methods) {
		sYy <- kernel2d( Y, levels[ level], W=levelW, X=fftY, xdim=xdim, Nxy=bigN)
		sYx <- kernel2d( X, levels[ level], W=levelW, X=fftX, xdim=xdim, Nxy=bigN)
		tmp <- upscale2dfun(sYy=sYy, sYx=sYx, threshold=c(thresholds[threshold,]), subset=sub)
		if( threshold==1) out$upscaling$rmse[level] <- tmp$rmse
		out$upscaling$bias[level,threshold] <- tmp$bias
		out$upscaling$ts[level,threshold] <- tmp$ts
		out$upscaling$ets[level,threshold] <- tmp$ets
	   }
	} # end of for 'level' loop.
   } # end of for 'threshold' loop.
   if( verbose) print( Sys.time() - begin.time)
   class( out) <- "hoods2d"
   return( out)
} # end of 'hoods2d' function.

hoods2dPrep <-
function( Fcst.name, Obs.name, thresholds=NULL, Pe=NULL, levels=NULL, max.n=NULL, subset=NULL,
                                loc=NULL, qs=NULL, units=NULL) {
   ##
   ## Function to set up an object that can be used with all neighborhood type spatial verification functions.
   ##
   ## Arguments:
   ##
   ## 'Fcst.name', 'Obs.name' characters giving the names of the forecast and observed 'k X m' numeric matrics giving the forecast and
   ##	observed fields, resp.
   ## 'thresholds' (optional) numeric vector of length 'q >= 1' (or 'q X 2' matrix) giving the threshold(s) to be used in calculating the
   ##   various statistics.  If different thresholds are to be used for the 'Fcst' field as for the 'Obs' field, then the 'Fcst'
   ##	thresholds are in the first column, and those for the 'Obs' field in the second column (in this case, you may want to provide
   ##	quantile names to 'qs').  If a vector, then the thresholds are taken to be the same for both fields.  If NULL, then thresholds
   ##	are taken to be the quantiles: 0, .10, .25, .33, .50, .66, .75, .90, .95 of each field, resp.  
   ## 'Pe', (optional) single numeric giving the frequency threshold over which to determine binary fields for methods such as minimum coverage.  If
   ##	it is NULL, then it is taken to be the most relaxed requirement (i.e., that an event occurs at least once in a neighborhood) of 'Pe'=1/(n^2),
   ##	where 'nlen' is the length of the neighborhood.
   ## 'levels' (optional) vector of odd integers giving the neighborhood lengths to use.
   ## 'max.n' (optional) single odd numeric giving the maximum neighborhood length to use.  Must be less than 2N-1, where N
   ##   is the length of the longer side of 'Fcst'.  If this is even, 1 is subtracted from it.  Only used if 'levels' is
   ##   NULL.  If both 'levels' and 'max.n' are NULL, then all odd integers from 1 to 2N-1 are used.
   ## 'subset' (optional) integer vector describing which points should be considered within a domain. May simply be a subset of
   ##   points. If NULL, no tiling is performed.  See Details section for more information.
   ## 'loc' (optional) 'k*m X 2' matrix giving the lon/lat coordinates for each data point.
   ## 'qs' (optional) character vector giving the names of the quantiles used for the thresholds.  Only used for the 'fss2dPlot' function.
   ## 'units' (optional) character giving the units.  Only used with 'fss2dPlot'.
   ##
   ## Value: A list object with components: All arguments input with the same names, and:
   ##
   ## 'xdim' numeric vector of length 2 giving the dimensions of each field in the verification set (i.e., these should be the same for both fields
   ##	or these functions won't work).
   ## 'Nxy' single numeric giving the total number of points in each field (i.e., prod( xdim)).
   ##
   out <- list()
   out$Fcst.name <- Fcst.name
   out$Obs.name <- Obs.name
   Fcst <- get( Fcst.name)
   Obs <- get( Obs.name)
   dimF <- dim( Fcst)
   out$xdim <- dimF
   if( any(dimF != dim( Obs))) stop("fss2d: dim of Obs must be the same as dim of Fcst.")
   Nxy <- prod( dimF)
   out$Nxy <- Nxy
   if( is.null( levels)) {
        if( is.null( max.n)) max.n <- 2*max(dimF)-1
        else {
           if( max.n %%2 == 0) max.n <- max.n-1
           if( max.n > 2*max( dimF)-1) stop(paste("fss2d: max.n must be less than 2N-1, where N is ", max(dimF), sep=""))
        } # end of if else 'max.n' stmts.
        if( max.n < 1) stop("fss2d: max.n must be a positive integer.")
        levels <- seq(1,max.n,2)
   } # end of if no 'levels' given.
   out$levels <- levels
   out$max.n <- max.n
   if( is.null( thresholds)) {
	thresholds <- cbind( quantile( c(Fcst), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95)),
                                quantile(c(Obs), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95)))
	qs <- as.character( c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.90, 0.95))
   } else if( !is.matrix( thresholds)) thresholds <- cbind( thresholds, thresholds)
   out$thresholds <- thresholds
   if( is.null( Pe)) out$Pe <- 1/(levels^2)
   if( length( Pe) == 1) out$Pe <- rep( Pe, length( levels))
   else  out$Pe <- Pe
   out$subset <- subset
   out$loc <- loc
   out$qs <- qs
   out$units <- units
   class( out) <- "hoods2dPrep"
   return( out)
} # end of 'hoods2dPrep' function.

plot.hoods2d <-
function(x, ...) {
   a <- get( x$prep.object)
   mets <- names( x)
   if( "upscaling" %in% mets) {
	try( upscale2dPlot( x$upscaling, args=a))
   } # end of if plot upscaling components.
   if( "mincvr" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$mincvr$ets, args=a, main=paste("Min. Coverage: Gilbert Skill Score (GSS)", sep="")))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$mincvr$far, args=a, main="Min. Coverage: False Alarm Ratio"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$mincvr$pod, args=a, main="Min. Coverage: Hit Rate"))
   } # end of if 'mincvr' stmts.
   if( "multi.event" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$multi.event$pod, args=a, main=paste("Multi-Event Contingency Table: Hit Rate", sep="")))
	a$ylab <- "False Alarm Rate"
        try( hoods2dPlot( x$multi.event$f, args=a, main="Multi-Event Contingency Table: False Alarm Rate"))
	a$ylab <- "Hanssen-Kuipers Score"
	try(hoods2dPlot( x$multi.event$hk, args=a, main="Multi-Event Contingency Table: HK"))
   } # end of if 'multi.event' stmts.
   if( "fuzzy" %in% mets) {
	par( mfrow=c(3,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$fuzzy$ets, args=a, main="Fuzzy Logic: Gilbert Skill Score (GSS)"))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$fuzzy$far, args=a, main="Fuzzy Logic: False Alarm Ratio (FAR)"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$fuzzy$pod, args=a, main="Fuzzy Logic: Hit Rate"))
   } # end of if 'fuzzy' stmts.
   if( "joint" %in% mets) {
	par( mfrow=c(4,2), bg="beige")
	a$ylab <- "Gilbert Skill Score"
        try( hoods2dPlot( x$joint$ets, args=a, main="Joint Probability: Gilbert Skill Score (GSS)"))
	a$ylab <- "False Alarm Ratio"
        try( hoods2dPlot( x$joint$far, args=a, main="Joint Probability: False Alarm Ratio (FAR)"))
	a$ylab <- "Hit Rate"
        try( hoods2dPlot( x$joint$pod, args=a, main="Joint Probability: Hit Rate"))
   } # end of if 'joint' stmts.
   if( "fss" %in% mets) {
	look <- a
	look$values <- x$fss$values
	look$fss.random <- x$fss$fss.random
	look$fss.uniform <- x$fss$fss.uniform
	try( fss2dPlot( look, main="Fractions Skill Score (FSS)"))
   } # end of if 'fss' stmts.
   if( "pragmatic" %in% mets) {
	par( mfrow=c(2,2), bg="beige")
	a$ylab <- "Brier Score"
	try( hoods2dPlot( x$pragmatic$bs, args=a, main="Pragmatic: Brier Score (BS)"))
	a$ylab <- "Brier Skill Score"
	try( hoods2dPlot( x$pragmatic$bss, args=a, main="Pragmatic: Brier Skill Score (BSS)"))
   } # end of if 'pragmatic' stmts.
   invisible()
} # end of 'plot.hoods2d' function.

fss2dPlot <-
function(x, ...) {
   ##
   ## Function to make useful plots for the output from 'fss2d'.
   ##
   ## Arguments:
   ##
   ## 'x' a list object of class "fss2d".
   ## '...' optional arguments to the 'image' and 'image.plot' functions.
   ##
   ## Value: No value is returned, but a two-panel plot is created: the first is an image plot (quilt plot) showing the FSS scores
   ##	for each threshold/level combination.  The second is a line plot with a different colored line for each threshold.  The second
   ##	plot also shows the FSS for a random forecast (dotted line) and based on a uniform distribution (dotted line).
   ##
   # image plot
   odim <- dim( x$values)
   q <- odim[2]
   l <- odim[1]
   if( is.null( odim)) stop("fss2d: values must be a matrix.")
   if( q == 1 & l == 1) stop("fss2d: values must be a matrix with at least one dimension > 1.")
   par( mfrow=c(1,2), mar=c(5.1, 4.1, 3.1, 4.1), bg="beige")
   # image plot
   if( is.null( x$qs) & all( x$thresholds[,1] == x$thresholds[,2])) {
        a1.labels <- as.character( x$thresholds[,1])
        xlb <- paste("Threshold (", x$units, ")", sep="")
   } else if( !is.null( x$qs)) {
        a1.labels <- x$qs
        xlb <- "Threshold (quantile)"
   } else {
        a1.labels <- as.character( 1:dim( x$thresholds)[1])
        xlb <- "Threshold index number"
   }
   image( t( x$values), xaxt="n", yaxt="n", xlab=xlb, ylab="Neighborhood size (grid squares)", col=c("grey", tim.colors(256)), ...)
   text( x=seq(0,1,,q)[ rep(1:q,l)], y=seq(0,1,,l)[ rep(1:l,each=q)], labels=round( t( x$values), digits=2))
   axis( 1, at=seq(0,1,,q), labels=a1.labels)
   if( any( x$thresholds[,1] != x$thresholds[,2])) warning("fss2dPlot: thresholds differ for the two fields, labels only for the first")
   axis( 2, at=seq(0,1,,length(x$levels)), labels=x$levels)
   # lines( seq( 0, 1,,q), x$fss.random, lty=1, col="darkorange")
   # lines( seq( 0, 1,,q), x$fss.uniform, lty=2, col="darkred")
   image.plot( x$values, legend.only=TRUE, col=c("grey", tim.colors(256)), ...)

   # line plot
   look <- x$values
   matplot( look, ylim = c(0,1), ylab = "FSS", xlab = "Neighborhood size (grid squares)", type = "l", lty = 1, axes = FALSE , lwd = 2)
   abline(h=c(x$fss.uniform), col=1:q, lty=2)
   abline(h=c(x$fss.random), col=1:q, lty=3)
   axis(2)
   box()
   axis(1, at = 1:l, lab = x$levels)
   grid()
   legend("topleft", legend = a1.labels, col = 1:q, title = xlb, inset = 0.02, lwd = 2 )
   invisible()
} # end of 'fss2dPlot' function.

Fourier2d <-
function(x, xdim=NULL) {
   ##
   ## Function to compute the Fast Fourier Transform (FFT) of a field 'x'.
   ##
   ## Arguments:
   ##
   ## 'x' 'n X m' numeric matrix.
   ## 'xdim' (optional) numeric vector of length 2.  If NULL, it will be computed.
   ##
   ## Details: This function creates a matrix of zeros with twice the dimension of the input field, 'x'.
   ##	'x' is then put in the upper left corner of the matrix of zeros, missing values are set to zero,
   ##	and the FFT is calculated for the large matrix.  This is what is returned.  The output can then be
   ##	used by the 'kernel2d' function in order to possibly save an FFT in computing time at some point.
   ##
   ## See also: 'fft', 'kernel2d', 'hoods2d'
   ##
   ## Value: '2n X 2m' numeric matrix giving the Fourier transformed field.
   ##
   if( is.null( xdim)) xdim <- dim( x)
   out <- matrix( 0, 2*xdim[1], 2*xdim[2])
   out[1:xdim[1],1:xdim[2]] <- x
   out[ is.na( out)] <- 0
   return( fft( out))
} # end of 'Fourier2d' function.

vxstats <-
function(Fcst, Obs, which.stats=c("bias", "ts", "ets", "pod", "far", "f", "hk", "mse"), subset=NULL) {
   ##
   ## Function to calculate various traditional verification statistics for a gridded
   ## verification set.
   ##
   ## Arguments:
   ##
   ## 'Fcst', 'Obs' 'k X m' logical or numeric matrices of forecast and observed values, resp.
   ## 'which.stats' character vector telling which verification statistics should be computed.
   ## 'subset' numeric vector indicating a subset of points over which to calculate the statistics.  If NULL, then the entire
   ##	fields are used.
   ##
   ## Details:
   ##
   ## The possible statistics that can be computed, as determined by 'which.stats' are:
   ##
   ##	"bias" the number of forecast events divided by the number of observed events (sometimes called frequency bias).
   ##	"ts" threat score, given by hits/(hits + misses + false alarms)
   ##	"ets" equitable threat score, given by (hits - hits.random)/(hits + misses + false alarms - hits.random), where
   ##		'hits.random' is the number of observed events times the number of forecast events divided by the total
   ##		number of forecasts.
   ##	"pod" probability of detecting an observed event (aka, hit rate).  It is given by hits/(hits + misses).
   ##	"far" false alarm ratio, given by (false alarms)/(hits + false alarms).
   ##	"f" false alarm rate (aka probability of false detection) is given by (false alarms)/(correct rejections + false alarms).
   ##	"hk" Hansen-Kuipers Score is given by the difference between the hit rate ("pod") and the false alarm rate ("f").
   ##	"mse" mean square error (not a contingency table statistic, but can be used with binary fields).  This is the only
   ##		statistic that can be calculated here that does not require binary fields for 'Fcst' and 'Obs'.
   ##
   ## Warnings: It is up to the user to provide the appropriate type of fields for the given statistics
   ##	to be computed.  For example, they must be binary for all types of 'which.stats' except "mse".
   ##
   ## Value: list object with component names the same as 'which.stats' giving the single numeric
   ##	value of each statistic for the two fields.
   ##
   out <- list()
   xdim <- dim( Fcst)
   if( "mse" %in% which.stats) {
        if( is.null( subset)) {
         Nxy <- sum( colSums( !is.na( Fcst) & !is.na( Obs), na.rm=TRUE), na.rm=TRUE)
         out$mse <- sum( colSums( (Fcst - Obs)^2, na.rm=TRUE), na.rm=TRUE)/Nxy
        } else {
           out$mse <- mean( (c(Fcst)[ subset] - c(Obs)[subset])^2, na.rm=TRUE)
        }
   } # end of if do MSE.
   if( any( c("bias", "ts", "ets", "pod", "far", "f", "hk") %in% which.stats)) {
	if( !is.logical( Fcst)) Fcst <- as.logical( Fcst)
	if( !is.logical( Obs)) Obs <- as.logical( Obs)
	if( is.null( subset)) {
	   hits <- sum( colSums( matrix( Fcst & Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   miss <- sum( colSums( matrix( !Fcst & Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   fa   <- sum( colSums( matrix( Fcst & !Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	   if( any( c("ets", "f", "hk") %in% which.stats)) cn <- sum( colSums( matrix( !Fcst & !Obs, xdim[1], xdim[2]), na.rm=TRUE), na.rm=TRUE)
	} else {
	   hits <- sum( c(Fcst)[subset] & c(Obs)[subset], na.rm=TRUE)
           miss <- sum( !c(Fcst)[subset] & c(Obs)[subset], na.rm=TRUE)
           fa   <- sum( c(Fcst)[subset] & !c(Obs)[subset], na.rm=TRUE)
           if( any( c("ets", "f", "hk") %in% which.stats)) cn <- sum( !c(Fcst)[subset] & !c(Obs)[subset], na.rm=TRUE)
	}
	if( "bias" %in% which.stats) {
	   if( (hits + fa == 0) & (hits + miss == 0)) out$bias <- 1
	   else if( hits + miss == 0) out$bias <- (hits + fa)/(1e-8)
	   else out$bias <- (hits + fa)/(hits + miss)
 	}
	if( "ts" %in% which.stats) {
	   if( hits == 0) out$ts <- 0
	   else out$ts <- hits/(hits + miss + fa)
	}
	if( "ets" %in% which.stats) {
	   if( (hits + miss == 0) | (hits + fa == 0)) hits.random <- 0
	   else hits.random <- (hits + miss)*(hits + fa)/(hits + miss + fa + cn)
	   if( hits + miss + fa == 0) out$ets <- 0
	   else out$ets <- (hits - hits.random)/(hits + miss + fa - hits.random)
	}
	if( any( c("pod", "hk") %in% which.stats)) {
	   if( hits + miss == 0) pod <- 0
	   else pod <- hits/(hits + miss)
	   if( "pod" %in% which.stats) out$pod <- pod
	} 
	if( "far" %in% which.stats) {
	   if( hits + fa == 0) out$far <- 0
	   else out$far <- fa/(hits + fa)
	}
	if( any( c("f", "hk") %in% which.stats)) {
	   if( cn + fa == 0) f <- 0
	   else f <- fa/(cn + fa)
	   if( "f" %in% which.stats) out$f <- f
	   if( "hk" %in% which.stats) out$hk <- pod - f
	}
    } # end of if any contingency table scores stmts.
   class( out) <- "vxstats"
   return( out)
} # end of 'vxstats' function.

hoods2dSetUpLists <-
function( which.methods, mat) {
   out <- list()
   if( "upscaling" %in% which.methods) {
      out$upscaling <- list()
      out$upscaling$rmse <- numeric( nrow( mat)) 
      out$upscaling$bias <- out$upscaling$ts <- out$upscaling$ets <- mat
   } # end of if set up 'upscaling' component.
   if( "mincvr" %in% which.methods) {
      out$mincvr <- list()
      out$mincvr$pod <- out$mincvr$far <- out$mincvr$ets <- mat
   } # end of if calculate minimum coverage stmts.
   if( "multi.event" %in% which.methods) {
      out$multi.event <- list()
      out$multi.event$pod <- out$multi.event$f <- out$multi.event$hk <- mat
   } # end of if calculate "multi.event" stmts.
   if( "fuzzy" %in% which.methods) {
      out$fuzzy <- list()
      out$fuzzy$pod <- out$fuzzy$far <- out$fuzzy$ets <- mat
   } # end of if do fuzzy logic method.
   if( "joint" %in% which.methods) {
      out$joint <- list()
      out$joint$pod <- out$joint$far <- out$joint$ets <- mat
   } # end of if do joint prob method.
   if( "fss" %in% which.methods) {
      out$fss <- list()
      out$fss$values <- mat
      out$fss$fss.random <- out$fss$fss.uniform <- numeric(0)
   } # end of if do FSS method.
   if( "pragmatic" %in% which.methods) {
      out$pragmatic <- list()
      out$pragmatic$bs <- out$pragmatic$bss <- mat
   } # end of do pragmatic method.
   return( out)
} # end of 'hoods2dSetUpLists' function.

MinCvg2dfun <-
function( sIy, sIx, subset=NULL) {
   ##
   ## Function to calculate the minimum coverage neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy', 'sIx' (optional) 'k X m' binary forecast and observed matrices, resp.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##	all of the points are used.
   ##
   ## Value: list object with components: "pod", "far" and "ets".
   ##
   out <- vxstats( sIy, sIx, which.stats=c("pod", "far", "ets"), subset=subset)
   return( out)
} # end of 'MinCvg2dfun' function.

multicon2dfun <-
function(sIy, Ix, subset=NULL) {
   ##
   ## Function to calculate the multi-event contingency table neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy' 'k X m' binary smoothed forecast matrix.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components "pod" and "f" and "hk".
   ##
   out <- vxstats( sIy, Ix, which.stats=c("pod", "f", "hk"), subset=subset)
   return( out)
} # end of 'multicon2dfun' function.

fuzzyjoint2dfun <-
function( sPy, sPx, subset=NULL) {
   ##
   ## Function to calculate the fuzzy logic and joint probability neighborhood methods.
   ##
   ## Arguments:
   ##
   ## 'sPy', 'sPx' smoothed 'k X m' forecast and observed matrices, resp., output from 'kernel2d'.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components: "fuzzy" and "joint", each of which are themselves
   ##	list objects each with components: "pod", "far" and "ets".
   ##
   out <- list()
   vxfun <- function(n11, n01, n10, n00) {
	pod <- n11/(n11 + n01)
	far <- n10/(n11 + n10)
	hits.random <- (n11 + n01)*(n11 + n10)/(n11 + n01 + n10 + n00)
	ets <- (n11 - hits.random)/(n11 + n01 + n10 - hits.random)
	return( list(pod=pod, far=far, ets=ets))
   } # end of internal 'vxfun' function.
   if( is.null( subset)) {
      hits <- sum( colSums( sPx < sPy, na.rm=TRUE), na.rm=TRUE)
      miss <- sum( colSums( sPx < (1-sPy), na.rm=TRUE), na.rm=TRUE)
      fa   <- sum( colSums( (1-sPx) < sPy, na.rm=TRUE), na.rm=TRUE)
      cn   <- sum( colSums( (1-sPx) < (1-sPy), na.rm=TRUE), na.rm=TRUE)
      out$fuzzy <- vxfun( hits, miss, fa, cn)
      hits <- sum( colSums( sPx*sPy, na.rm=TRUE), na.rm=TRUE)
      miss <- sum( colSums( sPx*(1-sPy), na.rm=TRUE), na.rm=TRUE)
      fa   <- sum( colSums( (1-sPx)*sPy, na.rm=TRUE), na.rm=TRUE)
      cn   <- sum( colSums( (1-sPx)*(1-sPy), na.rm=TRUE), na.rm=TRUE)
      out$joint <- vxfun( hits, miss, fa, cn)
   } else {
      hits <- sum( c(sPx)[subset] < c(sPy)[subset], na.rm=TRUE)
      miss <- sum( c( sPx)[subset] < (1-c(sPy)[subset]), na.rm=TRUE)
      fa   <- sum( (1-c(sPx)[subset]) < c(sPy)[subset], na.rm=TRUE)
      cn   <- sum( (1-c(sPx)[subset]) < (1-c(sPy)[subset]), na.rm=TRUE)
      out$fuzzy <- vxfun( hits, miss, fa, cn)
      hits <- sum( (c(sPx)[subset])*(c(sPy)[subset]), na.rm=TRUE)
      miss <- sum( (c( sPx)[subset])*(1-c(sPy)[subset]), na.rm=TRUE) 
      fa   <- sum( (1-c(sPx)[subset])*c(sPy)[subset], na.rm=TRUE)
      cn   <- sum( (1-c(sPx)[subset])*(1-c(sPy)[subset]), na.rm=TRUE)
      out$joint <- vxfun( hits, miss, fa, cn)
   } # end of if else 'subset' stmt.
   return( out)
} # end of 'fuzzyjoint2dfun' function.

pragmatic2dfun <-
function( sPy, Ix, mIx=NULL, subset=NULL) {
   ##
   ## Function to calculate the pragmatic neighborhood statistics.
   ##
   ## Argmunets:
   ##
   ## 'sPy' 'k X m' matrix of smoothed forecast event frequencies.
   ## 'Ix' 'k X m' binary matrix of observed events.
   ## 'mIx' single numeric giving the base rate.  If NULL, it will be computed here.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components 'bs' and 'bss' giving the Brier and Brier Skill Scores.
   ##
   bs <- vxstats( sPy, Ix, which.stats="mse", subset=subset)$mse
   if( is.null( subset)) {
      if( is.null( mIx)) mIx <- mean( Ix, na.rm=TRUE)
      denom <- sum( colSums( (mIx - Ix)^2, na.rm=TRUE), na.rm=TRUE)/sum( !is.na( Ix), na.rm=TRUE)
   } else {
      if( is.null( mIx)) mIx <- mean( c(Ix)[subset], na.rm=TRUE)
      denom <- mean( (mIx - c(Ix)[subset])^2, na.rm=TRUE)
   }
   bss <- 1 - bs/denom
   return( list( bs=bs, bss=bss))
} # end of 'pragmatic2dfun' function.

upscale2dfun <-
function( sYy, sYx, threshold=NULL, subset=NULL) {
   ##
   ## Function to calculate the upscaling neighborhood statistics.
   ##
   ## Arguments:
   ##
   ## 'sIy', 'sIx' 'k X m' upscaled forecast and observed matrices (e.g., as output from 'kernel2d'), resp.
   ## 'threshold' (optional) numeric of length 2 giving the value over which to calculate the statistics.  If NULL,
   ##	only RMSE is calculated.  The forecast threshold is first, and the observbed one second.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: list object with components "rmse", and if 'threshold' is not NULL, "bias", "gss", and "ts".
   ##
   out <- list()
   out$rmse <- sqrt( vxstats( sYy, sYx, which.stats="mse", subset=subset)$mse)
   if( !is.null( threshold)) {
	xdim <- dim( sYy)
	binmat <- matrix(0, xdim[1], xdim[2]) 
	sIx <- sIy <- binmat
	sIx[ sYx >= threshold[ 2]] <- 1
	sIy[ sIy >= threshold[ 1]] <- 1
	tmp <- vxstats( sIy, sIx, which.stats=c("bias", "ts", "ets"), subset=subset)
	out$bias <- tmp$bias
	out$ts <- tmp$ts
	out$ets <- tmp$ets
   } # end of if 'threshold' stmts.
   return( out)
} # end of 'upscale2dfun'.

hoods2dPlot <-
function(x, args, ...) {
   require( fields)
   odim <- dim( x)
   q <- odim[2]
   l <- odim[1]
   if( is.null( odim)) stop("hoods2dPlot: values must be a matrix.")
   if( q == 1 & l == 1) stop("hoods2dPlot: values must be a matrix with at least one dimension > 1.")
   if( is.null( args$qs) & all( args$thresholds[,1] == args$thresholds[,2])) {
        a1.labels <- as.character( args$thresholds[,1])
        xlb <- paste("Threshold (", args$units, ")", sep="")
   } else if( !is.null( args$qs)) {
        a1.labels <- args$qs
        xlb <- "Threshold (quantile)"
   } else {
        a1.labels <- as.character( 1:dim( args$thresholds)[1])
        xlb <- "Threshold index number"
   }

   image( t( x), xaxt="n", yaxt="n", xlab=xlb, ylab="Neighborhood size (grid squares)", col=c("grey", heat.colors(12)), ...)
   text( x=seq(0,1,,q)[ rep(1:q,l)], y=seq(0,1,,l)[ rep(1:l,each=q)], labels=round( t( x), digits=2))
   axis( 1, at=seq(0,1,,q), labels=a1.labels)
   if( any( args$thresholds[,1] != args$thresholds[,2])) warning("hoods2dPlot: thresholds differ for the two fields, labels only for the first")
   axis( 2, at=seq(0,1,,length(args$levels)), labels=args$levels)
   image.plot( t( x), legend.only=TRUE, col=c("grey", heat.colors(12)), ...)

   matplot( x, ylab=args$ylab, xlab="Neighborhood size (grid squares)", type = "l", lty = 1, axes = FALSE , lwd = 2)
   axis(2)
   box()
   axis(1, at = 1:l, lab = args$levels)
   grid()
   legend("topleft", legend = a1.labels, col = 1:q, title = xlb, inset = 0.02, lwd = 2 )
   invisible()
} # end of 'hoods2dPlot' function.

upscale2dPlot <-
function(object, args, ...) {
   par( mfrow=c(4,2), bg="beige")
   try( hoods2dPlot( object$ets, args=args, main="Upscaling: Gilbert Skill Score (GSS)", ...))
   try( hoods2dPlot( object$ts, args=args, main="Upscaling: Threat Score (TS)", ...))
   try( hoods2dPlot( object$bias, args=args, main="Upscaling: Bias", ...))
   try( plot( args$levels, object$rmse, type="b", xlab="Neighborhood Length (grid squares)", ylab="RMSE", main="Upscaling: RMSE", col="darkblue"))
} # end of 'upscale2dPlot' function.

fss2dfun <- function(sPy, sPx, subset=NULL, verbose=FALSE) {
   ## 
   ## Function to calculate the FSS neighborhood statistics.  Not to be confused
   ## with 'fss2d', which does the same thing, but this function begins with the smoothed
   ## fields, and only calculates the final score.  To be used internally by 'fss2d'.
   ## 
   ## Arguments:
   ## 
   ## 'sPy', 'sPx' smoothed 'k X m' forecast and observed matrices, resp., output from 'kernel2d'.
   ## 'subset' numeric indicating over which points the summary scores should be calculated.  If NULL,
   ##   all of the points are used.
   ##
   ## Value: single numeric giving the FSS value.
   ##
   id1 <- !is.na( sPy)
   id2 <- !is.na( sPx)
   if( verbose) cat("Finding the numbers of non-missing values for each field.\n")
   if( is.null( subset)) {
      N1 <- sum( colSums( id1, na.rm=TRUE), na.rm=TRUE)
      N2 <- sum( colSums( id2, na.rm=TRUE), na.rm=TRUE)
      Nxy <- sum( colSums( id1 & id2, na.rm=TRUE), na.rm=TRUE)
      if( verbose) cat(Nxy, " total number of non-missing points.\n")
      num <- sum( colSums( (sPy - sPx)^2, na.rm=TRUE), na.rm=TRUE)/Nxy
      if( verbose) cat("MSE = ", num, "\n")
      denom <- sum( colSums( sPy^2, na.rm=TRUE), na.rm=TRUE)/N1 + sum( colSums( sPx^2, na.rm=TRUE), na.rm=TRUE)/N2
   } else {
      num <- mean( (c(sPy)[subset]-c(sPx)[subset])^2, na.rm=TRUE)
      denom <- mean( (c(sPy)[subset])^2+(c(sPx)[subset])^2, na.rm=TRUE)
   }
   if( verbose) cat("Reference MSE = ", denom, "\n")
   return(1-num/denom)
} # end of 'fss2dfun' function.

pphindcast2d <- function(obj, which.score="ets", verbose=FALSE, ...) {
   Fcst <- get( obj$Fcst.name)
   Obs <- get( obj$Obs.name)
   xdim <- obj$xdim
   Nxy <- obj$Nxy
   subset <- obj$subset
   thresholds <- obj$thresholds
   levels <- obj$levels
   q <- dim( thresholds)[1]
   l <- length( levels)
   outmat <- Pthresh <- matrix( NA, l, q)
   out <- list()
   out$which.score <- which.score
   findthresh <- function( p, Ix, sPx, binmat, which.score, subset=NULL) {
	sIx <- binmat
	sIx[ sPx >= p] <- 1
	return( -vxstats( sIx, Ix, which.stats=which.score, subset=subset)[[which.score]])
   } # end of 'findthresh' internal function.
   binmat <- matrix(0, xdim[1], xdim[2])
   for( threshold in 1:q) {
	Ix <- Iy <- binmat
	Ix[ Obs >= thresholds[threshold,2]] <- 1
	Iy[ Fcst >= thresholds[threshold,1]] <- 1
	IxFFT <- Fourier2d( Ix, xdim=xdim)
	IyFFT <- Fourier2d( Iy, xdim=xdim)
	for( level in 1:l) {
	   Wlvl <- kernel2d( Ix, levels[ level], setup=TRUE)
	   sPy <- kernel2d( Iy, levels[level], W=Wlvl, X=IyFFT, xdim=xdim, Nxy=Nxy)
	   sPx <- kernel2d( Ix, levels[level], W=Wlvl, X=IxFFT, xdim=xdim, Nxy=Nxy)
   	   Pthresh[ level, threshold] <- optim( 0, findthresh, Ix=Ix, sPx=sPx, binmat=binmat, which.score=which.score, subset=subset, 
							lower=0, upper=1, method="L-BFGS-B", ...)$par
	   if( verbose) cat("Pthresh = ", Pthresh[level,threshold], " for obs threshold = ", thresholds[ threshold,2], " and level = ", levels[level], "\n")
	   sIy <- binmat
	   sIy[ sPy >= Pthresh[ level, threshold]] <- 1
	   outmat[ level, threshold] <- vxstats( sIy, Ix, which.stats=which.score)[[which.score]]
	} # end of for 'level' loop.
   } # end of for 'threshold' loop.
   out$values <- outmat
   out$Pthresh <- Pthresh
   return( out)
} # end of 'pphindcast' function.

deltamm <- function(x, p = 2, max.delta = Inf, const = Inf, type = c( "sqcen", "original" ),
    N = NULL, verbose = FALSE, ...) {

    if (verbose) begin.tiid <- Sys.time()

    if (is.null(x$X.feats) && is.null(x$Y.feats)) stop("deltamm: no features to merge/match!")
    if (is.null(x$X.feats)) stop("deltamm: no verification features present.")
    if (is.null(x$Y.feats)) stop("deltamm: no model features present.")

    out <- x

    out$match.type <- "deltamm"
    out$match.message <- "Objects merged/matched based on the Baddeley Delta metric via the deltamm function."

    a <- attributes( x )

    type <- match.arg( type )

    if( type == "original" ) res <- deltammOrig( x = x, p = p, max.delta = max.delta, const = const, verbose = verbose, ... )
    else if( type == "sqcen" ) res <- deltammSqCen( x = x, p = p, max.delta = max.delta, const = const, N = N, verbose = verbose, ... )

    out$Q <- res$Q
    out$unmatched <- res$unmatched
    out$matches   <- res$matches
    out$merges    <- res$merges

    out$MergeForced <- TRUE

    if( verbose ) print( Sys.time() - begin.tiid )

    class( out ) <- "matched"

    return(out)

} # end of 'deltamm' function.

deltammSqCen <- function(x, p = 2, max.delta = Inf, const = Inf, N = NULL, verbose = FALSE, ...) {

    a <- attributes( x )

    out <- list()

    X <- x$X.feats
    Xhat <- x$Y.feats

    xdim <- a$xdim

    if( is.null( N ) ) {

	N <- max( xdim )
	if( N %% 2 == 0 ) N <- N + 1

    } # end of if 'N' not given stmt.

    bdfun <- function( id, OB, FC, p, const, N, verbose, ... ) {

	j <- id[ 1 ]
	k <- id[ 2 ]

	if( verbose ) cat("Calculating Baddeley Delta Metric between forecast feature ", k, " and observed feature ", j, "\n")

	# if( class( OB[[ 1 ]] ) == "im" )
	A <- as.matrix( OB[[ j ]] )
	# else A <- OB[[ k ]][[ j ]]

	# if( class( FC[[1]] ) == "im" )
	B <- as.matrix( FC[[k]] )
        # else B <- FC[[ j ]][[ k ]]

	res <- censqdelta( x = A, y = B, N = N, p = p, const = const, ... )

    } # end of internal 'bdfun' function.

    n <- length( X )
    m <- length( Xhat )

    no.X <- is.null( X ) || n == 0
    no.Xhat <- is.null( Xhat ) || m == 0

    if( no.X || no.Xhat ) {

        if( no.X && no.Xhat && verbose ) cat("No identified features in either field.  Therefore, no matches/merges made.\n")
        else if( no.X && verbose ) cat("No identified observed features.  Therefore, no matches/merges made.\n")
        else if( no.Xhat && verbose ) cat("No identified model features.  Therefore, no matches/merges made.\n")

        funmatched <- 1:m
        vxunmatched <- 1:n

        matches <- cbind(integer(0), integer(0))
        merges <- NULL
        out$unmatched <- list(X = vxunmatched, Xhat = funmatched)
        out$matches <- matches
        out$merges <- merges

        class(out) <- "matched"

        return(out)

    } # end of if no features in one or both fields stmts.

    if( verbose ) cat( "Step 1: Finding Upsilon matrix containing the Baddeley delta between each individual feature across fields.\n" )
    ind <- cbind( rep( 1:n, m ), rep( 1:m, each = n ) )

    Upsilon <- apply( ind, 1, bdfun, OB = X, FC = Xhat, p = p, const = const, N = N, verbose = verbose, ... )
    Upsilon <- matrix( Upsilon, n, m )

    if( verbose ) {

        cat( "Step 1 completed.  Upsilon matrix given by:\n" )
        print( Upsilon )

    }

    Psi <- Ksi <- matrix( NA, n, m )
    o.Ksi <- t( apply( Upsilon, 1, order, na.last = TRUE ) )
    o.Psi <- apply( Upsilon, 2, order, na.last = TRUE )
    Ksi[, 1] <- apply( Upsilon, 1, min, na.rm = TRUE )
    Psi[1, ] <- apply( Upsilon, 2, min, na.rm = TRUE )

    if( verbose ) cat( " Finding ", m - 1, " potential forecast merges for each of the ", n, " observed features.\n" )

    for( j in 1:n ) {

        o <- ( 1:m )[ o.Ksi[ j, ] ]

        if( m >= 2 ) {

	    newobj <- list()
	    newobj[[ 1 ]] <- union.owin( Xhat[[ o[ 1 ] ]], Xhat[[ o[ 2 ] ]] )

        } else newobj <- Xhat

        if( m >= 3 ) for( i in 3:m ) newobj[[ i - 1 ]] <- union.owin( newobj[[ i - 2 ]], Xhat[[ o[ i ] ]] )

    } # end of for 'j' loop.

    if( m > 1 ) {

        if( verbose ) cat( "Calculating delta metrics for forecast merges.\n" )

        ind <- cbind( rep( 1:n, m - 1 ), rep( 1:( m - 1 ), each = n ) )
        look <- paste( ind[, 1 ], ind[, 2 ], sep = "-" )
        ind <- ind[ !duplicated( look ), ] 

        tmp <- apply( ind, 1, bdfun, OB = X, FC = newobj, p = p, const = const, N = N, verbose = verbose, ... )
        Ksi[, 2:m ] <- tmp

    } else Ksi <- Upsilon

    if( verbose ) cat( "\nFinding ", n - 1, " potential observed merges for each of the ", m, " forecast features.\n" )

    for( k in 1:m ) {

        if( verbose ) cat( k, " " )

        o <- ( 1:n )[ o.Psi[, k ] ]

        if( n >= 2 ) {

	    newobj <- list()
	    newobj[[ 1 ]] <- union.owin( X[[ o[ 1 ] ]], X[[ o[ 2 ] ]] )
	
        } else newobj <- X

        if( n >= 3 ) for( i in 3:n ) newobj[[ i - 1 ]] <- union.owin( newobj[[ i - 2 ]], X[[ o[ i ] ]] )

    } # end of for 'k' loop.

    if (n > 1) {

        if( verbose ) cat( "\nAll metrics for forecast merges found.  Calculating delta metrics for observed merges.\n" )

        ind <- cbind( rep( 1:( n - 1 ), m ), rep( 1:m, each = n - 1 ) )
        look <- paste( ind[, 1], ind[, 2], sep = "-" )
        ind <- ind[ !duplicated( look ), ]

        tmp <- apply( ind, 1, bdfun, OB = newobj, FC = Xhat, p = p, const = const, N = N, verbose = verbose, ... )

        Psi <- t(Psi)
        Psi[, 2:n ] <- tmp
        Psi <- t( Psi )

    } else Psi <- Upsilon

    bigQ <- array( NA, dim = c( n, m, 3 ) )
    bigQ[, , 1] <- Upsilon
    bigQ[, , 2] <- Psi
    bigQ[, , 3] <- Ksi

    out$Q <- bigQ

    if (all(bigQ > max.delta)) {

        if( verbose ) cat( "\nAll delta metrics are larger than max.delta, no merges/matches.\n" )

        funmatched <- 1:m
        vxunmatched <- 1:n

        matches <- cbind( integer(0), integer(0) )

        merges <- NULL
        out$unmatched <- list( X = vxunmatched, Xhat = funmatched )
        out$matches <- matches
        out$merges <- merges

        return(out)

    } else {

	if( verbose ) {

            cat("Psi:\n")
            print(Psi)
            cat("Ksi:\n")
            print(Ksi)
            cat("\nAll Baddeley metrics found.  Book keeping ...\n")

        } 

        J <- K <- list()

        ind <- cbind( rep( 1:n, 3 * m ), rep( rep( 1:m, each = n ), 3 ), rep( 1:3, each = n * m ), 1:( 3 * n * m ) )

        for( jk in 1:( 3 * n * m ) ) {

            if( jk <= n * m ) {

                J[[ jk ]] <- ind[ jk, 1 ]
                K[[ jk ]] <- ind[ jk, 2 ]

            } else if( ( jk > n * m ) && ( jk <= 2 * n * m ) ) {

                K[[ jk ]] <- ind[ jk, 2 ]
                jj <- 1:( ind[ jk, 1 ] )
                J[[ jk ]] <- ( 1:n )[ o.Psi[ jj, ind[ jk, 2 ] ]]

            } else {

                J[[ jk ]] <- ind[ jk, 1 ]
                kk <- 1:( ind[ jk, 2 ] )
                K[[ jk ]] <- ( 1:m )[ o.Ksi[ ind[ jk, 1 ], kk ] ]

            }

        } # end of for 'jk' loop.

        iter <- 1
        matches <- cbind( integer(0), integer(0) )
        nn <- 1:n
        mm <- 1:m

        bigQ[ bigQ > max.delta ] <- NA
        bigQ <- c( bigQ )
        efun <- function( x, ftr ) return( any( is.element( ftr, x ) ) )

        while( length( nn ) > 0 && length( mm ) > 0 && iter < 3 * n * m + 1 ) {

            if( any( !is.na( bigQ ) ) ) {

                minQ <- min( bigQ, na.rm = TRUE )

                id <- ( ind[, 4] )[ bigQ == minQ ]
                id <- id[ !is.na( id ) ]
                id <- id[ 1 ]

                vx <- J[[ id ]]
                fc <- K[[ id ]]

                newmatch <- cbind( fc, vx )
                matches <- rbind( matches, newmatch )

                bigQ[ id ] <- NA

                id1 <- unlist( lapply( J, efun, ftr = vx ) )
                id2 <- unlist( lapply( K, efun, ftr = fc ) )
                bigQ[ id1 | id2 ] <- NA

                nn <- nn[ !is.element( nn, vx ) ]
                mm <- mm[ !is.element( mm, fc ) ]
                iter <- iter + 1

            } else iter <- Inf

        } # end of while loop.

        out$unmatched <- list( X = nn, Xhat = mm )
        colnames( matches ) <- c( "Forecast", "Observed" )
        out$matches <- matches
        merges <- MergeIdentifier( matches )

        out$merges <- merges

    } # end of if else all Q values larger than maximum allowed delta stmts.

    return( out )

} # end of 'deltammSqCen' function.

deltammOrig <- function( x, p = 2, max.delta = Inf, const = Inf, verbose = FALSE, ... ) {

    begin.tiid <- Sys.time()

    a <- attributes(x)
    out <- list()

    bfun <- function(id, OB, FC, p, const, verbose) {
        j <- id[1]
        k <- id[2]

	if( verbose ) cat("Calculating Baddeley Delta Metric between forecast feature ", k, " and observed feature ", j, "\n")

        if( class( OB[[ 1 ]] ) == "im" ) A <- OB[[ j ]]
        else A <- OB[[ k ]][[ j ]]

        if( class( FC[[ 1 ]] ) == "im") B <- FC[[k]] 
        else B <- FC[[ j ]][[ k ]]

        if (!is.infinite(const)) {

            A <- eval.im(pmin.int(A, const))
            B <- eval.im(pmin.int(B, const))

        }

        if (is.infinite(p)) {

            res <- eval.im(abs(A - B))
            res <- summary(res)$max

        } else {

            res <- eval.im(abs(A - B)^p)
            res <- summary(res)$mean
            res <- res^(1/p)

        }

        return(res)

    } # end of internal 'bfun' function.

    X <- x$X.feats
    Xhat <- x$Y.feats

    # xdim <- dim(X[[1]][["m"]])

    xdim <- a$xdim

    n <- length(X)
    m <- length(Xhat)
    no.X <- is.null(X) || n == 0
    no.Xhat <- is.null(Xhat) || m == 0

    if( no.X || no.Xhat ) {

        if( no.X && no.Xhat && verbose ) cat( "No identified features in either field.  Therefore, no matches/merges made.\n" )
        else if( no.X && verbose ) cat( "No identified observed features.  Therefore, no matches/merges made.\n" )
        else if( no.Xhat && verbose ) cat( "No identified model features.  Therefore, no matches/merges made.\n" )

        funmatched <- 1:m
        vxunmatched <- 1:n
        matches <- cbind(integer(0), integer(0))
        merges <- NULL

        out$unmatched <- list(X = vxunmatched, Xhat = funmatched)
        out$matches <- matches
        out$merges <- merges

        class(out) <- "matched"

        return(out)

    } # end of if no features in one or both fields stmts.

    if( verbose ) cat( "Finding a bounding box for the entire set of features.\n" )
    Xbb <- do.call( boundingbox, X )
    Ybb <- do.call( boundingbox, Xhat )

    bb <- boundingbox( Xbb, Ybb )

    if( verbose ) {

        cat( "Finished finding overall bounding box for the entire set of features.\n" )
        print( bb )

    }

    Xbox <- lapply( X, rebound, rect = bb )
    Ybox <- lapply( Xhat, rebound, rect = bb )

    if( verbose ) cat( "Finding the distance maps for each feature in each field.\n" )

    dX <- lapply( Xbox, distmap, ... )
    dXhat <- lapply( Ybox, distmap, ... )

    if( verbose ) cat( "Step 1: Finding Upsilon matrix containing the Baddeley delta between each individual feature across fields.\n" )
    ind <- cbind( rep( 1:n, m ), rep( 1:m, each = n ) )

    Upsilon <- apply( ind, 1, bfun, OB = dX, FC = dXhat, p = p, const = const, verbose = verbose )
    Upsilon <- matrix( Upsilon, n, m )

    if (verbose) {
        cat("Step 1 completed.  Upsilon matrix given by:\n")
        print(Upsilon)
    }
    Psi <- Ksi <- matrix(NA, n, m)
    o.Ksi <- t(apply(Upsilon, 1, order, na.last = TRUE))
    o.Psi <- apply(Upsilon, 2, order, na.last = TRUE)
    Ksi[, 1] <- apply(Upsilon, 1, min, na.rm = TRUE)
    Psi[1, ] <- apply(Upsilon, 2, min, na.rm = TRUE)
    if (verbose) 
        cat("Finding distance maps for potential merges.\n")
    if (verbose) 
        cat(" Finding ", m - 1, " potential forecast merges for each of the ", 
            n, " observed features.\n")
    Ksi.dmap <- Psi.dmap <- list()
    for (j in 1:n) {
        if (verbose) 
            cat(j, " ")
        o <- (1:m)[o.Ksi[j, ]]
        newobj <- list()
        if (m >= 2) 
            newobj[[1]] <- union.owin(Xhat[[o[1]]], Xhat[[o[2]]])
        else newobj <- Xhat
        if (m >= 3) 
            for (i in 3:m) newobj[[i - 1]] <- union.owin(newobj[[i - 
                2]], Xhat[[o[i]]])
        newobj2 <- lapply(newobj, rebound, rect = bb)
        Ksi.dmap[[j]] <- lapply(newobj2, distmap, ...)
    }
    if (verbose) 
        cat("\nFinding ", n - 1, " potential observed merges for each of the ", 
            m, " forecast features.\n")
    for (k in 1:m) {
        if (verbose) 
            cat(k, " ")
        o <- (1:n)[o.Psi[, k]]
        newobj <- list()
        if (n >= 2) 
            newobj[[1]] <- union.owin(X[[o[1]]], X[[o[2]]])
        else newobj <- X
        if (n >= 3) 
            for (i in 3:n) newobj[[i - 1]] <- union.owin(newobj[[i - 
                2]], X[[o[i]]])
        newobj2 <- lapply(newobj, rebound, rect = bb)
        Psi.dmap[[k]] <- lapply(newobj2, distmap, ...)
    }
    if (verbose) 
        cat("\nAll necessary distance maps have been found.\n")
    if (m > 1) {
        if (verbose) 
            cat("Calculating delta metrics for forecast merges.\n")
        ind <- cbind(rep(1:n, m - 1), rep(1:(m - 1), each = n))
        look <- paste(ind[, 1], ind[, 2], sep = "-")
        ind <- ind[!duplicated(look), ]
        tmp <- apply(ind, 1, bfun, OB = dX, FC = Ksi.dmap, p = p, 
            const = const, verbose = verbose)
        Ksi[, 2:m] <- tmp
    }
    else Ksi <- Upsilon
    if (n > 1) {

        if( verbose ) cat( "\nAll metrics for forecast merges found.  Calculating delta metrics for observed merges.\n" )
        ind <- cbind(rep(1:(n - 1), m), rep(1:m, each = n - 1))
        look <- paste(ind[, 1], ind[, 2], sep = "-")
        ind <- ind[!duplicated(look), ]

        tmp <- apply( ind, 1, bfun, OB = Psi.dmap, FC = dXhat, p = p, const = const, verbose = verbose )

        Psi <- t( Psi )
        Psi[, 2:n ] <- tmp
        Psi <- t( Psi )

    } else Psi <- Upsilon

    bigQ <- array(NA, dim = c(n, m, 3))
    bigQ[, , 1] <- Upsilon
    bigQ[, , 2] <- Psi
    bigQ[, , 3] <- Ksi
    out$Q <- bigQ

    if( all( bigQ > max.delta ) ) {

        if( verbose ) cat( "\nAll delta metrics are larger than max.delta, no merges/matches.\n" )

        funmatched <- 1:m
        vxunmatched <- 1:n
        matches <- cbind(integer(0), integer(0))
        merges <- NULL
        out$unmatched <- list(X = vxunmatched, Xhat = funmatched)
        out$matches <- matches
        out$merges <- merges
        return(out)

    } else {

        if (verbose) {
            cat("Psi:\n")
            print(Psi)
            cat("Ksi:\n")
            print(Ksi)
            cat("\nAll Baddeley metrics found.  Book keeping ...\n")
        }
        J <- K <- list()
        ind <- cbind(rep(1:n, 3 * m), rep(rep(1:m, each = n), 
            3), rep(1:3, each = n * m), 1:(3 * n * m))
        for (jk in 1:(3 * n * m)) {
            if (jk <= n * m) {
                J[[jk]] <- ind[jk, 1]
                K[[jk]] <- ind[jk, 2]
            }
            else if ((jk > n * m) && (jk <= 2 * n * m)) {
                K[[jk]] <- ind[jk, 2]
                jj <- 1:(ind[jk, 1])
                J[[jk]] <- (1:n)[o.Psi[jj, ind[jk, 2]]]
            }
            else {
                J[[jk]] <- ind[jk, 1]
                kk <- 1:(ind[jk, 2])
                K[[jk]] <- (1:m)[o.Ksi[ind[jk, 1], kk]]
            }
        }
        iter <- 1
        matches <- cbind(integer(0), integer(0))
        nn <- 1:n
        mm <- 1:m
        bigQ[bigQ > max.delta] <- NA
        bigQ <- c(bigQ)
        efun <- function(x, ftr) return(any(is.element(ftr, x)))
        while (length(nn) > 0 && length(mm) > 0 && iter < 3 * 
            n * m + 1) {
            if (any(!is.na(bigQ))) {
                minQ <- min(bigQ, na.rm = TRUE)
                id <- (ind[, 4])[bigQ == minQ]
                id <- id[!is.na(id)]
                id <- id[1]
                vx <- J[[id]]
                fc <- K[[id]]
                newmatch <- cbind(fc, vx)
                matches <- rbind(matches, newmatch)
                bigQ[id] <- NA
                id1 <- unlist(lapply(J, efun, ftr = vx))
                id2 <- unlist(lapply(K, efun, ftr = fc))
                bigQ[id1 | id2] <- NA
                nn <- nn[!is.element(nn, vx)]
                mm <- mm[!is.element(mm, fc)]
                iter <- iter + 1
            }
            else iter <- Inf
        }
        out$unmatched <- list(X = nn, Xhat = mm)
        colnames(matches) <- c("Forecast", "Observed")
        out$matches <- matches
        merges <- MergeIdentifier(matches)

        out$merges <- merges

    }

    if (verbose) print(Sys.time() - begin.tiid)

    # class(out) <- "matched"

    return(out)

} # end of 'deltammOrig' function.

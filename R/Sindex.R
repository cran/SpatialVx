Sindex <- function(x, thresh=NULL, ...) {

    UseMethod("Sindex", x)

} # end of 'Sindex' function.

Sindex.SpatialVx <- function(x, thresh=NULL, ..., time.point = 1, obs = 1, model = 1) {

    a <- attributes(x)

    ## Begin: Get the data sets
    dat <- datagrabber(x, time.point = time.point, obs = obs, model = model)


    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Sindex.default(x=X, thresh=thresh, ..., loc=a$loc)
    res2 <- Sindex.default(x=Xhat, thresh=thresh, ..., loc=a$loc)

    res <- rbind(res1, res2)

    rownames(res) <- c(a$obs.name[ obs ], a$model.name[ model ] )

    return(res)

} # end of 'Sindex.SpatialVx' function.

Sindex.default <- function(x, thresh=NULL, ..., loc=NULL) {

    if(!is.null(thresh)) x[x < thresh] <- 0
    n <- sum(colSums(x>0,na.rm=TRUE),na.rm=TRUE)
    n2 <- sqrt(n)

    if( floor(n2) == n2 ) Pmin <- 4 * n2
    else Pmin <- 2 * ( floor( 2 * n2 ) + 1 )

    if(is.null(loc)) {

	xdim <- dim(x)
	loc <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    }

    id <- c(x) != 0
    corners <- apply(loc[id,], 2, range, finite = TRUE)
    P <- 2 * ( diff( corners[,1] ) + 1 ) + 2 * ( diff( corners[,2] ) + 1 )

    res <- c(Pmin / P, Pmin, P)
    names(res) <- c("Sindex", "Pmin", "P")

    return(res)

} # end of 'Sindex.default' function.
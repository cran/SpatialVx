Cindex <- function(x, thresh=NULL, connect.method="C", ...) {

    UseMethod("Cindex", x)

} # end of 'Cindex' function.

Cindex.SpatialVx <- function(x, thresh=NULL, connect.method="C", ..., time.point=1, obs = 1, model=1) {

    a <- attributes(x)

    ## Begin: Get the data sets
    dat <- datagrabber(x, time.point = time.point, obs = obs, model = model)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Cindex.default(x=X, thresh=thresh, connect.method=connect.method, ...)
    res2 <- Cindex.default(x=Xhat, thresh=thresh, connect.method=connect.method, ...)

    res <- c(res1, res2)

    names(res) <- c(a$obs.name[ obs ], a$model.name[ model ] )

    return(res)

} # end of 'Cindex.SpatialVx' function.

Cindex.default <- function(x, thresh=NULL, connect.method="C", ...) {

   if(!is.null(thresh)) x[x < thresh] <- 0

   NP <- sum( colSums( x > 0, na.rm = TRUE ), na.rm = TRUE )

   x[x==0] <- NA
   x <- as.im(x)
   x <- connected(x, method = connect.method)

   NC <- max( as.numeric( as.matrix( x ) ), na.rm = TRUE )

   return( 1 - ( NC - 1 ) / ( sqrt( NP ) + NC ) )

} # end of 'Cindex.default' function.
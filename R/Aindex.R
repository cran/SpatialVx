



Aindex <- function(x, thresh=NULL, dx=1, dy=1, ...) {

    UseMethod("Aindex", x)

} # end of 'Aindex' function.

Aindex.SpatialVx <- function(x, thresh=NULL, dx=1, dy=1, ..., time.point = 1, obs = 1, model = 1) {

    a <- attributes(x)

    ## Begin: Get the data sets
    dat <- datagrabber(x, time.point = time.point, obs = obs, model = model)
   
    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Aindex.default(x=X, thresh=thresh, dx=dx, dy=dy, ...)
    res2 <- Aindex.default(x=Xhat, thresh=thresh, dx=dx, dy=dy, ...)

    res <- rbind(res1, res2)

    rownames(res) <- c(a$obs.name[ obs ], a$model.name[ model ])

    return(res)
} # end of 'Aindex.SpatialVx' function.

Aindex.default <- function(x, thresh=NULL, dx=1, dy=1, ...) {
    if(is.null(thresh)) thresh <- 1e-8 
    x <- as.im(x)
    x <- solutionset(x>=thresh)
    ch <- convexhull(x)
    A <- sum( colSums( as.matrix( x ), na.rm = TRUE ), na.rm = TRUE ) * dx * dy
    Aconvex <- area.owin(ch)*dx*dy
    Aindex <- A/Aconvex

    res <- c(Aindex, A, Aconvex, dx, dy)
    names(res) <- c("Aindex", "A", "Aconvex", "dx", "dy")

    return(res)
} # end of 'Aindex.default' function.

griddedVgram <- function(object, zero.in=TRUE, zero.out=TRUE, ...) {
   out <- object
   X <- get(object$Vx.name)
   Y <- get(object$Fcst.name)
   out$zero.in <- zero.in
   out$zero.out <- zero.out
   out$Vx.vgram.matrix <- list()
   out$Fcst.vgram.matrix <- list()
   if(zero.in) out$Vx.vgram.matrix[[1]] <- vgram.matrix(dat=X, ...)
   else out$Vx.vgram.matrix[[1]] <- NULL
   if(zero.in) out$Fcst.vgram.matrix[[1]] <- vgram.matrix(dat=Y, ...)
   else out$Fcst.vgram.matrix[[1]] <- NULL
   if(zero.out) out$Vx.vgram.matrix[[2]] <- variogram.matrix(dat=X, zero.out=TRUE, ...)
   else out$Vx.vgram.matrix[[2]] <- NULL
   if(zero.out) out$Fcst.vgram.matrix[[2]] <- variogram.matrix(dat=Y, zero.out=TRUE, ...)
   else out$Fcst.vgram.matrix[[2]] <- NULL
   class(out) <- c(class(object), "griddedVgram")
   return(out)
} # end of 'griddedVgram' function.

plot.griddedVgram <- function(x, ...) {
   mainX <- x$Vx.name
   mainY <- x$Fcst.name
   if(x$zero.in) {
      vgX <- x$Vx.vgram.matrix[[1]]
      vgY <- x$Fcst.vgram.matrix[[1]]
   }
   if(x$zero.out) {
      vgX.zero <- x$Vx.vgram.matrix[[2]]
      vgY.zero <- x$Fcst.vgram.matrix[[2]]
   }
   if(((x$zero.in) & !(x$zero.out)) | (!(x$zero.in) & (x$zero.out))) par(mfrow=c(2,2), mar=rep(4.1,4), bg="beige")
   if(x$zero.in & x$zero.out) par(mfrow=c(4,2), mar=rep(4.1,4), bg="beige")
   if(x$zero.in) {
      plot(vgX$d, vgX$vgram, xlab="separation distance", ylab="variogram", main=mainX, col="darkblue")
      points(vgX$d.full, vgX$vgram.full, pch=".", cex=1.25, col="darkblue")
      plot(vgY$d, vgY$vgram, xlab="separation distance", ylab="variogram", main=mainY, col="darkblue") 
      points(vgY$d.full, vgY$vgram.full, pch=".", cex=1.25, col="darkblue")
      class(vgX) <- "vgram.matrix"
      class(vgY) <- "vgram.matrix"
      zl <- range(c(c(vgX$vgram.full), c(vgY$vgram.full)), finite=TRUE)
      plot(vgX, xlab="x separations", ylab="y separations", zlim=zl, ...)
      plot(vgY, xlab="x separations", ylab="y separations", zlim=zl, ...)
   }
   if(x$zero.out) {
      plot(vgX.zero$d, vgX.zero$vgram, xlab="separation distance", ylab="variogram", main=paste(mainX, " (non-zeros only)", sep=""), col="darkblue")
      points(vgX.zero$d.full, vgX.zero$vgram.full, pch=".", cex=1.25, col="darkblue")
      plot(vgY.zero$d, vgY.zero$vgram, xlab="separation distance", ylab="variogram", main=paste(mainY, " (non-zeros only)", sep=""), col="darkblue")
      points(vgY.zero$d.full, vgY.zero$vgram.full, pch=".", cex=1.25, col="darkblue")
      class(vgX.zero) <- "vgram.matrix"
      class(vgY.zero) <- "vgram.matrix"
      zl <- range(c(c(vgX.zero$vgram.full), c(vgY.zero$vgram.full)), finite=TRUE)
      plot(vgX.zero, xlab="x separations", ylab="y separations", zlim=zl, ...)
      plot(vgY.zero, xlab="x separations", ylab="y separations", zlim=zl, ...)
   }
   invisible()
} # end of 'plot.griddedVgram' function.

corrskill <- function(x,y,...) {
   good <- !is.na(x) & !is.na(y)
   n <- sum(good, na.rm=TRUE)
   s1 <- sd(c(x), ...)
   s2 <- sd(c(y), ...)
   m1 <- sum(colSums(x,na.rm=TRUE),na.rm=TRUE)/n
   m2 <- sum(colSums(y,na.rm=TRUE),na.rm=TRUE)/n
   out <- ((n/(n-1))/(s1*s2))*(x-m1)*(y-m2)
   return(out)
} # end of 'corrskill' function.

sqerrloss <- function(x,y,...) return((x-y)^2)
abserrloss <- function(x,y,...) return(abs(x-y))
distmaploss <- function(x,y, threshold=0, const=Inf, ...) {
   if(length(threshold)==1) threshold <- rep(threshold,2)
   xdim <- dim(x)
   if(is.finite(const)) {
	x <- cbind(matrix(0, xdim[1], const), x, matrix(0, xdim[1], const))
	x <- rbind(matrix(0, const, xdim[2]+2*const), x, matrix(0, const, xdim[2] + 2*const))
	y <- cbind(matrix(0, xdim[1], const), y, matrix(0, xdim[1], const))
        y <- rbind(matrix(0, const, xdim[2]+2*const), y, matrix(0, const, xdim[2] + 2*const))
   }
   x <- im(x)
   y <- im(y)
   x <- solutionset(x >= threshold[2])
   y <- solutionset(y >= threshold[1])
   bb <- bounding.box(as.rectangle(x), as.rectangle(y))
   x <- rebound(x, bb)
   y <- rebound(y, bb)
   dy <- distmap(y, ...)$v
   dx <- distmap(x, ...)$v
   if(is.finite(const)) {
	dy[dy>const] <- const
	dx[dx>const] <- const
	dy <- dy[(const+1):(xdim[1]+const),(const+1):(xdim[2]+const)]
	dx <- dx[(const+1):(xdim[1]+const),(const+1):(xdim[2]+const)]
   }
   return(abs(dx - dy))
} # end of 'distmaploss' function.


spatMLD <- function(x,y1,y2,lossfun="corrskill", trend="ols", loc=NULL, maxrad=20, dx=1, dy=1, zero.out=FALSE, ...) {
   out <- list()
   out$Vx.name <- as.character(substitute(x))
   out$Mod1.name <- as.character(substitute(y1))
   out$Mod2.name <- as.character(substitute(y2))
   out$trend <- trend
   xdim <- dim(x)
   out$lossfun <- lossfun
   out$lossfun.args <- list(...)
   out$vgram.args <- list(maxrad=maxrad,dx=dx,dy=dy)
   g1 <- do.call(lossfun, args=c(list(x=x,y=y1),list(...)))
   g2 <- do.call(lossfun, args=c(list(x=x,y=y2),list(...)))
   d <- matrix(g1-g2, xdim[1], xdim[2])
   if(zero.out) {
	# zeros <- d == 0
	zeros <- (x==0) & (y1==0) & (y2==0)
	beta <- mean(!zeros, na.rm=TRUE)
	out$zeros <- zeros
	out$beta <- beta
   }
   if(trend=="ols") {
	if(is.null(loc)) loc <- cbind(rep(1:xdim[1],xdim[2]), rep(1:xdim[2],each=xdim[1]))
	else out$loc <- as.character(substitute(loc))
	dat <- data.frame(y=c(d), x1=loc[,1], x2=loc[,2])
	fit <- lm(y~x1+x2, data=dat)
	tr <- matrix(predict(fit),xdim[1],xdim[2])
	if(zero.out) tr[zeros] <- 0
	d <- d - tr
	out$trend <- fit
   } else if(is.numeric(trend)) {
	if(zeros) warning("spatMLD: zero.out is TRUE, but trend provided.  Not setting original zeros back after removing trend.")
	d <- d - trend
	out$trend <- trend
   }
   out$d <- d
   if(!zero.out) vg <- vgram.matrix(dat=d, R=maxrad, dx=dx, dy=dy)
   else vg <- variogram.matrix(dat=d, R=maxrad, dx=dx, dy=dy, zero.out=zero.out)
   out$lossdiff.vgram <- vg
   # tmp <- data.frame(h=vg$d, v=vg$vgram)
   # Exp.vgram <- function(h, sigma2=1, theta=1) return(sigma2*(1-exp(-h/theta)))
   # vgmodel <- nls(v~Exp.vgram(h=h,sigma2=s^2,theta=r), data=tmp, start=list(s=sqrt(vg$vgram[1]), r=maxrad))
   # out$vgmodel <- vgmodel
   class(out) <- "spatMLD"
   return(out)
} # end of 'spatMLD' function.

fit.spatMLD <- function(object) {
   out <- object
   maxrad <- out$vgram.args$maxrad
   vgdat <- data.frame(h=out$lossdiff.vgram$d, v=out$lossdiff.vgram$vgram)
   v1 <- sqrt(vgdat$v[1])
   Exp.vgram <- function(h, sigma2=1, theta=1) return(sigma2*(1-exp(-h/theta)))
   vgmodel <- nls(v~Exp.vgram(h=h,sigma2=s^2,theta=r), data=vgdat, start=list(s=v1, r=maxrad))
   out$vgmodel <- vgmodel
   return(out)
} # end of 'fit.spatMLD' function.

plot.spatMLD <- function(x, ...) {
   msg <- paste(x$lossfun, ": ", x$Mod1.name, " vs ", x$Mod2.name, " (", x$Vx.name, ")", sep="")
   par(mfrow=c(2,2), mar=rep(4.1,4), bg="beige")
   if(is.null(x$zeros)) image.plot(x$d, main=msg, axes=FALSE)
   else {
	Im <- x$d
	Im[x$zeros] <- NA
	image.plot(Im, main=msg, axes=FALSE)
   }
   hist(x$d, breaks="FD", xlab="Mean Loss Differential", col="darkblue", freq=FALSE, main=msg)
   a <- x$lossdiff.vgram
   plot(a$d, a$vgram, col="darkblue", xlab="separation distance", ylab="variogram")
   if(!is.null(x$vgmodel)) {
	# tmp <- coef(x$vgmodel)
   	# sig2 <- tmp[1]^2
   	# r <- tmp[2]
   	# b <- sig2*(1 - exp(-a$d/r))
	# lines(a$d, b, col="darkorange", lwd=1.5)
	lines(a$d, predict(x$vgmodel), col="darkorange", lwd=1.5)
   	legend("bottomright", legend=c("Empirical", "Model"), pch=c("o", ""), col=c("darkblue","darkorange"), lty=c(0,1), lwd=1.5, bty="n")
   }
   plot.vgram.matrix(a,main="variogram by direction")
   invisible()
} # end of 'plot.spatMLD' function.

summary.spatMLD <- function(object, ...) {
   out <- object
   msg <- paste(object$lossfun, ": ", object$Mod1.name, " vs ", object$Mod2.name, " (against verification: ", object$Vx.name, ")", sep="")
   print(msg)
   d <- object$d
   if(zero.out <- !is.null(object$beta)) {
	cat("\n", "Estimate of beta present, so calculating Dbar and test statistic over non-zero entries only.\n")
	cat("\n", "Frequency of non-zero loss differential (beta) is: ", object$beta, "\n")
	d[d==0] <- NA
   }
   good <- !is.na(d)
   n <- length(d[good])
   cat("\n", "number of non-zero loss differential points is: ", n, "\n")
   cat("\n", "Mean Loss Differential: \n")
   Dbar <- sum(colSums(d,na.rm=TRUE),na.rm=TRUE)/n
   out$Dbar <- Dbar
   cat("\n", Dbar, "\n", "\n")
   a <- object$lossdiff.vgram
   cat("\n", "Summary of empirical variogram values:\n")
   print(stats(a))
   if(!is.null(object$vgmodel)) {
      # co <- coef(object$vgmodel)
      # sig2 <- co[1]^2
      # r <- co[2]
      # b <- sig2*exp(-a$d.full/r)
      # if(zero.out) denom <- sqrt(mean(object$beta^2*b,na.rm=TRUE))
      # else 
      b <- predict(object$vgmodel)
      denom <- sqrt(mean(b,na.rm=TRUE))
      SV <- Dbar/denom
      out$test.statistic <- SV
      cat("Test Statistic for null hypothesis of equal predictive ability on average\n")
      print(SV)
      pval <- pnorm(abs(SV), lower.tail=FALSE)
      cat("p-value for (two-sided) test is: ", pval, "\n")
      out$p.value <- pval
   }
   invisible(out)
} # end of 'summary.spatMLD' function.

variogram.matrix <- function (dat, R = 5, dx = 1, dy = 1, zero.out=FALSE) 
{
    SI <- function(ntemp, delta) {
        n1 <- 1:ntemp
        n2 <- n1 + delta
        good <- (n2 >= 1) & (n2 <= ntemp)
        cbind(n1[good], n2[good])
    }
    if(zero.out) dat[dat==0] <- NA
    N <- ncol(dat)
    M <- nrow(dat)
    m <- min(c(round(R/dx), M))
    n <- min(c(round(R/dy), N))
    ind <- rbind(as.matrix(expand.grid(0, 1:n)), as.matrix(expand.grid(1:m, 
        0)), as.matrix(expand.grid(c(-(m:1), 1:m), 1:n)))
    d <- sqrt((dx * ind[, 1])^2 + (dy * ind[, 2])^2)
    good <- (d > 0) & (d <= R)
    ind <- ind[good, ]
    d <- d[good]
    ind <- ind[order(d), ]
    d <- sort(d)
    nbin <- nrow(ind)
    holdVG <- rep(NA, nbin)
    holdN <- rep(NA, nbin)
    for (k in 1:nbin) {
        MM <- SI(M, ind[k, 1])
        NN <- SI(N, ind[k, 2])
	numNA <- sum(is.na(dat[MM[,1],NN[,1]]) | is.na(dat[MM[,2],NN[,2]]),na.rm=TRUE)
        holdN[k] <- length(MM) * length(NN) - numNA
        BigDiff <- (dat[MM[, 1], NN[, 1]] - dat[MM[, 2], NN[,2]])
        holdVG[k] <- mean(0.5 * (BigDiff)^2, na.rm=TRUE)
    }
    top <- tapply(holdVG * holdN, d, FUN = "sum")
    bottom <- tapply(holdN, d, FUN = "sum")
    dcollapsed <- as.numeric(names(bottom))
    vgram <- top/bottom
    dimnames(vgram) <- NULL
    out <- list(vgram = vgram, d = dcollapsed, ind = ind, d.full = d, 
        vgram.full = holdVG, N = holdN, dx = dx, dy = dy)
    class(out) <- "vgram.matrix"
    return(out)
} # end of 'variogram.matrix' function.

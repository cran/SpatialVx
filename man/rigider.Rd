\name{rigider}
\alias{rigider}
\alias{plot.rigided}
\alias{print.rigided}
\alias{summary.rigided}
\alias{rigidTransform}

\title{

Rigid Transformation

}

\description{

Find the optimal rigid transformation for a spatial field (e.g. an image).

}

\usage{

rigider(x1, x0, p0, init = c(0, 0, 0), type = c("regular", "fast"),
    translate = TRUE, rotate = FALSE, loss, loss.args = NULL,
    interp = "bicubic", stages = TRUE,
    verbose = FALSE, ...)

\method{plot}{rigided}(x, ...)

\method{print}{rigided}(x, ...)

\method{summary}{rigided}(object, ...)

rigidTransform(theta, p0, N, cen)

}

\arguments{

  \item{x1, x0}{
    matrices of same dimensions giving the forecast (or 1-energy) and observation (or 0-energy) fields, resp.
}
  \item{x, object}{list object of class \dQuote{rigided} as output by \code{rigider}.}
  \item{N}{(optional) the dimension of the fields (i.e., if \code{x1} and \code{x0} are n by m, then \code{N} is the product m * n).}
  \item{cen}{N by 2 matrix whoes rows are all the same giving the center of the field (used to subtract before determining rotations, etc.).}
  \item{p0}{
    N by 2 matrix giving the coordinates for the 0-energy (observed) field.
}
  \item{init}{
    (optional) numeric vector of length equal to the number of parameters (e.g., 2 for translation only, 3 for both, and 1 for rotation only).  If missing, then these will be estimated by taking the difference in centroids (translation) and the difference in orientation angles (rotation) as determined using image moments by way of \code{imomenter}.
}
  \item{theta}{numeric vector of length 1, 2 or 3 (depending on whether you want to translate only (2), rotate only (1) or both (3)) giving the rigid transformation parameters.}
  \item{type}{character stating whether to optimize a loss function or just find the centroid (and possibly orientation angle) difference(s).}
  \item{translate, rotate}{
logical, should the optimal translation/rotation be found?
}
  \item{loss}{
character naming a loss function (see details) to use in optimizing the rigid transformation (defaults to square error loss.
}
  \item{loss.args}{
named list giving any optional arguments to \code{loss}.
}
  \item{interp}{
character naming the 2-d interpolation method to use in calls to \code{Fint2d}.  Must be one of \dQuote{round} (default), \dQuote{bilinear} or \dQuote{bicubic}.
}
  \item{stages}{
logical.  Should the optimal translation be found before finding both the optimal tranlsation and rotation?
}
  \item{verbose}{
logical.  Should progress information be printed to the screen?
}
  \item{\dots}{
optional arguments to \code{nlminb}.
}

}

\details{

A rigid transformation translates coordinates of values in a matrix and/or rotates them.  That is, if (r, s) are coordinates in a field with center (c1, c2), then the rigid transformation with parameters (x, y) and theta is given by:

(r, s) + (x, y) + Phi ((r, s) - (c1, c2)),

where Phi is the matrix with first column given by (cos( theta ), - sin( theta)) and second column given by (sin( theta), cos( theta )).

The optimal transformation is found by way of numerical optimization using the \code{nlminb} function on the loss function given by \code{loss}.  If no value is given for \code{loss}, then square error loss is assumed.  In this case, the loss function is based on an assumption of Gaussian errors, but this assumption is only important if you try to make inferences based on this model, in which case you should probably think much harder about what you are doing.  In particular, the default objective function, Q, is given by:

Q = - sum( ( F(W(s)) - O(s) )^2 / (2 * sigma^2) - (N / 2) * log( sigma^2 ),

where s are the coordinates, W(s) are the rigidly transformed coordinates, F(W(s)) is the value of the 1-enegy field (forecast) evaluated at W(s) (which is interpolated as the translations typically do not give integer translations), O(s) is the 0-energy (observed) field evaluated at coordinate s, and sigma^2 is the estimated variance of the error field.  A good alternative is to use \dQuote{QcorrRigid}, which calculates the correlation between F and O instead, and has been found by some to give better performance.

The function \code{rigidTransform} performs a rigid transform for given parameter values.  It is intended as an internal function, but may be of use to some users.

}

\value{
A list object of class \dQuote{rigided} is returned with components:
  \item{call }{the function call.}
  \item{translation.only }{If stages argument is true, this part is the optimal translation before rotation.}
  \item{rotate}{ optimal translation and rotation together, if stages argument is true.}
  \item{initial}{initial values used.}
  \item{interp.method}{ same as input argument interp.}
  \item{optim.args}{optional arguments passed to nlminb.}
  \item{loss, loss.args}{same as input arguments.}
  \item{par}{optimal parameter values found.}
  \item{value}{value of loss function at optimal parameters.}
  \item{x0, x1, p0}{same as input arguments.}
  \item{p1}{transformed p0 coordinates.}
  \item{x1.transformed}{The field F(W(s)).}
}
\author{

Eric Gilleland

}

\note{

Finding the optimal rigid transformation can be very tricky when applying both rotatons and translations.  This function helps, but for some fields may require more user input than is ideal, and should be considered experimental for the time being; as the examples will demonstrate.  It does seem to work well for translations only, which has been the recommended course of action for the CRA method.

}

\seealso{

\code{\link{nlminb}}, \code{\link{Fint2d}}

}

\examples{

# Simple uninteresting example for the R robots.
x <- y <- matrix(0, 20, 40)

x[ 12:18, 2:3 ] <- 1

y[ 13:19, 5:6 ] <- 1

xycoords <- cbind(rep(1:20, 40), rep(1:40, each = 20))

tmp <- rigider(x1 = x, x0 = y, p0 = xycoords)
tmp
plot(tmp)

# Rotate a coordinate system.
data( "geom000" )

loc <- cbind(rep(1:601, 501), rep(1:501, each = 601))

# Rotate the coordinates by pi / 4.
th <- c(0, 0, pi / 4)
names(th) <- c("x", "y", "rotation")
cen <- colMeans(loc[ geom000 > 0, ])
loc2 <- rigidTransform(theta = th, p0 = loc, cen = cen)

geom101 <- Fint2d(X = geom000, Ws = loc2, s = loc, method = "round")

\dontrun{

image.plot(geom101)

# Try to find the optimal rigid transformation.
# First, allow a translation as well as rotation.

tmp <- rigider(x1 = geom101, x0 = geom000, p0 = loc,
    rotate = TRUE, verbose = TRUE)
tmp
plot(tmp)

# Now, only allow rotation, which does not work as
# well as one would hope.
tmp <- rigider(x1 = geom101, x0 = geom000, p0 = loc,
    translate = FALSE, rotate = TRUE, verbose = TRUE)
tmp
plot(tmp)

# Using correlation.
tmp <- rigider(x1 = geom101, x0 = geom000, p0 = loc,
    rotate = TRUE, loss = "QcorrRigid", verbose = TRUE)
tmp
summary(tmp)
plot(tmp)

##
## Examples from ICP phase 1.
##
## Geometric cases.
##

data( "geom001" )
data( "geom002" )
data( "geom003" )
data( "geom004" )
data( "geom005" )


tmp <- rigider(x1 = geom001, x0 = geom000, p0 = loc, verbose = TRUE)
tmp
plot(tmp)

tmp <- rigider(x1 = geom002, x0 = geom000, p0 = loc, verbose = TRUE)
tmp
plot(tmp)

tmp <- rigider(x1 = geom003, x0 = geom000, p0 = loc, verbose = TRUE)
tmp
plot(tmp)

tmp <- rigider(x1 = geom004, x0 = geom000, p0 = loc, verbose = TRUE)
tmp
plot(tmp)

# Note: Above is a scale error rather than a rotation, but can we
# approximate it with a rotation?
tmp <- rigider(x1 = geom004, x0 = geom000, p0 = loc, rotate = TRUE,
    verbose = TRUE)
tmp
plot(tmp)
# Ok, maybe need to give it better starting values?  Or, run it again
# with just the translation.

tmp <- rigider(x1 = geom005, x0 = geom000, p0 = loc, verbose = TRUE)
tmp
plot(tmp)


}
}
\keyword{ manip }

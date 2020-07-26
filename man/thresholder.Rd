\name{thresholder}
\alias{thresholder}
\alias{thresholder.default}
\alias{thresholder.SpatialVx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a Threshold to a Field
}
\description{
App;y a threshold to a field and return either a binary field or a field with replace.width everywhere the rule is not true.
}
\usage{
thresholder(x, type = c("binary", "replace.below"), th, rule = ">=",
    replace.with = 0, ...)

\method{thresholder}{default}(x, type = c("binary", "replace.below"), th, rule = ">=", 
    replace.with = 0, ... )

\method{thresholder}{SpatialVx}(x, type = c("binary", "replace.below"), th, rule = ">=",
    replace.with = 0, ..., time.point = 1, obs = 1, model = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A field or \dQuote{SpatialVx} object to which to apply the thresholds.
}
  \item{type}{
character describing which type of field(s) to return: binary or replace.
}
  \item{rule}{
If \code{type} is \dQuote{binary}, return 0 when the rule applied to a grid point's value is not true in relation to the threshold value, and 1 elsewhere.  If \code{type} is \dQuote{replace.below}, then return \code{replace.with} wherever the rule is not true and return the original value otherwise.  By default, it replaces values below the threshold with zero (hence its name), but if \code{rule} is, e.g., \dQuote{<=}, then it will replace values above with zero; or whatever value is chosen for \code{replace.with}.
}
  \item{replace.with}{
Only used if \code{type} is \dQuote{replace.below}.  The value with which to replace values that are below (default) the threshold.
}
  \item{th}{Value of the threshold (default) or index to which row of threshold matrices in \code{thresholds} attribute of \dQuote{SpatialVx} object.  Must be a single number.  }
  \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{\dots}{ Not used.
}
}
\details{
At each point, p, in the field, the expression: p \code{rule} \code{threshold} is applied.  If \code{type} is \dQuote{binary}, then if the expression is false, zero is returned for that grid point, and if it is true, then one is returned.  If \code{type} is \dQuote{replace.below}, then if the expression is false, \code{replace.with} is returned for that grid point, and if true, then the original value is returned.  By default, the original field is returned, but with values below the threshold set to zero.  If \code{rule} is \dQuote{<=}, then \code{replace.below} will actually replace values above the threshold with \dQuote{replace.with} instead.

If applied to a \dQuote{SpatialVx} class object, then observation \code{obs} and model \code{model} at time point \code{time.point} will each be thresholded using the respective \code{th} threshold value for the observed and modeled fields as taken from the thresholds attribute of the object (see the help file for \code{make.SpatialVx}).
}
\value{
A field of the same dimension as \code{x} if a matrix.  If \code{x} is a \dQuote{SpatialVx} class object, then a list is returned with components:
  \item{X, Xhat }{The matrices giving the respective thresholded fields for the observation and forecast.}
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{make.SpatialVx}}
}
\examples{

x <- matrix( 12 + rnorm( 100, 10, 10 ), 10, 10 )

par( mfrow = c(2, 2) )
image.plot( thresholder( x, th = 12 ), main = "binary" )

image.plot( thresholder( x, type = "replace.below", th = 12 ),
    main = "replace.below" ) 

image.plot( thresholder( x, th = 12, rule = "<=" ),
    main = "binary with rule <=" )

image.plot( thresholder( x, type = "replace.below", th = 12, rule = "<=" ),
    main = "replace.below with rule <=" )

par( mfrow = c(1,1) )
\dontrun{
data("geom000")
data("geom004")
data("ICPg240Locs")

hold <- make.SpatialVx( geom000, geom004, thresholds = c(0.01, 50.01),
    projection = TRUE, map = TRUE, loc = ICPg240Locs, loc.byrow = TRUE,
    field.type = "Geometric Objects Pretending to be Precipitation",
    units = "mm/h", data.name = "ICP Geometric Cases", obs.name = "geom000",
    model.name = "geom004" )

# Note: th = 1 means threshold = 0.01.
look <- thresholder( hold, th = 1 )

image.plot( look$X )
contour( look$Xhat, add = TRUE, col = "white" )

# Note: th = 2, means threshold = 50.01
look <- thresholder( hold, th = 2 )

image.plot( look$X )
contour( look$Xhat, add = TRUE, col = "white" )

look <- thresholder( hold, th = 1, rule = "<" )

image.plot( look$X )
contour( look$Xhat, add = TRUE, col = "white" )
}
}
\keyword{ manip }

\name{pphindcast2d}
\alias{pphindcast2d}
\title{
Function to perform the practically perfect hindcast neighborhood verification method
}
\description{
Function to perform the practically perfect hindcast neighborhood verification method.  Finds the
optimal threhsold, Pthresh, and calculates the desired statistic for that threshold.
}
\usage{
pphindcast2d(obj, which.score = "ets", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A list object returned by the \code{hoods2dPrep} function.
}
  \item{which.score}{
character stating which verification score is to be used.  Must be one that is accepted by \code{vxstats}.
}
  \item{verbose}{
logical, should progress information be printed ot the screen?
}
  \item{\dots}{
optional arguments to the \code{optim} function.  May not include lower, upper or method as these are hard
coded into the function.
}
}
\details{
The practically perfect hindcast method is described in Ebert (2008).  Using a similar notation as that described therein (and in the help
page for \code{hoods2d}), the method is a SO-NF approach that first compares the observed binary field (obtained from the trheshold(s)
provided by obj from \code{hoods2dPrep}), Ix, with the smoothed binary field, <Px>s.  This smoothed binary field is thresholded by
Pthresh to obtain a new binary field.  The value of Pthresh that maximizes the verification score (provided by the which.score argument)
is then used to compare Ix with <Iy>s, the binary forecast field obtained by thresholding the smoothed binary forecast field Iy using
the value of Pthresh found above.  The verification statistic determined by which.score is calculated between Ix and <Iy>s.
}
\value{
 A list object is returned with components:
    \item{which.score }{value of which.score, same as the argument passed in.}
    \item{Pthresh }{l X q matrix giving the value of Pthresh applied at each level (rows) and threshold (columns).}
    \item{values}{l X q matrix giving the value of which.score found for each level (rows) and threshold (columns).}
}
\references{
Ebert EE, 2008. Fuzzy verification of high resolution gridded
     forecasts: A review and proposed framework.  \emph{Meteorol. Appl.},
     \bold{15}, 51-64. DOI: 10.1002/met.25 Available at
     http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf
}
\author{
Eric Gilleland
}

\section{Warning }{The value Pthresh is optimized under the assumption that larger values of which.score are better.} 

\seealso{
\code{\link{hoods2d}}, \code{\link{kernel2dsmooth}}, \code{\link{vxstats}}, \code{\link{hoods2dPlot}}, \code{\link{optim}}
}
\examples{
x <- y <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y[ sample(1:50,20), sample(1:50,20)] <- rexp( 400)
hold <- hoods2dPrep("y", "x", thresholds=c(0.1, 0.5), levels=c(1, 3))
look <- pphindcast2d( hold)
look
\dontrun{
data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01), levels=c(1, 3, 65), loc=ICPg240Locs, units="in/100")
look <- pphindcast2d( hold, verbose=TRUE)
par( mfrow=c(1,2))
hoods2dPlot( look$values, args=list( thresholds=hold$thresholds, units="in/100", levels=hold$levels), main="Gilbert Skill Score")
}
}
\keyword{ math }

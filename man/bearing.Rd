\name{bearing}
\alias{bearing}
\title{
Bearing from One Spatial Location to Another
}
\description{
Find the bearing from one spatial location to another.
}
\usage{
bearing(point1, point2, deg = TRUE, aty = "compass")
}
\arguments{
  \item{point1, point2}{ two-column numeric matrices giving lon/lat coordinates for the origin point(s) (\code{point1}) and the destination point(s) (\code{point2}).}
  \item{deg}{
logical, should the output be converted from radians to degrees?
}
  \item{aty}{
character stating either \dQuote{compass} (default) or \dQuote{radial}.  The former gives the standard compass bearing angle (0 is north, increase clockwise), and the latter is for polar coordinates (0 is East, increase counter-clockwise).
}
}
\details{
The bearing, beta, of a point B as seen from a point A is given by

     beta = atan2(S,T)

where

     S = cos(phi_B) * sin(L_A - L_B), and

     T = cos(phi_A)*sin(phi_B) - sin(phi_A)*cos(phi_B)*cos(L_A - L_B)

where phi_A (phi_B) is the latitude of point A (B), and L_A (L_B) is the longitude of point A (B).

Note that there is no simple relationship between the bearing of A to B vs. the bearing of B to A.  The bearing given here is in the usual R convention for lon/lat information, which gives points east of Greenwich as negative longitude, and south of the equator as negative latitude.
}
\value{
numeric giving the bearing angle.
}
\references{
Keay, W. (1995) Land Navigation: Routefinding with Map & Compass, Coventry, UK: Clifford Press Ltd., ISBN 0319008452, 978-0319008454
}
\author{
Eric Gilleland and Randy Bullock, bullock \dQuote{at} ucar.edu
}

\seealso{
\code{\link{atan2}}, \code{\link{FeatureAxis}}, \code{\link{rdist.earth}}
}
\examples{

# Boulder, Colorado and Wallaroo, Australia.
A <- rbind(c(-105.2833, 40.0167), c(137.65, -33.9333))

# Wallaroo, Australia and Boulder, Colorado.
B <- rbind(c(137.65, -33.9333), c(-105.2833, 40.0167))

bearing(A,B)
bearing(A,B,aty="radial")

plot(A, type="n", xlab="", ylab="")
points(A[,1], A[,2], pch="*", col="darkblue")

# Boulder, Colorado to Wallaroo, Australia.
arrows(A[1,1], A[1,2], A[2,1], A[2,2], col="red", lwd=1.5)
}
\keyword{ math }

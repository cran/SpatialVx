\name{CircleHistogram}
\alias{CircleHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Circular Histogram
}
\description{
	Make a circle histogram, also known as a wind-rose diagram.
}
\usage{
CircleHistogram(wspd, wdir, numPetals = 12, radians = FALSE,
		COLS = NULL, scale.factor = 3, varwidth = TRUE,
		minW = NULL, maxW = NULL, circFr = 10,
		main = "Wind Rose", cir.ind = 0.05,
		max.perc = NULL, leg = FALSE, units = "units",
		verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wspd}{numeric vector of length \code{n} giving the magnitudes for the histograms (i.e., the petal colors).  For example, with MODE output, one might use the centroid distances between the paired objects, or the Baddeley Delta metric, or whatever other attribute is of interest.  The argument is called \code{wspd} because in the original context of windrose graphs, the windspeed was used for the petal colors.
}
  \item{wdir}{ numeric vector of length \code{n} giving the angles (from the north) between two matched objects.  In the original context of a windrose diagram, the direction from the north of the wind.  See \code{bearing} for a way to obtain such vectors.
}
  \item{numPetals}{ numeric giving the number of petals to use.  }
  \item{radians}{ logical if TRUE the angles displayed are radians, if FALSE degrees.  }
  \item{COLS}{ vector defining the colors to be used for the petals.  See, for example, \code{topo.colors}.  }
  \item{scale.factor}{numeric determining the line widths (scaled against the bin sizes), only used if \code{varwidth} is TRUE. }
  \item{varwidth}{ logical determining whether to vary the widths of the petals or not.  }
  \item{minW, maxW}{single numerics giving the minimum and maximum break ranges for the histogram of each petal.  If NULL, it will be computed as min( wspd) and max( wspd).}
  \item{circFr}{numeric giving the bin width. }
  \item{main}{character string giving the title to add to the plot.}
  \item{cir.ind}{numeric only used if \code{max.perc} is NULL.}
  \item{max.perc}{numeric giving the maximum percentage to show}
  \item{leg}{logical determining whether or not to add a legend to the plot.}
  \item{units}{character string giving the units for use with the legend.  Not used if \code{leg} is FALSE.}
  \item{verbose}{logical telling whether or not to print information to the screen.}
  \item{\dots}{optional arguments to the \code{plot} or \code{text} functions.}
}
\details{
	The windrose diagram, or circle histogram, is similar to a regular histogram, but adds more information when direction is important.  Binned directions are placed in a full circle with frequencies of for each angle shown via the lengths of each pedal.  Colors along the pedal show the conditional histogram of another variable along the pedal.
}
\value{
	A plot is produced.  If assigned to an object, then a list object is returned with the components:
  \item{summary }{a list object giving the histogram information for each petal.}
  \item{number.obs }{numeric giving the number of observations.}
  \item{number.calm}{numeric giving the number of zero wspd and missing values.}
}
% \references{
%% ~put references to the literature/web site here ~
% }
\author{
	Matt Pocernich
}
\seealso{
\code{\link{bearing}}
}
\examples{
set.seed( 1001 )
wdir <- runif( 1000, 0, 360 )
set.seed( 2002 )
wspd <- rgamma( 1000, 15 )
CircleHistogram( wspd = wspd, wdir = wdir, leg = TRUE )

}
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")

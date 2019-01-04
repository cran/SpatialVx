\name{clusterer}
\alias{clusterer}
\alias{clusterer.default}
\alias{clusterer.SpatialVx}
\alias{summary.clusterer}
\alias{plot.clusterer}
\alias{plot.summary.clusterer}
\alias{print.clusterer}
\title{
Cluster Analysis Verification
}
\description{
Perform Cluster Analysis (CA) verifcation per Marzban and Sandgathe (2006).
}
\usage{
clusterer(X, Y = NULL, ...)

\method{clusterer}{default}(X, Y = NULL, ..., xloc = NULL, xyp = TRUE, threshold = 1e-08, 
    linkage.method = "complete", stand = TRUE, trans = "identity", 
    a = NULL, verbose = FALSE)

\method{clusterer}{SpatialVx}(X, Y = NULL, ..., time.point = 1, obs = 1, model = 1, xyp = TRUE, 
    threshold = 1e-08, linkage.method = "complete", stand = TRUE, 
    trans = "identity", verbose = FALSE)

\method{plot}{clusterer}(x, ..., mfrow = c(1, 2), col = c("gray", tim.colors(64)), 
    horizontal = FALSE)

\method{plot}{summary.clusterer}(x, ...)

\method{print}{clusterer}(x, ...)

\method{summary}{clusterer}(object, ...)
}
\arguments{
  \item{X,Y}{\code{clusterer} default method, these are m by n matrices giving the verification and forecast fields, resp.

\dQuote{SpatialVx} method function, \code{X} is an object of class \dQuote{SpatialVx} and \code{Y} is not used (a warning is given if it is not missing and not NULL).}
  \item{object,x}{list object of class \dQuote{clusterer} as returned by \code{clusterer} (or \code{summary.clusterer} in the case of \code{plot.summary.clusterer}).}
  \item{xloc}{(optional) numeric mn by 2  matrix giving the gridpoint locations.  If NULL, this will be created using 1:m and 1:n.}
  \item{xyp}{logical, should the cluster analysis be performed on the locations and intensities (TRUE) or only the locations (FALSE)?}
  \item{threshold}{numeric of length one or two giving the threshold to apply to each field (>=).  If length is two, the first value corresponds to the threshold for the verification field, and the second to the foreast field.}
  \item{linkage.method}{character naming a valid linkage method accepted by \code{hclust}.}
  \item{stand}{logical, should the data matrices consisting of \code{xloc} and each field first be standardized before performing cluster analysis?}
  \item{trans}{character naming a function to be applied to the field intensities before performing the CA.  Only used if \code{xyp} is TRUE.  Default applies no transformation.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
    \item{a}{(optional) list giving object attributes associated with a \dQuote{SpatialVx} class object.  The \code{clusterer} method for \dQuote{SpatialVx} objects calls the default method function, and uses this argument to pass the attributes through to the final returned object, as well as to grab location information.}
    \item{mfrow}{mfrow parameter (see help file for \code{par}).  If NULL, then the parameter is not re-set.  }
    \item{col}{color vector for image plots of fields after applying the threshold(s).}
    \item{horizontal}{logical, should the image plot color legend be placed horizontally or vertically?  Only for image plot sof the fields.}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{\dots}{optional arguments to the \code{hclust} function.  In the case of the \code{summary} method function, \code{z} and/or \code{sigma} giving a numeric value used to find the cut-off given by median + z*sigma for detemining matched objects (see Marzban and Sandgathe 2006) where defaults of 1 and the standard deviation of minimum inter-cluster distances are used, and \code{silent} (logical should information be printed to the screen (FALSE) or not (TRUE); default is to print to the screen.  In the case of the \code{plot} method functions, these are optional arguments to the \code{summary} method function.}
}
\details{
  This function performs cluster analysis (CA) on positive values from each of two fields in a verification set using the hclust function from package fastcluster.  Inter-cluster distances are computed between each cluster of each field at every level of the CA.  The function clusterer performs CA on both fields, and finds the inter-cluster distances across fields for every possible combination of objects at each iteration of each CA.  The summary method function finishes the analysis by determining hits, misses and false alarms as well as the numbers of clusters.  It also computes CSI for each number of cluster combinations.  This is the verification approach described in Marzban and Sandgathe (2006).

The \code{plot} method function creates a 4 by 2 panel of plots.  The top two plots give image plots of the verification and forecast fields with grid points below the threshold(s) showing zero.  The next two plots are dendrograms as performed by the plot method function for \code{hclust} (\code{dendrogram}) objects.  The next row gives a histogram of the minimum inter-cluster distances, then box plots showing the hits, misses and false alarms for every possible combination of levels of each CA.  Finally, the bottom two plots show, for each combination of CA level (i.e., numbers of clusters), the CSI and average error (inter-cluster distance) for all matched objects.  These last three plots are the ones made by the plot method for values returned from the \code{summary} method function.

\code{print} is currently not very useful here, but it prevents printing a big mess to the screen.
}
\value{
A list object of class \dQuote{clusterer} is returned with components:
\item{linkage.method}{character vector of length one or two giving the linkage method as passed into the function.  The length is two only if the McQuitty method is chosen in which case this method is used for the CA, but not for the inter-cluster differencs across fields (average is used for that instead).}
\item{trans}{character naming the transformation function applied to the intensities.}
\item{N}{numeric giving the size of the fields.}
\item{threshold}{numeric of length two giving the threshold applied to each field.}
\item{NCo,NCf}{numeric vectors giving the number of clusters at each iteration of the CA for the verification and forecast fields, resp.}
\item{cluster.identifiers}{a list with components X and Y giving lists of lists identifying specific CA components at each level of the CA for both fields.}
\item{idX,idY}{logical vectors describing which grid points were included in the CA for each field (i.e., which grid points were >= threshold and had non-missing values).}
\item{cluster.objects}{a list with components X and Y giving the objects returned by hclust for each field.}
\item{inter.cluster.dist}{a list of list objects with NCf by NCo matrix components giving the inter-cluster distances (between verification and forecast fields) for each iteration of CA for each field.}
\item{min.intercluster.dists}{numeric vector givng the minimum values inter.cluster.dist at each iteration.  Used to determine the cut-off for matched objects.}

The summary method function returns a list with the same components as above, but also the components:
\item{cutoff}{The cut-off value used for determining matches.}
\item{csi,AvgErr}{NCo by NCf numeric matrix giving the critical success index (CSI) and average intercluster error (distance) based on matched/un-matched objects.}
\item{HMF}{NCo by NCf by 3 array giving the hits, misses and false alarms based on matched/un-matched objects.}

If the argument a is not NULL, then these are returned as attributes of the returned object.  In the case of \dQuote{SpatialVx} objects, the attributes are preserved.

plot and print methods do not return anything.
}
\references{
Marzban, C. and Sandgathe, S. (2006) Cluster analysis for verification of precipitation fields.  \emph{Wea. Forecasting}, \bold{21}, 824--838.
}
\author{
Eric Gilleland
}
\note{Contact Caren Marzban, marzban \dQuote{at}  u.washington.edu, for questions about the method, and Eric Gilleland, ericg \dQuote{at} ucar.edu, for problems with the code.}
\section{Warning}{
Although some effort has been put into making the functions in this package as computationally efficient as possible, there is a lot of bookeeping involved with this approach, and the current functions are probably not as efficient as they could be.  In any case, they will likely be slow for large data sets.  The function can work quickly on large fields if an adequately high threshold is used (e.g., if threshold=10 is replaced for 16 in the not run example below, the function is VERY slow).  Performing the actual cluster analysis on each field is fast because the hclust function from the fastcluster package is used, which works very well.  However, bookeeping after the CA is done employs a lot of loops within loops, which possibly can be made more efficient (and maybe someday will be), but for now...

If it is desired to simply look at the CA for the two fields, the function hclust from fastcluster can be used, which essentially replaces the hclust function from the stats package with a faster version, but otherwise operates the same as far as what is returned, etc., and the same method functions can be employed.
}

\seealso{
\code{\link[stats]{hclust}}, \code{\link[fastcluster]{hclust}}, \code{\link{as.dendrogram}}, \code{\link{cutree}}, \code{\link{make.SpatialVx}}, \code{\link{CSIsamples}}
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )
look <- clusterer(X=UKobs6, Y=UKfcst6, threshold=16, trans="log", verbose=TRUE)
plot( look )

\dontrun{
data( "UKloc" )

# Now, do the same thing, but using a "SpatialVx" object.
hold <- make.SpatialVx( UKobs6, UKfcst6, loc = UKloc, map = TRUE,
    field.type = "Rainfall", units = "mm/h",
    data.name = "Nimrod", obs.name = "obs 6", model.name = "fcst 6" )

look2 <- clusterer(hold, threshold=16, trans="log", verbose=TRUE)
plot( look2 )
# Note that values differ because now we're using the
# actual locations instead of integer indicators of
# positions.
}
}
\keyword{ cluster }
\keyword{ hplot }

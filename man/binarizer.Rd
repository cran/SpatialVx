\name{binarizer}
\alias{binarizer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create Binary Fields
}
\description{
	Convert a spatial field to a binary spatial field via thresholding.
}
\usage{
binarizer(X, Xhat, threshold = NULL,
	  rule = c(">", "<", ">=", "<=", "<>", "><", "=<>",
		   "<>=", "=<>=", "=><", "><=", "=><="),
	  value = c("matrix", "owin"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X, Xhat}{
	  matrix or \dQuote{owin} objects.
}
  \item{threshold}{
	single number, numeric vector of two numbers, or two-by-two matrix; depending on the value of rule.  May be missing or null if both \code{X} and \code{Xhat} are \dQuote{owin} class objects, in which case binary fields are made wherever these fields are greater than zero, and other arguments are ignored.
}
  \item{rule}{
	  character giving the rule for identifying 1-valued grid squares.  For example, if rule is the default (\dQuote{>}), then threshold should either be a single numeric or a vector of length two.  If the latter, it specifies a different threshold for \code{X} and \code{Xhat}, the result is that everywhere \code{X} > \code{threshold} will have a value of 1 and zero otherwise, etc.  The rule \dQuote{<>} means 1-values whenever \code{X} (\code{Xhat}) are less than the lower threshold and higher than the higher threshold, etc.  For rules requiring two threshold values, if a two-by-two matrix is given, the first column is associated with \code{X} and the second with \code{Xhat}.}
  \item{value}{
	  character telling whether the returned object be a list with two matrices or a list with two \dQuote{owin} class objects.
}
  \item{\dots}{
}
}
\details{

	The binary fields are created by assigning ones according to the rule:
	1. ">": if X > threshold, assign 1, zero otherwise.
        2. "<": if X < threshold, assign 1, zero otherwise.
        3. ">=": if X >= threshold, assign 1, zero otherwise.
        4. "<=": if X <= threshold, assign 1, zero otherwise.
        5. "<>": if X < threshold[ 1 ] or X > threshold[ 2 ], assign 1, zero otherwise.
        6. "><": if threshold[ 1 ] < X < threshold[ 2 ], assign 1, zero otherwise.
	7. "=<>": if threshold[ 1 ] <= X or X > threshold[ 2 ], assign 1, zero otherwise.
	8. "<>=": if threshold[ 1 ] < X or X >= threshold[ 2 ], assign 1, zero otherwise.
        9. "=<>=": if X <= threshold[ 1 ] or X >= threshold[ 2 ], assign 1, zero otherwise.
        10. "=><": if threshold[ 1 ] <= X < threshold[ 2 ], assign 1, zero otherwise.
        11. "><=": if threshold[ 1 ] < X <= threshold[ 2 ], assign 1, zero otherwise.
	12. "=><=": if threshold[ 1 ] <= X <= threshold[ 2 ], assign 1, zero otherwise.
}
\value{
	A list object with two components is returned with the first component being the binary version of X and the second that for Xhat.  These fields will either be matrices of the same dimension as X and Xhat or they will be owin objects depending on the value argument.
}
\author{
	Eric Gilleland
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{thresholder}}, \code{\link{im}}, \code{\link{solutionset}}
}
\examples{
data( "obs0601" )
data( "wrf4ncar0531" )
bin <- binarizer( X = obs0601, Xhat = wrf4ncar0531, threshold = 2.1 )

image.plot( bin[[ 1 ]] )
image.plot( bin[[ 2 ]] )

bin2 <- binarizer( X = obs0601, Xhat = wrf4ncar0531,
		  threshold = 2.1, value = "owin" )
plot( bin2[[ 1 ]] )
plot( bin2[[ 2 ]] )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")

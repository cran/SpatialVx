
\name{calculate_dFSS}
\alias{calculate_dFSS}
\title{
%%  ~~function to do ... ~~
Calculate dFSS Score Value
}
\description{
Calculates the value for the dFSS binary distance metric. The dFSS uses the Fraction Skill Score to provide a measure of spatial displacement of precipitation in two precipitation fields.   
}
\usage{
calculate_dFSS(fbin1, fbin2)
}
\arguments{
  \item{fbin1}{
A numeric matrix representing the first binary field. Only values 0 and 1 are allowed in the matrix.   
}
  \item{fbin2}{
A numeric matrix representing the second binary field.  Only values 0 and 1 are allowed in the matrix. The matrix needs to have the same dimensions as \code{fbin1}.   
}
}
\details{
The dFSS uses the Fraction Skill Score to provide a measure of spatial displacement of precipitation in two precipitation fields. 

The function requires two binary fields as input. A binary field can only have values of 0 or 1 and can be obtained through a thresholding process of the original continuous precipitation field (e.g., by setting all values below a selected precipitation threshold to zero, and all values above the threshold to one). 

The dFSS has a requirement that the frequency bias of precipitation needs to be small in order for the metric to work properly (i.e. the number of non-zero grid points has to be similar in both binary fields). The unbiased fields can be obtained from the original continuous precipitation fields via the use of a frequency (percentile) threshold. For example, instead of using a predefined physical threshold (e.g. 1 mm/h), which might produce binary fields with a different number of non-zero points, a frequency threshold (e.g. 5 \%) can be used which guarantees that both fields will have the same number of non-zero grid-points and will thus be unbiased (provided that enough grid points in the domain contain non-zero precipitation). Function \code{quantile} can be used to determine the value of a physical threshold that corresponds to a prescribed frequency threshold. 

If the frequency bias is larger than 1.5 the function will work but produce a warning. If the frequency bias is larger than 2 the function will produce an error. The dFSS value can only be calculated if both fields contain at least one non-zero grid point. For correct interpretation of the results and some other considerations please look at the "Recipe" in the Conclusions section of Skok and Roberts (2018). 

The code utilizes the fast method for computing fractions (Faggian et al., 2015) and the Bisection method to arrive more quickly at the correct displacement. Optionally, a significantly faster R code that requires significantly less memory and uses some embedded C++ code is available upon request from the author.
}
\value{
The function returns a single numeric value representing the size of the estimated spatial displacement (expressed as a number of grid points - see the example below).
}
\references{
Skok, G. and Roberts, N. (2018), Estimating the displacement in precipitation forecasts using the Fractions Skill Score. Q.J.R. Meteorol. Soc. doi:10.1002/qj.3212. 

Faggian N., Roux B., Steinle P., Ebert B., 2015: Fast calculation of the Fractions Skill Score, MAUSAM, 66 (3), 457-466.
}
\author{
Gregor Skok (Gregor.Skok@fmf.uni-lj.si)
}
\seealso{
\code{\link{hoods2d}}
}
\examples{

# ---------------------------------------------
# A simple example with two 500 x 500 fields
# ---------------------------------------------

# generate two empty 500 x 500 fields where all values are 0 
fbin1=matrix(0, 500, 500, byrow = FALSE)
fbin2=fbin1

# in the fields define a single 20x20 non-zero region of precipitation
# that is horizontally displaced in the second field by 100 grid points   
fbin1[200:220,200:220]=1
fbin2[200:220,300:320]=1

# calulate dFSS value
dFSS=calculate_dFSS(fbin1, fbin2)

# print dFSS value 
print(dFSS)

# The example should output 97 which means that the spatial displacement
# estimated by dFSS is 97 grid points.
}

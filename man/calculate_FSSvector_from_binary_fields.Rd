\name{calculate_FSSvector_from_binary_fields}
\alias{calculate_FSSvector_from_binary_fields}
\title{
Calculate FSS Values
}
\description{
Calculates the value of Fraction Skill Score (FSS) for multiple neighborhood sizes.}
\usage{
calculate_FSSvector_from_binary_fields(fbin1, fbin2, nvector)
}
\arguments{
  \item{fbin1}{
A numeric matrix representing the first binary field. Only values 0 and 1 are allowed in the matrix.   
}
  \item{fbin2}{
A numeric matrix representing the second binary field.  Only values 0 and 1 are allowed in the matrix. The matrix needs to have the same dimensions as \code{fbin1}.   
}
  \item{nvector}{
A numeric vector containing neighborhood sizes for which the FSS values are to be calculated. Only positive odd values are allowed in the vector. A square neighborhood shape is assumed and the specified value represents the length of square side.
}
}
\details{
Fractions Skill Score is a neighborhood-based spatial verification metric frequently used for verifying precipitation (see Roberts and Lean, 2008, for details).

The function requires two binary fields as input. A binary field can only have values of 0 or 1 and can be obtained through a thresholding process of the original continuous precipitation field (e.g., by setting all values below a selected precipitation threshold to zero, and all values above the threshold to one). Either a predefined physical threshold (e.g. 1 mm/h) or a frequency threshold (e.g. 5 \%) can be used to produce the binary fields from the original continuous precipitation fields. If a frequency threshold is used the binary fields will be unbiased and the FSS value will asymptote to 1 at large neighborhoods. Function \code{quantile} can be used to determine the value of a physical threshold that corresponds to a prescribed frequency threshold. 

The code utilizes the fast method for computing fractions (Faggian et al., 2015) that enables fast computation of FSS values at multiple neighborhood sizes. Optionally, a significantly faster R code that requires significantly less memory and uses some embedded C++ code is available upon request from the author.
}
\value{
A numeric vector of the same dimension as \command{nvector} that contains the FSS values at corresponding neighborhood sizes.
}
\references{
Roberts, N.M., Lean, H.W., 2008. Scale-Selective Verification of Rainfall Accumulations from High-Resolution Forecasts of Convective Events. Mon. Wea. Rev. 136, 78-97.

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

# generate two empty 500 x 500 binary fields where all values are 0 
fbin1=matrix(0, 500, 500, byrow = FALSE)
fbin2=fbin1

# in the fields define a single 20x20 non-zero region of precipitation that
# is horizontally displaced in the second field by 100 grid points 
fbin1[200:220,200:220]=1
fbin2[200:220,300:320]=1

# specify a vector of neighborhood sizes for which the FSS values are to be calculated
nvector = c(1,51,101,201,301,601,901,1501)

# calulate FSS values
FSSvector=calculate_FSSvector_from_binary_fields(fbin1, fbin2, nvector)

# print FSS values 
print(FSSvector)

# The example should output:
#  0.00000000 0.00000000 0.04271484 0.52057596 0.68363656 0.99432823 1.00000000 1.00000000
}

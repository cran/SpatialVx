\name{calculate_FSSwind}
\alias{calculate_FSSwind}
\title{
Calculate FSSwind Score Value
}
\description{
Calculates the value for the FSSwind metric that can be used for spatial verification of 2D wind fields. 
}
\usage{
calculate_FSSwind(findex1, findex2, nvector)
}
\arguments{
  \item{findex1}{
A numeric matrix representing the first wind class index field. Only integer values larger than 0 are allowed in the matrix. Each integer value corresponds to a certain wind class.    
}
  \item{findex2}{
A numeric matrix representing the second wind class index field. Only integer values larger than 0 are allowed in the matrix. Each integer value corresponds to a certain wind class. The matrix needs to have the same dimensions as \code{findex1}.   
}
  \item{nvector}{
A numeric vector containing neighborhood sizes for which the FSSwind score value is to be calculated. Only positive odd values are allowed in the vector.
}
}
\details{
The FSSwind is based on the idea of the Fractions Skill Score, a neighborhood-based spatial verification metric frequently used for verifying precipitation. The FSSwind avoids some of the problems of traditional non-spatial verification metrics (the "double penalty" problem and the failure to distinguish between a "near miss" and much poorer forecasts) and can distinguish forecasts even when the spatial displacement of wind patterns is large. Moreover, the time-averaged score value in combination with a statistical significance test enables different wind forecasts to be ranked by their performance (see Skok and Hladnik, 2018, for details).

The score can be used to spatially compare two 2D wind vector fields. In order to calculate the FSSwind value, wind classes have first to be defined. The choice of classes will define how the score behaves and influence the results. The definition of classes should reflect what a user wants to verify. The score will evaluate the spatial matching of the areas of the wind classes. The class definitions should cover the whole phase space of possible wind values (i.e. to make sure every wind vector can be assigned a wind class). It does not make sense choosing a class definition with an overly large number of classes or a definition where a single class would totally dominate over all the other classes. Once the class definition is chosen, the wind vector fields need to be converted to wind class index fields with each class assigned an unique integer value. These wind class index fields serve as input to the calculate_FSSwind function. The FSSwind can have values between 0 and 1 with 0 indicating the worst possible forecast and 1 indicating a perfect forecast.  For guidance on correct interpretation of the results and other important considerations please refer to Skok and Hladnik (2018). 

The code utilizes the fast method for computing fractions (Faggian et al., 2015). Optionally, a significantly faster R code that requires significantly less memory and uses some embedded C++ code is available upon request from the author.
}
\value{
A numeric vector of the same dimension as \code{nvector} that contains the FSSwind values at corresponding neighborhood sizes.
}
\references{
Skok, G. and V. Hladnik, 2018: Verification of Gridded Wind Forecasts in Complex Alpine Terrain: A New Wind Verification Methodology Based on the Neighborhood Approach. Mon. Wea. Rev., 146, 63-75, https://doi.org/10.1175/MWR-D-16-0471.1

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

# generate two 500 x 500 wind class index fields where all values are 1 (wind class 1)
findex1=matrix(1, 500, 500, byrow = FALSE)
findex2=findex1

# in the fields generate some rectangular areas with other wind classes (classes 2,3 and 4)
findex1[001:220,200:220]=2
findex1[100:220,300:220]=3
findex1[300:500,100:200]=4

findex2[050:220,100:220]=2
findex2[200:320,300:220]=3
findex2[300:500,300:500]=4

# specify a vector of neighborhood sizes for which the FSSwind values are to be calculated
nvector = c(1,51,101,201,301,601,901,1501)

# calulate FSSwind values
FSSwindvector=calculate_FSSwind(findex1, findex2, nvector)

# print FSSwind values 
print(FSSwindvector)

# The example should output:
# 0.6199600 0.6580598 0.7056385 0.8029494 0.8838075 0.9700274 0.9754587 0.9756134
}

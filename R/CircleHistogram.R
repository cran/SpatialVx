CircleHistogram <- function(wspd, wdir, numPetals = 12, radians=FALSE,
                         COLS = NULL, scale.factor = 3, varwidth = TRUE,
			 minW=NULL, maxW=NULL,
                         circFr = 10, main='Wind Rose', cir.ind =
                         0.05, max.perc = NULL, leg  = FALSE, units
                         = "units", verbose=FALSE, ...){
##
## Function to create a circle histogram.  Originally called 'windrose.alt', and written by
## Matt Pocernich.
##
## 'wspd' numeric vector of length 'n' giving the magnitudes for the histograms (i.e., the petal colors).
##	For example, with MODE output, one might use the centroid distances between the paired objects, or
##	the Baddeley Delta metric, or whatever other attribute is of interest.  The argument is called 'wspd'
##	because in the original context of the 'windrose.alt' function, one used the windspeed for the petal
##	colors.
## 'wdir' numeric vector of length 'n' giving the angles (from the north) between two matched objects.  In
##	the original context of 'windrose.alt' this would be the direction from the north of the wind.
## 'numPetals' numeric giving the number of petals to use.
## 'radians' logical if TRUE the angles displayed are radians, if FALSE degrees.
## 'COLS' vector defining the colors to be used for the petals.  See, for example, 'topo.colors'.
## 'scale.factor' numeric determining the line widths (scaled against the bin sizes), only used if 'varwidth'
##	is TRUE.
## 'varwidth' logical determining whether to vary the widths of the petals or not.
## 'minW', 'maxW', single numerics giving the minimum and maximum break ranges for each petal's histogram.  If NULL,
##	this will be computed as the min( wspd) and the max( wspd).
## 'circFr' numeric giving the bin width.
## 'main' character string giving the title to add to the plot.
## 'cir.ind' numeric only used if 'max.perc' is NULL.
## 'max.perc' numeric
## 'leg' logical determining whether or not to add a legend to the plot.
## 'units' character string giving the units for use with the legend.  Not used if 'leg' is FALSE.
## '...' optional arguments to the 'plot' or 'text' functions.
##
## Value:
## 
## A plot is produced.  If assigned to an object, then a list object is returned with the components:
##
## 'summary' a list object giving the histogram information for each petal.
## 'number.obs' numeric giving the number of observations (i.e., 'n' from above).
## 'number.calm' numeric giving the number of zero 'wspd' and NA values.
##

## sample data
#wdir      <- runif(1000, 0, 360)
#wspd      <-  rgamma(1000, 15)
#numPetals <- 12
#radians   <- FALSE
#circFr    <- 5  ## bin width
#max.bins  <- NULL ## allow user to specify the max number of bins
#cir.ind   <- 0.05
#max.perc <- 0.3
#scale.factor <- 3
#main <- "test"
#COLS <- NULL
#leg <- TRUE
#units <- "units"
  
 #Checks to see if dir in Degrees or Radians

if (radians){
  wdir <- wdir * 180/pi
}

## assumes data has 2 columns names dir and mag  
op<- par(mar = c(1,1,2,1))
on.exit(par(op))
### rm NA ~ calms ###


#### remove na's and count calms
#### following NWS protocol, when dir == 0 weather calm

no   <- length(wspd)  ## 
ind  <- wspd == 0 | is.na(wdir) 
calm <- sum( ind)
wspd <- wspd[!ind]
wdir <- wdir[!ind]

###### internal function used to plot circles
circles<- function(rad, fill = FALSE){ ## draws circle, radius in units of plot
  ## internal function for windrose
x <- rad* cos(seq(0, 2 * pi, length = 1000) )
y <- rad*sin(seq(0, 2 * pi, length = 1000) )
lines(x,y,col = 1, lty = 2)
if(fill){polygon(x,y, xpd= FALSE, col = "white")}
}
####
step <- 360/numPetals

ddir <- (wdir+step/2)%%360# values rounded up go to Sector 0
ddir[ddir == 0] <- 360  ## anything exactly 0, becomes 360

##########################

aa<- seq(0,360, step)
theta <- aa[-length(aa)]
ind1 <-  as.numeric(cut(ddir, aa)) ## dir indicator

## make bb large enought to include max.bin

N <- length(wdir) ## number of records


## calculate breaks
# bb  <- seq(0, max(wspd) + circFr, by = circFr)
if( is.null( minW) & is.null( maxW)) bb <- seq( min( wspd), max( wspd) + circFr, by=circFr)
else if( is.null( minW))  bb <- seq( min( wspd), maxW + circFr, by=circFr)
else if( is.null( maxW)) bb <- seq( minW, max( wspd) + circFr, by=circFr)
else bb <- seq( minW, maxW + circFr, by=circFr)

max.bins<- length(bb)- 1  ## number of bins

DD <- matrix(NA, nrow = numPetals, ncol = max.bins)
## a matrix with a row for each  petal, column for each bin
for(i in 1:numPetals){ 
if( verbose) cat(i, " ")
 sub    <- wspd[ind1 == i]
DD[i,]  <- hist(sub, breaks = bb, plot = FALSE)$counts/ length(wspd)
}

L       <- apply(DD, 1, sum) ## overall length of each petal

PLPT <- matrix(0, nrow = numPetals, ncol = 4 + 2 * (max.bins))

PLPT[,1] <- theta ## plot angle
PLPT[,2] <- L/sum(L)    ##
indx <- 5
for(i in 1:(max.bins)){
  ll <- apply(matrix(DD[,1:i], ncol = i),1 ,sum)
  PLPT[, indx]    <- sin(theta/360*2*pi)* ll # x plot point for each bin
  PLPT[, indx +1] <- cos(theta/360*2*pi)* ll # y plot point for each bin
  indx            <- indx +2
}

############

if(is.null(max.perc)){
cc  <- max(abs(range(PLPT[, -c(1:4)] ) ) )
cc <- ceiling(cc*(1/cir.ind))* cir.ind
}else
{ cc <- max.perc}

par( lend = 1 )  ## butted line ends
if(leg){par(mar = c(3,3,3,10)) }

plot(0,0, type = "n", ylim = c(-cc, cc),  xlim = c(-cc, cc ),
     main = main, axes = FALSE, xlab = "", ylab = "", ... )
box()

indx <- 3
indy <- 4

if(is.null(COLS) ){COLS <- rainbow(max.bins)}

for(k in 1:(max.bins) ){

  for(j in 1:nrow(PLPT)){
if(varwidth){X<- k}else {X<- 1}
    lines(x = PLPT[j,c(indx, indx + 2) ], y = PLPT[j,c(indy, indy+2)],
      lwd = X*scale.factor, col = COLS[k])
} ## close j

indx <- indx + 2
indy <- indy + 2

} ## close k

CC <- seq(0, cc, cir.ind) 
for(i in 1:length(CC) ) {
  circles(CC[i])
text(0, CC[i], paste(CC[i]*100, "%") )
}

 text(x = c(cc, 0, -cc, 0), y = c(0, cc, 0, -cc ),
      c("90", "0", "270", "180"),
     font = 1, pos = c(3, 4, 3, 4) , offset = 0, ...)

############### legend
if(leg){
  if(!is.null(max.perc)){ cc <- max.perc}
TXT <- paste( round( bb[-length(bb)], digits = 2 ), "to", round( bb[-1], digits = 2 ), units )
par(xpd = NA)
legend(x = cc+ cc*0.15, y = 0, yjust = 0.5,  legend = TXT, pch = 16, col = COLS)

} ## close legend
##############

OUT<- DD # round(new/sum(OUT[m,]),3)

rownames(OUT) <- paste(aa[-length(aa)], "deg")

colnames(OUT) <- paste(  bb[-length(bb)], "to", bb[-1] )
                       
r <- list(summary = OUT, number.obs = no, number.calm = calm)
invisible(r)

} ## close function

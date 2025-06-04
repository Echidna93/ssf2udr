library(terra)
library(raster)
# casts one distribution to another
#' @param n size of dist
#' @param dist initial distribution
#' @param method string for method used by base::rank function
#' @export
cast <- function(n, dist, method = "average"){
  # generate a vector of places "ranks" based on value in the dist
  ranks <- rank(dist, ties.method = method)
  # tranform to a normal distribution
  qnorm(ranks / (n + 1))
}


#' initiates a continous landscape matrix with habitat values between 0 and 1
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrix <- function(nrow, ncol, seed){
  set.seed(seed)
  n <- nrow*ncol
  uVals <- runif(n, 0, 1) # uniform distribution
  normVals <- cast(n, uVals)
  # convert to ncol * nrow matrix
  land <- matrix(normVals, ncol = ncol, nrow = nrow)
  # normalize by the range to keep dist of values the same
  rangeNormalize(land)
}
#' initiates a continuous landscape matrix with increasing habitat values
#' (from left to right) between zero and one
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrixIncreasing <- function(nrow, ncol){
  # make landscape of size nrow*ncol with values of 0
  m <- matrix(sample(0, nrow*ncol, replace = TRUE), nrow = nrow, ncol = ncol)
  j <- 1 # iterator
  for(i in seq(ncol,((ncol*nrow) - ncol), by = ncol)){
    m[(i + 1):(i+ncol)] <- j * 0.02
    j <- j + 1
  }
  m
}

rangeNormalize <- function(land){
  return(((land - min(land))/(max(land) - min(land))))
}

# smooths by the padded matrix by a factor sf using the terra::focal function
#' @param pad padded matrix
#' @param sf smoothing factor; size of focal window
#' @param land matrix to smooth
#' @export
smooth_pad_terra <- function(pad, sf, land){
  w <- raster::focalWeight(rast(pad), sf, type = "circle")
  smooth <- terra::focal(rast(pad), w = w, fun = "mean")
  
  smooth <- smooth[(sf + 1):(nrow(smooth) - sf),(sf + 1):(nrow(smooth) - sf)]
  smooth <- matrix(smooth$focal_mean, nrow = nrow(land), ncol = ncol(land))
  # smooth[smooth > mean(smooth)] <- 1
  # smooth[smooth <= mean(smooth)] <- 0
  rangeNormalize(matrix(as.vector(unlist(smooth)), nrow = nrow(land), ncol = ncol(land)))
}

createBootstrap <- function(smooths, thinVals, nboot, m){
  pops <- data.frame(matrix(NA, 0, 4))
  names(pops) <- c("betaISSF", "smoothingFctr", "nThin", "moransI")
  for(i in 1:length(smooths)){
    for(j in 1:length(thinVals)){
      for(k in 1:nboot){
        population <- metaDat[which(metaDat$smoothingFctr == smooths[i] & metaDat$nThin == thinVals[j]),]$betaISSF
        pop.morans <- metaDat[which(metaDat$smoothingFctr == smooths[i] & metaDat$nThin == thinVals[j]),]$moransI
        nobs <- length(population)
        bootdat <- population[sample(1:nobs, nobs, replace=TRUE)]
        boot.moran <- pop.morans[sample(1:nobs, nobs, replace=TRUE)] 
        pops <- rbind(pops, cbind(mean(bootdat), smooths[i], thinVals[j], mean(boot.moran)))
      }
    }
  }
  return(pops)
}
# creates a matrix which is padded by 2x the smoothing function 
# on each side of the original domain
#' @param land the original matrix
#' @param sf the number of rows to pad by
#' @export
createPaddedMatrix <- function(land, sf){
  pad_land <- matrix(NA, nrow(land) + 2*sf, ncol(land) + 2*sf)
  pad_land[(sf + 1):(nrow(land) + sf), (sf + 1):(ncol(land) + sf)] <- land
  # top left
  pad_land[1:sf, 1:sf] <- land[(nrow(land)-sf + 1):nrow(land), (ncol(land)-(sf ) + 1):ncol(land)]
  # top right
  pad_land[1:sf, (nrow(pad_land) - sf + 1):(nrow(pad_land))] <- land[(nrow(land)-sf + 1):nrow(land), 1:sf] 
  # bottom left
  pad_land[(nrow(pad_land) - sf + 1):(nrow(pad_land)), 1:sf] <- land[1:sf, (ncol(land) - sf + 1):ncol(land)]
  # bottom right
  pad_land[(nrow(pad_land)- sf + 1):(nrow(pad_land)), (nrow(pad_land) - sf + 1):(nrow(pad_land))] <- land[1:sf, 1:sf]
  # bottom row
  pad_land[(nrow(pad_land)-sf + 1):nrow(pad_land), (sf + 1):(ncol(pad_land) - sf)] <- land[1:sf,]
  # top row
  pad_land[1:sf, (sf + 1):(ncol(pad_land) - sf)] <- land[(nrow(land) - sf + 1):nrow(land),]
  # left
  pad_land[(sf + 1):(ncol(pad_land) - sf),1:sf] <- land[,(nrow(land) - sf + 1):nrow(land)]
  # right
  pad_land[(sf + 1):(ncol(pad_land) - sf),(ncol(pad_land) - sf + 1):(ncol(pad_land))]  <- land[,1:sf]
  pad_land
}

clampToSize <- function(dir, ncol){
  # if(dir > ncol){
  #   return(dir - ncol)
  # }
  # else if(dir <= 0){
  #   return(dir + ncol)
  # }
  # return(dir)
  dir_int <- as.integer(dir)
  ((dir_int - 1) %% ncol) + 1
}
checkTrackUD <- function(realizations, nrow, ncol){
  ud <- matrix(0, nrow = nrow, ncol = ncol)
  for(i in 1:nrow(realizations)){
    ud[realizations[i,]$x.proj, realizations[i,]$y.proj] <- ud[realizations[i,]$x.proj, realizations[i,]$y.proj] + 1
  }
  return(ud)
}

iterativeSmooth <- function(land, sf, rho, nrow, ncol){
  smooth <- land
  for(i in seq_len(rho)){
    pad <- createPaddedMatrix(smooth, sf)
    smooth <- smooth_pad_terra(pad,sf,smooth)
    smooth <- matrix(cast(ncell(smooth), as.vector(smooth)), nrow, ncol)
    smooth <- rangeNormalize(smooth)
  }
  return(smooth)
}

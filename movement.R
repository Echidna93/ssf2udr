library(terra)
library(raster) # for ncell function
library(amt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel) # for running foreach in parallel
library(MASS) # for mvrnorm
library(survival) # for clogit
#' initiates a continous landscape matrix with habitat values between 0 and 1
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrix <- function(nrow, ncol, binary=TRUE){
  matrix(runif(nrow*ncol, 0, 1), ncol = ncol, nrow = nrow)
}
#' initiates a continuous landscape matrix with increasing habitat values
#' (from left to right) between zero and one
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrixIncreasing <- function(nrow, ncol, binary=TRUE){
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
  return((land - min(land))/(max(land) - min(land)))
}

# smooths by the padded matrix by a factor sf using the terra::focal function
#' @param pad padded matrix
#' @param sf smoothing factor; size of focal window
#' @param land matrix to smooth
#' @export
smooth_pad_terra <- function(pad, sf, land){
  smooth <- terra::focal(rast(pad), w = sf, fun = "mean",
                         NAonly = TRUE, padValue = 0)
  smooth <- smooth[(sf + 1):(nrow(smooth) - sf),(sf + 1):(nrow(smooth) - sf)]
  # smooth <- matrix(smooth$focal_mean, nrow = nrow(land), ncol = ncol(land))
  # smooth[smooth > mean(smooth)] <- 1
  # smooth[smooth <= mean(smooth)] <- 0
  rangeNormalize(matrix(as.vector(unlist(smooth)), nrow = nrow(land), ncol = ncol(land)))
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

squishToSize <- function(dir){
  if(dir > ncol){
    return(dir - ncol)
  }
  else if(dir <= 0){
    return(dir + ncol)
  }
  return(dir)
}
checkTrackUD <- function(track){
  mat <- matrix(0, nrow = nrow, ncol = ncol)
  #realSteps <- cbind(track$x, track$y)
  for(i in 1:nrow(track)){
    mat[track[i,]$x.proj, track[i,]$y.proj] <- mat[track[i,]$x.proj, track[i,]$y.proj] + 1
  }
  return(mat)
}
# CONSTANTS --------------------------------------------------------------------

nrow <- 100
ncol <- 100
nsims <- nrow*ncol
startTime <- as.POSIXct("2016-11-07 00:00:00 UTC")
smoothingFactorL <- c(1,3,5,7)
nreps = 100
ntraj <- 10
lvars <- 5
sigmaSqEta <- 0.2
nburnin <- 5000
drctnlPers <- 2
drctnPrev <- 'stay'
out.dat <- data.frame(matrix(nrow = 0, ncol = 4))
# create ID for the replicate
# replicate - smoothingFactor - beta
names(out.dat) <- c("t",
                    "cell",
                    "xMod",
                    "yMod")
# metaDat holds data on regression coefficents
metaDat  <- data.frame(matrix(nrow = 0, ncol = 14))

names(metaDat) <- c("rep",
                    "beta", # beta1 is the assigned coeff
                    "slctnCoff", # selection coeff is the retrieved from regression
                    "sl_obs", 
                   "smoothngFctr",
                    "sl_",
                    "shape",
                    "scale",
                    "var_log_sl_",
                    "var_sl_log_sl_",
                    "moransI",
                    "movePenalty",
                    "nThin"
)
landscape <- makeLandscapeMatrix(nrow = nrow, ncol = ncol)
xyTrackDat <- data.frame(matrix(NA, nrow = 0, ncol = 2))
names(xyTrackDat) <- c("x", "y")
# set up file directory
path <- "./data/output" # main file path
# make the domain dir
domName <- paste("domain", ncol, "by", nrow, sep = "-")
if(!dir.exists(paste0(path, "/",domName))){
  dir.create(paste0(path, "/", domName)) # create dir
}
path <- paste0(path, "/", domName)
# create and register worker nodes
# cl <- parallel::makeCluster(4)
# registerDoParallel(cl)

# SIMULATION -------------------------------------------------------------------
for(h in 1:1){
  if(h == 1){
    betaOne <- c(0,1,2,3,4)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- rep(3, lvars)
    thinVals <- rep(150, lvars)
  }
  else if(h == 2){
    betaOne <- rep(0, lvars)
    movePenalties <- c(0.5, 0.25, 1)
    smoothingFactorL <- rep(3, lvars)
    thinVals <- rep(150, lvars)
  }
  else if(h == 3){
    betaOne <- rep(0.5, lvars)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- c(1,3,5)
    thinVals <- rep(150, lvars)
  }
  else if(h == 4){
    betaOne <- rep(0.5, lvars)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- rep(3, lvars)
    thinVals <- rep(100, 150, 200)
  }
for(i in 1:500){
  p <- sample(c(1:lvars), 1, replace = TRUE)
  smoothingFactor <- smoothingFactorL[p] 
  movePen <- movePenalties[p] 
  theta <- betaOne[p]
  nThin <- thinVals[p] # grab thinning value for this iteration
  print(theta)
  if(smoothingFactor == 1){
    landscape_smooth <- landscape
  }else{
    pad <- createPaddedMatrix(landscape, smoothingFactor)
    landscape_smooth <- smooth_pad_terra(pad, smoothingFactor, landscape)
  }
  for(k in 1:1){
    # reset xyTrackDat
    xyTrackDat <- data.frame(matrix(NA, nrow = 0, ncol = 5))
    names(xyTrackDat) <- c("x.proj", "y.proj","x.real","y.real","t")
    startTime <- as.POSIXct("2016-11-07 00:00:00 UTC")
    x.init <- sample(1:ncol, 1, replace = TRUE)
    y.init <- sample(1:ncol, 1, replace = TRUE)
    # get init location
    xyTrackDat[1,]$x.proj = x.init
    xyTrackDat[1,]$y.proj = y.init
    
    xyTrackDat[1,]$x.real = x.init
    xyTrackDat[1,]$y.real = y.init
    
    xyTrackDat[1,]$t <- startTime
    
    up.x <- xyTrackDat[1,]$x.proj- 1
    up.y <- xyTrackDat[1,]$y.proj
    left.x <- xyTrackDat[1,]$x.proj
    left.y <- xyTrackDat[1,]$y.proj- 1
    right.x <- xyTrackDat[1,]$x.proj
    right.y <- xyTrackDat[1,]$y.proj+ 1 
    down.x <- xyTrackDat[1,]$x.proj+ 1
    down.y <- xyTrackDat[1,]$y.proj 
    stay.x <- xyTrackDat[1,]$x.proj
    stay.y <- xyTrackDat[1,]$y.proj
    # real coordinates
    up.real.x <- up.x
    up.real.y <- up.y
    left.real.x <- left.x
    left.real.y <- left.y
    right.real.x <- right.x
    right.real.y <- right.y 
    down.real.x <- down.x
    down.real.y <- down.y 
    stay.real.x <-stay.x
    stay.real.y <-stay.y 
    
    for(iter in 1:(10000)){
      # make decision
      currTime <- startTime
      up.x <- squishToSize(up.x)
      up.y <- squishToSize(up.y)
      down.x <- squishToSize(down.x)
      down.y <- squishToSize(down.y)
      left.x <- squishToSize(left.x)
      left.y <- squishToSize(left.y)
      right.x <- squishToSize(right.x)
      right.y <- squishToSize(right.y)
      stay.x <- squishToSize(stay.x)
      stay.y <- squishToSize(stay.y)
      
      up.val <- landscape_smooth[up.x, up.y]
      down.val <- landscape_smooth[down.x, down.y]
      left.val <- landscape_smooth[left.x, left.y]
      right.val <- landscape_smooth[right.x, right.y]
      stay.val <- landscape_smooth[stay.x, stay.y]
      
      weights <- c(exp(-movePen + (up.val * theta)), # up
                   exp(-movePen + (down.val * theta)), # down
                   exp(-movePen + (left.val * theta)), # left
                   exp(-movePen + (right.val * theta)), # right
                   exp(0 + (stay.val * theta))) # stay
      
      locs.proj <- list('up' = list("x" = up.x, "y" = up.y),
                'down' = list('x' = down.x, 'y' = down.y),
                'left' = list('x' = left.x, 'y' = left.y),
                'right' = list('x' = right.x, 'y' = right.y),
                'stay' = list('x' = stay.x, 'y' = stay.y))
      
      locs.real <- list('up' = list("x" = up.real.x, "y" = up.real.y),
                        'down' = list('x' = down.real.x, 'y' = down.real.y),
                        'left' = list('x' = left.real.x, 'y' = left.real.y),
                        'right' = list('x' = right.real.x, 'y' = right.real.y),
                        'stay' = list('x' = stay.real.x, 'y' = stay.real.y))
      
      decision <- sample(c(1:5), 1, prob = weights, replace = TRUE)
      # print(decision)
      # print(weights)
      # print()
      xyTrackDat[iter + 1, ]$x.proj<- locs.proj[[decision]]$x
      xyTrackDat[iter + 1, ]$y.proj<- locs.proj[[decision]]$y
      
      xyTrackDat[iter + 1, ]$x.real<- locs.real[[decision]]$x
      xyTrackDat[iter + 1, ]$y.real<- locs.real[[decision]]$y
      
      # grab new directions
      up.x <- xyTrackDat[iter + 1,]$x.proj- 1
      up.y <- xyTrackDat[iter + 1,]$y.proj
      left.x <- xyTrackDat[iter + 1,]$x.proj
      left.y <- xyTrackDat[iter + 1,]$y.proj- 1
      right.x <- xyTrackDat[iter + 1,]$x.proj
      right.y <- xyTrackDat[iter + 1,]$y.proj + 1 
      down.x <- xyTrackDat[iter + 1,]$x.proj+ 1
      down.y <- xyTrackDat[iter + 1,]$y.proj 
      stay.x <- xyTrackDat[iter + 1,]$x.proj
      stay.y <- xyTrackDat[iter + 1,]$y.proj
      
      up.real.x <- xyTrackDat[iter + 1,]$x.real- 1
      up.real.y <- xyTrackDat[iter + 1,]$y.real
      left.real.x <- xyTrackDat[iter + 1,]$x.real
      left.real.y <- xyTrackDat[iter + 1,]$y.real- 1
      right.real.x <- xyTrackDat[iter + 1,]$x.real
      right.real.y <- xyTrackDat[iter + 1,]$y.real+ 1 
      down.real.x <- xyTrackDat[iter + 1,]$x.real+ 1
      down.real.y <- xyTrackDat[iter + 1,]$y.real 
      stay.real.x <- xyTrackDat[iter + 1,]$x.real
      stay.real.y <- xyTrackDat[iter + 1,]$y.real

      currTime <- as.POSIXct(currTime) + lubridate::minutes(1)
      xyTrackDat[iter + 1, ]$t <- currTime # update time
    }
    # ANALYSIS ---------------------------------------------------------------------
    
    # thin the movement dataset
    ## remove burnin
    xyTrackDat <- xyTrackDat[nburnin:nrow(xyTrackDat),]
    # ## remove every nth sample
    out.dat <- out.dat[seq(1:nrow(out.dat)) %% nThin == 0,]   
    xyTrackDat2 <- xyTrackDat
    xyTrackDat <- xyTrackDat[seq(1:nrow(xyTrackDat)) %% nThin == 0,]  
    # project cells to e-space
    
    # REAL TRACK
    ## get sls from real track
    trk.real <- make_track(as_tibble(xyTrackDat), .x = x.real,
                       .y = y.real,
                       .t = t) 
    
    trk.real$x_ <- trk.real$x_ + rnorm(nrow(trk.real), 0, sigmaSqEta)
    trk.real$y_ <- trk.real$y_ + rnorm(nrow(trk.real), 0, sigmaSqEta)
    
    stps.real <- trk.real %>% steps() 
    
    
    # PROJECTED TRACK
    # turn x,y proj into tracks object
    trk.proj <- make_track(as_tibble(xyTrackDat), .x = x.proj,
                           .y = y.proj,
                           .t = t) 
    # add gaussian noise
    trk.proj$x_ <- trk.proj$x_ + rnorm(nrow(trk.proj), 0, sigmaSqEta)
    trk.proj$y_ <- trk.proj$y_ + rnorm(nrow(trk.proj), 0, sigmaSqEta)
    
    stps.proj <- trk.proj %>% steps()
    stps.proj$sl_ <- stps.real$sl_ # add real sls back
    
    stps <- stps.proj %>% random_steps(n_control = 30)
  
    # fit ISSF
    # 
    # squash steps outside the domain back in
    stps$x2_ <- stps$x2_ %% ncol # row
    stps$y2_ <- stps$y2_ %% ncol # row
    stps$x1_ <- stps$x1_ %% ncol # row
    stps$y1_ <- stps$y1_ %% ncol # row
    
    stps <- stps %>% amt::extract_covariates(rast(landscape_smooth)) # extract covars
    
    stps <- stps %>% mutate(log_sl_ = log(sl_),
                            land = lyr.1)
    
    #hist(landscape_smooth[cbind(stps[which(stps$case_),]$x1_, stps[which(stps$case_),]$y1_)])
    mod <- amt::fit_issf(stps, case_ ~ log_sl_ + sl_ + land + strata(step_id_))
    print(mod$model$coefficients[3])
    # grab sl_ and log_sl_ distr
    norms <- mvrnorm(ntraj, cbind(mod$model$coefficients['log_sl_'],
                                  mod$model$coefficients['sl_']),
                     vcov(mod$model)[1:2, 1:2])
    scale <- update_sl_distr(mod, log_sl_ = norms[1:nrow(norms), 1], sl_ = norms[1:nrow(norms), 2])$params$scale
    shape <- update_sl_distr(mod, log_sl_ = norms[1:nrow(norms), 1], sl_ = norms[1:nrow(norms), 2])$params$shape
    # plot(rast(matrix(transDat$num, ncol = ncol, nrow = nrow)))
    metaDat <- rbind(metaDat,cbind(
      rep = i,
      beta1 = theta,
      slctnCoff = unlist(mod$model$coefficients[3]), # grab LS regression coefficient
      var_slctn = unlist(vcov(mod$model)[3,3]), # grab variance
      sl_ = unlist(mod$model$coefficients['sl_']),
      sl_obs = 1, 
      scale = scale, # grab variance
      shape = shape,
      var_log_sl_ = unlist(vcov(mod$model)[1,1]), # grab variance
      var_sl_log_sl_ = unlist(vcov(mod$model)[1,2]),
      smoothingFctr = unlist(smoothingFactor),
      moransI = unlist(Moran(raster(landscape_smooth))),
      movePenalty = unlist(movePen),
      nThin = unlist(nThin)))
    }
  }
}
# PLOTTING ---------------------------------------------------------------------

# Below thins the trajectories exclude trajectories that cross the boundary
out.dat$xMod <- if_else(is.na(out.dat$xMod), 0, out.dat$xMod)

out.dat$yMod <- if_else(is.na(out.dat$yMod), 0, out.dat$yMod)

# create lagged vector for previous statement
out.dat$xModPrev <- data.table::shift(out.dat$xMod, 1, fill = 0)
out.dat$yModPrev <- data.table::shift(out.dat$yMod, 1, fill = 0)
# filter out points where xModPrev doesn't equal xMod
out.dat.filter <- out.dat %>% filter(xModPrev == xMod, yModPrev == yMod)

trk <- make_track(as_tibble(out.dat.filter), .x = x,
                  .y = y,
                  .t = t)
plot(rast(landscape_smooth), xlim = c(1,ncol-1), ylim = c(1,nrow-1),
     col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2,
                       alpha = NULL))
lines(xyTrackDat, col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))
#points(trk[which(trk$case_ == FALSE),],col = "yellow", lwd=1, xlim = c(0,50), ylim=c(0,50))


# PLOT PATH WITH GGPLOT --------------------------------------------------------
# qplot(x, y, data = out.dat, aes(color="red"))+ 
#   geom_path(aes(color="red")) + scale_x_continuous(limits = c(1,ncol)) + scale_y_continuous(limits = c(1,nrow))

# pivot on the data.frame 
j <- landscape_smooth %>% as.data.frame() %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>% mutate(
Var1 = factor(Var1, levels = 1:ncol),
Var2 = factor(gsub("V", "", Var2), levels = 1:nrow))

# plot
ggplot(j, aes(Var1, Var2,)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high="black") +
  geom_path(data = out.dat.filter, aes(x = x, y = y, color = "red")) +
  geom_point(data = out.dat, aes(x = x, y = y, color = "pink"))

# plot on continuous domain
out.dat$coord  <-  paste0("(", out.dat$xMod, ",", out.dat$yMod, ")")
ggplot(j, aes(x = Var1, y = Var2,)) + 
  facet_grid(~ coord) +   
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high="black") +
  geom_path(data = out.dat, aes(x = x, y = y, color = "red")) +
  geom_point(data = out.dat, aes(x = x, y = y, color = "pink"))



# out.dat.long <- out.dat %>% mutate(uid = paste0(rep, "-", smoothingFactor, "-", beta1)) %>%
#   mutate(logRSS = loghttp://127.0.0.1:17001/graphics/e7beda11-ce12-411c-adbf-ce0c4c7660ec.png((hb1*b0)/(hb0*b1))) %>%
#   group_by(uid) %>% mutate(meanlogRSS = mean(logRSS))
# 
ggplot(metaDat, aes(y = unlist(slctnCoff), x = unlist(smoothingFctr), color = factor(unlist(beta1)))) + geom_boxplot()
# ggplot(out.dat.long, aes(x = factor(smoothingFactor), y = logRSS, color = factor(beta1))) + geom_boxplot()
# 
# 
# write.table(out.dat, "./outdat-sf-1:7-50-50", sep = ",")
# hist(pts$lyr.1)


# # PLOT PATH WITH BASE PLOT ---------------------------------------------------

# plot(rast(landscape_smooth), xlim = c(1,ncol), ylim = c(1,nrow), col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
# lines(make_track(as_tibble(out.dat), .x = x, .y = y, .t = t), col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))


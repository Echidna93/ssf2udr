library(terra)
library(raster) # for ncell function
library(amt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel) # for running foreach in parallel
library(MASS) # for mvrnorm
library(survival) # for clogit
#' initiates a landscape matrix of vectors of random 0's and 1's
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrix <- function(nrow, ncol, binary=TRUE){
  #matrix(runif(nrow*ncol, 0, 1), ncol = ncol, nrow = nrow)
  m <- matrix(sample(0, nrow*ncol, replace = TRUE), nrow = nrow, ncol = ncol)
  m[1:ncell(m)/2] <- 0.5
  m
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an makeDecision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param 
#' @export
makeDecision<-function(landscape, betas_l, cell,
                       cells, mods, transition_prob, drctnPrev,
                       drctnlPers){
  # TODO check out sample int
  prob <- transition_prob[cell,] # grab row of transDat
  # check previous direction for drctnl persistance
  # prob[,drctnPrev] <- if_else(drctnPrev!='stay',
  #                                  drctnlPers*prob[,drctnPrev],
  #                                  prob[,drctnPrev])
  selection <- sample(c(2:6), 1, prob = prob[,3:7])
  # selection <- selection + 1 # add one so we don't get the cell reference
  drctn <- colnames(cells)[selection]
  if(drctn == "stay"){
    modVec <- list(ymod = 0, xmod = 0)
  }
  if(drctn == "up" || drctn == "down"){
    modVec <-  list(ymod = mods[cell, selection], xmod = 0)
  }else{
    modVec <- list(ymod = 0, xmod = mods[cell, selection]) 
  }

  c(cell = cells[cell, selection], modVec, drctnPrev)
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an makeDecision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param landscape matrix to smooth
#' @param smoothingFactor size of moving window to smooth
#' @export
generateSteps <- function(smoothingFactor){
  steps <- seq(from = -1, to = 1, by = 1)
  smooth <- seq(1, smoothingFactor, by = 1)
  steps <- expand.grid(steps, steps) %>% as.data.frame() %>% dplyr::rename(x = Var1, y = Var2) %>% filter(x != -y) %>% filter(x != y)
  steps <- rbind(c(0,0),steps)
  for(i in 1:nrow(steps)){
    for(k in 1:length(smooth)){
      steps <- rbind(steps, smooth[k] * steps[i,])
    }
  }
  steps %>% unique()
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
  matrix(as.vector(unlist(smooth)), nrow = nrow(land), ncol = ncol(land))
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

# getMods creates a dataframe that tracks when each step 
# leaves the original domain
#' @param nrow number of rows
#' @param ncol number of columns
#' @param landscape matrix to smooth
#' @export
getMods <- function(nrow, ncol, landscape){
  cells <- data.frame(cell = 1:ncell(landscape), stay = 0,
                      right = 0, left = 0, up = 0, down = 0)
  for(i in 1:nrow(cells)){
    # top left
    if(i == 1){
      cells[i,]$left <- -1
      cells[i,]$up <- 1
      cells[i,]$down <- 0
      cells[i,]$right <- 0
    }
    # bottom left
    else if(i == ncol){
      cells[i,]$left <- -1
      cells[i,]$up <- 0
      cells[i,]$down <- 1
      cells[i,]$right <- 0
    }
    # bottom right
    else if(i == ncol*nrow){
      cells[i,]$left <- 0
      cells[i,]$up <- 0
      cells[i,]$down <- 1
      cells[i,]$right <- 1
    }
    # top right
    else if(i == (nrow*ncol - ncol) + 1){
      cells[i,]$left <- 0 
      cells[i,]$up <- -1
      cells[i,]$down <- 0
      cells[i,]$right <- 1
    }
    # bottom row
    else if(i %% nrow == 0 & (i != ncol & i != ncol*nrow)){
      cells[i,]$left <-  0
      cells[i,]$up <- 0
      cells[i,]$down <- 1
      cells[i,]$right <- 0
    }
    # right edge
    else if((i > (nrow * ncol - ncol) + 1) & (i < nrow*ncol)){
      cells[i,]$left <- 0
      cells[i,]$up <- 0
      cells[i,]$down <- 0
      cells[i,]$right <- 1
    }
    # left edge
    else if((i > 1) & (i < ncol)){
      cells[i,]$left <- -1
      cells[i,]$up <- 0
      cells[i,]$down <- 0
      cells[i,]$right <- 0
    }
    # top edge
    else if((i %% nrow == 1) & (i != 1 & i != (nrow*ncol - ncol) + 1)){
      cells[i,]$left <- 0
      cells[i,]$up <- -1
      cells[i,]$down <- 0
      cells[i,]$right <- 0
    }
    else{
      cells[i,]$left <- 0
      cells[i,]$up <- 0 
      cells[i,]$down <- 0
      cells[i,]$right <- 0
    }
  }
  return(cells)
}

# getSteps creates a dataframe that tracks the index (cell) of the next 
# step for each cell in the domain
#' @param nrow number of rows
#' @param ncol number of columns
#' @param land matrix to smooth
#' @export
getSteps <- function(nrow, ncol, land){
  cells <- data.frame(cell = 1:ncell(land), stay = 1:ncell(land),
                      right = 0, left = 0, up = 0, down = 0)
  for(i in 1:nrow(cells)){
    # top left
    if(i == 1){
      cells[i,]$left <- (ncol*nrow - ncol) + 1
      cells[i,]$up <- ncol
      cells[i,]$down <- i + 1
      cells[i,]$right <- i + ncol
    }
    # bottom left
    else if(i == ncol){
      cells[i,]$left <- ncol * nrow
      cells[i,]$up <- i - 1
      cells[i,]$down <- 1
      cells[i,]$right <- i + ncol
    }
    # bottom right
    else if(i == ncol*nrow){
      cells[i,]$left <- i - ncol
      cells[i,]$up <- i - 1
      cells[i,]$down <- ncol*nrow - ncol + 1
      cells[i,]$right <- ncol
    }
    # top right
    else if(i == (nrow*ncol - ncol) + 1){
      cells[i,]$left <- i - nrow 
      cells[i,]$up <- ncol*nrow
      cells[i,]$down <- i + 1
      cells[i,]$right <- 1
    }
    # bottom edge
    else if(i %% nrow == 0 & (i != ncol & i != ncol*nrow)){
      cells[i,]$left <-  i - nrow
      cells[i,]$up <- i - 1
      cells[i,]$down <- i - ncol + 1
      cells[i,]$right <- i + nrow
    }
    # right edge
    else if((i > (nrow * ncol - ncol) + 1) & (i < nrow*ncol)){
      cells[i,]$left <- i - nrow
      cells[i,]$up <- i - 1
      cells[i,]$down <- i + 1
      cells[i,]$right <- i %% ncol
    }
    # left edge
    else if((i > 1) & (i < ncol)){
      cells[i,]$left <- (ncol*nrow - ncol) + i
      cells[i,]$up <- i - 1
      cells[i,]$down <- i + 1
      cells[i,]$right <- i + ncol
    }
    # top edge
    else if((i %% nrow == 1) & (i != 1 & i != (nrow*ncol - ncol) + 1)){
      cells[i,]$left <- i - ncol
      cells[i,]$up <- i + (ncol - 1)
      cells[i,]$down <- i + 1
      cells[i,]$right <- i + nrow
    }
    else{
      cells[i,]$left <- i - nrow
      cells[i,]$up <- i - 1 
      cells[i,]$down <- i + 1
      cells[i,]$right <- i + nrow
      }
    }
  return(cells)
}

# createTransDat creates a dataframe that the probability of movement for 
# each step on the map
#' @param cells the dataframe containing the cells # of the possible steps
#' @param landscape_smooth smoothed domain from smooth_pad_terra
#' @param betas_l list of beta values for landscape types
#' @param move_penalty the propensity for the agent to move or not move
#' @export
createTransDat <- function(cells, landscape_smooth, theta, move_penalty){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- exp(0 + (as.vector(landscape_smooth) * theta))
  transMat$left <- exp(-move_penalty + (landscape_smooth[cells$left]) * theta)
  transMat$up <- exp(-move_penalty + (landscape_smooth[cells$up] * theta))
  transMat$down <- exp(-move_penalty + (landscape_smooth[cells$down] * theta))
  transMat$right <- exp(-move_penalty + (landscape_smooth[cells$right] * theta))
  return(transMat)
}
createCellDat <- function(cells, landscape_smooth, betas_l, move_penalty){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- as.integer(betas_l[as.character(landscape_smooth)])
  transMat$left <- as.integer(betas_l[as.character(landscape_smooth[cells$left])])
  transMat$up <- as.integer(betas_l[as.character(landscape_smooth[cells$up])])
  transMat$down <- as.integer(betas_l[as.character(landscape_smooth[cells$down])])
  transMat$right <- as.integer(betas_l[as.character(landscape_smooth[cells$right])])
  return(transMat)
}
# helper function for extrat_covariates
# takes columns of points
getBadRows <- function(cols, ncol){
  # coerce cols into dataframe
  cols <- as.data.frame(cols)
  names(cols) <- c("x", "y")
  # check which columns are outside of the range of the matrix
  return(unique(which(cols$y < 1 | cols$y > ncol | cols$x < 1 | cols$x > ncol)))
}

# wrapper around extract covariates
# takes in a list of indicies with points that fall outside
# and coerces them back into the domain
projectPoints <- function(pts){
  
}

# CONSTANTS --------------------------------------------------------------------

nrow <- 50
ncol <- 50
betaOne <- c(2, 2.5, 3, 3.5)
movePenalties <- c(0,0.25, 0.5, 1)
thinVals <- c(100, 150, 200, 250) # number of samples to throw out before writing
betas_l <- list('0' = 0,
                '1' = 1)
nsims <- nrow*ncol
startTime <- as.POSIXct("2016-11-07 00:00:00 UTC")
smoothingFactorL <- c(1,3,5,7)
nreps = 100
ntraj <- 10
lvars <- 4
nburnin <- 10000
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
landscape <- makeLandscapeMatrix(nrow, ncol, TRUE)
cells <- getSteps(nrow, ncol, landscape)
mods <- getMods(nrow, ncol, landscape)

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
    betaOne <- c(0)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- rep(1, lvars)
    thinVals <- rep(100, lvars)
  }
  else if(h == 2){
    betaOne <- rep(0, lvars)
    movePenalties <- c(0.5, 0.25, 0.5, 1)
    smoothingFactorL <- rep(3, lvars)
    thinVals <- rep(150, lvars)
  }
  else if(h == 3){
    betaOne <- rep(0.5, lvars)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- c(1,3,5,7)
    thinVals <- rep(150, lvars)
  }
  else if(h == 4){
    betaOne <- rep(0.5, lvars)
    movePenalties <- rep(0, lvars)
    smoothingFactorL <- rep(3, lvars)
    thinVals <- rep(100, 150, 200, 250)
  }
for(i in 1:1){
  smoothingFactor <- smoothingFactorL[i] 
  movePen <- movePenalties[i] 
  theta <- betaOne[i]
  nThin <- thinVals[i] # grab thinning value for this iteration
  print(theta)
  if(smoothingFactor == 1){
    landscape_smooth <- landscape
  }else{
    pad <- createPaddedMatrix(landscape, smoothingFactor)
    landscape_smooth <- smooth_pad_terra(pad, smoothingFactor, landscape)
  }
  transDat <- createTransDat(cells, landscape_smooth, theta, movePen)
  transDat$num <- 0
  # doParallel routine
  # k is the number of trajectories
  for(k in 1:(nrow*ncol)){
    print(k)
    nThin <- sample(c(10, 50, 100), 1, replace = TRUE)
    out.dat <- data.frame(matrix(nrow = 0, ncol = 4))
    # create ID for the replicate
    # replicate - smoothingFactor - beta
    names(out.dat) <- c("t",
                        "cell",
                        "xMod",
                        "yMod")
    # metaDat holds data on regression
    currTime <- startTime
    sampIter <- 1
    xmod <- 0 # x-mod and y mod need reset for each simulation
    ymod <- 0
    # for randomly sampling the landscape--testing only
    # loc <- data.frame(cell = sample(1:ncell(landscape_smooth), 1, replace = TRUE), ymod = 0, xmod = 0)
    # iteratively stepping thru landscape
    loc <- data.frame(cell = k, ymod = 0, xmod = 0)
    transDat[k,]$num <- transDat[k,]$num + 1
    for(iter in 1:(nrow*ncol)){
      loc <- makeDecision(landscape_smooth, theta, loc$cell, cells, mods, transDat,
                          drctnPrev, drctnlPers)
      drctnPrev <- loc$drctnPrev
      xmod <- xmod + loc$xmod # one of mods is always zero
      ymod <- ymod + loc$ymod
      currTime <- as.POSIXct(currTime) + lubridate::minutes(1)
      # for main sim
      # transDat[loc$cell,]$num <- if_else(iter >= nburnin, transDat[loc$cell,]$num + 1,
      #                                    transDat[loc$cell,]$num)
      
      transDat[loc$cell,]$num <- transDat[loc$cell,]$num + 1
      out.dat <- rbind(out.dat, cbind(
        t = currTime,
        cell = loc$cell,
        xMod = xmod, # keeps track of the number of times agent has passed the x-axis
        yMod = ymod  # keeps track of the number of times agent has passed the y-axis
        ))
    }
    # ANALYSIS ---------------------------------------------------------------------
    # out.dat <- out.dat[nburnin:nrow(out.dat),]
    
    # thin the movement dataset   
    # out.dat <- out.dat[seq(1:nrow(out.dat)) %% nThin == 0,]   

    # make x and y column
    out.dat$x = if_else(out.dat$cell %% ncol != 0, out.dat$cell %% ncol, ncol)
    out.dat$y = if_else(out.dat$cell %% ncol != 0, ceiling(out.dat$cell/ncol), ceiling(out.dat$cell / ncol))
    out.dat$t = as.POSIXct(out.dat$t)

    trk <- make_track(as_tibble(out.dat), .x = x,
                       .y = y,
                       .t = t) 
    stps  <- steps(trk) # %>% random_steps(n_control = 30)
    sl_obs <- stps$sl_
    stps <- random_steps(stps, n_control = 30)
    # # 
    # # # round the decimal places off
    # # # clamp values to size
    stps[which(!stps$case_),]$x2_ <- ceiling(stps[which(!stps$case_),]$x2_ %% ncol) # row
    stps[which(!stps$case_),]$y2_ <- ceiling(stps[which(!stps$case_),]$y2_ %% ncol) # row

    # # # extract landscape values
    stps <- stps %>% amt::extract_covariates(rast(landscape_smooth))
    # 
    # # remove incomplete strata
    # # randStps <- randStps %>% amt::remove_incomplete_strata() # might not need this?
    # # add integer to avoid sl of zero
    # # stps$sl_ <- stps$sl_ + 1
    # 
    stps <- stps %>% mutate(sl_ = sl_ + 1,
                            log_sl_ = log(sl_),
                            land = lyr.1)

    # fit ISSF
    # mod <- amt::fit_issf(stps, case_ ~ log_sl_ + sl_ + land + strata(step_id_))
    # print(mod.surv)
    
    # Below thins the trajectories exclude trajectories that cross the boundary
    # out.dat$xMod <- if_else(is.na(out.dat$xMod), 0, out.dat$xMod)
    # out.dat$yMod <- if_else(is.na(out.dat$yMod), 0, out.dat$yMod)
    # 
    # # create lagged vector for previous statement
    # out.dat$xModPrev <- data.table::shift(out.dat$xMod, 1, fill = 0)
    # out.dat$yModPrev <- data.table::shift(out.dat$yMod, 1, fill = 0)
    # # filter out points where xModPrev doesn't equal xMod
    # out.dat.filter <- out.dat %>% filter(xModPrev == xMod, yModPrev == yMod)
    # 
    # trk <- make_track(as_tibble(out.dat), .x = x,
    #                   .y = y,
    #                   .t = t)
    # plot(rast(landscape_smooth), xlim = c(1,ncol), ylim = c(1,nrow),
    #      col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2,
    #                        alpha = NULL))
    # lines(make_track(as_tibble(out.dat.filter), .x = y, .y = x, .t = t),
    #       col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))
    # 
    # 
    #print(mod$model$coefficients[3])
    
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
      sl_obs = sl_obs, 
      scale = scale, # grab variance
      shape = shape,
      var_log_sl_ = unlist(vcov(mod$model)[1,1]), # grab variance
      var_sl_log_sl_ = unlist(vcov(mod$model)[1,2]),
      smoothingFctr = unlist(smoothingFactor),
      moransI = unlist(Moran(raster(landscape_smooth))),
      movePenalty = unlist(movePen),
      nThin = unlist(nThin)))
      # subDomName <- paste("smoothingFactor", smoothingFactor, sep = "-") # smoothing Factor
      # ssubDomName <- paste("beta1", betas_l$'1', sep = "-") # betas
      # sssubDomName <- paste("movement-penalty", movePen, sep = "-") # movement penalty
      # ssssubDomName <- paste("thinning", nThin, sep = "-")
      # fname <- paste0("-", subDomName, "-", ssubDomName,
      #                 "-", sssubDomName, "-", ssssubDomName, "-movement-data-rep-", i, "-realization-", k, sep = "")
      # #dir.create(fpath, recursive = TRUE) # create dir
      # #ame <- paste0("movement-data-rep-", i, "-realization-", k, sep = "")
      # write.table(out.dat,
      #             file = paste0(path, "/", fname), sep = ",")
      # 
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
lines(trk, col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))
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


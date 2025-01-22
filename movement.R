library(terra)
library(raster) # for ncell function
library(amt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel) # for running foreach in parallel

#' initiates a landscape matrix of vectors of random 0's and 1's
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
makeLandscapeMatrix <- function(numrow, numcol, binary=TRUE){
  x = seq(1, numrow, by = 1)
  y = seq(1, numcol, by = 1)
  coords <- list(expand.grid(x,y))
  vals <- sample(x = c(0,1),size = numrow*numcol, replace = TRUE)
  matrix(vals, ncol = numcol, nrow=numrow)
  
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an makeDecision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param 
#' @export
makeDecision<-function(landscape, betas_l, cell, cells, mods, transition_prob){
  selection <- sample(c(2:6), 1, prob = transition_prob[cell,3:7])
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

  c(cell = cells[cell, selection], modVec)
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
  smooth <- matrix(smooth$focal_mean, nrow = nrow(land), ncol = ncol(land))
  smooth[smooth > mean(smooth)] <- 1
  smooth[smooth <= mean(smooth)] <- 0
  smooth
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

# createTransDate creates a dataframe that the probability of movement for 
# each step on the map
#' @param cells the dataframe containing the cells # of the possible steps
#' @param landscape_smooth smoothed domain from smooth_pad_terra
#' @param betas_l list of beta values for landscape types
#' @param move_penalty the propensity for the agent to move or not move
#' @export
createTransDat <- function(cells, landscape_smooth, betas_l, move_penalty){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- exp(0 + c((landscape_smooth) * as.integer(betas_l[as.character(landscape_smooth)])))
  transMat$left <- exp(-move_penalty + c((landscape_smooth[cells$left]) * as.integer(betas_l[as.character(landscape_smooth[cells$left])])))
  transMat$up <- exp(-move_penalty + c((landscape_smooth[cells$up]) * as.integer(betas_l[as.character(landscape_smooth[cells$up])])))
  transMat$down <- exp(-move_penalty + c((landscape_smooth[cells$down]) * as.integer(betas_l[as.character(landscape_smooth[cells$down])])))
  transMat$right <- exp(-move_penalty + c((landscape_smooth[cells$right]) * as.integer(betas_l[as.character(landscape_smooth[cells$right])])))
  return(transMat)
}

# CONSTANTS --------------------------------------------------------------------

nrow <- 100
ncol <- 100
betaOne <- c(1,2, 2.5)
movePenalties <- c(0, 0.5, 1)
thinVals <- c(100, 150, 200) # number of samples to throw out before writing
betas_l <- list('0' = 0,
                '1' = 1)
nsims <- nrow*ncol
startTime <- as.POSIXct("2016-11-07 00:00:00 UTC")
smoothingFactorL <- c(1,3,5,7)
vars <- expand.grid(betaOne, movePenalties, thinVals, smoothingFactorL)
vars <- vars[sample(1:nrow(vars)), ]
nreps = 100
out.dat <- data.frame(matrix(nrow = 0, ncol = 4))
# create ID for the replicate
# replicate - smoothingFactor - beta
names(out.dat) <- c("t",
                    "cell",
                    "xMod",
                    "yMod")
# metaDat holds data on regression coefficents
metaDat  <- data.frame(matrix(nrow = 0, ncol = 7))
names(metaDat) <- c("rep",
                    "beta", # beta1 is the assigned coeff
                    "slctnCoff", # selection coeff is the retrieved from regression
                    "smoothngFctr",
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
for(i in 69:nrow(vars)){
  landscape <- makeLandscapeMatrix(nrow, ncol, TRUE)
  smoothingFactor <- vars[i,4] #smoothingFactorL[sample(1:length(smoothingFactorL), 1)]
  movePen <- vars[i,2] #movePenalties[sample(1:length(movePenalties), 1)]
  betas_l$'1' <- vars[i,1] #betaOne[sample(1:length(betaOne), 1)]
  nThin <-  vars[i,3] # thinVals[sample(1:length(thinVals), 1)] # grabbing the thinning value for this iteration
  if(smoothingFactor == 1){
    landscape_smooth <- landscape
  }else{
    pad <- createPaddedMatrix(landscape, smoothingFactor)
    landscape_smooth <- smooth_pad_terra(pad, smoothingFactor, landscape)
  }
  transDat <- createTransDat(cells, landscape_smooth, betas_l, movePen)
  transDat$num <- 0
  # doParallel routine
  for(k in 1:15){
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
    loc <- data.frame(cell = sample(1:ncell(landscape_smooth), 1, replace = TRUE), ymod = 0, xmod = 0)
    transDat[k,]$num <- transDat[k,]$num + 1
    for(iter in 1:(nsims)){
      loc <- makeDecision(landscape_smooth, betas_l, loc$cell, cells, mods, transDat)
      xmod <- xmod + loc$xmod # one of mods is always zero
      ymod <- ymod + loc$ymod
      currTime <- as.POSIXct(currTime) + lubridate::minutes(1)
      transDat[loc$cell,]$num <- transDat[loc$cell,]$num + 1
      out.dat <- rbind(out.dat, cbind(
        t = currTime,
        cell = loc$cell,
        xMod = xmod, # keeps track of the number of times agent has passed the x-axis
        yMod = ymod  # keeps track of the number of times agent has passed the y-axis
      ))
    }
    # ANALYSIS ---------------------------------------------------------------------
    
    # thin the movement dataset   
    out.dat <- out.dat[seq(1:nrow(out.dat)) %% nThin == 0,]   
    
    # make x and y column
    out.dat$x = out.dat$cell %% ncol
    out.dat$y = out.dat$cell/ncol
    out.dat$t = as.POSIXct(out.dat$t)
    trk <- make_track(as_tibble(out.dat), .x = x,
                      .y = y,
                       .t = t) #%>% random_points() %>% extract_covariates(rast(landscape_smooth))

    # mod <- stats::glm(case_ ~ lyr.1, family = "binomial", data = trk)
    # mod$model$coefficients[2]
    # mod <- trk %>% random_points() %>% extract_covariates(rast(landscape_smooth)) %>% amt::fit_rsf(case_ ~ lyr.1)
    # resample track
    # trkRes <- track_resample(trk, rate = lubridate::hours(1), tolerance = lubridate::minutes(5))
    hist(steps(trk)$sl_)
    stps  <- steps(trk) %>% random_steps(n_control = 30)

    # round the decimal places off
    # clamp values to size
    stps[which(!stps$case_),]$x2_ <- stps[which(!stps$case_),]$x2_ %% ncol # row
    stps[which(!stps$case_),]$y2_ <- stps[which(!stps$case_),]$y2_ %% ncol # row
    
    # extract landscape values
    stps <- stps %>% amt::extract_covariates(rast(landscape_smooth))

    # remove incomplete strata
    # randStps <- randStps %>% amt::remove_incomplete_strata() # might not need this?
    # add integer to avoid sl of zero
    stps$sl_ <- stps$sl_ + 1

    # fit ISSF
    mod <- amt::fit_issf(stps, case_ ~ log(sl_) + sl_ + factor(lyr.1) + strata(step_id_))
    
    k <- k + 1
    
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
    print(exp(mod$model$coefficients[3]))
    metaDat <- rbind(metaDat,cbind(
      rep = i,
      beta1 = unlist(betas_l$'1'),
      slctnCoff = unlist(exp(mod$model$coefficients[3])), # grab LS regression coefficient
      smoothingFctr = unlist(smoothingFactor),
      moransI = unlist(Moran(raster(landscape_smooth))),
      movePenalty = unlist(movePen),
      nThin = unlist(nThin)))
    # if(l %% 100 == 0){
    #   subDomName <- paste("smoothingFactor", smoothingFactor, sep = "-") # smoothing Factor
    #   ssubDomName <- paste("beta1", betas_l$'1', sep = "-") # betas
    #   sssubDomName <- paste("movement-penalty", movePen, sep = "-") # movement penalty
    #   ssssubDomName <- paste("thinning", nThin, sep = "-")
    #   fname <- paste0("-", subDomName, "-", ssubDomName,
    #                   "-", sssubDomName, "-", ssssubDomName, "-movement-data-rep-", i, "-realization-", k, sep = "")
    #   #dir.create(fpath, recursive = TRUE) # create dir
    #   #ame <- paste0("movement-data-rep-", i, "-realization-", k, sep = "")
    #   write.table(out.dat,
    #               file = paste0(path, "/", fname), sep = ",")
    #   
    write.table(data.frame(metaDat), file = paste0(path, "/", "metaData","/", "summaryData-2"), sep = ",") # write summary data
    # }
    
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

trk <- make_track(as_tibble(out.dat), .x = x,
                  .y = y,
                  .t = t)
plot(rast(landscape_smooth), xlim = c(0,ncol-1), ylim = c(0,nrow-1),
     col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2,
                       alpha = NULL))
lines(make_track(out.dat.filter, .x = x, .y =y, .t = t),
      col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))
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


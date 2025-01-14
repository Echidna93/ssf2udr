library(terra)
library(raster)
library(amt)
library(dplyr)
library(tidyr)
library(ggplot2)

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

getNeighborsPad <- function(land, pad_rast, sf, i, j){
  pad <- sf
  l <- 1
  nbrs <- c()
  for(k in 1:sf){
    right <- pad_rast[(pad + i), (pad + j + k)]
    left <- pad_rast[(pad + i), (pad + j - k)]
    up <- pad_rast[(pad + i - k), (pad + j)]
    down <- pad_rast[(pad + i + k),(pad + j)]
    # diagonal right up
    dru <- pad_rast[(pad + i - k), (pad + j + k)]
    # diagonal right down
    drd <- pad_rast[(pad + i + k), (pad + j + k)]
    # diagonal left down
    dld <- pad_rast[(pad + i + k), (pad + j - k)]
    # diagonal left up
    dlu <- pad_rast[(pad + i - k), (pad + j - k)]
    nbrs[l] <- right
    l <- l + 1
    nbrs[l] <- left
    l <- l + 1
    nbrs[l] <- up
    l <- l + 1
    nbrs[l] <- down
    l <- l + 1
    nbrs[l] <- dru
    l <- l + 1
    nbrs[l] <- drd
    l <- l + 1
    nbrs[l] <- dld
    l <- l + 1
    nbrs[l] <- dlu
    l <- l + 1
  }
  nbrs[length(nbrs) + 1] <- pad_rast[(pad + i), (pad + j)]
  nbrs
}

smooth_pad_terra <- function(pad, sf, land){
  smooth <- terra::focal(rast(pad), w = sf, fun = "mean",
                         NAonly = TRUE, padValue = 0)
  smooth <- smooth[(sf + 1):(nrow(smooth) - sf),(sf + 1):(nrow(smooth) - sf)]
  smooth <- matrix(smooth$focal_mean, nrow = nrow(land), ncol = ncol(land))
  smooth[smooth > mean(smooth)] <- 1
  smooth[smooth <= mean(smooth)] <- 0
  smooth
  }

smooth_pad_rast <- function(land, pad_rast, sf){
  smooth_mat <- matrix(0, nrow(land), ncol(land))
  for(i in 1:nrow(land)){
    for(j in 1:ncol(land)){
      smooth_mat[i,j] <- mean(getNeighborsPad(ladn, pad_rast, sf, i, j))
    }
  }
  med <- median(smooth_mat)
  smooth_mat[smooth_mat > med] <- 1
  smooth_mat[smooth_mat <= med] <- 0
  smooth_mat
}

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



getSteps <- function(nrow, ncol, landscape){
  cells <- data.frame(cell = 1:ncell(landscape), stay = 1:ncell(landscape),
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
createTransDat <- function(cells, landscape_smooth, betas_l, move_penalty){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- exp(0 + c((landscape_smooth) * as.integer(betas_l[as.character(landscape_smooth)])))
  transMat$left <- exp(-move_penalty + c((landscape_smooth[cells$left]) * as.integer(betas_l[as.character(landscape_smooth[cells$left])])))
  transMat$up <- exp(-move_penalty + c((landscape_smooth[cells$up]) * as.integer(betas_l[as.character(landscape_smooth[cells$up])])))
  transMat$down <- exp(-move_penalty + c((landscape_smooth[cells$left]) * as.integer(betas_l[as.character(landscape_smooth[cells$down])])))
  transMat$right <- exp(-move_penalty + c((landscape_smooth[cells$right]) * as.integer(betas_l[as.character(landscape_smooth[cells$right])])))
  return(transMat)
}

# CONSTANTS --------------------------------------------------------------------

nrow <- 100
ncol <- 100
betaOne <- c(1, 1.5, 2, 2.5, 3)
movePenalties <- c(0, 0.25, 0.5,1)
thinVals <- c(100, 150, 200) # number of samples to throw out before writing
betas_l <- list('0' = 0,
                '1' = 1)
nsims <- nrow*ncol
startTime <- as.POSIXct("2016-11-07 00:00:00 UTC")
smoothingFactorL <- c(1,3,5)
nreps = 2
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


# SIMULATION -------------------------------------------------------------------
rbenchmark::benchmark(for(i in 1:1){
  landscape <- makeLandscapeMatrix(nrow, ncol, TRUE)
  smoothingFactor <- smoothingFactorL[sample(1:length(smoothingFactorL), 1)]
  movePen <- movePenalties[sample(1:length(movePenalties), 1)]
  betas_l$'1' <- betaOne[sample(1:length(betaOne), 1)]
  nThin <- thinVals[sample(1:length(thinVals), 1)] # grabbing the thinning value for this iteration
  if(smoothingFactor == 1){
    landscape_smooth <- landscape
  }else{
    pad <- createPaddedMatrix(landscape, smoothingFactor)
    landscape_smooth <- smooth_pad_terra(pad, smoothingFactor, landscape)
  }
    transDat <- createTransDat(cells, landscape_smooth, betas_l, movePen)
    transDat$num <- 0  
    for(k in 1:1){
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
     loc <- data.frame(cell = k, ymod = 0, xmod = 0)
     transDat[k,]$num <- transDat[k,]$num + 1
        for(iter in 1:(1000)){
            loc <- makeDecision(landscape_smooth, betas_l, loc$cell, cells, mods, transDat)
            xmod <- xmod + loc$xmod # one of mods is always zero
              ymod <- ymod + loc$ymod
              currTime <- as.POSIXct(currTime) + lubridate::minutes(1)
              transDat[loc$cell,]$num <- transDat[loc$cell,]$num + 1
              sampIter <- sampIter + 1
              out.dat <- rbind(out.dat, cbind(
                t = currTime,
                cell = loc$cell,
                xMod = xmod, # keeps track of the number of times agent has passed the x-axis
                yMod = ymod  # keeps track of the number of times agent has passed the y-axis
              ))
              #print(out.dat)

          }
        # ANALYSIS ---------------------------------------------------------------------
        
        # thin the movement dataset   
        out.dat <- out.dat[seq(1:nrow(out.dat)) %% nThin == 0,]   
     
        # make x and y column
        out.dat$x = ceiling(out.dat$cell/ncol)
        out.dat$y = ifelse(out.dat$cell %% ncol == 0, ncol, out.dat$cell %% ncol)
        out.dat$t = as.POSIXct(out.dat$t)
        trk <- make_track(as_tibble(out.dat), .x = x,
                          .y = y,
                          .t = t)
        # resample track
        trkRes <- track_resample(trk, rate = lubridate::hours(1), tolerance = lubridate::minutes(15))
        stps  <- steps(trkRes) %>% amt::random_steps(n_control = 30)

        # round the decimal places off
        stps$x2_ <- round(stps$x2_)
        stps$y2_ <- round(stps$y2_)

        # clamp values to size
        stps$x2_ <- ceiling(stps$x2_/ncol) # column
        stps$y2_ <- stps$y2_ %% ncol # row

        # extract landscape values
        stps <- stps %>% amt::extract_covariates(rast(landscape_smooth))

        # remove incomplete strata
        # randStps <- randStps %>% amt::remove_incomplete_strata() # might not need this?
        # add integer to avoid sl of zero
        stps$sl_ <- stps$sl_ + 1
        
        # fit ISSF
        mod <- amt::fit_issf(stps, case_ ~ as.factor(lyr.1) + log(sl_) + sl_ + strata(step_id_))

        metaDat <- rbind(metaDat, cbind(
          rep = rep,
          beta1 = betas_l$'1',
          slctnCoff = mod$model$coefficients[1], # grab LS regression coefficient
          smoothingFctr = smoothingFactor,
          moransI = Moran(raster(landscape)),
          movePenalty = movePen,
          nThin = nThin))
    }
      
    
    subDomName <- paste("smoothingFactor", smoothingFactor, sep = "-") # smoothing Factor
    ssubDomName <- paste("beta1", betas_l$'1', sep = "-") # betas
    sssubDomName <- paste("movement-penalty", movePen, sep = "-") # movement penalty
    ssssubDomName <- paste("thinning", nThin, sep = "-")
    fpath <- paste0(path, "/", subDomName, "/", ssubDomName,
                    "/", sssubDomName, "/", ssssubDomName)
    dir.create(fpath, recursive = TRUE) # create dir
    write.table(out.dat,
                file = paste0(fpath, "/", "movement-data"), sep = ",") # write movement data to file
    }, replications = 10) 
write.table(path, "/", metaDat, file = "summaryData", sep = ",") # write summary data



# PLOTTING ---------------------------------------------------------------------

# Below thins the trajectories exclude trajectories that cross the boundary

# create lagged vector for previous statement
out.dat$xModPrev <- data.table::shift(out.dat$xMod, 1, fill = 0)
out.dat$yModPrev <- data.table::shift(out.dat$yMod, 1, fill = 0)
# filter out points where xModPrev doesn't equal xMod
out.dat.filter <- out.dat %>% filter(xModPrev == xMod, yModPrev == yMod)
trk <- make_track(as_tibble(out.dat), .x = x,
                  .y = y,
                  .t = t)
plot(rast(landscape_smooth), xlim = c(1,ncol), ylim = c(1,nrow),
     col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2,
                       alpha = NULL))
lines(make_track(as_tibble(out.dat.filter), .x = x, .y = y, .t = t),
      col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))


# PLOT PATH WITH GGPLOT --------------------------------------------------------
# qplot(x, y, data = out.dat, aes(color="red"))+ 
#   geom_path(aes(color="red")) + scale_x_continuous(limits = c(1,ncol)) + scale_y_continuous(limits = c(1,nrow))

# pivot on the data.frame 
#j <- landscape_smooth %>% as.data.frame() %>% rownames_to_column("Var1") %>% pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>% mutate(
# Var1 = factor(Var1, levels = 1:ncol),
# Var2 = factor(gsub("V", "", Var2), levels = 1:nrow)) 

# plot
# ggplot(j, aes(x, y,)) + geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = "white", high="black") +
#   geom_path(data = out.dat, aes(x = x, y = y, color = "red")) +
#   geom_point(data = out.dat, aes(x = x, y = y, color = "pink"))

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
ggplot(metaDat, aes(y = unlist(slctnCoff), x = unlist(smoothingFctr), color = factor(unlist(beta1)))) + geom_smooth()
# ggplot(out.dat.long, aes(x = factor(smoothingFactor), y = logRSS, color = factor(beta1))) + geom_boxplot()
# 
# 
# write.table(out.dat, "./outdat-sf-1:7-50-50", sep = ",")
# hist(pts$lyr.1)


# # PLOT PATH WITH BASE PLOT ---------------------------------------------------

# plot(rast(landscape_smooth), xlim = c(1,ncol), ylim = c(1,nrow), col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
# lines(make_track(as_tibble(out.dat), .x = x, .y = y, .t = t), col = "red", lwd=2, xlim = c(0,50), ylim=c(0,50))


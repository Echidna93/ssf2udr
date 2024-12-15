# library(tidyverse)
library(plotrix)
library('plot.matrix')

#' initiates a landscape matrix of vectors of random 0's and 1's
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
make_landscape_matrix <- function(numrow, numcol, binary=TRUE){
  x = seq(1, numrow, by = 1)
  y = seq(1, numcol, by = 1)
  coords <- list(expand.grid(x,y))
  vals <- sample(x = c(0,1),size = numrow*numcol, replace = TRUE)
  matrix(vals, ncol = numcol, nrow=numrow)
  
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an make_decision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param 
#' @export
make_decision<-function(landscape,R, l, betas_l, cell, cells, transition_prob){
  selection <- sample(c(1:5), 1, prob = transition_prob[cell,3:7])
  selection <- selection + 1
  cells[cell, selection]
}
#' Chooses best possible landscape component to move to
#' TODO alter in the case of an make_decision being fed an empty list, make 
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

smooth_rast <- function(landscape, smoothingFactor){
  matrix <- landscape
  for(i in 1:length(landscape)){
    nbrs <- getNeighbors(nrow(landscape), ncol(landscape), i, landscape, smoothingFactor)
    matrix[i] <- mean(landscape[nbrs])
  }
  med <- mean(matrix)
  matrix[matrix > med] <- 1
  matrix[matrix <= med] <- 0
  matrix
}


#' initiates a data frame of inds
#' @param n.initial # inds
#' @export
makeInds <- function(nInds, cellInit, dim){
  return(data.frame(1, cell = cellInit, t=1))
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
  smooth <- terra::focal(rast(pad), w = sf, fun = "median", NAonly = TRUE, padValue = 0)
  smooth <- smooth[(sf + 1):(nrow(smooth) - sf),(sf + 1):(nrow(smooth) - sf)]
  smooth <- matrix(smooth$focal_median, nrow = nrow(land), ncol = ncol(land))
  # med <- median(smooth)
  # smooth[smooth > med] <- 1
  # smooth[smooth <= med] <- 0
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

getNeighbors <- function(nrow, ncol, loc, landscape, sl){
  nbrs <- c()
  n <- 1
  for(s in 1:sl){
    nbrs[n] <- n + s*nrow # right
    if(loc + s*nrow > ncol*nrow){
      nbrs[n] <- (loc + s*nrow) %% (ncol*nrow)
      n <- n + 1
    }else{
      nbrs[n] <- loc + s*nrow  
      n <- n + 1
    }
    if(loc - s*nrow < 1){
        nbrs[n] <- ncol*nrow + (loc - s*nrow)
        n <- n + 1
    }else{
      nbrs[n] <- loc - s*nrow  
      n <- n + 1
    }
    
    # up
    if(loc - s < (ceiling(loc/ncol)*ncol - (ncol-1))){
      nbrs[n] <- (((loc - (loc %% ncol)) + 1) + ncol - 1) - (s - (loc %% ncol))
      n <- n + 1
    }else{
      nbrs[n] <- loc - s  
      n <- n + 1
    }
    # down
    if(loc + s > ceiling(loc/ncol)*ncol){
      nbrs[n] <- abs((ceiling(loc/ncol) * ncol) - loc - s) + ((ceiling(loc/ncol)-1) * ncol)
      n <- n + 1
    }else{
      nbrs[n] <- loc + s
      n <- n + 1
    }
    if(loc - (s*ncol + s) < 1){
      nbrs[n] <- (nrow*ncol - ((((s - ceiling(loc/ncol)) * ncol)) + abs(((loc - ((s*ncol + s))))) %% ncol))
      n <- n + 1
      }else{
     nbrs[n] <- loc - (s*ncol + s)
     n <- n + 1
    }
    # if((loc + (s*ncol + s)) > (ncol * nrow)){
    #   nbrs[n] <- (loc + ((s - 1)*ncol + s)) %% ncol*nrow
    #   print((loc + ((s - 1)*ncol + s)) %% ncol*nrow)
    #   n <- n + 1
    # }else{
    #   nbrs[n] <- (loc + (s*ncol + s))
    #   n <- n + 1
    # }
  }
  nbrs[length(nbrs) + 1] <- loc
  return(nbrs)
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


getSteps <- function(nrow, ncol, landscape){
cells <- data.frame(cell = 1:ncell(landscape), stay = 1:ncell(landscape), right = 0, left = 0,
                      up = 0, down = 0)
for(i in 1:nrow(cells)){
  if(i == 1){
    cells[i,]$left <- (ncol*nrow - ncol) + 1
    cells[i,]$up <- ncol
    cells[i,]$down <- i + 1
    cells[i,]$right <- i + ncol
  }
  else if(i == ncol){
    cells[i,]$left <- ncol * nrow
    cells[i,]$up <- i - 1
    cells[i,]$down <- 1
    cells[i,]$right <- i + ncol
  }
  else if(i == ncol*nrow){
    cells[i,]$left <- i - ncol
    cells[i,]$up <- i - 1
    cells[i,]$down <- ncol*nrow - ncol + 1
    cells[i,]$right <- ncol 
  }
  else if(i == (nrow*ncol - ncol) + 1){
    cells[i,]$left <-i - nrow 
    cells[i,]$up <- ncol*nrow
    cells[i,]$down <- i + 1
    cells[i,]$right <- 1
  }
  # bottom row
  
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
createTransDat <- function(cells, landscape_smooth, betas_l){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- exp(0 + c((landscape_smooth) * as.integer(betas_l[as.character(landscape_smooth)])))
  for(i in 1:nrow(transMat)){
    transMat[i,]$left <- exp(-1 + c((landscape_smooth[[cells[i,]$left]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$left]])])))
    transMat[i,]$up <- exp(-1 + c((landscape_smooth[[cells[i,]$up]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$up]])])))
    transMat[i,]$down <- exp(-1 + c((landscape_smooth[[cells[i,]$down]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$down]])])))
    transMat[i,]$right <- exp(-1 + c((landscape_smooth[[cells[i,]$right]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$right]])])))
  }
  return(transMat)
}

nrow <- 50
ncol <- 50
betas_l <- list('0' = -1,
                '1' = 1)
nsims <- nrow*ncol
R <- 1
l <- 1
smoothingFactorL <- c(1,3,5)
nreps = 20
out.dat <- data.frame(matrix(nrow = 0, ncol = 10))
names(out.dat) <- c("rep",
                    "smoothingFactor",
                    "moran",
                    "beta1",
                    "beta0",
                    "beta0",
                    "hb1",
                    "hb0",
                    "b1",
                    "b0")
locs <- data.frame(matrix(0, nrow = 0, ncol = 4))
names(locs) <- c("x", "y", "t")
betaOne <- c(1, 2, 3)
# betaZero <- c(0, -0.1, -0.3)
landscape <- make_landscape_matrix(nrow, ncol, TRUE)
cells <- getSteps(nrow, ncol, landscape)
for(t in 1:nreps){
  landscape <- make_landscape_matrix(nrow, ncol, TRUE)
  #print(t)
  for(p in 1:length(smoothingFactorL)){
    if(p == 1){
      landscape_smooth <- landscape
      smoothingFactor <- 1
    }else{
     smoothingFactor <- smoothingFactorL[p]
     pad <- createPaddedMatrix(landscape, smoothingFactor)
     landscape_smooth <- smooth_pad_terra(pad, smoothingFactor, landscape)
     #landscape_smooth <- smooth_rast(landscape, smoothingFactor)
     # 
     # landscape_smooth <-matrix(focal(rast(landscape), w = 3, mean, pad = TRUE,
     #                                        padValue = NA, na.rm = TRUE, wrap = TRUE),
     #                                    nrow = nrow, ncol = ncol)
     #  med <- median(landscape_smooth)
     #  landscape_smooth[landscape_smooth > med] <- 1
     #  landscape_smooth[landscape_smooth <= med] <- 0
     }
    # print(betas_l)
    for(l in 1:length(betaOne)){
      betas_l$'1' <- betaOne[l]
  #    for(d in 1:length(betaZero)){
  #      betas_l$'0' <- betaZero[d]
      transDat <- createTransDat(cells, landscape_smooth, betas_l)
      transDat$num <- 0
      for(k in 1:nsims){
        loc <- makeInds(1,k,nrow) # start pos
        transDat[k,]$num <- transDat[k,]$num + 1
        # locs <- rbind(locs, new_loc)
        for(iter in 1:(nrow*ncol)){
          loc$cell <- make_decision(landscape_smooth, R, l, betas_l, loc$cell, cells, transDat)
          # loc$t <- new_loc$t + 1
          transDat[loc$cell,]$num <- transDat[loc$cell,]$num + 1
          }
      }
      out.dat <- rbind(out.dat, cbind(
                              rep = t,
                              smoothingFactor,
                              Moran(raster(landscape_smooth)),
                              beta1 = betas_l$'1',
                              beta0 = betas_l$'0',
                              hb1 = sum(transDat[which(landscape_smooth==1),]$num),
                              hb0 = sum(transDat[which(landscape_smooth==0),]$num),
                              b1 = length(which(landscape_smooth == 1)),
                              b0 = length(which(landscape_smooth ==0 ))))
      }
    }
  }
#}
out.dat.long <- out.dat %>% mutate(uid = paste0(rep, "-", smoothingFactor, "-", beta1)) %>%
  mutate(logRSS = log((hb1*b0)/(hb0*b1))) %>%
  group_by(uid) %>% mutate(meanlogRSS = mean(logRSS))

ggplot(out.dat.long, aes(x = smoothingFactor, y = logRSS, color = factor(beta1))) + geom_smooth()
ggplot(out.dat.long, aes(x = factor(smoothingFactor), y = logRSS, color = factor(beta1))) + geom_boxplot()


write.table(out.dat, "./outdat-sf-1:7-50-50", sep = ",")
hist(pts$lyr.1)

plot(rast(landscape_smooth))
points(make_track(as_tibble(locs), .x = x, .y = y, .t = t), col = "white", pch=20)

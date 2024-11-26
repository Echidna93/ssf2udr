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

#' initiates a data frame of inds
#' @param n.initial # inds
#' @export
makeInds <- function(nInds, cellInit, dim){
  return(data.frame(1, cell = cellInit, t=1))
}

getSteps <- function(nrow, ncol, landscape){
cells <- data.frame(cell = 1:ncell(landscape), stay = 1:ncell(landscape), right = 0, left = 0,
                      up = 0, down = 0)
for(i in 1:nrow(cells)){
  if(i == 1){
    cells[i,]$left <- ncol
    cells[i,]$up <- (ncol*nrow - ncol) + 1
    cells[i,]$down <- i + ncol 
    cells[i,]$right <- i + 1
  }
  else if(i == ncol){
    cells[i,]$left <- i - 1
    cells[i,]$up <- ncol * nrow
    cells[i,]$down <- i + ncol
    cells[i,]$right <- 1
  }
  else if(i == ncol*nrow){
    cells[i,]$left <- i - 1
    cells[i,]$up <- i - ncol
    cells[i,]$down <- ncol
    cells[i,]$right <- ncol*nrow - ncol + 1
  }
  else if(i == (nrow*ncol - ncol) + 1){
    cells[i,]$left <- ncol*nrow
    cells[i,]$up <- i - nrow
    cells[i,]$down <- 1
    cells[i,]$right <- i + 1
  }
  # bottom row
  else if((i > (nrow * ncol - ncol) + 1) & (i < nrow*ncol)){
    cells[i,]$left <- i - 1
    cells[i,]$up <- i - nrow
    cells[i,]$down <- i %% nrow
    cells[i,]$right <- i + 1
  }
  # right edge
  else if(i %% nrow == 0 & (i != ncol & i != ncol*nrow)){
    cells[i,]$left <- i - 1
    cells[i,]$up <- i - nrow
    cells[i,]$down <- i + nrow
    cells[i,]$right <- i - ncol + 1
  }
  # left edge
  else if((i %% nrow == 1) & (i != 1 & i != (nrow*ncol - ncol) + 1)){
    cells[i,]$left <- i + (ncol - 1) 
    cells[i,]$up <- i - ncol
    cells[i,]$down <- i + nrow
    cells[i,]$right <- i + 1
  }
  # top edge
  else if((i > 1) & (i < ncol)){
    cells[i,]$left <- i - 1
    cells[i,]$up <- (nrow*ncol - nrow) + (i %% nrow)
    cells[i,]$down <- i + nrow
    cells[i,]$right <- i + 1
  }
  else{
    cells[i,]$left <- i -1
    cells[i,]$up <- i - nrow 
    cells[i,]$down <- (i + nrow)
    cells[i,]$right <- i + 1
    }
  }
  return(cells)
}
createTransDat <- function(cells, landscape_smooth){
  transMat <- data.frame(cell = 1:ncell(landscape),num = 0, stay = 0, right = 0, left = 0,
                         up = 0, down = 0)
  transMat$stay <- exp(0 + c(exp((landscape_smooth) * as.integer(betas_l[as.character(landscape_smooth)]))))
  for(i in 1:nrow(transMat)){
    transMat[i,]$left <- exp(-1 + c(exp((landscape_smooth[[cells[i,]$left]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$left]])]))))
    transMat[i,]$up <- exp(-1 + c(exp((landscape_smooth[[cells[i,]$up]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$up]])]))))
    transMat[i,]$down <- exp(-1 + c(exp((landscape_smooth[[cells[i,]$down]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$down]])]))))
    transMat[i,]$right <- exp(-1+ c(exp((landscape_smooth[[cells[i,]$right]]) * as.integer(betas_l[as.character(landscape_smooth[[cells[i,]$right]])]))))
  }
  return(transMat)
}

nrow <- 20
ncol <- 20
betas_l <- list('0' = -0.1,
                '1' = 1)
landscape <- make_landscape_matrix(nrow, ncol, TRUE)
nsims <- nrow*ncol
R <- 1
l <- 1
smoothingFactorL <- c(1, 3, 5)
nreps = 10
out.dat <- data.frame(matrix(nrow = 0, ncol = 8))
names(out.dat) <- c("rep",
                    "smoothingFactor",
                    "beta1",
                    "beta0",
                    "hb1",
                    "hb0",
                    "b1",
                    "b0")
locs <- data.frame(matrix(0, nrow = 0, ncol = 4))
names(locs) <- c("x", "y", "t")
betaOne <- c(1, 2, 3)
cells <- getSteps(nrow, ncol, landscape)
for(t in 1:nreps){
  print(t)
  for(p in 1:length(smoothingFactorL)){
    if(p == 1){
      landscape_smooth <- landscape
      smoothingFactor <- 1
    }else{
     smoothingFactor <- smoothingFactorL[p]
     landscape_smooth <-matrix(focal(rast(landscape), w = smoothingFactor, mean, pad = TRUE,
                                            padValue = NA, na.rm = TRUE, wrap = TRUE),
                                        nrow = nrow, ncol = ncol)
      med <- median(landscape_smooth)
      landscape_smooth[landscape_smooth > med] <- 1
      landscape_smooth[landscape_smooth <= med] <- 0
     }
    # print(betas_l)
    for(l in 1:length(betaOne)){
      transDat$num <- 0
      betas_l$'1' <- betaOne[l]
      transDat <- createTransDat(cells, landscape_smooth)
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
                              beta1 = betas_l$'1',
                              beta0 = betas_l$'0',
                              hb1 = transDat[which(landscape_smooth==1),]$num,
                              hb0 = transDat[which(landscape_smooth==0),]$num,
                              b1 = length(which(landscape_smooth == 1)),
                              b0 = length(which(landscape_smooth ==0 ))))
    }
}
  write.table(matrix(transDat$num, nrow = nrow, ncol = ncol),
              file = paste0("./", "sf", smoothingFactor,"_", "size", nrow, "-", ncol), sep = ",")
}
write.table(out.dat, "./outdat-sf-1:7-50-50", sep = ",")
hist(pts$lyr.1)

plot(rast(landscape_smooth))
points(make_track(as_tibble(locs), .x = x, .y = y, .t = t), col = "white", pch=20)

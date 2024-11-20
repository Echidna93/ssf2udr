# library(tidyverse)
library(plotrix)
library('plot.matrix')

betas_l <- list('0' = 0, 
                '1' = 1)

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

# helper function
# divides by max value of matrix
divide_by_max<-function(x, max){
  x / max
}
#' Updates individual locations
#' @param data_frame holds data about inds
#' @export
move<-function(data_frame, landscape, nrow, ncol, R, l, betas_l){
  nbrs<-get_neighbors(loc = c(data_frame$x,data_frame$y), nrow, ncol, landscape)
  new_loc<-make_decision(landscape=landscape, nbrs=nbrs, R, l, betas_l)
  return(data.frame(x = new_loc[[1]], y = new_loc[[2]], t = data_frame[1,]$t + 1))
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an make_decision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param 
#' @export
make_decision<-function(landscape, nbrs, R, l, betas_l){
  # Random movement
  # indx <- round(runif(1, 1, length(nbrs)))
  # decision_vec<-c(nbrs[[indx]][1],nbrs[[indx]][2])
  # #printnbrs)
  probs <- mapply(function(x, y) exp(-y) * exp(betas_l[as.character(landscape[x[1],][x[2]])][[1]]), nbrs, c(0, 1, 1, 1, 1))
  probs <- probs / sum(probs)
  # print(probs)
  selection <- rmultinom(1, c(1:5), prob = probs)
  indx <- which(selection==1)
  decision_vec<-c(nbrs[[indx]][1],nbrs[[indx]][2])
  decision_vec
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
  steps <- expand.grid(steps, steps) %>% as.data.frame() %>% dplyr::rename(x = Var1, y = Var2) %>%
    filter(x != y) %>% filter(x != -y)
  steps <- rbind(c(0,0),steps)
  for(i in 1:nrow(steps)){
    for(k in 1:length(smooth)){
      steps <- rbind(steps, smooth[k] * steps[i,])
    }
  }
  steps <- steps %>% unique()
}

#' Chooses best possible landscape component to move to
#' TODO alter in the case of an make_decision being fed an empty list, make 
#' else case
#' TODO implement sorting function
#' @param landscape matrix to smooth
#' @param smoothingFactor size of moving window to smooth
#' @export
smooth_rast <- function(landscape, smoothingFactor, ncol, nrow){
  matrix <- landscape
  steps <- generateSteps(smoothingFactor)
  l <- replicate(nrow(steps), c(0,0), simplify = FALSE)
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      loc <- c(i,j)
      for(k in 1:nrow(steps)){
        if(!((loc[1] + steps[k,]$x) > ncol & (loc[1] + steps[k,]$y) < 1)){ # right edge
          l[[k]][1] <- loc[1] + steps[k,]$x
        }
        if(!((loc[2] + steps[k,]$y) < 1 & ((loc[2] + steps[k,]$y) > nrow))){
          l[[k]][2] <-  loc[2] + steps[k,]$y
        }
        if((loc[1] + steps[k,]$x) < 1){
          l[[k]][1] <- ncol + (steps[k,]$x - 1)
        }
        if((loc[1] + steps[k,]$x) > ncol){ # right edge
          l[[k]][1] <- 1 + (steps[k,]$x - 1)
        }
        if((loc[2] + steps[k,]$y) < 1){
          l[[k]][2] <-  nrow + (steps[k,]$y - 1)
        }
        if((loc[2] + steps[k,]$y) > nrow){
          l[[k]][2] <- 1 + (steps[k,]$y - 1)
        }
      }
      nbr_vals_ <- lapply(l, function(x) landscape[x[1],][x[2]])
      if(any(is.na(nbr_vals_))){
        #printl)
      }
      matrix[i,j] <- round(mean(unlist(nbr_vals_)))
    }
  }
  matrix
}

makeIndiceList <- function(nrow, ncol){
  as.data.frame(expand.grid(seq(1:nrow), seq(1:ncol))) %>% rename(x = Var1, y = Var2)
}

#' initiates a data frame of inds
#' @param n.initial # inds
#' @export
makeInds <- function(nInds, xInit, yInit, dim){
  return(data.frame(x=xInit, y=yInit, t=1))
}

#' Helper function
#' returns neighborhood of cells as a list of vectors
#' @param loc current coordinates of individual
#' @param nrow # of rows in landscape matrix
#' @param ncol # columns in landscape matrix
#' @export
get_neighbors<-function(loc, nrow, ncol, landscape){
  l<-list(c(0,0),
          c(0,0),
          c(0,0),
          c(0,0),
          c(0,0))
  # check if either x,y element of loc is greater than
  # the dimension of the landscape matrix
  # case 1 on left edge
  steps <- list(c(x = 0, y = 0),
                c(x = 1, y = 0),
                c(x = 0, y = 1),
                c(x = 0, y = -1),
                c(x = -1, y = 0))
  for(i in 1:length(steps)){
    if(!((loc[1] + steps[[i]][1]) > ncol) & !((loc[1] + steps[[i]][1]) < 1)){ # right edge
      l[[i]][1] <- loc[1] + steps[[i]][1]
    }
    if(!((loc[2] + steps[[i]][2]) < 1) & !((loc[2] + steps[[i]][2]) > nrow)){
      l[[i]][2] <-  loc[2] + steps[[i]][2]
    }
    if((loc[1] + steps[[i]][1]) < 1){
      l[[i]][1] <- ncol
    }
    if((loc[1] + steps[[i]][1]) > ncol){ # right edge
      l[[i]][1] <- 1
    }
    if((loc[2] + steps[[i]][2]) < 1){
      l[[i]][2] <-  nrow
    }
    if((loc[2] + steps[[i]][2]) > nrow){
      l[[i]][2] <- 1
    }
  }
  #printl)
  return(l)
  
}
# Main simulation
# TODO make main simulation loop, using functional? (lapply)
check_inds_locs<-function(ind, data_frame){
  for(i in 1:nrow(data_frame)){
    if(is_same_location(ind, data_frame[i,]) & is_not_same_id(ind, data_frame[i,])){
      
    }
  }
}

#helper function is location the same
is_same_location<-function(ind1, ind2){
  (ind1$xloc == ind1$yloc) & (ind1$yloc == ind2$yloc)
}

runSim <- function(nrow, ncol, landscape_smooth, betas_l, out.dat,
                     R,l, smoothingFactor){

}


#helper checks if two ids are the same
is_not_same_id<-function(ind1, ind2){
  !(ind1$id == ind2$id)
}
# CONSTANTS
# TODO look into making these changeable by a user in an x11() window?
# calling the random number generator beforehand
# long format; sparse matrix
# generate random numbers a priori
# clear environment every few steps
nrow <- 10
ncol <- 10
landscape <- make_landscape_matrix(nrow, ncol, TRUE)
nsims <- nrow*ncol
R <- 1
l <- 1
smoothingFactorL <- c(1:7)
nreps = 10
out.dat <- data.frame(matrix(nrow = 0, ncol = 3))
names(out.dat) <- c("smoothingFactor",
                    "meanIntsty1",
                    "meanIntsty0")
locs <- data.frame(matrix(0, nrow = 0, ncol = 3))
names(locs) <- c("x", "y", "t")
indxList <- makeIndiceList(nrow, ncol)
for(i in 1:length(smoothingFactorL)){
  smoothingFactor <- smoothingFactorL[i]
  landscape_smooth <- smooth_rast(landscape, smoothingFactor, nrow, ncol)
  for(j in 1:nreps){
    intensity_mat <- matrix(0, nrow = nrow, ncol = ncol)
    for(k in 1:nsims){
      new_loc <- makeInds(1,indxList[k,]$x,indxList[k,]$y,nrow) # start pos
      intensity_mat[new_loc[1,]$x, new_loc[1,]$y] <- intensity_mat[new_loc[1,]$x, new_loc[1,]$y] + 1    
      # locs <- rbind(locs, new_loc)
      for(iter in 1:(nrow*ncol)){
        new_loc <- move(new_loc, landscape_smooth, nrow, ncol, R, l, betas_l)
        # locs <- rbind(locs, new_loc)
        intensity_mat[new_loc[1,]$x, new_loc[1,]$y] <- intensity_mat[new_loc[1,]$x, new_loc[1,]$y] + 1    
        }
      }
    out.dat <- rbind(out.dat, cbind(smoothingFactor,
                            meanIntsty1 = mean(intensity_mat[which(landscape_smooth==1)]/((nsims)*(nrow*ncol))),
                            meanIntsty0 = mean(intensity_mat[which(landscape_smooth==0)])/((nsims)*(nrow*ncol))))
  }
}
write.table(out.dat, "./outdat-sf-4-50-50", sep = ",")
hist(pts$lyr.1)

plot(rast(landscape_smooth))
points(make_track(as_tibble(locs), .x = x, .y = y, .t = t), col = "white", pch=20)

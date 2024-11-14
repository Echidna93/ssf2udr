# library(tidyverse)
library(plotrix)
library('plot.matrix')

betas_l <- list('0' = 0, 
                '1' = 1)

#' initiates a landscape matrix of vectors of random 0's and 1's
#' @param nrow number of rows in matrix
#' @param ncol number of columns in matrix
#' @export
make_landscape_matrix <- function(numrow, numcol, R, binary=TRUE){
    x = seq(1, numrow, by = 1)
    y = seq(1, numcol, by = 1)
    coords <- list(expand.grid(x,y))
    vals <- sample(x = c(0,1),size = numrow*numcol, replace = TRUE)
    vals = unlist(sapply(coords, function(x) round(cos((pi * x[1])/R) * cos((pi * x[2])/R))))
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
    nbrs<-get_neighbors(loc = c(data_frame[1,]$x,data_frame[1,]$y), nrow, ncol, landscape)
    #print(nbrs)
    # new_loc<-nbrs[[round(runif(1,1,length(nbrs)))]]
    new_loc<-make_decision(landscape=landscape, nbrs=nbrs, R, l, betas_l)
    #new_loc
    data_frame[1,]$x<-new_loc[[1]]
    data_frame[1,]$y<-new_loc[[2]]
    data_frame[1,]$t <- data_frame[1,]$t + 1
    data_frame
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
  # print(nbrs)
  probs <- mapply(function(x, y) exp(-y)*exp(-2*betas_l[as.character(landscape[x[1],][x[2]])][[1]]), nbrs, c(0, 1, 1, 1, 1))
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
        print(l)
      }
      matrix[i,j] <- round(mean(unlist(nbr_vals_)))
      }
   }
matrix
}

#' initiates a data frame of inds
#' @param n.initial # inds
#' @param dim dimension of a vector for random assignment of location
#' @param nI # infected individuals to start
#' @export
make_inds <- function(n.initial,dim){
  id<-1:n.initial
  x<-round(runif(1, min=0, max=dim))
  y<-round(runif(1, min=0, max=dim))
  inds <- data.frame(id = id, x=x, y=y,t=1, stringsAsFactors=FALSE) 
  inds
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
      print(l)
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

#helper checks if two ids are the same
is_not_same_id<-function(ind1, ind2){
  !(ind1$id == ind2$id)
}

# CONSTANTS
# TODO look into making these changeable by a user in an x11() window?
numrow=100
numcol=100
R <- exp(1*1 + 0 * 0)
landscape<-make_landscape_matrix(numrow, numcol, R, TRUE)
landscape_smooth <- smooth_rast(landscape, 2, 100, 100)
landscape_smooth3 <- smooth_rast(landscape, 3, 100, 100)
nsims <- 100
for(k in 1:nsims){
  l = 1
  start<-make_inds(1,numrow)
  locs_smooth <- start
  nsims <- 100
  for(iter in 1:1000){
    if(iter == 1){
      new_loc<-move(start, landscape_smooth, numrow, numcol, R, l, betas_l)  
    }else{
      old_loc <- new_loc
      new_loc <- move(old_loc, landscape_smooth, numrow, numcol,R, l, betas_l)
    }
    locs_smooth<-rbind(locs_smooth, new_loc)
  }
  ls_rast <- rast(landscape)
  pts <- terra::extract(ls_rast, data.frame(x = locs$x, y = locs$y))
  prop.1[k] <- length(which(pts ==1))/nrow(pts)
}


hist(pts$lyr.1)

plot(rast(landscape))
points(make_track(as_tibble(locs_smooth), .x = x, .y = y, .t = t), col = "white", pch=20)

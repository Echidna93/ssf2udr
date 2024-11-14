# use the amt data
fisher_dat <- amt::amt_fisher
nrows = 100
ncols = 100
# generate map
land <- terra::rast(nrows=nrows, ncols=ncols, xmin=min(fisher_dat$x_),
                    xmax=max(fisher_dat$x_), ymin=min(fisher_dat$y_), ymax = max(fisher_dat$y_))
make_landscape_matrix <- function(numrow, numcol, binary=TRUE){
    matrix(sample(c(1,2,3), replace=TRUE, size=numrow*numcol), ncol = numcol, nrow=numrow)
}
mat <- make_landscape_matrix(100, 100) # make landscape matrix
land <- rast(mat)

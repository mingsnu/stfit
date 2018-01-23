###########
library(raster)
library(rasterVis)
files = list.files("../data/2014", "*.tif$", full.names = TRUE)
## smoothing spline vs fourier basis spline regression
r = brick(stack(files[1:10]))
levelplot(r, par.settings = RdBuTheme)
mat = t(values(r)) ## note that the image value is stacked row-wise
wmat = matrix(c(0.375, 0.562, 0.375, 0.562, 0.75, 0.562, 0.375, 0.562, 0.375), 3)
tmp = mean_est(mat, 31, 31, wmat)
image(matrix(tmp, 31))
levelplot(raster(matrix(tmp, 31)))

#### Validating mean_est function
mat = t(matrix(1:120, 30)) # 4 images, each is 5x6, row stacking
matrix(apply(mat,2, mean), byrow = TRUE, nrow = 5)
tmp = mean_est(mat, 5, 6, wmat)
resmat = matrix(tmp, byrow=TRUE, nrow=5)


## check
nrow = 5; ncol=6
## (1,1)
sum(t(mat[,c(1,2, c(1,2)+ncol)]) * c(t(wmat[2:3, 2:3]))) / sum(c(t(wmat[2:3, 2:3]))) / 4
## (1,2)
sum(t(mat[,c(1,2,3, c(1:3)+ncol)]) * c(t(wmat[2:3, 1:3]))) / sum(c(t(wmat[2:3, 1:3]))) / 4
## (2,1)
sum(t(mat[,c(1,2,c(1,2)+ncol, c(1,2)+2*ncol)]) * c(t(wmat[1:3, 2:3]))) / sum(c(t(wmat[1:3, 2:3]))) / 4
## (2,2)
sum(t(mat[,c(1,2,3, c(1:3)+ncol, c(1:3)+2*ncol)]) * c(t(wmat))) / sum(c(t(wmat))) / 4 == resmat[2,2]

#### Validating mean_est function on data with NA
mat = t(matrix(1:120, 30))
mat[1,1] = NA
tmp = mean_est(mat, 5, 6, wmat)
matrix(tmp, byrow=TRUE, nrow=5)

## check
## (1,1)
sum(t(mat[,c(1,2, c(1,2)+ncol)]) * c(t(wmat[2:3, 2:3])), na.rm = TRUE) / (sum(c(t(wmat[2:3, 2:3]))) *4 - wmat[2,2])
## (1,2)
sum(t(mat[,c(1,2,3, c(1:3)+ncol)]) * c(t(wmat[2:3, 1:3])), na.rm = TRUE) / (sum(c(t(wmat[2:3, 1:3]))) *4 - wmat[2,1])
## (2,1)
sum(t(mat[,c(1,2,c(1,2)+ncol, c(1,2)+2*ncol)]) * c(t(wmat[1:3, 2:3])), na.rm = TRUE) / (sum(c(t(wmat[1:3, 2:3])))*4 - wmat[1,2])
## (2,2)
sum(t(mat[,c(1,2,3, c(1:3)+ncol, c(1:3)+2*ncol)]) * c(t(wmat)), na.rm = TRUE) / (sum(c(t(wmat)))*4 - wmat[1,1])



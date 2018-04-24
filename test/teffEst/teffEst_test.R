##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)


###############################
#### I. Landsat data test #####
###############################
df = read_feather("../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

###################################
#### 1. Overall mean estimaton ####
###################################
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest = meanEst(doy, mat, doyeval = 1:365)
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 31)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm)

###################################
#### 2. Time effect estimation ####
###################################
## remove outlier pixels
for(i in 1:length(meanest$outlier$outidx)){
  mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
}
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

teffarray = teffEst(year, doy, rmat,
                    doyeval = meanest$doyeval,
                    h.cov = 100, h.sigma2 = 300)

### visualization 1
op <- par(ask=TRUE)
yeareval = as.numeric(dimnames(teffarray)[[1]])
doyeval = as.numeric(dimnames(teffarray)[[2]])
for(i in 1:length(yeareval)){
  ind = which(year==yeareval[i])
  plot(doy[ind], rmat[ind,1], col = i)
  lines(doyeval, teffarray[i,,1], col = i)
}
par(op)
### visualization 2
op <- par(ask=FALSE)
yeareval = as.numeric(dimnames(teffarray)[[1]])
doyeval = as.numeric(dimnames(teffarray)[[2]])
plot(0, xlim = range(doyeval), ylim = range(rmat, na.rm = TRUE), type = "n")
for(i in 1:length(yeareval)){
  ind = which(year==yeareval[i])
  points(doy[ind], rmat[ind,1], col = i, pch=19)
  lines(doyeval, teffarray[i,,1], col = i, lwd=3)
  locator(1)
}
par(op)

### visualization 3
## sequence of images for one specific year
teffmat_stack = mat2stack(teffmat[12,,], 31)
levelplot(teffmat_stack[[seq(1, 100, by = 5)]], par.settings = colthm)



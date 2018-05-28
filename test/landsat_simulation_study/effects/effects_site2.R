##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

df = read_feather("../../data/features_2_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 3000] = NA

registerDoParallel(6)
stfit::opts$set(temporal_mean_est = spreg)
mat = mat0
#### Detect and remove outliers by first fitting a mean estimation ======
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                  clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
## remove outlier pixels
for(i in 1:length(meanest$outlier$outidx)){
  mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
}
## remove outlier images
outlier.img.idx = meanest$idx$idx.outlier
for(i in outlier.img.idx){
  mat[outlier.img.idx,] = NA
}

#########################
#### mean estimation ####
#########################
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                  clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
## saveRDS(meanest, "effects_output/meanest.rds")
mat_mean_imp = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

############################
#### teffect estimation ####
############################
## calculate the residuals
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
teffarray = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = 100, h.sigma2 = 300)
## saveRDS(teffarray, "effects_output/teffarray.rds")

#### RMSE using mean+teffect
mat_teff_imp = mat_mean_imp
yearidx = unlist(lapply(year, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
doyidx = unlist(lapply(doy, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
for(i in 1:nrow(mat_teff_imp)){
  mat_teff_imp[i,] = mat_teff_imp[i,] + teffarray[yearidx[i], doyidx[i],]
}

############################
#### seffect estimation ####
############################
yearidx = unlist(lapply(year, function(x, y)
  which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
doyidx = unlist(lapply(doy, function(x, y)
  which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
for (i in 1:nrow(rmat)) {
  rmat[i, ] = rmat[i, ] - teffarray[yearidx[i], doyidx[i], ]
}
seffmat = seffEst(rmat, 31, 31, nnr = 30, h.cov = 2, h.sigma2 = 2)$seffmat
mat_seff_imp = mat_teff_imp + seffmat
## saveRDS(seffmat, "effects_output/seffmat.rds")

##################
##### RMSE #######
##################
RMSE(mat0, mat_mean_imp, na.rm=TRUE)
## [1] 121.072
RMSE(mat0, mat_teff_imp, na.rm=TRUE)
## 112.8716
RMSE(mat0, mat_seff_imp, na.rm=TRUE)
## 91.84515


RMSE(mat, mat_mean_imp, na.rm=TRUE)
## [1] 74.98979
RMSE(mat, mat_teff_imp, na.rm=TRUE)
## 61.66152
RMSE(mat, mat_seff_imp, na.rm=TRUE)
## 25.29047

####################################
##### NO teffect; ONLY seffect #####
####################################
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
seffmat1 = seffEst(rmat, 31, 31, nnr = 30, h.cov = 2, h.sigma2 = 2)$seffmat
mat_seff_imp1 = mat_mean_imp + seffmat
RMSE(mat0, mat_seff_imp1, na.rm=TRUE)
## 99.29852
RMSE(mat, mat_seff_imp1, na.rm=TRUE)
## 44.44821



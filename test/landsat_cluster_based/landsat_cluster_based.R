##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA
hist(mat)

colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)
idx = doy >=1 & doy <=365
# idx = doy >=100 & doy <=200

#####################################
#### visualize the original data ####
#####################################
ss = mat2stack(mat, 31, seq(1,100, by=5))
levelplot(ss, par.settings = colthm)

##########################################################
#### gapfill based on  pixel level temporal smoothing ####
##########################################################
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
res1 = gapfill(year[idx], doy[idx], mat[idx,], 31,31, h = 0, nnr=5)
# res1 = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=5)
## temporal trend visulization
ssmean1 = mat2stack(res1$temporal.mean, 31)
levelplot(ssmean1[[seq(1, 100, by = 5)]], par.settings = colthm)

###############################
##### clustering analysis #####
###############################
cmat = values(t(ssmean1[[seq(100, 200, by = 10)]]))
cres = kmeans(cmat, centers = 10, nstart = 10)
levelplot(raster(matrix(cres$cluster, 31)), margin=FALSE)

###########################################################
#### gapfill based on cluster level temporal smoothing ####
###########################################################
Gapfill::opts$set(temporal_mean_est = lpreg)
res2 = gapfill(year[idx], doy[idx], mat[idx,], 31,31, h = 0, nnr=5, cluster = cres$cluster)
#res2 = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=5, cluster = cres$cluster)

## temporal trend visulization
ssmean2 = mat2stack(res2$temporal.mean, 31)
levelplot(ssmean2[[seq(1, 100, by = 5)]], par.settings = colthm)


##########################################
#### Imputed partial missing images ######
##########################################
res = res1
pdf("../output/landsat_partial_nnr5_pixel_wise_temporal_smoothing.pdf")
for(i in res$idx$idx.partialmissing){
  r1 = raster(matrix(mat[idx,][i,], 31))
  r2 = raster(matrix(res$imputed.mat[i,], 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = colthm))
  # Sys.sleep(2)
}
dev.off()
pdf("../output/landsat_partial_resid_nnr5_pixel_wise_temporal_smoothing.pdf")
for(i in res$idx$idx.partialmissing){
  m = res$temporal.mean[which(seq(min(doy[idx]), max(doy[idx])) == doy[idx][i]),]
  r1 = raster(matrix(mat[idx,][i,] - m, 31))
  r2 = raster(matrix(res$imputed.mat[i,] - m, 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = colthm))
  # Sys.sleep(2)
}
dev.off()

res = res2
pdf("../output/landsat_partial_nnr5_cluster_wise_temporal_smoothing.pdf")
for(i in res$idx$idx.partialmissing){
  r1 = raster(matrix(mat[idx,][i,], 31))
  r2 = raster(matrix(res$imputed.mat[i,], 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = colthm))
  # Sys.sleep(2)
}
dev.off()
pdf("../output/landsat_partial_resid_nnr5_cluster_wise_temporal_smoothing.pdf")
for(i in res$idx$idx.partialmissing){
  m = res$temporal.mean[which(seq(min(doy[idx]), max(doy[idx])) == doy[idx][i]),]
  r1 = raster(matrix(mat[idx,][i,] - m, 31))
  r2 = raster(matrix(res$imputed.mat[i,] - m, 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = colthm))
  # Sys.sleep(2)
}
dev.off()
##########################
#### Outlier images ######
##########################
for(i in res$idx$idx.outlier){
  r1 = raster(matrix(mat[idx,][i,], 31))
  r2 = raster(matrix(res$imputed.mat[i,], 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = colthm))
}

############################
##### Simulation study #####
############################
res = res1
df = read_feather("../data/features_106_wide.feather")
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
set.seed(20180124)
n = 6
fidx = sample(res$idx$idx.fullyobserved, n) ## full observed image index
pidx = sample(res$idx$idx.partialmissing, n) ## partial observed image index

r.list = list()
for(i in 1:n){
  r.list[[i]] = raster(matrix(mat[fidx[i], ], 31))
}
## apply missing patterns to fully observed images
for(i in 1:n){
  mat[fidx[i],][is.na(mat[pidx[i],])] = NA
}
for(i in 1:n){
  r.list[[i+n]] = raster(matrix(mat[fidx[i], ], 31))
}
s = stack(r.list)
levelplot(s, par.settings = colthm)

MSE = matrix(0, 5, n)
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
Gapfill::opts$set(temporal_mean_est = llreg)
Gapfill::opts$set(temporal_mean_est = lpreg)
for(k in 1:5){
  ## Gapfill
  res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "lc")
  ## res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k)
  
  ## imputed images
  for(i in 1:n){
    r.list[[i+2*n]] = raster(matrix(res$imputed.mat[fidx[i],], 31))
  }
  s1 = stack(r.list)
  # levelplot(s1, par.settings = RdBuTheme, layout=c(3,3))
  
  ## Calculate MSE
  mat1 = values(s1)
  MSE[k,] = apply((mat1[,1:n] - mat1[,1:n+2*n])^2, 2, sum) / apply(mat1[,1:n+n], 2, function(x) sum(is.na(x)))
}
MSE
## smooth_spline
# [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
# [1,] 20681.42 18325.96 3070.127  7323.574 22642.91 38805.01
# [2,] 19504.34 16617.45 3660.064  6972.476 22503.04 33981.75
# [3,] 17182.69 15309.92 3672.560  7505.769 21400.69 25220.93
# [4,] 13656.42 21943.22 3689.017  7132.400 24719.33 17279.50
# [5,] 23519.12 25722.28 5217.264 10329.878 24733.64 26478.23
## llreg
# [,1]      [,2]     [,3]      [,4]     [,5]     [,6]
# [1,] 51802.88  54273.53 10870.61  9565.664 29231.35 225468.1
# [2,] 43594.40  53876.50 10748.46  9378.665 29231.35 176766.8
# [3,] 38085.79  61828.33 10452.57  9136.471 29231.35 113454.9
# [4,] 36313.46  66794.68 10884.84 11405.680 29231.35 121261.6
# [5,] 33762.50 120534.08 11015.85 10315.146 29231.35 142574.7
## lpreg
# [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
# [1,] 16453.52 15971.24 6374.898  9430.331 27905.23 41756.94
# [2,] 16850.97 15398.08 6582.681  9118.598 28576.34 33616.25
# [3,] 19420.43 17224.70 6515.007  9353.449 24157.97 21228.70
# [4,] 13630.67 20450.86 6752.756  8707.836 27580.86 15351.84
# [5,] 19949.50 21488.99 6845.936 10856.245 24460.00 25001.04
### MSE using "lc" method, loess
MSE = matrix(0, 5, 3)
Gapfill::opts$set(temporal_mean_est = llreg)
Gapfill::opts$set(temporal_mean_est = lpreg)
for(k in 1:5){
  ## Gapfill
  res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "lc", cluster = cres$cluster)
  ## res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k)
  
  ## imputed images
  for(i in 1:n){
    r.list[[i+2*n]] = raster(matrix(res$imputed.mat[fidx[i],], 31))
  }
  s1 = stack(r.list)
  # levelplot(s1, par.settings = RdBuTheme, layout=c(3,3))
  
  ## Calculate MSE
  mat1 = values(s1)
  MSE[k,] = apply((mat1[,1:n] - mat1[,1:n+2*n])^2, 2, sum) / apply(mat1[,1:n+n], 2, function(x) sum(is.na(x)))
}
MSE
#################################################
######### MSE using year >= 2000 data ###########
#################################################
### MSE using "lc" method, smooth_spline
# [,1]     [,2]     [,3]
# [1,] 19572.41 14424.70 3648.579
# [2,] 17742.11 14815.96 3834.272
# [3,] 17701.20 16178.37 3769.953
# [4,] 12887.35 13593.87 4760.713
# [5,] 20251.94 13413.21 4224.410
### MSE using "lc" method, loess
# [,1]     [,2]     [,3]
# [1,] 11170.162 14737.80 5050.032
# [2,] 11358.858 15714.91 5347.452
# [3,] 13660.099 17226.81 4853.116
# [4,]  9989.565 14173.28 5070.335
# [5,] 19086.853 13555.67 5033.158
### MSE using "lc" method, llreg
# [,1]     [,2]     [,3]
# [1,] 15697.920 39660.83 5284.211
# [2,] 13930.911 36024.50 5680.658
# [3,] 12298.597 32968.12 5212.175
# [4,]  9860.581 27667.26 5600.273
# [5,] 11381.210 24923.09 5286.580
#################################################
######### MSE using year >= 2010 data ###########
#################################################
### MSE using "lc" method, meanCurve
# [,1]     [,2]     [,3]
# [1,] 13694.649 11645.54 53027.74
# [2,] 12535.964 13725.02 51095.97
# [3,] 10471.302 12854.14 48750.09
# [4,]  8758.225 12741.40 45267.84
# [5,]  9737.190 10047.46 43043.83
### MSE using "lc" method, loess
# [,1]      [,2]     [,3]
# [1,] 14409.435 10308.731 27758.80
# [2,] 13322.852 11857.228 27470.17
# [3,] 10486.803 10911.912 26615.26
# [4,]  7624.877 11191.374 26120.23
# [5,]  7734.626  8776.976 26648.26

res = res2
MSE.mean = rep(0,n)
for(i in 1:n){
  idx = is.na(mat1[,i+n])
  MSE.mean[i] = sum((mat1[idx,i] - res$temporal.mean[doy[fidx[i]],idx])^2) / sum(idx)
}
MSE.mean
# [1] 32966.95 29772.31 40613.33

## add temporal mean for comparison
for(i in 1:n){
  r.list[[i+3*n]] = raster(matrix(res$temporal.mean[doy[fidx[i]],], 31))
}
s2 = stack(r.list)
levelplot(s2, par.settings = RdBuTheme, layout=c(3,4))

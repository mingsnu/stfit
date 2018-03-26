##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
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
Gapfill::opts$set(temporal_mean_est = function(x, y, x.eval, plot = FALSE){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > 4){
    x = x[nonna.idx]
    y = y[nonna.idx]
    loessfit <- loess(y~x, span = 0.1, control = loess.control(surface = "direct"))
    res = predict(loessfit, data.frame(x = x.eval))
    # res[x.eval < min(x) | x.eval > max(x)] = NA
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
})
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
n = 3
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

MSE = matrix(0, 5, 3)
Gapfill::opts$set(temporal_mean_est = Gapfill::meanCurve)
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

### MSE using "lc" method, loess
MSE = matrix(0, 5, 3)
Gapfill::opts$set(temporal_mean_est = function(x, y, x.eval, plot = FALSE){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > 4){
    x = x[nonna.idx]
    y = y[nonna.idx]
    loessfit <- loess(y~x, span = 0.1, control = loess.control(surface = "direct"))
    res = predict(loessfit, data.frame(x = x.eval))
    # res[x.eval < min(x) | x.eval > max(x)] = NA
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
})
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
### MSE using "lc" method, meanCurve
# [,1]      [,2]     [,3]
# [1,]  8675.732  7249.692 5678.767
# [2,]  9877.618  7732.764 5227.063
# [3,]  9482.812  9341.619 5073.503
# [4,] 15129.848  9185.436 5433.924
# [5,] 19338.603 11023.234 4661.543
### MSE using "lc" method, loess
# [,1]     [,2]     [,3]
# [1,] 12446.05  9834.94 6719.890
# [2,] 13709.62 11136.49 6647.363
# [3,] 17178.59 12981.21 6399.325
# [4,] 22545.50 16104.49 6472.858
# [5,] 30258.83 16777.37 6219.564

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

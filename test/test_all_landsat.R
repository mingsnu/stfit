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
## res = gapfill(year, doy, mat, 31,31,h=0, doyrange=1:90, nnr=5)
res1 = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=5)
# res = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=15, method="emp")

res = res1
## temporal trend visulization
r.list = list()
for(i in 1:365){
  r.list[[i]] = raster(matrix(res$temporal.mean[i,], 31))
}
s = stack(r.list)
levelplot(s[[seq(1, 365, 10)]], par.settings = RdBuTheme)

##########################################
#### Imputed partial missing images ######
##########################################
pdf("landsat_partial_nnr5.pdf")
for(i in 1:length(res$idx$idx.partialmissing)){
  r1 = raster(matrix(mat[res$idx$idx.partialmissing[i],], 31))
  r2 = raster(matrix(res$imputed.partial[i,], 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  # Sys.sleep(2)
}
dev.off()
pdf("landsat_partial_resid_nnr5.pdf")
for(i in 1:length(res$idx$idx.partialmissing)){
  idx = res$idx$idx.partialmissing[i]
  m = res$temporal.mean[doy[idx],]
  r1 = raster(matrix(mat[idx,] - m, 31))
  r2 = raster(matrix(res$imputed.partial[i,] - m, 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  # Sys.sleep(2)
}
dev.off()

##########################
#### Outlier images ######
##########################
for(i in 1:length(res$idx$idx.outlier)){
  r = raster(matrix(mat[res$idx$idx.outlier[i],], 31))
  print(levelplot(r, par.settings = RdBuTheme))
  Sys.sleep(2)
}

############################
##### Simulation study #####
############################
df = read_feather("../data/features_106_wide.feather")
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
set.seed(20180124)
n = 3
fidx = sample(res$idx$idx.fullyobserved, n)
pidx = sample(res$idx$idx.partialmissing, n)

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
levelplot(s, par.settings = RdBuTheme)

MSE = matrix(0, 30, 3)
for(k in 1:30){
  ## Gapfill
  res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "emp")
  ## res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k)
  
  ## imputed images
  for(i in 1:n){
    r.list[[i+2*n]] = raster(matrix(res$imputed.partial[which(res$idx$idx.partialmissing == fidx[i]),], 31))
  }
  s1 = stack(r.list)
  # levelplot(s1, par.settings = RdBuTheme, layout=c(3,3))

  ## Calculate MSE
  mat1 = values(s1)
  MSE[k,] = apply((mat1[,1:n] - mat1[,1:n+2*n])^2, 2, sum) / apply(mat1[,1:n+n], 2, function(x) sum(is.na(x)))
}
MSE
### MSE using "lc" method
# [,1]      [,2]      [,3]
# [1,]  9857.966  8919.988 12922.965
# [2,] 10898.452  9472.533 11730.300
# [3,] 10837.929 10915.312 10100.070
# [4,] 11298.320 11982.844  9127.858
# [5,] 18344.951 12745.338 17264.901
# [6,] 23521.472  8782.628 11554.994
# [7,] 17230.563 13701.762 13852.066
# [8,] 13209.412 11896.244 10808.330
# [9,] 14583.408 10844.089  9408.221
# [10,] 15805.717  8619.229  8760.538

### MSE using "emp" method
# [,1]       [,2]        [,3]
# [1,] 1059093.976 1499916.44  1932776.14
# [2,] 1512780.798 3285803.15  8329764.07
# [3,]  371528.130 3742518.94  4236515.93
# [4,]   24909.644 2464875.73 17967230.85
# [5,]   11764.987 1766265.07 11547098.86
# [6,]    9279.183  928640.34  6315799.16
# [7,]    7705.863  769264.61  1172488.15
# [8,]    7817.138  111440.42   162621.98
# [9,]    8146.553   69236.23    28290.70
# [10,]    5809.718   68367.26    66391.90
# [11,]    8131.858   24121.19    39850.45
# [12,]   11900.851   34731.78    20551.71
# [13,]   11340.432   31836.98    14483.31
# [14,]   10675.170   32904.58    33094.63
# [15,]    5115.963   54290.53    41002.96
# [16,]    6903.540   60914.18    35420.60
# [17,]    6337.571   38705.16    31425.13
# [18,]    6906.187   58873.11    26409.38
# [19,]    9101.781   29018.42    25340.30
# [20,]    8778.127   21361.35    25246.29
# [21,]    8617.795   24125.12    24233.26
# [22,]    8327.677   26035.61    22131.47
# [23,]    9047.249   20744.10    20616.19
# [24,]    8524.894   21885.15    19142.74
# [25,]    6407.207   22345.43    17318.38
# [26,]    6521.879   22580.42    15544.33
# [27,]    5177.164   22482.76    13415.80
# [28,]    5113.894   18767.85    12779.87
# [29,]    4265.739   17599.01    12448.30
# [30,]    4127.897   17715.10    12002.93
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

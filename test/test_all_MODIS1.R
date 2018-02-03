##### Test for MODIS data
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)

files = list.files("../data/MODIS2003", "*.tif$", full.names = TRUE)
## brick.list contains a list of bricks at different 
brick.list = list()
for(i in 1:length(files)){
  brick.list[[i]] = brick(files[i])
}
brick.list
## Only use 1 and 3 for test purpose
#############################################
###### Temporal difference estimation #######
#############################################
## the difference between brick one and brick three
br13 = brick.list[[1]] - brick.list[[3]]
br13 = rmOutlier(br13)
## smooth the temporal difference
br13.tm = doyMeanEst(br13)

#########################################################
###### Impute image one with three and vice versa #######
#########################################################
br1 = cover(brick.list[[1]], brick.list[[3]] + br13.tm)
br3 = cover(brick.list[[3]], brick.list[[1]] - br13.tm)

cols <- terrain.colors(30) # heat.colors(20)
rnames = gsub("MOD11A1Day2003_100x100.", "", names(br1))
levelplot(br1[[1:21]], col.regions = cols,
          main="MOD11A1Day2003_100x100",
          names.attr = rnames[1:21])
cols <- terrain.colors(30) # heat.colors(20)
rnames = gsub("MOD11A1Day2003_100x100.", "", names(br3))
levelplot(br3[[1:21]], col.regions = cols,
          main="MOD11A1Day2003_100x100",
          names.attr = rnames[1:21])


## 100 x 100 was too slow, try 30 x 30 first
nrow = ncol = 30
br0 <- crop(br1, extent(br1, 1, nrow, 1, ncol))

mat = t(values(br0))
year = rep(2003, 365)
doy = 1:365
idx = 100:250
res = gapfill1(year[idx], doy[idx], mat[idx,], nrow, ncol, h = 1, doyrange = doy[idx], nnr=5, method="lc")
res1 = gapfill1(year[idx], doy[idx], mat[idx,], nrow, ncol, h = 1, doyrange = doy[idx], nnr=5, method="emp")
## the i-th partial missing images
for(i in 1:length(res$idx$idx.partialmissing)){
  r1 = raster(matrix(mat[idx,][res$idx$idx.partialmissing[i],], nrow))
  r2 = raster(matrix(res$imputed.partial[i,], nrow))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  Sys.sleep(2)
}

for(i in 1:length(res1$idx$idx.partialmissing)){
  r1 = raster(matrix(mat[idx,][res$idx$idx.partialmissing[i],], nrow))
  r2 = raster(matrix(res$imputed.partial[i,], nrow))
  r3 = raster(matrix(mat[idx,][res1$idx$idx.partialmissing[i],], nrow))
  r4 = raster(matrix(res1$imputed.partial[i,], nrow))
  s = stack(r1, r2, r3, r4)
  print(levelplot(s, par.settings = RdBuTheme))
  Sys.sleep(2)
}

## mean plot
r.list = list()
for(i in 1:nrow(res$temporal.mean)){
  r.list[[i]] = raster(matrix(res$temporal.mean[i,], 30, byrow = TRUE))
}
s = stack(r.list)
levelplot(s[[1:30]])

## original image
r1.list = list()
for(i in 1:nrow(res$temporal.mean)){
  r1.list[[i]] = raster(matrix(mat[idx,][i,], 30, byrow = TRUE))
}
s1 = stack(r1.list)
levelplot(s1[[1:30]])

##################
###### MSE #######
##################
set.seed(20180124)
n = 3
fidx = sample(res$idx$idx.fullyobserved, n)
pidx = sample(res$idx$idx.partialmissing, n)

r.list = list()
for(i in 1:n){
  r.list[[i]] = raster(matrix(mat[idx[fidx[i]], ], nrow))
}
## apply missing patterns to fully observed images
for(i in 1:n){
  mat[idx[fidx[i]],][is.na(mat[idx[pidx[i]],])] = NA
}

for(i in 1:n){
  r.list[[i+n]] = raster(matrix(mat[idx[fidx[i]], ], nrow))
}
s = stack(r.list)
levelplot(s, par.settings = RdBuTheme)

MSE = matrix(0, 30, 3)
for(k in 1:30){
  ## Gapfill
  res = gapfill(year[idx], doy[idx], mat[idx,], nrow, ncol, h = 0, doyrange = doy[idx], nnr=k, method="emp")

  ## imputed images
  for(i in 1:n){
    r.list[[i+2*n]] = raster(matrix(res$imputed.partial[which(res$idx$idx.partialmissing == fidx[i]),], nrow))
  }
  s1 = stack(r.list)
  # levelplot(s1, par.settings = RdBuTheme, layout=c(3,3))
  
  ## Calculate MSE
  mat1 = values(s1)
  MSE[k,] = apply((mat1[,1:n] - mat1[,1:n+2*n])^2, 2, sum) / apply(mat1[,1:n+n], 2, function(x) sum(is.na(x)))
}
MSE
### "lc" method
# [,1]     [,2]     [,3]
# [1,] 49471.96 37017.93 15537.55
# [2,] 43651.48 30422.16 17162.28
# [3,] 42320.59 47624.28 19594.98
# [4,] 40163.13 37315.82 37726.06
# [5,] 33607.74 41932.51 28313.98
# [6,] 33006.29 25840.91 26837.80
# [7,] 33254.36 36096.35 25218.17
# [8,] 31045.35 86280.89 18900.57
# [9,] 24870.94 95465.85 17229.44
# [10,] 22009.78 66883.97 20777.57
### "emp" method
# [,1]      [,2]        [,3]
# [1,] 1474940.61  22893.90 16719694.60
# [2,] 1951338.69  19409.87 37555487.48
# [3,]  851314.20  41565.54  2178897.08
# [4,]  230035.33  23423.64   216848.61
# [5,]   53219.57  36575.84    71641.02
# [6,]   45441.42  16237.61    36224.45
# [7,]   44792.49  42020.50    26147.12
# [8,]   36499.00  80912.57    20578.62
# [9,]   61484.15 109107.14    21853.73
# [10,]   25826.73  52440.66    23268.74
# [11,]   18492.56  36653.20    21146.66
# [12,]   26383.64  25810.77    28730.25
# [13,]   41236.50  37586.51    24002.24
# [14,]   56004.72  23488.10    27257.53
# [15,]   58667.56  30581.73    33929.65
# [16,]   49761.69  30807.77    29358.23
# [17,]   54602.78  43147.45    29308.87
# [18,]   85860.00  53590.00    25530.34
# [19,]  347475.13  58261.02    20431.22
# [20,]  281696.67  64298.91    15833.61
# [21,]  162930.13  68310.30    16091.98
# [22,]   82517.64  76276.78    13139.48
# [23,]   50303.89  87730.67    12592.78
# [24,]   45434.25  86559.30    13313.40
# [25,]   38112.32  90910.49    12976.17
# [26,]   31947.03  82829.66    12091.38
# [27,]   24003.61  78318.98    12460.13
# [28,]   20656.90  68886.83    12004.88
# [29,]   17409.86  59897.71    10341.90
# [30,]   17409.86  59897.71    10341.90



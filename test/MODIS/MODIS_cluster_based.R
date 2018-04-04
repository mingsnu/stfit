library(raster)
library(rasterVis)
library(doParallel)
library(Matrix)
library(Gapfill)
library(foreach)

#############################
#### Clustering analysis ####
#############################
dat = readRDS("output_pixel_nnr_1/lvl3_impu.rds")
msk = getMask(dat)
## sum(msk) ##[1] 130101
dat[dat<23000] = NA
dat[dat>35000] = NA

cmat = t(dat[seq(76,290, by = 10), !msk])
dim(cmat)
apply(cmat, 2, function(x) sum(is.na(x)))
##  MYD11A1Day2010.76  MYD11A1Day2010.86  MYD11A1Day2010.96 MYD11A1Day2010.106 
##                 86                671                403                588 
## MYD11A1Day2010.116 MYD11A1Day2010.126 MYD11A1Day2010.136 MYD11A1Day2010.146 
##                688                231                258                413 
## MYD11A1Day2010.156 MYD11A1Day2010.166 MYD11A1Day2010.176 MYD11A1Day2010.186 
##                671                688                230                357 
## MYD11A1Day2010.196 MYD11A1Day2010.206 MYD11A1Day2010.216 MYD11A1Day2010.226 
##                622                 14                393                534 
## MYD11A1Day2010.236 MYD11A1Day2010.246 MYD11A1Day2010.256 MYD11A1Day2010.266 
##                163                688                 27                437 
## MYD11A1Day2010.276 MYD11A1Day2010.286 
##                688                358

## using mean to impute missing values
cmat = apply(cmat, 2, function(x) {
    x[is.na(x)] = mean(x, na.rm=TRUE)
    x
})
ncl = 100
cres = kmeans(cmat, centers = ncl, nstart = 10)
cluster = rep(0, ncol(dat))
cluster[!msk] = cres$cluster
saveRDS(cluster, "output_cluster_nnr_1/cluster.rds")

pdf("output_cluster_nnr_1/cluster_100.pdf")
print(levelplot(raster(matrix(cluster, 1200))))
dev.off()

###############################
### Level three imputation ####
###############################
dat = readRDS("output_pixel_nnr_1/lvl2_impu.rds")
dat[dat<23000] = NA
dat[dat>35000] = NA

year = rep(2010, 365)
doy = 1:365
nrow = 30; ncol=30;
k = 1
registerDoParallel(cores=14)

res3.list1 = foreach(n=1:800) %dopar% {
  ii = floor((n-1)/40) + 1
  jj = (n-1) %% 40 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  Gapfill::opts$set(temporal_mean_est = Gapfill::llreg)
  res = gapfill(year, doy, dat[,bIdx], nrow, ncol, h = 1, doyrange = 1:365, nnr=k, method="lc",
                cluster = cluster[bIdx])
  list(res$temporal.mean, res$imputed.mat)
}

saveRDS(res3.list1, "output_cluster_nnr_1/res3_list1.rds")

for(n in 1:800){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    dat[,bIdx] = res3.list1[[n]][[2]]
}

rm(res3.list1);gc()

res3.list2 = foreach(n=801:1600) %dopar% {
  ii = floor((n-1)/40) + 1
  jj = (n-1) %% 40 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
    Gapfill::opts$set(temporal_mean_est = Gapfill::llreg)
  res = gapfill(year, doy, dat[,bIdx], nrow, ncol, h = 1, doyrange = 1:365, nnr=k, method="lc",
                cluster = cluster[bIdx])
  list(res$temporal.mean, res$imputed.mat)
}

saveRDS(res3.list2, "output_cluster_nnr_1/res3_list2.rds")

for(n in 801:1600){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    dat[,bIdx] = res3.list2[[n-800]][[2]]
}

rm(res3.list2);gc()

saveRDS(dat, "output_cluster_nnr_1/lvl3_impu.rds")
## > range(dat, na.rm=TRUE)
## [1] -17746870  27830180

##############################
########## plots #############
##############################
dat0 = readRDS("../data/MYD11A1Day2010.rds")
dat[dat < min(dat0, na.rm = TRUE) - 1000] = NA
dat[dat > max(dat0, na.rm = TRUE) + 1000] = NA

pdf("output_cluster_nnr_1/lvl3_1200x_1200_partial_imputed_1-100.pdf")
for(i in 1:100){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output_cluster_nnr_1/lvl3_1200x_1200_partial_imputed_101-200.pdf")
for(i in 101:200){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output_cluster_nnr_1/lvl3_1200x_1200_partial_imputed_201-300.pdf")
for(i in 201:300){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output_cluster_nnr_1/lvl3_1200x_1200_partial_imputed_301-365.pdf")
for(i in 301:365){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()


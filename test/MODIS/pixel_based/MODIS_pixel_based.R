library(raster)
library(rasterVis)
library(doParallel)
library(Matrix)
library(Gapfill)
library(foreach)

##############################
#### Multi-lvl imputation ####
##############################
## 1. Divide 1200x1200 image into 300x300 images
## 2. Divide 300x300 image into 30x30
## Also take care of the case where there are "water body" in the image.
## if(!file.exists("../../data/MYD11A1Day2010.rds")){
## n  MODIS = stack("../../../MODIS_LST_new/2010/MYD11A1Day2010.tif")
##   dat0 = t(values(MODIS))
##   saveRDS(dat0, "../data/MYD11A1Day2010.rds")
## }
dat0 = readRDS("../data/MYD11A1Day2010.rds")

## can use getMask function to get the 'black holes'
msk0 = getMask(dat0)
pdf("output/mas0.pdf")
levelplot(raster(matrix(msk0, 1200, 1200, byrow=TRUE)))
dev.off()

year = rep(2010, 365)
doy = 1:365
idx = 1:365 ## doy index on which to to imputation
### dat is the imputation target
dat = dat0[idx,]
year = year[idx]
doy = doy[idx]

#############################
### Level one imputation ####
#############################
## systematic sampling every 50 pixels
## Using a 24 x 24 image to represent the 1200x1200 image
idx1 = c(t(outer(seq(25, 1200, by = 50), seq(25, 1200, by = 50),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
dat1 = dat[,idx1]

## visualization of the coarse lattice
rlist = list()
for(i in 1:365){
  rlist[[i]] = raster(matrix(dat1[i,], 24, 24, byrow=TRUE))
}
s = stack(rlist)

pdf("output/lvl1_24x24.pdf")
for(i in 1:12){
    print(levelplot(s[[((i-1)*30 +1):(i*30)]]))
}
levelplot(s[[361:365]])
dev.off()

## visualization of imputed results
nrow = 24; ncol=24;
res1 = gapfill(year, doy, dat1, nrow, ncol, h = 1, doyrange = doy, nnr=1, method="lc")

pdf("output/lvl1_24x_24_partial_imputed.pdf")
for(l in res1$idx$idx.partialmissing){
  r1 = raster(matrix(dat1[l,], 24))
  r2 = raster(matrix(res1$imputed.mat[l,], 24))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()

na.idx1 = is.na(dat1)
dat[,idx1][na.idx1] = res1$imputed.mat[na.idx1]
saveRDS(dat, "output/lvl1_impu.rds")

#############################
### Level two imputation ####
#############################
dat = readRDS("output/lvl1_impu.rds")
## range(dat, na.rm=TRUE)
## [1] -1638664.8   369162.9
## dat[dat<23900] = NA
## dat[dat>35000] = NA

## break image into 300 x 300 and do gapfilling
## systematic sampling every 10 pixels
idx2 = c(t(outer(seq(5, 300, by = 10), seq(5, 300, by = 10),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 300 + cidx
                 })))
nrow = 30; ncol=30;

registerDoParallel(cores=9)
res2.list = foreach(n=1:16) %dopar% {
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat2 = dat[,bIdx[idx2]]
  ## msk1 = getMask(dat2)
  ## plot(raster(matrix(msk1, 300, 300, byrow=TRUE)))
  gapfill(year, doy, dat2, nrow, ncol, h = 1, doyrange = doy, nnr=1, method="lc")
}

## impute missing values at the second level 'grid points'
for(n in 1:length(res2.list)){
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat2 = dat[,bIdx[idx2]]
  na.idx2 = is.na(dat2)
  dat[,bIdx[idx2]][na.idx2] = res2.list[[n]]$imputed.mat[na.idx2]
}
saveRDS(res2.list, "output/res2.list.rds")
saveRDS(dat, "output/lvl2_impu.rds")

###############################
### Level three imputation ####
###############################
dat = readRDS("output/lvl2_impu.rds")
## range(dat, na.rm=TRUE)
## [1] -1638664.8   369162.9
dat[dat<23900] = NA
## dat[dat>35000] = NA

nrow = 30; ncol=30;
k = 10
t1=proc.time()
registerDoParallel(cores=14)

res3.list1 = foreach(n=1:800) %dopar% {
  ii = floor((n-1)/40) + 1
  jj = (n-1) %% 40 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat3 = dat[,bIdx]
  gapfill(year, doy, dat3, nrow, ncol, h = 1, doyrange = doy, nnr=k, method="lc")$imputed.mat
}

saveRDS(res3.list1, "output/res3_list1.rds")

for(n in 1:800){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    dat[,bIdx] = res3.list1[[n]]
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
  dat3 = dat[,bIdx]
  gapfill(year, doy, dat3, nrow, ncol, h = 1, doyrange = doy, nnr=k, method="lc")$imputed.mat
}

saveRDS(res3.list2, "output/res3_list2.rds")

for(n in 801:1600){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    dat[,bIdx] = res3.list2[[n-800]]
}

rm(res3.list2);gc()

t2=proc.time()
t2-t1

sum(is.na(dat))
saveRDS(dat, "output/lvl3_impu.rds")
## > range(dat, na.rm=TRUE)
## [1] -17746870  27830180

dat[dat < min(dat0, na.rm = TRUE) - 1000] = NA
dat[dat > max(dat0, na.rm = TRUE) + 1000] = NA

pdf("output/lvl3_1200x_1200_partial_imputed_1-100.pdf")
for(i in 1:100){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_partial_imputed_101-200.pdf")
for(i in 101:200){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_partial_imputed_201-300.pdf")
for(i in 201:300){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_partial_imputed_301-365.pdf")
for(i in 301:365){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(dat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()



rm(dat)
meanmat = matrix(0, 365, 1440000)
for(n in 1:800){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    meanmat[,bIdx] = res3.list1[[n]]$temporal.mean
}
for(n in 801:1600){
    ii = floor((n-1)/40) + 1
    jj = (n-1) %% 40 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                     FUN = function(ridx, cidx){
                         (ridx-1) * 1200 + cidx
                     })))
    meanmat[,bIdx] = res3.list2[[n-800]]$temporal.mean
}

sum(is.na(meanmat))
saveRDS(meanmat, "meanest.rds")
## range(meanmat, na.rm=TRUE)
## [1] -17746870  27830180

meanmat[meanmat<23900] = NA
meanmat[meanmat>35000] = NA

pdf("output/lvl3_1200x_1200_partial_imputed_1-100.pdf")
for(i in 1:100){
    r1 = raster(matrix(dat0[i,], 1200))
    r2 = raster(matrix(meanmat[i,], 1200))
    s = stack(r1, r2)
    print(levelplot(s))
}
dev.off()



## classification
dat= readRDS("lvl3_impu.rds")
dat[dat<23900] = NA
dat[dat>35000] = NA

dim(dat)
cdat = dat[seq(100,200, by=5),]
cdat = round(cdat/1000,2)
cl = kmeans(cdat, 10, nstart=10)
pdf("plot/cluster.pdf")
rasterVis::levelplot(raster(matrix(cl$cluster, 1200, byrow = TRUE)))
dev.off()

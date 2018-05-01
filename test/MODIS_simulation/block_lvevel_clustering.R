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

dat0 = readRDS("./output/dat1.rds")

## can use getMask function to get the mask
msk0 = getMask(dat0)
pdf("output/mask0.pdf")
levelplot(raster(matrix(msk0, 1200, 1200, byrow=TRUE)))
dev.off()

dat = dat0
year = rep(2010, 365)
doy = 1:365
idx = 1:365 ## doy index on which to to imputation

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

pdf("output/lvl1/lvl1_24x24.pdf")
for(i in 1:12){
  print(levelplot(s[[((i-1)*30 +1):(i*30)]]))
}
levelplot(s[[361:365]])
dev.off()

## visualization of imputed results
nrow = 24; ncol=24;
registerDoParallel(cores = 14)
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
## res1 = gapfill_modis(doy, dat1, nrow, ncol, ncluster = 0, breaks=NULL, intermediate.dir = "./output/lvl1/")

## pdf("output/lvl1/lvl1_24x_24_partial_imputed.pdf")
## for(l in res1$idx$idx.partialmissing){
##   r1 = raster(matrix(dat1[l,], 24))
##   r2 = raster(matrix(res1$imat[l,], 24))
##   s = stack(r1, r2)
##   print(levelplot(s))
## }
## dev.off()

## na.idx1 = is.na(dat1)
## dat[,idx1][na.idx1] = res1$imat[na.idx1]
## saveRDS(dat, "output/lvl1_impu.rds")

res1 = gapfill_modis(doy, dat1, nrow, ncol, ncluster = 0, breaks=NULL, intermediate.dir = "./output/lvl1/", outlier.action = "remove")

pdf("output/lvl1/lvl1_24x_24_partial_imputed_outlier_removed.pdf")
for(l in res1$idx$idx.partialmissing){
  r1 = raster(matrix(dat1[l,], 24))
  r2 = raster(matrix(res1$imat[l,], 24))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()

na.idx1 = is.na(dat1)
dat[,idx1][na.idx1] = res1$imat[na.idx1]
saveRDS(dat, "output/lvl1_impu_outlier_removed.rds")

#############################
### Level two imputation ####
#############################
dat = readRDS("output/lvl1_impu_outlier_removed.rds")
## range(dat, na.rm=TRUE)
## [1] 23724 32939

## break image into 300 x 300 and do gapfilling
## systematic sampling every 10 pixels on each 300x300 images
idx2 = c(t(outer(seq(5, 300, by = 10), seq(5, 300, by = 10),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 300 + cidx
                 })))
nrow = 30; ncol=30;

registerDoParallel(cores=16)
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
  gapfill_modis(doy, dat2, nrow, ncol, ncluster = 0, breaks=NULL, intermediate.save = FALSE,
                intermediate.dir = "./output/lvl2/", outlier.action = "remove")
}

## visualize the level two imputation results.
for(n in 1:16){
  pdf(paste0("output/lvl2/lvl2_30x_30_partial_imputed_outlier_removed_n", n, ".pdf"))
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  for(l in 1:365){
    r1 = raster(matrix(dat[,bIdx[idx2]][l,], 30))
    r2 = raster(matrix(res2.list[[n]]$imat[l,], 30))
    s = stack(r1, r2)
    print(levelplot(s))
  }
  dev.off()    
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
  dat[,bIdx[idx2]][na.idx2] = res2.list[[n]]$imat[na.idx2]
}

saveRDS(res2.list, "output/res2.list.rds")
saveRDS(dat, "output/lvl2_impu.rds")


###############################
### Level three imputation ####
###############################
dat = readRDS("output/lvl2_impu.rds")
## range(dat, na.rm=TRUE)
## [1] -46500.62 111518.68
## [1] -1638664.8   369162.9
dat[dat<23000] = NA
dat[dat>33000] = NA

.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X %*% lmfit$coefficient)
}
Gapfill::opts$set(temporal_mean_est = customfun)

registerDoParallel(cores=10)
res3.list1 = foreach(n=1:8) %dopar% {
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat2 = dat[,bIdx]
  ## msk1 = getMask(dat2)
  ## plot(raster(matrix(msk1, 300, 300, byrow=TRUE)))
  gapfill_modis(doy, dat2, 300, 300, ncluster = 500,
                intermediate.dir = paste0("output/lvl3/block", n, "/"))[c("imat","idx")]
}
saveRDS(res3.list1, "./output/res3.list1.rds")
## res3.list1 = readRDS("./output/res3.list1.rds")

for(n in 1:8){
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat[,bIdx] = res3.list1[[n]]$imat
}

## rm(res3.list1)
gc()
res3.list2 = foreach(n=9:16) %dopar% {
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat2 = dat[,bIdx]
  ## msk1 = getMask(dat2)
  ## plot(raster(matrix(msk1, 300, 300, byrow=TRUE)))
  gapfill_modis(doy, dat2, 300, 300, ncluster = 500,
                intermediate.dir = paste0("output/lvl3/block", n, "/"))[c("imat","idx")]
}
saveRDS(res3.list2, "./output/res3.list2.rds")

for(n in 9:16){
  ii = floor((n-1)/4) + 1
  jj = (n-1) %% 4 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*300+1, ii*300), seq((jj-1)*300+1, jj*300),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 1200 + cidx
                   })))
  dat[,bIdx] = res3.list2[[n-8]]$imat
}

saveRDS(dat, "./output/dat_imputed.rds")

#### visualization
pdf("output/lvl3_1200x_1200_imputed_1-100.pdf")
for(i in 1:100){
  r1 = raster(matrix(dat0[i,], 1200))
  r2 = raster(matrix(dat[i,], 1200))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_imputed_101-200.pdf")
for(i in 101:200){
  r1 = raster(matrix(dat0[i,], 1200))
  r2 = raster(matrix(dat[i,], 1200))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_imputed_201-300.pdf")
for(i in 201:300){
  r1 = raster(matrix(dat0[i,], 1200))
  r2 = raster(matrix(dat[i,], 1200))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()

pdf("output/lvl3_1200x_1200_imputed_301-365.pdf")
for(i in 301:365){
  r1 = raster(matrix(dat0[i,], 1200))
  r2 = raster(matrix(dat[i,], 1200))
  s = stack(r1, r2)
  print(levelplot(s))
}
dev.off()


## nrow = 30; ncol=30;
## k = 10
## t1=proc.time()
## registerDoParallel(cores=16)

## res3.list1 = foreach(n=1:800) %dopar% {
##   ii = floor((n-1)/40) + 1
##   jj = (n-1) %% 40 + 1
##   ## block index
##   bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                    FUN = function(ridx, cidx){
##                      (ridx-1) * 1200 + cidx
##                    })))
##   dat3 = dat[,bIdx]
##   gapfill_modis(doy, dat3, nrow, ncol, ncluster = 0, breaks=NULL, nnr = k,
##                 intermediate.save = FALSE,
##                 intermediate.dir = "./output/lvl3/", outlier.action = "remove")
## }

## saveRDS(res3.list1, "output/res3_list1.rds")

## for(n in 1:800){
##     ii = floor((n-1)/40) + 1
##     jj = (n-1) %% 40 + 1
##     ## block index
##     bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                      FUN = function(ridx, cidx){
##                          (ridx-1) * 1200 + cidx
##                      })))
##     dat[,bIdx] = res3.list1[[n]]$imat
## }

## rm(res3.list1);gc()

## res3.list2 = foreach(n=801:1600) %dopar% {
##   ii = floor((n-1)/40) + 1
##   jj = (n-1) %% 40 + 1
##   ## block index
##   bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                    FUN = function(ridx, cidx){
##                      (ridx-1) * 1200 + cidx
##                    })))
##   dat3 = dat[,bIdx]
##   gapfill_modis(doy, dat2, nrow, ncol, ncluster = 0, breaks=NULL, nnr = k,
##                 intermediate.save = FALSE,
##                 intermediate.dir = "./output/lvl3/", outlier.action = "remove")
## }

## saveRDS(res3.list2, "output/res3_list2.rds")

## for(n in 801:1600){
##     ii = floor((n-1)/40) + 1
##     jj = (n-1) %% 40 + 1
##     ## block index
##     bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                      FUN = function(ridx, cidx){
##                          (ridx-1) * 1200 + cidx
##                      })))
##     dat[,bIdx] = res3.list2[[n-800]]$imat
## }
## rm(res3.list2);gc()

## t2=proc.time()
## t2-t1

## sum(is.na(dat))
## saveRDS(dat, "output/lvl3_impu.rds")
## ## > range(dat, na.rm=TRUE)
## ## [1] -17746870  27830180

## dat[dat < min(dat0, na.rm = TRUE) - 1000] = NA
## dat[dat > max(dat0, na.rm = TRUE) + 1000] = NA

## pdf("output/lvl3_1200x_1200_partial_imputed_1-100.pdf")
## for(i in 1:100){
##     r1 = raster(matrix(dat0[i,], 1200))
##     r2 = raster(matrix(dat[i,], 1200))
##     s = stack(r1, r2)
##     print(levelplot(s))
## }
## dev.off()

## pdf("output/lvl3_1200x_1200_partial_imputed_101-200.pdf")
## for(i in 101:200){
##     r1 = raster(matrix(dat0[i,], 1200))
##     r2 = raster(matrix(dat[i,], 1200))
##     s = stack(r1, r2)
##     print(levelplot(s))
## }
## dev.off()

## pdf("output/lvl3_1200x_1200_partial_imputed_201-300.pdf")
## for(i in 201:300){
##     r1 = raster(matrix(dat0[i,], 1200))
##     r2 = raster(matrix(dat[i,], 1200))
##     s = stack(r1, r2)
##     print(levelplot(s))
## }
## dev.off()

## pdf("output/lvl3_1200x_1200_partial_imputed_301-365.pdf")
## for(i in 301:365){
##     r1 = raster(matrix(dat0[i,], 1200))
##     r2 = raster(matrix(dat[i,], 1200))
##     s = stack(r1, r2)
##     print(levelplot(s))
## }
## dev.off()



## rm(dat)
## meanmat = matrix(0, 365, 1440000)
## for(n in 1:800){
##     ii = floor((n-1)/40) + 1
##     jj = (n-1) %% 40 + 1
##     ## block index
##     bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                      FUN = function(ridx, cidx){
##                          (ridx-1) * 1200 + cidx
##                      })))
##     meanmat[,bIdx] = res3.list1[[n]]$temporal.mean
## }
## for(n in 801:1600){
##     ii = floor((n-1)/40) + 1
##     jj = (n-1) %% 40 + 1
##     ## block index
##     bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
##                      FUN = function(ridx, cidx){
##                          (ridx-1) * 1200 + cidx
##                      })))
##     meanmat[,bIdx] = res3.list2[[n-800]]$temporal.mean
## }

## sum(is.na(meanmat))
## saveRDS(meanmat, "meanest.rds")
## ## range(meanmat, na.rm=TRUE)
## ## [1] -17746870  27830180

## meanmat[meanmat<23900] = NA
## meanmat[meanmat>35000] = NA

## pdf("output/lvl3_1200x_1200_partial_imputed_1-100.pdf")
## for(i in 1:100){
##     r1 = raster(matrix(dat0[i,], 1200))
##     r2 = raster(matrix(meanmat[i,], 1200))
##     s = stack(r1, r2)
##     print(levelplot(s))
## }
## dev.off()


library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
##### MYD Day plot
dat0 = readRDS("./data/MYD11A1Day2010.rds")
# colthm = RdBuTheme()
# colthm$regions$col = rev(colthm$regions$col)
msk = readRDS("./data/msk.rds")
dat0[,msk] = NA
dat0 = dat0/100 ## the unit after transform is 'K'
## stackRaster
r.list = list()
dd = seq(1,365, 30)
for(i in 1:length(dd)){
  r.list[[i]] = raster(matrix(dat0[dd[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
png("./output/MYD11A1Day2010.png", width=700, height=550)
print(levelplot(s[[1:12]], names.attr=as.character(dd[1:12]), main = "",
                layout=c(4,3)), at = seq(240, 325, length.out = 20))
dev.off()

##### MOD Day plot
dat0 = readRDS("./data/MOD11A1Day2010.rds")
# colthm = RdBuTheme()
# colthm$regions$col = rev(colthm$regions$col)
msk = readRDS("./data/msk.rds")
dat0[,msk] = NA
dat0 = dat0/100 ## the unit after transform is 'K'
## stackRaster
r.list = list()
dd = seq(1,365, 30)
for(i in 1:length(dd)){
  r.list[[i]] = raster(matrix(dat0[dd[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
png("./output/MOD11A1Day2010.png", width=700, height=550)
print(levelplot(s[[1:12]], names.attr=as.character(dd[1:12]), main = "",
                layout=c(4,3)), at = seq(240, 325, length.out = 20))
dev.off()

## percentage of missing values
N = sum(!msk)
pctmissing = apply(dat0[,!msk], 1, function(x) sum(is.na(x))/N)
summary(pctmissing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2035  0.6150  0.7471  0.7247  0.8734  1.0000
## the overall percentage of missing
sum(is.na(dat0[,!msk]))/(N*nrow(dat0))
# [1] 0.7246504
saveRDS(pctmissing, "./output/pctmissing.rds")
pdf("./output/pctmissing.pdf", width=6.9, height=5.5)
par(mar=c(4.2,4.2,0.2,0.2))
hist(pctmissing, breaks = 20, xlab = "Percentage of missing data over land areas",  main="",
     xlim=c(0.2,1))
dev.off()

## percentage of missing values after daily merging
dat1 = readRDS("./data/MYD11A1Day2010_daily_imputed_lm.rds")
pctmissing = apply(dat1[,!msk], 1, function(x) sum(is.na(x))/N)
summary(pctmissing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1441  0.4456  0.6087  0.6120  0.7946  0.9960 
sum(is.na(dat1[,!msk]))/(N*nrow(dat1))
# 0.6120424

## mask
png("./output/msk.png", width=620, height=550)
myTheme <- rasterTheme(region=c('white', gray.colors(1))) 
r = raster(matrix(msk, 1200, byrow = TRUE))
r = ratify(r)
rat <- levels(r)[[1]]
rat$landcover <- c('Land', 'Water')
rat$class <- c('a', 'b')
levels(r) <- rat
levelplot(r, par.settings=myTheme, margin = FALSE)
dev.off()


######## day imputation plots
### level 3 plot
dat = readRDS("./data/MYD11A1Day2010_daily_imputed_lm.rds")
msk = readRDS("./data/msk.rds")
dat[, msk] = NA
idx1 = c(t(outer(seq(25, 1200, by = 50), seq(25, 1200, by = 50),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
res1 = readRDS("./output_day/res1.rds")
l =  res1$idx$idx.partialmissing[170]
r1 = raster(matrix(dat[,idx1][l,], 24, byrow = TRUE))
r2 = raster(matrix(res1$imat[l,], 24, byrow = TRUE))
s1 = stack(r1, r2)
print(levelplot(s1, names.attr=c("Level 3 image", "Imputed level 3 image")))

### level 2 plot
dat = readRDS("./output_day/lvl1_impu_outlier_removed.rds")
## visualize the level two imputation results.
bIdx = c(t(outer(seq(5, 1200, by = 10), seq(5, 1200, by = 10),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
r3 = raster(matrix(dat[,bIdx][l,], 120, byrow = TRUE))
dat1 = readRDS("output_day/lvl2_impu.rds")
r4 = raster(matrix(dat1[,bIdx][l,], 120, byrow = TRUE))
s2 = stack(r3, r4)
print(levelplot(s2, names.attr=c("Level 2 image", "Imputed level 2 image")))

### level 1 plot
dat = readRDS("./output_day/lvl2_impu.rds")
r5 = raster(matrix(dat[l,], 1200, byrow = TRUE))
dat1 = readRDS("output_day/dat_imputed.rds")
r6 = raster(matrix(dat1[l,], 1200, byrow = TRUE))
s3 = stack(r5, r6)
print(levelplot(s3, names.attr=c("Original image", "Imputed image")))

### individual plots
r.list = list(r1, r2, r3, r4, r5, r6)
ll = seq(279,317, length.out = 16)
for(i in 1:6){
  pdf(paste0("./output_day/lst_day", l, "_", i, ".pdf"), width=5.7, height=5.2)
  print(levelplot(r.list[[i]]/100, margin = FALSE, at = ll))
  dev.off()
}

### plot T2 vs lm imputed T2
dat0 = readRDS("./data/MYD11A1Day2010.rds")
msk = readRDS("./data/msk.rds")
dat0[,msk] = NA
dat0 = dat0/100 ## the unit after transform is 'K'
dat1 = readRDS("./data/MYD11A1Day2010_daily_imputed_lm.rds")
dat1 = dat1/100
r.list = list()
idx = c(1, 201, 280)
for(i in 1:length(idx)){
  r.list[[i]] = raster(matrix(dat0[idx[i],], 1200, byrow = TRUE))
  r.list[[3 + i]] = raster(matrix(dat1[idx[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
levelplot(s, names.attr=c("DOY 1", "DOY 201", "DOY 280", "Imputed DOY 1", "Imputed DOY 201", "Imputed DOY 280"))

## lm imiputed T2 MYD plot
dat1 = readRDS("./data/MYD11A1Day2010_daily_imputed_lm.rds")
dat1 = dat1/100
dd = seq(1,365, 30)
r.list=list()
for(i in 1:length(dd)){
  r.list[[i]] = raster(matrix(dat1[dd[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
png("./output/MYD11A1Day2010_lm_imputed.png", width=700, height=550)
print(levelplot(s[[1:12]], names.attr=as.character(dd[1:12]), main = "",
                layout=c(4,3)), at = seq(240, 325, length.out = 20))
dev.off()




dat1 = readRDS("./data/MYD11A1Day2010_daily_imputed_lm.rds")
dat1 = dat1/100
## cluster based mean estimation
bIdx = c(t(outer(1:300, 1:300, 
                 FUN = function(ridx, cidx) {
                   (ridx - 1) * 1200 + cidx
                 })))
dat_b1 = dat1[,bIdx]
levelplot(mat2stack(dat_b1, 300, byrow=TRUE)[[1:3]])

registerDoParallel(cores = 8)
## stfit::opts$set(temporal_mean_est = stfit::spreg)
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)
msk1 = getMask(dat_b1)
meanest1 = meanEst(1:365, dat_b1, doyeval = 1:365, redo = FALSE, msk = msk1,
                   clipRange = c(230, 350), clipMethod = "truncate", img.nrow = 300, 
                   img.ncol = 300)
saveRDS(meanest1, "./output/meanest1.rds")
##levelplot(mat2stack(meanest1$meanmat, 300, byrow=TRUE)[[1:3]], at = seq(253, 273, 1))

cmat = t(meanest1$meanmat[seq(10, 350, by = 10), !msk1])
cres = kmeans(cmat, centers = 500, nstart = 10, iter.max = 50)
cluster = rep(0, ncol(dat_b1))
cluster[!msk1] = cres$cluster
saveRDS(cluster, "./output/dat_b1_cluster.rds")
meanest_cl = meanEst(1:365, dat_b1, doyeval = 1:365, redo = FALSE,
                     cluster = cluster, msk = msk1, clipRange = c(230, 350), 
                     clipMethod = "truncate", img.nrow = 300, 
                     img.ncol = 300)
## levelplot(mat2stack(meanest_cl$meanmat, 300, byrow=TRUE)[[1:3]], at = seq(253, 273, 1))
levelplot(mat2stack(rbind(meanest1$meanmat[1,], meanest_cl$meanmat[1,]), 300, byrow=TRUE), 
          at = seq(253, 273, 1), names.attr=c("Mean(DOY1/block1)", "Cluster-based mean(DOY1/block1)"))
levelplot(raster(matrix(cluster, 300, byrow=TRUE)), margin=FALSE)

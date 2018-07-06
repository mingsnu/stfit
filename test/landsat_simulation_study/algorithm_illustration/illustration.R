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

df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

registerDoParallel(6)
stfit::opts$set(temporal_mean_est = spreg)
mat = mat0

pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)
pmat = readRDS("../missing_pattern/output/missing_pattern.rds")
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(132, 261, 265, 615, 657) 
fidx = 198
year[fidx] ## 2004
doy[fidx] ## 228
fmat = mat0[fidx, ]
missing.idx = is.na(pmat[8,])
mat[fidx, missing.idx] = NA

r.list=list()
#####
## F6; 
r.list[[1]] = raster(matrix(fmat, 31)) ## original 
r.list[[2]] = raster(matrix(mat[fidx,], 31)) ## artifical missing value

levelplot(raster(matrix(fmat, 31)), par.settings = colthm, margin=FALSE)
levelplot(raster(matrix(mat[fidx,], 31)), par.settings = colthm, margin=FALSE)

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
## saveRDS(meanest, "output/meanest.rds")
mat_mean_imp = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

r.list[[3]] = raster(matrix(mat_mean_imp[fidx,], 31))
levelplot(raster(matrix(mat_mean_imp[fidx,], 31)), par.settings = colthm, margin=FALSE)

############################
#### teffect estimation ####
############################
## calculate the residuals
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

r.list[[4]] = raster(matrix(rmat[fidx,], 31))
levelplot(raster(matrix(rmat[fidx,], 31)), par.settings = colthm, margin=FALSE)

teffarray = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = 100, h.sigma2 = 300)
## saveRDS(teffarray, "output/teffarray.rds")

#### RMSE using mean+teffect
mat_teff_imp = mat_mean_imp
yearidx = unlist(lapply(year, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
doyidx = unlist(lapply(doy, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
for(i in 1:nrow(mat_teff_imp)){
  mat_teff_imp[i,] = mat_teff_imp[i,] + teffarray[yearidx[i], doyidx[i],]
}


### visualization 1
pdf("output/teff_pixel157.pdf", width =6, height=5.5)
par(mar=c(4,4,0.1,0.1))
yeareval = as.numeric(dimnames(teffarray)[[1]])
doyeval = as.numeric(dimnames(teffarray)[[2]])
## year 2004
plot(doy, rmat[,157], col = gray(0.8), pch=19, ylab = "Residuals", xlab= "DOY")
abline(h=0, lty=2, col=gray(0.5))
ind = which(year==2004)
points(doy[ind], rmat[ind,157], col = 2, pch=19)
lines(doyeval, teffarray[5,,157], col = 2, lwd = 2)
dev.off()

r.list[[5]] = raster(matrix(teffarray[yearidx[fidx], doyidx[fidx],], 31))
levelplot(raster(matrix(teffarray[yearidx[fidx], doyidx[fidx],], 31)), par.settings = colthm, margin=FALSE)
levelplot(raster(matrix(mat_teff_imp[fidx,], 31)), par.settings = colthm, margin=FALSE)


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
#saveRDS(seffmat, "output/seffmat30.rds")
mat_seff_imp[fidx, !missing.idx] = mat0[fidx, !missing.idx]

r.list[[6]] = raster(matrix(rmat[fidx,], 31))
r.list[[7]] = raster(matrix(seffmat[fidx,], 31))
r.list[[8]] = raster(matrix(mat_seff_imp[fidx,], 31))

levelplot(raster(matrix(rmat[fidx,], 31)), par.settings = colthm, margin=FALSE)
levelplot(raster(matrix(seffmat[fidx,], 31)), par.settings = colthm, margin=FALSE)
levelplot(raster(matrix(mat_seff_imp[fidx,], 31)), par.settings = colthm, margin=FALSE)

pdf("output/fullobs.pdf", width =6, height=5.5)
levelplot(r.list[[1]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
pdf("output/partialmissing.pdf", width =6, height=5.5)
levelplot(r.list[[2]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
pdf("output/mean.pdf", width =6, height=5.5)
levelplot(r.list[[3]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
pdf("output/mean_resid.pdf", width =6, height=5.5)
levelplot(r.list[[4]], par.settings = colthm, margin = FALSE)
dev.off()
pdf("output/teff.pdf", width =6, height=5.5)
levelplot(r.list[[5]], par.settings = colthm, margin = FALSE)
dev.off()
pdf("output/teff_resid.pdf", width =6, height=5.5)
levelplot(r.list[[6]], par.settings = colthm, margin = FALSE)
dev.off()
pdf("output/seff.pdf", width =6, height=5.5)
levelplot(r.list[[7]], par.settings = colthm, margin = FALSE)
dev.off()
pdf("output/imputed_theoretical.pdf", width =6, height=5.5)
levelplot(r.list[[8]], par.settings = colthm, margin = FALSE)
dev.off()



          
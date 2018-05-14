library(raster)
library(rasterVis)
library(stfit)
library(doParallel)
library(Matrix)
library(foreach)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)
smapdat = readRDS("../data/smapdat.rds")
doy = 1:365
idx = 1:365 ## doy index on which to to imputation
## dimention 406 964
#############################
### Level one imputation ####
#############################
## systematic sampling every 20 pixels
## Using a 21 x 49 image to represent the 406x964 image
idx1 = c(t(outer(seq(3, 406, by = 20), seq(2, 964, by = 20),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 964 + cidx
                 })))
dat1 = smapdat[,idx1]

## visualization of the coarse lattice ============================
rlist = list()
for(i in 1:365){
  rlist[[i]] = raster(matrix(dat1[i,], 21, 49))
}
s = stack(rlist)

pdf("../output/lvl1/lvl1_21x49.pdf")
for(i in 1:12){
  print(levelplot(s[[((i-1)*30 +1):(i*30)]], par.settings=colthm))
}
levelplot(s[[361:365]], par.settings=colthm, layout=c(6,5))
dev.off()

## visualization of imputed results ====================
nrow = 21; ncol=49;
registerDoParallel(cores = 6)
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
res1 = gapfill_modis(doy, dat1, nrow, ncol, ncluster = 0, breaks=NULL, intermediate.dir = "./output/lvl1/")

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

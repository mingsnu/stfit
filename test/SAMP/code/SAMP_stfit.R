library(raster)
library(rasterVis)
library(stfit)
library(doParallel)
library(Matrix)
library(foreach)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)
smapdat = readRDS("../data/smapdat_rowstack.rds")
doy = 1:365
doyeval = 1:365
# msk0 = getMask(smapdat)
# levelplot(raster(matrix(msk0, 406, 964, byrow=TRUE), xmn=0, xmx=49, ymn=0, ymx=21))

# msk1 = msk0
# msk1[0:405*964+400] = FALSE
# msk1[0:405*964+401] = FALSE
# msk1[0:405*964+650] = FALSE
# msk1[0:405*964+651] = FALSE
# msk1[80*964 + 1:964] = FALSE
# msk1[81*964 + 1:964] = FALSE
# msk1[161*964 + 1:964] = FALSE
# msk1[162*964 + 1:964] = FALSE
# msk1[242*964 + 1:964] = FALSE
# msk1[243*964 + 1:964] = FALSE
# msk1[323*964 + 1:964] = FALSE
# msk1[324*964 + 1:964] = FALSE
# levelplot(raster(matrix(msk1, 406, 964, byrow=TRUE), xmn=0, xmx=49, ymn=0, ymx=21))
# 

## dimention 406 964
#############################
### Level one imputation ####
#############################
## systematic sampling every 20 pixels
## Using a 21 x 49 image to represent the 406x964 image
grid1.row = seq(3, 406, by = 10)
grid1.col = seq(2, 964, by = 10)
idx1 = c(t(outer(grid1.row, grid1.col,
                 FUN = function(ridx, cidx){
                   (ridx-1) * 964 + cidx
                 })))
dat1 = smapdat[,idx1]
nrow1 = length(grid1.row)
ncol1 = length(grid1.col)
## visualization of the coarse lattice ============================
rlist = list()
for(i in 1:365){
  rlist[[i]] = raster(matrix(dat1[i,], nrow1, ncol1, byrow = TRUE), xmn=0, xmx=49, ymn=0, ymx=21)
}
s = stack(rlist)

pdf("./output/lvl1/lvl1_", nrow1, "x", ncol1, ".pdf")
for(i in 1:12){
  print(levelplot(s[[((i-1)*30 +1):(i*30)]], par.settings=colthm))
}
levelplot(s[[361:365]], par.settings=colthm, layout=c(5,6))
dev.off()
pdf("./output/lvl1/lvl1_original.pdf")
levelplot(s[[seq(1,365,by=30)]], par.settings=colthm)
dev.off()

## visualization of imputed results ====================
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
mat = dat1
msk = getMask(mat)
pdf("./output/lvl1/lvl1_msk.pdf")
levelplot(raster(matrix(msk, nrow1, byrow = TRUE), xmn=0, xmx=49, ymn=0, ymx=21))
dev.off()


clipRange = c(0, 0.8)
clipMethod = "nnr"
meanest = meanEst(doy, mat, doyeval = doyeval, msk = msk, clipRange=clipRange, clipMethod = clipMethod,
                  img.nrow = nrow1, img.ncol = ncol1)
ms = mat2stack(meanest$meanmat, nrow1, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
pdf("./output/lvl1/lvl1_mean.pdf")
levelplot(ms[[seq(1,365,by=30)]], par.settings=colthm)
dev.off()
#### if using cluster
# ncluster = 50
# cmat = t(meanest$meanmat[seq(10, 350, by = 10), !msk])
# cres = kmeans(cmat, centers = ncluster, nstart = 10, 
#               iter.max = 50)
# cluster = rep(0, ncol(mat))
# cluster[!msk] = cres$cluster
# meanest_cl = meanEst(doy, mat, doyeval = doyeval, 
#                      cluster = cluster, msk = msk, clipRange = clipRange, 
#                      clipMethod = clipMethod, img.nrow = nrow1, 
#                      img.ncol = ncol1)

meanest_cl = meanest
outlier.img.idx = meanest_cl$idx$idx.outlier
for(i in outlier.img.idx){
  mat[outlier.img.idx,] = NA
}
## remove outlier pixels
for(i in 1:length(meanest_cl$outlier$outidx)){
  mat[meanest_cl$outlier$outidx[i], meanest_cl$outlier$outlst[[i]]] = NA
}
## calculate the residuals
rmat = mat - meanest_cl$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest_cl$doyeval)),]
rs = mat2stack(rmat, nrow1, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
pdf("./output/lvl1/lvl1_resid.pdf")
levelplot(rs[[seq(1,365,by=30)]], par.settings=colthm, 
          at=seq(-0.1, 0.15, length.out = 20))
dev.off()

## there are some pixels never have overlapping observations!!
## rmat[,c(583, 881)]
seffest = seffEst(rmat, nrow1, ncol1, nnr=30, h.cov = 2, h.sigma2 = 2, msk = msk)
rmat.imp = seffest$seffmat
rmat.imp[rmat.imp==0] = NA
res = mat2stack(rmat.imp, nrow1, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
pdf("./output/lvl1/lvl1_residest.pdf")
levelplot(res[[seq(1,365,by=30)]], par.settings=colthm, 
          at=seq(-0.1, 0.15, length.out = 20))
dev.off()
# cor(c(rmat), c(rmat.imp), use= "complete.obs")
# 0.833
# plot(c(rmat), c(rmat.imp)); abline(a=0,b=1)
# rdiff = rmat - rmat.imp
# sum(rdiff^2, na.rm=TRUE) / sum(!is.na(rdiff))
# ## 0.0005179786
dat1.imp = meanest$meanmat + seffest$seffmat
dat1[is.na(dat1)] = dat1.imp[is.na(dat1)]
dat1s = mat2stack(dat1, nrow1, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
levelplot(dat1s[[seq(1,365,by=30)]], par.settings=colthm)

smapdat[,idx1][is.na(dat1)] = dat1.imp[is.na(dat1)]
# tmpr = raster(matrix(smapdat[1,], 406, 964, byrow = TRUE), 
#               xmn=0, xmx=964, ymn=0, ymx=406)
# levelplot(tmpr)
saveRDS(smapdat, "./output/smapdat_lvl1_imputed.rds")


smapdat = readRDS("./output/smapdat_lvl1_imputed.rds")
## visualization of imputed results ====================
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
msk0 = getMask(smapdat)
clipRange = c(0, 0.8)
clipMethod = "nnr"
meanest2 = meanEst(doy, smapdat, doyeval = doyeval, msk = msk0, clipRange=clipRange,
                   clipMethod = clipMethod, img.nrow = 406, img.ncol = 964)
saveRDS(meanest2, "./output/smapdat_meanest.rds")
ms = mat2stack(meanest2$meanmat, 406, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
pdf("./output/smapdat_meanest.pdf")
levelplot(ms[[seq(1,365,by=30)]], par.settings=colthm)
dev.off()
# levelplot(raster(matrix(smapdat[1,], 406, byrow=TRUE), xmn=0, xmx=49, ymn=0, ymx=21),
#           par.settings=colthm, at = seq(0,0.8, length.out = 20))
# levelplot(ms[[1]], par.settings=colthm, at = seq(0,0.8, length.out = 20))


bIdx = c(outer(81:162*964, 400:650, "+"))
bIdx = c(outer(0:80*964, 1:400, "+"))
brmat = msk0[bIdx]
sum(!brmat)

seffest2 = seffEst(rmat[, bIdx], breaks$img.nrow, breaks$img.ncol, nnr = nnr, 
        h.cov = h.scov, h.sigma2 = h.ssigma2, msk = msk[bIdx])




breaks = list(block.nrow = 10, block.ncol = 10, img.nrow = 30, img.ncol = 30)
nblocks = breaks$block.nrow*breaks$block.ncol
res3.list = foreach(n=1:nblocks) %dopar% {
  ii = floor((n-1)/breaks$block.ncol) + 1
  jj = (n-1) %% breaks$block.ncol + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*breaks$img.nrow+1, ii*breaks$img.nrow), 
                   seq((jj-1)*breaks$img.ncol+1, jj*breaks$img.ncol),
                   FUN = function(ridx, cidx){
                     (ridx-1) * img.ncol + cidx
                   })))
  seffEst(rmat[, bIdx], breaks$img.nrow, breaks$img.ncol, nnr = nnr, 
          h.cov = h.scov, h.sigma2 = h.ssigma2, msk = msk[bIdx])
}

# pdf("output/lvl1/lvl1_21x_49_partial_imputed.pdf")
# for(l in res1$idx$idx.partialmissing){
#   r1 = raster(matrix(dat1[l,], 24))
#   r2 = raster(matrix(res1$imat[l,], 24))
#   s = stack(r1, r2)
#   print(levelplot(s))
# }
# dev.off()

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

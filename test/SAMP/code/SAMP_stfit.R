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
  rlist[[i]] = raster(matrix(dat1[i,], 21, 49, byrow = TRUE), xmn=0, xmx=49, ymn=0, ymx=21)
}
s = stack(rlist)

pdf("./output/lvl1/lvl1_21x49.pdf")
for(i in 1:12){
  print(levelplot(s[[((i-1)*30 +1):(i*30)]], par.settings=colthm))
}
levelplot(s[[361:365]], par.settings=colthm, layout=c(5,6))
dev.off()
levelplot(s[[seq(1,365,by=30)]], par.settings=colthm)

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
#res1 = gapfill_modis(doy, dat1, nrow, ncol, ncluster = 0, breaks=NULL, intermediate.dir = "../output/lvl1/")
mat = dat1
msk = getMask(mat)
levelplot(raster(matrix(msk,21, byrow = TRUE), xmn=0, xmx=49, ymn=0, ymx=21))
doyeval = 1:365
clipRange = c(0, 0.8)
clipMethod = "nnr"
img.nrow = 21
img.ncol = 49
meanest = meanEst(doy, mat, doyeval = doyeval, msk = msk, clipRange=clipRange, clipMethod = clipMethod,
                  img.nrow = img.nrow, img.ncol = img.ncol)
ms = mat2stack(meanest$meanmat, 21, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
levelplot(ms[[seq(1,365,by=30)]], par.settings=colthm)
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

rs = mat2stack(rmat, 21, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
levelplot(rs[[seq(1,365,by=30)]], par.settings=colthm, at=seq(-0.1, 0.15, length.out = 20))

## there are some pixels never have overlapping observations!!
## rmat[,c(583, 881)]

seffest = seffEst(rmat, 21, 49, nnr=30, h.cov = 2, h.sigma2 = 2, msk = msk)
res = mat2stack(seffest$seffmat, 21, byrow=TRUE, xmn=0, xmx=49, ymn=0, ymx=21)
levelplot(res[[seq(1,365,by=30)]], par.settings=colthm, at=seq(-0.1, 0.15, length.out = 20))


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

dat = readRDS("./data/MOD11A1Day2010.rds")
msk = readRDS("./data/msk.rds")
dat[, msk] = NA
year = rep(2010, 365)
doy = 1:365
idx = 1:365 ## doy index on which to to imputation
registerDoParallel(cores = 8)
## stfit::opts$set(temporal_mean_est = stfit::spreg)
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=27))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)
## Using a 24 x 24 image to represent the 1200x1200 image
idx1 = c(t(outer(1:30, 1:30,
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
dat1 = dat[,idx1]
nrow = 30; ncol=30;


res1 = gapfill_modis(doy, dat1, nrow, ncol, nnr = 30,
                     ncluster = 0, breaks=NULL,
                     outlier.action = "remove", intermediate.save = FALSE)
rlist = list()
for(i in 1:10){
  rlist[[i]] = raster(matrix(res1$imat[i,], nrow, ncol, byrow=TRUE))
}
simpu= stack(rlist)
levelplot(simpu)


meanest1 = meanEst(doy, dat1, doyeval = 1:365, 
                  clipRange = c(23000, 35000), clipMethod = "nnr", img.nrow = 30, 
                  img.ncol = 30)
rlist = list()
for(i in 1:10){
  rlist[[i]] = raster(matrix(meanest1$meanmat[i,], nrow, ncol, byrow=TRUE))
}
smean = stack(rlist)
levelplot(smean)

par(mar=c(4.3,4.3,0.3,0.3))
pixidx = which.max(apply(dat1, 2, function(x) min(which(!is.na(x)))))

plot(1:365, dat1[,pixidx]/100, ylim = c(24500,32000)/100, pch = 19, col=gray(0.5, alpha=0.8), 
     xlab = "Day of the year", ylab = "LST")
lines(1:365, customfun(1:365, dat1[,pixidx])/100, lwd = 2)
#lines(1:365, stfit::smooth_spline(1:365, dat1[,337], 1:365))

cmat = t(meanest1$meanmat[seq(10, 350, by = 10),])
cres = kmeans(cmat, centers = 50, nstart = 10, 
              iter.max = 50)
cluster = cres$cluster
cidx = which(cluster == cluster[pixidx])
plot(rep(1:365, length(cidx)), c(dat1[,cidx])/100, ylim = c(24500,32000)/100, pch = 19, col=gray(0.5, alpha=0.2), 
     xlab = "Day of the year", ylab = "LST")
lines(1:365, customfun(rep(1:365, length(cidx)), c(dat1[,cidx]))/100, lwd = 2)

      
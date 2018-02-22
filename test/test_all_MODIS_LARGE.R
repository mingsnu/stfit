library(raster)
library(rasterVis)
## Divide 300x300 image into 30x30 images and combine together after imputation seperately
## Also take care of the case where there are "black holes" in the image.
dat = readRDS("../data/dat_300x300.rds")
## can use getMask function to get the 'black holes'
msk0 = getMask(dat)
plot(raster(matrix(msk0, 300, 300, byrow=TRUE)))

################################################
### Break down into 30x30 and do gapfilling ####
################################################
year = rep(2003, 365)
doy = 1:365
idx = 100:250

idx.mat = matrix(1:90000, 300, byrow = TRUE)
nrow = 30; ncol=30;
# nrow = 50; ncol=50;

res.list = list()
k=0
for(i in 1:6){
  cat("i = ", i, "\n")
  for(j in 1:6){
    cat("j = ", j, "\n")
    k = k+1
    if(i + j ==2) next
    mat = dat[idx, c(t(idx.mat[seq(1+(i-1)*nrow, i*nrow), seq(1+(j-1)*ncol, j*ncol)]))]
    res.list[[k]] = gapfill(year[idx], doy[idx], mat, nrow, ncol, h = 1, doyrange = doy[idx], nnr=1, method="lc")
  }
}
# saveRDS(res.list, file = "../data/dat_300x300_30x30block_result.rds")


#######################################
### Imputed images: block by block ####
#######################################
res.list = readRDS("../data/dat_300x300_30x30block_result.rds")
## imputation results, block by block
k=0
for(i in 1:5){
  for(j in 1:5){
    k = k + 1
    for(l in res.list[[k]]$idx$idx.partialmissing){
      r1 = raster(matrix(dat[idx, c(t(idx.mat[seq(1+(i-1)*nrow, i*nrow), seq(1+(j-1)*ncol, j*ncol)]))][l,], 30))
      r2 = raster(matrix(res.list[[k]]$imputed.mat[l,], 30))
      s = stack(r1, r2)
      print(levelplot(s))
      Sys.sleep(2)
    }
  }
}

#################################
### Imputed images: combined ####
#################################
dat.imputed = dat
k=0
for(i in 1:10){
  for(j in 1:10){
    k = k+1
    dat.imputed[idx, c(t(idx.mat[seq(1+(i-1)*nrow, i*nrow), seq(1+(j-1)*ncol, j*ncol)]))] =
      res.list[[k]]$imputed.mat
    }
}
dat.imputed[dat.imputed <0 | dat.imputed>40000] = NA

pdf("MODIS_30x30_plot.pdf")
for(t in 1:length(idx)){
  if(all(is.na(dat[idx[t],]))) next
  r1 = raster(matrix(dat[idx[t],], 300, byrow = TRUE))
  r2 = raster(matrix(dat.imputed[idx[t], ], 300, byrow = TRUE))
  plot(levelplot(stack(r1, r2)))
  print(t)
}
dev.off()



## the i-th outlier images
for(i in 1:2){
  r1 = raster(matrix(dat[idx,][i,], 300))
  r2 = raster(matrix(imp.dat[idx,][i,], 300))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  Sys.sleep(2)
}
plot(raster(matrix(imp.dat[i,], 300, 300, byrow=TRUE)))



i = j = 1
mat = dat[, c(t(idx.mat[seq(1+(i-1)*nrow, i*nrow), seq(1+(j-1)*ncol, j*ncol)]))]
levelplot(raster(matrix(dat1[2,], nrow)))
msk = getMask(mat)
plot(raster(matrix(msk, nrow, ncol, byrow=TRUE)))

# r = raster(matrix(dat[2,], 300))
# levelplot(r)
# levelplot(crop(r, extent(r, 1+(i-1)*nrow, i*nrow, 1+(j-1)*ncol, j*ncol)))

res = gapfill1(year[idx], doy[idx], mat[idx,], nrow, ncol, h = 1, doyrange = doy[idx], nnr=2, method="lc")


## the i-th partial missing images
for(i in 1:length(res$idx$idx.partialmissing)){
  r1 = raster(matrix(mat[idx,][res$idx$idx.partialmissing[i],], nrow))
  r2 = raster(matrix(res$imputed.partial[i,], nrow))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  Sys.sleep(2)
}



################################
######## Merge together ########
################################
MODIS_merge = mosaic(MODIS1, MODIS2, fun = mean)
plot(MODIS_merge[[1:40]])





## mask
msk = raster("../data/bf.h11v04.tif")
plot(msk)
msk1 = crop(msk, extent(msk, 601, 900, 601, 900))
plot(msk1)



##
aa = cov(matrix(rnorm(100), 25))
aa=cbind(aa,0)
aa = rbind(aa, 0)
aa[,5] = aa[,4]
aa[5,] = aa[4,]
aa[4,] = 0
aa[,4] = 0
ev = eigen(aa)
lambda = ev$values[1:4]
phi = ev$vectors[,1:4]
phi %*% diag(lambda) %*% t(phi)
aa

tmp = t(phi) * lambda
phi %*% tmp


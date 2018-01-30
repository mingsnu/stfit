## Divide 300x300 image into 30x30 images and combine together after imputation seperately
## Also take care of the case where there are "black holes" in the image.
dat = readRDS("../data/dat_300x300.rds")
## can use getMask function to get the 'black holes'
msk = getMask(dat)
plot(raster(matrix(msk, 300, 300, byrow=TRUE)))



## Test mosaic function
nrow = 25; ncol=25;
i = j = 1
MODIS1 = crop(MODIS, extent(MODIS, 1+(i-1)*nrow, i*nrow, 1+(j-1)*ncol, j*ncol))
plot(MODIS1[[1:40]])
i = 1; j = 2
MODIS2 = crop(MODIS, extent(MODIS, 1+(i-1)*nrow, i*nrow, 1+(j-1)*ncol, j*ncol))
plot(MODIS2[[1:40]])

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


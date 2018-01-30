# MODIS = stack("../data/MOD11A1Day2003_300x300.tif")
# dat = t(values(MODIS))
# saveRDS(dat, "../data/dat_300x300.rds")
dat = readRDS("../data/dat_300x300.rds")
## Test getMask function
msk = getMask(dat)
plot(raster(matrix(msk, 300, 300, byrow=TRUE)))

# tmp = stack("../data/MOD11A1Day2003_50x50.tif")
# brk = brick(tmp) ## in memory
# msk1 = getMask(brk)
# plot(msk1)

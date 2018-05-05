tstset=c( 8,47,56,76,85, 99,117,139,147,170,184,197,221,239,256,282,295,313,327,346)
mskset=c(12,42,58,67,94,110,126,129,154,174,177,199,222,232,252,285,304,306,332,342)

dat = readRDS("../data/MYD11A1Day2010.rds")
r1 = raster(matrix(dat[8,], 1200))
r2 = raster(matrix(dat[12,], 1200))
s = stack(r1,r2)
levelplot(s)

## cover tstset with mskset
dat[tstset, ][is.na(dat[mskset,])] = NA
r1 = raster(matrix(dat[8,], 1200))
r2 = raster(matrix(dat[12,], 1200))
s = stack(r1,r2)
levelplot(s)

saveRDS(dat, "./data/MYD11A1Day2010_simulated.rds")

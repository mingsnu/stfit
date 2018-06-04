library(raster)
tstset=c( 8,47,56,76,85, 99,117,139,147,170,184,197,221,239,256,282,295,313,327,346)
mskset=c(12,42,58,67,94,110,126,129,154,174,177,199,222,232,252,285,304,306,332,342)
dat = readRDS("./data/MYD11A1Day2010.rds")
idat = readRDS("./output/dat_imputed.rds")
evalidx = !is.na(values(raster("./data/bf.h11v04.tif")))
idx = is.na(dat[mskset, evalidx])

tdat = dat[tstset, evalidx]
idat = idat[tstset, evalidx]
d = c(tdat[idx])-c(idat[idx])
sqrt(sum(d^2, na.rm=TRUE)/sum(!is.na(d)))
cor(c(tdat[idx]), c(idat[idx]), use = "complete.obs")
## [1] 267.6731
## [1] 0.9773728

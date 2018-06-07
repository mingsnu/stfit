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
## nc1000
## [1] 268.0228
## [1] 0.977315
## nbasis39 nocluster
## [1] 258.4132
## [1] 0.9788721
## nbasis 27 cluster800
## [1] 277.9573
## [1] 0.975626
## nbasis 27 for first two step and nbasis 11 for the third step; ncluster500
## [1] 269.1594
## [1] 0.9771795


## daily lm impu
## [1] 208.1006
## [1] 0.9855059
## daily shift impu
## [1] 210.5348
## [1] 0.98489

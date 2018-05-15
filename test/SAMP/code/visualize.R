library(raster)
library(rasterVis)
library(stfit)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

smapdat = readRDS("../data/smapdat.rds")
s = mat2stack(smapdat, 406, idx=1:365, xmn=-179.8133, xmx=179.8133, ymn=-83.63197, ymx=83.63197)
levelplot(s[[1:6]], par.settings = colthm, layout=c(2,3))

## mask ===
msk = getMask(smapdat)
levelplot(raster(matrix(msk,406), xmn=-179.8133, xmx=179.8133, ymn=-83.63197, ymx=83.63197))

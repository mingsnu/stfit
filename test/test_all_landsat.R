##### Test for Landsat data
library(feather)
library(dplyr)
df = read_feather("../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
## res = gapfill(year, doy, mat, 31,31,h=0, doyrange=1:90, nnr=5)
res = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=30)

## temporal trend visulization
r.list = list()
for(i in 1:365){
  r.list[[i]] = raster(matrix(res$temporal.mean[i,], 31))
}
s = stack(r.list)
levelplot(s[[seq(1, 365, 10)]], par.settings = RdBuTheme)

s = raster(matrix(res$temporal.mean[1,], 31))
for(i in seq(2,365,length.out = 30)){
  r = raster(matrix(res$temporal.mean[i,], 31))
  s = stack(s,r)
}
levelplot(s, par.settings = RdBuTheme)

## the i-th partial missing images
for(i in 1:length(res$ids$ids.partialmissing)){
  r1 = raster(matrix(mat[res$ids$ids.partialmissing[i],], 31))
  r2 = raster(matrix(res$imputed.partial[i,], 31))
  s = stack(r1, r2)
  print(levelplot(s, par.settings = RdBuTheme))
  Sys.sleep(2)
}

# df1_pm = readRDS("df1_partial_missing.rds")
# partial_imputed = PACE(resid.mat[1:length(idx2),], ev.vec, sigma2, ev.val, mc.cores)

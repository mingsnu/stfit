##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../../data/features_106_wide.feather")
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

## focus on year >= 2000 for test purpose

for(yy in c(1995,2000,2010, 2015)){
  tmpdf = df[df$year == yy,]
  year = tmpdf$year
  doy = tmpdf$doy
  mat = as.matrix(tmpdf[,-c(1:2)])
  mat[mat > 2000] = NA
  r.list = list()
  for(i in 1:nrow(mat)){
    r.list[[i]] = raster(matrix(mat[i,], 31))
  }
  s = stack(r.list)
  pdf(paste0("year", yy, ".pdf"))
  print(levelplot(s, par.settings = colthm, names.attr=as.character(doy), main = paste0("Landsat scene p026r031 Site 106, ", yy),
            layout=c(9,5)))
  dev.off()
}


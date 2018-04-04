##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

##########################################################
#### gapfill based on  pixel level temporal smoothing ####
##########################################################
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
res1 = gapfill(year[idx], doy[idx], mat[idx,], 31,31, h = 0, nnr=5)
# res1 = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=5)
## temporal trend visulization
ssmean1 = mat2stack(res1$temporal.mean, 31)
levelplot(ssmean1[[seq(1, 100, by = 5)]], par.settings = colthm)

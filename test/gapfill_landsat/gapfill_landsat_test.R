##########################################################
######## GAPFILLING PROCEDURE FOR LANDSAT DATA ###########
##########################################################
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(stfit)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

###############################
#### I. Landsat data test #####
###############################
df = read_feather("../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

#### visualize mat before imputation
mat_stack = mat2stack(mat, 31)
levelplot(mat_stack[[seq(1, 365, by = 10)]], par.settings = colthm, at=seq(0, 1500, 100))
registerDoParallel(cores = 8)
stfit::opts$set(temporal_mean_est = stfit::spreg)
mat_imputed = gapfill_landsat(year, doy, mat, 31, 31)
## saveRDS(mat_imputed, "./output/mat_imputed.rds")
mat_imputed_stack = mat2stack(mat_imputed, 31)
levelplot(mat_imputed_stack[[seq(1, 365, by = 10)]], par.settings = colthm, at=seq(0, 1500, 100))


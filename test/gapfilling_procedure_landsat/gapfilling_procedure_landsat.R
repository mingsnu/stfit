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
mat1 = mat
mat1[is.na(mat1)] = 0
mat1_stack0 = mat2stack(mat1, 31)
levelplot(mat1_stack0[[seq(1, 365, by = 20)]], par.settings = colthm, at=seq(0, 1500, 100))

## initialize mat_imputed
imat = mat

###################################
#### 1. Overall mean estimaton ####
###################################
registerDoParallel(cores = 8)
stfit::opts$set(temporal_mean_est = stfit::spreg)
system.time({
  meanest1 = meanEst(doy, mat, doyeval = 1:365)
})

#### comparing results not using redo.
# system.time({
#   meanest2 = meanEst(doy, mat, doyeval = 1:365, redo = FALSE)
# })
# hist(meanest1$meanmat-meanest2$meanmat)
# range(meanest2$meanmat, na.rm = TRUE)
# range(mat, na.rm=TRUE)
# ## mean visulization
# mean_stack = mat2stack(meanest$meanmat, 31)
# levelplot(mean_stack[[seq(1, 365, by = 20)]], par.settings = colthm)

###################################
#### 2. Time effect estimation ####
###################################
## remove outlier pixels
for(i in 1:length(meanest$outlier$outidx)){
  mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
}
## remove outlier images
outlier.img.idx = meanest$idx$idx.outlier
for(i in outlier.img.idx){
  mat[outlier.img.idx,] = NA
}
## calculate the residuals
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

## estimate the temporal effect using residuals
## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
teffarray = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = 100, h.sigma2 = 300)
## saveRDS(teffarray, "./output/teffarray.rds")

######################################
#### 3. Spatial effect estimation ####
######################################
## claculate residuals after removing temporal effect
yearidx = unlist(lapply(year, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
doyidx = unlist(lapply(doy, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
for(i in 1:nrow(rmat)){
  rmat[i,] = rmat[i,] - teffarray[yearidx[i], doyidx[i],]
}
## saveRDS(rmat, "./output/rmat4seff.rds")
## rmat = readRDS("./output/rmat4seff.rds")

#### visualize residials before doing spatial imputation
rmat1 = rmat
rmat1[is.na(rmat1)] = 0
seff_stack0 = mat2stack(rmat1, 31)
levelplot(seff_stack0[[seq(1, 365, by = 20)]], par.settings = colthm)

## estimate the spatial effect using residuals
## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
registerDoParallel(cores = 8)
seffest = seffEst(rmat, 31, 31, nnr = 1)

#### visualize residials after doing spatial imputation
seff_stack = mat2stack(seffest$seffmat, 31)
levelplot(seff_stack[[seq(1, 365, by = 20)]], par.settings = colthm)
seffest = seffEst(rmat, 31, 31, nnr = 15)
seff_stack = mat2stack(seffest$seffmat, 31)
levelplot(seff_stack[[seq(1, 365, by = 20)]], par.settings = colthm)

#######################
#### 4. Gapfilling ####
#######################
## partially missing images: mean + time effect + spatial effect 
## all missing images: mean + time effect
## first calculate the theoretically imputed mat.
mat_imputed = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
for(i in 1:nrow(mat_imputed)){
  mat_imputed[i,] = mat_imputed[i,] + teffarray[yearidx[i], doyidx[i],]
}
mat_imputed = mat_imputed + seffest$seffmat

#### visualize theoretically imputated mat
mat_imputed_stack = mat2stack(mat_imputed, 31)
levelplot(mat_imputed_stack[[seq(1, 365, by = 20)]], par.settings = colthm, at=seq(0, 1500, 100))

#### final imputation
## if keep originally observed values
imat[is.na(imat)] = mat_imputed[is.na(imat)]
## if remove outliers
imat = mat
imat[is.na(imat)] = mat_imputed[is.na(imat)]




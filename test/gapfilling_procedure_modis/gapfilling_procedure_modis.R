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
registerDoParallel(cores = 8)

####################################
#### MODIS data test with mask #####
#####################################
mat = readRDS("../data/dat_300x300.rds")
doy = 1:365
msk = getMask(mat)
levelplot(raster(matrix(msk, 300)), margin = FALSE)
###################################
#### 1. Overall mean estimaton ####
###################################
## when there is no repeated measures, spreg seems to work better
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest = meanEst(1:365, mat, doyeval = 1:365, msk = msk)
# saveRDS(meanest, "./output/MODIS_with_mask_mean_est.rds")
# meanest = readRDS("./output/MODIS_with_mask_mean_est.rds")
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 300)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm)

#####################################################
#### 2. Cluster analysis based on inital results ####
#####################################################
cmat = t(meanest$meanmat[seq(10, 350, by = 10),!msk])
cres = kmeans(cmat, centers = 50, nstart = 10, iter.max = 50)
cluster = rep(0, ncol(mat))
cluster[!msk] = cres$cluster
# saveRDS(cluster, "./output/cluster.rds")
# cluster = readRDS("./output/cluster.rds")
levelplot(raster(matrix(cluster, 300)), margin=FALSE)

################################################
#### 3. Overall mean estimaton with cluster ####
################################################
## NEVER USE smooth_spline when using clusters
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest_cl = meanEst(1:365, mat, doyeval = 1:365, cluster = cluster, msk = msk)
## mean visulization
mean_stack_cl = mat2stack(meanest_cl$meanmat, 300)
levelplot(mean_stack_cl[[seq(10, 350, by = 30)]], par.settings = colthm)
## saveRDS(meanest_cl, "./output/MODIS_with_mask_mean_est_cl.rds")
## meanest_cl = readRDS("./output/MODIS_with_mask_mean_est_cl.rds")
#### visualize mat before imputation
mat1 = mat
mat1[is.na(mat1)] = 0
mat1_stack0 = mat2stack(mat1, 300)
levelplot(mat1_stack0[[seq(1, 365, by = 30)]], par.settings = colthm, at=seq(20000, 33000, 1000))

## initialize mat_imputed
imat = mat

#########################################
#### 2. 'temporal' effect estimation ####
#########################################
## remove outlier images
outlier.img.idx = meanest$idx$idx.outlier
# ## outlier image visulization
# outlier_img_stack = mat2stack(mat[outlier.img.idx,], 300)
# levelplot(outlier_img_stack, par.settings = colthm)

for(i in outlier.img.idx){
  mat[outlier.img.idx,] = NA
}
## remove outlier pixels
for(i in 1:length(meanest_cl$outlier$outidx)){
  mat[meanest_cl$outlier$outidx[i], meanest_cl$outlier$outlst[[i]]] = NA
}

## calculate the residuals
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
## saveRDS(rmat, "./output/rmat4teff.rds")
## rmat = readRDS("./output/rmat4teff.rds")

####### the results are very small.
####### no need to do this step.
# ## estimate the 'temporal' effect using residuals
# ## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
# registerDoParallel(cores = 8)
# cefflist = ceffEst(doy, rmat, cluster,
#                     doyeval = 1:365, h.cov = 100, h.sigma2 = 300, max.per.cluster = 30,
#                     t.grid.num = 50)
# range(cefflist[[1]]$mat)
# ## saveRDS(cefflist, "./output/cefflist.rds")

######################################
#### 3. Spatial effect estimation ####
######################################
# ## claculate residuals after removing temporal effect
# #### visualize residials before doing spatial imputation
seff_stack0 = mat2stack(rmat, 300)
levelplot(seff_stack0[[seq(1, 365, by = 30)]], par.settings = colthm)
# levelplot(seff_stack0[[c(1,21)]], par.settings = colthm)
# cc = seff_stack0[[c(1,21)]]
# levelplot(crop(cc, extent(cc, 1,30,1,30)), par.settings = colthm)

## estimate the spatial effect using residuals
registerDoParallel(cores = 8)
res3.list = foreach(n=1:100) %dopar% {
  ii = floor((n-1)/10) + 1
  jj = (n-1) %% 10 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 300 + cidx
                   })))
  seffEst(rmat[, bIdx], 30, 30, nnr = 6, msk = msk[bIdx])
}
## saveRDS(res3.list, "./output/res3_list.rds")
#### visualize residials after doing spatial imputation
# aa = seffEst(rmat[, bIdx], 30, 30, nnr = 3, msk = msk[bIdx])
# seff_stack1 = mat2stack(aa$seffmat, 30)
# levelplot(seff_stack1[[seq(1, 365, by = 20)]], par.settings = colthm)
# tmp = mat2stack(rmat[, bIdx], 30)
# levelplot(tmp[[seq(1, 365, by = 20)]], par.settings = colthm)

seffmat = matrix(NA, nrow(rmat), ncol(rmat))
for(n in 1:100){
  ii = floor((n-1)/10) + 1
  jj = (n-1) %% 10 + 1
  ## block index
  bIdx = c(t(outer(seq((ii-1)*30+1, ii*30), seq((jj-1)*30+1, jj*30),
                   FUN = function(ridx, cidx){
                     (ridx-1) * 300 + cidx
                   })))
  seffmat[,bIdx] = res3.list[[n]]$seffmat
}
## saveRDS(seffmat, "./output/seffmat.rds")

#### visualize residials after doing spatial imputation
seff_stack = mat2stack(seffmat, 300)
levelplot(seff_stack[[seq(1, 365, by = 30)]], par.settings = colthm)


#######################
#### 4. Gapfilling ####
#######################
## partially missing images: mean + time effect + spatial effect 
## all missing images: mean + time effect
## first calculate the theoretically imputed mat.
mat_imputed = meanest$meanmat
mat_imputed = mat_imputed + seffmat

#### visualize theoretically imputated mat
mat_imputed_stack = mat2stack(mat_imputed, 300)
levelplot(mat_imputed_stack[[seq(1, 365, by = 30)]], par.settings = colthm)

#### final imputation
## if keep originally observed values
imat[is.na(imat)] = mat_imputed[is.na(imat)]
## if remove outliers
imat = mat
imat[is.na(imat)] = mat_imputed[is.na(imat)]




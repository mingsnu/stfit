##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
library(fda)
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
res0 = gapfill(year, doy, mat, 31,31, h = 1, nnr=1, outlier.tol = 0.2)

############################
##### Simulation study #####
############################
res = res0
set.seed(20180124)
n = 6
fidx = sample(res$idx$idx.fullyobserved, n) ## full observed image index
pidx = sample(res$idx$idx.partialmissing, n) ## partial observed image index
## >fidx
## [1] 485 587 107  82 475 599
## > pidx
## [1] 282 311 156 298 590 390

fmat = mat[fidx[1:n], ]
## saveRDS(fmat, "output/fmat.rds")
## apply missing patterns to fully observed images
for(i in 1:n){
  mat[fidx[i],][is.na(mat[pidx[i],])] = NA
}
## artificial partial missing images
pmat = mat[fidx[1:n], ]

##################################################################################
#### MSE based on different different nnr size and temporal smoothing methods ####
##################################################################################
nnr_size = 30
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)

MSE1 = foreach(k = 1:nnr_size, .combine = "rbind") %dopar%{
    ## Gapfill
    res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "lc")
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output_pixel_based_nnr/nnr_", k, "_smooth_spline.rds"))
    ## Calculate MSE
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}

saveRDS(MSE1, paste0("output_pixel_based_nnr/MSE_smooth_spline.rds"))
## v8.6
##              [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
## result.1 24882.47 18531.31 2963.482 8807.417 20024.72 30655.18
## result.2 22453.14 16890.18 3509.072 8174.786 19991.95 26155.96
## result.3 19873.31 15879.84 3441.894 8621.155 19365.10 18854.01

Gapfill::opts$set(temporal_mean_est = llreg)
MSE2 = foreach(k = 1:nnr_size, .combine = "rbind") %dopar%{
    ## Gapfill
    res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "lc")
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output_pixel_based_nnr/nnr_", k, "_llreg.rds"))
    ## Calculate MSE
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}
saveRDS(MSE2, paste0("output_pixel_based_nnr/MSE_llreg.rds"))

Gapfill::opts$set(temporal_mean_est = lpreg)
MSE3 = foreach(k = 1:nnr_size, .combine = "rbind") %dopar%{
    ## Gapfill
    res = gapfill(year, doy, mat, 31,31, h=0, doyrange=1:365, nnr=k, method = "lc")
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output_pixel_based_nnr/nnr_", k, "_lpreg.rds"))
    ## Calculate MSE
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}
saveRDS(MSE3, paste0("output_pixel_based_nnr/MSE_lpreg.rds"))

## sum(MSE1-MSE2<0)/length(MSE1)
## [1] 0.35
## sum(MSE1-MSE3<0)/length(MSE1)
## [1] 0.6611111
## sum(MSE2-MSE3<0)/length(MSE1)
## [1] 0.7166667

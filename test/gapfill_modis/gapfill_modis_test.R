
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
## when there is no repeated measures, spreg seems to work better
## stfit::opts$set(temporal_mean_est = stfit::spreg)
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)

registerDoParallel(cores = 8)
mat_imputed = gapfill_modis(doy, mat, 300, 300, doyeval = 1:365, msk = msk,
         breaks = list(block.nrow = 10, block.ncol = 10, img.nrow = 30, img.ncol = 30),
         outlier.action = "keep",
         teff = FALSE, seff = TRUE)
saveRDS(mat_imputed, "./output/mat_imputed.rds")

mat_imputed_stack = mat2stack(mat_imputed, 300)
levelplot(mat_imputed_stack[[seq(1, 365, by = 30)]], par.settings = colthm)
mat_stack = mat2stack(mat, 300)
levelplot(mat_stack[[seq(1, 365, by = 30)]], par.settings = colthm)


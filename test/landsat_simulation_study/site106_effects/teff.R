library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(stfit)
library(abind)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

#######################################
##### Simulation study for summer #####
#######################################
#### partial missing image indexes with different missing percentage
pmat = readRDS("../missing_pattern/output/missing_pattern.rds")
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(88, 132, 261, 265, 615)
fidx = fidx2
fmat = mat0[fidx, ]
if(!dir.exists("output"))
  dir.create("output")

###### define function used for mean estimation ######
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

## collapse the partial and full index for parallel
## matrix of MxN, column stacking
N = nrow(fmat)
M = nrow(pmat)
registerDoParallel(10)
res = foreach(n = 1:(M*N)) %dopar% {
  i = (n - 1) %% M + 1 ## ROW INDEX
  j = (n - 1) %/% M + 1 ## COLUMN INDEX
  mat = mat0
  ## apply missing patterns to fully observed images
  missing.idx = is.na(pmat[i,])
  mat[fidx[j], missing.idx] = NA
  
  if(file.exists(paste0("./output/P", i, "_F_", j, ".rds"))){
    res1 <- readRDS(paste0("./output/mean_P", i, "_F_", j, ".rds"))
  } else {
    res1 <- gapfill_landsat(year, doy, mat, 31, 31, teff = TRUE, seff = FALSE,
                            use.intermediate.result = TRUE, intermediate.save = TRUE,
                            intermediate.dir = paste0("./intermediate_results/P_", i,
                                                      "_F_", j, "/"))
    saveRDS(res1, paste0("./output/mean_P", i, "_F_", j, ".rds"))
  }
  imat = res1$imat[fidx[j],]
  c(RMSE(fmat[j, missing.idx], imat[missing.idx]),
    NMSE(fmat[j, missing.idx], imat[missing.idx]),
    ARE(fmat[j, missing.idx], imat[missing.idx]),
    cor(fmat[j, missing.idx], imat[missing.idx]))
}
saveRDS(res, "./output/teff_res.rds")


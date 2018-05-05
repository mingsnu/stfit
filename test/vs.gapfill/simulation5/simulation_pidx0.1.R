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

#### partial missing images with different missing percentage
res = readRDS("../our_output/our_res.rds")
missingpct = apply(mat0[res$idx$idx.partialmissing,], 1, function(x) sum(is.na(x))/length(x))
pidx0.1 = res$idx$idx.partialmissing[which(missingpct < 0.1)]
## tmpstack = mat2stack(mat0[pidx0.1,], 31)
## levelplot(tmpstack, par.settings = colthm)
pidx0.1 = pidx0.1[c(1,2,7,8,9)]

#######################################
##### Simulation study for pidx0.1 #####
#######################################

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
Gapfill::opts$set(temporal_mean_est = customfun)


###### variables used for Gapfill package ######
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))

registerDoParallel(10)


######## fully observed image indexes ########
fidx = sort(c(481, 587, 107, 82, 460, 599))
fmat = mat0[fidx, ]

## initialize error metrics matrices
RMSEmat1 =  matrix(NA, length(fidx), length(pidx0.1))
NMSEmat1 = matrix(NA, length(fidx), length(pidx0.1))
AREmat1 = matrix(NA, length(fidx), length(pidx0.1))

for(i in 1:length(pidx0.1)){
  for(j in 1:length(fidx)){
    mat = mat0
    ## apply missing patterns to fully observed images
    missing.idx = is.na(mat[pidx0.1[i],])
    mat[fidx[j], missing.idx] = NA

    #### proposed method
    res1 <- gapfill_landsat(year, doy, mat, 31, 31, teff = TRUE, seff = FALSE,
                            use.intermediate.result = TRUE, intermediate.save = TRUE,
                            intermediate.dir = paste0("./pidx0.1/res1_pidx_", pidx0.1[i],
                                                       "_fidx_", fidx[j], "/"))
    saveRDS(res1, paste0("./pidx0.1/res1_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
    imat = res1$imat[fidx[j],]
    RMSEmat1[j, i] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat1[j, i] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat1[j, i] = ARE(fmat[j, missing.idx], imat[missing.idx])
  }
}

saveRDS(RMSEmat1, "./pidx0.1/RMSEmat1.rds")
saveRDS(NMSEmat1, "./pidx0.1/NMSEmat1.rds")
saveRDS(AREmat1, "./pidx0.1/AREmat1.rds")

RMSEmat1.1 = readRDS("../simulation2/pidx0.1/RMSEmat1.rds") ## without doy.break
RMSEmat1.2 = readRDS("../simulation3/pidx0.1/RMSEmat1.rds") ## with doy.break
RMSEmat1.3 = readRDS("../simulation4/pidx0.1/RMSEmat1.rds") ## overall mean
RMSEmat2 = readRDS("../simulation2/pidx0.1/RMSEmat2.rds")
## overall mean vs. +temporal effect: with temporal effect better
table(c(RMSEmat1 < RMSEmat1.3))
## FALSE  TRUE 
##     7    23
table(c(RMSEmat1 < RMSEmat2)) ## + temporal effect is even better than gapfill package result
## FALSE  TRUE 
##    13    17 
table(c(RMSEmat1 < RMSEmat1.2)) ## adding temporal effect is indeed better
## FALSE  TRUE 
##    20    10
table(c(RMSEmat1 < RMSEmat1.1)) ## adding temporal effect is indeed better
## FALSE  TRUE 
##    23     7
table(c(RMSEmat1.1 < RMSEmat1.2)) ## no obvious advantage using doy.break
## FALSE  TRUE 
##    15    15 


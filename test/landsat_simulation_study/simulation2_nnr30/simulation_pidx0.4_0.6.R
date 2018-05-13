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

#### partial missing image indexes with different missing percentage
pidx0.1 = c(66, 75, 348, 573, 605)
pidx0.4_0.6 = c(74, 156, 273, 285, 326)
pidx0.8_0.95 = c(112, 184, 318, 448, 508)
#### fully observed image indexes from different seasons
fidx1 = c(261, 265, 312, 387, 581)
fidx2 = c(145, 276, 444, 481, 587)
fidx3 = c(198, 202, 493, 549, 557)
fidx4 = c(82, 293, 505, 609, 615)
fidx = c(fidx1, fidx2, fidx3, fidx4)
fmat = mat0[fidx, ]

############################################
##### Simulation study for pidx0.4_0.6 #####
############################################

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

registerDoParallel(16)
## initialize error metrics matrices
RMSEmat1 = RMSEmat2 = matrix(NA, length(fidx), length(pidx0.4_0.6))
NMSEmat1 = NMSEmat2 = matrix(NA, length(fidx), length(pidx0.4_0.6))
AREmat1 = AREmat2 = matrix(NA, length(fidx), length(pidx0.4_0.6))
CORmat1 = CORmat2 = matrix(NA, length(fidx), length(pidx0.4_0.6))

N = length(fidx)
M = length(pidx0.1)
res = foreach(n = 1:(M*N)) %dopar% {
  i = (n - 1) %/% N + 1 ## COLUMN INDEX 
  j = (n - 1) %% N + 1 ## ROW INDEX
    mat = mat0
    ## apply missing patterns to fully observed images
    missing.idx = is.na(mat[pidx0.4_0.6[i],])
    mat[fidx[j], missing.idx] = NA

    #### proposed method
    if(file.exists(paste0("./pidx0.4_0.6/res1_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))){
        res1 <- readRDS(paste0("./pidx0.4_0.6/res1_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
    } else {
        res1 <- gapfill_landsat(year, doy, mat, 31, 31, nnr =30,
                            use.intermediate.result = FALSE, intermediate.save = FALSE)
        saveRDS(res1, paste0("./pidx0.4_0.6/res1_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
    }    
    imat = res1$imat[fidx[j],]
    c(RMSE(fmat[j, missing.idx], imat[missing.idx]),
      NMSE(fmat[j, missing.idx], imat[missing.idx]),
      ARE(fmat[j, missing.idx], imat[missing.idx]),
      cor(fmat[j, missing.idx], imat[missing.idx]))
}

saveRDS(res, "./pidx0.4_0.6/res.rds")
## saveRDS(RMSEmat1, "./pidx0.4_0.6/RMSEmat1.rds")
## saveRDS(NMSEmat1, "./pidx0.4_0.6/NMSEmat1.rds")
## saveRDS(AREmat1, "./pidx0.4_0.6/AREmat1.rds")
## saveRDS(CORmat1, "./pidx0.4_0.6/CORmat1.rds")

## RMSEmat1 = readRDS("./pidx0.4_0.6/RMSEmat1.rds")
## RMSEmat2 = readRDS("./pidx0.4_0.6/RMSEmat2.rds")
## res = readRDS("./pidx0.4_0.6/res.rds")
## RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
## RMSEmat2 = readRDS("../simulation2/pidx0.4_0.6/RMSEmat2.rds")
## RMSEmat1 > RMSEmat2
## table(c(RMSEmat1 < RMSEmat2))
## RMSEmat1.0 = readRDS("../simulation2/pidx0.4_0.6/RMSEmat1.rds")
## RMSEmat1 > RMSEmat1.0

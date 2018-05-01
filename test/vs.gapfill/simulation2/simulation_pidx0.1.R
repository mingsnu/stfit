library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(Gapfill)
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
tmpstack = mat2stack(mat0[pidx0.1,], 31)
levelplot(tmpstack, par.settings = colthm)
pidx0.1 = pidx0.1[c(1,2,7,8,9)]

idx0.4_0.6 = res$idx$idx.partialmissing[which(missingpct <= 0.6 & missingpct > 0.4)]
tmpstack = mat2stack(mat0[idx0.4_0.6,], 31)
levelplot(tmpstack, par.settings = colthm)

idx0.8_0.9 = res$idx$idx.partialmissing[which(missingpct <= 0.9 & missingpct > 0.8)]
tmpstack = mat2stack(mat0[idx0.8_0.9,], 31)
levelplot(tmpstack, par.settings = colthm)

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
RMSEmat1 = RMSEmat2 = matrix(NA, length(fidx), length(pidx0.1))
NMSEmat1 = NMSEmat2 = matrix(NA, length(fidx), length(pidx0.1))
AREmat1 = AREmat2 = matrix(NA, length(fidx), length(pidx0.1))


for(i in 1:length(pidx0.1)){
  for(j in 1:length(fidx)){
    mat = mat0
    ## apply missing patterns to fully observed images
    missing.idx = is.na(mat[pidx0.1[i],])
    mat[fidx[j], missing.idx] = NA

    #### proposed method
    res1 <- gapfill_landsat(year, doy, mat, 31, 31,
                            use.intermediate.result = FALSE, intermediate.save = FALSE)
    saveRDS(res1, paste0("./pidx0.1/res1_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
    imat = res1$imat[fidx[j],]
    RMSEmat1[j, i] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat1[j, i] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat1[j, i] = ARE(fmat[j, missing.idx], imat[missing.idx])
    
    #### Gapfill package
    datarray = array(NA, dim = c(31, 31, 46, 16), dimnames = list(1:31, 1:31, doybinuni, yearuni))
    for(ii in 1:16){
      for(jj in 1:46){
        idx = year == yearuni[ii] & doybin == doybinuni[jj]
        if(sum(idx) == 1)
          datarray[,,jj,ii] = matrix(mat[year == yearuni[ii] & doybin == doybinuni[jj],], 31) else
            if(sum(idx) > 1)
              warning("Multiple matches.")
      }
    }
    
    yidx = which(year[fidx[j]] == yearuni)
    didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
    didxinterval = max(1,didx-6):min(46, didx + 6)
    yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
    tmpmat = datarray[,,didxinterval, yidxinterval]
    ## gapfill::Image(tmpmat)
    res2 = gapfill::Gapfill(tmpmat, clipRange = c(0, 1800), dopar = TRUE)
    saveRDS(res2, paste0("./pidx0.1/res2_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
    imat = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    
    RMSEmat2[j, i] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat2[j, i] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat2[j, i] = ARE(fmat[j, missing.idx], imat[missing.idx])
  }
}

saveRDS(RMSEmat1, "./pidx0.1/RMSEmat1.rds")
saveRDS(RMSEmat2, "./pidx0.1/RMSEmat2.rds")
saveRDS(NMSEmat1, "./pidx0.1/NMSEmat1.rds")
saveRDS(NMSEmat2, "./pidx0.1/NMSEmat2.rds")
saveRDS(AREmat1, "./pidx0.1/AREmat1.rds")
saveRDS(AREmat2, "./pidx0.1/AREmat2.rds")
## RMSEmat1 > RMSEmat2
##       [,1]  [,2]  [,3]  [,4]  [,5]
## [1,] FALSE FALSE  TRUE FALSE FALSE
## [2,]  TRUE FALSE  TRUE FALSE FALSE
## [3,] FALSE FALSE  TRUE FALSE FALSE
## [4,] FALSE FALSE FALSE FALSE FALSE
## [5,] FALSE  TRUE  TRUE FALSE  TRUE
## [6,]  TRUE  TRUE  TRUE FALSE  TRUE
## > table(c(RMSEmat1 < RMSEmat2))

## FALSE  TRUE 
##    11    19 

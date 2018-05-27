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
##### Simulation study for spring #####
#######################################
#### partial missing image indexes with different missing percentage
pmat = readRDS("../missing_pattern/output/missing_pattern.rds")
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(132, 261, 265, 615, 338) 
fidx = fidx1
fmat = mat0[fidx, ]
if(!dir.exists("gapfill_spring"))
  dir.create("gapfill_spring")


###### variables used for Gapfill package ######
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))

## collapse the partial and full index for parallel
## matrix of MxN, column stacking
N = nrow(fmat)
M = nrow(pmat)
registerDoParallel(16)
res = foreach(n = 1:(M*N)) %dopar% {
  i = (n - 1) %% M + 1 ## ROW INDEX
  j = (n - 1) %/% M + 1 ## COLUMN INDEX
  mat = mat0
  ## apply missing patterns to fully observed images
  missing.idx = is.na(pmat[i,])
  mat[fidx[j], missing.idx] = NA
  
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
  
  if(file.exists(paste0("./gapfill_spring/gapfill_spring_P", i, "_F_", j, ".rds"))){
    res1 <- readRDS(paste0("./gapfill_spring/gapfill_spring_P", i, "_F_", j, ".rds"))
  } else {
    res1 = gapfill::Gapfill(tmpmat, clipRange = c(0, 1800), dopar = TRUE)
    saveRDS(res1, paste0("./gapfill_spring/gapfill_spring_P", i, "_F_", j, ".rds"))
  }
  imat = c(res1$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
  c(RMSE(fmat[j, missing.idx], imat[missing.idx]),
    NMSE(fmat[j, missing.idx], imat[missing.idx]),
    ARE(fmat[j, missing.idx], imat[missing.idx]),
    cor(fmat[j, missing.idx], imat[missing.idx]))
}
saveRDS(res, "./gapfill_spring/res.rds")

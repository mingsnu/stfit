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

df = read_feather("../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA


##################################
#### gapfill::Gapfill method #####
##################################
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))
datarray = array(NA, dim = c(31, 31, 46, 16), dimnames = list(1:31, 1:31, doybinuni, yearuni))
for(i in 1:16){
  for(j in 1:46){
    idx = year == yearuni[i] & doybin == doybinuni[j]
    if(sum(idx) == 1)
      datarray[,,j,i] = matrix(mat[year == yearuni[i] & doybin == doybinuni[j],], 31) else
        if(sum(idx) > 1)
          warning("Multiple matches.")
  }
}
str(datarray)
gapfill::Image(datarray[,,1:16, 1:8])

system.time({
  registerDoParallel(8)
  gapfill_res <- gapfill::Gapfill(datarray, clipRange = c(0, 2000), dopar = TRUE)
})
gapfill::Image(gapfill_res$fill[,,1:16, 1:8])
saveRDS(gapfill_res, "./gapfill_output/gapfill_res.rds")
# data has 707296 values: 197625 (28%) observed
# 509671 (72%) missing
# 509671 (72%) to predict
# started at 2018-04-26 09:58:02.
# elapsed time is 22.754 mins (0.00268 secs per NA).
# user   system  elapsed 
# 3154.042   87.175 1365.242 

##################################
#### Gapfill::gapfill method #####
##################################
system.time({
  registerDoParallel(8)
  our_res <- gapfill_landsat(year, doy, mat, 31, 31, intermediate.dir = "./our_output/",
                             intermediate.save = FALSE, use.intermediate.result = FALSE)
})
# user   system  elapsed 
# 2465.583   95.854  355.525 
saveRDS(our_res, "./our_output/our_res.rds")

############################
##### Simulation study #####
############################
res = our_res
set.seed(20180124)
n = 6
fidx = sample(res$idx$idx.fullyobserved, n) ## full observed image index
pidx = sample(res$idx$idx.partialmissing, n) ## partial observed image index
# > fidx
# [1] 481 587 107  82 460 599
# > pidx
# [1] 253 287 148 650 579 366

fmat = mat[fidx[1:n], ]
## saveRDS(fmat, "simulation1/fmat.rds")
## apply missing patterns to fully observed images
for(i in 1:n){
  mat[fidx[i],][is.na(mat[pidx[i],])] = NA
}
## artificial partial missing images
pmat = mat[fidx[1:n], ]

##################################################################################
#### MSE based on different different nnr size and temporal smoothing methods ####
##################################################################################
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

####### Gapfill with no temporal effect
system.time({
  registerDoParallel(8)
  our_res1 <- gapfill_landsat(year, doy, mat, 31, 31, intermediate.dir = "./simulation1/",
                              teff = FALSE)
})
# user  system elapsed 
# 75.680   3.745  34.746 
saveRDS(our_res1, "./simulation1/our_res1.rds")
imat1 = our_res1$imat[fidx[1:n],]
apply((fmat - imat1)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
# [1]  9401.989 48888.104  2401.824 10851.620 38626.121  3129.779

####### Gapfill with temporal effect
system.time({
  registerDoParallel(8)
  our_res2 <- gapfill_landsat(year, doy, mat, 31, 31, intermediate.dir = "./simulation1/",
                              teff = TRUE)
})
# user   system  elapsed 
# 2849.805  352.474  464.637
saveRDS(our_res2, "./simulation1/our_res2.rds")
imat2 = our_res2$imat[fidx[1:n],]
apply((fmat - imat2)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
# [1]  6793.560 37849.273  2475.201  8611.613 33930.876  2282.639
# 
# mat_imputed_stack = mat2stack(mat_imputed, 31)
# levelplot(mat_imputed_stack[[seq(1, 365, by = 10)]], par.settings = colthm, at=seq(0, 1500, 100))
# 
i=3
levelplot(raster(matrix(fmat[i,], 31)))
levelplot(raster(matrix(mat[fidx[i],], 31)))
levelplot(raster(matrix(imat2[i,], 31)))


############ gapfill ############
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))
datarray = array(NA, dim = c(31, 31, 46, 16), dimnames = list(1:31, 1:31, doybinuni, yearuni))

for(i in 1:16){
  for(j in 1:46){
    idx = year == yearuni[i] & doybin == doybinuni[j]
    if(sum(idx) == 1)
      datarray[,,j,i] = matrix(mat[year == yearuni[i] & doybin == doybinuni[j],], 31) else
        if(sum(idx) > 1)
          warning("Multiple matches.")
  }
}

system.time({
  registerDoParallel(8)
  gapfill_res1 <- gapfill::Gapfill(datarray, clipRange = c(0, 2000), dopar = TRUE)
})
saveRDS(gapfill_res1, "./simulation1/gapfill_res1.rds")

MSE1 = rep(0, n)
gapfill_res1 = list()
registerDoParallel(8)
for(i in 1:n){
  yidx = which(year[fidx[i]] == yearuni)
  didx = which(findInterval(doy[fidx[i]], seq(1,365, by=8)) == doybinuni)
  ## levelplot(raster(datarray[,,didx,yidx]))
  didxinterval = max(1,didx-4):min(46, didx + 4)
  yidxinterval = max(1, yidx - 3):min(16, yidx + 3)
  tmpmat = datarray[,,didxinterval, yidxinterval]
  ## gapfill::Image(tmpmat)
  gapfill_res1[[i]] = gapfill::Gapfill(tmpmat, clipRange = c(0, 2000), dopar = TRUE)
  MSE1[i] = sum((c(gapfill_res1[[i]]$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)]) - fmat[i,])^2) / sum(is.na(pmat[i,]))
  ## levelplot(raster(gapfill_res1[[i]]$fill[,,5,4]))
}
saveRDS(gapfill_res1, "./simulation1/gapfill_res1.rds")
## > MSE1
## [1]  9772.7309 17885.3934  4405.5507 11604.0390         NA   930.0504


MSE2 = rep(0, n)
gapfill_res2 = list()
registerDoParallel(8)
for(i in 1:n){
  yidx = which(year[fidx[i]] == yearuni)
  didx = which(findInterval(doy[fidx[i]], seq(1,365, by=8)) == doybinuni)
  ## levelplot(raster(datarray[,,didx,yidx]))
  didxinterval = max(1,didx-5):min(46, didx + 5)
  yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
  tmpmat = datarray[,,didxinterval, yidxinterval]
  ## gapfill::Image(tmpmat)
  gapfill_res2[[i]] = gapfill::Gapfill(tmpmat, clipRange = c(0, 2000), dopar = TRUE)
  MSE2[i] = sum((c(gapfill_res2[[i]]$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)]) - fmat[i,])^2) / sum(is.na(pmat[i,]))
  ## levelplot(raster(gapfill_res2[[i]]$fill[,,5,4]))
}
saveRDS(gapfill_res2, "./simulation1/gapfill_res2.rds")
# > MSE2
# [1]  9254.771 16482.933  3370.354 16233.064        NA  1008.186

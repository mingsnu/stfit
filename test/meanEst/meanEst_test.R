##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
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

###################################
#### 1. Overall mean estimaton ####
###################################
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)

#### can also use customized function. Compared with spreg function, the customfun is faster
#### since it avoids repeated evaluating the design matrix.
# .X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
# customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
#         nonna.idx = !is.na(y)
#         if(sum(nonna.idx) < minimum.num.obs)
#           return(rep(NA, 365))
#         ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
#         lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
#         return(.X %*% lmfit$coefficient)
# }
# Gapfill::opts$set(temporal_mean_est = customfun)

meanest = meanEst(doy, mat, doyeval = 1:365)
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 31)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm, at=seq(0, 1600, 100))

## outlier images
levelplot(stack(raster(matrix(mat[meanest$idx$idx.outlier[15],], 31)),
                raster(matrix(meanest$meanmat[doy[meanest$idx$idx.outlier[15]],], 31))))

## mean estimation by specifying clipRange and clipMethod.
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0, 1500), clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 31)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm, at=seq(0, 1600, 100))

#####################################################
#### 2. Cluster analysis based on inital results ####
#####################################################
cmat = t(meanest$meanmat[seq(10, 250, by = 10),])
cres = kmeans(cmat, centers = 10, nstart = 10)
levelplot(raster(matrix(cres$cluster, 31)), margin=FALSE)

################################################
#### 3. Overall mean estimaton with cluster ####
################################################
## NEVER USE smooth_spline when using clusters
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest_cl = meanEst(doy, mat, doyeval = 1:365, cluster = cres$cluster)
## mean visulization
mean_stack_cl = mat2stack(meanest_cl$meanmat, 31)
levelplot(mean_stack_cl[[seq(1, 100, by = 5)]], par.settings = colthm)

##############################
#### II. MODIS data test #####
##############################
dat = stack("../data/MOD11A1Day2003_50x50.tif")
dat = t(values(dat))

###################################
#### 1. Overall mean estimaton ####
###################################
## when there is no repeated measures, spreg seems to work better
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest = meanEst(1:365, dat, doyeval = 1:365)
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 50)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm)

#####################################################
#### 2. Cluster analysis based on inital results ####
#####################################################
cmat = t(meanest$meanmat[seq(10, 250, by = 10),])
cres = kmeans(cmat, centers = 10, nstart = 20, iter.max = 50)
levelplot(raster(matrix(cres$cluster, 50)), margin=FALSE)


#########################################
#### III. MODIS data test with mask #####
#########################################
dat = readRDS("../data/dat_300x300.rds")
msk = getMask(dat)
levelplot(raster(matrix(msk, 300)), margin = FALSE)
###################################
#### 1. Overall mean estimaton ####
###################################
## when there is no repeated measures, spreg seems to work better
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest = meanEst(1:365, dat, doyeval = 1:365, msk = msk)
## mean visulization
mean_stack = mat2stack(meanest$meanmat, 300)
levelplot(mean_stack[[seq(1, 100, by = 5)]], par.settings = colthm)

#####################################################
#### 2. Cluster analysis based on inital results ####
#####################################################
cmat = t(meanest$meanmat[seq(10, 350, by = 10),!msk])
cres = kmeans(cmat, centers = 50, nstart = 10, iter.max = 50)
cluster = rep(0, ncol(dat))
cluster[!msk] = cres$cluster
levelplot(raster(matrix(cluster, 300)), margin=FALSE)

################################################
#### 3. Overall mean estimaton with cluster ####
################################################
## NEVER USE smooth_spline when using clusters
Gapfill::opts$set(temporal_mean_est = Gapfill::spreg)
meanest_cl = meanEst(1:365, dat, doyeval = 1:365, cluster = cluster, msk = msk)
## mean visulization
mean_stack_cl = mat2stack(meanest_cl$meanmat, 300)
levelplot(mean_stack_cl[[seq(1, 100, by = 20)]], par.settings = colthm)




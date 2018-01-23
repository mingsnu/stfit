##### Test for MODIS data
library(raster)
library(rasterVis)
library(devtools)
library(refund)
load_all()

files = list.files("../data/MODIS2003", "*.tif$", full.names = TRUE)
## brick.list contains a list of bricks at different 
brick.list = list()
for(i in 1:length(files)){
  brick.list[[i]] = brick(files[i])
}
brick.list

#############################################
###### Temporal difference estimation #######
#############################################
## the difference between brick one and brick three
br13 = brick.list[[1]] - brick.list[[3]]
br13 = rmOutlier(br13)
## smooth the temporal difference
br13.tm = doyMeanEst(br13)

#########################################################
###### Impute image one with three and vice versa #######
#########################################################
br1 = cover(brick.list[[1]], brick.list[[3]] + br13.tm)
br3 = cover(brick.list[[3]], brick.list[[1]] - br13.tm)

cols <- terrain.colors(30) # heat.colors(20)
rnames = gsub("MOD11A1Day2003_100x100.", "", names(br1))
levelplot(br1[[1:21]], col.regions = cols,
          main="MOD11A1Day2003_100x100",
          names.attr = rnames[1:21])
cols <- terrain.colors(30) # heat.colors(20)
rnames = gsub("MOD11A1Day2003_100x100.", "", names(br3))
levelplot(br3[[1:21]], col.regions = cols,
          main="MOD11A1Day2003_100x100",
          names.attr = rnames[1:21])

#############################################################
###### Combine br1 and br3 for the following analysis #######
#############################################################
df = as.data.frame(rbind(t(values(br1)), t(values(br3))))
colnames(df) = paste0("pixel", 1:10000)
df$year = rep(2003, nrow(df))
df$doy = rep(1:365, 2)

#######################################
###### Temporal mean estimation #######
#######################################
## 1. Remove outliers
brick.list = lapply(brick.list, function(rst) rmOutlier(rst))
## 2. Temporal mean estimation (smoothing on temporal domain)
brick.list.tm <- lapply(brick.list, function(rst) doyMeanEst(rst))

##################################
###### Residual estimation #######
##################################
# br1.resid = br1 - brick.list.tm[[1]]
# br3.resid = br3 - brick.list.tm[[3]]
br1.resid = brick.list[[1]] - brick.list.tm[[1]]
br3.resid = brick.list[[3]] - brick.list.tm[[3]]

#######################################
###### Residual mean estimation #######
#######################################

#############################################
###### Residual covariance estimation #######
#############################################
##################################
###### Residual prediction #######
##################################
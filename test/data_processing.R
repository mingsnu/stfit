#######################################
###### Ladsat: long data to wide ######
#######################################
## the wide version uses less storage space and easier to program
library(feather)
df <- read_feather("../data/features_106.feather")
colnames(df) = c("lat", "lon", paste0("f", 1:11), "year", "doy", "pixel")
df[df == -9999] = NA
## f11 is Fmask: 0-clear land; 1-clear water; 2- cloud shadow; 
## 3-snow; 4-cloud; 255-no observation;
## Define missing values for `f1` by using `f11`
df$f1[df$f11 %in% c(2,3,4,255)] = NA

df = df[, c("year", "doy", "pixel", "f1")]
widedf = dcast(df, year + doy ~ pixel, value.var="f1")
colnames(widedf) = c("year", "doy", paste0("pixel", 1:961))
write_feather(widedf, "../data/features_106_wide.feather")

#######################################
###### MODIS: 300x300 image crop ######
#######################################
## save 300x300 images for use later. This images contain constatnt missing part
MODIS = stack("../data/MOD11A1Day2003.tif")
MODIS3 = crop(MODIS, extent(MODIS, 601, 900, 601, 900))
plot(MODIS3[[1:9]])
writeRaster(MODIS3, "../data/MOD11A1Day2003_300x300.tif")
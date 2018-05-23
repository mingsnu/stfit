library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)
df = read_feather("../data/features_106_wide.feather")
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

##########################################
#### fully observed images selection #####
##########################################
##### Four seasons ====================================
##### Astronomical
## Spring: March 20 ~ June 20:  DOY 80 ~ 171
## Summer: June 21 ~ Sep 22: 172 ~ 265
## Fall: Sep 23 ~ Dec 20: 266 ~ 354
## Winter: Dec 21 ~ March 20: DOY 1 ~ 79; 354 ~ 365
##### Meteorological season
## Spring:March 1 ~ May 31: DOY 60 ~ 151
## Summer: June 1 ~ August 31: DOY 152 ~ 243
## Fall: Sep 1 ~ Nov 30: DOY 244 ~ 334
## Winter: Dec 1 ~ Feb 28: DOY 335 ~ 365; 1 ~ 59
full.idx = which(apply(mat, 1, function(x) all(!is.na(x))))
spring.idx = full.idx[which(doy[full.idx] %in% 60:151)]
summer.idx = full.idx[which(doy[full.idx] %in% 152:243)]
fall.idx = full.idx[which(doy[full.idx] %in% 244:334)]
winter.idx = full.idx[which(doy[full.idx] %in% c(1:59, 335:365))]
spring.idx
summer.idx
fall.idx
winter.idx

##### selected fully observed images ==================
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(88, 132, 261, 265, 615) 

##### image information talbe =========================
df = data.frame(year = year[c(fidx1, fidx2, fidx3, fidx4)],
           DOY = doy[c(fidx1, fidx2, fidx3, fidx4)],
           idx = c(fidx1, fidx2, fidx3, fidx4))
print(xtable::xtable(df[,1:2]), include.rownames=FALSE)

##### Plot ============================================
pdf("./plot/site106_fully_observed.pdf", width=6.9, height=5.5)
s = mat2stack(mat[c(fidx1, fidx2, fidx3, fidx4),], 31)
levelplot(s, names.attr = c(paste0("F", 1:20)), par.settings = colthm, layout = c(5,4))
dev.off()




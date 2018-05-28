library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)
df = read_feather("../../data/features_2_wide.feather")
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 3000] = NA

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

levelplot(mat2stack(mat[spring.idx,],31), at = seq(0,1800, length=20))
year[spring.idx[c(2, 10, 18, 26, 31)]]
doy[spring.idx[c(2, 10, 18, 26, 31)]]
findInterval(doy[spring.idx[c(2, 10, 18, 26, 31)]], seq(1,365, by=8))

levelplot(mat2stack(mat[summer.idx,],31), at = seq(0,1800, length=20))
year[summer.idx[c(2, 12, 18, 25, 32)]]
doy[summer.idx[c(2, 12, 18, 25, 32)]]
findInterval(doy[summer.idx[c(2, 12, 18, 25, 32)]], seq(1,365, by=8))

levelplot(mat2stack(mat[fall.idx,],31), at = seq(0,1800, length=20))
year[fall.idx[c(2, 12, 19, 26, 32)]]
doy[fall.idx[c(2, 12, 19, 26, 32)]]
findInterval(doy[fall.idx[c(2, 10, 18, 26, 32)]], seq(1,365, by=8))

levelplot(mat2stack(mat[winter.idx,],31), at = seq(0,1800, length=20))
year[winter.idx]
doy[winter.idx]
findInterval(doy[winter.idx], seq(1,365, by=8))

##### selected fully observed images ==================
fidx1 = c(13, 101, 267, 432, 485)
fidx2 = c(21, 110, 192, 280, 493)
fidx3 = c(33, 121, 295, 458, 563) #609 
fidx4 = c(95, 128, 222, 261) 

##### image information talbe =========================
df = data.frame(year = year[c(fidx1, fidx2, fidx3, fidx4)],
                DOY = doy[c(fidx1, fidx2, fidx3, fidx4)],
                idx = c(fidx1, fidx2, fidx3, fidx4))
print(xtable::xtable(df[,1:2]), include.rownames=FALSE)

##### Plot ============================================
pdf("./plot/site2_fully_observed.pdf", width=6.9, height=5.5)
s = mat2stack(mat[c(fidx1, fidx2, fidx3, fidx4),], 31)
levelplot(s, names.attr = c(paste0("F", 1:19)), par.settings = colthm, layout = c(5,4))
dev.off()




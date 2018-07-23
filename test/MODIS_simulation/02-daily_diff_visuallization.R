library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)

dat1 = readRDS("./data/MYD11A1Day2010_simulated.rds")
dat2 = readRDS("../data/MOD11A1Day2010.rds")

msk = getMask(dat2)

## Only use 1 and 3 for test purpose
#############################################
###### Temporal difference estimation #######
#############################################
## the difference between brick one and brick three
dat12diff = dat1 - dat2

pdf("./output/datadiff_visualization.pdf")
par(mfrow=c(3,1))
for(i in sample(1:ncol(dat12diff), 300)){
  if(!msk[i]){
    plot(dat12diff[,i], ylim = c(-1000, 1000))
  }
}
dev.off()
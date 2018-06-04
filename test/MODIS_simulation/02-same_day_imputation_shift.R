library(dplyr)
library(doParallel)
library(Matrix)
library(stfit)
library(foreach)

dat1 = readRDS("./data/MYD11A1Day2010_simulated.rds")
dat2 = readRDS("./data/MOD11A1Day2010.rds")

if(!file.exists("./data/msk.rds")){
  msk = getMask(dat2)
  saveRDS(msk, "./data/msk.rds")
}

diff_median = sapply(1:365, function(i) median(dat1[i,] - dat2[i,], na.rm=TRUE))
idx = is.na(dat1)
dat1[idx] = (dat2+diff_median)[idx]
saveRDS(dat1, "./output/MYD11A1Day2010_simulated_daily_imputed_shift.rds")

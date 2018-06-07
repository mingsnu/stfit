library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)

dat1 = readRDS("./data/MYD11A1Day2010_simulated.rds")
dat2 = readRDS("../data/MOD11A1Day2010.rds")

## Only use 1 and 3 for test purpose
#############################################
###### Temporal difference estimation #######
#############################################
## the difference between brick one and brick three
dat12diff = dat1 - dat2
boxstat = boxplot(c(dat12diff), plot = FALSE)
dat12diff[dat12diff > boxstat$stats[5,] | dat12diff < boxstat$stats[1,]] = NA
## saveRDS(dat12diff, "./output/dat12diff.rds")
## dat12diff = readRDS("./output/dat12diff.rds")

#########################################
#### whole image into 12 x 12 blocks ####
#########################################
registerDoParallel(cores=5)
midx = rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31))

diff_month_median = foreach(i=1:12, .combine = rbind) %dopar%{
  foreach(n = 1:144, .combine=c) %do%{
    ii = floor((n-1)/12) + 1
    jj = (n-1) %% 12 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*100+1, ii*100), seq((jj-1)*100+1, jj*100),
                     FUN = function(ridx, cidx){
                       (ridx-1) * 1200 + cidx
                     })))
    median(dat12diff[midx==i, bIdx], na.rm=TRUE)
  }
}


## every Ten Day Index
tdidx = rep(1:36, c(rep(11, 3), rep(10, 31), rep(11, 2)))
diff_td_median = foreach(i=1:36, .combine = rbind) %dopar%{
  foreach(n = 1:144, .combine=c) %do%{
    ii = floor((n-1)/12) + 1
    jj = (n-1) %% 12 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*100+1, ii*100), seq((jj-1)*100+1, jj*100),
                     FUN = function(ridx, cidx){
                       (ridx-1) * 1200 + cidx
                     })))
    median(dat12diff[tdidx==i, bIdx], na.rm=TRUE)
  }
}
td_na_idx = which(is.na(diff_td_median), arr.ind = TRUE)
## impute missing valus of ten days median with monthly median.
for(i in 1:nrow(td_na_idx)){
  diff_td_median[td_na_idx[i,1], td_na_idx[i,2]] =
    diff_month_median[findInterval(td_na_idx[i,1]*10, cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31)), rightmost.closed = TRUE) + 1, td_na_idx[i,2]]
}


## every Day Index
didx = 1:365
diff_day_median = foreach(i=1:365, .combine = rbind) %dopar%{
  foreach(n = 1:144, .combine=c) %do%{
    ii = floor((n-1)/12) + 1
    jj = (n-1) %% 12 + 1
    ## block index
    bIdx = c(t(outer(seq((ii-1)*100+1, ii*100), seq((jj-1)*100+1, jj*100),
                     FUN = function(ridx, cidx){
                       (ridx-1) * 1200 + cidx
                     })))
    median(dat12diff[didx==i, bIdx], na.rm=TRUE)
  }
}
day_na_idx = which(is.na(diff_day_median), arr.ind = TRUE)

## impute missing valus of everyday median with ten days median.
for(i in 1:nrow(day_na_idx)){
  diff_day_median[day_na_idx[i,1], day_na_idx[i,2]] =
    diff_td_median[findInterval(day_na_idx[i,1], cumsum(c(rep(11, 3), rep(10, 31), rep(11, 2))),
                                rightmost.closed = TRUE) + 1, day_na_idx[i,2]]
}
saveRDS(diff_day_median, "./output/diff_day_median.rds")

####################################
###### Impute with each other ######
####################################
whichblock <- function(x){
  ii = (x-1)%% 1200 + 1
  jj = (x-1) %/% 1200 + 1
  ((jj -1)%/% 100 + 1) + ((ii -1) %/% 100)*12
}

registerDoParallel(cores = 6)
dat1_imputed = foreach(i = 1:365, .combine=rbind) %dopar%{
  dat2[i,] + c(outer(1:1200, 1:1200, function(jj,ii) diff_day_median[i, ((jj-1) %/% 100 + 1) + ((ii - 1) %/% 100)*12]))
}

tmpidx = which(is.na(dat1))
dat1[tmpidx] = dat1_imputed[tmpidx]
saveRDS(dat1, "./output/dat1.rds")

## ###################################
## #### scatterplot dat1 vs. dat2 ####
## ###################################
## dat0= readRDS("./data/MYD11A1Day2010_simulated.rds")

## i=120000
## lm.fit = lm(dat0[,i]~dat2[,i])
## summary(lm.fit)
## ypred =  predict(lm.fit, newdata = data.frame(x=dat2[,i]))
## cor(dat0[,i], ypred, use="complete.obs")
## cor(dat0[,i], dat1_imputed[,i], use="complete.obs")

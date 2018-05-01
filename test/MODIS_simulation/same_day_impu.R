##### Test for MODIS data
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)

## if(!file.exists("../../data/MYD11A1Day2010.rds")){
## n  MODIS = stack("../../../MODIS_LST_new/2010/MYD11A1Day2010.tif")
##   dat0 = t(values(MODIS))
##   saveRDS(dat0, "../data/MYD11A1Day2010.rds")
## }
dat1 = readRDS("./data/MYD11A1Day2010_simulated.rds")
dat2 = readRDS("../data/MYD11A1Day2010.rds")

msk = getMask(dat1)

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

## diff mean estimation
diff.median = apply(dat12diff, 1, median,  na.rm = TRUE)
registerDoParallel(cores = 18)
M = nrow(dat12diff)
bs = fda::create.fourier.basis(rangeval=c(0,365), nbasis=11)
X = fda::eval.basis(1:365, bs)

system.time({
  diff.mat = foreach(i = 1:ncol(dat12diff)) %dopar% {
    if(msk[i])
      return(rep(NA, M)) else{
        y = dat12diff[,i]
        nonna.idx = !is.na(y)
        if(sum(nonna.idx) < 10)
          return(diff.median)
        lmfit = lm.fit(X[nonna.idx,], y[nonna.idx])
        return(X %*% lmfit$coefficient)
      }
  }
})
##     user   system  elapsed 
## 9412.599   64.403 9288.976

tmp = apply(diff.mat, 2, min, na.rm=TRUE)
idx = which(tmp < -1000)
pdf("./output/outlier1439961.pdf")
par(mfrow=c(2,1))
x.eval = 1:365
for(i in idx){
  y = dat12diff[,i]
  tmp1 = X %*% lmfit$coefficient
  tmp2 = smooth_spline(x.eval[nonna.idx], y[nonna.idx], x.eval)
  plot(x.eval[nonna.idx], y[nonna.idx], xlim = c(1,365), ylim = range(c(y, tmp1, tmp2), na.rm=TRUE))
  lines(x.eval, tmp1)
  lines(x.eval, tmp2)
}
dev.off()

rangebydoy = apply(dat12diff, 1, range, na.rm=TRUE)
## MOD11A1Day2010.1   -652  275
## MOD11A1Day2010.2   -653  275
## MOD11A1Day2010.3   -653  275
## MOD11A1Day2010.4    Inf -Inf
## MOD11A1Day2010.5   -653  274
## MOD11A1Day2010.361 -653  275
## MOD11A1Day2010.362 -651  275
## MOD11A1Day2010.363 -653  273
## MOD11A1Day2010.364  -72  274
## MOD11A1Day2010.365 -646  273
range(dat1, na.rm=TRUE)
## [1] 23724 32746

## mean.mat: columns are pixel index, rows are doy index (ex. 365 x 961)
diff.mat = do.call("cbind", diff.mat)
## saveRDS(diff.mat, "./output/diff_mat.rds")
## diff.mat = readRDS("./output/diff_mat.rds")

pdf("random30lines.pdf")
idx1 = seq(1, ncol(dat12diff), length = 30)
plot(0, xlim=c(1,365), ylim= range(dat12diff[,idx1], na.rm=TRUE), type = "n")
for(i in idx1){
  yy = dat12diff[,i]
  yy.idx = !is.na(yy)
  lines(x.eval[yy.idx], yy[yy.idx], col=rgb(0, 0, 1, 0.5))
}
dev.off()


##################################
##### median diff by month #######
##################################
midx = rep(1:12, c(31,28,31,30,31,30,31,31,30,31,30,31))
mdiff.mat = matrix(0, 12, ncol(dat12diff))

for(i in 1:12){
  mdiff.mat[i,] = apply(dat12diff[midx==i,],2, median, na.rm = TRUE)
}


registerDoParallel(cores = 12)
mdiff.mat = foreach(i=1:12, .combine = rbind) %dopar%{
  apply(dat12diff[midx==i,],2, median, na.rm = TRUE)
}

###################################
#### scatterplot dat1 vs. dat2 ####
###################################
pdf("./output/scatterplot_dat1_dat2.pdf")
pidx = sample(1:ncol(dat1), 50)
for(i in pidx){
  if(sum(!is.na(dat1[,i] - dat2[,i])) >0)
    plot(dat1[,i],dat2[,i])
}
dev.off()


i=2
lm.fit = lm(dat1[,i]~dat2[,i])
summary(lm.fit)
ypred =  predict(lm.fit, newdata = data.frame(x=dat2[,i]))
range(dat2[,i] - dat1[,i], na.rm=TRUE)
range(ypred  - dat1[,i], na.rm=TRUE)
## [1] -653.2162  804.4651
sd(ypred - dat1[,i], na.rm=TRUE)
## [1] 253.2291
sd(dat2[,i] - dat1[,i], na.rm=TRUE)

medianbydoy = apply(dat12diff, 1, median, na.rm=TRUE)
lm.fit1 = lm(dat1[,i]~dat2[,i] + medianbydoy)
summary(lm.fit1)
ypred1 =  predict(lm.fit1, newdata = data.frame(x=dat2[,i]))
sd(ypred1 - dat1[,i], na.rm=TRUE)
## [1] 250.4942
range(ypred1  - dat1[,i], na.rm=TRUE)
## [1] -629.6109  739.7980


#########################################
#### whole image into 12 x 12 blocks ####
#########################################
registerDoParallel(cores=16)
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

registerDoParallel(cores = 14)
dat1_imputed = foreach(i = 1:365, .combine=rbind) %dopar%{
  dat2[i,] + c(outer(1:1200, 1:1200, function(jj,ii) diff_day_median[i, ((jj-1) %/% 100 + 1) + ((ii - 1) %/% 100)*12]))
}

dat2_imputed = foreach(i = 1:365, .combine=rbind) %dopar%{
  dat1[i,] - c(outer(1:1200, 1:1200, function(jj,ii) diff_day_median[i, ((jj-1) %/% 100 + 1) + ((ii - 1) %/% 100)*12]))
}

tmpidx = which(is.na(dat1))
dat1[tmpidx] = dat1_imputed[tmpidx]
tmpidx = is.na(dat2)
dat2[tmpidx] = dat2_imputed[tmpidx]

saveRDS(dat1, "./output/dat1.rds")
saveRDS(dat2, "./output/dat2.rds")

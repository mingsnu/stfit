##########################################################
##### Local linear covariance estiamtion (R version) #####
##########################################################
cov.est = function(ids, time, resid, h1, h2 = h1, tt, Kern) {
  nonna.idx = !is.na(resid)
  ids = ids[nonna.idx]
  time = time[nonna.idx]
  resid = resid[nonna.idx]
  N = length(ids)
  if(length(time) != N  || length(resid) != N)
    stop("cov.est: length of ids, TT and resid should be the same.")
  id.uni = unique(ids)
  n = length(id.uni)
  ## loop through all individuals
  for (i in 1:n) {
    ind = which(ids == id.uni[i])
    m = length(ind)
    e.temp = resid[ind]
    t.temp = time[ind]
    EE.temp = e.temp %*% t(e.temp)
    T.temp = matrix(t.temp, m, m)
    
    SS.temp = as.vector(T.temp)
    TT.temp = as.vector(t(T.temp))
    EE.temp = as.vector(EE.temp)
    
    ## index for the diagnal elements
    dindex = seq(1, m ^ 2, m + 1)
    EE.temp = EE.temp[-dindex]
    SS.temp = SS.temp[-dindex]
    TT.temp = TT.temp[-dindex]
    
    if (i == 1) {
      EE = EE.temp
      SS = SS.temp
      TT = TT.temp
    } else{
      EE = c(EE, EE.temp)
      SS = c(SS, SS.temp)
      TT = c(TT, TT.temp)
    }
  }
  ## local linear covariance function estimation
  R0.hat = matrix(0, length(tt), length(tt))
  for (ii in 1:length(tt)) {
    for (jj in 1:ii) {
      X.des = cbind(rep(1, length(SS)), SS - tt[ii], TT - tt[jj])
      KW = Kern((SS - tt[ii]) / h1) * Kern((TT - tt[jj]) / h1)
      temp = X.des * KW
      a = solve(t(temp) %*% X.des) %*% t(temp) %*% EE
      R0.hat[ii, jj] = a[1, 1]
      if (jj < ii)
        R0.hat[jj, ii] = R0.hat[ii, jj]
    }
  }
  ## variance function estimation
  sigma2.hat = llreg(x = time, y = resid^2, x.eval=tt, h = h2, Kern=Kern)
  rlist = list(R0.hat = R0.hat, sigma2.hat = sigma2.hat, tt = tt)
  return(rlist)
}

###################################
##### Test using Landsat data #####
###################################
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA
mat[mat < 0] = NA
hist(mat)

#######################################
## temporal smoothing for one pixel ###
#######################################
#pdf("output/smoothing_methods_plot.pdf")
i=3
y = mat[,i]
nnaidx = !is.na(y)
ids = year[nnaidx]
x = doy[nnaidx]
y = y[nnaidx]
uni.ids = unique(ids)

## mean estimation
yfit = smooth_spline(x, y)
resid = y - yfit
.outlier <- function(y){
  whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
  which(y < whisker[1] | y > whisker[2])
}
outid = .outlier(resid)

## remove outliers
ids = ids[-outid]
x = x[-outid]
y = y[-outid]

## individual curves
plot(x, y, ylim = range(y), type = "n")
##for(i in 1:length(uni.ids)){
for(i in 1:length(uni.ids)){
  idx = ids == uni.ids[i]
  lines(x[idx], y[idx])
}
## redo mean estimation
yfit = smooth_spline(x, y)
lines(sort(x), yfit[order(x)], lwd = 4)

## residual plots
resid = y - yfit
plot(x, resid, ylim = range(resid), type = "n")
##for(i in 1:length(uni.ids)){
for(i in 1:length(uni.ids)){
  idx = ids == uni.ids[i]
  lines(x[idx], resid[idx])
}
###############################################
#### covariance estimation based on R code ####
###############################################
tt = 1:365
covf = cov.est(ids, x, resid, 100, 50, tt, epan)
persp(tt,tt,covf$R0.hat,theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
range(covf$R0.hat)
plot(covf$sigma2.hat)
##########################################
############ direct calculation ##########
##########################################
h = 100
W = epan(seq(-h, h, by = 1)/h)
tt = 1:365
R0.hat = lc_cov_1d_est(ids, x, resid, W, tt)
persp(tt,tt, R0.hat,theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
range(R0.hat)
plot(covf$sigma2.hat)
##################################################
## to avoid boundary effect, circulate the data ##
##################################################
ids1 = c(ids + 1, ids, ids - 1)
x1 = c(-(365 - x), x, x + 365)
resid1 = rep(resid, 3)
## truncate it to smaller interval
idx1 = which(x1 > -100 & x1 < 465)
ids1 = ids1[idx1]
x1 = x1[idx1]
resid1 = resid1[idx1]
tt = 1:365
R0.hat1 = lc_cov_1d_est(ids1, x1, resid1, W, tt)
persp(tt,tt, R0.hat1,theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
sigma2 = llreg(x, resid^2, h = 200, x.eval = tt)

ev = eigen(R0.hat1)
ev$values[ev$values < 0] = 0
ev.idx = which.min(cumsum(ev$values)/sum(ev$values) < 0.99)
cat("The first ", ev.idx, " eigen values are used...\n")
ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
ev.val = ev$values[1:ev.idx]
tmp = sigma2 - diag(R0.hat1)
nugg = mean(tmp)
nugg.fun = rep(nugg, length(tt))


#### only look at certen period
idx2 = 50:300
persp(tt[idx2],tt[idx2], covf$R0.hat[idx2, idx2],theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
persp(tt[idx2],tt[idx2], R0.hat[idx2, idx2],theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')
persp(tt[idx2],tt[idx2], R0.hat1[idx2, idx2],theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='s',ylab='t',zlab='R0(s,t)',ticktype='detailed')

covf$R0.hat[150:160, 150:160]
R0.hat1[150:160, 150:160]
R0.hat[150:160, 150:160]

lc_cov_1d(ids1, x1, resid1, W, 150, 150)
lc_cov_1d(ids, x, resid, W, 150, 150)
lc_cov_1d1(ids, x, resid, W,1,1)


#############################
## R version of lc_cov_1d_ ##
#############################
lc_cov_1d1 <- function(ids, time, resid, W, t1,t2){
  W_size = length(W)
  sumEEKK = 0.0
  sumKK = 0.0
  N = length(ids)
  time_min = vecmin(time);
  time_max = vecmax(time);
  
  k1_start = max(t1 - W_size%/%2+1, time_min);
  k2_start = max(t2 - W_size%/%2+1, time_min);
  k1_stop = min(t1 + W_size%/%2+1, time_max);
  k2_stop = min(t2 + W_size%/%2+1, time_max);
      
  for (i in 1:N){
    # if(i %in% c(32, 48, 63, 87, 158, 184))
    #   browser()
    if(time[i] >= k1_start & time[i] < k1_stop){
      for(j in 1:N){
        # if(j %in% c(32, 48, 63, 87, 158, 184))
        #   browser()
        if(i == j)
          next
        if(ids[i] == ids[j]){
          if(time[j] >= k2_start & time[j] < k2_stop){
            sumEEKK = sumEEKK + resid[i]*resid[j]*W[time[i] - t1 + W_size%/%2+1]*W[time[j] - t2 + W_size%/%2+1];
            sumKK = sumKK +  W[time[i] - t1 + W_size%/%2+1]*W[time[j] - t2 + W_size%/%2+1];
          }
        }
      }
    }
  }
  # browser()
  if(sumKK == 0.0){
    return(NA);
  } else{
    return(sumEEKK/sumKK);
  }
}

## test
lc_cov_1d(ids, x, resid, W, 1,1)
lc_cov_1d1(ids, x, resid, W,1,1)

###################################################
######### bandwidth selection for llreg ###########
###################################################
Proj.wi <- function(tt, time, h, sigma2=rep(1, length(time)), Kern){
  if(length(tt) != 1)
    stop("Proj.wi: tt should be a single number")
  vecX = time - tt
  vecK = Kern(vecX/h)/sigma2
  S0 = mean(vecK)
  S1 = mean(vecX*vecK) ## higher order term than S0 and S2, keep this term for accuracy 
  S2 = mean(vecX^2*vecK)
  denom = S2*S0 - S1^2
  if(denom != 0){
    return((S2-S1*vecX)*vecK/denom/length(time)) ## local linear estimator
  } else
    return(vecK/sum(vecK)) ## local constant estimator
}
llreg.gcv <- function(h.list, x, y, Kern){
  if(!ident(length(x), length(y)))
    stop("TT, y don't have the same dimension!")
  score=rep(0,length(h.list))
  for(i in 1:length(h.list)){
    h = h.list[i]
    Proj.mat = t(sapply(x, FUN = Proj.wi, time = x, h = h, Kern = Kern))
    resid = y - Proj.mat %*% y
    smoothmat = mean(diag(Proj.mat))
    score[i] = mean(resid^2)/(1-smoothmat)^2
  }
  return(data.frame(h=h.list, score=score))
}
y = resid^2
sigma2 = llreg.gcv(seq(200, 365, by = 10), x, y, epan)
plot(sigma2)
sigma2 = llreg(x,resid^2, h = 100)
mean((y-sigma2)^2)
plot(x,y)
lines(sort(x), sigma2[order(x)])
sigma2 = llreg(x,resid^2, h = 50)
mean((y-sigma2)^2)

library(locpol)
cvBwSel <- regCVBwSelC(x,resid^2, 1, EpaK, interval= c(20, 365))
thBwSel <- thumbBw(x,resid^2, 1, EpaK)
piBwSel <- pluginBw(x,resid^2, 1, EpaK)


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

registerDoParallel(16)
dat1_imputed = foreach(i = 1:ncol(dat1), .combine = cbind) %dopar%{
    if(msk[i])
        return(NA) else {
                       lmfit = lm(dat1[,i]~dat2[,i])
                       return(predict(lmfit, newdata = data.frame(cbind(1,dat2[,i]))))
                   }
}
tmpidx = which(is.na(dat1))
dat1[tmpidx] = dat1_imputed[tmpidx]
saveRDS(dat1, "./output/MYD11A1Day2010_simulated_daily_imputed.rds")

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

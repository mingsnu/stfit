##dat1 = readRDS("./data/MYD11A1Day2010_simulated.rds")
##dat2 = readRDS("./data/MOD11A1Day2010.rds")
dat3 = readRDS("./data/MYD11A1Nit2010_simulated.rds")
dat4 = readRDS("./data/MOD11A1Nit2010.rds")
if(!file.exists("./data/msk.rds")){
  msk = getMask(dat2)
  saveRDS(msk, "./data/msk.rds")
} else {
    msk = readRDS("./data/msk.rds")
}
##dat1[,msk] = NA
##dat2[,msk] = NA
dat3[,msk] = NA
dat4[,msk] = NA

DMlm=function(x, y){
    ##daily merge using linear regression, y is to be predict, x is predictor
    ##x and y should be the maxtrix of the raster stack
    Sxy=x*y
    tmp = Sxy*0
    Sx=x + tmp
    Sy=y + tmp
    Sx2=x*x + tmp
    n = 1 + tmp
    Sxy=rowSums(Sxy,na.rm=T)
    Sx=rowSums(Sx,na.rm=T)
    Sy=rowSums(Sy,na.rm=T)
    Sx2=rowSums(Sx2,na.rm=T)
    n=rowSums(n,na.rm=T)
    n[n<10]=NA ##if smaller than 10 samples, not processed
    bb=((Sxy-Sx*Sy/n)/(Sx2-Sx*Sx/n)) #slope
    aa=(Sy/n-bb*Sx/n)##intercept
    return(aa+bb*x)##predicted
}

## dat1_imputed = t(DMlm(t(dat2), t(dat1)))
## tmpidx = which(is.na(dat1))
## dat1[tmpidx] = dat1_imputed[tmpidx]
## saveRDS(dat1, "./data/MYD11A1Day2010_simulated_daily_imputed_lm.rds")

dat3_imputed = t(DMlm(t(dat4), t(dat3)))
tmpidx = which(is.na(dat3))
dat3[tmpidx] = dat3_imputed[tmpidx]
saveRDS(dat3, "./data/MYD11A1Nit2010_simulated_daily_imputed_lm.rds")

library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(stfit)
library(abind)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)


df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

#### partial missing image indexes with different missing percentage
pidx0.1 = c(66, 75, 348, 573, 605)
pidx0.4_0.6 = c(74, 156, 273, 285, 326)
pidx0.8_0.95 = c(112, 184, 318, 448, 508)
#### fully observed image indexes from different seasons
fidx1 = c(261, 265, 312, 387, 581)
fidx2 = c(145, 276, 444, 481, 587)
fidx3 = c(198, 202, 493, 549, 557)
fidx4 = c(82, 293, 505, 609, 615)
fidx = c(fidx1, fidx2, fidx3, fidx4)
fmat = mat0[fidx, ]

############################################
##### Simulation study for pidx0.4_0.6 #####
############################################
registerDoParallel(10)
## initialize error metrics matrices
RMSEmat3 = matrix(NA, length(fidx), length(pidx0.4_0.6))
NMSEmat3 = matrix(NA, length(fidx), length(pidx0.4_0.6))
AREmat3 = matrix(NA, length(fidx), length(pidx0.4_0.6))
CORmat3 = matrix(NA, length(fidx), length(pidx0.4_0.6))

for(i in 1:length(pidx0.4_0.6)){
  for(j in 1:length(fidx)){
    pmat = fmat[j,]
    ## apply missing patterns to fully observed images
    missing.idx = is.na(mat[pidx0.4_0.6[i],])
    pmat[missing.idx] = NA
    gdata = as.geodata(cbind(expand.grid(seq(1,31), seq(1, 31)), pmat))
    vario = variog(gdata)
    wls = variofit(vario, ini = c(40000, 40))
    ## wls = variofit(vario, ini = c(40000, 30))
    loci <- expand.grid(seq(0,1,l=31), seq(0,1,l=31))
    # predicting by ordinary kriging
    kc <- krige.conv(gdata, loc=expand.grid(seq(1,31), seq(1, 31)),
                     krige=krige.control(cov.pars= wls$cov.pars))
    imat = kc$predict
    image(kc, main="kriging estimates")
    image(matrix(pmat, 31))
    
    #### proposed method
    RMSEmat3[j, i] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat3[j, i] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat3[j, i] = ARE(fmat[j, missing.idx], imat[missing.idx])
    CORmat3[j, i] = cor(fmat[j, missing.idx], imat[missing.idx])
  }
}

RMSEmat1 = readRDS("./pidx0.4_0.6/RMSEmat1.rds")
RMSEmat2 = readRDS("./pidx0.4_0.6/RMSEmat2.rds")

RMSEmat1 < RMSEmat3
table(c(RMSEmat1 < RMSEmat3))
RMSEmat2 < RMSEmat3
table(c(RMSEmat2 < RMSEmat3))




library(feather)
library(dplyr)
library(geoR)
library(raster)
library(rasterVis)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])
mat0[mat0 > 2000] = NA

#######################################
##### Simulation study for spring #####
#######################################
#### partial missing image indexes with different missing percentage
pmat = readRDS("../missing_pattern/output/missing_pattern.rds")
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(88, 132, 261, 265, 615)
fidx = c(fidx1, fidx2, fidx3, fidx4)
fmat = mat0[fidx, ]
if(!dir.exists("output"))
  dir.create("output")

RMSEmat = matrix(NA, nrow(pmat), length(fidx))
NMSEmat = matrix(NA, nrow(pmat), length(fidx))
AREmat = matrix(NA, nrow(pmat), length(fidx))
CORmat = matrix(NA, nrow(pmat), length(fidx))
for(i in 1:nrow(pmat)){
  for(j in 1:length(fidx)){
    fmatj = fmat[j,]
    ## apply missing patterns to fully observed images
    missing.idx = is.na(pmat[i,])
    fmatj[missing.idx] = NA
    gdata = as.geodata(cbind(expand.grid(seq(1,31), seq(1, 31)), fmatj))
    vario = variog(gdata)
    wls = variofit(vario, ini = c(40000, 40))
    loci <- expand.grid(seq(0,1,l=31), seq(0,1,l=31))
    # predicting by ordinary kriging
    kc <- krige.conv(gdata, loc=expand.grid(seq(1,31), seq(1, 31)),
                     krige=krige.control(cov.pars= wls$cov.pars))
    imat = kc$predict
    # image(kc, main="kriging estimates")
    # image(matrix(fmatj, 31))
    #### proposed method
    RMSEmat[i, j] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat[i, j] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat[i, j] = ARE(fmat[j, missing.idx], imat[missing.idx])
    CORmat[i, j] = cor(fmat[j, missing.idx], imat[missing.idx])
  }
}

saveRDS(list(RMSEmat, NMSEmat, AREmat, CORmat), "./output/res.rds")


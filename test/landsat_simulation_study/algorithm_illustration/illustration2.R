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

#######################################
##### Simulation study for summer #####
#######################################
#### partial missing image indexes with different missing percentage
##### selected partially observed images indexes
pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)[c(6,8,14,15)]
pmat = readRDS("../missing_pattern/output/missing_pattern.rds")[c(6,8,14,15),]
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(132, 261, 265, 615, 657) 
fidx = c(fidx1, fidx2, fidx3, fidx4)[c(3,7,14,18)]
fmat = mat0[fidx, ]
if(!dir.exists("output2"))
  dir.create("output2")

###### define function used for mean estimation ######
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)

## collapse the partial and full index for parallel
## matrix of MxN, column stacking
N = nrow(fmat)
M = nrow(pmat)
registerDoParallel(10)
r.list = list()
for(i in 1:4){
  mat = mat0
  ## apply missing patterns to fully observed images
  missing.idx = is.na(pmat[i,])
  mat[fidx[i], missing.idx] = NA
  
  if(file.exists(paste0("./output2/output2_P", pidx[i], "_F", fidx[i], ".rds"))){
    res1 <- readRDS(paste0("./output2/output2_P", pidx[i], "_F", fidx[i], ".rds"))
  } else {
    res1 <- gapfill_landsat(year, doy, mat, 31, 31, nnr=30,
                            use.intermediate.result = FALSE, intermediate.save = FALSE)
    saveRDS(res1, paste0("./output2/output2_P", pidx[i], "_F", fidx[i], ".rds"))
  }
  imat = res1$imat[fidx[i],]
  r.list[[(i-1)*3+1]] = raster(matrix(mat0[fidx[i],], 31))
  r.list[[(i-1)*3+2]] = raster(matrix(mat[fidx[i],], 31))
  r.list[[(i-1)*3+3]] = raster(matrix(imat, 31))
}

s = stack(r.list)
pdf("output2/impu_res_example.pdf", width =6, height=5.5)
levelplot(s,par.settings = colthm, index.cond=list(c(1,4,7,10,2,5,8,11,3,6,9,12)),
          names.attr = c(rbind(paste0("F", c(3,7,14,18)), paste0("F", c(3,7,14,18), "P", c(6,8,14,15)),
                         paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " impu"))))

dev.off()





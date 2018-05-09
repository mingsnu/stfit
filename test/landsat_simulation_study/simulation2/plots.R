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


###### variables used for Gapfill package ######
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))


#### partial missing images with different missing percentage
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

########## pidx0.1
tmpdat = mat0[pidx0.1,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.1.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()

for(i in 1:length(pidx0.1)){
    ## artificial partial missing images applying pidx0.1[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.1[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using stfit and gapfill
    imat0.1stfit = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.1/res1_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
        imat0.1stfit[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.1/res2_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1stfit, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.1_", pidx0.1[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(length(fidx),4)))
    dev.off()
}


########## pidx0.4_0.6
tmpdat = mat0[pidx0.4_0.6,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.4_0.6.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()

for(i in 1:length(pidx0.4_0.6)){
    ## artificial partial missing images applying pidx0.4_0.6[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.4_0.6[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using stfit and gapfill
    imat0.1stfit = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.4_0.6/res1_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
        imat0.1stfit[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.4_0.6/res2_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1stfit, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.4_0.6_", pidx0.4_0.6[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(length(fidx),4)))
    dev.off()
}


########## pidx0.8_0.95
tmpdat = mat0[pidx0.8_0.95,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.8_0.95.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()

for(i in 1:length(pidx0.8_0.95)){
    ## artificial partial missing images applying pidx0.8_0.95[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.8_0.95[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using stfit and gapfill
    imat0.1stfit = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.8_0.95/res1_pidx_", pidx0.8_0.95[i], "_fidx_", fidx[j], ".rds"))
        imat0.1stfit[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.8_0.95/res2_pidx_", pidx0.8_0.95[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1stfit, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.8_0.95_", pidx0.8_0.95[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(length(fidx),4)))
    dev.off()
}



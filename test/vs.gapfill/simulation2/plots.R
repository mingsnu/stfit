library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(Gapfill)
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
res = readRDS("../our_output/our_res.rds")
missingpct = apply(mat0[res$idx$idx.partialmissing,], 1, function(x) sum(is.na(x))/length(x))
pidx0.1 = res$idx$idx.partialmissing[which(missingpct < 0.1)]
pidx0.1 = pidx0.1[c(1,2,7,8,9)]
pidx0.4_0.6 = res$idx$idx.partialmissing[which(missingpct <= 0.6 & missingpct > 0.4)]
pidx0.4_0.6 = pidx0.4_0.6[1:5]
pidx0.8_0.95 = res$idx$idx.partialmissing[which(missingpct <= 0.95 & missingpct > 0.8)]

########## pidx0.1
tmpdat = mat0[pidx0.1,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.1.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()

fidx = sort(c(481, 587, 107, 82, 460, 599))
## original fully observed images
fmat = mat0[fidx, ]

for(i in 1:length(pidx0.1)){
    ## artificial partial missing images applying pidx0.1[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.1[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using gastif and gapfill
    imat0.1gas = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.1/res1_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gas[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.1/res2_pidx_", pidx0.1[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1gas, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.1_", pidx0.1[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(6,4)))
    dev.off()
}


########## pidx0.4_0.6
tmpdat = mat0[pidx0.4_0.6,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.4_0.6.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()
fidx = sort(c(481, 587, 107, 82, 460, 599))
## original fully observed images
fmat = mat0[fidx, ]

for(i in 1:length(pidx0.4_0.6)){
    ## artificial partial missing images applying pidx0.4_0.6[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.4_0.6[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using gastif and gapfill
    imat0.1gas = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.4_0.6/res1_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gas[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.4_0.6/res2_pidx_", pidx0.4_0.6[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1gas, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.4_0.6_", pidx0.4_0.6[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(6,4)))
    dev.off()
}


########## pidx0.8_0.95
tmpdat = mat0[pidx0.8_0.95,]
tmpdat[!is.na(tmpdat)] = 1
pdf("./plot/missing_pattern_pidx0.8_0.95.pdf")
tmpstack = mat2stack(tmpdat, 31)
levelplot(tmpstack, colorkey=FALSE, names.attr = c(paste0("P", 1:5)))
dev.off()
fidx = sort(c(481, 587, 107, 82, 460, 599))
## original fully observed images
fmat = mat0[fidx, ]

for(i in 1:length(pidx0.8_0.95)){
    ## artificial partial missing images applying pidx0.8_0.95[i]
    pmat0.1 = fmat
    missing.idx = is.na(mat0[pidx0.8_0.95[i],])
    pmat0.1[, missing.idx] = NA
    ## initialize imputed partial missing images using gastif and gapfill
    imat0.1gas = imat0.1gap = pmat0.1
    for(j in 1:length(fidx)){
        res1 = readRDS(paste0("./pidx0.8_0.95/res1_pidx_", pidx0.8_0.95[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gas[j,] = res1$imat[fidx[j],]
        yidx = which(year[fidx[j]] == yearuni)
        didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
        didxinterval = max(1,didx-6):min(46, didx + 6)
        yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
        res2 = readRDS(paste0("./pidx0.8_0.95/res2_pidx_", pidx0.8_0.95[i], "_fidx_", fidx[j], ".rds"))
        imat0.1gap[j,] = c(res2$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
    }
    s = mat2stack(rbind(fmat, pmat0.1, imat0.1gas, imat0.1gap), 31)
    pdf(paste0("./plot/pidx0.8_0.95_", pidx0.8_0.95[i], ".pdf"))
    print(rasterVis::levelplot(s, par.settings = colthm, layout=c(6,4)))
    dev.off()
}



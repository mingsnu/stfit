library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(stfit)
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

df = read_feather("../../data/features_106_wide.feather")
## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

res = readRDS("../core/our_output/our_res.rds")
##########################################
#### partial missing image selection #####
##########################################
## 0~30; 30~70; 70~100
missingpct = apply(mat[res$idx$idx.partialmissing,], 1, function(x) sum(is.na(x))/length(x))
pidx0.3 = res$idx$idx.partialmissing[which(missingpct < 0.3)]
pidx0.3 = pidx0.3[c(2,12,25,31,32)]
pidx0.3.missingpct = apply(mat[pidx0.3,], 1, function(x) sum(is.na(x))/length(x))
pidx0.3_0.7 = res$idx$idx.partialmissing[which(missingpct <= 0.7 & missingpct > 0.3)]
pidx0.3_0.7 = pidx0.3_0.7[c(1,3,19, 21, 53)]
pidx0.3_0.7.missingpct = apply(mat[pidx0.3_0.7,], 1, function(x) sum(is.na(x))/length(x))
pidx0.7_0.99 = res$idx$idx.partialmissing[which(missingpct <= 0.99 & missingpct > 0.7)]
pidx0.7_0.99 = pidx0.7_0.99[c(2, 5, 9, 10, 12)]
pidx0.7_0.99.missingpct = apply(mat[pidx0.7_0.99,], 1, function(x) sum(is.na(x))/length(x))
# s = mat2stack(mat[pidx0.3_0.7,], nrow = 31)
# levelplot(s)
# data.frame(year = year[c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99)],
#            DOY = doy[c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99)],
#            idx = c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99),
#            pct = apply(mat[c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99),], 1, function(x) sum(is.na(x))/length(x)))
##### missing pattern image information =================
df = data.frame(pn = paste0("P", 1:15),
                pct=apply(mat[c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99),], 1,
                          function(x) sum(is.na(x))/length(x)))
print(xtable::xtable(df), include.rownames=FALSE)

#### save missing pattern data to file  =================
saveRDS(mat[c(pidx0.3, pidx0.3_0.7, pidx0.7_0.99),], "output/missing_pattern.rds")


###########################################
###########################################
##### Plot ============================================
pmat = readRDS("output/missing_pattern.rds")
pmat[!is.na(pmat)] = 1
pdf("./plot/missing_pattern.pdf", width=6.9, height=5.5)
s = mat2stack(pmat, 31)
myTheme <- rasterTheme(region=gray.colors(1)) 
levelplot(s, colorkey=FALSE, names.attr = c(paste0("P", 1:15)), par.settings = myTheme)
dev.off()




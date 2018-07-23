library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)

dat0 = readRDS("./data/MYD11A1Nit2010.rds")
# colthm = RdBuTheme()
# colthm$regions$col = rev(colthm$regions$col)
msk = readRDS("./data/msk.rds")
dat0[,msk] = NA
dat0 = dat0/100 ## the unit after transform is 'K'
## stackRaster
r.list = list()
dd = seq(1,365, 30)
for(i in 1:length(dd)){
  r.list[[i]] = raster(matrix(dat0[dd[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
png("./output/MYD11A1Nit2010.png", width=700, height=550)
print(levelplot(s[[1:12]], names.attr=as.character(dd[1:12]), main = "",
                layout=c(4,3)))
dev.off()

## percentage of missing values
N = sum(!msk)
pctmissing = apply(dat0[,!msk], 1, function(x) sum(is.na(x))/N)
summary(pctmissing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2769  0.5851  0.7122  0.7009  0.8360  0.9940 
## the overall percentage of missing
sum(is.na(dat0[,!msk]))/(N*nrow(dat0))
## 0.7008861
saveRDS(pctmissing, "./output/pctmissing_nit.rds")
pdf("./output/pctmissing_nit.pdf", width=6.9, height=5.5)
par(mar=c(4.2,4.2,0.2,0.2))
hist(pctmissing, breaks = 20, xlab = "Percentage of missing data over land areas",  main="")
dev.off()

## percentage of missing values after daily merging
dat1 = readRDS("./data/MYD11A1Nit2010_daily_imputed_lm.rds")
pctmissing = apply(dat1[,!msk], 1, function(x) sum(is.na(x))/N)
summary(pctmissing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1549  0.4636  0.6025  0.5961  0.7265  0.9824 
sum(is.na(dat1[,!msk]))/(N*nrow(dat1))
# 0.5960919


######## day imputation plots
### level 3 plot
dat = readRDS("./data/MYD11A1Nit2010_daily_imputed_lm.rds")
msk = readRDS("./data/msk.rds")
dat[, msk] = NA
idx1 = c(t(outer(seq(25, 1200, by = 50), seq(25, 1200, by = 50),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
res1 = readRDS("./output_nit/res1.rds")
l =  res1$idx$idx.partialmissing[1]
r1 = raster(matrix(dat[,idx1][l,], 24))
r2 = raster(matrix(res1$imat[l,], 24))
s1 = stack(r1, r2)
print(levelplot(s1, names.attr=c("Level 3 image", "Imputed level 3 image")))

### level 2 plot
dat = readRDS("./output_day/lvl1_impu_outlier_removed.rds")
## visualize the level two imputation results.
bIdx = c(t(outer(seq(5, 1200, by = 10), seq(5, 1200, by = 10),
                 FUN = function(ridx, cidx){
                   (ridx-1) * 1200 + cidx
                 })))
r3 = raster(matrix(dat[,bIdx][l,], 120))
dat1 = readRDS("output_day/lvl2_impu.rds")
r4 = raster(matrix(dat1[,bIdx][l,], 120))
s2 = stack(r3, r4)
print(levelplot(s2, names.attr=c("Level 2 image", "Imputed level 2 image")))

### level 1 plot
dat = readRDS("./output_day/lvl2_impu.rds")
r5 = raster(matrix(dat[l,], 1200))
dat1 = readRDS("output_day/dat_imputed.rds")
r6 = raster(matrix(dat1[l,], 1200))
s3 = stack(r5, r6)
print(levelplot(s3, names.attr=c("Original image", "Imputed image")))

### individual plots
r.list = list(r1, r2, r3, r4, r5, r6)
ll = seq(279,317, length.out = 16)
for(i in 1:6){
  pdf(paste0("./output_day/lst_day", l, "_", i, ".pdf"), width=5.7, height=5.2)
  print(levelplot(r.list[[i]]/100, margin = FALSE, at = ll))
  dev.off()
}










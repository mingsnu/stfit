library(raster)
tstset=c(14,27,47,63,76,95,107,129,152,169,183, 202,218,225,255,273,289,308,316,353)
mskset=c(10,31,44,52,74,94,112,138,154,170,185, 193,216,229,260,275,284,307,325,357)
dat = readRDS("./data/MYD11A1Nit2010.rds")
idat = readRDS("./data/MYD11A1Nit2010_simulated_daily_imputed_lm.rds")
idat1 = readRDS("./output_nnr50_nit_sim_clip23500/dat_imputed.rds")


#### calculate the accuracy on testset
tdat = dat[tstset, ]
idat = idat[tstset, ]
idat1 = idat1[tstset,]
idx = is.na(dat[mskset, ])

## the overall accuracy
d = c(tdat[idx])-c(idat1[idx])
sqrt(sum(d^2, na.rm=TRUE)/sum(!is.na(d)))
cor(c(tdat[idx]), c(idat1[idx]), use = "complete.obs")

## accuracy for daily merging
d = c(tdat[idx])-c(idat[idx])
sqrt(sum(d^2, na.rm=TRUE)/sum(!is.na(d)))
cor(c(tdat[idx]), c(idat[idx]), use = "complete.obs")

## accuracy for our algorithm (excluding daily merging part)
idx1 = is.na(idat[idx])
d = c(tdat[idx][idx1]) - c(idat1[idx][idx1])
sqrt(sum(d^2, na.rm=TRUE)/sum(!is.na(d)))
cor(c(tdat[idx][idx1]), c(idat1[idx][idx1]), use = "complete.obs")

## percentage daily merging imputed
pct1 = sum(!idx1)/sum(idx)
pct1
## percentage gapfilling merging imputed
pct2 = sum(!is.na(idat1[idx]))/sum(idx)
pct2

## Overall accuracy
## [1] 216.153
## [1] 0.9814926
## daily lm impu
## [1] 176.2332
## [1] 0.9880501
## our method only
## [1] 304.0221
## [1] 0.9658377

## [1] 216.153
## [1] 0.9814926
## [1] 176.2332
## [1] 0.9880501
## [1] 304.0221
## [1] 0.9658377
## [1] 0.4165464
## [1] 0.9103658


## [1] 215.9183
## [1] 0.9815328
## [1] 176.2332
## [1] 0.9880501
## [1] 303.3679
## [1] 0.9659797
## [1] 0.4165464
## [1] 0.9103658

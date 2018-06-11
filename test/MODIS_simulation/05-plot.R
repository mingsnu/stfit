library(raster)
library(stfit)
tstset=c( 8,47,56,76,85, 99,117,139,147,170,184,197,221,239,256,282,295,313,327,346)
mskset=c(12,42,58,67,94,110,126,129,154,174,177,199,222,232,252,285,304,306,332,342)
dat = readRDS("./data/MYD11A1Day2010.rds")
## idat = readRDS("./data/MYD11A1Day2010_simulated_daily_imputed_lm.rds")
idat1 = readRDS("./output_nnr50_sim/dat_imputed.rds")

#### calculate the accuracy on testset
tdat = c(dat[tstset, ])/100
## idat = c(idat[tstset, ])/100
idat1 = c(idat1[tstset,])/100
idx = c(is.na(dat[mskset, ]))


pdf("./plots/T2_scatter.pdf")
par(mar = c(4.2, 4.2, 0.2, 0.2))
smoothScatter(tdat[idx], idat1[idx], xlab = "Observed T2 (K)", ylab = "Predicted T2 (K)",
               colramp = colorRampPalette(c("white", "blue", "green", "yellow", "orange", "red")),
              xlim = c(230, 340), ylim = c(230, 340))
rmspe = round(stfit::RMSE(tdat[idx], idat1[idx]), 2)
r2 = round(cor(tdat[idx], idat1[idx], use = "complete.obs"), 2)
mape = round(stfit::ARE(tdat[idx], idat1[idx]), 2)
text(228, 340, paste0("RMSPE: ", rmspe), adj = c(0,0.5))
text(228, 335, paste0("MAPE: ", mape), adj = c(0,0.5))
text(228, 330, bquote(R^2: .(r2)), adj = c(0,0.5))
dev.off()

## pdf("./plots/day_lm.pdf")
## par(mar = c(4.2, 4.2, 0.2, 0.2))
## smoothScatter(tdat[idx], idat[idx])
## dev.off()

## pdf("./plots/day_hmri.pdf")
## par(mar = c(4.2, 4.2, 0.2, 0.2))
## idx1 = is.na(idat[idx])
## smoothScatter(tdat[idx][idx1], idat1[idx][idx1])
## dev.off()

tstset=c(14,27,47,63,76,95,107,129,152,169,183, 202,218,225,255,273,289,308,316,353)
mskset=c(10,31,44,52,74,94,112,138,154,170,185, 193,216,229,260,275,284,307,325,357)
dat = readRDS("./data/MYD11A1Nit2010.rds")
## idat = readRDS("./data/MYD11A1Nit2010_simulated_daily_imputed_lm.rds")
idat1 = readRDS("./output_nnr50_nit_sim_clip23500/dat_imputed.rds")

tdat = c(dat[tstset, ])/100
## idat = c(idat[tstset, ])/100
idat1 = c(idat1[tstset,])/100
idx = c(is.na(dat[mskset, ]))

pdf("./plots/T4_scatter.pdf")
par(mar = c(4.2, 4.2, 0.2, 0.2))
smoothScatter(tdat[idx], idat1[idx], xlab = "Observed T4 (K)", ylab = "Predicted T4 (K)",
              colramp = colorRampPalette(c("white", "blue", "green", "yellow", "orange", "red")),
              xlim = c(230, 340), ylim = c(230, 340))
rmspe = round(stfit::RMSE(tdat[idx], idat1[idx]), 2)
r2 = round(cor(tdat[idx], idat1[idx], use = "complete.obs"), 2)
mape = round(stfit::ARE(tdat[idx], idat1[idx]), 2)
text(228, 340, paste0("RMSPE: ", rmspe), adj = c(0,0.5))
text(228, 335, paste0("MAPE: ", mape), adj = c(0,0.5))
text(228, 330, bquote(R^2: .(r2)), adj = c(0,0.5))
dev.off()


library(raster)
library(stfit)
smoothScatter1 <- function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", 
                                                                                             blues9)), nrpoints = 100, ret.selection = FALSE, pch = ".", 
                            cex = 1, col = "black", transformation = function(x) x^0.25, 
                            postPlotHook = box, xlab = NULL, ylab = NULL, xlim, ylim, 
                            xaxs = par("xaxs"), yaxs = par("yaxs"), ...) 
{
  if (!is.numeric(nrpoints) || nrpoints < 0 || length(nrpoints) != 
      1) 
    stop("'nrpoints' should be numeric scalar with value >= 0.")
  nrpoints <- round(nrpoints)
  ret.selection <- ret.selection && nrpoints > 0
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab
  x <- cbind(xy$x, xy$y)[I <- is.finite(xy$x) & is.finite(xy$y), 
                         , drop = FALSE]
  if (ret.selection) 
    iS <- which(I)
  if (!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
    x <- x[I <- min(xlim) <= x[, 1] & x[, 1] <= max(xlim), 
           , drop = FALSE]
    if (ret.selection) 
      iS <- iS[I]
  }
  else {
    xlim <- range(x[, 1])
  }
  if (!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
    x <- x[I <- min(ylim) <= x[, 2] & x[, 2] <= max(ylim), 
           , drop = FALSE]
    if (ret.selection) 
      iS <- iS[I]
  }
  else {
    ylim <- range(x[, 2])
  }
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth,list(xlim, ylim)) 
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens[] <- transformation(dens)
  image(xm, ym, z = dens, col = colramp(256), xlab = xlab, 
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs, 
        ...)
  if (!is.null(postPlotHook)) 
    postPlotHook()
  if (nrpoints > 0) {
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) == 
                ncol(dens))
    ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] - 
                                                          xm[1]))
    iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] - 
                                                          ym[1]))
    sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
    x <- x[sel, , drop = FALSE]
    points(x, pch = pch, cex = cex, col = col)
    if (ret.selection) 
      iS[sel]
  }
}


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
smoothScatter1(tdat[idx], idat1[idx], xlab = "Observed T2 (K)", ylab = "Predicted T2 (K)",
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
smoothScatter1(tdat[idx], idat1[idx], xlab = "Observed T4 (K)", ylab = "Predicted T4 (K)",
              colramp = colorRampPalette(c("white", "blue", "green", "yellow", "orange", "red")),
              xlim = c(230, 340), ylim = c(230, 340))
rmspe = round(stfit::RMSE(tdat[idx], idat1[idx]), 2)
r2 = round(cor(tdat[idx], idat1[idx], use = "complete.obs"), 2)
mape = round(stfit::ARE(tdat[idx], idat1[idx]), 2)
text(228, 340, paste0("RMSPE: ", rmspe), adj = c(0,0.5))
text(228, 335, paste0("MAPE: ", mape), adj = c(0,0.5))
text(228, 330, bquote(R^2: .(r2)), adj = c(0,0.5))
dev.off()


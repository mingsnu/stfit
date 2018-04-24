##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA
mat[mat < 0] = NA
hist(mat)

#######################################
## temporal smoothing for one pixel ###
#######################################
pdf("output/smoothing_methods_plot.pdf")
x = doy
x.eval = 1:365
for(i in seq(1, 961, 10)){
  y = mat[, i]
  ypred1 = smooth_spline(x, y, x.eval)
  ypred2 = llreg(x, y, x.eval, 60, epan)
  ypred3 = lpreg(x, y, x.eval,0.3)
  ypred4 = spreg(x, y, x.eval, nbasis=11)
  
  plot(x, y, pch = 20, col=rgb(0.2,0.2,0.2, alpha=0.5), ylim=c(0, 2000))
  lines(x.eval, ypred1, col=2, lwd = 2)
  lines(x.eval, ypred2, col=3, lwd = 2)
  lines(x.eval, ypred3, col=4, lwd = 2)
  lines(x.eval, ypred4, col=5, lwd = 2)
  legend("topright", legend = c("ss", "ll.60", "lp.0.3", "sp.11"), col = 2:5,lty=1, lwd =2)
}
dev.off()

##########################################
### temporal smoothing for one cluster ###
##########################################
cluster=readRDS("cluster.rds")
x.eval = 1:365
pdf("output/smoothing_methods_cluster_plot.pdf")
for(i in unique(sort(cluster))){
  clidx = cluster == i
  x = rep(doy, sum(clidx))
  y = c(mat[,clidx])
  idx = which(!is.na(y))
  y = y[idx]
  x = x[idx]
  
  ypred1 = smooth_spline(x, y, x.eval)
  ypred2 = llreg(x, y, x.eval, 60, epan)
  ypred3 = lpreg(x, y, x.eval,0.3)
  ypred4 = spreg(x, y, x.eval, nbasis=11)
  
  plot(x, y, pch = 20, col=rgb(0.2,0.2,0.2, alpha=0.05), ylim=c(0, 2000))
  lines(x.eval, ypred1, col=2, lwd = 2)
  lines(x.eval, ypred2, col=3, lwd = 2)
  lines(x.eval, ypred3, col=4, lwd = 2)
  lines(x.eval, ypred4, col=5, lwd = 2)
  legend("topright", legend = c("ss", "ll.60", "lp.0.3", "sp.11"), col = 2:5,lty=1, lwd =2)
}
dev.off()

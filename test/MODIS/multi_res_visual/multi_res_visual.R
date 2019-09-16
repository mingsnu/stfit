# ## level 1 
# par(mar=rep(0,4))
# x = seq(0,1,length.out = 26)[-c(1,26)]
# y = expand.grid(x,x)
# plot(y$Var1, y$Var2, xlim=c(0,1), ylim=c(0,1), axes=F, ylab="", xlab="", pch=15, col = "lightgray", cex=0.5)
# rect(0,0,1,1)
# ## level 2
# xx = seq(0,1, length.out = 5)[-c(1,5)]
# for(i in 1:length(xx)){
#   lines(c(0,1), c(xx[i], xx[i]))
#   lines(c(xx[i], xx[i]), c(0,1))
# }
# x = seq(0,0.25,length.out = 32)[-c(1,32)]
# y = expand.grid(x,x)
# points(y$Var1, y$Var2, pch=".", col = 1)
# rect(0,0,1,1)
# 
# ## level 3
# xx = seq(0,0.25, length.out = 12)[-c(1,12)]
# for(i in 1:length(xx)){
#   lines(c(0,0.25), c(xx[i], xx[i]), col="blue")
#   lines(c(xx[i], xx[i]), c(0,0.25), col = "blue")
# }
# 

# # ## 300x300
# myTheme <- rasterTheme(region=c('white', gray.colors(1)))
# library(rasterVis)
# mat = matrix(0, 1200, 1200)
# mat[seq(5, 1200, by = 10), seq(5,1200, by = 10)] = 0.6
# r = raster(mat)
# #pdf("multires_illustration1.pdf", width =7, height = 7)
# levelplot(r, margin = FALSE, par.settings= myTheme, colorkey=FALSE)
# #dev.off()
# 
# mat = matrix(0, 120, 120)
# mat[seq(5,120, by=10), seq(5,120, by=10)] = 0.6
# r = raster(mat)
# #pdf("multires_illustration2.pdf", width =7, height = 7)
# levelplot(r, margin = FALSE, par.settings= myTheme, at=seq(0,1,0.1), colorkey=FALSE)
# #dev.off()
# 
# mat = matrix(0.6, 30, 30)
# r = raster(mat)
# #pdf("multires_illustration2.pdf", width =7, height = 7)
# levelplot(r, margin = FALSE, par.settings= myTheme, at=seq(0,1,0.1), colorkey=FALSE)

####### design plots
par(mar=rep(0,4), bg=NA)
x = seq(5,120,by=10)
y = expand.grid(x,x)
##plot(y$Var1, y$Var2, xlim=c(0,120), ylim=c(0,120), axes=F, ylab="", xlab="", pch=15, col = "white", cex=0.5, type="n")
plot(y$Var1, y$Var2, xlim=c(0,120), ylim=c(0,120), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5)
rect(0,0,120,120)

xx = seq(30,90, by=30)
for(i in 1:length(xx)){
  lines(c(0, 120), c(xx[i], xx[i]), lty=2, col="gray")
  lines(c(xx[i], xx[i]), c(0,120), lty = 2, col="gray")
}


####### illustration of design plots for large images
## original image
# par(mar=rep(0,4), bg=NA)
x = seq(1,79,by=3)
y = expand.grid(x,x)
plot(y$Var1, y$Var2, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y$Var1, y$Var2, pch=19, col = gray(0.4), cex=0.8)
rect(0,0,80,80, lwd=2)
xx = seq(26.5,80, by=27)
xx = seq(8.5,80, by=9)
for(i in 1:length(xx)){
  lines(c(0, 80), c(xx[i], xx[i]), lty=2, col="red", lwd=3)
  lines(c(xx[i], xx[i]), c(0,80), lty = 2, col="red", lwd=3)
}

## lvl 2 sampling points
x = seq(1,79,by=3)
y = expand.grid(x,x)
x1 = x[seq(2,length(x), by =3)]
y1 = expand.grid(x1,x1)
idx = paste(y[,1], y[,2], sep=".")
idx1 = paste(y1[,1], y1[,2], sep=".")
y = y[!idx %in% idx1,]
plot(y$Var1, y$Var2, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y$Var1, y$Var2, pch=19, col = gray(0.8), cex=0.8)
points(y1$Var1, y1$Var2, pch=15, col = gray(0.4), cex=1)
rect(0,0,80,80, lwd=2)


## lvl 2 image
x = seq(1,79,by=3)
y = expand.grid(x,x)
x1 = x[seq(2,length(x), by =3)]
y1 = expand.grid(x1,x1)
idx = paste0(y[,1], y[,2])
idx1 = paste0(y1[,1], y1[,2])
y = y[!idx %in% idx1,]
plot(y$Var1, y$Var2, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y1$Var1, y1$Var2, pch=15, col = gray(0.4), cex=1)
rect(0,0,80,80, lwd=2)
xx = seq(26.5,80, by=27)
for(i in 1:length(xx)){
  lines(c(0, 80), c(xx[i], xx[i]), lty= 2, col="blue", lwd=3)
  lines(c(xx[i], xx[i]), c(0,80), lty = 2, col="blue", lwd=3)
}

## lvl 3 sampling points
x = seq(1,79,by=3)
y = expand.grid(x,x)
x1 = x[seq(2,length(x), by =3)]
y1 = expand.grid(x1,x1)
x2 = x1[seq(2, length(x1), by = 3)]
y2 = expand.grid(x2, x2)
idx = paste0(y[,1], y[,2])
idx1 = paste0(y1[,1], y1[,2])
idx2 = paste0(y2[,1], y2[,2])
y = y[!idx %in% idx1,]
y1 = y1[!idx1 %in% idx2,]
plot(y$Var1, y$Var2, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y1$Var1, y1$Var2, pch=15, col = gray(0.8), cex=1)
points(y2$Var1, y2$Var2, pch=17, col = gray(0.4), cex=1.2)
rect(0,0,80,80, lwd=2)

######## illustration of multi resolution imputation
## lvl 3 imputation
x2 = x1[seq(2, length(x1), by = 3)]
y2 = expand.grid(x2, x2)
plot(0, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y2$Var1, y2$Var2, pch=17, col = gray(0.4), cex=1.2)
rect(0,0,80,80, lwd=2)

## lvl 2 imputation
x1 = x[seq(2,length(x), by =3)]
y1 = expand.grid(x1,x1)
x2 = x1[seq(2, length(x1), by = 3)]
y2 = expand.grid(x2, x2)
idx1 = paste0(y1[,1], y1[,2])
idx2 = paste0(y2[,1], y2[,2])
y1 = y1[!idx1 %in% idx2,]
plot(0, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y1$Var1, y1$Var2, pch=15, col = gray(0.5), cex=1)
points(y2$Var1, y2$Var2, pch=17, col = gray(0.4), cex=1.2)
rect(0,0,80,80, lwd=2)
xx = seq(26.5,80, by=27)
for(i in 1:length(xx)){
  lines(c(0, 80), c(xx[i], xx[i]), lty= 2, col="blue", lwd=3)
  lines(c(xx[i], xx[i]), c(0,80), lty = 2, col="blue", lwd=3)
}

## lvl 1 imputation
# par(mar=rep(0,4), bg=NA)
x = seq(1,79,by=3)
y = expand.grid(x,x)
x1 = x[seq(2,length(x), by =3)]
y1 = expand.grid(x1,x1)
x2 = x1[seq(2, length(x1), by = 3)]
y2 = expand.grid(x2, x2)
idx = paste0(y[,1], y[,2])
idx1 = paste0(y1[,1], y1[,2])
idx2 = paste0(y2[,1], y2[,2])
y = y[!idx %in% idx1,]
y1 = y1[!idx1 %in% idx2,]

plot(0, xlim=c(0,79), ylim=c(0,79), axes=F, ylab="", xlab="", pch=15, col = "gray", cex=0.5, type = "n")
points(y$Var1, y$Var2, pch=19, col = gray(0.5), cex=0.8)
rect(0,0,80,80, lwd=2)
points(y1$Var1, y1$Var2, pch=15, col = gray(0.5), cex=1)
points(y2$Var1, y2$Var2, pch=17, col = gray(0.5), cex=1.2)
xx = seq(26.5,80, by=27)
xx = seq(8.5,80, by=9)
for(i in 1:length(xx)){
  lines(c(0, 80), c(xx[i], xx[i]), lty=2, col="red", lwd=2)
  lines(c(xx[i], xx[i]), c(0,80), lty = 2, col="red", lwd=2)
}


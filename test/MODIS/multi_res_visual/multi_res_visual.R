## level 1 
par(mar=rep(0,4))
x = seq(0,1,length.out = 26)[-c(1,26)]
y = expand.grid(x,x)
plot(y$Var1, y$Var2, xlim=c(0,1), ylim=c(0,1), axes=F, ylab="", xlab="", pch=20, col = "lightgray")
rect(0,0,1,1)
## level 2
xx = seq(0,1, length.out = 5)[-c(1,5)]
for(i in 1:length(xx)){
  lines(c(0,1), c(xx[i], xx[i]))
  lines(c(xx[i], xx[i]), c(0,1))
}
x = seq(0,0.25,length.out = 32)[-c(1,32)]
y = expand.grid(x,x)
points(y$Var1, y$Var2, pch=".", col = 1)
rect(0,0,1,1)

## level 3
xx = seq(0,0.25, length.out = 32)[-c(1,32)]
for(i in 1:length(xx)){
  lines(c(0,0.25), c(xx[i], xx[i]), col="blue")
  lines(c(xx[i], xx[i]), c(0,0.25), col = "blue")
}
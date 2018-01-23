###########
library(raster)
library(rasterVis)
files = list.files("../data/2014", "*.tif$", full.names = TRUE)
## smoothing spline vs fourier basis spline regression
r = brick(stack(files))
levelplot(r, par.settings = RdBuTheme)
mat = t(values(r)) ## note that the image value is stacked row-wise
missIdx = apply(mat, 1, function(x) all(is.na(x)))
mat = mat[!missIdx, ]

covest = sparse_cov_est(mat, 31,31,1)
scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, dims = c(961, 961), symmetric = TRUE)
image(scovest)
image(scovest[1:50, 1:50])
scovest.eigen = eigen(scovest)
str(scovest.eigen)
plot(scovest.eigen$values)
image(matrix(scovest.eigen$vectors[,1], 31))
image(matrix(scovest.eigen$vectors[,2], 31))
image(matrix(scovest.eigen$vectors[,4], 31))



#### Validating sparse_cov_est function
mat = matrix(rnorm(100), 25)
mat = t(mat - apply(mat,1, mean))
tmp = sparse_cov_est(mat, 5, 5, 2)
stmp = sparseMatrix(tmp$ridx, tmp$cidx, x = tmp$value, dims = c(25,25), symmetric = TRUE)
stmp
image(stmp)
stmp[1:5,1:5]
aa = cov(mat)
aa[1:5,1:5]

##
df1 = readRDS("df1.rds")
mean_cov_est = meanCovEst(df1, 0, matrix(0, 961,961))
mean_cov_est$R.matrix[1:10,1:10]
mce = readRDS("test/mean_cov_est.rds")
mce$R.matrix[1:10,1:10]
mat = matrix(df1$yresid, 961)
mean = apply(mat, 1, mean, na.rm = TRUE)
head(mean)
head(mce$mean.curve)

mat = mat - mean
covest = sparse_cov_est(mat, 31, 31, 5)
scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, dims = c(961, 961), symmetric = TRUE)
image(scovest)
scovest[1:10, 1:10]

scovest.eigen = eigen(scovest)
image(matrix(scovest.eigen$vectors[,1], 31))
image(matrix(scovest.eigen$vectors[,2], 31))
image(matrix(scovest.eigen$vectors[,4], 31))
sev = scovest.eigen$values
sev[sev<0] = 0
which.min(cumsum(sev)/sum(sev) < 0.99)
plot(sev)

dcovest.eigen = eigen(mce$R.matrix)
image(matrix(dcovest.eigen$vectors[,1], 31))
image(matrix(dcovest.eigen$vectors[,2], 31))
ev = dcovest.eigen$values
ev[ev<0] = 0
which.min(cumsum(ev)/sum(ev) < 0.99)
plot(ev)


10 x 10 sparse Matrix of class "dsCMatrix"

[1,] 15287.192 11528.618  9949.501  8640.314  8605.827  8155.825     .         .        .        .   
[2,] 11528.618 21418.839  6996.810 10050.090 10342.738  9789.687  9428.786     .        .        .   
[3,]  9949.501  6996.810 25543.845 10882.351 11404.107 10796.227 10298.889 10334.90     .        .   
[4,]  8640.314 10050.090 10882.351 18148.897 17375.807 17253.214 16534.326 16024.77 15857.49     .   
[5,]  8605.827 10342.738 11404.107 17375.807 20466.958 19981.791 19346.404 18760.51 18020.95 18298.81
[6,]  8155.825  9789.687 10796.227 17253.214 19981.791 20661.165 19819.846 19012.31 18091.70 18532.20
[7,]     .      9428.786 10298.889 16534.326 19346.404 19819.846 20904.765 19508.50 18314.66 18368.21
[8,]     .         .     10334.903 16024.771 18760.505 19012.313 19508.504 19421.92 17761.69 17746.44
[9,]     .         .         .     15857.494 18020.948 18091.702 18314.663 17761.69 17595.06 17049.94
[10,]     .         .         .         .     18298.811 18532.200 18368.212 17746.44 17049.94 18031.37
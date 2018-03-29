##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
df = read_feather("../data/features_106_wide.feather")

## focus on year >= 2000 for test purpose
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA
#hist(mat)

colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)

#####################################
#### visualize the original data ####
#####################################
##ss = mat2stack(mat, 31, seq(1,100, by=5))
#levelplot(ss, par.settings = colthm)

##########################################################
#### gapfill based on  pixel level temporal smoothing ####
##########################################################
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
res1 = gapfill(year, doy, mat, 31,31, h = 1, doyrange = 1:365, nnr=1)
# res1 = gapfill(year, doy, mat, 31,31, h = 0, doyrange = 1:365, nnr=5)
## temporal trend visulization
ssmean1 = mat2stack(res1$temporal.mean, 31)
## pdf("output/ssmean1_plot.pdf")
## for(i in 1:365)
##     print(levelplot(ssmean1[[i]], par.settings = colthm))
## dev.off()

############################
##### Simulation study #####
############################
##### Create artificial missing partens and apply to mat
res = res1
set.seed(20180124)
## number of artifical missing images to make
n = 6
fidx = sample(res$idx$idx.fullyobserved, n) ## full observed image index
pidx = sample(res$idx$idx.partialmissing, n) ## partial observed image index
## The fully observed images
fmat = mat[fidx[1:n], ]
## saveRDS(fmat, "output/fmat.rds")
## apply missing patterns to fully observed images
for(i in 1:n){
  mat[fidx[i],][is.na(mat[pidx[i],])] = NA
}
## artificial partial missing images
pmat = mat[fidx[1:n], ]
## saveRDS(pmat, "output/pmat.rds")

#########################################################
#### gapfill based on pixel level temporal smoothing ####
#########################################################
## number of nnr to use
k = 10
Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
## Gapfill
res = gapfill(year, doy, mat, 31,31, h=1, doyrange=1:365, nnr=k, method = "lc")
## imputed images
imat = res$imputed.mat[fidx[1:n],]
## saveRDS(imat, "output/imat.rds")
## Calculate MSE
MSE = apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
MSE
### MSE using "lc" method, smooth_spline, k = 1
## [1] 24882.475 18531.305  2963.482  8807.417 20024.720 30655.182
### MSE using "lc" method, smooth_spline, k = 10
## [1]  9623.177 13413.383  6728.366 11174.873 14530.421  4111.339

###########################################################
#### gapfill based on cluster level temporal smoothing ####
###########################################################
Gapfill::opts$set(temporal_mean_est = Gapfill::llreg)
cmat = t(mat[setdiff(res1$idx$idx.fullyobserved, fidx),])

MSE1 = foreach(ncl= seq(5, 80, by = 5), .combine = "rbind") %dopar%{
    ## clustering analysis
    cres = kmeans(cmat, centers = ncl, nstart = 10)
    res = gapfill(year, doy, mat, 31,31, h=1, doyrange=1:365, nnr=k, method = "lc", cluster = cres$cluster)
    ## imputed images
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output/imat_ncl_", ncl, "_k_10.rds")) 
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}
MSE1

## msedf = as.data.frame(rbind(c(0, MSE), cbind(seq(5, 80, by = 5), MSE1)))
## colnames(msedf) = c("ncl", paste0("s", 1:6))
## write.csv(msedf, "output/msedf_llreg_k_10.csv", row.names=FALSE)

## MSE using year >= 2000; llreg; k = 1
##               [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
## result.1  16198.73 15421.72 14297.70 13725.171 20696.79 45858.71
## result.2  12232.78 15548.45 12619.47 11118.220 19845.20 49448.34
## result.3  13515.63 15042.21 12925.24  9810.496 16841.64 47682.13
## result.4  14202.37 17496.78 10880.10  8066.347 17831.96 53301.16
## result.5  11821.93 17657.36 14770.08 10282.297 18375.04 57667.79
## result.6  15533.89 16442.65 14740.61 10440.577 16482.80 36896.87
## result.7  13541.47 18427.80 14541.41  9194.375 18379.44 58051.50
## result.8  15386.52 16922.65 11922.80 11073.048 15975.21 62333.80
## result.9  18948.07 17489.70 11756.86 10764.650 15949.45 65701.47
## result.10 26666.83 17562.36 13037.14 11001.421 16413.72 63856.89
## result.11 24350.91 18774.96 15797.43  9957.700 19706.25 44648.64
## result.12 15134.89 19595.54 12680.55 11262.802 19022.66 48645.86
## result.13 24385.38 21299.37 14418.77 10899.599 16947.58 47805.34
## result.14 24000.60 17005.08 14368.29 11255.277 19587.49 37590.54
## result.15 26208.88 15470.06 16508.67  9928.019 19007.25 58919.29
## result.16 14755.05 20108.27 14349.27 10716.696 16601.38 49601.04
## MSE using year >= 2000; llreg; k = 10
##                [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
## result.1   8245.605 11773.97 8630.903 11918.396 30573.28 5264.601
## result.2   6463.642 14814.15 6945.797 10142.214 20919.51 5424.867
## result.3   6827.795 13511.04 9010.810  7925.588 19219.88 4806.355
## result.4   5570.625 14381.02 8134.396  7195.371 18111.47 4129.821
## result.5   5136.005 14001.08 9591.737  8042.173 15173.31 4619.273
## result.6   5401.299 13300.15 8602.546  7585.966 17852.89 4546.469
## result.7   5436.913 13979.15 8560.096  7725.382 17324.55 5011.835
## result.8   4942.892 14052.11 7510.976  7145.756 18639.90 4844.128
## result.9   6208.032 16930.46 8216.332  7385.954 19420.55 4341.705
## result.10  5556.950 17006.35 8473.495  6974.657 17744.40 4186.392
## result.11  6181.524 17643.31 9160.925  7701.398 17023.22 5237.194
## result.12 15882.313 17911.49 8468.834  6721.060 19086.58 4917.100
## result.13 14080.411 17274.81 9387.145  7148.523 20212.72 4405.456
## result.14  7269.319 17844.98 7886.821  7285.983 16799.77 5097.575
## result.15 16133.326 16155.61 9128.313  6865.969 14781.29 5014.883
## result.16 19519.816 15473.38 9051.305  6645.482 20147.97 6773.992

Gapfill::opts$set(temporal_mean_est = Gapfill::lpreg)
cmat = t(mat[setdiff(res1$idx$idx.fullyobserved, fidx),])

MSE2 = foreach(ncl= seq(5, 80, by = 5), .combine = "rbind") %dopar%{
    ## clustering analysis
    cres = kmeans(cmat, centers = ncl, nstart = 10)
    res = gapfill(year, doy, mat, 31,31, h=1, doyrange=1:365, nnr=k, method = "lc", cluster = cres$cluster)
    ## imputed images
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output/imat_ncl_", ncl, "_k_10.rds")) 
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}
MSE2
msedf = as.data.frame(rbind(c(0, MSE), cbind(seq(5, 80, by = 5), MSE2)))
colnames(msedf) = c("ncl", paste0("s", 1:6))
write.csv(msedf, "output/msedf_lpreg_k_1.csv", row.names=FALSE)
## MSE using year >= 2000; lpreg; k = 1
##               [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
## result.1  26915.42 16721.58 6362.292 12424.341 19371.28 13777.46
## result.2  24628.82 17350.46 5907.520 11259.212 17687.16 10812.05
## result.3  24515.02 17156.24 4772.492  9470.697 16352.59 10317.59
## result.4  26189.44 15372.46 5098.637  9555.741 15011.25 10328.49
## result.5  24736.53 15303.66 4702.539 10241.526 15271.31 10070.49
## result.6  27411.38 15603.18 4530.325  9310.303 15174.77 10864.95
## result.7  25852.92 15210.73 4507.911  9807.086 15168.26 10272.92
## result.8  26890.62 15348.49 4724.754  9835.321 15021.95 10232.99
## result.9  25784.51 15140.32 4387.357 10152.568 15391.72 10207.29
## result.10 26004.51 15207.12 4432.073  9996.293 14595.89 10596.54
## result.11 25889.90 15309.16 4693.095  9204.485 15309.96 10057.13
## result.12 26147.65 14936.75 4437.108 10084.463 15148.05 10589.86
## result.13 26850.40 15930.85 4392.584  9307.836 14830.13 10684.54
## result.14 26328.06 15235.44 4547.591  9637.563 15228.11 10460.59
## result.15 27827.98 15494.99 4651.199  9940.806 14935.52 19273.61
## result.16 26754.46 15671.74 4566.384  9436.868 14919.27 10554.95
## MSE using year >= 2000; lpreg; k = 10
##               [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
## result.1  14144.77 13492.77 8016.635 15544.27 27972.95 4422.302
## result.2  11172.00 15856.46 6283.631 13899.09 20114.64 4930.919
## result.3  10202.10 15777.54 8045.401 16618.11 20759.18 4448.222
## result.4  10114.02 15665.71 8510.993 16081.43 15265.45 4760.123
## result.5  11031.10 15809.42 7463.264 15451.71 19150.57 4368.442
## result.6  10913.47 14868.26 7358.232 16035.41 14398.58 4484.235
## result.7  11331.96 14894.43 7177.993 14672.71 14455.08 4272.042
## result.8  10928.55 14269.26 7867.067 13552.75 15149.59 4209.243
## result.9  10257.07 14756.46 7590.749 15660.52 13372.90 4409.658
## result.10 10182.44 15089.98 7689.482 16401.12 13539.27 3843.292
## result.11 10875.35 14685.77 7656.356 14794.61 14715.45 4120.988
## result.12 10999.33 14012.93 7616.298 14949.82 12026.35 4058.814
## result.13 10372.23 14504.94 7810.006 14477.89 13721.53 4089.087
## result.14 10858.61 15116.26 7897.648 15601.98 13279.76 4225.349
## result.15 10266.58 14328.45 7104.950 14176.69 12283.81 4159.677
## result.16 10234.69 14466.30 7316.254 14647.33 12651.14 3934.677


Gapfill::opts$set(temporal_mean_est = Gapfill::smooth_spline)
cmat = t(mat[setdiff(res1$idx$idx.fullyobserved, fidx),])

MSE3 = foreach(ncl= seq(5, 80, by = 5), .combine = "rbind") %dopar%{
    ## clustering analysis
    cres = kmeans(cmat, centers = ncl, nstart = 10)
    res = gapfill(year, doy, mat, 31,31, h=1, doyrange=1:365, nnr=k, method = "lc", cluster = cres$cluster)
    ## imputed images
    imat = res$imputed.mat[fidx[1:n],]
    saveRDS(imat, paste0("output/smooth_spline_k_10/imat_ncl_", ncl, "_k_10.rds")) 
    apply((fmat - imat)^2, 1, sum) / apply(pmat, 1, function(x) sum(is.na(x)))
}
MSE3
msedf = as.data.frame(rbind(c(0, MSE), cbind(seq(5, 80, by = 5), MSE3)))
colnames(msedf) = c("ncl", paste0("s", 1:6))
write.csv(msedf, "output/smooth_spline_k_10/msedf_ss_k_10.csv", row.names=FALSE)
## MSE using year >= 2000; smooth_spline; k = 1
##               [,1]     [,2]     [,3]      [,4]     [,5]     [,6]
## result.1  16294.94 15793.05 14410.17 13367.879 20801.68 46187.78
## result.2  12812.87 16248.97 12276.39 10131.361 18891.85 46081.97
## result.3  12761.75 16413.15 12521.37 10023.702 17424.02 53862.87
## result.4  11652.35 15619.10 11256.23 11480.706 17702.38 53120.88
## result.5  12610.36 15291.17 15230.47 10369.332 17762.75 61107.24
## result.6  13229.49 17207.23 14329.77 10311.944 18036.95 57604.15
## result.7  15166.84 18743.82 11548.54 10638.030 16278.47 60928.71
## result.8  11437.18 16588.38 11626.76  9944.384 17409.88 52697.77
## result.9  16523.16 15823.38 16573.77 10074.957 18520.42 45311.26
## result.10 16580.91 16571.05 15125.51 11199.353 18159.43 50610.08
## result.11 23968.31 16244.28 13788.51  9374.754 16743.87 49884.80
## result.12 16239.11 15620.62 15955.72  9657.683 17850.06 38894.76
## result.13 24423.10 17976.41 14315.97 10645.051 15696.90 48006.20
## result.14 27837.98 15514.41 15956.17 10023.294 19018.29 54295.49
## result.15 26270.80 16459.49 17002.55  9914.851 19316.37 52147.39
## result.16 24010.43 21253.19 13309.60 10910.295 16148.64 44669.56

## MSE using year >= 2000; smooth_spline; k = 10
##                [,1]     [,2]     [,3]      [,4]     [,5]      [,6]
## result.1   8051.682 14680.58 12553.70 14084.549 27803.14  7883.912
## result.2   7743.998 16227.43 13975.71 12421.600 27553.16  6071.320
## result.3   7061.001 20615.36 15155.63 11631.766 19974.29 10667.085
## result.4   8211.807 17496.25 10435.31  9769.020 19565.07 10228.898
## result.5   7670.967 20478.73 12188.25  8160.009 14096.84 12062.802
## result.6   8240.559 14375.27 17742.19 10004.260 14720.93 16000.769
## result.7   7809.395 15145.42 14124.34  8816.226 14367.78 15302.313
## result.8   8171.820 19337.65 11471.32 10414.426 13414.94 12792.248
## result.9   7410.336 20738.42 11846.79 11153.152 14864.81 14869.832
## result.10 11897.330 20155.25 11207.27  9304.357 12486.34 13397.398
## result.11 11017.602 18529.60 18812.77 10433.722 13178.63 16289.656
## result.12  9182.264 17794.50 16244.10 11836.222 13304.84 15660.891
## result.13 10591.497 21979.86 13308.61 10568.872 12657.48 17247.019
## result.14  9054.259 19784.45 19963.35 12004.338 14718.60 16067.532
## result.15  6966.334 23032.39 13756.41 11473.244 12600.67 12679.259
## result.16  7885.073 20953.57 16040.43 12432.591 13229.86 16958.191


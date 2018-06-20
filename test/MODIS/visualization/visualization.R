##### Test for Landsat data
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
dat0 = readRDS("../../data/MYD11A1Day2010.rds")
colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)
i = 100
r = raster(matrix(dat0[100,], 1200, byrow = TRUE))
pdf(paste0("year2010_doy_", i, ".pdf"))
print(levelplot(r, par.settings = colthm, main = paste0("Landsat scene p026r031 Site 106, ", i),
                margin = FALSE))
dev.off()

## stackRaster
r.list = list()
dd = seq(1,365, 30)
for(i in 1:length(dd)){
  r.list[[i]] = raster(matrix(dat0[dd[i],], 1200, byrow = TRUE))
}
s = stack(r.list)
pdf(paste0("MODIS_2010.pdf"))
print(levelplot(s, par.settings = colthm, names.attr=as.character(dd), main = "MODIS MYD11A1Day 2010",
                layout=c(5,3)))
dev.off()

## MYD VS MOD
dat0 = readRDS("../../data/MYD11A1Day2010.rds")
dat1 = readRDS("../../data/MOD11A1Day2010.rds")
msk = getMask(dat1)
dat0[,msk] = NA
i = 100
r0 = raster(matrix(dat0[i,], 1200, byrow = TRUE))
r1 = raster(matrix(dat1[i,], 1200, byrow = TRUE))
s = stack(r0, r1)
pdf("year2010_doy_100_MYD_vs_MOD.pdf")
print(levelplot(s, names.attr=c("MYD", "MOD"), main = paste0("MODIS 2010 DOY 100 MYD vs. MOD")))
dev.off()

## msk
pdf(paste0("MODIS_mask.pdf"))
print(levelplot(raster(matrix(msk, 1200, byrow=TRUE)), main = paste0("MODIS mask"),
                margin = FALSE))
dev.off()
##### Test for Landsat data
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(Gapfill)
df = read_feather("../data/features_106_wide.feather")
df = df %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat = as.matrix(df[,-c(1:2)])
mat[mat > 2000] = NA

colthm = RdBuTheme()
colthm$regions$col = rev(colthm$regions$col)
## focus on year >= 2000 for test purpose

########## bin days ############
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
## yearuni = yearuni[seq(1, length(yearuni), by = 2)]
doybinuni = sort(unique(doybin))
N = length(yearuni)

mat4plot = matrix(NA, length(doybinuni)*N, ncol(mat))
for(j in 1:length(doybinuni)){
  for(i in 1:N){
    idx = year == yearuni[i] & doybin == doybinuni[j]
    if(sum(idx) == 1)
      mat4plot[(j-1)*N + i, ] = mat[year == yearuni[i] & doybin == doybinuni[j],] else
        if(sum(idx) > 1)
          warning("Multiple matches.")
  }
}
str(mat4plot)
## gapfill::Image(datarray[,,1:16, 1:8])
s = mat2stack(mat4plot, 31)

# pdf(paste0("output/data_overview_16_le1200.pdf"))
# for(i in 1:7){
#   print(rasterVis::levelplot(s[[seq((i-1)*6*N+1, i*6*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 6)),
#                        layout = c(N,6), at = seq(0,1200, 20)))
# }
# print(rasterVis::levelplot(s[[(7*6*N+1):(length(doybinuni)*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 4)),
#                            layout = c(N,6), at = seq(0,1200, 20)))
# dev.off()
# 
# pdf(paste0("output/data_overview_16_le1200.pdf"))
# for(i in 1:4){
#   print(rasterVis::levelplot(s[[seq((i-1)*10*N+1, i*10*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 10)),
#                              layout = c(N,10), at = seq(0,1200, 20)))
# }
# print(rasterVis::levelplot(s[[(i*10*N+1):(length(doybinuni)*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 6)),
#                            layout = c(N,10), at = seq(0,1200, 20)))
# dev.off()

pdf(paste0("output/data_overview_four_seasons_le1200.pdf"))
print(rasterVis::levelplot(s[[seq(1, 12*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 12)),
                             layout = c(N,12), at = seq(0,1200, 20)))
print(rasterVis::levelplot(s[[seq(12*N+1, 23*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 11)),
                           layout = c(N,11), at = seq(0,1200, 20)))
print(rasterVis::levelplot(s[[seq(23*N+1, 34*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 11)),
                           layout = c(N,11), at = seq(0,1200, 20)))
print(rasterVis::levelplot(s[[(34*N+1):(length(doybinuni)*N)]], par.settings = colthm, names.attr=as.character(rep(yearuni, 12)),
                           layout = c(N,12), at = seq(0,1200, 20)))
dev.off()
# 
# for(yy in c(1995,2000,2010, 2015)){
#   tmpdf = df[df$year == yy,]
#   year = tmpdf$year
#   doy = tmpdf$doy
#   mat = as.matrix(tmpdf[,-c(1:2)])
#   mat[mat > 2000] = NA
#   r.list = list()
#   for(i in 1:nrow(mat)){
#     r.list[[i]] = raster(matrix(mat[i,], 31))
#   }
#   s = stack(r.list)
#   pdf(paste0("year", yy, ".pdf"))
#   print(levelplot(s, par.settings = colthm, names.attr=as.character(doy), main = paste0("Landsat scene p026r031 Site 106, ", yy),
#             layout=c(9,5)))
#   dev.off()
# }

## Based on the images, five fully observed images are chosen from each season: doy 1~96; 97~184; 185~273; 274~365
fidx1 = c(261, 265, 312, 387, 581)
fidx2 = c(145, 276, 444, 481, 587)
fidx3 = c(198, 202, 493, 549, 557)
fidx4 = c(82, 293, 505, 609, 615)
s = mat2stack(mat[c(fidx1, fidx2, fidx3, fidx4),], 31)
levelplot(s, par.settings = colthm, layout = c(5,4))



